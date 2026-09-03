#!/usr/bin/env node
/**
 * 宣言インデックス —— 在庫検索の対策(2026-08-21、2026-09-03 に拡張)
 *
 * ★動機(実測): 宣言が 9,937 個になり、「使える補題がすでにあるか」を探すのに
 * 木全体(198,000 行)を `grep` する往復が増えていた。
 * ★1 本のテキストにしておけば **1 回の `grep` で済む**。
 *
 * ★★2026-09-03 の拡張(案 I / 案 D)。動機は「事後的に既存が判明した」15 件の実測:
 *   - 自前の在庫の取りこぼし 11 件(73%)。索引はあったのに**名前しか入っていない**ので、
 *     `tools/lean-idioms.md` 自身の指針「探すべきは道具の名前ではなく**結論の形**」
 *     「結論に現れるリテラル(`1728`)で grep する」が**機械的に実行できなかった**。
 *     → 案 I: **statement(型)の本文を索引に入れる**。
 *   - mathlib の取りこぼし 4 件(27%)。`mathlib` は索引に 1 行も入っていなかった
 *     (「mathlib は ℘ を丸ごと持っていた」第 597 / 「EGA IV §8 を持っていた」など)。
 *     → 案 D: **mathlib も索引する**(別ファイル。ABC3 の grep を遅くしないため)。
 *
 * 出力(既定 `.cache/`、gitignore 済み):
 *   `decl-index.txt`    1 行 1 宣言 —— `種別 <TAB> 名前 <TAB> ファイル:行 <TAB> statement`
 *   `src-index.txt`     1 行 1 locator —— `論文 pPDF頁 <TAB> item <TAB> ファイル:行 <TAB> 宣言`
 *   `mathlib-index.txt` 同上(mathlib)。★別ファイルなのは、ABC3 だけを引きたいときに
 *                       10 倍の行を舐めさせないためである。
 *
 * 使い方:
 *   node tools/decl-index.mjs            # ABC3 のみ(既定、速い)
 *   node tools/decl-index.mjs --mathlib  # mathlib も作り直す(数十秒)
 *   grep -i "rootMap" .cache/decl-index.txt
 *   grep "1728" .cache/decl-index.txt          # ★結論のリテラルで引く(案 I)
 *   grep -i "weierstrass" .cache/mathlib-index.txt   # ★mathlib の在庫を引く(案 D)
 *   grep "Proposition 5.5" .cache/src-index.txt
 */

import { readFileSync, readdirSync, statSync, mkdirSync, writeFileSync, existsSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const LEAN_SRC = join(ROOT, 'lean', 'ABC3');
const MATHLIB_SRC = join(ROOT, 'lean', '.lake', 'packages', 'mathlib', 'Mathlib');
const OUT_DIR = join(ROOT, '.cache');
const WANT_MATHLIB = process.argv.includes('--mathlib');

function walk(dir, acc = []) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else if (p.endsWith('.lean')) acc.push(p);
  }
  return acc;
}

const DECL_RE =
  /^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(theorem|lemma|def|abbrev|instance|structure|inductive|class)\s+([A-Za-z_][\w'!?₀-₉]*(?:\.[\w'!?₀-₉]+)*)/;
const NS_RE = /^\s*namespace\s+([A-Za-z_][\w'.₀-₉]*)/;
const END_RE = /^\s*end\s+([A-Za-z_][\w'.₀-₉]*)\s*$/;

/** 宣言の statement(型)を 1 行に畳んで返す。
 *
 * ★どこまでを statement とするか: 宣言行から `:=` / `where` / 次の宣言 の手前まで。
 *   証明本体は入れない——入れると索引が本文と同じ大きさになり、`grep` の意味が消える。
 * ★上限 400 文字。長い型は先頭だけで十分に効く(結論のリテラルは前の方に出る)。
 */
const STATEMENT_MAX = 400;
function statementOf(lines, start) {
  const buf = [];
  for (let j = start; j < Math.min(start + 24, lines.length); j++) {
    let s = lines[j];
    if (j > start && (DECL_RE.test(s) || NS_RE.test(s) || END_RE.test(s))) break;
    // 行コメントを落とす(`--` の前だけ残す)
    s = s.replace(/--.*$/, ' ');
    const cut = s.search(/:=|\bwhere\b/);
    if (cut >= 0) { buf.push(s.slice(0, cut)); break; }
    buf.push(s);
    if (buf.join(' ').length > STATEMENT_MAX * 2) break;
  }
  return buf.join(' ').replace(/\s+/g, ' ').trim().slice(0, STATEMENT_MAX);
}

/** 1 本の木を舐めて宣言と locator を集める。 */
function scan(srcRoot, { collectSrc }) {
  const decls = [];
  const srcs = [];
  for (const file of walk(srcRoot)) {
    const rel = relative(srcRoot, file).replace(/\\/g, '/');
    const lines = readFileSync(file, 'utf8').split(/\r?\n/);
    const nsStack = [];
    let pendingSrcDecl = null;
    for (let i = 0; i < lines.length; i++) {
      const line = lines[i];
      const ns = NS_RE.exec(line);
      if (ns) { nsStack.push(ns[1]); continue; }
      const en = END_RE.exec(line);
      if (en && nsStack.length && nsStack[nsStack.length - 1] === en[1]) { nsStack.pop(); continue; }
      const m = DECL_RE.exec(line);
      if (!m) continue;
      const full = [...nsStack, m[2]].join('.');
      // ★4 列目が statement(案 I)。名前で引いて外した失敗形を潰すための列である。
      decls.push(`${m[1].padEnd(9)}\t${full}\t${rel}:${i + 1}\t${statementOf(lines, i)}`);
      if (!collectSrc) continue;
      if (m[2].endsWith('.src')) pendingSrcDecl = { name: full, rel, line: i + 1, at: i };
      // `.src` の中身(paper / pdfPage / item)を後続 6 行から拾う
      if (pendingSrcDecl && i - pendingSrcDecl.at <= 6) {
        const blob = lines.slice(pendingSrcDecl.at, Math.min(pendingSrcDecl.at + 7, lines.length)).join(' ');
        const paper = /paper\s*:=\s*"([^"]*)"/.exec(blob);
        const page = /pdfPage\s*:=\s*(\d+)/.exec(blob);
        const item = /item\s*:=\s*"([^"]*)"/.exec(blob);
        if (paper && item) {
          srcs.push(`${paper[1]} p${page ? page[1] : '?'}\t${item[1]}\t${pendingSrcDecl.rel}:${pendingSrcDecl.line}\t${pendingSrcDecl.name}`);
          pendingSrcDecl = null;
        }
      }
    }
  }
  decls.sort();
  srcs.sort();
  return { decls, srcs };
}

mkdirSync(OUT_DIR, { recursive: true });

const abc3 = scan(LEAN_SRC, { collectSrc: true });
writeFileSync(join(OUT_DIR, 'decl-index.txt'), abc3.decls.join('\n') + '\n', 'utf8');
writeFileSync(join(OUT_DIR, 'src-index.txt'), abc3.srcs.join('\n') + '\n', 'utf8');
console.log(`decl-index.txt: ${abc3.decls.length} 宣言(statement つき)`);
console.log(`src-index.txt:  ${abc3.srcs.length} locator`);

const mlPath = join(OUT_DIR, 'mathlib-index.txt');
if (WANT_MATHLIB) {
  if (!existsSync(MATHLIB_SRC)) {
    console.log(`mathlib-index.txt: ${MATHLIB_SRC} が無いので作らなかった`);
  } else {
    const ml = scan(MATHLIB_SRC, { collectSrc: false });
    writeFileSync(mlPath, ml.decls.join('\n') + '\n', 'utf8');
    console.log(`mathlib-index.txt: ${ml.decls.length} 宣言`);
  }
} else if (existsSync(mlPath)) {
  console.log('mathlib-index.txt: 既にある(作り直すなら --mathlib)');
} else {
  console.log('mathlib-index.txt: ★まだ無い。`node tools/decl-index.mjs --mathlib` で作る');
}
