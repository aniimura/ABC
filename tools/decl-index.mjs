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
 *   grep "CompactSpace" .cache/mathlib-index.txt     # ★無名 instance も入る(2026-09-05)
 *
 * ★★2026-09-05 の拡張(backlog M5)。動機は「mathlib に無い」の誤判定 4 件のうち 1 件が
 *   **索引の穴**だったこと: `instance : CompactSpace Gal(K/k)` は名前を持たないので
 *   索引に 1 行も入らず、名前で grep しても出なかった。mathlib の `instance` 行の
 *   **58.8%(15,798 本)がこれ**である。→ 無名 instance を**型を鍵として**索引に入れる。
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

const MODS = '(?:private\\s+|protected\\s+|noncomputable\\s+|scoped\\s+|local\\s+)*';
/** ★★★★★**名前の文字類**(2026-09-06 に拡張。第 1036)。
 *
 * ★旧版は `[A-Za-z_][\w'!?₀-₉]*` だったが、JS の `\w` は **ASCII だけ**である。
 * そのため **名前が非 ASCII 文字で切れていた**。実測(2026-09-06):
 *   `WeierstrassCurve.preΨ₄` → 索引上 `WeierstrassCurve.pre`
 *   `coeff_Ψ`               → 索引上 `coeff_`
 * したがって `Ψ|Φ` で名前欄を grep すると **0 件**と出るが、
 * 実体は存在する(statement 欄には 133 行あった)。
 *
 * ★★これはこの木で 5 件出ている「不在の誤判定」と**同じ回路**である
 * (`ULift.field` / `continuousCohomology` / `Ẑ` / `CompactSpace Gal` / `IsArithFrobAt`)。
 * 無名 instance を入れた 2026-09-05 の拡張(backlog M5)と同じ動機である。
 *
 * ★入れた範囲(Lean の識別子で実際に使われるものだけ):
 *   `\u00C0-\u024F` Latin Extended / `\u0370-\u03FF` ギリシャ(Ψ Φ σ μ π …)
 *   `\u1D00-\u1DBF` 音声拡張(ᵥ ᵢ …) / `\u2070-\u209F` 上付き・下付き
 *   `\u2100-\u214F` 文字風記号(ℓ ℝ ℤ ℕ ℚ ℂ ℘ …) / `\uD800-\uDFFF` 代理対(𝓞 𝔓 𝕜 …)
 * ★★**矢印(`\u2190-\u21FF`)や論理記号(`\u2200-`)は入れていない** ——
 * 入れると `def f : A → B` の型まで名前に食い込むからである。 */
const IDX = "\\w'!?\\u00C0-\\u024F\\u0370-\\u03FF\\u1D00-\\u1DBF\\u2070-\\u209F\\u2100-\\u214F\\uD800-\\uDFFF";
const ID1 = "A-Za-z_\\u00C0-\\u024F\\u0370-\\u03FF\\u2100-\\u214F\\uD800-\\uDFFF";
const NAME = `[${ID1}][${IDX}]*(?:\\.[${IDX}]+)*`;
const DECL_RE = new RegExp(
  `^\\s*(?:@\\[[^\\]]*\\]\\s*)?${MODS}(theorem|lemma|def|abbrev|instance|structure|inductive|class)\\s+(${NAME})`);

/** ★**無名 instance**(2026-09-05、メタ第 2 回 / backlog M5)。
 *
 * `DECL_RE` は種別の後ろに**名前を要求する**ので、`instance : CompactSpace G := …` のように
 * 名前を持たない instance は索引に 1 行も入っていなかった。
 * 実測: mathlib の `instance` 行 26,846 本のうち **15,798 本(58.8%)が無名**である。
 *
 * ★これが「mathlib に X は無い」の誤判定を生んでいた(backlog M4 の 4 件目):
 * `CompactSpace Gal(Ω/F)` は `FieldTheory/Galois/Profinite.lean:329` に
 * `instance [IsGalois k K] : CompactSpace Gal(K/k)` として**在る**のに、
 * 名前が無いので索引を grep しても出ず、ソース直読で初めて見つかった。
 *
 * 無名なので**鍵は型**である——名前の列には `⟨無名⟩` を置き、
 * 4 列目の statement(型)で引けるようにする。
 * `instance (priority := low) Foo.bar : …` のように名前の前に括弧が来る形も
 * ここで拾う(`DECL_RE` は括弧で止まるため)。 */
const ANON_INST_RE = new RegExp(`^\\s*(?:@\\[[^\\]]*\\]\\s*)?${MODS}instance\\b(.*)$`);
const INST_NAME_RE = new RegExp(`^\\s*(?:\\(\\s*priority\\s*:=[^)]*\\)\\s*)?(${NAME})(?=\\s|:|$)`);
const ANON = '⟨無名⟩';

const NS_RE = /^\s*namespace\s+([A-Za-z_][\w'.₀-₉]*)/;
const END_RE = /^\s*end\s+([A-Za-z_][\w'.₀-₉]*)\s*(?:--.*)?$/;
const SEC_RE = /^\s*section\b\s*([A-Za-z_][\w'.₀-₉]*)?\s*(?:--.*)?$/;
const END_BARE_RE = /^\s*end\s*(?:--.*)?$/;

/** 開いている scope の積み。**成分 1 つ = 1 段**で積む。
 *
 * ★★2026-09-06(backlog M22 の後半)。旧実装は `namespace A.B.C` を**1 段**で積み、
 * `end NAME` は**頂と完全一致**したときだけ落としていた。しかし Lean は
 * ドット付きの名前空間を**部分的に閉じてよい**:
 *
 *   namespace Order.Frame.MinimalAxioms   -- 旧: 1 段 "Order.Frame.MinimalAxioms"
 *   end MinimalAxioms                     -- 旧: 一致しないので落ちない ★
 *   end Order.Frame                       -- 旧: 一致しないので落ちない ★
 *
 * その結果、後続の宣言に**閉じたはずの名前空間が被り続ける**。実測(2026-09-06):
 * `Order/CompleteBooleanAlgebra.lean:913` の 1 件は索引上
 * `Order.Frame.MinimalAxioms.Order.Coframe.MinimalAxioms.CompleteDistribLattice.
 *  MinimalAxioms.CompletelyDistribLattice.MinimalAxioms.Equiv.coframe` と 9 段になっていた
 * (実体は `Equiv.coframe`)。
 *
 * ★台帳 M22 の見立て(「`end Foo -- コメント` を pop できないため」)は**外れ**である。
 * `end NAME -- コメント` は mathlib 全体で **45 行**しかなく、
 * ドット付き namespace は **1,486 本**ある。数えて分かった。
 *
 * ★`section` も積む。積まないと `end 名前つき section` が
 * 同名の名前空間成分を**誤って落とす**(旧実装にあった危険)。
 * section は名前に寄与しないので `ns:false` で持つ。 */
const pushNamespace = (stack, name) => { for (const c of name.split('.')) stack.push({ ns: true, name: c }); };
const popEnd = (stack, name) => {
  if (name === null) { // 裸の `end` = 無名 section を閉じる
    const top = stack[stack.length - 1];
    if (top && !top.ns && top.name === null) stack.pop();
    return;
  }
  const segs = name.split('.');
  if (stack.length >= segs.length) { // 末尾 segs.length 段が成分ごとに一致すれば落とす
    let ok = true;
    for (let k = 0; k < segs.length; k++) if (stack[stack.length - segs.length + k].name !== segs[k]) { ok = false; break; }
    if (ok) { stack.length -= segs.length; return; }
  }
  const top = stack[stack.length - 1]; // ドット付きの名前を持つ section を 1 段で閉じる形
  if (top && top.name === name) stack.pop();
};
const nsPath = (stack) => stack.filter((e) => e.ns).map((e) => e.name);

/** 宣言名に名前空間を被せる。
 *
 * ★★2026-09-06(backlog M22)。**`_root_.` は名前空間を抜ける指示**であって
 * 名前の一部ではない。無条件に `[...nsStack, name].join('.')` すると
 *   `namespace Algebra` の中の `theorem _root_.RingHom.commSemiringToCommRing` が
 *   索引上 `Algebra._root_.RingHom.commSemiringToCommRing` になる。
 * 実測(2026-09-06): `.cache/mathlib-index.txt` の名前欄に `._root_.` が **3,041 件**、
 * 先頭が `_root_.` のものが 52 件(ABC3 側は 0 件——この書き方をしていない)。
 * ★M14(名前欄が非 ASCII で切れる 24,659 件)・M20(statement 欄が autoParam で切れる
 * 6,271 件)と**同じ回路の 3 例目**である——「mathlib に無い」の誤判定を生む。 */
const qualify = (nsStack, name) =>
  (name.startsWith('_root_.') ? name.slice(7) : [...nsPath(nsStack), name].join('.'));

/** 文字列リテラルを**長さを保ったまま**空白へ伏せる。
 *
 * ★長さが 1 文字も変わらないのが要点である。伏せた側で見つけた**添字を
 * そのまま原文の `slice` に使える**ので、切り出しの位置がずれない。
 * `tools/check.mjs` の `maskLeanStrings` と同じ型(backlog M15)。 */
const maskStrings = (s) => s.replace(/"(?:[^"\\\n]|\\.)*"/g, (m) => ' '.repeat(m.length));

/** 宣言の statement(型)を 1 行に畳んで返す。
 *
 * ★どこまでを statement とするか: 宣言行から `:=` / `where` / 次の宣言 の手前まで。
 *   証明本体は入れない——入れると索引が本文と同じ大きさになり、`grep` の意味が消える。
 * ★上限 400 文字。長い型は先頭だけで十分に効く(結論のリテラルは前の方に出る)。
 *
 * ★★★2026-09-06(backlog M15/M20)。**素朴な `search(/:=|\bwhere\b/)` は結論を壊す**。
 * 壊れ方は 2 つあり、どちらも「構文で切る前に文字列と括弧を見ていない」ことに由来する:
 *
 *   1. **括弧の中の `:=`**。名前つき引数 `f (p := p)`・`(priority := 100)`・
 *      `{ x := 0 }` で切れてしまう。第 7 回(M17)が `unwired.mjs` を書いたとき
 *      **5 件で結論が壊れる**と実測した。→ **深さ 0 の `:=` でだけ切る**。
 *   2. **文字列の中の `:=` / `where` / `--`**。docstring ではなく
 *      `def foo.src := { paper := "…" }` のような**台帳の本文**に出る。
 *      → **伏せてから切る**(`maskStrings`)。
 *
 * ★深さは `( [ { ⟨ ⦃` で数える。行をまたいで持ち越す(複数行の署名があるため)。
 * ★`--` の行コメント落としも**伏せた側**で探す。伏せないと
 *   文字列に `--` を含む行が途中で切れる。
 */
const STATEMENT_MAX = 400;
const OPEN = '([{⟨⦃';
const CLOSE = ')]}⟩⦄';
function statementOf(lines, start) {
  const buf = [];
  let depth = 0;
  for (let j = start; j < Math.min(start + 24, lines.length); j++) {
    const raw = lines[j];
    if (j > start && (DECL_RE.test(raw) || ANON_INST_RE.test(raw) || NS_RE.test(raw) || END_RE.test(raw))) break;
    // 行コメントを落とす(`--` の前だけ残す)。★伏せた側で探す。
    let s = raw;
    let m = maskStrings(raw);
    const c = m.indexOf('--');
    if (c >= 0) { s = `${s.slice(0, c)} `; m = `${m.slice(0, c)} `; }
    // ★深さ 0 の `:=` / `where` でだけ切る。
    let cut = -1;
    for (let k = 0; k < m.length; k++) {
      const ch = m[k];
      if (OPEN.includes(ch)) { depth++; continue; }
      if (CLOSE.includes(ch)) { if (depth > 0) depth--; continue; }
      if (depth !== 0) continue;
      if (ch === ':' && m[k + 1] === '=') { cut = k; break; }
      if (ch === 'w' && m.startsWith('where', k) && !/[\w'!?]/.test(m[k - 1] ?? ' ')
          && !/[\w'!?]/.test(m[k + 5] ?? ' ')) { cut = k; break; }
    }
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
      if (ns) { pushNamespace(nsStack, ns[1]); continue; }
      const sec = SEC_RE.exec(line);
      if (sec) { nsStack.push({ ns: false, name: sec[1] ?? null }); continue; }
      const en = END_RE.exec(line);
      if (en) { popEnd(nsStack, en[1]); continue; }
      if (END_BARE_RE.test(line)) { popEnd(nsStack, null); continue; }
      let m = DECL_RE.exec(line);
      if (!m) {
        // ★無名 instance(または `instance (priority := …) 名前` のように
        //   `DECL_RE` が括弧で止まる形)。名前が取れなければ `⟨無名⟩` を置く。
        const a = ANON_INST_RE.exec(line);
        if (!a) continue;
        const nm = INST_NAME_RE.exec(a[1]);
        m = [line, 'instance', nm ? nm[1] : ANON];
      }
      const full = qualify(nsStack, m[2]);
      // ★4 列目が statement(案 I)。名前で引いて外した失敗形を潰すための列である。
      //   ★無名 instance では**ここだけが手がかり**なので、型が入っていることが要である。
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
