#!/usr/bin/env node
/**
 * 「mathlib に X は無い」を**もう一度測る**道具(2026-09-05、メタ第 2 回 / backlog M4)
 *
 * ## ★動機(実測)
 *
 * 2026-09-05 の 1 日で、「mathlib に無い」という記録が **4 件覆った**:
 *   `ULift.field` / `continuousCohomology` / `ProfiniteGrp.…completion` /
 *   `CompactSpace Gal(Ω/F)`(`FieldTheory/Galois/Profinite.lean:329` に無名 instance)。
 * いずれも「その数学を書き直しかけた」段階で見つかっている。**不在の誤りは高い**
 * ——既にある数学を数百行書き直させるからである。
 *
 * 不在の主張は `.absent` 49 件 + 自由文 355 件 = 404 件あり、
 * **再実行できる検索パターンを持つのは 30 件(7%)** だった。
 * `Meta/Claim.lean` は `absent (searched : String)` で探索範囲を要求しているのに、
 * 実際に書けているのはごく一部で、**93% は機械で再検査できない**。
 *
 * ## 規約(`check.mjs` の G11 が要求する形)
 *
 *     .absent "mathlib 全体を re:`WeilPairing|weil_pairing` で検索して 0 件(2026-08-20 実測)"
 *                                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 再実行できるパターン
 *
 * * `re:` の直後にバッククォートで囲んだ**正規表現**を書く。複数書いてよい。
 * * 測定時のヒット件数が 0 でないときは `re:`pat`→2` のように**そのときの件数**を書く
 *   (「2 件ヒットするがいずれも別物」という記録が実在する。件数を書いておけば
 *   **増えたときだけ**鳴らせる)。省略時は 0 とみなす。
 * * 測定日 `YYYY-MM-DD` も書く。日付の無い測定は古びたことが分からない。
 *
 * ## 何と照合するか
 *
 * `.cache/mathlib-index.txt`(`node tools/decl-index.mjs --mathlib`)の**全行**。
 * 索引は `種別 <TAB> 名前 <TAB> ファイル:行 <TAB> statement` なので、
 * **名前でも型でも**引ける。★2026-09-05 から**無名 instance も入る**(backlog M5)
 * ——これが入るまで、上の 4 件目は名前で grep しても出なかった。
 *
 * ## 使い方
 *
 *     node tools/absent-recheck.mjs           # 規約つき + 推定の両方を報告
 *     node tools/absent-recheck.mjs --strict  # 規約つきだけ(推定は出さない)
 *     node tools/absent-recheck.mjs --json
 *     node tools/absent-recheck.mjs --try 'CompactSpace Gal'   # ★書く前に件数を測る
 *     node tools/absent-recheck.mjs --selftest                 # 器具自身の較正
 *
 * 終了コード 1 は「**規約つきの主張が覆っている可能性がある**」ことだけを意味する
 * (推定側は参考であってゲートではない)。
 */

import { readFileSync, readdirSync, statSync, existsSync, mkdtempSync, mkdirSync, writeFileSync, rmSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { tmpdir } from 'node:os';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const argOf = (flag, dflt) => {
  const i = process.argv.indexOf(flag);
  return i >= 0 && process.argv[i + 1] ? process.argv[i + 1] : dflt;
};
const LEAN_SRC = argOf('--root', join(ROOT, 'lean', 'ABC3'));
const INDEX = argOf('--index', join(ROOT, '.cache', 'mathlib-index.txt'));
const JSON_OUT = process.argv.includes('--json');
const STRICT = process.argv.includes('--strict');
const CASE = process.argv.includes('--case');   // 既定は大文字小文字を無視する
const TRY = argOf('--try', null);
const SELFTEST = process.argv.includes('--selftest');
const SAMPLES = 3;

// ────────────────────────────────────────────────────────────────
// `.absent "…"` の実引数を拾う
// ────────────────────────────────────────────────────────────────

function walk(dir, acc = []) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else if (p.endsWith('.lean')) acc.push(p);
  }
  return acc;
}

/** コメントだけを空白へ潰す(**文字列は残す**——中身がこの道具の対象だから)。 */
export const stripCommentsKeepStrings = (src) => src
  .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
  .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length));

/** Lean の文字列リテラルを読む(改行を含みうる)。 */
function readString(src, i) {
  let out = '';
  for (let j = i + 1; j < src.length; j++) {
    const c = src[j];
    if (c === '\\') { out += src[j] + src[j + 1]; j++; continue; }
    if (c === '"') return { s: out, end: j + 1 };
    out += c;
  }
  return null;
}

/** 木の中の `.absent "…"` を全部返す。`{ file, line, arg }` */
export function collectAbsent(dir) {
  const rows = [];
  for (const f of walk(dir)) {
    const src = stripCommentsKeepStrings(readFileSync(f, 'utf8'));
    const re = /\.absent\s*(?=")/g;
    let m;
    while ((m = re.exec(src)) !== null) {
      const q = src.indexOf('"', m.index);
      const r = readString(src, q);
      if (!r) continue;
      rows.push({
        file: relative(dir, f).replace(/\\/g, '/'),
        line: src.slice(0, m.index).split('\n').length,
        arg: r.s,
      });
      re.lastIndex = r.end;
    }
  }
  return rows;
}

// ────────────────────────────────────────────────────────────────
// パターンの取り出し
// ────────────────────────────────────────────────────────────────

/** 規約 `re:`pat`→n` を全部取り出す。`n` 省略時は 0。 */
export function declaredPatterns(arg) {
  const out = [];
  const re = /re:\s*`([^`]+)`\s*(?:(?:→|->)\s*(\d+))?/g;
  let m;
  while ((m = re.exec(arg)) !== null) out.push({ pat: m[1], was: m[2] === undefined ? 0 : Number(m[2]) });
  return out;
}

/** 規約が無い記録から、それらしいパターンを**推定**する(参考であってゲートではない)。 */
const TOO_GENERIC = new Set([
  'mathlib', 'Mathlib', 'lean', 'Lean', 'sorry', 'true', 'false', 'Type', 'Prop',
  'grep', 'ls', 'papers.json', 'absent', 'unmeasured',
]);
export function guessPatterns(arg) {
  const out = new Set();
  // ★推定したパターンは**語境界で囲む**。囲まないと `valence` が `Equivalence` に
  //   当たって 1,435 件になり(実測)、報告が読めなくなる。
  //   規約つきのパターンは著者が書いた正規表現そのままを使う——囲まない。
  const bounded = (t) => `\\b(?:${t})\\b`;
  // (1) バッククォートで囲まれた語(既存の良い記録が実際に使っている形)
  for (const m of arg.matchAll(/`([^`]+)`/g)) {
    const t = m[1].trim();
    if (t.length < 3 || t.includes('/') || /\.lean$/.test(t) || TOO_GENERIC.has(t)) continue;
    if (!/^[A-Za-z_][\w'.|*\\+?^$()\[\]\- ]*$/.test(t)) continue;
    out.add(bounded(t));
  }
  // (2) 裸の選言(`A|B|C`)——「A|B|C で grep、0 件」と書いてある記録が実在する
  for (const m of arg.matchAll(/[A-Za-z_][\w'-]*(?:\|[A-Za-z_][\w'-]*)+/g)) {
    if (TOO_GENERIC.has(m[0])) continue;
    out.add(bounded(m[0]));
  }
  return [...out];
}

// ────────────────────────────────────────────────────────────────
// 索引に当てる
// ────────────────────────────────────────────────────────────────

function runPatterns(lines, pats) {
  const res = new Map();
  for (const p of pats) {
    let re;
    try { re = new RegExp(p, CASE ? '' : 'i'); } catch { res.set(p, { bad: true }); continue; }
    let n = 0; const hits = [];
    for (const l of lines) {
      if (!re.test(l)) continue;
      n++;
      if (hits.length < SAMPLES) hits.push(l);
    }
    res.set(p, { n, hits });
  }
  return res;
}

// ────────────────────────────────────────────────────────────────

const t0 = Date.now();
if (!existsSync(INDEX)) {
  console.error(`${relative(ROOT, INDEX)} が無い。\`node tools/decl-index.mjs --mathlib\` で作ること`);
  process.exit(2);
}
const lines = readFileSync(INDEX, 'utf8').split('\n').filter(Boolean);

// ★`--try`: `.absent` を**書く前に**件数を測るための口。
//   `re:`pat`→N` の `N` はここで得た数字を書く。
if (TRY !== null) {
  const v = runPatterns(lines, [TRY]).get(TRY);
  if (v.bad) { console.error(`正規表現として読めない: ${TRY}`); process.exit(2); }
  console.log(`\`${TRY}\` → ${v.n} 件(${relative(ROOT, INDEX)}、${CASE ? '大小区別' : '大小無視'})`);
  for (const h of v.hits) console.log(`   ${h.replace(/\t/g, '  ').slice(0, 160)}`);
  process.exit(0);
}

// ★`--selftest`: 器具が「覆っている」を実際に鳴らせるかを較正する。
//   ゲート(`check.mjs`)と同じ規律——**壊れた入力を落とせることを確かめてから使う**。
if (SELFTEST) {
  const tmp = mkdtempSync(join(tmpdir(), 'abc3-absent-'));
  mkdirSync(join(tmp, 'Found'), { recursive: true });
  // 索引に必ず在るもの / 必ず無いものを使う(前者は Lean の基本型、後者は無意味語)。
  writeFileSync(join(tmp, 'Found', 'a.lean'),
    'def x.status := LeanStatus.absent "mathlib 全体を re:`CommMonoidWithZero`→0 で検索して 0 件(2026-01-01)"\n' +
    'def y.status := LeanStatus.absent "mathlib 全体を re:`ZzzNoSuchDeclZzz`→0 で検索して 0 件(2026-01-01)"\n' +
    'def z.status := LeanStatus.absent "実測して無かった(2026-01-01)"\n', 'utf8');
  const rs = collectAbsent(tmp);
  for (const r of rs) { r.declared = declaredPatterns(r.arg); r.guessed = r.declared.length ? [] : guessPatterns(r.arg); }
  const rr = runPatterns(lines, rs.flatMap((r) => r.declared.map((d) => d.pat)));
  const over = rs.filter((r) => r.declared.some((d) => (rr.get(d.pat)?.n ?? 0) > d.was));
  rmSync(tmp, { recursive: true, force: true });
  const checks = [
    ['E1 在るものを「無い」と書いた記録を鳴らせる', over.some((r) => /CommMonoidWithZero/.test(r.arg))],
    ['E2 本当に無いものは鳴らさない', !over.some((r) => /ZzzNoSuchDecl/.test(r.arg))],
    ['E3 パターンの無い記録を「引けない」と数える', rs.some((r) => !r.declared.length && !r.guessed.length)],
    ['E4 3 件すべて拾えている', rs.length === 3],
  ];
  for (const [label, ok] of checks) console.log(`  ${ok ? 'ok ' : 'NG '} ${label}`);
  const bad = checks.filter(([, ok]) => !ok).length;
  console.log(`\n  absent-recheck selftest: ${checks.length - bad}/${checks.length} PASS`);
  process.exit(bad ? 1 : 0);
}

const rows = collectAbsent(LEAN_SRC);

// 走らせるパターンを一度に集めてから 1 パスで当てる(索引は 46MB ある)。
const all = new Set();
for (const r of rows) {
  r.declared = declaredPatterns(r.arg);
  r.guessed = r.declared.length ? [] : guessPatterns(r.arg);
  for (const d of r.declared) all.add(d.pat);
  for (const g of r.guessed) all.add(g);
}
const res = runPatterns(lines, [...all]);

const overturned = [];   // 規約つきで件数が増えたもの
const suspicious = [];   // 推定でヒットが出たもの
const unmeasurable = []; // パターンを取り出せなかったもの
for (const r of rows) {
  if (r.declared.length) {
    for (const d of r.declared) {
      const v = res.get(d.pat);
      if (v && !v.bad && v.n > d.was) overturned.push({ r, d, v });
    }
  } else if (r.guessed.length) {
    for (const g of r.guessed) {
      const v = res.get(g);
      if (v && !v.bad && v.n > 0) suspicious.push({ r, g, v });
    }
  } else {
    unmeasurable.push(r);
  }
}

const nDeclared = rows.filter((r) => r.declared.length).length;
const nGuessed = rows.filter((r) => !r.declared.length && r.guessed.length).length;

if (JSON_OUT) {
  console.log(JSON.stringify({
    index: relative(ROOT, INDEX), indexLines: lines.length,
    total: rows.length, declared: nDeclared, guessed: nGuessed, unmeasurable: unmeasurable.length,
    overturned: overturned.map(({ r, d, v }) => ({ file: r.file, line: r.line, pat: d.pat, was: d.was, now: v.n, hits: v.hits })),
    suspicious: suspicious.map(({ r, g, v }) => ({ file: r.file, line: r.line, pat: g, now: v.n, hits: v.hits })),
    unmeasurableRows: unmeasurable.map((r) => ({ file: r.file, line: r.line, arg: r.arg.replace(/\s+/g, ' ').slice(0, 120) })),
    ms: Date.now() - t0,
  }, null, 1));
} else {
  console.log(`★不在の主張の再検査 —— 索引 ${relative(ROOT, INDEX)}(${lines.length} 宣言、` +
              `${CASE ? '大小区別' : '大小無視'})`);
  console.log(`  .absent ${rows.length} 件 / 規約つき ${nDeclared} / 推定で引ける ${nGuessed} / ` +
              `引けない ${unmeasurable.length}\n`);

  console.log(`— 規約つき(re:\`…\`)で覆っている可能性: ${overturned.length} 件`);
  for (const { r, d, v } of overturned) {
    console.log(`  ★ ${r.file}:${r.line}  \`${d.pat}\`  記録 ${d.was} → いま ${v.n} 件`);
    for (const h of v.hits) console.log(`       ${h.replace(/\t/g, '  ').slice(0, 150)}`);
  }
  if (!overturned.length) console.log('  (なし)');

  if (!STRICT) {
    console.log(`\n— 推定パターンでヒットが出たもの(★参考。当たっているとは限らない): ${suspicious.length} 件`);
    for (const { r, g, v } of suspicious) {
      console.log(`  ? ${r.file}:${r.line}  \`${g}\` → ${v.n} 件`);
      for (const h of v.hits) console.log(`       ${h.replace(/\t/g, '  ').slice(0, 150)}`);
    }
    if (!suspicious.length) console.log('  (なし)');

    console.log(`\n— パターンを取り出せない(= 機械で再検査できない)記録: ${unmeasurable.length} 件`);
    for (const r of unmeasurable) {
      console.log(`  - ${r.file}:${r.line}  ${r.arg.replace(/\s+/g, ' ').slice(0, 110)}`);
    }
  }
  console.log(`\n${Math.round((Date.now() - t0) / 100) / 10} 秒`);
}
process.exit(overturned.length ? 1 : 0);
