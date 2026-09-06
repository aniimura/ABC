#!/usr/bin/env node
/**
 * 配線されていない既存の数学を探す —— `sorry` の隣に**もう証明がある**かを機械で当てる。
 *
 * ★動機(実測、2026-09-06 メタ第 7 回)
 *   2026-09-05〜06 の 2 日で、**未解決の数学ではなく配線漏れ**だった節点が 4 件出た:
 *     Divisor/weil ・ Divisor/arith ・ NumberField/cheb ・ GaloisRep/WeilPairingDef。
 *   4 件とも「たまたま在庫調査をした agent が気づいた」だけで、**機械は 1 件も見つけていない**。
 *   `tools/frontier.mjs` は「`sorry` があるノード」を出すが、
 *   **その `sorry` が既に `Found` で解けているか**を見ていない。
 *
 * ★この道具が答える 2 つの問い(既存の道具の挙動は変えない。新規 1 本)
 *   (1) `--wired`(既定): `sorry` を持つ宣言の**結論**に近い、`sorry` 無しの宣言を 3 件出す。
 *       ★**完全一致は狙わない。** 4 件目も `[CharZero F]`・`[IsAlgClosed F]` を足す差があった。
 *       人が 5 分で確かめられる候補が出れば勝ちである。
 *   (2) `--dead`: **空撃ち** —— import はされているのに、
 *       import 元が**その語を 1 つも使っていない**モジュール(3 件目の型)。
 *
 * 使い方:
 *   node tools/unwired.mjs                    # sorry ごとに候補 3 件(既定)
 *   node tools/unwired.mjs --min 60           # 一致率の下限で絞る(既定 0)
 *   node tools/unwired.mjs --node Skeleton/Divisor/SchemeWeil.lean
 *   node tools/unwired.mjs --dead             # 空撃ち
 *   node tools/unwired.mjs --json
 *   node tools/unwired.mjs --marks            # ★候補と空撃ちを 1 回の走査で返す(frontier.mjs 用)
 *   node tools/unwired.mjs --selftest         # 較正(既知 6 組が上位 3 件に入るか)
 *
 * ★読み方(実測に基づく)
 *   - 一致率 80% 以上 かつ 情報量 15 以上 … **ほぼ当たり**。まず開いて見る。
 *   - 結論が `ℝ` / `WeilDiv X` のように短いものは 100% が出るが**意味がない**。
 *     情報量(idf の絶対量)を併記してあるので、そこで捨てる。
 *   - `Check/` は既定で在庫から外す —— 退化例・反証は**わざと同じ形**をしており、
 *     答えではないのに上位に来る(実測で雑音の主因だった)。`--check` で入る。
 */

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const SRC = join(ROOT, 'lean', 'ABC3');
const args = process.argv.slice(2);
const flag = (f) => args.includes(f);
const opt = (f) => { const i = args.indexOf(f); return i >= 0 ? args[i + 1] : null; };

function walk(dir, acc = []) {
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else if (p.endsWith('.lean')) acc.push(p);
  }
  return acc;
}

/* ------------------------------------------------------------------ 字句 */

/** 名前の文字類。`tools/decl-index.mjs` と同じもの(2026-09-06 の M14 修正込み)。
 *  ★ASCII だけの `\w` にすると `preΨ₄` が `pre` で切れる。 */
const IDX = "\\w'!?\\u00C0-\\u024F\\u0370-\\u03FF\\u1D00-\\u1DBF\\u2070-\\u209F\\u2100-\\u214F\\uD800-\\uDFFF";
const ID1 = "A-Za-z_\\u00C0-\\u024F\\u0370-\\u03FF\\u2100-\\u214F\\uD800-\\uDFFF";
const NAME = `[${ID1}][${IDX}]*(?:\\.[${IDX}]+)*`;
const MODS = '(?:private\\s+|protected\\s+|noncomputable\\s+|scoped\\s+|local\\s+)*';
const DECL_RE = new RegExp(
  `^(?:@\\[[^\\]]*\\]\\s*)?${MODS}(theorem|lemma|def|abbrev|instance|structure|inductive|class)\\s+(${NAME})`);
const TOKEN_RE = new RegExp(`[${ID1}][${IDX}]*(?:\\.[${IDX}]+)*|[0-9]+`, 'g');
const BARE_RE = new RegExp(`[${ID1}][${IDX}]*`, 'g');
const NS_RE = /^\s*namespace\s+([A-Za-z_][\w'.₀-₉]*)/;
const END_RE = /^\s*end\s+([A-Za-z_][\w'.₀-₉]*)\s*$/;
/** `.src` / `.needs` などの記録用 def。数学ではないので在庫にも問いにも入れない。 */
const META_DECL = /\.(src|needs|deviation|deviations|note|notes)$/;

/** ブロックコメント・行コメントを空白に潰す。 */
function decomment(text) {
  let out = '', i = 0, depth = 0;
  while (i < text.length) {
    if (text.startsWith('/-', i)) { depth++; i += 2; out += '  '; continue; }
    if (depth && text.startsWith('-/', i)) { depth--; i += 2; out += '  '; continue; }
    if (depth) { out += text[i] === '\n' ? '\n' : ' '; i++; continue; }
    if (text.startsWith('--', i)) { while (i < text.length && text[i] !== '\n') { out += ' '; i++; } continue; }
    out += text[i]; i++;
  }
  return out;
}

/** 文字列リテラルを潰す。
 *  ★これが無いと `.needs` の**本文に「sorry」と書いてあるだけ**の 33 件を
 *  sorry として数えてしまう(実測: 90 件 → 56 件)。 */
function stripStrings(text) {
  let out = '', i = 0, inStr = false;
  while (i < text.length) {
    const c = text[i];
    if (!inStr && c === '"') { inStr = true; out += ' '; i++; continue; }
    if (inStr) {
      if (c === '\\') { out += '  '; i += 2; continue; }
      if (c === '"') { inStr = false; out += ' '; i++; continue; }
      out += c === '\n' ? '\n' : ' '; i++; continue;
    }
    out += c; i++;
  }
  return out;
}

const OPEN = '([{⟨⦃', CLOSE = ')]}⟩⦄';

/** 署名(`:=` / `where` の手前まで)を 1 行に畳む。
 *  ★**深さ 0 の `:=` でだけ切る** —— `f (p := p)` の名前つき引数で切ると
 *  結論が `(cyclotomicCharacterObject (p` のように途中で消える(実測 5 件)。 */
function signatureOf(lines, start, end) {
  const buf = [];
  let d = 0;
  outer: for (let j = start; j < end; j++) {
    const s = lines[j];
    for (let i = 0; i < s.length; i++) {
      const c = s[i];
      if (OPEN.includes(c)) d++;
      else if (CLOSE.includes(c)) d--;
      else if (d === 0 && (s.startsWith(':=', i) || /^\bwhere\b/.test(s.slice(i)))) {
        buf.push(s.slice(0, i)); break outer;
      }
    }
    buf.push(s);
    if (buf.join(' ').length > 4000) break;
  }
  return buf.join(' ').replace(/\s+/g, ' ').trim();
}

/** 署名から**結論**を切り出す(深さ 0 の `:` の右)。 */
function conclusionOf(sig) {
  let d = 0;
  for (let i = 0; i < sig.length; i++) {
    const c = sig[i];
    if (OPEN.includes(c)) d++;
    else if (CLOSE.includes(c)) d--;
    else if (c === ':' && d === 0 && sig[i + 1] !== '=' && sig[i + 1] !== ':' && sig[i - 1] !== ':') {
      return sig.slice(i + 1).trim();
    }
  }
  return '';
}

const SYMS = ['≠', '≃', '≅', '↔', '⁻¹', '∀', '∃', '∧', '∨', '¬', '∈', '∉', '⊆', '⊂', '∣', '≤', '≥',
  '<', '>', '=', '→', '•', '∑', '∏', '∘', '⊗', '×', '∪', '∩', '^', '+', '*', '/', '%'];

/** 結論を鍵の集合にする。識別子・その末尾成分・部分語(camel/snake)・記号。
 *  ★部分語を入れているのは `weilPairing` と `weilPairingVal` を繋ぐためである
 *  (較正 1 番目がこれで当たる)。 */
function keysOf(text) {
  const out = [];
  for (const m of text.matchAll(TOKEN_RE)) {
    const t = m[0];
    out.push('id:' + t);
    const last = t.split('.').pop();
    if (last !== t) out.push('id:' + last);
    for (const w of last.split(/[_']/).flatMap((s) => s.split(/(?=[A-Z])/)))
      if (w.length >= 3) out.push('w:' + w.toLowerCase());
  }
  for (const s of SYMS) if (text.includes(s)) out.push('s:' + s);
  return out;
}

/* ------------------------------------------------------------ 木を舐める */

function scan() {
  const decls = [];
  const mods = new Map();
  for (const file of walk(SRC)) {
    const rel = relative(SRC, file).replace(/\\/g, '/');
    const mod = 'ABC3.' + rel.replace(/\.lean$/, '').split('/').join('.');
    const raw = readFileSync(file, 'utf8');
    const lines = decomment(raw).split(/\r?\n/);
    const starts = [];
    for (let i = 0; i < lines.length; i++) { const s = lines[i]; if (s.trim() && !/^\s/.test(s)) starts.push(i); }
    const nsAt = []; { const st = [];
      for (let i = 0; i < lines.length; i++) {
        const ns = NS_RE.exec(lines[i]); if (ns) { st.push(ns[1]); nsAt[i] = st.join('.'); continue; }
        const en = END_RE.exec(lines[i]);
        if (en && st.length && st[st.length - 1] === en[1]) { st.pop(); nsAt[i] = st.join('.'); continue; }
        nsAt[i] = st.join('.');
      } }
    const names = new Set(), imports = [];
    for (const l of raw.split(/\r?\n/)) { const im = /^import\s+(\S+)/.exec(l); if (im) imports.push(im[1]); }
    for (let k = 0; k < starts.length; k++) {
      const i = starts[k];
      const m = DECL_RE.exec(lines[i]); if (!m) continue;
      const end = k + 1 < starts.length ? starts[k + 1] : lines.length;
      const block = stripStrings(lines.slice(i, end).join('\n')).replace(/\s+/g, ' ').trim();
      const hasSorry = /\bsorry\b/.test(block);
      /** ★**裸の sorry** —— 証明本体が `:= sorry` / `:= by sorry` だけのもの。
       *  途中まで書いてある sorry は「残っている数学」、裸の sorry は
       *  「まだ誰も手を付けていない」= 在庫が既にある可能性が高い。実測でこれが効く
       *  (`semistableAt_veluQuotientFull` は本文が既に `Found` を 3 本呼んでいる)。 */
      const bare = hasSorry && /:=\s*(by\s+)?sorry$/.test(block);
      const sig = signatureOf(lines, i, end);
      const short = m[2].split('.').pop();
      if (!META_DECL.test(m[2])) names.add(short);
      decls.push({
        kind: m[1], name: (nsAt[i] ? nsAt[i] + '.' : '') + m[2], short, rel, mod,
        line: i + 1, sig, concl: conclusionOf(sig), hasSorry, bare, meta: META_DECL.test(m[2]),
      });
    }
    const tokens = new Set(); for (const t of raw.matchAll(BARE_RE)) tokens.add(t[0]);
    mods.set(mod, { rel, names, imports, tokens, decls: 0 });
  }
  for (const d of decls) if (!d.meta) mods.get(d.mod).decls++;
  return { decls, mods };
}

/* -------------------------------------------------------- (1) 候補を出す */

function buildIndex(decls, withCheck) {
  const pool = decls.filter((d) => !d.hasSorry && !d.meta && (d.concl || d.sig).length > 4
    && (withCheck || !d.rel.startsWith('Check/')));
  const df = new Map();
  const keys = pool.map((d) => {
    const s = new Set(keysOf(d.concl || d.sig));
    for (const k of s) df.set(k, (df.get(k) ?? 0) + 1);
    return s;
  });
  const N = pool.length;
  const idf = (k) => Math.log((N + 1) / ((df.get(k) ?? 0) + 1));
  const inv = new Map();
  keys.forEach((s, i) => { for (const k of s) { if (!inv.has(k)) inv.set(k, []); inv.get(k).push(i); } });
  return { pool, keys, idf, inv };
}

/** 問い(結論の文字列)に対して候補を返す。
 *  cov  … 問いの情報量のうち候補が覆う割合
 *  mass … 覆った情報量の絶対値。★短い結論(`ℝ` など)を捨てるのはここ。
 *  imported … 問いのファイルが直接 import しているモジュールの集合(印をつけるだけ)。 */
function rank(ix, concl, excludeRel, topK, imported) {
  const qk = [...new Set(keysOf(concl))];
  const qmass = qk.reduce((a, k) => a + ix.idf(k), 0) || 1;
  const acc = new Map();
  for (const k of qk) {
    const w = ix.idf(k);
    if (w < 0.5) continue;                       // ありふれた語は引かない
    for (const i of ix.inv.get(k) ?? []) acc.set(i, (acc.get(i) ?? 0) + w);
  }
  const rows = [...acc].map(([i, s]) => ({ i, cov: s / qmass, mass: s, size: ix.keys[i].size }));
  rows.sort((a, b) => b.cov - a.cov || b.mass - a.mass || a.size - b.size);
  const out = [];
  for (const r of rows) {
    const d = ix.pool[r.i];
    if (d.rel === excludeRel) continue;          // 同じファイルの兄弟は答えではない
    out.push({ name: d.name, short: d.short, rel: d.rel, line: d.line, kind: d.kind,
      cov: r.cov, mass: r.mass, qmass, imported: imported ? imported.has(d.mod) : false });
    if (out.length >= topK) break;
  }
  return out;
}

/* ---------------------------------------------------------- (2) 空撃ち */

/** 宣言を 1 つも持たないモジュール = 束ねるだけの `Skeleton.lean` 等。
 *  ★これを消費者に数えると全部が「空撃ち」になる(実測 2,594 辺のうち大半)。 */
function deadModules(mods, decls) {
  const sorryMods = new Set(decls.filter((d) => d.hasSorry && !d.meta).map((d) => d.mod));
  const inbound = new Map();
  for (const [mod, m] of mods) for (const im of m.imports) {
    if (!mods.has(im)) continue;
    if (!inbound.has(im)) inbound.set(im, []);
    inbound.get(im).push(mod);
  }
  const rows = [];
  for (const [mod, m] of mods) {
    if (!m.names.size) continue;
    const ins = (inbound.get(mod) ?? []).filter((n) => mods.get(n).decls > 0);
    const aggr = (inbound.get(mod) ?? []).length - ins.length;
    const alive = ins.filter((n) => [...m.names].some((x) => mods.get(n).tokens.has(x)));
    if (alive.length) continue;
    rows.push({ rel: m.rel, decls: m.names.size, consumers: ins.length, aggregators: aggr,
      sorry: sorryMods.has(mod),
      kind: ins.length === 0 ? '消費者なし' : '空撃ち', blanks: ins.map((n) => mods.get(n).rel) });
  }
  return rows;
}

/* ------------------------------------------------------------- 較正 */

/** ★既知の 6 組。1 番目だけは**修正前の署名**(第 1036 で配線して消えた)を持っている。 */
const FIXTURES = [
  { label: 'weilPairing_nondeg → exists_pairing_ne_one(第 1036、4 件目)',
    rel: 'Skeleton/GaloisRep/WeilPairingDef.lean',
    concl: '∃ Q : W.Point, n • Q = 0 ∧ weilPairing W n P Q ≠ 1',
    want: 'exists_pairing_ne_one' },
  { label: 'inertia_recoverable → inertia_recoverable_of_prop12(backlog M3 が手で見つけた辺)',
    rel: 'Skeleton/PGC/Section1Cor13.lean', find: 'inertia_recoverable',
    want: 'inertia_recoverable_of_prop12' },
  { label: 'isDiscreteValuationRing_stalk_of_codimOne(5 件目)',
    rel: 'Skeleton/Divisor/SchemeWeil.lean', find: 'isDiscreteValuationRing_stalk_of_codimOne',
    want: 'isDiscreteValuationRing_stalk_of_codimOne' },
  { label: 'ordAtDiv_mul → ordPt_mul(5 件目と同じファイル)',
    rel: 'Skeleton/Divisor/SchemeWeil.lean', find: 'ordAtDiv_mul', want: 'ordPt_mul' },
  { label: 'finite_support_ordAtDiv → finite_ordPt_ne_zero',
    rel: 'Skeleton/Divisor/SchemeWeil.lean', find: 'finite_support_ordAtDiv', want: 'finite_ordPt_ne_zero' },
  { label: 'alpha_in_modl_image → alpha_mem_map_of_galTate',
    rel: 'Skeleton/GenEll/GaloisLocal.lean', find: 'alpha_in_modl_image', want: 'alpha_mem_map_of_galTate' },
];

/* ------------------------------------------------------------- 本体 */

const t0 = Date.now();
const { decls, mods } = scan();

if (flag('--dead')) {
  const rows = deadModules(mods, decls);
  const only = opt('--owner');
  const sel = rows.filter((r) => !only || r.rel.startsWith(only));
  if (flag('--json')) { console.log(JSON.stringify({ generated: new Date().toISOString(), dead: sel }, null, 1)); process.exit(0); }
  const hit = sel.filter((r) => r.kind === '空撃ち');
  const shown = [...hit, ...sel.filter((x) => x.kind !== '空撃ち')]
    .filter((r) => flag('--all') || r.rel.startsWith('Skeleton/'))
    .sort((a, b) => (b.sorry - a.sorry) || (b.consumers - a.consumers) || (b.decls - a.decls) || a.rel.localeCompare(b.rel));
  console.log('★空撃ち —— import されているのに 1 語も使われていないモジュール');
  console.log(`  宣言を持つモジュール ${[...mods.values()].filter((m) => m.names.size).length} 本`);
  console.log(`  うち **実消費者がいるのに 1 語も使われていない**: ${hit.length} 本`);
  console.log(`  うち 実消費者が 0(束ねる ${'`Skeleton.lean`'} だけが import): ${sel.length - hit.length} 本`);
  console.log('  ★束ねるだけのモジュール(宣言 0)は消費者に数えていない。');
  console.log(`  ★★このうち **sorry を持つもの ${shown.filter((r) => r.sorry).length} 本** —— ` +
    '「誰も使わない sorry」である。配線するか畳むかの判断が要る。');
  console.log();
  for (const r of shown) {
    console.log(`  ${r.sorry ? '★' : ' '}[${r.kind}] ${r.rel}  (宣言 ${r.decls})${r.sorry ? '  ← sorry あり' : ''}`);
    for (const b of r.blanks.slice(0, 3)) console.log(`        ← ${b} が import しているが 1 語も使っていない`);
  }
  if (!flag('--all')) console.log('\n  （Skeleton 以外を隠している。--all で全部）');
  console.log(`\n所要 ${((Date.now() - t0) / 1000).toFixed(2)} 秒`);
  process.exit(0);
}

const ix = buildIndex(decls, flag('--check'));

if (flag('--selftest')) {
  let ok = 0;
  for (const f of FIXTURES) {
    let concl = f.concl;
    if (!concl) {
      const q = decls.find((d) => d.rel === f.rel && d.short === f.find);
      concl = q ? q.concl : '';
    }
    const r = concl ? rank(ix, concl, f.rel, 3) : [];
    const at = r.findIndex((x) => x.short === f.want);
    const pass = at >= 0;
    if (pass) ok++;
    console.log(`${pass ? 'OK  ' : 'NG  '} ${f.label}`);
    console.log(`      ${pass ? `${at + 1} 位 / 3` : '上位 3 件に無い'}` +
      (r.length ? `  ── ${r.map((x) => `${x.short}(${(x.cov * 100) | 0}%)`).join(' , ')}` : '  ── 候補なし'));
  }
  console.log(`\n較正 ${ok}/${FIXTURES.length}`);
  process.exit(ok === FIXTURES.length ? 0 : 1);
}

const minCov = Number(opt('--min') ?? 0) / 100;
const node = opt('--node');
const topK = Number(opt('--top') ?? 3);
const queries = decls.filter((d) => d.hasSorry && !d.meta && !d.rel.startsWith('Meta/')
  && (!node || d.rel === node || d.rel.startsWith(node)));

const report = [];
for (const q of queries.sort((a, b) => a.rel.localeCompare(b.rel) || a.line - b.line)) {
  const imported = new Set(mods.get(q.mod)?.imports ?? []);
  const cands = q.concl ? rank(ix, q.concl, q.rel, topK, imported).filter((c) => c.cov >= minCov) : [];
  report.push({ rel: q.rel, line: q.line, kind: q.kind, name: q.short, concl: q.concl, bare: q.bare, cands });
}

/** ★`--marks`(2026-09-06、メタ第 9 回。backlog M18): `--json` の中身に `dead` も入れて返す。
 *  `tools/frontier.mjs` が**印を付けるために 1 回だけ呼ぶ**口である。
 *  ★2 回呼ぶ(`--json` と `--dead --json`)と木を 2 度舐めることになる:
 *  実測 2.12 + 1.87 = 3.99 秒 → 1 回にまとめて **2.35 秒**。
 *  `--json` 単体の出力は 1 バイトも変えない(`dead` は `--marks` のときだけ足す)。 */
if (flag('--json') || flag('--marks')) {
  const out = { generated: new Date().toISOString(), pool: ix.pool.length, report };
  if (flag('--marks')) out.dead = deadModules(mods, decls);
  console.log(JSON.stringify(out, null, 1));
  process.exit(0);
}

const isStrong = (r) => r.bare && r.cands.some((c) => c.cov >= 0.8 && c.mass >= 15);
const strong = report.filter(isStrong);
console.log('★配線されていない既存の数学の候補');
console.log(`  在庫(sorry 無しの宣言) ${ix.pool.length} 件${flag('--check') ? '（Check/ 込み）' : '（Check/ を除く）'}`);
console.log(`  問い(sorry を持つ宣言) ${queries.length} 件 / うち **裸の sorry**(本体が \`by sorry\` だけ) ${queries.filter((q) => q.bare).length} 件`);
console.log(`  ★★裸 かつ 一致率 80% 以上 かつ 情報量 15 以上: **${strong.length} 件** ← まずここを見る`);
const thin = report.filter((r) => r.cands.length && r.cands[0].qmass < 15);
if (thin.length) console.log(`  ・結論の情報量が 15 未満で**引けない**問い ${thin.length} 件: ` +
  `${thin.map((r) => r.name).join(' , ')}（結論が \`ℝ\` などの \`def\`。名前で引くしかない）`);
console.log();
for (const r of report) {
  const mark = isStrong(r) ? '★' : ' ';
  console.log(`${mark} ${r.rel}:${r.line}  ${r.kind} ${r.name}  ${r.bare ? '[裸]' : '[途中まで書いてある]'}`);
  console.log(`    結論: ${r.concl.slice(0, 160)}${r.concl.length > 160 ? ' …' : ''}`);
  if (!r.cands.length) { console.log('    （候補なし）'); console.log(); continue; }
  for (const c of r.cands) {
    console.log(`    ${String((c.cov * 100) | 0).padStart(3)}%  情報量 ${c.mass.toFixed(0).padStart(3)}  ` +
      `${c.name}  (${c.rel}:${c.line})${c.imported ? '  [import 済]' : ''}`);
  }
  console.log();
}
console.log(`所要 ${((Date.now() - t0) / 1000).toFixed(2)} 秒`);
console.log('★読み方: 一致率は**結論の情報量のうち候補が覆った割合**。');
console.log('  情報量が小さい(結論が `ℝ` など)ときは 100% でも意味が無い。');
console.log('  当たっていたら、その `Found` を import して配線する（4 件目は 5 行だった）。');
