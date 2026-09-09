// Gained 族が出口(AxSenTate / AxLemmaGraded / AxLemma)まで「推移的に」繋がっているかを測る。
// v1 の欠陥を 2 つ直した:
//   (a) docstring(/- ... -/)の中の `theorem ...` を宣言として拾っていた（偽陽性）
//   (b) 1 ホップしか見ていなかった（推移閉包を取る）
import fs from 'node:fs';
import path from 'node:path';

const DIR = 'lean/ABC3/Found/PGC';
const files = fs.readdirSync(DIR).filter(f => f.endsWith('.lean'));

// --- ブロックコメントを空白に潰す（行数は保つ） ---
function stripComments(text) {
  const out = text.split('');
  let depth = 0;
  for (let i = 0; i < text.length - 1; i++) {
    if (depth === 0 && text[i] === '/' && text[i + 1] === '-') { depth = 1; }
    else if (depth > 0 && text[i] === '/' && text[i + 1] === '-') { depth++; }
    if (depth > 0) {
      if (out[i] !== '\n') out[i] = ' ';
      if (text[i] === '-' && text[i + 1] === '/') {
        depth--;
        if (out[i + 1] !== '\n') out[i + 1] = ' ';
        i++;
      }
    }
  }
  if (depth > 0 && out[text.length - 1] !== '\n') out[text.length - 1] = ' ';
  return out.join('');
}

const clean = new Map();
for (const f of files) clean.set(f, stripComments(fs.readFileSync(path.join(DIR, f), 'utf8')));

// --- 宣言の切り出し（コメント除去済みの本文から） ---
const decls = [];
for (const f of files) {
  const lines = clean.get(f).split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const m = /^(?:noncomputable\s+)?(?:private\s+)?(theorem|lemma|def)\s+([A-Za-z_][A-Za-z0-9_'!?]*)/.exec(lines[i]);
    if (!m) continue;
    decls.push({ kind: m[1], name: m[2], file: f, line: i + 1, lines });
  }
}
// 本文（次の宣言まで）
for (let k = 0; k < decls.length; k++) {
  const d = decls[k];
  const nxt = decls[k + 1] && decls[k + 1].file === d.file ? decls[k + 1].line - 1 : d.lines.length;
  d.body = d.lines.slice(d.line - 1, nxt).join('\n');
  d.sig = d.body.split(/:=/)[0];
}
const byName = new Map();
for (const d of decls) if (!byName.has(d.name)) byName.set(d.name, d);

// --- 呼び出しグラフ: d.body に出てくる宣言名 ---
const allNames = [...byName.keys()].filter(n => n.length >= 6 && !n.endsWith('.src'));
const nameSet = new Set(allNames);
for (const d of decls) {
  const toks = d.body.match(/[A-Za-z_][A-Za-z0-9_'!?]*/g) || [];
  const s = new Set();
  for (const t of toks) if (nameSet.has(t) && t !== d.name) s.add(t);
  d.calls = s;
}

// --- 出口: 結論に AxSenTate / AxLemmaGraded が出る theorem ---
const exitPat = /\b(AxSenTate|AxLemmaGraded)\b/;
const exits = decls.filter(d => {
  if (d.kind === 'def') return false;
  const idx = d.sig.lastIndexOf(') :');
  const tail = idx >= 0 ? d.sig.slice(idx) : d.sig;
  return exitPat.test(tail);
});

// --- Gained 族 ---
const FAM = /^exists_norm_sub_algebraMap_le_prod_axDecay/;
const fam = decls.filter(d => d.kind === 'theorem' && FAM.test(d.name));
const famSet = new Set(fam.map(d => d.name));

// --- 推移閉包（出口から下向きに何が呼ばれるか） ---
function reach(startName) {
  const seen = new Set([startName]);
  const stack = [startName];
  while (stack.length) {
    const n = stack.pop();
    const d = byName.get(n);
    if (!d) continue;
    for (const c of d.calls) if (!seen.has(c)) { seen.add(c); stack.push(c); }
  }
  return seen;
}

console.log('宣言総数(コメント除去後) = ' + decls.length);
console.log('Gained 族 theorem = ' + fam.length + ' / ファイル = ' + new Set(fam.map(d => d.file)).size);
console.log('出口(結論 AxSenTate/AxLemmaGraded)の theorem = ' + exits.length);
console.log('');

console.log('=== 出口ごとに、推移的に Gained 族へ届くか ===');
let hit = 0;
for (const e of exits) {
  const R = reach(e.name);
  const used = [...R].filter(n => famSet.has(n));
  if (used.length) { hit++; console.log(`★届く  ${e.file}:${e.line} ${e.name}  ← ${used.join(', ')}`); }
}
console.log(`届く出口 = ${hit} / ${exits.length}`);
console.log('');

console.log('=== 逆向き: Gained 族から上向きに、出口へ届くか（誰が呼んでいるか） ===');
// 逆辺
const callers = new Map();
for (const d of decls) for (const c of d.calls) {
  if (!callers.has(c)) callers.set(c, new Set());
  callers.get(c).add(d.name);
}
function reachUp(startName) {
  const seen = new Set([startName]);
  const stack = [startName];
  while (stack.length) {
    const n = stack.pop();
    for (const c of (callers.get(n) || [])) if (!seen.has(c)) { seen.add(c); stack.push(c); }
  }
  return seen;
}
const exitNames = new Set(exits.map(d => d.name));
let up = 0;
for (const d of fam) {
  const U = reachUp(d.name);
  const es = [...U].filter(n => exitNames.has(n));
  if (es.length) { up++; console.log(`★${d.name} → ${es.join(', ')}`); }
}
console.log(`出口へ届く Gained 族 = ${up} / ${fam.length}`);
console.log('');

console.log('=== 参考: AxLemmaGraded を結論に持つ theorem ===');
for (const d of decls) {
  if (d.kind !== 'theorem') continue;
  const idx = d.sig.lastIndexOf(') :');
  const tail = idx >= 0 ? d.sig.slice(idx) : d.sig;
  if (/\bAxLemmaGraded\b/.test(tail)) {
    console.log(`  ${d.file}:${d.line}\t${d.name}\t${tail.replace(/\s+/g, ' ').slice(0, 80)}`);
  }
}

// --- 健全性検査: 既知の辺が本当に見つかるか（0 件が script のバグでないことの確認） ---
console.log('');
console.log('=== 健全性検査（既知の道が見えるか）===');
const probes = [
  ['axSenTate_of_deep_bound', 'axLemma_of_wildDescent_Icc'],
  ['axSenTate_of_cyclotomic_loss', 'axDecay_of_loss_bound'],
  ['axSenTate_of_axDecay', 'prod_le_rpow_of_geometric_shift'],
  ['axSenTate_of_deep_bound', 'axWildDescent_prime'],
  ['exists_norm_sub_algebraMap_le_prod_axDecay', 'norm_natCast_le_of_dvd'],
];
for (const [a, b] of probes) {
  if (!byName.has(a)) { console.log(`  ${a} → (宣言が見つからない)`); continue; }
  const R = reach(a);
  console.log(`  ${a} ⇝ ${b} : ${R.has(b) ? 'YES' : 'NO'}   (到達数 ${R.size})`);
}
