// ゴール「[GenEll] を IUT が要る範囲で Found/ に完成させる」の進捗を機械で数える。
//
// ★この道具が答える問い: **需要のある項目のうち何件が実装済みか。**
// ★長期ゴールの本体は `ResearchPaper/genell-goal.md`。
//
// ★`tools/frdi-progress.mjs` と同じ規則で作った別ファイルである。共通化していない理由:
//   (1) FrdI 側は別セッションが稼働中で、共有ファイルを触ると衝突する。
//   (2) ★本道具は frdi-progress に**無い修正**を1つ入れている(下記)。
//       先に共通化すると、その修正が FrdI 側の数を予告なく動かす。
//   共通化はどちらも止まっているときに行うこと。
//
// ★★frdi-progress(および tools/intra-graph.mjs)との違い —— 節見出しでも本体を切る
//   両者は項目の本体を「次の**宣言行**まで」と近似する。節の最後の項目では
//   **次の節の導入文を巻き込む**。2026-08-16 に GenEll で実害を実測した:
//     §2 は `Theorem 2.1` 1 項目しか持たないため、その本体が `Section 3:` の
//     導入文(`… is "rather large" [cf. Theorem 3.8]`)まで伸び、
//     ★**Theorem 2.1 が Theorem 3.8 に依存する**という誤った辺が立っていた。
//     到達件数は 17 件と出るが、節見出しで切れば **9 件**である。
//   → 本道具は本体の終端を「次の宣言行 **または** 次の `Section N:` 見出し」とする。
//
// ★分母の限界(frdi-progress と同じ):
//   pdftotext 経由で目視していない(事実2)。番号の無い依存を数え落とし、
//   「cf.」型の案内を数え過ぎる。**桁を見るための数**である。
//
// ★分子の規則(frdi-progress と同じ。2値):
//   `.src` の `item` が **`Kind N.M` ちょうど**のものだけを実装済みとする。
//   条つき(`"Lemma 3.1, (i)"`)は **数えないが、見えるように印字する**。
//
// 使い方: node tools/genell-progress.mjs [--json]

import { readFileSync, writeFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const FILE_OF = JSON.parse(readFileSync(join(ROOT, 'tools', '_fileof.json'), 'utf8'));
const TAGRE = new RegExp(`\\[([A-Za-z][A-Za-z0-9]*)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');

const SECRE = /^Section \d+:/;

function load(tag) {
  const f = FILE_OF[tag];
  const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) return null;
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) { const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  const stops = [];
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(declRe);
    if (m) { list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] }); stops.push(i); }
    else if (SECRE.test(lines[i])) stops.push(i);      // ★ここが frdi-progress との違い
  }
  for (const d of list) d.end = stops.find((s) => s > d.line) ?? lines.length;
  const decls = new Map();
  for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  return { decls, lines };
}

const G = load('GenEll');
if (!G) { console.error('★[GenEll] の .txt が無い。0_Source を確認すること。'); process.exit(1); }

// ── 分母: 需要の推移閉包
const CONSUMERS = ['IUTchI', 'IUTchII', 'IUTchIII', 'IUTchIV', 'EtTh'];
const direct = new Map();
const missingConsumers = [];
for (const tag of CONSUMERS) {
  const P = load(tag);
  if (!P) { missingConsumers.push(tag); continue; }
  for (const m of P.lines.join('\n').matchAll(TAGRE)) {
    if (m[1] !== 'GenEll') continue;
    const k = `${m[2]} ${m[3]}`;
    if (!G.decls.has(k)) continue;
    if (!direct.has(k)) direct.set(k, new Set());
    direct.get(k).add(tag);
  }
}

const bodyOf = (name) => { const d = G.decls.get(name); return G.lines.slice(d.line, d.end).join('\n'); };
function intraEdges(name) {
  const body = bodyOf(name);
  const spans = []; const out = new Set();
  for (const m of body.matchAll(TAGRE)) spans.push([m.index, m.index + m[0].length]);
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && G.decls.has(to)) out.add(to);
  }
  return [...out];
}
const need = new Set(direct.keys());
{ const q = [...need]; while (q.length) { const v = q.shift(); for (const w of intraEdges(v)) if (!need.has(w)) { need.add(w); q.push(w); } } }

// ── 分子: Found/ の `.src`
const ITEM_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*$`);
const PART_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*[,—-]`);
const SRC_RE = /\.src[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?item\s*:=\s*"([^"]*)"/g;
function walk(d, out = []) { for (const e of readdirSync(d)) { const p = join(d, e); if (statSync(p).isDirectory()) walk(p, out); else out.push(p); } return out; }
const LEAN = join(ROOT, 'lean', 'ABC3');
const done = new Map();
const partial = new Map();
for (const f of walk(LEAN).filter((p) => p.endsWith('.lean'))) {
  const rel = relative(LEAN, f).split('\\').join('/');
  if (rel.split('/')[0] !== 'Found') continue;
  const text = readFileSync(f, 'utf8');
  for (const m of text.matchAll(SRC_RE)) {
    if (m[1] !== 'GenEll') continue;
    const im = ITEM_RE.exec(m[2]);
    if (!im) {
      const pm = PART_RE.exec(m[2]);
      if (pm) {
        const pk = `${pm[1]} ${pm[2]}`;
        if (!partial.has(pk)) partial.set(pk, new Set());
        partial.get(pk).add(m[2]);
      }
      continue;
    }
    const k = `${im[1]} ${im[2]}`;
    if (!done.has(k)) done.set(k, new Set());
    done.get(k).add(rel);
  }
}

const doneInNeed = [...done.keys()].filter((k) => need.has(k)).sort();
const doneOutside = [...done.keys()].filter((k) => !need.has(k)).sort();
const remaining = [...need].filter((k) => !done.has(k));
const bySec = (k) => k.split(' ')[1].split('.')[0];
const cmp = (a, b) => {
  const pa = a.split(' ')[1].split('.').map(Number), pb = b.split(' ')[1].split('.').map(Number);
  return pa[0] - pb[0] || pa[1] - pb[1] || (pa[2] ?? 0) - (pb[2] ?? 0);
};

if (process.argv.includes('--json')) {
  writeFileSync(join(ROOT, 'ResearchPaper', 'genell-needed.json'), JSON.stringify({
    generated: 'node tools/genell-progress.mjs --json',
    limits: [
      '分母は測定値であって定義ではない。番号の無い依存を数えていないので下界寄り。',
      '一方「cf.」の案内も引用として数えているので上に振れる方向もある。桁を見るための数である。',
      'pdftotext 経由なので項目名は壊れうる(事実2)。逐語には使えない。',
      '★本体の終端は「次の宣言行 または 次の Section N: 見出し」。frdi-progress は節見出しで切らない。',
      '分子は .src の形だけを見る。locator の検証は check.mjs の G1 が行う。',
    ],
    needed: [...need].sort(cmp).map((k) => ({
      item: k, page: G.decls.get(k).page, section: bySec(k),
      citedBy: [...(direct.get(k) ?? [])], done: done.has(k),
    })),
    total: need.size, done: doneInNeed.length,
  }, null, 1) + '\n');
  console.log('書き出し: ResearchPaper/genell-needed.json');
}

console.log(`★ゴール進捗: [GenEll] の必要分 ${doneInNeed.length} / ${need.size} 件 (${(doneInNeed.length / need.size * 100).toFixed(0)}%)`);
if (missingConsumers.length) console.log(`  ★.txt が無い需要側: ${missingConsumers.join(', ')} —— 分母はそのぶん小さく出ている`);
console.log(`  直接名指し ${direct.size} 件 → 推移閉包 ${need.size} 件 / GenEll 全 ${G.decls.size} 件`);
if (doneOutside.length) console.log(`  ★必要分の外に実装したもの: ${doneOutside.join(' / ')}`);
{
  const partInNeed = [...partial.keys()].filter((k) => need.has(k) && !done.has(k)).sort(cmp);
  if (partInNeed.length) {
    console.log(`  ★条つき \`.src\` はあるが命題全体の \`.src\` が無いもの ${partInNeed.length} 件`);
    console.log('     ★★これは**数に入らない**。命題全体が揃ったときに条なしの `.src` を 1 個足すこと。');
    for (const k of partInNeed) console.log(`     ${k} — 条つき ${partial.get(k).size} 個 (${[...partial.get(k)].join(' / ')})`);
  }
}

console.log('\n節別:');
for (const s of [...new Set([...need].map(bySec))].sort((a, b) => +a - +b)) {
  const all = [...need].filter((k) => bySec(k) === s);
  const d = all.filter((k) => done.has(k));
  const bar = '█'.repeat(d.length) + '·'.repeat(all.length - d.length);
  console.log(`  §${s}  ${String(d.length).padStart(2)}/${String(all.length).padStart(2)}  ${bar}`);
}

console.log('\n★次に着手する候補(必要分のうち未実装、物理ページ順):');
for (const k of remaining.sort((a, b) => G.decls.get(a).page - G.decls.get(b).page).slice(0, 12)) {
  const c = direct.get(k);
  console.log(`  p.${String(G.decls.get(k).page).padStart(3)}  ${k.padEnd(20)}${c ? `  ← ${[...c].join(',')} が直接引く` : '  (推移閉包で入ったもの)'}`);
}
if (remaining.length > 12) console.log(`  … 他 ${remaining.length - 12} 件`);
