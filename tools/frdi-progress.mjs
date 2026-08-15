// ゴール「[FrdI] を IUT が要る範囲で Found/ に完成させる」の進捗を機械で数える。
//
// ★この道具が答える問い: **54 件のうち何件が実装済みか。**
//
// ★分母(54)の出し方と、その限界(先に書く):
//   分子でも分母でもなく、まず **分母が測定値である**ことを明示する。
//   (1) IUTchI–IV と EtTh の本文から `[FrdI], Kind N.M` の形の引用を集める(直接需要)
//   (2) それを FrdI 内の同一論文内参照で推移閉包を取る
//   ★限界1: **番号の無い依存を数えていない。** §0 の語彙、原文が「immediate」と書く段、
//            本文中で名前だけで参照される概念は辺にならない。→ **分母は下界寄り**である。
//   ★限界2: `pdftotext` 経由なので項目名が壊れうる(事実2)。**逐語には使えない。**
//   ★限界3: 「cf.」の案内も引用として数えている。→ **分母を上に振らせる方向**。
//   限界1 と 限界3 は逆向きなので、分母は「桁を見るための数」である。
//
// ★分子の出し方:
//   `Found/` の `.src` を読み、`item` の先頭が `Kind N.M` の形のものを実装済みとする。
//   ★`check.mjs` の G1 が locator(論文登記・ページ範囲・sectionId 実在)を検証しているので、
//     ここでは形だけ見る。**G1 を通っていない `.src` はゲートが落とす。**
//
// ★★分子の規則(2026-08-15 明文化。子からの指摘を受けて親が決めた):
//   **`.src` は「その原典項目を**完全に**実装した」という主張である。**
//   この道具は部分実装と完全実装を区別しない。だから区別は `.src` を付けるかどうかで表す。
//   - 項目の一部だけ実装した段階では `.src` を **付けない**。
//   - 項目へ向けた**構成**(例: `MonoidOn.charOn`)にも付けない。原典の項目ではないから。
//   ★この規則だと、大きい項目(例: Proposition 1.5 は (i)(ii)(iii) と 22 条)は
//     全部終わるまで数に出ない。**それでよい。指標は完了を測るものであって労力ではない。**
//   ★「部分実装」の状態を器具に持たせない理由: 状態を増やすと、どこまでで「部分」を
//     名乗ってよいかの線引きが要り、そこが緩む。**付けるか付けないかの2値のほうが硬い。**
//   ★逆インセンティブは働かない: 数を上げる唯一の方法は項目を完成させることである。
//
// 使い方: node tools/frdi-progress.mjs [--json]

import { readFileSync, writeFileSync, readdirSync, statSync, existsSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const SRC = join(ROOT, 'ResearchPaper', '0_Source');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';
const FILE_OF = JSON.parse(readFileSync(join(ROOT, 'tools', '_fileof.json'), 'utf8'));
// ★タグに数字を許すこと。`[Mzk15]` を masking し損ねると偽の同一論文内参照になる(2026-08-15 の実測)。
const TAGRE = new RegExp(`\\[([A-Za-z][A-Za-z0-9]*)\\],?\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');

function load(tag) {
  const f = FILE_OF[tag];
  const p = f ? join(SRC, `${f}.txt`) : null;
  if (!p || !existsSync(p)) return null;
  const lines = readFileSync(p, 'utf8').split(/\r?\n/);
  const pageOf = new Array(lines.length).fill(0);
  { let pg = 0; for (let i = 0; i < lines.length; i++) { const m = lines[i].match(/^===== \[page (\d+)\]/); if (m) pg = Number(m[1]); pageOf[i] = pg; } }
  // ★宣言行 = 行頭 + 番号 + ピリオド + (行末 または 括弧つき表題)。緩めると継続行を誤検出する。
  const declRe = new RegExp(`^[ \\t]*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\.[ \\t]*(?:$|\\()`);
  const list = [];
  for (let i = 0; i < lines.length; i++) {
    const m = lines[i].match(declRe);
    if (m) list.push({ name: `${m[1]} ${m[2]}`, line: i, page: pageOf[i] });
  }
  for (let i = 0; i < list.length; i++) list[i].end = i + 1 < list.length ? list[i + 1].line : lines.length;
  const decls = new Map();
  for (const d of list) if (!decls.has(d.name)) decls.set(d.name, d);
  return { decls, lines };
}

const F = load('FrdI');
if (!F) { console.error('★[FrdI] の .txt が無い。0_Source を確認すること。'); process.exit(1); }

// ── 分母: 需要の推移閉包
const CONSUMERS = ['IUTchI', 'IUTchII', 'IUTchIII', 'IUTchIV', 'EtTh'];
const direct = new Map();   // "Kind N.M" -> Set(引いた論文)
const missingConsumers = [];
for (const tag of CONSUMERS) {
  const P = load(tag);
  if (!P) { missingConsumers.push(tag); continue; }
  for (const m of P.lines.join('\n').matchAll(TAGRE)) {
    if (m[1] !== 'FrdI' && m[1] !== 'Mzk14') continue;   // FrdI の引用タグは2通り
    const k = `${m[2]} ${m[3]}`;
    if (!F.decls.has(k)) continue;                        // FrdI に宣言が無いものは辿れない
    if (!direct.has(k)) direct.set(k, new Set());
    direct.get(k).add(tag);
  }
}

const bodyOf = (name) => { const d = F.decls.get(name); return F.lines.slice(d.line, d.end).join('\n'); };
function intraEdges(name) {
  const body = bodyOf(name);
  const spans = []; const out = new Set();
  for (const m of body.matchAll(TAGRE)) spans.push([m.index, m.index + m[0].length]);
  const ire = new RegExp(`(${KIND})\\s+(\\d+(?:\\.\\d+)+)`, 'g');
  for (const m of body.matchAll(ire)) {
    if (spans.some(([a, b]) => m.index >= a && m.index < b)) continue;
    const to = `${m[1]} ${m[2]}`;
    if (to !== name && F.decls.has(to)) out.add(to);
  }
  return [...out];
}
const need = new Set(direct.keys());
{ const q = [...need]; while (q.length) { const v = q.shift(); for (const w of intraEdges(v)) if (!need.has(w)) { need.add(w); q.push(w); } } }

// ── 分子: Found/ の `.src` が指す原典項目
// ★★2026-08-16 修正: **項目名が丸ごと `Kind N.M` のものだけを数える。**
//   以前は先頭一致(`^…`)だったので `item := "Proposition 1.9, (i)"` が
//   キー "Proposition 1.9" に潰れ、★**条が 1 つでもあれば命題全体が実装済みに数えられた。**
//   これは上の docstring が宣言している規則(「全部終わるまで数に出ない」)と**食い違っていた**。
//   ★文脈を持たない検証役の指摘で判明した。器具を docstring の意図に戻す。
const ITEM_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*$`);
// 条つき `.src`(`"Proposition 1.9, (i)"` 等)を拾う——**数えないが、見えるようにする。**
const PART_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*[,—-]`);
const SRC_RE = /\.src[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?item\s*:=\s*"([^"]*)"/g;
function walk(d, out = []) { for (const e of readdirSync(d)) { const p = join(d, e); if (statSync(p).isDirectory()) walk(p, out); else out.push(p); } return out; }
const LEAN = join(ROOT, 'lean', 'ABC3');
const done = new Map();     // "Kind N.M" -> Set(ファイル)  ★完全実装の主張
const partial = new Map();  // "Kind N.M" -> Set(条つき item 文字列)  ★数えない
for (const f of walk(LEAN).filter((p) => p.endsWith('.lean'))) {
  const rel = relative(LEAN, f).split('\\').join('/');
  if (rel.split('/')[0] !== 'Found') continue;
  const text = readFileSync(f, 'utf8');
  for (const m of text.matchAll(SRC_RE)) {
    if (m[1] !== 'FrdI') continue;
    const im = ITEM_RE.exec(m[2]);
    if (!im) {
      const pm = PART_RE.exec(m[2]);                     // 条つき —— 記録だけ
      if (pm) {
        const pk = `${pm[1]} ${pm[2]}`;
        if (!partial.has(pk)) partial.set(pk, new Set());
        partial.get(pk).add(m[2]);
      }
      continue;                                          // §0 の語彙などは番号項目でない
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
  writeFileSync(join(ROOT, 'ResearchPaper', 'frdi-needed.json'), JSON.stringify({
    generated: 'node tools/frdi-progress.mjs --json',
    limits: [
      '分母は測定値であって定義ではない。番号の無い依存(§0 の語彙、原文が immediate と書く段)を数えていないので下界寄り。',
      '一方「cf.」の案内も引用として数えているので上に振れる方向もある。両者は逆向きで、分母は桁を見るための数である。',
      'pdftotext 経由なので項目名は壊れうる(事実2)。逐語には使えない。',
      '分子は .src の形だけを見る。locator の検証は check.mjs の G1 が行う。',
    ],
    needed: [...need].sort(cmp).map((k) => ({
      item: k, page: F.decls.get(k).page, section: bySec(k),
      citedBy: [...(direct.get(k) ?? [])], done: done.has(k),
    })),
    total: need.size, done: doneInNeed.length,
  }, null, 1) + '\n');
  console.log(`書き出し: ResearchPaper/frdi-needed.json`);
}

console.log(`★ゴール進捗: [FrdI] の必要分 ${doneInNeed.length} / ${need.size} 件 (${(doneInNeed.length / need.size * 100).toFixed(0)}%)`);
if (missingConsumers.length) console.log(`  ★.txt が無い需要側: ${missingConsumers.join(', ')} —— 分母はそのぶん小さく出ている`);
console.log(`  直接名指し ${direct.size} 件 → 推移閉包 ${need.size} 件 / FrdI 全 ${F.decls.size} 件`);
if (doneOutside.length) console.log(`  ★必要分の外に実装したもの: ${doneOutside.join(' / ')}`);
{
  const partInNeed = [...partial.keys()].filter((k) => need.has(k) && !done.has(k)).sort(cmp);
  if (partInNeed.length) {
    console.log(`  ★条つき \`.src\` はあるが命題全体の \`.src\` が無いもの ${partInNeed.length} 件`);
    console.log(`     ★★これは**数に入らない**。`);
    console.log(`     ★命題全体が揃ったときに \`item := "Kind N.M"\`(条なし)の \`.src\` を 1 個足すこと。`);
    for (const k of partInNeed) console.log(`     ${k} — 条つき ${partial.get(k).size} 個`);
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
for (const k of remaining.sort((a, b) => F.decls.get(a).page - F.decls.get(b).page).slice(0, 12)) {
  const c = direct.get(k);
  console.log(`  p.${String(F.decls.get(k).page).padStart(3)}  ${k.padEnd(20)}${c ? `  ← ${[...c].join(',')} が直接引く` : '  (推移閉包で入ったもの)'}`);
}
if (remaining.length > 12) console.log(`  … 他 ${remaining.length - 12} 件`);
