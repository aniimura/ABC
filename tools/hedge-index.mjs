// 原文の「省略の合図」(immediately / formally / routine / one verifies …)を数え、
// 項目ごとに地図にする道具。
//
// ★動機(2026-08-20)
//   原典は証明を「it follows immediately」「a routine verification」で畳むことが多い。
//   我々はそれを 1 つずつ開いてきたが、**どこにいくつ残っているかを見ていなかった**。
//   CLAUDE.md の「進め方」は「必要なものが出たらスケルトンを足す」だが、
//   ★**省略の合図は「これから必要になるものの予告」である**——先に数えておけば、
//   分解(chain)を書く前に規模の当たりが付く。
//
// ★出力
//   1) 合図の語ごとの件数(論文全体)
//   2) 項目(Proposition N.M 等)ごとの合図の数と、我々の実装状況(条なし .src / 条つき / 未)
//   3) 未着手の項目のうち合図が多いもの = 「畳まれた量が多い」順の作業候補
//
// ★限界
//   * 合図の語が必ずしも省略とは限らない(定義の中の "clearly" など)。
//   * 「項目に属する」の判定は**直前の見出し**で行う近似である。
//   * 原文 txt は gitignore 下にあるので、無ければ静かに終わる(CI で落とさない)。
//
// 使い方: node tools/hedge-index.mjs [--paper FrdI] [--json] [--item "Proposition 1.10"]

import { readFileSync, existsSync, readdirSync, statSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..');

const PAPERS = {
  FrdI: 'The Geometry of Frobenioids I.txt',
  FrdII: 'The Geometry of Frobenioids II.txt',
};

const args = process.argv.slice(2);
const asJson = args.includes('--json');
const paperIdx = args.indexOf('--paper');
const paper = paperIdx >= 0 ? args[paperIdx + 1] : 'FrdI';
const itemIdx = args.indexOf('--item');
const onlyItem = itemIdx >= 0 ? args[itemIdx + 1] : null;

/** 合図の語。★分類は「実際に開けてみた経験」に基づく(README を見よ)。 */
const HEDGES = [
  { key: 'routine', re: /\broutine\b/i, note: '★最大級。畳まれた量がいちばん多い' },
  { key: 'formally', re: /\bformal(?:ly)?\b/i, note: '新しい材料は要らない。手数はある' },
  { key: 'immediately', re: /\bimmediat(?:e|ely)\b/i, note: '玉石混交。要点検' },
  { key: 'one verifies', re: /\b(?:one|One)\s+(?:verifies|checks|computes)\b/, note: '短いことが多い' },
  { key: 'easily', re: /\beasil(?:y|ier)\b/i, note: '短いことが多い' },
  { key: 'clearly', re: /\bclearl(?:y)\b/i, note: '短いことが多い' },
  { key: 'well-known', re: /\bwell-known\b/i, note: '★外部の在庫を指す。mathlib を測ること' },
  { key: 'similarly', re: /\b(?:similarly|in a similar (?:way|manner))\b/i, note: '前の議論の再演' },
];

// ★見出しの判定: 番号のあとに '.' が来て、その行がそこで終わる(または '(' で題が続く)もの
// だけを見出しとする。'Proposition 5.5, (iii), below]' のような**前方参照を拾わない**ため。
const ITEM_RE = /^(Proposition|Theorem|Definition|Corollary|Lemma|Example|Remark)\s+([0-9]+\.[0-9]+(?:\.[0-9]+)?)\.(?:\s*$|\s*\()/;
const PAGE_RE = /^=====\s*\[page\s+(\d+)\]\s*=====/;

function readSource() {
  const f = PAPERS[paper];
  if (!f) return null;
  const p = path.join(REPO, 'ResearchPaper', '0_Source', f);
  if (!existsSync(p)) return null;
  return readFileSync(p, 'utf8');
}

/** Lean 側の `.src` を集める: item 文字列 → 条なしか条つきか */
function scanLean() {
  const roots = [path.join(REPO, 'lean', 'ABC3')];
  const items = new Map(); // "Proposition 1.10" → {bare:bool, qualified:number}
  const walk = (d) => {
    for (const e of readdirSync(d)) {
      const p = path.join(d, e);
      const st = statSync(p);
      if (st.isDirectory()) walk(p);
      else if (e.endsWith('.lean')) {
        const s = readFileSync(p, 'utf8');
        const re = /item\s*:=\s*"([^"]*)"/g;
        let m;
        while ((m = re.exec(s)) !== null) {
          const raw = m[1];
          const head = raw.match(/^(Proposition|Theorem|Definition|Corollary|Lemma|Example|Remark)\s+[0-9]+\.[0-9]+(?:\.[0-9]+)?/);
          if (!head) continue;
          const key = head[0];
          const rec = items.get(key) ?? { bare: false, qualified: 0 };
          // 条なし = item がちょうど見出しだけ(「, (ii)」や「 — …」が付かない)
          if (raw.trim() === key) rec.bare = true;
          else rec.qualified += 1;
          items.set(key, rec);
        }
      }
    }
  };
  for (const r of roots) if (existsSync(r)) walk(r);
  return items;
}

const src = readSource();
if (src === null) {
  console.log(`(原文 txt が無いので測れない: ResearchPaper/0_Source/${PAPERS[paper] ?? paper})`);
  process.exit(0);
}

const lines = src.split(/\r?\n/);
let page = null;
let item = null;
const perItem = new Map(); // key → {page, hits: [{key, line, page}], total}
const perHedge = new Map();

for (let i = 0; i < lines.length; i++) {
  const L = lines[i];
  const pm = L.match(PAGE_RE);
  if (pm) { page = Number(pm[1]); continue; }
  const im = L.match(ITEM_RE);
  if (im) {
    item = `${im[1]} ${im[2]}`;
    if (!perItem.has(item)) perItem.set(item, { page, hits: [], total: 0 });
    continue;
  }
  for (const h of HEDGES) {
    if (h.re.test(L)) {
      perHedge.set(h.key, (perHedge.get(h.key) ?? 0) + 1);
      if (item) {
        const rec = perItem.get(item);
        rec.hits.push({ key: h.key, line: i + 1, page });
        rec.total += 1;
      }
    }
  }
}

const leanItems = scanLean();

const rows = [...perItem.entries()]
  .map(([k, v]) => {
    const li = leanItems.get(k);
    const state = li?.bare ? '済' : li?.qualified ? `条つき${li.qualified}` : '未';
    return { item: k, page: v.page, total: v.total, state, hits: v.hits };
  })
  .filter((r) => r.total > 0)
  .filter((r) => !onlyItem || r.item === onlyItem);

if (asJson) {
  console.log(JSON.stringify({ paper, perHedge: Object.fromEntries(perHedge), rows }, null, 1));
  process.exit(0);
}

console.log(`★[${paper}] 原文の「省略の合図」の地図\n`);
console.log('-- 語ごとの件数(論文全体)');
for (const h of HEDGES) {
  const n = perHedge.get(h.key) ?? 0;
  if (n) console.log(`   ${h.key.padEnd(13)} ${String(n).padStart(4)}   ${h.note}`);
}
console.log('\n-- 項目ごと(合図の多い順、上位 25)');
console.log('   状態  合図  物理p  項目                        内訳');
for (const r of rows.slice().sort((a, b) => b.total - a.total).slice(0, 25)) {
  const brk = Object.entries(
    r.hits.reduce((m, h) => ({ ...m, [h.key]: (m[h.key] ?? 0) + 1 }), {})
  ).map(([k, n]) => `${k}×${n}`).join(' ');
  console.log(
    `   ${r.state.padEnd(6)}${String(r.total).padStart(3)}  ${String(r.page ?? '?').padStart(5)}  ${r.item.padEnd(26)}${brk}`
  );
}

const undone = rows.filter((r) => r.state === '未' || r.state.startsWith('条つき'));
const sum = undone.reduce((a, r) => a + r.total, 0);
console.log(`\n-- 未実装(条つき含む)の項目に残る合図: ${sum} 件 / ${undone.length} 項目`);
console.log('   ★これが「まだ開いていない省略」の下界である(合図の無い省略は写らない)。');
console.log('\n★使い方: 項目に着手する前に `--item "Proposition 1.10"` で内訳を見て、');
console.log('  合図 1 つを分解(frdi-decomposition.json)の節点 1 つに対応させると、');
console.log('  「畳まれた量」を最初に見積もれる。');
