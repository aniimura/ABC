// [CorrHyp] Found/ 完成の進捗を機械で数える。
//
// ★この道具が答える問い: **24項目のうち何件が Found/ で条なし `.src` を持つか。**
// ★長期ゴールの本体は `ResearchPaper/corrhyp-goal.md`。分母は `tools/paper-items.mjs CorrHyp`
//   と同じ 24 件(§1 5・§2 6・§3 3・§4 2・§5 7・§6 1)。
//
// ★`tools/genell-progress.mjs`/`tools/frdi-progress.mjs` と同じ分子の規則:
//   `.src` の `item` が **`Kind N.M` ちょうど**のものだけを実装済みとする。
//   条つき(`"Definition 1.1, (i)"`)は数えないが、見えるように印字する。
//
// ★分母は CorrHyp が単一論文の自己完結ゴールなので、GenEll のような
//   「下流論文からの需要の推移閉包」は不要——`paper-items.mjs` の 24 件全部が分母。
//
// 使い方: node tools/corrhyp-progress.mjs [--json]

import { readFileSync, readdirSync, statSync } from 'node:fs';
import { join, dirname, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const LEAN = join(ROOT, 'lean', 'ABC3');
const KIND = 'Definition|Proposition|Theorem|Corollary|Remark|Lemma|Example';

const SECTIONS = [
  { n: 1, items: ['Definition 1.1', 'Definition 1.2', 'Definition 1.3', 'Definition 1.4', 'Definition 1.5'] },
  { n: 2, items: ['Definition 2.1', 'Definition 2.2', 'Definition 2.3', 'Proposition 2.4', 'Theorem 2.5', 'Theorem 2.6'] },
  { n: 3, items: ['Definition 3.1', 'Proposition 3.2', 'Theorem 3.3'] },
  { n: 4, items: ['Lemma 4.1', 'Theorem 4.2'] },
  { n: 5, items: ['Lemma 5.1', 'Definition 5.2', 'Theorem 5.3', 'Lemma 5.4', 'Lemma 5.5', 'Lemma 5.6', 'Theorem 5.7'] },
  { n: 6, items: ['Theorem 6.1'] },
];
const ALL_ITEMS = SECTIONS.flatMap((s) => s.items);

const ITEM_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*$`);
const PART_RE = new RegExp(`^\\s*(${KIND})\\s+(\\d+(?:\\.\\d+)+)\\s*[,(（—-]`);
const SRC_RE = /\.src[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?item\s*:=\s*"([^"]*)"/g;

function walk(d, out = []) {
  for (const e of readdirSync(d)) {
    const p = join(d, e);
    if (statSync(p).isDirectory()) walk(p, out);
    else out.push(p);
  }
  return out;
}

const done = new Map(); // "Kind N.M" -> Set(ファイル)
const partial = new Map();
for (const f of walk(LEAN).filter((p) => p.endsWith('.lean'))) {
  const rel = relative(LEAN, f).split('\\').join('/');
  if (rel.split('/')[0] !== 'Found') continue;
  const text = readFileSync(f, 'utf8');
  for (const m of text.matchAll(SRC_RE)) {
    if (m[1] !== 'CorrHyp') continue;
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

const asJson = process.argv.includes('--json');
const result = { sections: [], totalDone: 0, totalItems: ALL_ITEMS.length };
for (const s of SECTIONS) {
  const items = s.items.map((it) => ({
    item: it,
    done: done.has(it),
    files: done.has(it) ? [...done.get(it)] : [],
    partialFiles: partial.has(it) ? [...partial.get(it)] : [],
  }));
  const nDone = items.filter((i) => i.done).length;
  result.totalDone += nDone;
  result.sections.push({ section: s.n, done: nDone, total: s.items.length, items });
}

if (asJson) {
  console.log(JSON.stringify(result, null, 2));
} else {
  console.log('★[CorrHyp] Found/ 完成の進捗(node tools/corrhyp-progress.mjs)');
  console.log();
  for (const s of result.sections) {
    console.log(`§${s.section}  —— ${s.done}/${s.total}`);
    for (const it of s.items) {
      const mark = it.done ? '✅' : (it.partialFiles.length ? '△' : '  ');
      console.log(`  ${mark} ${it.item}${it.done ? '  (' + it.files.join(', ') + ')' : ''}${it.partialFiles.length ? '  [条つき: ' + it.partialFiles.join('; ') + ']' : ''}`);
    }
  }
  console.log();
  const parts = result.sections.map((s) => `§${s.section} ${s.done}/${s.total}`).join(', ');
  console.log(`★ゴール表記: CorrHyp ${parts} —— 合計 ${result.totalDone}/${result.totalItems}`);
}
