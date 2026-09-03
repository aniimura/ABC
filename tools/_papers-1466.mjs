// 次の並列作業のために LocProP・CorrHyp を papers.json へ登記する(第 1466)。
// 両方とも tools/_fileof.json には既にあった(依存グラフの引用解決に使われていた)が、
// ResearchPaper/papers.json には無く、check.mjs の G1(pdfPage 範囲検査)が通らない状態だった。
import fs from 'fs';

const P = 'D:/Math_ABC3/ResearchPaper/papers.json';
const j = JSON.parse(fs.readFileSync(P, 'utf8'));
fs.writeFileSync(P + '.bak-2026-09-04', JSON.stringify(j, null, 1) + '\n', 'utf8');

const UNMEASURED_NOTE =
  '★★★まだ目視していない。テキスト層は取れているが、どの記号が壊れるかは未測定。' +
  '逐語を使う前に必ず 260dpi で目視すること。';

const ADD = {
  LocProP: {
    title: 'The Local Pro-p Anabelian Geometry of Curves',
    file: 'The Local Pro-p Anabelian Geometry of Curves', author: 'Shinichi Mochizuki',
    pdfPages: 103, pageOffset: 0,
  },
  CorrHyp: {
    title: 'Correspondences on Hyperbolic Curves',
    file: 'Correspondences on Hyperbolic Curves', author: 'Shinichi Mochizuki',
    pdfPages: 18, pageOffset: 0,
  },
};

let n = 0;
for (const [tag, e] of Object.entries(ADD)) {
  if (j.papers[tag]) { console.log(`  ★既に登記済み: ${tag}`); continue; }
  const src = `D:/Math_ABC3/ResearchPaper/0_Source/${e.file}`;
  if (!fs.existsSync(`${src}.txt`) || !fs.existsSync(`${src}.pdf`)) {
    console.log(`  ★PDF/txt が揃っていない: ${tag}`); continue;
  }
  j.papers[tag] = { ...e, notationRisk: 'unmeasured', notationNotes: UNMEASURED_NOTE, verifiedPages: [] };
  n++;
  console.log(`  + ${tag.padEnd(9)} pages=${e.pdfPages} offset=${e.pageOffset} year=(未測定)`);
}
fs.writeFileSync(P, JSON.stringify(j, null, 1) + '\n', 'utf8');
console.log(`papers.json: ${n} 件登記した(合計 ${Object.keys(j.papers).length})`);
