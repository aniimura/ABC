// 13 件の論文を papers.json へ登記する(第 1465)。
//
// これらは tools/_fileof.json には既にあり(依存グラフの引用解決に使われていた)、
// ResearchPaper/papers.json には無かった。0_Source に PDF/txt があることを確認し、
// 物理 p.1-3 を pdftotext で実測して pdfPages・pageOffset・author・(測れた場合のみ)year を埋める。
//
// ★notationRisk は "unmeasured" にする——ABSLean と同じ扱いである。
//   IUT 系の他論文(IUTchI-III・EtTh・AbsTopIII 等)はどれも high(オーバーバー・script・
//   下線の脱落)なので、これらも同種の危険を抱えている可能性が高いが、
//   実際に 260dpi で目視するまでは「high」と断定しない。
//
// ★year を書かなかったもの(AbsAnab・CanLift・CombCusp・AbsSect・MT・NodNon)は、
//   物理 p.1-3 に著者名に隣接する年が無かった、または見つかった年が
//   「2000 Mathematical Subject Classification」(MSC のスキーム年であって発表年ではない)
//   の誤検出だったため。推測で埋めない。
import fs from 'fs';

const P = 'D:/Math_ABC3/ResearchPaper/papers.json';
const j = JSON.parse(fs.readFileSync(P, 'utf8'));
fs.writeFileSync(P + '.bak-2026-09-03', JSON.stringify(j, null, 1) + '\n', 'utf8');

const UNMEASURED_NOTE =
  '★★★まだ目視していない。テキスト層は取れているが、どの記号が壊れるかは未測定' +
  '(GenEll では行列の順序が入れ替わる実害、IUTchI-III・EtTh・AbsTopIII ではオーバーバー・' +
  'script 文字・下線の脱落が実害になった)。逐語を使う前に必ず 260dpi で目視すること。';

const ADD = {
  IUTchIV: {
    title: 'Inter-universal Teichmuller Theory IV: Log-volume Computations and Set-theoretic Foundations',
    file: 'Inter-universal Teichmuller Theory IV', author: 'Shinichi Mochizuki', year: 2020,
    pdfPages: 87, pageOffset: 0,
  },
  AbsTopI: {
    title: 'Topics in Absolute Anabelian Geometry I: Generalities',
    file: 'Topics in Absolute Anabelian Geometry I', author: 'Shinichi Mochizuki', year: 2012,
    pdfPages: 83, pageOffset: 0,
  },
  AbsTopII: {
    title: 'Topics in Absolute Anabelian Geometry II: Decomposition Groups and Endomorphisms',
    file: 'Topics in Absolute Anabelian Geometry II', author: 'Shinichi Mochizuki', year: 2013,
    pdfPages: 76, pageOffset: 0,
  },
  FrdII: {
    title: 'The Geometry of Frobenioids II: Poly-Frobenioids',
    file: 'The Geometry of Frobenioids II', author: 'Shinichi Mochizuki', year: 2008,
    pdfPages: 72, pageOffset: 0,
  },
  AbsAnab: {
    title: 'Absolute Anabelian Geometry',
    file: 'Absolute Anabelian Geometry', author: 'Shinichi Mochizuki',
    pdfPages: 45, pageOffset: 0,
  },
  CanLift: {
    title: 'Canonical Liftings',
    file: 'Canonical Liftings', author: 'Shinichi Mochizuki',
    pdfPages: 34, pageOffset: 0,
  },
  CombGC: {
    title: 'A Combinatorial Version of the Grothendieck Conjecture',
    file: 'Combinatorial Grothendieck Conjecture', author: 'Shinichi Mochizuki', year: 2007,
    pdfPages: 29, pageOffset: 0,
  },
  CombCusp: {
    title: 'Combinatorial Cuspidalization of Hyperbolic Curves',
    file: 'Combinatorial Cuspidalization', author: 'Shinichi Mochizuki',
    pdfPages: 65, pageOffset: 0,
  },
  AbsSect: {
    title: 'Galois Sections in Absolute Anabelian Geometry',
    file: 'Galois Sections', author: 'Shinichi Mochizuki',
    pdfPages: 25, pageOffset: 0,
  },
  MT: {
    title: 'An Introduction to p-adic Teichmuller Theory',
    file: 'An Introduction to p-adic Teichmuller Theory', author: 'Shinichi Mochizuki',
    pdfPages: 51, pageOffset: 0,
  },
  HASurII: {
    title: 'A Survey of the Hodge-Arakelov Theory of Elliptic Curves II',
    file: 'A Survey of the Hodge-Arakelov Theory of Elliptic Curves II',
    author: 'Shinichi Mochizuki', year: 2000,
    pdfPages: 34, pageOffset: 0,
  },
  NodNon: {
    title: 'Topics Surrounding the Combinatorial Anabelian Geometry of Hyperbolic Curves I: '
         + 'Inertia Groups and Profinite Dehn Twists',
    file: 'Combinatorial Anabelian Topics I', author: 'Yuichiro Hoshi and Shinichi Mochizuki',
    pdfPages: 154, pageOffset: 0,
  },
  Config: {
    title: 'The Algebraic and Anabelian Geometry of Configuration Spaces',
    file: 'The Algebraic and Anabelian Geometry of Configuration Spaces',
    author: 'Shinichi Mochizuki and Akio Tamagawa', year: 2007,
    pdfPages: 49, pageOffset: 0,
  },
};

let n = 0;
for (const [tag, e] of Object.entries(ADD)) {
  if (j.papers[tag]) { console.log(`  ★既に登記済み: ${tag}。飛ばす`); continue; }
  const src = `D:/Math_ABC3/ResearchPaper/0_Source/${e.file}`;
  if (!fs.existsSync(`${src}.txt`) || !fs.existsSync(`${src}.pdf`)) {
    console.log(`  ★0_Source に PDF/txt が揃っていない: ${tag} (${e.file})。飛ばす`);
    continue;
  }
  j.papers[tag] = { ...e, notationRisk: 'unmeasured', notationNotes: UNMEASURED_NOTE, verifiedPages: [] };
  n++;
  console.log(`  + ${tag.padEnd(9)} pages=${e.pdfPages} offset=${e.pageOffset}` +
              (e.year ? ` year=${e.year}` : ' year=(未測定)'));
}
fs.writeFileSync(P, JSON.stringify(j, null, 1) + '\n', 'utf8');
console.log(`papers.json: ${n} 件登記した(合計 ${Object.keys(j.papers).length})`);
