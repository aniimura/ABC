import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const e = j.blocked.find(b => b.key.startsWith('[GenEll] Theorem 3.8(悪い素点の惰性'));
if (!e) { console.log('not found'); process.exit(1); }
e.closed2026_09_02 = '★★★★★★★★★★★★★★★★★★★★**この葉は閉じた（2026-09-02、第 1384）**。'
  + '`exists_h2_h1_unipotent` は `SSCurve.exists_h2_h1_unipotent_of_multRed`'
  + '（`Found/GenEll/AlphaUnipotent.lean`、無条件）で埋まり、'
  + '`Skeleton/GenEll/AlphaLocalData.lean` は sorry 0 になった。'
  + '☆第 1354 で「mathlib に完備離散付値環の有限拡大の理論が無い」と記録したブロッカーは'
  + '第 1357-1372 で自前で作って塞ぎ、第 1373-1383 で配管を通した。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded');
