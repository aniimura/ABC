// 2026-09-01（第 978）—— hcop の出どころ。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Hcop20260901 = {
  measuredAt: '2026-09-01',
  block: '第 978',
  supersedes: 'lemma35SmallHyps20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`vAdd_tateParam_eq_neg_jExp`（978）——' +
      'Tate 母数の付値は `−jExp`。第 932 を `q = tateParamR (E ⊗ Lv) h` に当てた形。',
    '★★★★★★★★★★★★`not_dvd_vAdd_tateParam_of_not_dvd_jExp`（978）——' +
      '原文の「`l` は局所高さと互いに素」を第 972 の `hcop` の形に直す。',
  ],
  howItWasDone:
    '★第 932 の仮説はすべて出た: ' +
    'Tate モデルの楕円性・`c₄ ≠ 0`・`j` の一致 → `tateModel_baseChange`（944）、' +
    '`1/j` の評価が `0` でない → `evalAdic_tateJinvSeries_eq_mul_unit`（`q·単元` だから）、' +
    '`E.j ≠ 0`・`E.c₄ ≠ 0` → `jExp p E < 0`（`j = 0` なら `jExp = 0`）。',
  remaining:
    '★残るのは **分裂性だけ**である（963 の二者択一＋整モデルの `IsCharNeTwoNF` 正規化、' +
    '非分裂は 938 → 929）。それが付けば 972 の仮説はすべて揃い、' +
    '`minDeltaExp_variableChange` で `E`・`E′` に戻して `IsMuAtBadPrimes` が閉じる。',
  note:
    '★944-978 の 35 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Hcop20260901 を書いた');
