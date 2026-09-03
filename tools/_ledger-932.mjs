// 2026-09-01（第 930-932）—— jExp の橋は分裂性を要求しないことが分かった。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35JBridge20260901 = {
  measuredAt: '2026-09-01',
  block: '第 930-932',
  supersedes: 'lemma35Final20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  discovery:
    '★★★★**`jExp` の橋は分裂性を要求しない**（第 932 ' +
    '`jExp_eq_neg_vAdd_of_j_tateCurveAt`）。' +
    '「`j(E ⊗ Lv) = j(E_q ⊗ Lv)`」だけを受けて `jExp p E = −v(q)` が出る。' +
    '☆したがって非分裂の場合も、捻りを経由せずに扱える可能性がある。',
  newTools: [
    'evalAdic_sub_constantCoeff_mem（930）——定数項を引けば I の元',
    'evalAdic_tateJinvSeries_eq_mul_unit（930）——1/j の評価は q × 単元' +
      '（tateJinvSeries = X·h で h の定数項が 1 だから）',
    'vAdd_evalAdic_tateJinvSeries（931）——v(1/j) = v(q)',
    'j_tateCurveAt_inv（932）——Tate 曲線の j の逆は 1/j の評価そのもの',
    'jExp_eq_neg_vAdd_of_j_tateCurveAt（932）——★v_p(j) = −v(q)',
  ],
  twoRoutesNow: {
    twist:
      '☆捻り経由（919-929）: 材料は全部揃っている。残るのは ' +
      '(1) d を非平方単数に取れば捻りは分裂、(2) 捻った対の hv・hw、(3) 組み立て。',
    jOnly:
      '☆★**j だけの経路（930-932）**: 非分裂の `E` に対しても ' +
      '`j(E ⊗ Lv) = j(E_q ⊗ Lv)` なる `q ∈ 𝔪` を作れれば、捻りは要らない。' +
      '★要るのは `q ↦ evalAdic tateJinvSeries q` の**全射性**' +
      '（`1/j ∈ I` を与えて `q` を解く。逐次近似／Hensel）。' +
      '☆単射性は第 875（`evalAdic_injective_of_coeff_one`）で済んでいる。',
  },
  recommendation:
    '★次のセッションは **j だけの経路**を先に測ることを勧める——' +
    '全射性 1 本で非分裂の場合が丸ごと消える見込みがある。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35JBridge20260901 を書いた');
