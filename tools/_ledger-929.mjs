// 2026-09-01（第 926-929）—— 非分裂の降下が「捻った対に分裂の連鎖を当てる」だけになった。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Final20260901 = {
  measuredAt: '2026-09-01',
  block: '第 926-929',
  supersedes: 'lemma35Twist20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  splitChain:
    '★★★★分裂の場合: stableLine_is_mu_of_coprime(906) → ' +
    'tateParam_quot_velu_j_dvr(914) / tateParam_quot_veluCurve_dvr(927) → ' +
    'minDeltaExp_eq_mul_of_veluMu(904) → lemma_3_5_velu_bad_delta(903)。',
  nonsplitChain:
    '★★★★非分裂の場合: quadTwist_veluQuotientFull(925) → ' +
    'j_veluCurve_variableChange(926) → tateParam_quot_veluCurve_dvr(927) → ' +
    'minDeltaExp_eq_mul_of_nonsplit(929)。' +
    '要石は 929——**捻った対で関係が成り立てばもとの対でも成り立つ**。',
  toolsReady: [
    'quadTwist / quadTwist_j / quadTwist_map / veluCurve_quadTwist / ' +
      'quadTwist_veluQuotientFull（920・923-925）',
    'splits_quadratic_of_root（921）——分裂性は剰余体に根が 1 つ',
    'isSquare_mul_of_not_isSquare（922）——非平方×非平方 = 平方',
    'minDeltaExp_eq_of_j_eq / minDeltaExp_eq_mul_of_twist（919）',
    'jExp_quadTwist / minDeltaExp_quadTwist（928）',
    'minDeltaExp_eq_mul_of_nonsplit（929）',
  ],
  remaining: [
    '☆(1) 「d を非平方単数に取れば quadTwist は p で分裂乗法還元」——' +
      '921 + 922 を並べる。要る配管は `integralModel R (quadTwist W d)` と ' +
      '`quadTwist (integralModel R W) d` の同定だけである。',
    '☆(2) 捻った対の `hv`・`hw`（μ_l の和との同定）——' +
      '分裂の場合とまったく同じ帳簿。925 が曲線の等式を、926 が変換則を与えている。',
    '☆(3) 上を並べて `lemma_3_5_velu_bad_delta`(903) に流す。',
  ],
  beyondLemma35:
    '☆Lemma 3.5 が閉じれば Lemma 3.7 → Theorem 3.8 → Corollary 4.3 / 4.4 は' +
    '既存の導出（第 367、Theorem 3.8 を仮説として受けた形で証明済み）が動く。' +
    '★Theorem 2.1（§2）だけは独立で、[NCBelyi]（non-critical Belyi maps）という' +
    '別論文の形式化を要する——本プロジェクトにはまだ無い。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Final20260901 を書いた');
