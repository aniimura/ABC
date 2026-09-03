// 2026-09-01（第 920-923）—— 非分裂の降下の材料が揃ったことを記録する。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Twist20260901 = {
  measuredAt: '2026-09-01',
  block: '第 920-923',
  supersedes: 'lemma35Remaining20260901b',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  materialsReady: [
    '`quadTwist` / `quadTwist_c₄` / `_c₆` / `_Δ` / `quadTwist_j`（920）——' +
      '捻りは c₄ ↦ d²c₄・c₆ ↦ d³c₆・Δ ↦ d⁶Δ なので j は不変。',
    '`splits_quadratic_of_root`（921）——分裂性は剰余体に根が 1 つあることに落ちる。',
    '`isSquare_mul_of_not_isSquare`（922）——有限体で非平方×非平方 = 平方。',
    '`veluCurve_quadTwist`（923）——veluCurve (W^d) (d²v) (d³w) = (veluCurve W v w)^d。',
    '`minDeltaExp_eq_mul_of_twist`（919）——j が同じ半安定な対で降りる。',
  ],
  remainingAssembly: [
    '☆(i) 非分裂の `E ⊗ Lv` に対し `d` を非平方単数に取り、`quadTwist` が ' +
      '`p` で**分裂**乗法還元をもつことを示す（921 + 922 を並べる）。',
    '☆(ii) 捻った点集合の Vélu の和が `d²v`・`d³w` になること' +
      '（`√d` の奇数冪は反転で安定な和の中で相殺する）。' +
      '★923 が曲線の水準の等式を与えているので、残るのは和の側の帳簿である。',
    '☆(iii) 上を並べて `minDeltaExp p E′ = l·minDeltaExp p E` を非分裂の場合にも出す。',
  ],
  splitChainComplete:
    '★★★★分裂の場合の連鎖は完結: stableLine_is_mu_of_coprime(906) → ' +
    'tateParam_quot_velu_j_dvr(914) → minDeltaExp_eq_mul_of_veluMu(904) → ' +
    'lemma_3_5_velu_bad_delta(903)。',
  beyondLemma35:
    '☆Lemma 3.5 が閉じたあとの残り: Lemma 3.7 → Theorem 3.8 → Corollary 4.3 / 4.4 は' +
    'Lemma 3.5 を消費する鎖であり、Theorem 2.1（§2）だけが独立で' +
    '[NCBelyi]（non-critical Belyi maps）という別論文を待っている。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Twist20260901 を書いた');
