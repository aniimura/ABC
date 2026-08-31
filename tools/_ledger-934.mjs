// 2026-09-01（第 933-934）—— 形式逆関数定理の存在の側が入った。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Inverse20260901 = {
  measuredAt: '2026-09-01',
  block: '第 933-934',
  supersedes: 'lemma35JBridge20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  newTools: [
    '★`evalAdic_surjective_of_coeff_one`（933）——形式逆関数定理の**存在の側**。' +
      '`f = X + 高次`、`t ∈ I` なら `f(q) = t` を満たす `q ∈ I` がある。' +
      '逐次近似（`q_{n+1} = q_n + (t − f(q_n))`）＋ IsPrecomplete ＋ IsHausdorff。' +
      '☆単射の側は第 875 で済んでいたので、これで両側が揃った。',
    '★`tateCurveAt_map_isElliptic` / `exists_tateParam_of_inv_j`（934）——' +
      '`1/j ∈ 𝔪` なら `j(E_q) = j(W)` なる Tate 母数 `q` が存在する。**分裂性を問わない**。',
  ],
  whereThisLeavesUs: {
    jLevel:
      '★`jExp` の水準は分裂性から完全に解放された: ' +
      '934 で `q` を作り、932 で `jExp p E = −v(q)` を出せる。',
    curveLevel:
      '☆ただし**曲線の水準**（V\u00e9lu の商の計算）は依然 `E ⊗ Lv = E_q`（変数変換を除いて）' +
      'を要求し、これは分裂の場合にしか成り立たない。' +
      '★したがって非分裂は捻り経由（919-929）で分裂の場合に帰着させる。',
  },
  remaining: [
    '☆(1) 「`d` を非平方単数に取れば `quadTwist` は `p` で分裂乗法還元」——' +
      '921 + 922 を並べる。要る配管は `integralModel R (quadTwist W d)` と ' +
      '`quadTwist (integralModel R W) d` の同定。',
    '☆(2) 捻った対の `hv`・`hw`（μ_l の和との同定）——分裂の場合と同じ帳簿。' +
      '925 が曲線の等式を、926 が変換則を与えている。',
    '☆(3) 組み立て: 925 → 926 → 927 → 934/932 → 929 → 903。',
  ],
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Inverse20260901 を書いた');
