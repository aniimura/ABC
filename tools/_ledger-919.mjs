// 2026-09-01（第 913-919）—— Lemma 3.5 の残りを最終形に書き直す。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Remaining20260901b = {
  measuredAt: '2026-09-01',
  block: '第 913-919',
  supersedes: 'lemma35Remaining20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  splitChainComplete:
    '★★★★**分裂乗法還元の場合の連鎖は完結している**: ' +
    'stableLine_is_mu_of_coprime(906) → tateParam_quot_velu_j_dvr(914) ' +
    '→ minDeltaExp_eq_mul_of_veluMu(904) → lemma_3_5_velu_bad_delta(903)。',
  closedThisStretch: [
    '★(a) Vélu の商を Tate モデルへ移す段（915 `j_veluQuotientFull_variableChange`）。',
    '★(b) `Q ⊗ Lv` を `Φ(ζ)` と同定する段（916-917 `exists_mu_point_of_stable`）。',
    '★(c) の降下部分（919 `minDeltaExp_eq_mul_of_twist`）——' +
      '半安定なら `minDeltaExp = max(0, −v_p(j))` なので `j` だけで決まり、捻りで変わらない。',
    '☆配管（913）: `j = Δ⁻¹c₄³` に直せば `IsElliptic` の `rw` motive 問題を避けられる。',
  ],
  remaining: [
    '☆**残るのは 1 つだけ**——`j` が同じ**分裂**乗法還元の対を 1 つ作ること。' +
      '不分岐な 2 次の捻り（`c₄ ↦ d²c₄`・`c₆ ↦ d³c₆` なので `j` は不変）がそれを与える。' +
      '★要るのは (i) 捻りの構成（mathlib の ModelsWithJ が使えるか要調査）、' +
      '(ii) 捻りが `p` で分裂乗法還元をもつこと（`d` を非平方単数に取る）、' +
      '(iii) 捻りが半安定であること。',
  ],
  note:
    '★外部文献への未証明の依存は無い。ht^Falt は第 704 で内製、' +
    'Néron–Ogg–Shafarevich は第 902 で不要になった。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Remaining20260901b を書いた');
