// 2026-09-01（第 907-912）—— Lemma 3.5 の残りを書き直す。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Remaining20260901 = {
  measuredAt: '2026-09-01',
  block: '第 907-912',
  supersedes: 'lemma35Remaining20260831',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★極小性の移行（第 908 `isMinimal_baseChange_of_c4`）——' +
      '一般形は弱近似を要するが、Lemma 3.5 が使う乗法還元の場合は ' +
      '`v_p(c₄) = 0` なので `isMinimal_of_c4_vAdd_eq_zero` が直接与える。',
    '★乗法還元の移行（第 909 `hasMultiplicativeReduction_baseChange`）——' +
      '`v(Δ) < 1` と `v(c₄) = 1` を完備化の付値に移すだけ。',
    '★分裂乗法還元（第 910）——`hsplit : True` を mathlib の構造体のフィールド ' +
      '（剰余体上の 2 次式の `Splits`）に置き換えて証明した。',
    '★大域→局所の Vélu（第 911 `veluQuotientFull_baseChange`）——' +
      '`E′ = E/⟨Q⟩` を底変換すると `E′ ⊗ A = (E ⊗ A)/⟨Q ⊗ A⟩`。',
    '★`hB`・`hBx`（第 912）——反転で安定な点集合なら自動的に成り立つ。',
  ],
  remaining: [
    '☆(a) `E ⊗ Lv` の Vélu の商を Tate モデルへ移す段: ' +
      '`tateParamR_spec` の変数変換 `D` を `veluQuotientFull_variableChange`（証明済み、' +
      '仮説 hB・hBx は第 912 が与える）で通す。',
    '☆(b) `Q ⊗ Lv` を `Φ(ζ)` と同定する段: `Φ` は加法同型（全射）なので ' +
      '`Q ⊗ Lv = Φ([x])` と書け、位数 l から `x^l = Q_E^k`、' +
      '`stableLine_is_mu_of_coprime`（第 906、証明済み）で `[x] = [ζ]`。',
    '☆(c) 非分裂乗法還元の降下（不分岐 2 次拡大へ上げて局所高さが変わらないことを使う）。' +
      '★原文は semi-stable としか言わないので必要。',
  ],
  chainNow:
    'stableLine_is_mu_of_coprime（906）→ tateParam_quot_velu_dvr（900）' +
    '→ minDeltaExp_eq_mul_of_veluMu（904）→ lemma_3_5_velu_bad_delta（903）' +
    '——★局所の仮説から Lemma 3.5 の結論まで一本で通っている。',
  note:
    '★外部文献への未証明の依存は無い。ht^Falt は第 704 で内製、' +
    'Néron–Ogg–Shafarevich は第 902 で不要になった。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Remaining20260901 を書いた');
