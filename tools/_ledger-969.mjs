// 2026-09-01（第 969）—— 商の楕円性は E′ から来る（循環を断つ）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluEllipticFromEprime20260901 = {
  measuredAt: '2026-09-01',
  block: '第 969',
  supersedes: 'lemma35PointImageEq20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`isElliptic_veluQuotient_vcPoint`（969）——' +
      'Vélu の商の楕円性は `E′` から来る。',
  ],
  whyItMatters:
    '★★第 967 は `veluQuotientFull (C • (E ⊗ K)) (…)` の楕円性を**インスタンスとして受けて**いた。' +
    'その出どころを第 962（`Δ = l¹²Δ`）に求めると `ζ` が要り、`ζ` は点から来るので**循環する**。' +
    '☆循環しない道: `veluQuotientFull_variableChange` により' +
    '`veluQuotientFull (C • (E ⊗ K)) (捻った点集合) = C • veluQuotientFull (E ⊗ K) (点集合)`、' +
    '右辺は `veluQuotientFull_baseChange` により `C • (E′ ⊗ K)`。' +
    '`E′ ⊗ K` は楕円だから変数変換しても楕円である。' +
    '★仮説 `hB`・`hBx` は点集合が反転で安定（949）なので自動的に満たされる。',
  remainingLayers: [
    '☆残るのは (1)-(9) を 1 つの証明として書き下す作業のみ。' +
      '楕円性インスタンスの出どころ（969）が付いたので、循環はもう無い。',
  ],
  note:
    '★944-969 の 26 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluEllipticFromEprime20260901 を書いた');
