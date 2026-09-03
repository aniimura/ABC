// 2026-09-01（第 993）—— 悪い素点での分裂／捻りの二者択一。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SplitOrTwistBad20260901 = {
  measuredAt: '2026-09-01',
  block: '第 993',
  supersedes: 'lemma35SplitAtCompletion20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`hasSplit_or_twist_at_bad_prime`（993）——' +
      '悪い素点では、完備化で**分裂乗法還元をもつ**か、' +
      '**ある捻りで 2 次式が分裂する**かのどちらかである。' +
      '☆992（完備化での二者択一）を ' +
      '`hasSplitMultiplicativeReduction_baseChange` に流しただけ。',
  ],
  chainComplete: [
    '★分裂性の連鎖が終点まで繋がった:',
    '976（乗法還元）→ 983（剰余体で `c₄ ≠ 0`）→ 981（整モデルの `IsCharNeTwoNF`）',
    '→ 989（剰余体の有限性・内製）→ 982（二者択一）→ 992（完備化で当てる）',
    '→ **993（`HasSplitMultiplicativeReduction` に流す）**',
    '→ 990（捻り `d` を `𝓞 L` に降ろす）→ 929（`Δ_min` を降ろす）',
  ],
  remaining:
    '☆残るのは (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段' +
    '（`p ∤ 2`。`p ∣ 2` は第 991 のとおり不分岐 2 次拡大＋943＋944）、' +
    '(ii) 993 の 2 つの枝を 972 に流し、`minDeltaExp_variableChange` で `E`・`E′` に戻す段。',
  note:
    '★944-993 の 45 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SplitOrTwistBad20260901 を書いた');
