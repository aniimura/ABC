// 2026-09-01（第 992）—— 完備化での分裂の二者択一。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35SplitAtCompletion20260901 = {
  measuredAt: '2026-09-01',
  block: '第 992',
  supersedes: 'lemma35CharTwoRoute20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`splits_or_twist_at_completion`（992）——' +
      '第 982 の二者択一を**数体の素点の完備化にそのまま当てた**。' +
      '☆第 983 では `[Fintype k]` が無くて当てられなかったが、' +
      '第 989（剰余体の有限性）を内製したので `Fintype.ofFinite` で通る。' +
      '★仮説は `W.a₁ = 0`・`W.a₃ = 0`（体の側の正規化）と乗法還元だけ。',
  ],
  chain: [
    '☆分裂性の連鎖はこれで繋がった:',
    '976（乗法還元）→ 983（剰余体で `c₄ ≠ 0`）→ 981（整モデルの `IsCharNeTwoNF`）',
    '→ 989（剰余体の有限性）→ 982（二者択一）→ **992（完備化で当てる）**',
    '→ 990（捻り `d` を `𝓞 L` に降ろす）→ 929（非分裂の降下）',
  ],
  remaining:
    '☆残るのは (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段' +
    '（`p ∤ 2`。`p ∣ 2` は第 991 のとおり不分岐 2 次拡大＋943＋944）、' +
    '(ii) 992 を `hasSplitMultiplicativeReduction_baseChange` に流す段、' +
    '(iii) 972 に並べて `minDeltaExp_variableChange` で戻す段。',
  note:
    '★944-992 の 44 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35SplitAtCompletion20260901 を書いた');
