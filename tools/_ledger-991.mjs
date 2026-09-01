// 2026-09-01（第 991）—— 剰余標数 2 の非分裂は 943＋944 で扱う（測定と訂正）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35CharTwoRoute20260901 = {
  measuredAt: '2026-09-01',
  block: '第 991（測定）',
  supersedes: 'lemma35ResidueRep20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  measurement:
    '★★★★**剰余標数 2 の非分裂は捻りでは直らない**。' +
    '☆有限体の標数 2 では `X² + aX + b`（`a ≠ 0`）は分裂しないことがある' +
    '（例: 𝔽₂ 上の `X² + X + 1`）。' +
    '一方その標数では非平方元が存在しないので、捻りという手が使えない。' +
    '★さらに整モデルを `IsCharNeTwoNF`（`a₁ = a₃ = 0`）に正規化するには `2` の可逆性が要る' +
    '（第 980・981）ので、`p ∣ 2` ではそこも通らない。',
  theFix:
    '★★★★**その場合は不分岐 2 次拡大に上がるのが正しい**——剰余体が `𝔽_q` から `𝔽_{q²}` に' +
    '広がるので `X² + X + 1` の類は分裂する。' +
    '☆そしてその道具は**本セッションで既に建ててある**:' +
    ' `isAdicComplete_valuationSubring`（第 943、任意の完備離散付値体の整数環）と' +
    ' `tateParamR_map` / `tateParam_descend`（第 944、Tate 母数は体の拡大で変わらない）。' +
    '★第 945 で「拡大の層は不要」と書いたのは **`μ_l`・`ζ` の層について**であり、' +
    '**分裂性の層については `p ∣ 2` で拡大が要る**（訂正）。',
  routes: [
    '☆したがって分裂性の手当ては素点で 2 通りに分かれる:',
    '(A) `p ∤ 2`: 981（体の側で `a₁ = a₃ = 0` に正規化）＋ 982（二者択一）＋ ' +
      '990（捻りを `L` に降ろす）＋ 929（非分裂の降下）。',
    '(B) `p ∣ 2`: 不分岐 2 次拡大 `Lv′` に上がって分裂させ、' +
      '943（`IsAdicComplete`）で台を作り、944（`tateParam_descend`）で `Lv` に降ろす。',
  ],
  remaining:
    '☆残るのは (A) の配管、(B) の不分岐 2 次拡大の構成、' +
    'および 972 に並べて `minDeltaExp_variableChange` で戻す段。' +
    '★mathlib の穴はもう無い（第 989 で最後の 1 つを埋めた）。',
  note:
    '★944-990 の 43 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35CharTwoRoute20260901 を書いた');
