// 2026-09-01（第 981）—— 整モデルの IsCharNeTwoNF は体の側で決まる。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35IntegralNF20260901 = {
  measuredAt: '2026-09-01',
  block: '第 981',
  supersedes: 'lemma35CharNeTwoNF20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`isCharNeTwoNF_integralModel`（981）——' +
      '体の側が `a₁ = a₃ = 0` なら整モデルもそうである。' +
      '☆mathlib の `integralModel_a₁_eq`・`integralModel_a₃_eq` と' +
      '`algebraMap R K` の単射性だけ。**剰余標数の条件は要らない**。',
  ],
  measurement:
    '★★これで正規化の場所が変わった: 第 980 は `R` の中で正規化する道' +
    '（`2` が `R` で単元＝`p ∤ 2`）だったが、981 により**体 `L` の側で正規化すればよい**——' +
    '標数 0 なら `2` は必ず可逆なので剰余標数の条件は消える。' +
    '☆残る問題は「体の側で正規化した曲線が `p` で整のままか」に移った。' +
    '`p ∤ 2` なら正規化の変数変換（`u = 1`, `r = 0`, `s = −a₁/2`, `t = −a₃/2`）は' +
    '`p` で整だからそのまま通る。`p ∣ 2` ではそこが崩れうる。',
  remaining:
    '☆残るのは (i) `p ∣ 2` の手当て（不分岐 2 次拡大か別の議論）、' +
    '(ii) 979 ＋ 981 を `hasSplitMultiplicativeReduction_baseChange` に流す段、' +
    '(iii) 972 に全部並べて `minDeltaExp_variableChange` で `E`・`E′` に戻す段。',
  note:
    '★944-981 の 38 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35IntegralNF20260901 を書いた');
