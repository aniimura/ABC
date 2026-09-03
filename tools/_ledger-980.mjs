// 2026-09-01（第 980）—— 整モデルを IsCharNeTwoNF に正規化する。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35CharNeTwoNF20260901 = {
  measuredAt: '2026-09-01',
  block: '第 980',
  supersedes: 'lemma35SplitCharTwo20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`exists_variableChange_isCharNeTwoNF_of_isUnit_two`（980）——' +
      '`2` が単元なら変数変換 1 回で `a₁ = a₃ = 0` にできる。' +
      '☆mathlib の `exists_variableChange_isCharNeTwoNF`（`[Invertible 2]`）に' +
      '`IsUnit (2 : R)` から作った `Invertible` を渡すだけ。' +
      '★完備化の整数環で `2` が単元になるのは `p ∤ 2` のとき（第 953 に `n = 2`）。',
  ],
  remaining:
    '☆残るのは (i) 979 ＋ 980 を `hasSplitMultiplicativeReduction_baseChange` に流す段' +
    '（`p ∣ 2` は第 979 の測定どおり別の手当てが要る）、' +
    '(ii) 972 に全部並べて `minDeltaExp_variableChange` で `E`・`E′` に戻す段。',
  suppliedSummary: [
    '★第 972 の仮説の供給表（2026-09-01 現在）:',
    '`C`・`hC`・`hc4ne`・`hc4` → 954',
    '極小性・乗法還元 → 976（908＋909＋973）',
    '`hp` → 964、`hlu` → 953、`hcop` → 978',
    '`h2`・`h2K`・`hql`・`hΔ` → 977',
    '`hvw` → 961＋962',
    '`P`・`hW′j`・`hellQ` → 972 が内製（970＋947＋971）',
    '分裂性 → 979＋980（`p ∤ 2`）／`p ∣ 2` は未',
  ],
  note:
    '★944-980 の 37 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35CharNeTwoNF20260901 を書いた');
