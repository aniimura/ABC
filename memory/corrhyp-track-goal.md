---
name: corrhyp-track-goal
description: CorrHyp(Correspondences on Hyperbolic Curves)のLean形式化——Skeleton(24/24)を2026-09-04に完了した。
metadata:
  type: project
---

`ResearchPaper/corrhyp-goal.md` に CorrHyp(Mochizuki, *Correspondences on Hyperbolic Curves*, 18頁)の
Skeleton を完了した: `§1 5/5, §2 6/6, §3 3/3, §4 2/2, §5 7/7, §6 1/1`(合計 24/24)。
§0(Theorem A/B/C)は本文中で後続定理の再掲と明記されているため対象外——[[genell-track-b]] の §1 起点と同じ扱い。

`Interface/CorrHyp/HyperbolicCurve.lean`(`HyperbolicCurveData`/`StackType` を posit)+
`Skeleton/CorrHyp/Section1.lean`〜`Section6.lean`。`lake build` 0 エラー、
`tools/check.mjs` で G1(出典・逐語照合)を全項目パス、G9(非空虚性の対照)14件は
プロジェクト全体の既知 debt として未着手のまま残した(ブロッキングではない)。

**Why:** G1 の逐語照合で `tools/check.mjs`(live pdftotext 呼び出し)と手打ちの
引用が何度も食い違った——`étale`(合成済み)対 `´etale`(pdftotextの分解形)、
`(g′,r′)` のプライムと `−`(U+2212)がこの1箇所だけ pdftotext で脱落、等。
`0_Source/*.txt` の静的キャッシュと check.mjs のライブ呼び出しで結果が違うことも
実測——[[gate-shell-pdftotext-differs]] と同型の罠が **同じツールの2つの実行**の
間でも起きる。

**How to apply:** 逐語引用を書いたら必ず `node tools/check.mjs`(または該当ページだけ
`pdftotext` を直接呼んで)照合すること。「.txt を読んで手で写す」だけでは信用できない。
[[genell-track-b]](ABC3b=このセッション)が GenEll と CorrHyp の両方を持つことになった。

★★**2026-09-04、Found 層(Track B)に着手・16 コミット**。`FuchsianGroup`
(`SL(2,ℝ)` の離散部分群、`ℂ` 上の解析的モデル)で §1(`Corr` の圏構造・
isogeny の同値関係)・§2(`self_le_commensurator`)・§3(`Proposition 3.2`
**完全証明**、`fg_of_fg_finiteIndex` = mathlib に無かった逆向き
Reidemeister–Schreier を自作)を sorry 無しで実装、純粋な群論で閉じる
範囲を掘り尽くした。§4 は `FieldLimit.lean` で着手(スキーム論トラック、
`FuchsianGroup` とは別建て)。残る節点は3つの独立した未構築の数学
(①双曲幾何=Gauss–Bonnet 0件、②スキーム論=`AffineTransitionLimit.lean`
に直撃する道具あり、③代数群論=非可換 Galois コホモロジー)に分岐する
——詳細と mathlib 補題名は `corrhyp-goal.md` §3。

★この過程で `/goal` の条件文が「0/N…達成」の自己矛盾形になり Stop hook が
延々と「未達成」判定を繰り返した——[[goal-condition-zero-numerator-trap]] に
教訓を分離して記録した。
