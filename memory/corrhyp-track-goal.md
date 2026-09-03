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
教訓を分離して記録した。ユーザーが後に正しい形式(分子=分母、
「§1 5/5, §2 6/6, §3 3/3, §4 2/2, §5 7/7, §6 1/1達成」= 24/24 全項目 Found
達成)で `/goal` を再設定——これは自己矛盾ではない、文字通りの最終目標。

★★**2026-09-04続報、23コミット目でLemma 4.1の核心を確立**。
`FieldLimit.lean`に`Spec K = lim Spec R`(`R`は`K`の有限生成`k`-部分環、
`AffineTransitionLimit.lean`が要求する余濾過的な極限錐)をsorry無しで
証明した(`isLimit_specKCone`)。環側の余極限(`isColimitToRingCatCocone`、
`Subalgebra.iSupLift`のRingHom版がmathlibに無かったので手で構成)を
`IsColimit.op`+`Scheme.Spec`の極限保存(`Γ⊣Spec`の右随伴)でSchemeへ運んだ。
配管の教訓: `FgSubalgebra k K`が`def`(非簡約)のため、Preorder-as-category
の文脈でCategory.comp_id等をsimp/rwすると transparency エラーで止まる——
mathlib自身が使う`set_option backward.isDefEq.respectTransparency false`
が直接効いた(`tools/lean-idioms.md`に追記予定)。

★正直な現状(24/24には程遠い): numbered itemとして完全に証明が閉じたのは
`Proposition 3.2`(§3)のみ。§1のCorr/isogenyは`FuchsianGroup`モデル上で
定義として機能する形まで到達。§4はLemma 4.1の**scaffolding**(Spec K=極限)
が完成し本体まであと一歩。§2(Margulis/Shimura本体・Thm2.5-2.6)・
§5(Gauss-Bonnet、mathlibに0件)・§6(Roydenの定理)は未着手で、
それぞれ独立した大きな未構築の数学を要する。

★★★★2026-09-04さらに続報(第25-29件、Interface修正・具体的instance・
道具による集計)。原文p.5の`Comm(Γ)`定義(離散性不問)に合わせて
`Interface`を修正(`IsDiscrete`/`Gamma_isDiscrete`を新設、純追加)。
`Found/CorrHyp/Instance.lean`に`HyperbolicCurveData`の具体的な項
`corrHypInstance`(`FuchsianGroup`で構成)を作り、その上で
**`Skeleton.CorrHyp.prop_3_2`のsorryを文字通り埋める**
(`prop_3_2_at_instance`、`funext+rfl`でSkeletonの文と関数として完全一致を
確認済み)ことに成功——「関連する具体モデル」ではなく「Skeletonの主張
そのものの実装」。さらに`Found/CorrHyp/ModularExample.lean`でモジュラー群
`SL(2,ℤ)`と主合同部分群`Γ(2)`(mathlibの`discreteSpecialLinearGroupIntRangeSL`・
`CongruenceSubgroup.Gamma`)という教科書的な例を使い、§1の
`Definition 1.1/1.3/1.4/1.5`を`.src`で登記した(`Definition 1.2`は
`FEt`がProp包みのため`IsTrivial`が退化することを確認し、正直に見送った)。

`tools/corrhyp-progress.mjs`(genell-progress.mjsと同じ分子規則で新設)の
機械集計: **`CorrHyp §1 4/5, §2 0/6, §3 1/3, §4 0/2, §5 0/7, §6 0/1
—— 合計5/24`**。§2の`Proposition 2.4`/`Theorem 2.5`は
`MargulisArithmetic`/`ShimuraArithmetic`のplaceholderで「閉じる」ことは
明示的に見送った(退化した証明になるため、`corrhyp-goal.md`に歯止めを記録)。
