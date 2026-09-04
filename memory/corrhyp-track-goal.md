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
機械集計(当初): **`CorrHyp §1 4/5, §2 0/6, §3 1/3, §4 0/2, §5 0/7, §6 0/1
—— 合計5/24`**。§2の`Proposition 2.4`/`Theorem 2.5`は
`MargulisArithmetic`/`ShimuraArithmetic`のplaceholderで「閉じる」ことは
明示的に見送った(退化した証明になるため、`corrhyp-goal.md`に歯止めを記録)。

★★★★★2026-09-04さらに続報(第30件、Hecke型対角共役でDefinition 2.1
を達成)。`Found/CorrHyp/HeckeExample.lean`: `g:=diag(2,1/2)∈SL(2,ℝ)`を
使い、`SL(2,ℤ)`が(古典的にarithmeticなので)自身のcommensuratorの中で
無限指数を持つことを完全にsorry無しで証明した——`Γ(4)`(mathlibに
`FiniteIndex`済み)を両方向の有限指数下界に使って`g2∈Comm(Γ_SL2Z)`を、
`g2^k`(`k≥1`)の`(1,1)`成分`2⁻ᵏ`が整数になりえないことから無限個の
相異なる剰余類の存在を、それぞれ確立した。配管の詰まり
(`Matrix.SpecialLinearGroup`のinstances-transparency問題)は
「`set A := (M:Matrix...)`で台のMatrixを先に1回だけ取り出す」という
具体的な回避策で解決——`tools/lean-idioms.md`に記録済み。

★★★★★★2026-09-04さらに続報(第31-33件、Definition 2.3を完成)。
ユーザーが「§2本体(Margulis/Shimura)に取り組む」を選択。調査の結果
`Definition 2.2`(Margulis-arithmetic、代数群の部分群スキーム分類・実点の
リー群構造がmathlibに丸ごと不在)は依然として人年規模と確認したが、
`Definition 2.3`(Shimura-arithmetic)は`F:=ℚ`(無限素点1個、「他で
非自明」が空虚に真)・`A:=M_2(ℚ)`(`Algebra.IsCentral.matrix`・
`IsSimpleRing.matrix`・`Module.finrank_matrix`で「4次元中心的単純多元環」
という四元数環の特徴づけを満たす)・`matrixEquivTensor`(唯一の無限
素点で自明)・`O_A:=M_2(ℤ)`(`Γ_SL2Z=O_A∩SL_2(ℝ)`を成分レベルで証明)
という4データすべてを`Found/CorrHyp/ShimuraArithmeticData.lean`で構成し、
`corrHypInstance`の`ShimuraArithmetic:=fun _↦False`placeholderを本物の
述語に差し替えた`corrHypInstance2`(`Instance2.lean`)へ接続、`.src`を
正当に登記した——`Proposition 3.2`が壊れないことも確認済み(安全な差し替え)。

**最終集計: `CorrHyp §1 4/5, §2 2/6, §3 1/3, §4 0/2, §5 0/7, §6 0/1
—— 合計7/24`**(このセッション開始時の1/24から7倍)。

★★★★★★★2026-09-04さらに続報(第34-36件、§4のScheme基盤を完成)。
ユーザーが「§4が使用している理論があり工数が多いという問題か」「著者が証明を
畳んでおり後追いが困難か。工数が多いことで踏みとどまる必要はない」と明示的に
指示。それまで関数体案(`F⊗_k K`が体にならず頓挫)で止まっていた`Space`/`FEt`/
`Ext`の設計を、実際にScheme案(`Space:=Over BaseK`、`BaseK:=Spec ℚ`)へ転換して
着手したところ、mathlibの`AlgebraicGeometry.Etale`/`IsFinite`
(`MorphismProperty`基盤、合成・恒等・base change安定性が`instance`で自動)が
見積もりよりずっと厚く、`FEt`一式(`idFEt`/`comp`/`pullback`/`pbFst`/`pbSnd`)が
1セッションで完成。`Ext`(係数拡大、`K:=ℝ`)は手作りの`Limits.pullback.map`では
なくmathlib自身の`Over.pullback f ⋙ Over.map f`の合成として組むと
`CategoryTheory.MorphismProperty.overPullbackMap`が base change 安定性からの
遺伝を一発で与える、という設計転換が鍵だった(`IsPullback.of_right`等の
手作りpasting は遥かに長くなり撤回)。

途中で`HyperbolicCurveData.Space:Type`(`Type 0`固定)と`Over BaseK:Type 1`
(Schemeは大きな圏なので宿命的に1universe上)が衝突する**設計上の欠陥**を
発見・修正——`Space`だけを`Type u`に universe 多相化した(`corrHypInstance`/
`corrHypInstance2`は`u=0`で無傷)。`corrHypInstance3`(`Instance3.lean`)として
具体化まで完成、`Lemma 4.1`のstatementが型検査を通ることも確認。

`Lemma 4.1`本体は未証明のまま(集計は7/24で変わらず)——
`AffineTransitionLimit.lean`の前提(`IsCofiltered`・`IsAffineHom`)は
`k:=ℚ,K:=ℝ`について`infer_instance`で自動充足することを確認したが、
「極限上の新しいスキームが有限段階のスキームのbase changeとして存在する」
という構成的降下(EGA IV 8.8-8.10相当)はmathlibに1つの補題としては
無く、組み立てが要る——`corrhyp-goal.md`§4に具体的な組み立て方の見通し
(`exists_isOpenCover_and_isAffine`→環の spreading→貼り合わせ→étale条件の
降下)を記録した。次に戻るときはここから。

コミット: `2af0967d`(SchemeFEt.lean・universe多相化)・`1626cc06`
(Instance3.lean)。[[corrhyp-scheme-universe-mismatch]] に universe の
教訓を分離して記録済み。

★★★★★★★★2026-09-04さらに続報(第37-45件、`Lemma 4.1`構成の核心部品が
完成)。同じセッション内で継続。`FieldLimit.lean`に多項式の係数降下
一式(`exists_fg_subalgebra_polynomial`/`_pair`/`_pair_monic`/`_family`)→
`StandardEtalePair`(有限エタール射の標準的表示)の降下
(`exists_fg_subalgebra_standardEtaleCond`/`_standardEtalePair`/
`_standardEtalePair_map`)→**`StandardEtalePair.Ring`のbase change**
(`Bivariate_equivMvPolynomial_map`・`standardEtalePair_ring_baseChange`・
`standardEtalePairRingBaseChange`:
`TensorProduct R K P₀.Ring ≃ₐ[K] (P₀.map (algebraMap R K)).Ring`)
まで、全9個の補題をsorry無しで積み上げた。

途中、「テンソル積は余極限を保存する」という一般圏論(`Under.pushout`が
left adjoint)を試みたが`Under.mk`のdefeq配管で頓挫し撤退——代わりに
mathlib既存の`MvPolynomial.algebraTensorAlgEquiv`・
`Algebra.TensorProduct.tensorQuotientEquiv`・`Ideal.quotientEquivAlg`を
`StandardEtalePair.Ring`(mathlibの明示的な二変数多項式商)に直接適用する、
より軽い道に乗り換えたのが鍵——一般論より具体構成のほうが速かった実例。

**これでLemma 4.1の構成的降下のうち最も技術的に重い部品(環レベルの
base change)が完成し、残るのはスキーム論の定型的な貼り合わせ
(アフィン開被覆のétale-locusでの細分→各片への適用→遷移射の一致の降下)
のみという段階に到達**。集計は7/24で変わらず(numbered itemとしての
`Lemma 4.1`自体はまだ)だが、`corrhyp-goal.md`§4に組み立て方の全体像を
記録済み。

コミット: `b45a6655`(map接続)・`8c5a010a`(Ring base change完成、
これが核心)・`57a9571c`(記録)。
