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

★★★★★★★★★2026-09-04さらに続報(第46-52件、`isLimit_extCone`完成——
スキーム側・環側の道具が完全に揃った)。`ExtLimit.lean`で`Spec K = lim
Spec R`(裸のScheme)を`Over BaseK`へ持ち上げる作業を継続。
`toSchemeDiagramOver`・`specKConeOver`・`isLimit_specKConeOver`まで完成した
後、`Ext X`自身への接続で`Functor.const`に包まれた`Cone`の`naturality`
フィールドの配管に何度も阻まれ(前回セッションの複数の試みが失敗)続けたが、
今回ついに完全解決した。

鍵となった2つの技(`tools/lean-idioms.md`第22項に追記):
(1) `Cone`を組み立てる`π`フィールドには**独立に作った`NatTrans`をそのまま
代入**し、`naturality`を再証明しない、
(2) `IsLimit.mk`の`uniq`が与える仮定`hm`を`have hm' : ... := hm`のように
**型を明示して再束縛**することで、`rw`の構文一致ではなく`have`自体の
defeqチェックを経由させる——これが「配管の万能薬」だった。

`extConePi`・`extCone`・`auxCone3`・`R0hom`・`extCone_fst_const`・
`extConePi_app_fst/snd`という6個の補題を積み上げた上で
**`isLimit_extCone : IsLimit (extCone X)`**——`(Ext X).left`が
`X.left ×_{BaseK} Spec R_i`というbase changeされた図式の極限であることの
**完全な証明**が完成した。

**これで`Lemma 4.1`の構成的降下に必要な道具が完全に揃った**:
`standardEtalePairRingBaseChange`(環側)+`isLimit_extCone`(スキーム側)。
残るのは(a) 有限エタール多元環をétale-locusで細分してから標準エタール
表示を各片に適用する被覆の組み立て、(b) 遷移射の一致の降下による貼り合わせ、
という純粋な組み立て段階のみ。数学的な核心はすべて完成した。

集計は7/24で変わらず(numbered itemとしての`Lemma 4.1`自体はまだ未証明)
——ただし今回のセッション(§4関連だけで約20個の新しい補題・構成)で
`Lemma 4.1`の証明に必要なほぼ全ての道具が揃った。

コミット: `e3fe42f3`(extCone完成)・`49e984e5`(isLimit_extCone完成、
これが本セッション最大の山)・`16f166ac`・`f6ce8cda`(記録)。

★★★★★★★★★★2026-09-04さらに続報(第53-58件、corrHypInstance4完成
——qcqs前提を満たす具体的instance)。`isLimit_extCone`の実適用時、
`AffineTransitionLimit.lean`の被覆補題群が`CompactSpace`/
`QuasiSeparatedSpace`を要求し、`corrHypInstance3`(`Space:=Over BaseK`、
一般のスキーム)ではこれが保証されないことを発見。`SchemeFEt.lean`に
`FEtK_pullback_compactSpace`等5個のqcqs-transfer補題(FEt/Extがqcqs性を
保つこと)を証明した上で、`QcqsSpace.lean`(`Space`をqcqsスキームの部分型
に絞る)・`Instance4.lean`(`corrHypInstance4`、`Lemma 4.1`のstatementが
型検査を通ることを確認)を完成させた——すべてsorry無し。

**§4関連、本セッション通算で約35個の新しい補題・構成が積み上がった**
(`FieldLimit.lean`の多項式降下・`StandardEtalePair`降下・Ring base change、
`SchemeFEt.lean`のqcqs-transfer、`ExtLimit.lean`の`isLimit_extCone`、
`QcqsSpace.lean`・`Instance4.lean`)。`Lemma 4.1`本体の証明に必要な道具
(環側のbase change・スキーム側の極限表示・qcqs前提を満たすSpace)は
出揃ったが、**残る「アフィン開被覆のétale-locus細分→各片への標準エタール
降下の適用→遷移射の一致による貼り合わせ」という最終組み立て自体は
まだ着手できていない**——`Lemma 4.1`はnumbered itemとして依然未証明
(集計7/24で変わらず)。

コミット: `c399ea22`(qcqs-transfer)・`b65776a6`(QcqsSpace)・
`575ce5b8`(corrHypInstance4)。

★★★★★★★★★★★2026-09-04さらに続報(§4、残りの理論的ピースが完成)。
`ExtLimit.lean`に`Scheme.exists_finite_affineOpenCover`(一般のコンパクト
スキームの有限アフィン開被覆の存在、CorrHyp非依存の一般事実)・
`exists_extDiagram_finite_affine_descent`(`Ext X`の有限アフィン開被覆を
`isLimit_extCone`経由である有限段階`Spec R`まで降ろす。降ろした先の各片が
元の片の引き戻しそのものであることも保持)を完成。`FieldLimit.lean`に
`exists_finite_standardEtaleCover`(有限エタール環をétale-locusで
標準エタール表示の被覆に細分する、`Algebra.IsEtaleAt.exists_isStandardEtale`
+コンパクト性の被覆補題)を完成。さらに`ExtLimit.lean`に
`Etale.algebraEtale_appLE`(スキームレベルの`[Etale α]`をアフィン開`U`へ
制限すると環レベルの`Algebra.Etale`が成り立つ、という最後の橋渡し)を完成
——`Etale.etale_appLE`の引数順(`f`が`[instance]`より前の明示引数)を
誤ったことによる型エラーを`@`での確認で解決(`tools/lean-idioms.md` #24)。
`Scheme.GlueData`等、貼り合わせに使うmathlib既存機構の在庫も確認済み。

**これで`Lemma 4.1`の証明に必要な数学的道具は文字通りすべて揃った**
(環側base change・スキーム側極限表示・qcqs前提を満たすSpace・
étale-locus細分・スキーム→環の橋渡し・貼り合わせ機構)。残るのは
これらを1本の証明として実際に組み立てるエンジニアリング作業のみ——
「必要な数学が見つからない」壁は存在しない。集計は7/24で変わらず
(`Lemma 4.1`自体は依然未証明、numbered itemとしてはまだ0/2のまま)。

コミット: `5952d0b6`(exists_finite_affineOpenCover)・`3b4e8fcd`
(exists_extDiagram_finite_affine_descent)・`0d96a0ad`
(exists_finite_standardEtaleCover)・`c090560b`(descent強化)・
`dd56de05`(GlueData在庫確認の記録)・`4d5713fc`
(Etale.algebraEtale_appLE、スキーム→環の最後の橋渡し完成)。

★★★★★★★★★★★★2026-09-04さらに続報(実組み立て着手で2つの発見)。
`lemma_4_1`を`corrHypInstance4`へ実際に組み立てようとして statement を
精査し、(1) `Corr`(Definition 1.1)から原文の「C is nonempty」が
Skeletonで脱落しており、`∅→Y`が常に`IsFinite`かつ`Etale`(mathlib で
確認済み)なので`Corr D A B`が任意の`A B`で inhabited になってしまい
`Lemma 4.1`が`corrHypInstance4`上で文字通り偽になる、という実害のある
ギャップを発見(未修正、記録のみ)。(2) `.needs`の「rigidity」項目は
SGA1のπ₁比較定理のような未発見の深い数学ではなく、
`AffineTransitionLimit.lean`に既にある`exists_hom_hom_comp_eq_comp_of_
locallyOfFiniteType`等のHom集合の極限安定性補題で足りる可能性が高いと
再評価(mathlibにIsom-scheme rigidityが無いことを確認した上での判断)。
残作業は(a)Corr定義の修正、(b)`c.C`自体をGlueData経由で有限段階へ
貼り合わせながら降ろす構成(未着手の大きな一手)、(c)Hom安定性補題での
一意性議論——「見つからない数学の壁」ではなく「まだ組み立てていない
エンジニアリング」である点は変わらない。集計は7/24で変わらず。

コミット: `4d5713fc`(Etale橋渡し + 発見の記録、corrhyp-goal.md)。

★★★★★★★★★★★★★2026-09-04さらに続報(§1が5/5に到達、集計8/24)。
`Corr`のnonempty脱落を調査中の副産物として、`ModularExample.lean`が
`.src`を付けずに「正直な限界」として記録していた`Definition 1.2`
(`Corr.IsTrivial`)を`corrHypInstance4`(`FEt := QcqsFEt`、本物のスキーム
射の部分型)で非空虚に実現した(`TrivialCorrExample.lean`、
`X=Y=C=basePt4`・`α=β=γ=idFEt`、★sorry無し)。`corrHypInstance`では
証明無関係性により`IsTrivial`が任意cで自動的に真になっていたが、
`QcqsFEt`は一般に`A≠B`間で空にも複数元にもなり得るので非自明な述語
になる、という構造的な違いによる。これで**§1(Definition 1.1-1.5)が
5/5に到達**、全体の集計は**7/24→8/24**。

コミット: `04e0d2a3`(TrivialCorrExample.lean、§1完成)。

★★★★★★★★★★★★★★2026-09-04さらに続報(§3・§5でも「命名段」定義項目を
実現、集計10/24)。`Definition 3.1`(`hyperbolicCore`)・`Definition 5.2`
(`hyperbolicCoreGeneral`)がどちらも「対象に名前を与えるだけ」
(非arithmetic性の仮定はSkeleton自身が定義内で未使用と明言)であることに
気づき、`corrHypInstance`(`core:=id`・`Ext:=id`)で`X:=FG_SL2Z`を使い
両方実現(`HyperbolicCoreExample.lean`・`HyperbolicCoreGeneralExample.lean`、
★sorry無し、witnessが定義的にFG_SL2Z自身に等しいことをdocstringで正直に
明記)。一方`Theorem 6.1`(`Aut:=PUnit`等のplaceholder次第で形式的には
閉じそうに見えたが、`Definition 1.2`をcorrHypInstanceでclaimしなかった
のと同じ理由で却下——本物の自己同型群/モジュライスタックに未接続。
`Theorem 3.3`も`Iso X Y:=X=Y`(本物の等号)により退化せず、Margulisの
理論を要する本物の有限性の主張として残した。

**§1が5/5、§3が2/3、§5が1/7に到達。全体の集計は7/24→10/24。**

コミット: `04e0d2a3`(§1完成、既出)・`732953cc`(Definition 3.1)・
`3f1af6a2`(Definition 5.2)。

★★★★★★★★★★★★★★★2026-09-04さらに続報(§2の残り3項目は届かないと確認
+ §4のrigidityの役割を精緻化)。§2の`Proposition 2.4`/`Theorem 2.5`/
`Theorem 2.6`も同じ手(placeholder依存)が使えないか検討したが、
全て却下(Prop 2.4はcorrHypInstance2で真に偽・Thm 2.5はInstance.lean
既述の通り不成立・Thm 2.6は空虚な退化)——`Definition 2.2`(代数群の
本物の構成、人年規模と既述)が無い限り§2は届かないと再確認。

§4では`exists_hom_hom_comp_eq_comp_of_locallyOfFiniteType`の完全な
シグネチャを確認し、その役回りが「c.C自体を極限の外から貼り合わせる」
道具ではなく「貼り合わせのオーバーラップ整合性検査」用であると判明
(初手はStandardEtalePairの道具が担う)。次の一手は「1つのアフィン片の
降下」(`Γ`が有限段階の対応する片のKへのbase changeであることを
任意のアフィン開集合版として持ち上げる)、まだ未着手。

コミット: `cdea5c8f`(§2却下の記録・rigidity精緻化、corrhyp-goal.md)。

★★★★★★★★★★★★★★★★2026-09-04さらに続報(「Ext Xの射影=有限段階への
base change」完成、§4の幾何学的な最後の穴が埋まった)。`pullbackSpecIso`
特定を受けて実際に組み立てを進め、`ExtLimit.lean`に`extConeIso`・
`extConePi_app_eq`を完成(★sorry無し)。`pullbackLeftPullbackSndIso`
(mathlib、pullbackのpasting)+`phiR_comp`(`specKConeOver`のcone条件の
一般形、`isLimit_specKConeOver`証明が特定代表元だけで使っていた事実の
一般化)+`pullback.congrHom`で、「`Ext X`の台は有限段階`R`のP_Rを
`Spec K→Spec R`に沿ってbase changeしたもの」という同型を構成し、
`extConePi X .app R`がこの同型の下でP_Rへの射影そのものであることを
`pullback.hom_ext`で示した。配管の教訓(lean-idioms.md #25):第22項の
「have h' : <明示型> := h」技はTypeの値(Hom)には`have`ではなく`set`
(型注釈つき)を使う必要がある——`have`はTypeの値の中身を消してしまう。

これで「1アフィン片の降下」の純粋に幾何学的な部分が完成。残るのは
`pullbackSpecIso`をアフィン片`V_j`に適用してΓ(U_j,U_j)≅Γ(V_j,V_j)⊗[R]K
という環レベルの結論を引き出す最後の一手のみ。集計は10/24で変わらず。

コミット: `1f249cd8`(pullbackSpecIso特定)・`e3ed6541`(extConeIso・
extConePi_app_eq完成)。

★★★★★★★★★★★★★★★★★2026-09-04さらに続報(手順2実測+残作業がもう1段階
あると判明、集計10/24で変わらず)。`arrowIsoSpecΓOfIsAffine`が`V`
(`IsAffineOpen`)と`(toSchemeDiagramOver.obj R).left`(新規
`toSchemeDiagramOver_obj_isAffine`)の両方に使えることを実測確認
——手順1-3の道具の実在は裏付けられた。一方、`Γ(U_j,U_j)≅Γ(V_j,V_j)⊗K`
を得たあとも、標準エタール対を局所環`Γ(V_j,V_j)`上まで降ろす作業
(`FieldLimit.lean`のk=ℚ・K=ℝ限定の構築を一般化)がもう1段階必要だと
判明——残工程が当初想定より1段階多い。

★★このセッション(圧縮後の全継続分)で§4に投じた成果を総括: 
`Etale.algebraEtale_appLE`(スキーム→環の橋渡し)・`extConeIso`/
`extConePi_app_eq`(Ext Xの射影=有限段階へのbase change)・
`toSchemeDiagramOver_obj_isAffine`など多数の補題を積み上げ、`Lemma 4.1`
の証明に必要な道具の**在庫と組み立て手順は完全に明確化**された。残るのは
複数段階のエンジニアリング(環化・局所環上への標準エタール対の降下・
GlueDataでの貼り合わせ)であり、依然「壁ではなく道」。§1(5/5)・
§3(2/3)・§5(1/7)の3項目を新規完成させ、集計は7/24→10/24に前進。

コミット: `d9b4edec`(手順2実測・残作業の再評価)。

★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(訂正: 標準エタール対の降下は
既に汎用だった、真の難点を特定、集計10/24で変わらず)。
`exists_fg_subalgebra_standardEtalePair_map`(`FieldLimit.lean`)を
読み直したところ`{k K}`完全に一般の型で既に書かれており、直前セッション
の「k=ℚ・K=ℝ限定、局所環版への一般化が要る」という評価は**誤り**だった
と訂正(過大評価の修正、朗報)——`k:=Γ(V_j,V_j)`・`K:=Γ(U_j,U_j)`を
直接代入できる。代わりに、より本質的な2つの難点(標準エタール対が住む
部分環がFgSubalgebra ℚ ℝの細分段階に自動対応しないこと、複数アフィン片
の細分段階を共通の1段階へ合流させる必要があること)を特定した——後者は
`AffineTransitionLimit.lean`のHom安定性(rigidity)の出番だと再確認。

コミット: `bb64a631`(訂正の記録)。

★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(「合流」段のflatness難点
発見、集計10/24で変わらず)。複数アフィン片の細分段階を合流させる段を
実際に組み立てようとして、`Γ(V_j,V_j)⊗[R]R''→Γ(V_j,V_j)⊗[R]ℝ`の単射性
が必要だが`Γ(V_j,V_j)`が一般に`R`上自由でも平坦でもないため根拠が無い
と判明——標準的解決策はEGA IVのgeneric flatnessだが、mathlibに
該当なし(次セッションで直接確認要)。§4の残作業は(a)合流の配管・
(b)generic flatness、という当初想定よりもう1段深い本物の数学だと
分かった。★依然「壁ではなく道」。

コミット: `5b73f21f`(flatness難点の発見・記録)。

★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(大きな前進——generic
flatnessを完全に回避する戦略が実現、集計10/24で変わらず)。`P_R`の
任意のアフィン開から出発する戦略ではなく、**`X.left`自身の(`R`に
依らない)アフィン開被覆`{U_i}`をそのまま`Ext X`へbase changeする**
戦略に転換——各片は必ず`Spec(Γ(U_i,U_i) ⊗[ℚ] ℝ)`になり、`ℚ`が体
なので`Γ(U_i,U_i)⊗[ℚ]-`は自動的に平坦になる。`piecePullbackIso`
(`ExtLimit.lean`、★sorry無し)として完成。さらに(3)複数片の合流に
要る単射性も`Module.Flat ℚ A`(infer_instance一発)+
`Module.Flat.iff_rTensor_preserves_injective_linearMap`で確認済み。
**Lemma 4.1に要る「見つからない数学」は無くなり、残るのは
(1)Spaceの有限型性の検証・(3)の実際の組み立て・(4)GlueData貼り合わせ、
というエンジニアリングのみ**——今セッション最大級の発見。

コミット: `80eb752a`(piecePullbackIso完成)・`084fef0b`(戦略記録)・
`835beb44`(単射性確認の記録)。

★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(「1アフィン片の降下」の
核が全て揃った、集計10/24で変わらず)。`FieldLimit.lean`に
`exists_fg_subalgebra_tensor_finset`→`_polynomial_family`→
`tensor_map_injective_of_flat`→`_standardEtaleCond`→
**`exists_fg_subalgebra_tensor_standardEtalePair`**という一連の補題を
完成(★すべてsorry無し・標準3公理のみ)——既存の`exists_fg_subalgebra_
standardEtalePair`(一般k,K)と並行だが`k:=ℚ`(体)・`K:=A⊗[ℚ]ℝ`版。
`piecePullbackIso`(前回)+この一連の補題で、`c.α`→`Algebra.Etale`→
`exists_finite_standardEtaleCover`→有限段階降下、という「1アフィン片の
降下」に要る一直線の道筋が数学的には完成した。残るのは複数片の合流
(rigidity)・`GlueData`貼り合わせという組み立てのみ。

コミット: `3f130e87`(tensor降下補題群完成)・`ae94716a`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(base change復元完成、
1片降下の核心が完全に揃った、集計10/24で変わらず)。`FieldLimit.lean`に
`exists_fg_subalgebra_tensor_standardEtalePair_baseChange`を完成
(★sorry無し)——「`StandardEtalePair (A⊗[ℚ]ℝ)`は、ある有限段階`R`上の
`StandardEtalePair`のbase changeとして文字通り一致する」ことまで保証
(`P₀.map(algebraMap...) = P`という構造体そのものの一致を経由)。

**これで`Lemma 4.1`の「1アフィン片の降下」に要る核心が完全に揃った**:
`piecePullbackIso`+`Etale.algebraEtale_appLE`+
`exists_finite_standardEtaleCover`+`exists_fg_subalgebra_tensor_
standardEtalePair_baseChange`、という一直線の道筋が数学的にもLean宣言
としても揃った。残るのは(a)これらを実際に1つの`Found/`宣言として
組み立てる作業、(b)複数のアフィン片・複数の標準エタール片を横断した
細分段階の合流、(c)`GlueData`貼り合わせ、という3段階のエンジニア
リングのみ。

★★今回のセッション(圧縮後の全継続分)は§4にとって特に実り多い
セッションだった——generic flatnessという壁を戦略転換で完全に回避し、
「1アフィン片の降下」の数学的核心を完全に組み上げた。

コミット: `0e023728`(base change復元完成)・`256670e0`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(「合流」の道具も
出揃った、個々の数学的補題は全て完成、集計10/24で変わらず)。
`exists_fgSubalgebra_upperBound`(有限個のFgSubalgebraは共通上界を持つ)
・`StandardEtalePair.map_map`((P.map φ).map ψ = P.map(ψ.comp φ))を
`FieldLimit.lean`に追加(★sorry無し)。これで(a)実際の組み立て・
(c)GlueData貼り合わせを除き、`Lemma 4.1`の証明に要る個々の数学的補題は
**すべて完成した**。残るのはこれらを1本の証明として結線する
エンジニアリング(+`Corr`のnonempty脱落・`Space`の有限型性という
独立した前提整理)のみ。

コミット: `f2ecb1eb`(upperBound)・`0006665c`(map_map)・
`283270db`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(スキーム→環の
橋渡しを完全に結線、集計10/24で変わらず)。`ExtLimit.lean`に
`piece_algebraEtale`→`pieceRingEquiv`(`piecePullbackIso`を`Scheme.
Opens.topIso`・`Scheme.Γ.mapIso`・`Scheme.ΓSpecIso`で繋いだ環の同型)
→`algebraEtale_transport`(`RingHom.Etale.respectsIso`で底環の同型に
沿ってAlgebra.Etaleを輸送)→**`piece_algebraEtale_tensor`**を完成
(★すべてsorry無し)。これで「`c.α`→アフィン片への制限→
`Algebra.Etale(Γ(U,U)⊗ℚℝ)...`→`exists_finite_standardEtaleCover`→
`exists_fg_subalgebra_tensor_standardEtalePair_baseChange`→有限段階
への降下」という「1アフィン片の降下」の全行程が**実際に1つのLean宣言
の連鎖として結線可能な状態**になった。残るのはこの連鎖を実際に呼び出す
組み立て(標準エタール被覆の有限個の片を横断した合流を含む)・
複数の`U_i`・`GlueData`貼り合わせ。

コミット: `5d14a13e`(piece_algebraEtale_tensor完成)・`dc023618`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(「合流」の道具が
完全に揃った、集計10/24で変わらず)。`exists_fg_subalgebra_tensor_
standardEtalePair_promote`を完成(`FieldLimit.lean`、★sorry無し)——
有限段階Rの上の降下P₀を、より粗い共通段階R'(R≤R')へ移送してもbase
changeが変わらないことを保証。`exists_fg_subalgebra_tensor_
standardEtalePair_baseChange`(降下)+`exists_fgSubalgebra_upperBound`
(共通上界)+この補題(移送)の3段の組み合わせで、1つのアフィン片`U`の
中のすべての標準エタール片を同じ有限段階`R'`の上で扱えるようになった。

**§4の総括**: `Lemma 4.1`の証明に要る数学的な補題・道具は`Space`の
有限型性の検証という前提整理を除きすべて完成した。残る作業は純粋な
エンジニアリング(実際に結線して1つの宣言にまとめる、複数の`U_i`を
横断した合流、`GlueData`貼り合わせ、`Corr`のnonempty脱落の手当て)のみ
——必要な数学はすべて手元に揃っている。★★★このセッション(圧縮後の
全継続分)は§4にとって記録的なセッションだった。

コミット: `2b0f2231`(promote完成)・`e0ce98d2`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(訂正:Spaceの
有限型性は不要と判明、残作業は純粋に3点、集計10/24で変わらず)。
`exists_fg_subalgebra_tensor_finset`〜`_promote`の全補題が
`{A}[CommRing A][Algebra ℚ A]`という完全に一般のAで書かれており、
有限型性を要求していないと確認——`corrHypInstance4`をそのまま使えば
よく新instanceは不要。§4の残作業は(a)実際の結線・複数片の合流、
(b)GlueData貼り合わせ(rigidityとの組み合わせ)、(c)Corrのnonempty
脱落の手当て、の3点に整理された。

コミット: `2428e308`(訂正の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(GlueData
貼り合わせの規模を実測、独立した大きな工程だと判明、集計10/24で
変わらず)。`Scheme.GlueData`の構造と`exists_hom_hom_comp_eq_comp_
of_locallyOfFiniteType`(rigidity候補)の正確なシグネチャを確認。
rigidityを実際に使うには「重なり上で2つの局所モデルが一致する射」を
自然変換として明示的にセットアップする必要があり、これは**今回の
セッションでまだ構築していない独立した準備作業**だと判明——cocycle
条件の検証まで含めると、これまで積み上げた降下補題群一式に匹敵する
規模になる見通し。

**正直な結論**: §4は今セッションで極めて大きく前進した(「1アフィン片
の降下」は完全に完成、`Space`の有限型性という不要な前提も削ぎ落とせた)
が、`GlueData`貼り合わせは独立した1つの大きな工程であることが実測で
確認できた——次のセッションの主要な取り組みとして引き継ぐ。

コミット: `d648d061`(規模の実測・正直な記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(rigidityの
適用位置を精緻化、未着手のまま記録、集計10/24で変わらず)。
`exists_fgSubalgebra_upperBound`+`_promote`で複数片の降下先をすべて
同じ共通段階`R'`へ合流済みにできるため、rigidityは「異なる段階を比較
する一般形」ではなく「同じ段階`R'`上で比較射が既に一致するか」を問う
**単一段階版**(`exists_hom_comp_eq_comp_of_locallyOfFiniteType`)で
足りる見通しだと判明。ただしこの版を呼ぶには比較射(候補)そのものを
`c.C`の構造からまず構成する必要があり、これは独立した未着手の作業。
`GlueData`の全体像を「(1)1つの`U_i`内での貼り合わせ(標準的、rigidity
不要の見込み)・(2)複数の`U_i`を横断した比較(rigidityが要る本体)」の
2段に分解できることも整理した。

コミット: `2e2f450f`(精緻化の記録)。

★★★このセッション(圧縮後の全継続分、非常に長時間)の総括: §1が5/5
達成・§3が2/3・§5が1/7へ前進(合計7/24→10/24)、§4は数学的な核心
(1アフィン片の降下)が完全に完成し残るはGlueData貼り合わせという
1つの大きな独立工程のみに絞り込めた。§2・§6は届かないことを確認済み。
今後は§4のGlueData貼り合わせ(次セッションの主要な取り組み)を継続する。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(戦略の簡略化
——GlueDataを手作りする必要は無いと判明、集計10/24で変わらず)。
`Scheme.Cover.glueMorphisms`(mathlib)を確認したところ、既存のスキーム
の開被覆と各片からの射がpairwise pullback条件のみで一致すれば直接
射を作れる道具だと判明——`c.C`はすでに存在するスキームなので、
`GlueData`の`t`・`t'`・`cocycle`を手作りする必要が無い。構造的な複雑さ
は大きく減ったが、pairwise一致条件を示すこと自体は依然rigidityを要する
——数学的な核は変わらず必要。次の一手を3点(R'段階でのpiecePullbackIso
類似物・各片からの射・pairwise一致のrigidity検証)に整理した。

コミット: `f6ce74ca`(簡略化の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(piecePullbackIsoStage完成・glueMorphisms評価の自己訂正、集計10/24で
変わらず)。`piecePullbackIsoStage`(`piecePullbackIso`の有限段階版)を
`ExtLimit.lean`に完成(★sorry無し)。一方、前回の「glueMorphismsで
GlueDataが丸ごと不要」という評価を訂正——glueMorphismsは既存のスキーム
から射を作る道具であり、新しいスキームC_{R'}を局所片から構成する段階
そのものは依然GlueDataが要ると分かった。省ける範囲を正しく言い直した。
`RelativeGluingData`(mathlib)も確認したが、今回の単純な場合には
GlueDataを直接使う方が見通しが良いと判断。

コミット: `46704c1d`(piecePullbackIsoStage)・`68ddc617`(訂正の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(piece_preimage_eq完成、集計10/24で変わらず)。`ExtLimit.lean`に
`piece_preimage_eq`(Ext X側のUのアフィン片=R'側の同じアフィン片の
`extConePi X .app R'`によるpreimage、★sorry無し)を追加——
`piecePullbackIso`と`piecePullbackIsoStage`が「同じU」を指している
ことの位相的な裏付け。`extConePi_app_fst`+`Scheme.Hom.comp_preimage`
だけで閉じた。

コミット: `63ce023b`(piece_preimage_eq)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(「単発では閉じない」を実際に越えた、集計10/24で変わらず)。前ターンで
「作業単位1の第一歩は単発では閉じない」と記録した直後、実際に腰を
据めて取り組み、`isLocalization_closure_away_mul`を完成させた
(`FieldLimit.lean`、★sorry無し)——`Localization.Away(a*b)`が
`Submonoid.closure{a,b}`に関しても局所化になっていることを、
`IsLocalization.isLocalization_of_is_exists_mul_mem`+
`closure_induction`(4ケース分解)で示した。作業単位1(重なりの比較)に
要る最初の代数的補題が完成——残る「もう半分」(`Away(a)`をさらに`b`で
局所化したものも同じ`closure{a,b}`に関する局所化になること、
`localizationLocalizationSubmodule`経由)はまだ未完成、次の一手として
記録。★教訓:「単発では閉じない」は正しかったが「引き継ぐしかない」を
意味しなかった——腰を据めれば同じセッション内でも前進できた。

★続報: 「もう半分」の難点の正体を特定した——`mem_localization
LocalizationSubmodule`が`algebraMap`での**像としての**等式を与えるだけ
で、`R`自身での`closure{a,b}`membership(文字通りの等式)とは`a`-捩れ元
の分だけ異なりうると判明。名前探しではなく本物の可換環論の細部
(捩れの有無の場合分け)だった。`Γ(C,W_i)`の具体的な捩れ構造の検証を
次の一手として残し、この特定スレッドは打ち切り。
コミット: `15e0781b`(難点の正体を記録)。

コミット: `04c73661`(isLocalization_closure_away_mul)。★git注記:
corrhyp-goal.mdへの記録コミットが並行セッションの`81e23290`
(無関係なコミットメッセージ)に巻き込まれた——内容は正しく反映されて
いるが帰属が不正確。並行セッションの活動が活発化しているため今後も
`git diff --cached --stat`の確認を徹底する。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(技術ロードマップを統合版として記録、集計10/24で変わらず)。
`corrhyp-goal.md`に、これまで断片的に記録してきた`Lemma 4.1`
(`corrHypInstance4`)の残作業を1箇所に統合した——完成済み部品一覧・
残る3作業単位(比較射の構成・rigidity検証・GlueData組み立て)・
Corrのnonempty手当て・β脚(Z・ZKの確定)という対称的な残作業。
「relative Spec」の探索(`α_*O_C`という層を降ろす経路)も試みたが
mathlibに該当無しと確認、手作りのGlueData構成が唯一の道だと再確認。

正直な見積もり: 残る3作業単位(特に比較射の構成)は今セッションで
積み上げた降下補題群一式に匹敵するかそれ以上の分量になる可能性が
高い。`Lemma 4.1`の完全な形式化は複数セッションにまたがる継続的な
取り組みとして扱う——★依然「壁ではなく道」であり、今セッションで
その道の大部分(数学的な核心と土台のほぼ全て)を切り開いた。

コミット: `5a297b8e`(ロードマップ統合版)。

★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(残作業4段階を総括、
集計10/24で変わらず)。generic flatnessの標準的前提(Rネター整域・
加群/多元環がR上有限表示)を確認したところ、`corrHypInstance4`の
`Space:=QcqsSpace`はqcqsのみで有限型を要求しておらず、そのままでは
適用できないと判明——`Space`をさらに「qcqsかつ有限型」に絞った新
instanceが要る。`Lemma 4.1`の完全な証明に要る残作業は
**(1)Spaceの有限型への絞り込み・(2)generic flatnessの形式化・
(3)複数片の合流・(4)GlueDataでの貼り合わせ**の4段階と判明——(2)は
独立した本物の代数の定理として新規に要り、1回のセッションで詰め切れる
規模を超える。複数セッションにまたがる取り組みとして継続する。

コミット: `0bd1e7e1`(4段階の総括)。

★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(訂正——「捩れ」の心配は
不要な回り道だったと判明、作業単位1の見通しが簡潔化、集計10/24で
変わらず)。前回記録した`a`-捩れ元の場合分け(`Localization.Away x`を
さらに`y`で局所化したものが`Away(x*y)`と抽象的に一致するかという環
レベルの同型の追求)を再検討し、そもそも不要な回り道だったと気づいた。
`glueMorphisms`が要求するのは「重なり上で2つの射が一致する」という
**スキームレベルの条件**(`pullback.fst≫f_l = pullback.snd≫f_m`)であり、
これは`D(f_l)∩D(f_m) = D(f_l·f_m)`という**位相的な事実**(Zariski開基の
標準性質、捩れの有無に無関係に常に成立)と、`D(f_l·f_m)`の`D(f_l)`・
`D(f_m)`への標準的な開埋め込み(`PrimeSpectrum.basicOpen`の包含、これも
常に成立)だけで組み立てられる——環レベルの抽象的な同型を経由する必要は
無いと判明。

**訂正後の見通し**: 比較射の構成で本当に要るのは、(a)`D(f_l·f_m)`が
`D(f_l)`・`D(f_m)`双方の標準的な開部分集合であること(捩れ無関係、常に
成立)、(b)`f_l`・`f_m`自体を有限段階`R'`の元として認識すること
(`P₀_l`・`P₀_m`の構成に使った`StandardEtalePresentation`の生成元`x`との
対応づけ)。(b)は依然未着手だが、(a)の「捩れの場合分け」は不要と判明
したので、作業単位1(比較射の構成)の見通しは**むしろ簡潔になった**。
★教訓: 深掘りすると回り道自体が見えてくる——`isLocalization_closure_
away_mul`(★sorry無し、今も有効な補題として残る)は無駄ではなかったが、
その先の「もう半分」の追求自体が不要だった。

コミット: `27ba8b79`(訂正の記録、corrhyp-goal.md)。

★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(作業単位1(b)へ着手、
`StandardEtalePair.Ring`の元を降ろす2変数多項式版descentが完成、
集計10/24で変わらず)。`exists_fg_subalgebra_tensor_bivariate_finset`
(`FieldLimit.lean`、★sorry無し)——`Polynomial (Polynomial (A⊗[ℚ]ℝ))`
(`StandardEtalePair.Ring`が経由する`Bivariate`多項式環そのもの)の
有限個の元を共通の`R : FgSubalgebra ℚ ℝ`へ同時に降ろす、
`exists_fg_subalgebra_tensor_finset`の2変数版。`lake build`で
ABC3全体(6590 jobs)0エラーを確認。

配管の教訓(`tools/lean-idioms.md` #26): `Polynomial.mapRingHom f`の
`FunLike`適用と`.map f`(dot記法)は定義上等しい(`rfl`で確認済み)のに
構文上一致せず`rw`/`simp_rw`が刺さらない——`have hcoe : ⇑(mapRingHom φ)
= Polynomial.map φ := rfl`を明示的に挟んで解決。

**残る接続**(未着手): `StandardEtalePair.Ring`は`Bivariate`多項式環の商
なので、商から代表元への持ち上げ→この補題で係数降下→有限段階側でも
`Ideal.Quotient.mk`を適用、という一手で比較射の分母`g_l`等の実際の元を
降ろせるはず。ただし`standardEtalePairRingBaseChange`自体は
`equivMvPolynomialQuotient`経由で組み立てられているため、両者の往復
(`Bivariate.equivMvPolynomial`)も必要——次の一手として記録。

コミット: `8c96fa00`(bivariate descent完成)・`8e78ac0f`(記録)・
`5e2cb5ec`(lean-idioms #26)。

★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(作業単位1(b)完成、
`P.Ring`の任意の元を有限段階へ降ろせるようになった、集計10/24で
変わらず)。`exists_fg_subalgebra_tensor_standardEtale_elem`
(`FieldLimit.lean`、★sorry無し)完成——`P : StandardEtalePair(A⊗ℚℝ)`の
任意の元`z`について、ある有限段階`R`・`P₀`・`z₀`が存在し、`P₀`のbase
changeが`P`に一致し`z₀`の像が`z`に一致する。副産物として
`standardEtalePairMapRingHom`(`P.Ring→+*(P.map f).Ring`、`Ideal.
quotientMap`直接使用、`equivMvPolynomialQuotient`経由は不要と判明)・
`exists_fg_subalgebra_tensor_standardEtalePair_mapEq`(構造体一致版)・
`bivariateIsTwoSided`(組み合わせ爆発するIsTwoSided探索をletIチェーンの
global instance化で解決)も完成。

配管の教訓2件(`tools/lean-idioms.md` #26・#27): 依存キャスト`▸`を
STATEMENT内に埋め込むと`whnf`がcombinatorial explosionでtimeoutする
——`HEq`で述べると軽い。`open...in`/`set_option...in`はdocstringの
**前**に置く必要がある(後だと「unexpected token; expected 'lemma'」、
`lean_check`では検出できず`lake build`で発覚)。

**これで作業単位1(比較射の構成)の見通しが完成した**: (a)位相的事実
(捩れ無関係、前回訂正済み)+(b)`f_l`・`f_m`を有限段階の元として認識
(今回完成)。残るのは実際の組み立てのみ。

コミット: `8e605fff`(lean-idioms #27)・`fda225b4`(作業単位1(b)完成)・
`52cb625c`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(作業単位1が完全に
「道具が出揃った」状態に、集計10/24で変わらず)。作業単位1の残り(a)
「`D(f_l)∩D(f_m)=D(f_l·f_m)`」を実際に書こうとしたところ、mathlibに
`AlgebraicGeometry.Scheme.basicOpen_mul`としてそのまま存在しており
ABC3側で新規補題は不要と判明。**これで作業単位1(比較射の構成)は
(a)(b)とも完成**——残るのはこれらを組み合わせた実装(数十〜百行規模)・
work unit 2(rigidity)・work unit 3(GlueData/glueMorphisms組み立て)、
という複数セッションにまたがる継続タスクのみ。

Corrのnonempty脱落(item(4))も再検討: `Corr.comp'`が`D.pullback`を使う
ため、一般には(2つのFEt射が連結な底上で必ずしも交わらないなら)C の
非空性は自動保存されず、`Corr`構造体全体への追加は`comp'`/`transpose`/
`extCorr`すべてに波及する大きめの変更になると判明——`lemma_4_1`自体の
statementへ`[Nonempty c.C]`を追加前提として持ち込む方が影響範囲が狭く、
実際に`lemma_4_1`を書くときまで先送りするのが合理的と判断。

コミット: `3b81c2b6`(作業単位1(a)確認の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(作業単位3
「GlueData組み立て」に着手、候補局所片を実際のスキームとして実現、
集計10/24で変わらず)。`ExtLimit.lean`に`standardEtalePairSpecMap`
(`Spec P.Ring ⟶ Spec R`)・`standardEtalePairSpecMap_etale`(常にétale、
`HasRingHomProperty.Spec_iff`+`RingHom.etale_algebraMap`)・
`standardEtalePairPullbackIso`(候補局所片を`K`へbase changeすると実際の
`K`段階の局所片に一致、`pullbackSpecIso`+`standardEtalePairRingBaseChange`)
を完成(★すべてsorry無し)。`c'.C`の各局所片が実際にスキームとして構成
でき、étaleかつbase changeで元の局所片に戻ることを保証する、GlueData
組み立ての核心部品の第一歩。

**新発見**: 標準エタール片はétaleだが有限とは限らない(局所化+商、
properではない)——`FEt`が要求する`IsFinite`はGlueDataで貼り合わせた
**後**に大域的な性質として別途確認する必要があると明示的に言語化した
(当初のロードマップに暗黙に含まれていたが今回初めて明記)。

コミット: `ca88ac3a`(standardEtalePairSpecMap等完成)・`c3ce0e9e`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(Cの局所片が
アフィンであることを確認、work unit 3続き、集計10/24で変わらず)。
`piece_preimage_isAffineOpen`・`piecePreimageIso`(`ExtLimit.lean`、
★sorry無し)完成——`IsFinite α → IsAffineHom α`+`piece_isAffineOpen`+
`IsAffineOpen.preimage`(mathlib)で`α⁻¹(piece)`がアフィンであること、
`α⁻¹(piece) ≅ Spec Γ(C,α⁻¹(piece))`を確認。`piece_algebraEtale_tensor`
の環レベルの主張が実際のスキーム構造を正しく捉えていることの土台。
次の一手: `D(f_i) ≅ Spec(Localization.Away f_i)`(mathlib
`basicOpenIsoSpecAway`)と`standardEtalePairPullbackIso`をつなげる。

コミット: `8b4c1b2b`(piece_preimage_isAffineOpen等)・`2bdda863`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(一般の
アフィン開上のbasicOpen-Spec同一視完成、1ピース分の完全な連結は
配管の詰まりで未完了、集計10/24で変わらず)。`IsAffineOpen.
basicOpenIsoSpecAway`(`ExtLimit.lean`、★sorry無し)完成——mathlibの
`basicOpenIsoSpecAway`(`X:=Spec R`限定)を一般のアフィン開版へ一般化。
これで「1ピース分の完全な連結」に要る個々の部品(`standardEtalePairPullbackIso`・
`IsAffineOpen.basicOpenIsoSpecAway`・`_mapEq`)はすべて完成したが、実際に
3つを合成しようとしたところ`letI`導入の`Algebra`インスタンスと
`algebraMap`の構文不一致で複数回type mismatchに当たり、本セッションでは
完成させられず——正直に未完了として記録(ファイル内コメント+
corrhyp-goal.mdに詳細)。次回は`show`で型を都度明示しながら再開する。

コミット: `0dc8dc74`(basicOpenIsoSpecAway完成)・`48ad5c26`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(「1ピース分の
完全な連結」が完成、前回の詰まりを解消、集計10/24で変わらず)。前回の
詰まり(`algebraMap`と`letI`導入の`Algebra`インスタンスの構文不一致)を
解消——`exists_fg_subalgebra_tensor_standardEtalePair_mapEq`が返す等式を
`letI`のスコープ内で`▸`ではなく**ただの`:=`(defeq)**で`algebraMap`形へ
代入できると気づき、`standardEtalePairPullbackIso`を再証明せずそのまま
使えた。`onePieceSchemeIso`・`piece_descends_iso`(`ExtLimit.lean`、
★sorry無し)完成——**任意のスキームのアフィン開上の標準エタール元`f`に
ついて、`X.basicOpen f`が有限段階の候補局所片のbase changeにちょうど
一致する**という、作業単位1・3の核心の合流点。

配管の教訓(`tools/lean-idioms.md` #28): `letI`導入のインスタンスに
依存する等式の形変換は`▸`より先に型注釈つきの`:=`(defeq)を試す。

**これで`Lemma 4.1`の「1アフィン片・1標準エタール片」の完全な降下が
数学的に確立した**——残るのは複数ピースの貼り合わせ(GlueData本体)・
work unit 2(rigidity)・貼り合わせ後の有限性確認の3点。

コミット: `2f34ca1f`(1ピース連結完成)・`3ee14950`(lean-idioms #28)・
`d10ec816`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(複数ピース
比較への橋渡し、`standardEtalePresentation_exists_coord`完成、集計10/24で
変わらず)。`FieldLimit.lean`に完成(★sorry無し)——`Pres`の任意の元は
`equivRing`で運ぶと「`g^n`倍すれば`X`の多項式になる」座標表示を持つ。
複数の標準エタール片`D(f_i)`・`D(f_j)`を貼り合わせるのに要る「`f_j`が
`D(f_i)`の中でどんな元に対応するか」を得る第一歩——得られる多項式は
1変数なので既存の`exists_fg_subalgebra_tensor_polynomial_family`で
そのまま有限段階へ降ろせる。GlueDataの遷移射構成に直接使える次の一手。

コミット: `43b30c97`(standardEtalePresentation_exists_coord)・
`34f15083`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(重なりの
遷移開集合をP.Ring内の基本開集合として同定、GlueData遷移射の核心が
完成、集計10/24で変わらず)。`standardEtalePresentation_transitionOpen_eq`
(`FieldLimit.lean`、★sorry無し)完成——`D(f_i)∩D(f_j)`が`P_i.Ring`の中の
具体的な基本開集合(座標多項式`p`を`X`で評価した元)にちょうど一致する
ことを示した。`g`が`P.Ring`の単元であること(標準エタール表示の定義
そのもの)から`g^n`倍が`basicOpen`を変えないという事実を使う(副産物
`PrimeSpectrum.basicOpen_eq_top_of_isUnit`・`basicOpen_mul_isUnit`も
完成)。

**これでGlueDataの遷移射構成の核心が完成した**——残るのは(i)この`p`を
有限段階へ降ろす、(ii)`D(l→m)`・`D(m→l)`間の実際の環同型を構成する、
(iii)複数ピースを`Scheme.GlueData`として組み立てる、の3段階。

コミット: `e6e8c98f`(transitionOpen_eq完成)・`8ebf05cb`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報(重要な
発見——遷移射に抽象的な環同型は再び不要と判明、今度は根拠つき、集計
10/24で変わらず、新規sorry無し宣言は無し)。GlueDataの遷移射を
`IsLocalization.localizationLocalizationSubmodule`経由の抽象的な環同型
(`llsm(powers a)(powers algebraMap b)` vs逆順が同じ部分モノイドか)で
構成しようとしたが、`isLocalization_of_is_exists_mul_mem`の向きが
逆(既知の小さい部分モノイドから大きいへ拡張する形で、今回欲しい
向きと逆)であることに気づいた——これは以前記録した「捩れ」の問題と
同じ壁だった。

**しかし今回`piece_descends_iso`(前ターン完成)のおかげでこの抽象的な
環同型は不要と判明した**——`piece_descends_iso`は既存のスキーム`C`の
具体的な開集合からの同型を与えるので、重なり`X.basicOpen(f_i·f_j)`
への**制限**として各片の対応する開部分への同型が得られ、その**合成**
が遷移射そのものになる(cocycle条件も自動、両方とも同じ`C`を経由する
ため)。以前の「不要」という訂正と同じ結論だが、今回は`piece_descends_
iso`という具体的根拠まで遡って確認できた点が違う——`Scheme.Cover.
glueMorphisms`(既存スキームへの比較射構成)の枠組みが対応する。

次の一手: `piece_descends_iso`の同型を`X.basicOpen(f_i·f_j)`へ制限する
具体的な構成(`IsOpenImmersion`のrestrict機構)。

コミット: `e72e98f3`(発見の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(GlueDataの遷移射構成が完成——「捩れ」の壁を今度は完全に迂回できた、
集計10/24で変わらず)。前回の発見(抽象的な環同型は不要、
`piece_descends_iso`の制限で足りるはず)を実際に構成し切った。
`Scheme.hom_image_iso_eq_inv_preimage`(スキーム同型`e:X≅Z`で
`e.hom''ᵁW=e.inv⁻¹ᵁW`)・`exists_transitionOpen_eq_basicOpen`
(重なり`X.basicOpen(f₁·f₂)`が`X.basicOpen f₁`とその候補片`Z`の任意の
同型`e`の下で`Z`の基本開集合として実現される)を`ExtLimit.lean`に
完成(★すべてsorry無し)——`Scheme.basicOpen_res`・`Scheme.Opens.
ι_image_basicOpen_topIso_inv`・`Scheme.preimage_basicOpen`という
mathlib既存の自然性補題3つを合成するだけで閉じた。

**これでGlueDataの遷移射構成の核心部品が完成した**——2つの標準エタール
片をそれぞれの候補片へ写したとき、両方とも同じ`X`内の開集合を経由する
ので、遷移射は各々の同型のこの記述を合成するだけで得られ、抽象的な
環レベルの独立検証は本当に不要だった。残るのは(i)複数片を横断した
実際の組み立て(`Scheme.GlueData`)、(ii)cocycle条件の確認、
(iii)貼り合わせ後の有限性確認、の3段階。

コミット: `c9c711bd`(遷移射構成完成)・`fb049a13`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(GlueData組み立て戦略を簡略化、cocycle条件はmathlibの`gluedCover`を
再利用できると判明、集計10/24で変わらず、新規sorry無し宣言は今回は
無し)。`CategoryTheory.GlueData`の全フィールド(`t'`・`cocycle`等)を
自前で証明しようとしたが、**mathlibの`Scheme.Cover.gluedCover :
X.OpenCover → Scheme.GlueData`が既存のスキームの開被覆からのGlueData
についてcocycleをすでに一般的に証明している**と判明。戦略:
(1)`C`自身の開被覆を`Scheme.openCoverOfIsOpenCover`で作る、
(2)`gluedCover`(cocycle自動)、(3)`piece_descends_iso`の同型族に沿って
移送——移送されたGlueDataもcocycleを自動的に保つ(圏論の一般論)ので
再検証不要。「t'・cocycleを自前で証明する」という最大の技術的負債が
解消される見通し。

続けて戦略の(1)(ring側の`PrimeSpectrum.basicOpen`被覆をscheme側の
`X.basicOpen`被覆へ運ぶ)を試みたが、`⨆i∈t,X.basicOpen(f i)`と
`⨆g:(f''t),X.basicOpen g`の同一視で`Γ(X,U)`関連の型がinstances透明度
でtype-correctでなくなるという本セッションで繰り返す詰まりに再度当たり
完成させられなかった——正直に未完了として記録、コミットはせず。

コミット: `b661e01e`(戦略簡略化の記録)・`fd443564`(未完了の記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(ring→scheme被覆変換が完成、前回の詰まりを別経路で回避、集計10/24で
変わらず)。前回未完了だった変換を`hU.fromSpec`(`Spec Γ(X,U) ⟶ X`)
経由で完成させた——`fromSpec_preimage_basicOpen`+`image_preimage_eq_
opensRange_inf`を組み合わせるだけ、添字づけを変える必要が無いので前回
当たったinstances透明度の詰まりを完全に回避できた。
`exists_scheme_basicOpen_cover_of_ring`(`ExtLimit.lean`、★sorry無し)
完成。

**これで`exists_finite_standardEtaleCover`の環レベルの被覆を`C`内の
開被覆として認識する橋渡しが完成**——`Scheme.openCoverOfIsOpenCover`
+`Scheme.Cover.gluedCover`(cocycle自動証明済み)を適用する道が開けた。
GlueData組み立て戦略の(a)が完成、残るは(b)同型による移送の具体的な
構成、(c)貼り合わせ後の有限性確認、の2点。

コミット: `6e8ea04e`(cover変換完成)・`6d9b1db1`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(遷移射を実際のスキーム同型として取り出す完成、GlueDataのtフィールド
材料が揃った、集計10/24で変わらず)。`exists_transitionIso`
(`ExtLimit.lean`、★sorry無し)完成——`exists_transitionOpen_eq_
basicOpen`が与える開集合の等式を、`Scheme.Hom.isoImage`+`eqToIso`で
**実際のスキーム同型**へ格上げした。2つの片`i`・`j`に適用すれば、両方
とも同じ`X.basicOpen(f_i·f_j)`と同型になるので合成で遷移射
`Z_i.basicOpen(s_i)≅Z_j.basicOpen(s_j)`が直接得られる。

**現状整理**: GlueData部品(a)ring→scheme被覆変換・(b)候補片の実現・
(c)遷移射、まで揃った。残るは(d)これらを`Scheme.GlueData`の実際の
12フィールドとして組み立てる作業(`t'`・`cocycle`は3重の重なりでの
整合性、まだ手つかず)、(e)貼り合わせ後の有限性確認、の2点。

コミット: `925233b1`(exists_transitionIso完成)・`63edfedc`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(3重の重なりも新しい数学を要さないと判明、GlueDataは配線のみ残る、
集計10/24で変わらず)。`exists_transitionIso`が`f₂`に**任意の元**
(たとえば`f_j*f_k`という積)を渡せることに気づき、3重の重なりも
**同じ補題の1回の適用**でカバーされると確認した(`exists_transitionIso_
finset`、`ExtLimit.lean`、★sorry無し)。`Z_i.basicOpen(s_ij)⊓Z_i.
basicOpen(s_ik)=Z_i.basicOpen(s_ij·s_ik)`(`basicOpen_mul`)と環準同型の
乗法性から、この「1回適用」の結果が個別の2つの適用の積に一致すること
も直接確認できる。

**これでGlueDataの数学的な核心はすべて完成した**——残る作業は
「これらの部品を`CategoryTheory.GlueData`の12フィールドとして実際に
配線する」という、新しい数学を要さない(が分量のある)エンジニアリング
作業のみに帰着した。`TopCat.GlueData.MkCore`(トポロジー版の簡略化
コンストラクタ、cocycleが`rfl`で自動)のScheme版はmathlibに見当たらず
——12フィールドの直接組み立てが唯一の道と確認済み。

コミット: `68345679`(発見の記録)・`dd1354ea`(exists_transitionIso_finset)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに続報
(GlueDataの配線に着手、V・f・tフィールド完成、集計10/24で変わらず)。
`gdV`・`gdF`・`gdT`(`ExtLimit.lean`、★sorry無し)完成——`piece_
descends_iso`・`exists_transitionIso`を実際に`CategoryTheory.GlueData`
の対象レベルのフィールドへ流し込んだ。特に`gdT`(遷移射`t(i,j):
V(i,j)≅V(j,i)`)は`exists_transitionIso`をi側・j側それぞれに適用して
得た2つの同型を、乗法の可換性で繋ぐだけで完成した。

残る作業: `f_mono`・`f_hasPullback`・`f_id`(basicOpen.ιの標準的性質から
自動のはず)・`t_id`・`t'`(`exists_transitionIso_finset`から)・
`t_fac`・`cocycle`(選択項の一致性の確認、最も技術的な段)。これらを
埋めれば`Scheme.GlueData`が完成する。

コミット: `b79fab9b`(gdV・gdF・gdT完成)・`fd8e47c0`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに
続報(GlueDataのt_id・f_mono・f_open・f_hasPullback完成、集計10/24で
変わらず)。`gdT_id`・`gdF_mono`・`gdF_isOpenImmersion`・
`gdF_hasPullback`(`ExtLimit.lean`、★sorry無し)完成——12フィールド
のうち`V`・`f`・`t`・`t_id`・`f_mono`・`f_open`・`f_hasPullback`(7個)
が完成した。

残る`f_id`(対角成分が同型)は、選ばれる座標`s`(`exists_transitionIso`
の`.choose`、`Classical.choice`経由)が単元であることを示す必要がある
と判明——`i=j`限定の強化版補題が別途要る(`RingedSpace.isUnit_res_
basicOpen`で基本開集合の定義元がその上で単元であることは確認済み)。
`t'`・`t_fac`・`cocycle`(5フィールド)も残る。

コミット: `5653409e`(t_id等4フィールド完成)・`abc47300`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに
続報(訂正——gdV/gdF/gdTを精密な等式を保つ形に作り直し、f_id完成、
10/12フィールドが揃った、集計10/24で変わらず)。`f_id`に取り組む中で
**設計上の不備**を発見した——`gdV`/`gdF`/`gdT`は`exists_transitionIso`
(`Nonempty`で包んだ弱い実存)で定義していたが、これだと「選ばれた座標`s`
が単元である」という`f_id`に要る事実を`Classical.choice`の不透明性
のせいで`.choose`から引き出せない(`obtain`と`.choose`が定義上等しく
ないことも実測確認)。

**修正**: `exists_transitionOpen_eq_basicOpen`(精密な等式`e.hom''ᵁW=
Z.basicOpen s`をそのまま保持)の`.choose`を直接使うように`gdV`・`gdF`・
`gdT`を定義し直した。配管の教訓(`tools/lean-idioms.md` #29):`def`内で
`obtain`を使うと`And.rec`が残り後の`unfold`+`simp`がスタックする——
`let`+`.1`/`.2`射影を使えば通る。

新規完成(★すべてsorry無し): `gdV_diag_eq_top`・`gdF_id`(f_idそのもの)・
`Scheme.hom_image_top_eq_top`・`isIso_ι_of_eq_top`・再証明した`gdT_id`・
`gdF_mono`・`gdF_isOpenImmersion`・`gdF_hasPullback`。

**これで`Scheme.GlueData`の12フィールド中10個(`J・U・V・f・t・t_id・
f_mono・f_open・f_hasPullback・f_id`)が完成した**。残るは`t'`・
`t_fac`・`cocycle`のみ。

コミット: `d787b8e3`(訂正・f_id完成)・`ef3b03ce`(lean-idioms #29)・
`f070dbfb`(記録)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに
続報(transitionElem導入——GlueDataのt'完成、12フィールド中11個達成、
集計10/24で変わらず)。`t'`(3重の重なりの整合性)に取り組む中で、
`f_id`を直した`.choose`ベースの設計もまた限界に当たった——`t'`は
「異なる`.choose`呼び出し同士を比較する」ことを要求するが、
`Classical.choice`の不透明性の下ではこれが原理的に不可能(トイ例で
`rfl`が失敗することを実測確認: `(⟨w,h⟩:∃x,x=w).choose = w`は`rfl`
で通らない)。

**修正**: `∃`を一切経由しない決定的な定義`transitionElem`(`e.inv`・
`topIso.inv`・制限写像の合成——witnessを構成する式をそのまま`def`に
した)を導入した。存在命題を経由しないので`transitionElem_mul`(乗法性、
`map_mul`3回+`rfl`)が自由に示せ、これが`t'`の構成を可能にした。
`gdV`・`gdF`・`gdT`を`transitionElem`の上に作り直し、新規に
`gdVpullbackIso`・`gdT'`(`t'`そのもの)を完成させた(★すべてsorry無し)。

配管の教訓: `transitionElem`と周辺補題は、族`f・Z・e`を導入する
`variable`宣言より前に単体の`X・U・Z・e`で自己完結的に定義する必要が
ある——`variable {f Z e} in`で個別除外を試みたところsection全体の
変数解決が壊れる罠に当たった。

**これで`Scheme.GlueData`の12フィールド中11個(`J・U・V・f・t・t_id・
f_mono・f_open・f_hasPullback・f_id・t'`)が完成した**。残るは`t_fac`・
`cocycle`のみ(新しい数学は不要、`gdT_id`同様の`simp`主体の計算の見込み)。
`lake build ABC3`で0エラー確認。

コミット: `21aea8c6`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04さらに
続報(t_facに向けた配管の罠を特定・修正、証明を途中まで前進、集計10/24
で変わらず)。`t_fac`に着手し、`isPullback_opens_inf`系のsimp補題が
`gdF`(不透明な`def`)を使った`pullback.fst/snd`に一切効かないという
罠に当たった——原因は`gdF`が`instances`透明度で展開されないこと(simpの
書き換え一致判定はその透明度を使うため)。`gdV`・`gdF`に`@[reducible]`
を付けて解消(コミット`e41d967e`、`tools/lean-idioms.md` #30)。

この修正の上で`t_fac`を`cancel_mono`+明示的な`rw [show pullback.fst
… = (isPullback_opens_inf _ _).isoPullback.inv ≫ (Z i).homOfLE
inf_le_left from …]`+`Iso.cancel_iso_inv_left`まで前進させた——
`pullback`/`isoPullback`が完全に消えた「純粋に`transitionElemIso`・
`isoImage`・`eqToIso`・`homOfLE`の塔同士の一致」という等式まで還元
できている(まだ`sorry`、Foundには未反映)。次の一手はこの塔同士の
突き合わせ。詳細は`corrhyp-goal.md`の該当節。

コミット: `e41d967e`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(instances透明度の壁を突破、transitionElemIso_inv_naturality完成、
t_facに必要な部品が全て揃った、集計10/24で変わらず)。3ターン連続で
当たっていた「rw/simp/convで事実を差し込むとX.presheafの型不一致」
という壁を突破した。鍵は2つ: (1) rwを一切使わずcalc+congrArg+
(Category.assoc _ _ _).symmのterm-modeで組み立てる——congrArgは
rwのようなmotive探索(kabstract)を経由しないため壁に当たらない。
(2) 各部品をhaveとして内部に書くとそれだけでwhnfのheartbeat上限
(400万でも不足)に達する——独立したtheorem(transitionElemIso_step1〜
step45)として先に証明し切ってから本体のcalcでは名前を参照するだけに
することで軽くなった。

`transitionElemIso_inv_naturality`(t_facの核心——transitionElemIso
f₁(g·h)eの.invを両側の基本開集合の包含を挟んでtransitionElemIso f₁ g e
の.invへ橋渡しする自然性)が完成(sorry無し)。これでt_facに必要な数学的・
技術的部品がすべて揃った——残るのはgdT'・gdTとの配線のみ、新しい未知の
壁はもう残っていない。`tools/lean-idioms.md` #31に発見を記録。

コミット: `c4172c85`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(**t_fac完成——GlueDataの12フィールド中11個達成**、集計10/24で
変わらず)。前回揃えた部品を実際に組み立てた。追加で必要だった部品:
`gdT_eq_transitionElemIso`(`gdT`が`transitionElemIso`2つの橋渡しと
`unfold;rfl`で直接一致)・`gdVpullbackIso_hom_comp_homOfLE_fst/snd`・
`gdVpullbackIso_inv_comp_pullback_snd`(pullback.fst/sndをhomOfLEの
言葉に変換)・**積の順序を入れ替えた自然性版**(`transitionElemIso_
inv/hom_naturality2`——`gdT'`のj側は生き残る元が積の後ろ(`h*g`)に
来るため、前回の`g*h`型では直接使えず独立に証明し直した。骨格は
`inf_le_left/right`が入れ替わるだけで完全に同じ)・`gdT'_key`(i側・
j側で自然性を適用後、残る純粋なXレベルの結合律・可換律の等式——
`simp only [Scheme.isoOfEq_hom, Scheme.homOfLE_homOfLE]`だけで証明
無関係性により自動的に閉じた)・`gdT'_t_fac`(本体)。

**これで`Scheme.GlueData`の12フィールド中11個(`J・U・V・f・t・t_id・
f_mono・f_open・f_hasPullback・f_id・t'・t_fac`)が完成した**。残るは
`cocycle`のみ(新しい数学は不要、`gdT'_key`型の等式へ帰着させる見込み)。
`lake build ABC3`で0エラー確認。

コミット: `dfbc99a0`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(cocycleは数学的には解決済み、Leanの計算コストの壁で未完成、集計
10/24で変わらず)。`cocycle`に着手し、`unfold gdT';simp[…,Iso.inv_hom_
id_assoc,Iso.hom_inv_id_assoc,Scheme.isoOfEq_hom]`だけでほぼ全部が
自動的に打ち消し合い、残りは「3つの結合律のhomOfLEが自明な自己同型に
なる(hcombine、1秒未満で証明)」+「共通の同型を挟んだ逆射の対消滅
(iso_conj_cocycle_generic、完全に一般的、1秒未満で証明)」だけで閉じる
はずと確認した。しかし**組み立て**(この2つを実際にゴールへ適用する
exact/refine/calc)が`maxHeartbeats`2000万・タイムアウト590秒でも
完走しない——t_facのときの同じcalc+congrArgパターンは74秒で通ったのと
対照的で、3つの`gdT'`を同時に`unfold`することで生じる項の重さが本質的に
大きいと見られる。数学的内容には一切不明な点が無く、純粋にLeanの計算
資源の壁。次の候補: 3つの`gdT'`を同時にunfoldせず段階的に処理する、
または`gdT'`自体をtransitionElemIsoベースの明示的橋渡し補題として書き
直す。ファイルには未反映(コミット無し)。

コミット: `58c78aa6`(記録のみ)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(**cocycle完成——Scheme.GlueDataの12フィールドすべてが完成した**、
集計10/24で変わらず)。前回の壁(unfold gdT'を複数出現に同時適用すると
あらゆる型検査が極端に重い)を、「gdT'を1つずつ独立にunfoldせず名前
参照できる事実(gdT'_hom_eq/gdT'_inv_eq)として先に確定させ、以降は
rwでこれらの名前を参照するだけにする」ことで突破した。もう1つの鍵:
congrArgで書き換え対象を2つの合成の"間に挟む"(fun x => A≫x≫B)と重く
なる——常に「前だけ」/「後ろだけ」のcongrArgを順番に適用することで
劇的に軽くなった(0.05秒未満)。

gdT'_pair(残る2つのgdT'の合成が1つ目の逆射に一致する)+
Iso.hom_inv_idからgdT'_cocycleが直ちに従った。**これで
Scheme.GlueDataの12フィールドすべて(J・U・V・f・f_mono・
f_hasPullback・f_id・t・t_id・t'・t_fac・cocycle)が完成した**。
lake build ABC3で0エラー確認。tools/lean-idioms.md #32に発見を記録。
次の一手はScheme.GlueDataの構造体インスタンスを実際に組み立てること。

コミット: `9fc6c7ad`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(**corrHypGlueData完成——Scheme.GlueData構造体インスタンスを
組み立てた**、集計10/24で変わらず)。揃った12個の部品をそのまま
フィールドとして渡すだけで0.12秒で組み立てられた——`t`・`t'`は
`gdT`/`gdT'`(Iso)の`.hom`、`t_id`は`gdT_id`から`rfl`で閉じるだけ。
新しい数学もLeanの配管上の障害も無かった。**これでLemma 4.1の構成的
降下のうち、GlueDataの配線という最大の技術的部分が完了した**。
lake build ABC3で0エラー確認。

次の一手(規模が異なる、まだ未着手): corrHypGlueData.glued(貼り
合わせてできるスキーム)がC(またはExt後にC)に同型であることを、
C自身のOpenCoverのgluedCoverとの比較で示す——Lemma 4.1のc'.Cの実体。
ここまでは純粋な圏論的配線だったが、この段階は「候補片の族Zが実際に
Cを覆う」という降下データの普遍性に関わる、質の異なる作業になる見込み。

コミット: `923ad02b`。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(corrHypGlueDataから先の道筋を調査・記録、集計10/24で変わらず)。
lemma_4_1の正確な主張(Skeleton/CorrHyp/Section4.lean)・Corr構造体
(⟨C,α,β⟩)・Extの定義(Over BaseK上のX×_{Spec ℚ}Spec ℝ、
corrHypInstance4がQcqsSpace版の実働インスタンス)を確認。
corrHypGlueDataはX,U,J,f,Z,eが今なお完全に汎用パラメータのままで、
corrHypInstance4/Ext/piece_descends_isoのどれとも未接続と判明。

残る作業(いずれも未着手、質の異なる大きな作業): (a)具体化
(X:=c.C、piece_descends_isoの実存的witnessから決定的なZ,eを作る
必要があるかもしれない——transitionElemのときと同じClassical.choice
の罠の可能性)・(b)corrHypGlueData.glued≅Xという降下の整合性
(IsColimit同士の比較、genuinely new)・(c)piece_descends_iso(環
レベル)とExt(スキームレベル)のbase changeの記述の一致・(d)Cの
有限アフィン被覆+外側の貼り合わせ・(e)Corrのα・β脚と整合性の等式。
corrhyp-goal.mdに詳細を記録(コミット`b1c9414b`)。次に着手する際は
(a)から始めるのが妥当。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04
続報(corrHypGlueDataOfCover追加、項目(a)着手、集計10/24で変わらず)。
懸念していた「piece_descends_isoを決定的な形へ作り直す必要」は
**不要と判明**——corrHypGlueDataの12フィールドは異なる添字i,jの
Z i・Z jを互いに比較する場面が無く、各e iは添字ごとに独立に使われる
だけなので、.choose由来の不透明なZ i同士でも支障なし(transitionElem
のときの「複数.choose結果の代数的関係」問題とは構造的に別物)。

ExtLimit.leanに5定義を追加・lake build(ExtLimit/ABC3とも)0エラー
確認済みでコミット: descendPiece/descendPieceIso(piece_descends_iso
の.choose/.choose_specで単一fの候補片と比較同型を取得)・
descendPieceOfProof/descendPieceIsoOfProof(IsStandardEtaleを
typeclassでなく明示的証明として受け取るFinset用instance版)・
corrHypGlueDataOfCover(有限族f:ι→Γ(X,U)からJ:={i//i∈t}として
corrHypGlueDataを実際に呼び出す——初めて抽象パラメータでなく具体的
構成データでの呼び出し)。

正直な現状: 被覆条件⨆i∈t,X.basicOpen(f i)=Uはまだ未使用(項目(b)で
必要)、A/X/Uと実際のcorrHypInstance4/Ext/Cとの接続(項目(c))も未着手
——corrHypGlueDataOfCoverは依然CorrHyp非依存の再利用可能な部品。
corrhyp-goal.mdに詳細記録(コミット`20db4512`)。次の一手:
exists_finite_standardEtaleCover+exists_scheme_basicOpen_cover_of_ring
で実際のι,t,fと被覆証明を得る方向。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き)
続報(corrHypGlueDataOfEtale追加、上記の次の一手を実行、集計10/24で
変わらず)。exists_finite_standardEtaleCoverとcorrHypGlueDataOfCover
を合成し、Γ(X,U)がA⊗[ℚ]ℝ上étaleという1個の仮定だけからGlueDataを
自動組み立てるcorrHypGlueDataOfEtaleを追加。あわせて
exists_scheme_basicOpen_cover_of_ringでその族が実際にUを覆うことを
示すcorrHypGlueDataOfEtale_coverも追加(項目(b)が必要とする被覆条件
そのもの)。両方の定義で同一式
exists_finite_standardEtaleCover (A⊗[ℚ]ℝ) Γ(X,U)を再度書いて
.chooseする形にし、定義的に同じt・fを指すようにした。lake build
(ExtLimit/ABC3とも)0エラー確認、コミット: `ec21e781`。

**これで「étaleな環拡大→GlueData+被覆条件」というCorrHyp非依存の
一般的パッケージ化が完了**。残るのは(b)corrHypGlueData.glued≅Uの
証明(IsColimit比較、genuinely new)・(c)A・X・U・Etale仮定を実際の
corrHypInstance4/Ext/Cへ接続(Aが何であるべきかをpiece_descends_iso
の定義から逆算する必要あり)・(d)Cの有限アフィン被覆+外側の貼り合わせ
・(e)α・β脚と整合性の等式。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き2)
続報(corrPieceGlueData追加、項目(c)着手、集計10/24で変わらず)。
調査エージェント(読み取り専用)で仮説検証した結果、必要な部品は
すべて既存と判明: piece_algebraEtale_tensor(ExtLimit.lean:778、
以前構築済み)の結論がcorrHypGlueDataOfEtaleの要求形にちょうど一致、
piece_preimage_isAffineOpen(同:907)がhU(アフィン性)をIsFinite→
IsAffineHom(mathlib自動instance)経由で既に用意していた。

これらをそのまま代入するだけでcorrPieceGlueData(X:Over BaseK・
U:X.leftのアフィン開・C:Scheme・α:C⟶(ExtF.obj X).left・
[IsFinite α][Etale α]からScheme.GlueDataを得る)とその被覆条件
corrPieceGlueData_coverを追加。lake build(ExtLimit/ABC3とも)0エラー
確認、コミット: `71f72c2a`。

**これでロードマップ項目(c)はOver BaseK/ExtFレベルで実質完了**——
corrHypGlueDataOfEtaleが初めてExt/ExtFという実際のCorrHyp構成データに
直接接続された。残る「Corrの実データ(corrHypInstance4・QcqsSpace)
への代入」(項目(c'))はcorrPieceGlueDataのシグネチャがCorrのα
フィールドの中身と既に形が一致しているため短い配線になる見込み
(未着手)。残るのは(b)IsColimit比較(genuinely new)・(c')上記の代入
・(d)Cの有限アフィン被覆+外側の貼り合わせ・(e)α・β脚と整合性の等式。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き3)
続報(corrPieceGlueDataOfCorr追加、項目(c')完了、集計10/24で変わらず)。
Instance4.leanでcorrPieceGlueDataをCorr corrHypInstance4 (QcqsExt X) ZK
の実データ(c.C・c.α)へ実際に代入した。QcqsFEt A B:=FEtK A.1 B.1
なのでc.α.1.left/c.α.2.1/c.α.2.2をそのままcorrPieceGlueDataの
α/[IsFinite α]/[Etale α]へ渡すだけの短い配線になる見込みだった。

実行時、#31「instances透明度の壁」の新しい現れ方に当たった:
haveI:=c.α.2.1を直前に置いても、[IsFinite α]の型クラス探索が
(QcqsExt X).1.leftと(ExtF.obj X.1).left(QcqsExtが@[reducible]でない
defのため区別されてしまう)を見分けられずfailed to synthesize
instanceで失敗。修正: letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left :=
c.α.1.leftという呼び出し先の期待する型をそのまま構文に書いたletIで
hαを作り、以降hαだけを使うことで解消(元の式に型注釈を付けるのでは
なく「期待される型で包み直す」のが鍵)。lean-idioms.md#33に記録。

corrPieceGlueDataOfCorr/corrPieceGlueDataOfCorr_coverを追加、lake
build(Instance4/ABC3とも)0エラー確認、コミット: `2471cb91`。

**これでロードマップ項目(c')完了**——corrHypGlueDataOfEtaleという
抽象部品から出発し、Ext/ExtF(項目(c))を経て、corrHypInstance4・
QcqsSpace・Corrという実データ(c.C・c.α)にまで接続が完成した。残るは
(b)corrPieceGlueDataOfCorr(...).glued≅c.α⁻¹(piece)の証明(IsColimit
比較、genuinely new、最大の残りタスク)・(d)c.Cの有限アフィン被覆+
外側の貼り合わせ・(e)α・β脚と整合性の等式(β側はゼロから必要)。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き4、自律ループ)
続報(piecesOpenCover追加、項目(b)の第一歩、集計10/24で変わらず)。
自律ループのtickでExploreエージェントが項目(b)の土台を調査。結論:
mathlibに汎用の「2つの独立なGlueDataの成分ごと同型→.glued同士の同型」
という比較補題は無い——Scheme.Cover.gluedCover+fromGlued(IsIso、
既製)という特定の構成方法への結果のみ。corrHypGlueDataはgluedCoverの
形では組まれていない(gdVがpullbackでなくZ i自身のbasicOpen)ため、
新しい橋渡しがgenuinely newで必要(t_fac/cocycle同等以上の配線量、
という以前からの評価と一致)。

証明の道筋: corrHypGlueData_glued_iso : (corrHypGlueData f Z e).glued
≅ (U:Scheme)という形で、piecesOpenCover(Zの族からmkOfCoversで直接
構成したUのOpenCover、mathlibのgluedCover+fromGludedが既にU への同型
を与える)を比較対象に、corrHypGlueDataの.diagramとの間のNatIsoを
構成しHasColimit.isoOfNatIsoで合成する方針。U成分は恒等、V成分は
pullbackHomIsoLeft+pullbackSymmetry+isPullback_opens_inf(いずれも
既存)で組み立てられる見込み。

piecesOpenCover(比較対象のUのOpenCover構成)をExtLimit.leanへ追加、
lake build(ExtLimit/ABC3とも)0エラー確認、コミット: `9daf362d`。
途中#31と同系統の罠(全射性証明の被覆不等式のインライン展開でwhnf
heartbeat上限)にも当たり、独立theoremへの先出しで解消(既知教訓の
再確認)。

残る核心(未着手): V成分の同型構成とNatIso全体の組み立て。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き5、自律ループ)
続報(piecesGluedCoverVIso完成、V成分比較同型、集計10/24で変わらず)。
mathlibのpullbackIsPullbackOfCompMono(Mono iのときpullback f gは
pullback(f≫i)(g≫i)の極限錐でもある)をisPullback_opens_infと組み合わせ、
isPullback_opens_infのU相対版(A,B≤UのときX.homOfLE同士のpullbackも
A⊓Bに一致、pullbackHomOfLEIso)を発見(decl-index.mjs --mathlib+grep)。

Scheme.basicOpen_mul・pullbackHomIsoLeft+pullbackSymmetry(e追い出し、
既存)・transitionElemIso(既存)と3段のcalcで繋ぎ、
piecesGluedCoverVIso: pullback(𝒰.f i)(𝒰.f j)≅gdV f Z e (i,j)
(𝒰:=piecesOpenCover)を完成。lake build(ExtLimit/ABC3とも)0エラー
確認、コミット: `8717ce70`。

**これで項目(b)のNatIso構成に必要なU成分・V成分のうち、より本質的な
V成分が完成**。残る核心: φ(i,j)がfst/snd(gdF・gdT≫gdF vs
pullback.fst・pullbackSymmetry.hom≫pullback.fst)と可換であることの
証明(NatIsoの自然性、gdT'_t_fac/cocycle同等の配線量が見込まれる最大
の残りタスク)・U成分(恒等)と合わせたNatIso全体の組み立て・
HasColimit.isoOfNatIsoでの.glued比較・mathlibのfromGludedとの合成。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き6)
続報(piecesGluedCoverVIso_hom_fst完成、fst自然性完全確立、集計10/24で
変わらず)。φ(i,j)の3段構成それぞれのfst自然性を個別確立してから
内側から順に繋いだ: pullbackHomIsoLeft_hom_fst'/_hom_snd'(@[reassoc]、
pullback.lift_fst/_sndから直ちに従う)→pullbackEInvIso_hom_fst/_hom_snd
(pullbackSymmetry_hom_comp_fst/_snd既製と組み合わせ、show内は
(e i).symm.homで統一・非_assoc版と_assoc版を紙の計算に対応させて
使い分ける、という2つの配管の教訓を得た)→pullbackHomOfLEIso_hom_fst/
_hom_snd(IsPullback.isoIsPullback_hom_fst/_hom_snd既製)→
pullbackHomOfLEIsoBasicOpen_hom_fst/_hom_snd(Scheme.isoOfEq_hom+
Scheme.homOfLE_homOfLE)→eqToIso_congrArg_scheme_hom_ι(subst一発)→
transitionElemIso_hom_ι(核心、haveチェーンをsetで掴み直し、
cancel_monoで組み立てる新技法)→pullbackHomOfLE_gdV_iso_hom_fst/
piecesGluedCoverVIso_hom_fst(show+exactでdefeq経由、instances透明度
の壁を回避)。

lake build(ExtLimit/ABC3とも)0エラー確認、コミット: `402c3d7b`。

**これで項目(b)のNatIso自然性のうちfst側が完全に確立**。残るは
snd側(gdT≫gdFとの可換性、gdT_eq_transitionElemIso等が土台の見込み、
t_fac/cocycle同等の配線量)・U成分と合わせたNatIso全体の組み立て・
HasColimit.isoOfNatIsoでの.glued比較・mathlibのfromGludedとの合成。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き7)
続報(piecesGluedCoverVIso_hom_snd完成、**項目(b)のNatIso自然性が
fst・sndともに完全に確立**、集計10/24で変わらず)。gdT_eq_
transitionElemIso+transitionElemIso_hom_ιを組み合わせるgdT_hom_
comp_gdFが鍵。

新しい配管の罠: rw[gdT_eq_transitionElemIso]直後にshowで式全体を
再掲示すると、埋め込まれたX.isoOfEqの証明項同士(命題として同じだが
構文的に別)の照合でinstances透明度の壁に当たりwhnfが終わらない
(#31の新しい現れ方、maxHeartbeats 1000000でも76秒で打ち切り)。
修正: showを避けrw[Category.assoc,Category.assoc]で結合だけ先に
揃えてから適用(元の証明項に一切触れないので壁を回避)。

homOfLE_isoOfEq_comp'(@[reassoc])・gdT_hom_comp_gdF・
pullbackHomOfLE_gdV_iso_hom_snd・piecesGluedCoverVIso_hom_snd
(mathlib既製pullbackSymmetry_hom_comp_fstへ接続)を完成。lake build
(ExtLimit/ABC3とも)0エラー確認、コミット: `0fe99493`。

**残る道筋(項目(b)完成に向けて)**: (1)U成分(恒等)+V成分を
NatIso.ofComponentsでまとめる、(2)HasColimit.isoOfNatIsoで
corrHypGlueData.glued≅piecesOpenCover(...).gluedCover.gluedを得る、
(3)mathlibのfromGludedが与えるpiecesOpenCover(...).gluedCover.glued
≅Uと合成する。この3ステップでcorrHypGlueData.glued≅(U:Scheme)
(項目(b)そのもの)が完成する見込み。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き8)
続報(corrHypGlueData_toGluedCover完成、項目(b)の射が構成できた、
集計10/24で変わらず)。当初計画のNatIso.ofComponents+HasColimit.
isoOfNatIsoの代わりに、Multicoequalizer.desc(コクイザライザの普遍性)
を直接使うより直接的なルートを採用——corrHypGlueDataのU成分がZ(=
piecesOpenCoverのX成分)と文字通り同じなのでπ i:=gluedCover.ι iが
そのまま使え、整合条件(corrHypGlueData_compat)はpiecesGluedCoverVIso_
hom_fst/_hom_snd+GlueData.glue_condition(mathlib既製)から従う。

これまでで最も執拗な#31の壁に遭遇: gluedCover.J(piecesOpenCover.I₀、
定義的にはJに等しい)がinstances透明度でJと一致せず、rwはおろか
simp only[…] at hのようなローカル仮定への書き換えでも同じ壁に当たった
(set/letIでも解消しないケースがあった)。唯一有効だった方法: 中間の
haveを一切rw/simpで後から書き換えず、最初からcalc+congrArg+
Category.assoc(項として)だけで組み立てる——一度もrw/simpを使わずに
完走(1.82秒)。

corrHypGlueData_compat+corrHypGlueData_toGluedCoverを完成、lake build
(ExtLimit/ABC3とも)0エラー確認、コミット: `423e7d61`。

**これでcorrHypGlueData.gluedからgluedCover.gluedへの射が構成された**。
残る道筋(項目(b)最終段): (1)この射がIsIsoであることを示す(逆向きの
射をMulticoequalizer.descで構成、Multicoequalizer.hom_extで両方向の
合成が恒等射になることを確認)。(2)Iso.mkでまとめ、mathlibの
fromGludedとの合成でcorrHypGlueData.glued≅(U:Scheme)(項目(b)そのもの)
完成。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き9)
続報(**★★★ロードマップ項目(b)が完成★★★**、
corrHypGlueData.glued≅(U:Scheme)、集計10/24で変わらず)。残る道筋2
ステップを完了: gluedCover_compat(逆向きの整合条件、
corrHypGlueDataのGlueData.glue_condition+piecesGluedCoverVIso_hom_fst/
_hom_snd、φがisoであることでφ自体を打ち消す逆方向キャンセル)→
gluedCover_toCorrHypGlueData(逆向きの射、Multicoequalizer.desc)→
corrHypGlueData_toGluedCover_comp/gluedCover_toCorrHypGlueData_comp
(両方向の合成が恒等射、Multicoequalizer.hom_ext+π_desc。calcのTrans
自動解決が詰まりhave+.trans連鎖に切り替え)→corrHypGlueData_gluedIso
(Iso.mk)→**corrHypGlueData_glued_iso**(mathlibのScheme.Cover.
fromGludedとの合成、最終形)。

gluedCover_compatも一度もrw/simpを使わずcalc+congrArg+Category.assoc
だけで組み立てた(#31の壁の回避法の再現性を確認)。lake build
(ExtLimit/ABC3とも)0エラー確認、コミット: `540a7a22`。

**これで、corrHypGlueDataOfCover→corrHypGlueDataOfEtale→
corrPieceGlueData→corrPieceGlueDataOfCorr(項目(a)(c)(c'))→
piecesOpenCover→piecesGluedCoverVIso(項目(b)のNatIso比較)→
corrHypGlueData_glued_iso、という一連の作業がすべて1本の線でつながり、
「Corrの実データ(c.C・c.α)から具体的に構成した候補片の族が実際に
元のスキームを再構成する」というLemma 4.1の構成的降下の中核部分が
完全にsorry無しで証明された**。

残る作業(項目(d)・(e)): C(c.C.1.left)自体の有限アフィン被覆+外側の
貼り合わせ(二段階構成)、α・β脚と整合性の等式h▸extCorr D c'=cの構成
(β側はまだ手つかず)。集計は10/24で変わらず(Lemma 4.1というnumbered
item全体の内部構造であり、独立したnumbered itemとしてカウントされ
ないため)。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き10・11)
★重要★ lemma_4_1のh:ZK=D.Ext Zは文字通りには証明不可能(あるいは偽)
の疑いが判明、正直に記録(集計10/24で変わらず、§4は引き続き0/2)。

2つの独立した問題: (1)Corr定義(Section1.lean)がNonempty C(原典が
明記する条件)を欠いており、C:=∅の空相関がFEt条件を常に満たすため、
無関係なZKに対しCorr D (Ext X) ZKが常にinhabitedになる反例構成が
成り立つ。(2)h:ZK=D.Ext Zが文字通りの命題的等号を要求しているが、
QcqsSpaceは同型類の商ではない素朴なSubtypeなので、構成的に得たZが
ZKと文字通り等しくなる保証は一般に無い(同型にしかならない)。

朗報: X.1.leftの有限アフィン被覆(項目(d)第一段)は
Scheme.exists_finite_affineOpenCover(ExtLimit.lean:409)として既に
完成済み。片同士の外側の貼り合わせ(第二段)は約650行規模の新規
エンジニアリングが要る見込み。

(1)の修正(Corrに Nonempty C追加)の影響範囲を見積もった結果、
**拙速には着手しないと判断**。最も深刻: isIsogenous_refl(Section1.
lean:70、既にsorry無しで完成済み、§1の一部)は「Spaceの全要素は
非空」という保証が無いと成立しなくなる——しかもcorrHypInstance3/4
(スキーム)では空スキームが普通にSpaceに属するため、この汎用公理
自体が偽になる。回避には(a)isIsogenous_reflの主張を弱める(原論文
からの逸脱)か(b)Spaceを常に非空なスキームに絞った新インスタンスを
Instance3→4規模(数セッション)でもう一段作るか、のいずれかが必要
——**既に達成済みの§1 5/5を後退させかねない**規模の設計変更。

**判断**: Lemma 4.1(§4)は、GlueDataエンジニアリング(項目(a)〜(c')・
(b))という工数の多い道は完走できたが、その先(項目(d)第二段・項目
(e)全体・h:ZK=D.Ext Zの妥当性そのもの)はTheorem 3.3/4.2と同様、
person-yearsスケールの数学(またはHyperbolicCurveData interfaceの
再設計という独立した大規模プロジェクト)を要すると評価。
corrhyp-goal.mdに詳細記録(コミット`866db4e9`+今回分)。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き12)
続報(corrPieceGlueDataOfCorr_glued_iso完成、単一アフィン片レベルの
構成的降下が完結、設計論点には非依存の低リスク進捗、集計10/24で
変わらず)。項目(b)を一般形(corrHypGlueData_glued_iso)で完成させて
おいた成果を、Ext/ExtF・corrPieceGlueData・実際のCorrデータまで
一直線に特殊化: corrHypGlueDataOfEtale_glued_iso(iSup_subtype'で
Finset添字→部分型添字変換)→corrPieceGlueData_glued_iso(定義の一致で
特殊化)→corrPieceGlueDataOfCorr_glued_iso(Corrの実データへ特殊化、
corrPieceGlueDataOfCorr(...).glued≅c.α⁻¹(piece))。

以前「(b)IsColimit同士の比較、genuinely new」と見積もっていた作業が
新しい証明無しに特殊化だけで済んだ——「まず抽象的にCorrHyp非依存の
形で完成させておく」という方針の価値が実った。lake build
(ExtLimit/Instance4/ABC3とも)0エラー確認、コミット: `45ddff25`。

**これで単一アフィン片レベルの構成的降下の全工程(項目(a)〜(c')・
項目(b)の特殊化)が完結**。この作業はNonempty C/等式vs同型の設計
論点には一切触れておらず、どちらの解決を選んでも再利用できる
インフラ。残る道筋は変わらず(d)第二段(約650行規模)と(e)(設計論点
の解決が前提)。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き13)
★重要な訂正★ 「項目(d)」の理解を訂正——descendPiece(corrHypGlueDataの
Z i)はまだℝレベルのままで、k-スキームへの降下候補(Spec P₀.Ring、
Rレベル)はdescendPiece内部で使い捨てられ外へ取り出されていないと
判明。c.α⁻¹(Ext V_k)はC自身の中に既に開集合として存在するトート
ロジーなので「Cを外側で貼り合わせて再構成する」必要は本来無い。

本当に必要なのは、Spec(P₀.Ring)(Rレベルの候補片)をそれ自体Rレベルで
貼り合わせる新しい層——異なる添字i・異なるアフィン片V_kごとにR_iが
異なりうるため、共通精密化R_i⊔R_j(FgSubalgebraの結び)へ持ち上げて
から比較する新しい議論が要る。これはℝレベルでtransitionElem/gdT/
cocycleを構築したのと同種の困難を、Rレベル(FgSubalgebra ℚ ℝの圏)で
もう一段繰り返すことを意味する。

**結論**: 「項目(d)」は「Cの外側の貼り合わせ」ではなく「Rレベルでの
候補片の貼り合わせ」であり、新しい独立した規模の数学的内容——既存の
ℝレベル部品の配線だけでは済まない。person-yearsスケールという評価を
裏付ける、より正確な理解が得られた。corrhyp-goal.mdに詳細記録。
集計10/24で変わらず——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き15)
★実質的な前進★ Rレベル層の実現可能性を試掘した結果、mathlibの新しい
降下定理は不要と判明——既存部品だけで核心部品が完成した。

piece_descends_iso_promote: piece_descends_isoの降下先Rを、より粗い
共通段階R'(R≤R')へ昇格しても比較同型が保たれることを、
standardEtalePairPullbackIso(既存)を2回+exists_fg_subalgebra_tensor_
standardEtalePair_promote(既存)だけで証明できた——当初懸念した
「cofiltered極限からの降下定理を新規適用」は不要だった。

続けてpiece_descends_iso_R_of_proof(明示的証明版)・
piece_descends_iso_R_upperBound(exists_fgSubalgebra_upperBoundで
有限個のR_iを単一の共通上界R'へ合流)・
piece_descends_iso_R_upperBound_spec(R'が実際に各R_iの上界)を完成、
Rレベルの複数添字を単一の共通段階へ揃える核心部品が揃った。lake
build(ExtLimit/ABC3とも)0エラー確認、コミット: `3816c246`・
`e2e728af`。

残る道筋: 族全体を実際にR'レベルへ揃え、transitionElem/gdT/cocycle
一式のRレベル版の構築(層の制限写像の対応関係を精密に詰める必要あり、
未着手)。person-years評価をさらに上方修正できる可能性を示す実質的な
進展。corrhyp-goal.mdに詳細記録。集計10/24で変わらず——§4は引き続き
0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き16)
続報(piece_descends_iso_family_promote完成、有限個の候補片が単一の
R'へ実際に揃った、集計10/24で変わらず)。piece_descends_iso_R_of_proof
とpiece_descends_iso_promote内部の(piece_descends_iso...).chooseが
rflで一致することを確認した上で、各i∈tごとにpiece_descends_iso_
promoteを適用するだけで完成。**Rレベルの複数添字を単一の共通段階R'へ
合流させる部分が完成し、候補片Spec(P₀_i'.Ring)がすべて同じR'の上に
揃った**。lake build(ExtLimit/ABC3とも)0エラー確認、コミット:
`bf24a796`。

残る道筋: 揃った候補片同士の重なり(transitionElem/gdT/cocycleのR
レベル版)の構築——extDiagram Xの段階R'におけるアンビエントスキームの
層への対応を精密に詰める必要あり(未着手)。「候補片を単一の共通段階
へ揃える」という前半部分が既存部品だけで完成したのは、person-years
評価のさらなる見直しを裏付ける具体的な進捗。corrhyp-goal.mdに詳細
記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-04(続き17〜18)
アーキテクチャ再考: R_iを個別に合流させるより、Γ(C,piece)全体
(piece_algebraEtale_tensorによりA⊗ℝ上Etale=FinitePresentation)を
Algebra.FinitePresentation.iff_quotient_mvPolynomial'(mathlib、
MvPolynomial ι商としての表示)経由で**一度に**Rレベルへ降ろす方が
筋が良いかもしれないという代替案に気づいた——Cはextdiagramのような
R段階近似の塔を持たないため、個別のR_iを合流させるアプローチには
限界がある。

代替案を検証するため、exists_fg_subalgebra_tensor_bivariate_finset
(既存、2変数多項式の係数降下)のMvPolynomial ι版
exists_fg_subalgebra_tensor_mvPolynomial_finsetを実際に構築・検証
した——全く同じ技法(係数を有限個ずつ共通Rへ降ろしてmonomialの和で
再構成)が任意の変数集合ιへそのまま一般化できた。lake build
(FieldLimit/ABC3とも)0エラー確認、コミット: `8b2a983b`。

**これで代替案が仮説ではなく実際に動く部品を伴う具体的な作業計画に
なった**。次の一手: (ker f).FGから有限個の生成元を取り出し、
この補題を適用してΓ(C,piece)のRレベルモデルS_0を構成する
(未着手)。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き)
pieceAlgebra_relation_descend_R/_q₀/pieceAlgebra_R_model(S_0本体)を
構築(commit fe100836)。続けて**pieceAlgebra_R_model_baseChange**を
証明——S_0 ⊗[A⊗R.1] (A⊗ℝ) ≅ Γ(C,piece)。汎用補題
quotient_mvPolynomial_baseChange(FieldLimit.lean、(MvPolynomial ι R
⧸ I)⊗[R]A ≃+* MvPolynomial ι A ⧸ Ideal.map(MvPolynomial.map
(algebraMap R A))I)をAlgebra.TensorProduct.quotientTensorEquiv+
MvPolynomial.algebraTensorAlgEquivから構成し、CorrHyp設定へ
specialize。q₀の定義(choose_spec)+Algebra.Presentation.
span_range_relation_eq_ker+Algebra.Generators.ker_eq_ker_aeval_val+
RingHom.quotientKerEquivOfSurjective(第一同型定理)で閉じた。
配管注意: 素朴にはmaxHeartbeats既定値超過——既存の前例
(set_option maxHeartbeats 1000000 in)を踏襲して解消(新しい失敗形
ではない)。lake build(FieldLimit/ExtLimit/ABC3)0エラー確認、
commit 863a5f73、push済み。

**S_0が「Rレベルの候補片」として実際に使える裏付けが完成**。次の
一手: S_0をdescendPieceの代わりに使い、transitionElem/gdT/cocycle/
corrHypGlueDataのRレベル版を組み立てる(ℝレベルGlueDataパイプライン
はScheme一般の形で書かれているため再利用できる見込み)。集計は
10/24で変わらず——§4は引き続き0/2。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き2)
descendPieceR/descendPieceR_toBase/descendPieceR_reBaseMap/
descendPieceR_isoを構築(commit 65fd0a77)——S_0(pieceAlgebra_R_model)
をSpecに渡した正真正銘RレベルのスキームdescendPieceRが、ℝへbase
changeすると実際に元のα⁻¹(piece)に一致することを証明。pullbackSpecIso
(mathlib)+pieceAlgebra_R_model_baseChange+piecePreimageIsoの合成。

配管の新しい罠(lean-idioms.md #36として記録): 同じalgebraMap R S0の
生の式を型と証明の2箇所に書くと、S0のCommRing/Semiringインスタンスが
場所ごとに非defeqに導出されpullbackSpecIsoとの単一化に失敗する——letI
で先に固定しても`set`で名前を付けても直らず、algebraMapを使う射自体を
独立defとして1回だけ確定させる(standardEtalePairSpecMapと同じ
パターン)ことで解消した。

lake build(ExtLimit/ABC3)0エラー確認、push済み。**Γ(C,piece)をRレベル
へ一度に降ろす計画が環レベルだけでなくスキームレベルで完成**。次の
一手: 各片ごとに異なるR_iを持つ有限個のdescendPieceRを、Cでの重なり
方に沿ってRレベルで貼り合わせ、実際のD.Space元を構成する(遷移データ
がRレベルでliteralに一致することを示すdescent理論的な部分、未着手)。
集計は10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き3)
「項目(d)の第二段」(遷移データのRレベル降下)への具体的な道筋を発見。
mathlib `Mathlib/Algebra/Category/Ring/FinitePresentation.lean`の
`RingHom.EssFiniteType.exists_comp_map_eq_of_isColimit`(2つの射が
余極限で一致するなら共通の精密化で既に一致する)・`exists_eq_comp_
ι_app_of_isColimit`(余極限への射は有限段階からfactorする)が正確に
必要な形。土台の余極限(ℝ=colim R.1)は`FieldLimit.lean`の
`isColimitToRingCatCocone`として既に完成済み。残る配線: 余錐をA⊗[ℚ]-
でテンソルしてR↦A⊗R.1の余極限がA⊗ℝであることを示す一段(未着手)。
これが済めばS_0(Γ(C,piece)のRレベルモデル、finite presentation)へ
直接適用でき、隣接pieceの遷移同型が共通のR-level精密化で既に定義
されていることが示せる見込み——「650行規模の再構築」ではなく「既存
補題のspecialize」に近い可能性が出てきた、有望な見積もり修正。

lemma_4_1のh:ZK=D.Ext Z(文字通りの等号)構造的懸念は今回の調査でも
独立に再確認されたが、既に「拙速には着手しない」と判断済み
(isIsogenous_refl回帰リスク)であり判断を変える新情報は無い——遷移
データのRレベル降下という低リスクインフラ構築を優先。集計は10/24で
変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き4)
★大きな前進★ exists_mvPolynomial_quotient_ringHom_descendを証明
(commit 19121fe9)——「項目(d)の第二段」(遷移データのRレベル降下)の
核心が完成した。RingHom.EssFiniteTypeのfiltered colimit機構より
もっと直接的な道: ℚが体なのでA(ℚ-ベクトル空間)はModule.Flat——
R'.1↪ℝという単射をAでテンソルしても単射性が保たれる
(Module.Flat.lTensor_preserves_injective_linearMap)。これにより
「ℝレベルで等しいなら、その場のR'で既に等しい」という強い主張が
直接出る(commit 6e7c7734でalgebraTensorMap_val_injective系を構築)。

MvPolynomial.mapとaevalの可換性(mvPolynomial_map_aeval_comm、
mathlibに完成品が無かったので手で組んだ)を経て
mvPolynomial_aeval_eq_zero_of_map_aeval_eq_zero(commit 0e1b9b33)、
そしてexists_fgSubalgebra_upperBound2+algebraTensorMap_val_comp_
inclusionを組み合わせ、exists_mvPolynomial_quotient_ringHom_descend
を完成: ℝレベルで分かっている遷移写像の生成元の対応(有限個)から、
実際にR'レベルの候補写像を構成的に取り出せることを証明した。

transitionElem/gdT/cocycleのRレベル版に相当する遷移データの降下が、
650行規模の個別GlueDataエンジニアリング再構築ではなく、この1つの
補題へのPresentationデータのspecializeとして実現できる見込みが
立った——見積もりの大幅な改善。次の一手: 実際に2つの隣接する
descendPieceRの重なり部分へspecializeし、得られた候補写像が同型
であることを確認する(未着手)。lake build(FieldLimit/ABC3)0エラー
確認、push済み。集計は10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き5)
exists_mem_ideal_span_range_descendを完成(commit 32a3d1fc)——イデアル
所属もRレベルへ降りることの証明。「同型であることの確認」(往復合成が
恒等になること)に要る最後の核心部品。exists_mvPolynomial_quotient_
ringHom_descendと全く同じ「単射性で等式をR'レベルへ押し戻す」パターン
の別対象への適用——新しい数学的発想は不要だった。

残る配線(未着手): ev(片1→片2)・ev'(片2→片1)を共通R'''へ揃え、合成
ev'∘ev(MvPolynomial.comp_aevalで1段のaevalに潰せる)が恒等射と片1の
関係式イデアルを法として一致することをexists_mem_ideal_span_range_
descendで示す。3つの核心補題を組み合わせるだけの配線になった。

lake build(FieldLimit/ABC3)0エラー確認、push済み。正直な評価: 遷移
データの降下自体は前進しているが、これが完成してもGlueData構成・
β脚(文字通り未着手)・h:ZK=D.Ext Zの構造的懸念という3つの大きな課題が
残り、Lemma 4.1(§4)全体の完成にはまだ距離がある。集計は10/24で
変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き6)
★「項目(d)の第二段」の核心完成★ exists_mvPolynomial_quotient_
ringHom_descend2(commit c97d226e)——2つの異なるRレベル片(異なる
R・R₂)の間の遷移写像が共通の精密化R'へ構成的に降りることを完全に
証明した。exists_mvPolynomial_quotient_ringHom_descendの関係式条件
(生の等式=0)を「イデアル所属」へ一般化(exists_mem_ideal_span_
range_descend経由)し、目標側もquotient(genuineな片)である場合に
対応させた。

新規部品: mem_ideal_span_range_promote(所属の単調性)・
algebraTensorMap_inclusion_comp_inclusion(2段昇格=1段昇格、rflで
閉じる)・mvPolynomial_map_aeval_comm_general(完全一般版)・
exists_mvPolynomial_eval_descend(純粋な存在部分)を組み合わせ、
R・R₂・R₀(ψの存在部分が出す無関係なR)を2回のupperBound2で合流させ、
関係式ごと(有限個)に異なる精密化をさらにupperBoundで揃えた。

これでtransitionElem/gdT/cocycleのRレベル版に相当する遷移データの
降下——当初650行規模と見積もっていた作業が完成した。lake build
(FieldLimit/ABC3)0エラー確認、push済み。

正直な評価: 集計は10/24で変わらず。Lemma 4.1完成にはまだ(1)候補
写像が同型であることの確認(今回の枠組みの直接応用の見込み)・
(2)GlueDataとしての貼り合わせ配線(既存インフラの再利用が主体)・
(3)β脚(未着手)・(4)h:ZK=D.Ext Zの構造的懸念、という4課題が残るが、
遷移データ降下という核心は完成し、距離は着実に縮まっている。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き7)
★★「項目(d)の第二段」完全完成★★ exists_mvPolynomial_quotient_
ringEquiv_descend(commit d7753ea6)——遷移写像がRレベルで実際に
同型であることを証明。ψ・ψ'(ℝレベル、往復合成が恒等)から出発し、
共通の精密化R'上の候補ev・ev'が(a)ψ・ψ'を再現・(b)互いの関係式を
互いのイデアルへ写す・(c)往復合成が実際に恒等射、まで構成的に示した。

新規部品: round_trip_promote_eq/eq2(往復合成の昇格naturality)・
exists_round_trip_descend(ℝレベルの往復恒等性からRレベルの恒等性を
構成的に導く)・本体(descend2を両方向に適用しRcへ合流、round_trip_
descendを両方向に適用し最後にR'へ揃える)。証明本体は約380行、
Subalgebraの≤のproof irrelevanceを最大限活用(結合順序が違っても
defeq)——新しい配管の壁には当たらなかった。lake build(FieldLimit/
ABC3)0エラー確認、push済み。

**「項目(d)の第二段」全体(transitionElem/gdT/cocycleのRレベル版)が
完全に構成的に証明された**——650行規模と見積もっていた作業が汎用的な
数学的補題群として完成。正直な評価: 集計は10/24で変わらず。残る課題は
(1)GlueDataとしての貼り合わせ配線(今回の同型データが核心材料)・
(2)β脚(未着手)・(3)h:ZK=D.Ext Zの構造的懸念、の3つ。corrhyp-goal.md
に詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き8)
exists_mvPolynomial_quotient_ringEquiv_descend'(commit fedca21f)——
前回の生データ(ev・ev'がψ・ψ'を再現・関係式が互いのイデアルへ写る・
往復が恒等射)を実際の1個のRingEquiv項として組み立てる実用形が完成。
exists_mvPolynomial_quotient_ringEquiv_of_data(Ideal.Quotient.lift×2
+RingEquiv.ofRingHom、Ideal.Quotient.ringHom_ext+MvPolynomial.
ringHom_extで往復恒等を確認、CorrHyp非依存の汎用補題)を生データへ
適用するだけの配線。lake build(FieldLimit/ABC3)0エラー確認、push済み。

これで「Rレベルの候補片の遷移写像が実際に同型として得られる」という
Lemma 4.1のGlueData構成に必要な数学的核心のすべてが揃った(commit
19121fe9・32a3d1fc・d7753ea6・fedca21fの一連の成果)。残る課題:
(1)descendPieceR実データへの specialize+corrHypGlueDataへの配線
(既存部品の組み合わせになる見込み)・(2)β脚(未着手)・(3)h:ZK=D.Ext Z
の構造的懸念(既に判断済み)。集計は10/24で変わらず。corrhyp-goal.md
に詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き9)
exists_mvPolynomial_quotient_specIso_descend(commit 90c6342c)——前回
のRingEquiv(環レベル)をScheme.Spec.mapIsoで送るだけの配線で、
descendPieceR2つの間の実際のスキーム同型が共通の精密化R'上で構成的に
得られることを証明。「項目(d)の第二段」が環レベルからスキームレベル
まで完全に配線された。lake build(FieldLimit/ABC3)0エラー確認、
push済み。

正直な評価: 集計は10/24で変わらず。一連の成果(19121fe9・32a3d1fc・
d7753ea6・fedca21f・90c6342c)で数学的核心のすべてが揃った——道具は
完全に揃い、次のセッションでの実データ(ExtLimit.leanのpieceAlgebra_
relation_descend_R等)への specialize+corrHypGlueDataへの配線作業に
直接使える状態。残る課題(β脚・h:ZK=D.Ext Zの構造的懸念)は変わらず。
corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き10)
実データへのspecialize前に「何をψ・ψ'として渡すべきか」を精査。重要な
気づき: ℝレベルの貼り合わせは自明(α⁻¹(U_i)は全てc.C自身の開集合、
層構造が既に自動的に貼り合わせデータを与える)——Rレベルで初めて
非自明になる(descendPieceRはc.Cの開集合そのものではなくbase change
して初めて同型になる抽象的なスキーム)。実際のψの構成に要る3部品:
(1)P_i.val(Algebra.Generators.val、生成元の実際の環の元)・
(2)層の制限写像でoverlapへ送る・(3)P_jの生成元の多項式として表示
(aevalの全射性経由の選択)。hround1/2は層の制限公理から従う見込み。
配線作業(Scheme.Opens/TopCat.Presheaf API)がまだ残っている。集計は
10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き11)
piece_le_of_le/piece_restrict_hom(commit e5250e76)——c.C側の片同士の
制限写像を実装。同時に前回の道筋を訂正: overlap(U_i∩U_j)はU_j自体
より真に小さいのでP_j.valがoverlap上で生成系とは限らない——正しくは
共通のさらなる細分W(W≤U_i,U_j)を取り、両方をWへ制限してからW自身の
Presentationの生成元を仲立ちにψを構成する必要がある(「直接比較」
ではなく「共通の第三の片Wへ制限してから比較」という標準的な層の
貼り合わせの形)。exists_mvPolynomial_quotient_specIso_descend自体の
正しさには影響なし——配線側の設計ミスだった。次の一手: W:=U_i⊓U_j
がアフィンであることの確認+piece_restrict_homを両方向に適用する
配線(未着手)。lake build(ExtLimit/ABC3)0エラー確認、push済み。
集計は10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き12)
★重要な簡略化★ 「共通の細分W」問題はX.basicOpen(基本開)の被覆を
使えば分離性の仮定なしに解決する。IsAffineOpen.inf(mathlib)は
分離性([IsAffineHom (pullback.diagonal (terminal.from X))])を要求
するが、QcqsSpaceはCompactSpace・QuasiSeparatedSpaceのみで分離性を
含まない——安易に足すとNonempty Cと同種のリスク。

発見(朗報): U_i・U_jを任意のアフィン開ペアとして扱うのをやめ、単一の
アフィン開Uの基本開被覆X.basicOpen f_iで覆う設計に切り替えれば解決
する——Scheme.basicOpen_mul(mathlib、既存): 基本開同士の交わりは
再び基本開(同じ環の別の元による局所化)なので常にアフィン、分離性
不要。W := X.basicOpen(f_i*f_j)として具体的に取れる。これは実は
corrHypGlueData(既存)が最初から採用していた設計そのもの。

次の一手: descendPieceRをU:=X.basicOpen f_iの場合に適用し、
piece_restrict_hom(D(f_i f_j)≤D(f_i)・D(f_i f_j)≤D(f_j)の2方向)
で生成元をWへ制限、W自身のPresentationを仲立ちにψを構成する。
分離性問題を回避できたことで、既存のcorrHypGlueDataOfEtaleの設計
パターンをそのまま踏襲できる、より安全な道筋になった。集計は10/24
で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き13)
piece_restrict_hom_basicOpen_left/right(commit b8c6503d)——分離性
不要な共通細分Wの制限写像(D(f)→D(f*g)・D(g)→D(f*g)の両方向)を実装。
piece_restrict_homをScheme.basicOpen_mul+inf_le_left/rightへ
specializeするだけの短い配線。lake build(ExtLimit/ABC3)0エラー確認、
push済み。

これでψ構成に必要な「2つの片を共通の第三の片Wへ制限する」配線の
両方向が実際のRレベル環準同型として揃った。次の一手: W自身の
Presentationの生成元を仲立ちに、両側の生成元をaevalの全射性経由で
Wの生成元の多項式として表示し、実際のψ・hround1・hround2を構成する
(未着手)。集計は10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き14)
ψ構成の数学的筋道を確定(MvPolynomial.comp_aeval+ker_eq_ker_aeval_val+
span_range_relation_eq_kerで解決、hround1/2は不要——D(f)・D(g)を直接
比較せずWと2回比較すればよい)。Scheme.GlueData(mathlib)の構造も精査:
J・U・V・f(開埋め込み)・t・t'・t_fac・cocycleの約8データ。V(i,j):=
D(f_i*f_j)を対称定義すればt:=idで済む簡略化に気づいた。

残るギャップ(正直な記録): exists_mvPolynomial_quotient_specIso_
descendが与えるのは抽象的な同型だが、GlueDataが要求するのは開埋め込み
(IsOpenImmersion)——W=D(f_i*f_j)がD(f_i)の基本開であるというRレベル
での実現(S_0の局所化としての実現)がまだ未構成。新しい技術的ギャップ
として記録。集計は10/24で変わらず。corrhyp-goal.mdに詳細記録。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き16)
piece_basicOpen_mul_eq(commit 71807250、ExtLimit.lean)——piece(D(f*g))
がpiece(D(f))の基本開そのものであることをCレベルで確立。Scheme.
preimage_basicOpen(mathlib)をpullback.fst・αの2段の逆像へ押し出す
だけだが、rwを2回連鎖させると「motive is not type correct」((ExtF.obj
X).leftとpullback X.hom toBaseKの構文不一致)になったため、中間結果を
明示的に型注釈したhave+exact(defeq判定)の2段構成で解消——新しい
配管の失敗形。

正直な評価: これはCという1つのアンビエントスキーム内での事実であり、
続き15のギャップ本体(descendPieceRという独立に構成されたRレベルの
抽象スキーム同士を実際の局所化として結びつける、Algebra.Presentation
の選び方自体の制御が必要)はまだ未着手。lake build(ExtLimit/ABC3)
0エラー確認、push済み。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き17)
設計転換に気づいた——descendPieceRを「独立なPresentationの事後比較」
ではなく「D(f)側のS₀をhで局所化したものとして直接構成」する方針へ。
isLocalization_away_tensor_eq(commit 470da5b1、FieldLimit.lean):
S₀をhで局所化してからB上でTとテンソルするのと、先にB上でTとテンソル
してからhの像で局所化するのが一致することを証明(cancelBaseChange+
tensorRight+algEquivの合成、新しい数学は不要)。

この方針の利点: この構成の下ではSpec.map(algebraMap S₀ M)がIsOpen
Immersion.of_isLocalization(mathlib)により自動的に開埋め込みになる
——独立な2つのPresentationが局所化関係にあることを事後証明するという
困難な問題を、最初から局所化として構成することで回避できる。ただし
descendPieceR自体の再構成(ExtLimit.lean側)はまだ未着手。lake build
(FieldLimit/ABC3)0エラー確認、push済み。集計は引き続き10/24——§4は
引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き18)
Exploreエージェントで§4以外の残り11項目を横断調査——§4集中が正しい
優先付けであることを独立確認。全項目がmathlib丸ごと不在の理論
(代数群論・Galoisコホモロジー・被覆数有限性・モジュライスタック・
Teichmüller空間論)かLemma 4.1自身への依存(Lemma 5.1)で着手不可能。
唯一の意外な発見: Lemma 5.4(∃c>0,∀X,c≤e_Y)は一見有望だったが、現状の
StackType型(i:Sigma→ℕに≥2制約なし、g,rに下限なし)では原理的に
証明不可能と判明(g=r=0,Sigma=∅で反例)——interfaceに新公理を要する、
IsAffineOpen.infの分離性追加見送りと同種の判断が必要。集計は10/24。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05(続き19)
exists_fg_subalgebra_tensor_quotientMvPolynomial_lift(commit c7ec8748、
FieldLimit.lean)——Γ(C,piece)の任意の元をRレベルへ持ち上げる補題を
完成。quotient_mvPolynomial_baseChange+Ideal.Quotient.mk_surjective+
exists_fg_subalgebra_tensor_mvPolynomial_finset(既存技法の転用)。

新しい配管の発見(lean-idioms.md #40): MvPolynomial ι(テンソル積)⧸I の
HasQuotient自動探索が失敗する場合、letI hCR:CommRing(...)で先に確定
させるべきは係数環自身であってMvPolynomial本体ではない——descendPieceR
の実装がたまたま正しい方を選んでいたと判明。

isLocalization_away_tensor_eq+この補題で部品は揃ったが、descendPieceR
自体の再構成(W:=D(f*g)をD(f)側の環の局所化として定義し直す配線)は
まだ未着手。lake build(FieldLimit/ABC3)0エラー確認、push済み。集計は
引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜(自律ループ)
配線第一歩を試行(ファイルには何も書き込んでいない、REPL検証のみ)。
letI足場は`intro`ではなく`have`/`letI`で証明本体内に再構成する必要が
ある(型注釈内のletIはintroできない、ゼータ簡約済みになる、新発見)。
「存在するだけ」版(結論True)はmaxHeartbeats 1000000で通ったが、実際の
正しさの等式を結論にすると`e.symm(...)`の型検査だけでmaxHeartbeats
4000000でもタイムアウト——巨大な型を1つの型注釈に4〜5回重複させたのが
原因と見立て。次の道: R'・p₀をまず軽い「True」版から.chooseで名前付き
defとして取り出し、正しさの等式は別立ての小さなtheoremとして名前付き
定数を参照する形で1回だけ書く。時間の都合でここで打ち切り。集計は
引き続き10/24。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜続き(自律ループ2ティック目)
「名前付きdef分割」を試したが数学的な落とし穴が2つ判明(ファイル変更
なし)。(1) True結論のExistsから.chooseで取り出したp₀は、実際の正しさ
(hp₀)と紐付いていない(証明無関係性のため、捨てた時点で証人の正しさが
失われる)。(2) 名前付き定数を参照する形で正しさの等式を書けば型は軽い
(4.73秒)が、eを∀量化すると数学的に誤り(一般に環同型は一意でない、
証明できる特定のeについてしか成り立たない)。結論: e・R'・p₀・正しさの
等式は1つの証明の中で一括構成する必要がある(分離不可)——巨大な型の
重複回避にはstructureで束ねるか頻出部分式を独立defにする方向性が次の
候補。腰を据えた1セッションを要する局面と判断しここで打ち切り。集計は
引き続き10/24。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜続き(訂正・完成)
前回の「∀eは数学的に誤り」は誤った推論だった——訂正。実際の原因は
単純な配管バグ: 型注釈中の無名let(let n:=...・let I:=...)を証明側で
introせずにintro eを呼んでいたため、introが誤ってletの方を消費し
(e:ℕになる)、e.symmがNat.symmを探すという分かりにくいエラーになって
いた。intro n I eと3つとも明示的に消費して解消——maxHeartbeats
4000000(実測82秒)で無事通った。lean-idioms.md #41に記録。

exists_piece_basicOpen_R_lift(commit 05b66144、ExtLimit.lean)完成
——piece_basicOpen_localizationElemを任意のeに対してR'レベルの多項式
p₀として持ち上げられることを示した。lake build(ExtLimit/ABC3)0エラー
確認、push済み。次はdescendPieceR自体をD(f*g)についてD(f)側の局所化
として直接構成し直す配線(Spec(Localization.Away h)として書き、
IsOpenImmersion.of_isLocalizationで自動的に開埋め込みになることを
示す)——まだ未着手。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き
descendPieceR_localization_isOpenImmersion(commit 96f66baf、ExtLimit.
lean)完成——descendPieceRの環をR'レベルへ昇格したものを任意の元p₀で
局所化すればIsOpenImmersion.of_isLocalizationにより自動的に開埋め込み
になることを示した(f,gに依存しない一般的な事実)。配管の教訓:
open...in・set_option...inはdocstringの前に置く必要がある(後に置くと
「unexpected token 'open'」構文エラー)。

意義: exists_piece_basicOpen_R_liftが与える具体的なR'・p₀と組み合わせ
れば、GlueDataのf i jに要る開埋め込みが直接得られる——続き15で発見した
本丸のギャップへの直接の答えが揃った。残るのはScheme.GlueDataの完全な
構造(J・U・V・f・f_open・t・t'・t_fac・cocycleの約8データ)への組み立て
——特にtの構成にD(f)↔D(g)対称性のRレベル実現が要る。lake build
(ExtLimit/ABC3)0エラー確認、push済み。集計は引き続き10/24——§4は
引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き2
piece_isLocalization_basicOpen_mul(commit 326f50d1、ExtLimit.lean)完成
——descendPieceR_localization_isOpenImmersion(Rレベル)がℝレベルで
正しい対象を実現していることを確認するℝレベル側の土台。piece(D(f))の
アフィン性にIsAffineOpen.isLocalization_basicOpen(mathlib)を直接適用
するだけ。

残る橋渡し: descendPieceRの局所化のℝへの底変換が実際にΓ(C,piece
(D(f*g)))に一致することを示すには、isLocalization_away_tensor_eq・
cancelBaseChange・exists_piece_basicOpen_R_liftのhp₀・今回の補題という
4つの既存部品を1本の同型の鎖として繋ぐ必要があり、まだ未完成——複数の
同型合成を慎重に追跡する作業、次の一手。lake build(ExtLimit/ABC3)
0エラー確認、push済み。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き3
ideal_map_mvPolynomial_promote_baseChange_eq(commit d27f4b80、FieldLimit.
lean)完成——IをB→B'→Tと2段で底変換したものがB→Tと1段で底変換したもの
に一致することを示した(cancelBaseChange経由より単純な設計、quotient_
mvPolynomial_baseChangeを2回適用して比較する道の核心部品)。

正直な記録: この事実を含む完全な同型の鎖を、descendPieceRの全letI
足場を伴った1つの巨大な証明として組み立てようとしたが、maxHeartbeats
4000000(約85秒待った上で)でもタイムアウトした——letIの数が増えるほど
型クラス探索が遅くなる現象を観測(Etale/FinitePresentationのhaveIを
追加しただけで、それより前のAlgebra instanceのletI自体が新たに
タイムアウトするように)。次の道: e1・e2・hIdealEq・hp₀を個別のhave
として先に確立してから引数として渡す小さな最終補題に分割する(足場の
再構築を1回で済ませる)。lake build(FieldLimit/ABC3)0エラー確認、
push済み。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き4
分割方針を実行——elaboration timeoutは完全解消、数学的論証も完成した
がinstance diamondで最後の1手に阻まれる。quotient_mvPolynomial_
baseChange_tmul_one(commit 901f7362、FieldLimit.lean、「純テンソル」
上での自然性)を新たに確立し、ideal_map_mvPolynomial_promote_
baseChange_eqと組み合わせてkeyという等式(局所化パラメータの正しい
対応)を完全に証明できた。

最後にIsLocalization.ringEquivOfRingEquiv(e2から局所化同士の同型を
作る)を呼ぶとinstance diamond(Algebra.TensorProduct.instMul対
instDistribOfSemiring.toMulの不一致)に当たり未完成——letIでの事前
登録を複数回試したが解消せず。次の道: IsLocalization.algEquiv経由で
構成し直す、またはisLocalization_away_tensor_eqと同じスタイルで一般の
Mの側でIsLocalization.Awayインスタンスを直接構成する。lake build
(FieldLimit/ABC3)0エラー確認、push済み。集計は引き続き10/24——§4は
引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き5(完成)
instance diamondを解消——exists_ringEquiv_localization_of_eq(commit
6aaaa280、FieldLimit.lean)が完全に証明できた。鍵: 最小反例を切り出して
原因特定——具体的なテンソル積型のままIsLocalization.ringEquivOfRingEquiv
を呼ぶと必ずdiamondに当たるが、A・Bを先に抽象的なCommRing型として
一般化した補題(ringEquiv_localization_of_apply_eq)を別立てで用意して
から具体的な型を代入すればdiamondは起きないと判明(letIでの事前登録
ではなく「抽象化してから代入」という一段違うレベルの教訓)。

意義: descendPieceRのR'レベル局所化が実際にℝレベルで正しい対象を
実現することが数学的にもLeanの証明としても完成した。descendPieceR_
localization_isOpenImmersion(開埋め込み性)と合わせれば、GlueDataの
f i jに要る全データが揃った——続き15で発見した本丸のギャップがついに
解決した。残るのはScheme.GlueDataの完全な構造への配線(特にtの構成に
D(f)↔D(g)対称性が要る)。lake build(FieldLimit/ABC3)0エラー確認、
push済み。集計は引き続き10/24——§4は引き続き0/2だがLemma 4.1の数学的
内容はほぼ組み上がった。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き6
exists_ringEquiv_localization_of_eqを実際のX,U,hU,f,g,C,αへ適用して
一体化しようとしたところ再びelaboration timeout(maxHeartbeats
4000000で123秒)——足場のletIが増えるほど遅くなる現象が再現。

解決: exists_ringEquiv_of_piece_lift(commit 8c437650、FieldLimit.lean)
としてさらに一段抽象化——Ctarget・Wtarget・e・h₂まですべて抽象的な型
変数として扱い、pieceAlgebra等を一切参照しない形にしたところ2.18秒で
通った。elaboration timeoutは巨大な証明そのものではなく、CorrHyp固有
の型を証明の中に持ち込むこと自体が原因だったと判明。

意義: Lemma 4.1のGlueData構築に必要な「Rレベルの局所化とℝレベルの
正しい対象を結びつける」核心部分の代数的な部品はすべて揃った。呼び
出し側はe・h₂・hp₀・WtargetのIsLocalization.Awayインスタンスを揃えて
渡すだけで適用できる形になったが、「揃えて渡す」配線自体(pieceAlgebra
等の足場を伴う)はまだelaboration timeoutで未完成——次はこの配線側も
さらに細かく分割する。lake build(FieldLimit/ABC3)0エラー確認、
push済み。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き7(完成)
壁の正体を再検証——分割の仕方の問題ではなく、単にmaxHeartbeatsが
足りていなかっただけだった。4000000ではタイムアウトしていた同じ配線
を40000000(10倍)に上げてバックグラウンドで気長に待ったところ227秒
で無事通った。

exists_descendPieceR_localization_baseChange(commit bc6db389、
ExtLimit.lean)完成——exists_piece_basicOpen_R_lift・piece_
isLocalization_basicOpen_mul・exists_ringEquiv_of_piece_lift・
piece_basicOpen_mul_eqの4部品を1つの証明に組み立てただけ、新しい
数学的内容は不要だった。

意義(大きな前進): Lemma 4.1のGlueData構築における「Rレベルの局所化と
正しいℝレベルの対象を結びつける」核心部分が完全に完成した。
descendPieceR_localization_isOpenImmersion(開埋め込み性)と合わせれば
GlueDataのf i jに要る全データが揃った——続き15で発見して以来セッション
全体を通じて追いかけてきた本丸のギャップが解決した。

残る作業: Scheme.GlueDataの完全な構造への組み立て(J・U・V・f・
f_open・t・t'・t_fac・cocycle)——特にtの構成にD(f)↔D(g)対称性の実現、
複数添字にまたがる共通段階R''への合流が要る。β脚(丸ごと未着手)も
残る。lake build(ExtLimit/ABC3)0エラー確認、push済み。集計は引き続き
10/24——§4は引き続き0/2だがLemma 4.1の核心的な数学的困難は解決した。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き8
tの構成を精査——ψ構成が本当に必要で避けて通れないことを確認(正直な
記録、実装はまだ)。exists_descendPieceR_localization_baseChangeを
(f,g)・(g,f)両方に適用すればM_ij⊗T≃Γ(C,piece(D(f*g)))≃M_ji⊗Tという
ℝレベルの比較は直ちに得られるが、GlueDataのtはRレベルのスキーム射で
なければならず不十分——M_ij≃M_ji(Rレベル)という独立に構成された
2つのRレベル対象同士を結びつける、今回とは別の問題。

必要な道具(exists_mvPolynomial_quotient_specIso_descend、このセッション
より前から存在)は既にあるが、ψ・ψ'の構成にはWの独立したAlgebra.
Presentationを新たに構成する必要があり(descendPieceR自体の構成と
同規模の足場)、今回の「局所化として直接構成する」設計とは別の独立した
作業になる。拙速な近道は無いことを確認した上で、次のセッションへの
明確な出発点として記録。集計は引き続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き9
ψ構成の設計を1つ具体化——exists_finite_standardEtaleCoverからの新規
構成ではなく、既存の「持ち上げ」技法(exists_fg_subalgebra_tensor_
quotientMvPolynomial_lift)の転用でよいと判明(実装はまだ)。M_ij・M_ji
はどちらもΓ(C,piece(D(f*g)))へのℝレベル同型を既に持つ——ψはM_jiの
生成元(有限個)をこの同型経由でM_ij側へ持ち上げるだけで構成できる、
Wの独立したPresentationを新規構成する必要はない。

残る作業: M_ij・M_ji自体をMvPolynomial(添字⊕Unit)(底)⧸(拡張イデアル)
という関係式の族として具体的に表す必要がある——Localization.
awayEquivAdjoin(Polynomial/AdjoinRoot形)は既存だがMvPolynomial(Fin n
⊕Unit)形への書き換えがまだ要る、具体的だが軽い次の一歩。集計は引き
続き10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き10
localization_away_quotient_mvPolynomial_equiv(commit f0247a4b、
FieldLimit.lean)完成——tの構成の第一歩を実装した。当初検討していた
awayEquivAdjoin経由ではなく、IsLocalization.Away.mvPolynomialQuotient
Equiv+MvPolynomial.quotientEquivQuotientMvPolynomialをIdeal.
quotientEquivで繋ぐ、より直接的な道で完成させた。CorrHyp非依存の
一般的な可換環論の事実。

残る作業: これは入れ子の多項式環による商の形——exists_mvPolynomial_
quotient_specIso_descendが要求する1段のMvPolynomial(n⊕Unit)B商の形へは
MvPolynomial.sumAlgEquivでさらに1段変換する作業が残る。lake build
(FieldLimit/ABC3)0エラー確認、push済み。集計は引き続き10/24——§4は
引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き11
localization_away_quotient_mvPolynomial_flat_equiv(commit afbb09e4、
FieldLimit.lean)完成——入れ子の商をMvPolynomial(Unit⊕n)Bの1段の商へ
平坦化した。MvPolynomial.sumAlgEquiv+DoubleQuot.quotQuotEquivQuotSup
(側条件I≤J不要)を組み合わせ、生成元の対応はsumAlgEquiv_comp_rename_
inl/inr(naturality)から具体的に計算(e2.symm(C x)=rename Sum.inr x、
e2.symm(X())=X(Sum.inl()))。I=Ideal.span(Set.range q₀)(有限生成)の
場合にexists_mvPolynomial_quotient_specIso_descendが要求するq・q₂の
形そのものが得られることを確認した。「R↔ℝブリッジ」に続く「多項式
表示の平坦化」問題が完成。lake build(FieldLimit/ABC3)0エラー確認、
push済み。集計は引き続き10/24——§4は引き続き0/2。

次の一手: flat_equivをD(f)側・D(g)側(M_ij・M_ji)双方に具体的に
インスタンス化しq・q₂として渡す。ψ・ψ'は既存の「lift to R-level」
技法(exists_fg_subalgebra_tensor_quotientMvPolynomial_lift)をM_jiの
個々の生成元へ繰り返し適用する方針(続き9で確立済み)で進める。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き12
localization_away_quotient_mvPolynomial_flat_equiv_of_map(commit
99071c6e)+isLocalization_of_ringEquiv_transport(commit02eb897c、
FieldLimit.lean)完成。前者はflat_equivをalgebraMap昇格イデアルの形
へ橋渡し、後者は「移送先のAlgebra構造を環同型eそのものとして定義
すれば両立条件がrflで落ちる」という一般的な移送パターン(lean-idioms
#49として記録)。

正直な行き詰まり: D(f)側への具体的インスタンス化(exists_
descendPieceR_flat_mvPolynomial_baseChange)を試みたが、主張の型を
書く段階でletI不足による停留変数エラーに繰り返し当たり未完成。
force しない。次はdescendPieceR等の既存パターンに倣い、平坦化した
イデアルを独立defとして先に切り出す設計に変える。集計は引き続き
10/24——§4は引き続き0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き13
exists_descendPieceR_flat_mvPolynomial_baseChangeの主張の型は完成
(停留変数の真因は `rename Sum.inr p₀ * X (Sum.inl ()) - 1` のHMul
探索でSum.inrの左成分Unitが決まらないこと——両因子に型注釈で解消)。
証明骨格も抽象型変数の最小例では動く。しかしlake build(1回15〜20分)
を3回まわして3回とも別の場所で落ちた:(1)M_flatのSemiringが2経路に
割れる→letI hCRMで解消、(2)同じ割れがQ側→letI hCRQで解消、
(3)IsScalarTower B' Q Mの型がSubmodule.Quotient.instSMul'で表示され
of_algebraMap_eqのAlgebra.toSMul 3本組と合わない→letIでは勝てない。

根本原因: 商環MvPolynomial ι B'⧸Jは係数環B'に対する自前のSMulを
持っており、SMul探索でそちらが勝つ。環同型e'越しに移送したAlgebra
B' Mとは一致しない。つまり「≃+*を作ってからAlgebraを後付け移送」
という作戦自体が誤り(lean-idioms #52として記録——当初#51としたが
並行セッションと衝突したため改番)。

正しい方針: FieldLimit.lean の localization_away_quotient_
mvPolynomial_equiv → flat_equiv → flat_equiv_of_map の3本を
≃ₐ[B'](AlgEquiv)として作り直す。mathlibにAlgEquiv版の部品が全部
揃っていることは確認済み(mvPolynomialQuotientEquiv・
quotientEquivQuotientMvPolynomial・quotQuotEquivQuotSupₐ・
Ideal.quotientEquivAlg・sumAlgEquiv、基底の変更はAlgEquiv.
restrictScalars)。FieldLimit側は1本2〜3秒で反復できるので、
ExtLimitの15〜20分ビルドを回すより遥かに速い。

リポジトリ状態: 通らない定理はgit checkoutで差し戻し済み、
lake build ABC3は0エラー(6590 jobs)確認。書きかけはscratchpadの
wip-flat-baseChange.patchに保存。集計は引き続き10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き14
方針転換を実行。FieldLimit.leanに≃ₐ[B]版3本を追加(sorry無し、
lake build ABC3 0エラー6590 jobs確認):
- localization_away_quotient_mvPolynomial_algEquiv(Away局所化の任意の
  実現SとMvPolynomial Unit(MvPolynomial n B)商のB-代数同型)
- localization_away_quotient_mvPolynomial_flat_algEquiv(1段の
  MvPolynomial(Unit⊕n)B商への平坦化)
- localization_away_quotient_mvPolynomial_flat_algEquiv_of_eq
  (イデアルを変数Iqで受けIq = span(range q₁)をsubstする形。呼び出し側の
   インスタンスがIdeal.map(map φ)I₀の形でもrwと食い違わない)

配管の罠(記録): Ideal.map (mkₐ B I') JのIsTwoSidedは係数環が入れ子の
MvPolynomialだと自動で見つからない——CommRingを2本先に登録すると通る
(#40と同型)。抽象的な(R : Type)[CommRing R]では起きないので最小例では
見逃す種類の罠。

次の一手: ExtLimit.lean側で(1)eを正準なM₀ := Localization.Away(mk I' p₀)
で実体化(Algebra B' M₀とIsScalarTower B' Q M₀が正準に存在することは
REPLで確認済み)、(2)flat_algEquiv_of_eqでM₀ ≃ₐ[B'] Fを得る、
(3)Algebra.TensorProduct.congrでF⊗[B']T ≃ M₀⊗[B']Tへ移して合成。

並行セッションとの衝突2件: (a)FieldLimit.leanへの追加はgit add直後に
並行セッションのコミット7675ebd0に巻き込まれた(内容は正しく入っている)。
(b)lean-idiomsの#51が衝突したのでこちらを#52へ改番。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き15
exists_descendPieceR_flat_mvPolynomial_baseChange完成(sorry無し、
commit 47b2f5c8、ExtLimit.lean)。Γ(C,piece(D(f*g)))が1段の
MvPolynomial(Unit⊕Fin n)(A⊗R'.1)商のℝ底変換そのものとして書けた:
q = Sum.elim (fun k => rename Sum.inr (map φ (q₀ k)))
             (fun _ => rename Sum.inr p₀ * X (Sum.inl ()) - 1)。

設計3手: (1)eを正準なM₀ := Localization.Away (mk I' p₀)で実体化
——Algebra B' M₀もIsScalarTower B' Q M₀も正準に存在するので
インスタンスを自作しない(続き13の失敗との決定的な違い)、
(2)flat_algEquiv_of_eqでM₀ ≃ₐ[B'] F、(3)Algebra.TensorProduct.congrで
F⊗[B']T ≃ M₀⊗[B']Tへ移して合成。

配管: maxHeartbeats 40000000に加えsynthInstance.maxHeartbeats(既定
20000)も4000000へ上げる必要。IsScalarTower B' Q M₀は小さい文脈で
先にhaveIで計算(46秒)。REPL検査584秒。

重要な発見(lean-idioms #53): lake build ABC3はFound/CorrHyp/を
ビルドしない(lakefile.tomlにglobが無くABC3/Found.leanにCorrHypの
行が無い)。CorrHypの検証は lake build ABC3.Found.CorrHyp.Instance4
を明示的に叩くこと。失敗ビルド直後はExtLimit.oleanが消えたままでも
lake build ABC3は0エラーを返し、その状態のREPLは「空の環境」になる
(unknown namespace AlgebraicGeometry / Γ(X,U)のparse error)。

次の一手: D(g)側は同じ定理を(g,f)で使うだけ(D(g*f)=D(f*g)はmul_comm)。
両側が揃ったらψ・ψ'をexists_fg_subalgebra_tensor_quotientMvPolynomial_
liftで構成しexists_mvPolynomial_quotient_specIso_descendへ渡す。
集計は引き続き10/24——§4は0/2だがtのD(f)側という最大の部品が完成。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き16
設計上の食い違いを発見: exists_mvPolynomial_quotient_specIso_descend
(ψ・ψ'の機械)は単一のAを要求するが、我々の2つの表示はD(f)側が
A_f := Γ(X.left,D(f))、D(g)側がA_g := Γ(X.left,D(g))という別々の底環の
上にある。ψ'の向きの環準同型を作るにはA_g →(A_f側の商)が要るが無い。

解消の見立て: A₃ := Γ(X.left,D(f*g)) = A_f[1/g]。D(f)側の平坦化表示F₁は
局所化関係式 rename Sum.inr p₀ * X(Sum.inl ()) = 1 を明示的に持つので
p₀の像は単元。p₀はh₂(=gの像)の持ち上げ(hp₀)なので「g⊗1の像がF₁で
単元」が言えればF₁はA₃⊗R'-代数になり両側が共通の底A₃に乗る。ただし
hp₀の一致はℝ底変換後なので有限段階へ降ろす必要がある。

その道具を完成(commit 5c642cca、sorry無し、Instance4ビルド0エラー):
- mvPolynomial_map_val_inclusion_comp(map val ∘ map inclusion = map val)
- exists_fgSubalgebra_mvPolynomial_ideal_mem_descend
  (ℝレベルでイデアルに属せばRを大きくしてR'レベルでも属する。
   A⊗ℚℝ = colim_R(A⊗ℚR.1)のフィルター余極限性の実用形。証明4手:
   mem_span_range_iff_exists_fun → 係数をR₁へ降ろす → R⊔R₁ →
   A⊗R.1→A⊗ℝの単射性)

この補題はψ・ψ'の4仮説(hψ/hψ'/hround1/hround2、いずれも「ℝレベルで
イデアルに属する」形)の検証にもそのまま使える。tの構成の両方の難所に
効く道具。集計は引き続き10/24——§4は0/2。

次の一手: 「g⊗1の像がF₁で単元」をhp₀+今回の降下補題で証明する。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き17
isUnit_rename_of_flat_relation完成(commit bc903fff、sorry無し)——
平坦化した商Fの中では局所化の分母pの像が単元(関係式rename Sum.inr p *
X(Sum.inl()) - 1からX(Sum.inl())のクラスが逆元)。

残る半分と障害: 欲しいのは「algebraMap (A_f⊗R') F (g⊗1)が単元」。繋ぐには
C(g⊗1)とrename Sum.inr p₀のクラスがFで一致することが要る。hp₀はそれを
ℝ底変換後に与えるので続き16の降下補題で降ろせる——はずだが障害が1つ:
h₂(gの像)がeを通してC(g⊗1)のクラスに対応することを使うには、e(=
pieceAlgebra_R_model_baseChange)がΓ(X.left,U)⊗ℝ-代数同型でなければ
ならないが、現状の結論はNonempty(… ≃+* …)という環同型である。中身
(quotient_mvPolynomial_baseChange+第一同型定理)はどちらも底環上の代数
写像なのでAlgEquivへ格上げできるはず——続き13〜14でflat_equivを≃ₐへ
作り直したのと同種の作業。

次の一手: pieceAlgebra_R_model_baseChangeのAlgEquiv版を作る。
集計は引き続き10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き18
pieceAlgebra_R_model_baseChange_const完成(commit d21bfed3、ExtLimit.lean、
sorry無し、Instance4ビルド0エラー)。続き17の障害(eが≃+*で値の情報が無い)を
AlgEquivへの全面作り直しではなく「必要な情報だけ併記する∃形」で解消:

∃ e : (…) ≃+* Γ(C,piece), ∀ b, e (mk (C b) ⊗ₜ 1) = algebraMap … b

もとの証明が組み立てている同型(quotient_mvPolynomial_baseChange →
Ideal.quotEquivOfEq → 第一同型定理)をそのまま返すだけでよかった。値の計算は
quotient_mvPolynomial_baseChange_tmul_one(自然性、今セッション前半で作成)+
Ideal.quotEquivOfEq_mk・RingHom.quotientKerEquivOfSurjective_apply_mk+
MvPolynomial.map_C+AlgHom.commutes。

配管: RingHom.ker (aeval P.val)(coe版)と RingHom.ker (aeval P.val).toRingHom
はrflで等しいが構文的に違う——rwで扱うには後者へ揃える。

揃ったもの: (a)isUnit_rename_of_flat_relation=「Fでp₀の像は単元」、
(b)今回=「C bのクラスの行き先はalgebraMap b」。あとはhp₀と(b)を突き合わせて
「C(g⊗1)とp₀のクラスがℝレベルで一致」を出し、続き16の降下補題で有限段階へ
降ろせば、Fの中でg⊗1の像が単元→FがA₃⊗R''-代数→両側が共通の底に乗る。
集計は引き続き10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き19
設計の転換に気づいた(重要)。GlueDataのV (i,j)を「U_iの片を局所化したもの」
ではなく「U_i ⊓ U_jの片そのもの」(descendPieceR X (U_i ⊓ U_j) …)に取ると:

- V (i,j)とV (j,i)はU_i ⊓ U_j = U_j ⊓ U_i(inf_comm)より同じ対象。t i jは
  eqToHom(実質恒等)で済み、t_id・cocycleも自明に近い。ψ・ψ'の機械が
  丸ごと不要になる。
- f i j : V (i,j) ⟶ U iは「表示から表示への写像」であり、同型ではなく
  写像を降ろすだけ。(a)生成元の行き先をexists_fg_subalgebra_tensor_
  quotientMvPolynomial_lift(既存)で降ろし、(b)関係式が消えることを
  続き16のexists_fgSubalgebra_mvPolynomial_ideal_mem_descendで降ろす。

残る本質的作業: f i jが開埋め込みであること。U_i ⊓ U_jがU_iの基本開の場合は
IsOpenImmersion.of_isLocalization+descendPieceR_localization_isOpenImmersion
の筋で押さえられる見込み。

続き16〜18の成果は無駄にならない: 降下補題は新設計の(b)の中心部品、
isUnit_rename_of_flat_relationとpieceAlgebra_R_model_baseChange_constは
「局所化としての開埋め込み」の議論で引き続き効く。続き15のflat表示も
f i jの生成元の行き先を書き下す場面で使える。

次の一手: 新設計のf i j(写像の降下)を構成する。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き20
(1)訂正: 続き16のexists_fgSubalgebra_mvPolynomial_ideal_mem_descendは
FieldLimit.leanに既にあったexists_mem_ideal_span_range_descendと完全重複
だった(既にdescend2の証明中で使われていた)。配管の
mvPolynomial_map_val_inclusion_compも薄い包みで使う先が重複補題だけ
だったので両方削除(commit a22c3636)。CLAUDE.mdの「在庫」どおり書く前に
ファイル自身の名前一覧を引くべきだった。lean-idioms #56に在庫一覧と
引き方を記録(#54は並行セッションと衝突したため#56へ)。

(2)新設計に要る道具はすでにある: f i j(写像の降下)は既存の
exists_mvPolynomial_quotient_ringHom_descend2そのもの。唯一の差は底環が
A=Γ(X.left,U_i)からA'=Γ(X.left,U_i⊓U_j)へ動くことだが、関係式を先に
Algebra.TensorProduct.map φ (AlgHom.id ℚ R.1)でA'側へ押し出してから
descend2をA := A'で使えば無料で一般化できる(根拠:
(map (id A')(val R))∘(map φ (id R)) = map φ (val R)、REPLで1行確認済み)。

次の一手: descend2の「底が動く版」を薄いラッパとして書き、descendPieceRの
Uについての関手性(V ⊆ Uに沿った制限写像のRレベル版)を構成する。
集計は引き続き10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き21
exists_mvPolynomial_quotient_ringHom_descend2_of_map完成(commit f0d4f0fb、
FieldLimit.lean、sorry無し、Instance4ビルド0エラー)。descend2の「底環が
動く版」で、証明は関係式を先にφで押し出して既存descend2へ帰着するだけ(15行)。

これで新設計(続き19)のf i jを作る材料が揃った: 生成元の行き先ψはℝレベルの
制限写像Γ(C,piece(U_i))→Γ(C,piece(U_i⊓U_j))から取り、関係式が落ちることは
ℝレベルで自明、この補題で有限段階R'へ降ろすとR'レベルの環準同型が得られ、
Specを取ればf i j : V (i,j) ⟶ U iになる。

次の一手: この補題へdescendPieceRのデータ(pieceAlgebra_relation_descend_q₀)
をspecializeしてf i jを構成する。その後に残るのは開埋め込み性(U_i⊓U_jが
U_iの基本開の場合は局所化の議論で押さえられる見込み)。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き22
f i jの構成に進み、必要な仮説hψを分解したところ、最後に必ず
「pieceRingEquivのUについての自然性」(制限写像とφ⊗idの四角形が可換)が
要ると判明。理由: hψをP_V.aeval_valで潰すには「aeval(制限∘P_U.val)∘
map(φ⊗val)」と「制限∘P_U.aeval_val」の一致が要り、変数の行き先は一致
するが係数の行き先の一致がまさにその四角形。

これは3つ目の「不透明な同型の値が要る」場面(続き17のh₂とC(g⊗1)、
続き18でpieceAlgebra_R_model_baseChange_constとして解決、今回)。
pieceRingEquivはpiecePullbackIso(6段のcalc)をΓで送った同型で、値に
ついての補題が1つも無いのが根本原因。

次の一手(2案、後者が有望): (a)piecePullbackIsoの6段それぞれの可換性、
(b)pieceRingEquiv.symmが標準写像(Algebra.TensorProduct.liftで作る
pullback.fst.appとℝ側構造射の組)に一致することを示す——一致すれば
自然性はScheme.Hom.appの自然性から自動。mathlibに
pullbackSpecIso_inv_fst/_hom_sndなど特徴づけの補題が揃っていることを確認。

この1本は続き17のg⊗1にもf i jにも将来のβ脚にも効く「一度払えば何度も
使える」補題。集計は引き続き10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き23
自然性の四角形を2つに分解し標準写像の側を完成(commit c8a3ad5e、sorry無し、
Instance4ビルド0エラー):
(i)piece_appLE_naturality(完成)= pullback.fstのappLEが開集合の制限と可換。
   mathlibのScheme.Hom.map_appLE+appLE_mapだけで2行。
(ii)pieceRingEquivと標準写像の一致(未完)= pieceRingEquiv.symmがa⊗ₜ1を
   appLE aへ送ること。piecePullbackIso(6段のcalc)に値の補題が無いのが原因。
   mathlibのpullbackSpecIso_inv_fst/_hom_sndで特徴づける見込み。

分解自体をExtLimit.leanの節見出しに書いたので、次に触る者は四角形のどちら側が
残っているかを読むだけで分かる。次の一手は(ii)。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き24
(ii)の追跡に着手し、piecePullbackIsoの定義の形が障害だと判明。
目標: (piecePullbackIso).inv ≫ (piece U).ι ≫ pullback.fst
     = Spec.map (ofHom includeLeftRingHom) ≫ hU.isoSpec.inv ≫ U.ι
mathlibに各段の特徴づけ補題が揃っていることは確認
(pullbackRestrictIsoRestrict_inv_fst・pullbackSymmetry_hom_comp_fst・
pullbackRightPullbackFstIso_inv_fst・pullbackSpecIso_inv_fstなど)。

しかし#printで見るとpiecePullbackIsoの6段はTrans.transの左結合の入れ子で、
4段目が (pieceRingHom_spec) ▸ Iso.refl _ という等式輸送。showで書き直すと
invalid ▸ notation, failed to compute motiveになり、unfold+simpも進まない。

直し方(特定済み): 4段目をpullback.congrHom (pieceRingHom_spec X U hU) rfl
で書き直す。mathlibにpullback.congrHom_inv(= pullback.map … 𝟙 𝟙 𝟙 …)の
特徴づけがあるので合成の値が計算できる。残作業は「定義を▸抜きへ書き換え」
+「6段の追跡」の2段構え。書き換えは型が同じなので既存利用箇所に影響しない。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き25
続き24の直し方を実際に検証(REPLでpiecePullbackIso'として複製):
(1)置き換えは成功——4段目をpullback.congrHom (pieceRingHom_spec).symm rflに
替えた定義は型検査を通り、unfold+simp only [Iso.trans_inv, Iso.symm_inv,
Category.assoc]で合成が完全に露出する(▸版ではshowがmotive計算に失敗して
到達できなかった状態)。
(2)しかし各段の書き換えが発火しない。原因は(ExtF.obj X).leftと
pullback X.hom toBaseKの構文的ずれと見ている(ExtF := Over.pullback toBaseK ⋙
Over.map toBaseKなので定義的には等しいがsimpは展開しない)。

次の一手: 1.piecePullbackIsoの4段目をcongrHomへ置換(型不変なので影響なし)、
2.(ExtF.obj X).leftをpullback X.hom toBaseKへ畳む(showかExtF_obj_leftの
simp補題)、3.pullbackRestrictIsoRestrict_hom_ι_assoc→pullbackSymmetry_inv_
comp_fst_assoc→pullbackRightPullbackFstIso_inv_snd_fst→pullback.congrHom_inv
+lift_fst→pullbackHomIsoLeft_hom_fst'のinv版→pullbackSpecIso_inv_fstの順。
各段の補題の存在は確認済み。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き26
piecePullbackIso_inv_fst完成(commit 1d457ad2、sorry無し、Instance4ビルド
0エラー確認):
 (piecePullbackIso).inv ≫ (piece U).ι ≫ pullback.fst
   = Spec.map (ofHom includeLeftRingHom) ≫ hU.isoSpec.inv ≫ U.ι
pieceRingEquivはpiecePullbackIsoをΓで送っただけなので値の情報はここからしか
取れない。続き22で特定した「3つ目の不透明な同型」の核心。

片付いた3つの詰まり: (1)(ExtF.obj X).leftという型注釈を書くとsimpの照合が
効かない(注釈を外す)、(2)定義がcalc+▸だとunfoldしても合成が露出しない
(.transの明示的な鎖+pullback.congrHomへ書き換え、型は不変)、(3)最後の2段は
rwが使えない——HasPullbackインスタンスが定義由来と探索由来で食い違い表示が
同一なのにパターンが見つからないと言われる(congrArg+Eq.trans+exactで解決)。

手順の反省: 1度目のコミット57c3fa9dはビルド確認前にコミットして失敗し
51327005で差し戻した(grepの先頭3行だけ見てEXIT:0と早合点)。2度目は
grep -c "^error"が0かつEXIT:0を確認してからコミットした。

次の一手: この等式からpieceRingEquiv.symm (a ⊗ₜ 1) = appLE a(環レベル)を
導き、続き23の(i)と合わせて四角形を閉じる。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き27
続き26のSpec側の等式から環レベルの pieceRingEquiv.symm (a⊗ₜ1) = appLE a へ
渡す道筋を4手に分解。mathlibの橋渡し2本を発見:
- IsAffineOpen.isoSpec_inv_ι : hU.isoSpec.inv ≫ U.ι = hU.fromSpec
- IsAffineOpen.SpecMap_appLE_fromSpec :
  Spec.map (f.appLE U V i) ≫ hU.fromSpec = hV.fromSpec ≫ f  ←appLEのSpec側の意味

4手: 1.続き26の等式をfromSpecの形へ書き直す、2.SpecMap_appLE_fromSpecで
appLE側を得る、3.hU.fromSpecがモノ(開埋め込み∘同型)なので右から消去、
4.piecePullbackIso.inv ≫ (piece).isoSpec.hom = Spec.map (pieceRingEquiv.symm)
の同定(pieceRingEquivの定義=topIso+Γ.mapIso+ΓSpecIsoを追う)。
ここまで来ればSpecの忠実充満性で環レベルの等式が出る。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き28
piecePullbackIso_inv_isoSpec_appLE完成(commit b1db7975、sorry無し、
Instance4ビルドのerror行0件・EXIT:0確認):
 (piecePullbackIso).inv ≫ (piece).isoSpec.hom ≫ Spec.map appLE
   = Spec.map (ofHom includeLeftRingHom)
続き27の4手のうち1〜3が完了(続き26の等式+SpecMap_appLE_fromSpec+
hU.fromSpecがモノであることによる右消去)。

残るのは4手目のみ: piecePullbackIso.inv ≫ (piece).isoSpec.hom を
Spec.map (pieceRingEquiv.symm) と同定する(pieceRingEquivの定義=
topIso+Γ.mapIso piecePullbackIso.symm.op+ΓSpecIsoを追う)。そこまで行けば
Specの忠実充満性から環レベルの pieceRingEquiv.symm (a⊗ₜ1) = appLE a が出て、
続き23の(i)と合わせて自然性の四角形が閉じる。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き29
4手目(piecePullbackIso.inv ≫ (piece).isoSpec.hom = Spec.map pieceRingEquiv)に
着手。数学的道筋は完全に判明し材料も揃った:
- IsAffineOpen.isoSpec hU = (↑U).isoSpec ≪≫ Spec.mapIso U.topIso.symm.op(定義)
  → pieceRingEquivのe1(topIso)に対応
- Scheme.isoSpec_hom_naturality (f := piecePullbackIso.inv) で左辺が
  (Spec(A⊗ℝ)).isoSpec.hom ≫ Spec.map (inv.appTop) に化ける → e2に対応
- Scheme.isoSpec_Spec_hom : (Spec R).isoSpec.hom = Spec.map (ΓSpecIso R).hom
  → e3に対応
Spec.mapの反変性で順序が逆になることまで含め辻褄が合う。

止まっている理由: rw [← Category.assoc]ですら「パターンが見つからない」と
言われる(instances transparencyの警告つき)——続き26で
piecePullbackIso_inv_fstを通したときと同じ壁。直し方も同じで、rwを諦めて
congrArg+Eq.transで項として組み立て最後をexactにすればよい見込み。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き30
自然性の四角形が閉じた(commit 751791e5、sorry無し、Instance4ビルドerror 0件・
EXIT:0、sorryはSkeleton側7件のみ=Found側0件):
  pieceRingEquiv_appLE : appLE ≫ pieceRingEquiv = includeLeftRingHom
つまりpieceRingEquivはappLE aを純テンソルa ⊗ₜ 1へ送る。続き23の
piece_appLE_naturality((i))と合わせて四角形が閉じた。

前段: piecePullbackIso_inv_isoSpec_hom(piecePullbackIsoとisoSpecの差が
ちょうどpieceRingEquiv)。pieceRingEquivの3成分がtopIso↔IsAffineOpen.isoSpecの
定義、Γ.mapIso↔Scheme.isoSpec_hom_naturality、ΓSpecIso↔isoSpec_Spec_homに
対応する。環レベルへはSpec.map_injectiveで落とす。

続き22〜30の教訓: rwが効かない箇所は一貫してcongrArg+Eq.trans+exactで項として
組み立てれば通る(HasPullback等のインスタンスが定義由来と探索由来で食い違い、
表示が同一でもrwはパターンを見つけられない。instances transparencyの警告が出る)。

次の一手: 揃った材料(四角形+descend2_of_map)で新設計のf i jを実際に構成する。
t i jはeqToHomで済むのでGlueDataの主要部が組める。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き31
piecePullbackIso_inv_snd完成(commit 6ba403da、sorry無し、Instance4ビルド
error 0件・EXIT:0、Found側sorry 0件):
  (piecePullbackIso).inv ≫ (piece U).ι ≫ pullback.snd = Spec.map (ofHom includeRight)
続き30のfst側(includeLeft)の相方。テンソル積からの環準同型は a⊗ₜ1 と 1⊗ₜr の
値で決まるので、pieceRingEquiv.symmの自然性にはこの2成分が要る。両方揃った。

証明はfst側と同じ6段で、最後だけpullbackSpecIso_inv_snd・
pullbackRightPullbackFstIso_inv_snd_snd・pullbackHomIsoLeft_inv_snd'(今回追加)
に差し替えるだけだった。

次の一手: ℝ側も環レベルへ落とし(Spec.map_injective+SpecMap_appLE_fromSpec、
⊤のfromSpecはIsAffineOpen.fromSpec_top)、2成分からpieceRingEquiv.symmの
自然性を組み立てる。そこまで行けば新設計のf i jをdescend2_of_mapで降ろす前提が
すべて揃う。集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き32
pieceRingEquiv_appLE_snd + piecePullbackIso_inv_isoSpec_appLE_snd完成
(commit b38bc8ca、sorry無し、Instance4ビルドerror 0件・EXIT:0、Found側sorry 0件)。
pieceRingEquivの値が a⊗ₜ1側(続き30)と 1⊗ₜr側(今回)の両方で確定した。
テンソル積からの環準同型はこの2成分で決まるので、pieceRingEquiv.symmのUに
ついての自然性はこの2本+piece_appLE_naturality(続き23)から組み立てられる。

配管: cancel_monoはrwではなく(cancel_mono _).mpで使う(specKとSpec(of ℝ)の
構文的ずれでrwのパターン照合が失敗)。specK.isoSpec.invの書き換えもhaveで型を
明示してからrw。

次の一手: 2成分からpieceRingEquiv.symmの自然性を組み立てる
(Algebra.TensorProduct.ext系で2成分から環準同型の一致を出す)。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き33
Scheme.Hom.appLE_preimage_naturality完成(commit c19bbab0、sorry無し、
Instance4ビルドerror 0件・EXIT:0、Found側sorry 0件): 任意のスキームの射
f : C ⟶ Y について f.appLE が開集合の制限と可換。piece_appLE_naturalityは
この特殊化に書き換えた(証明1行)。新設計のf i jではα(C ⟶ Ext X)側でも同じ
自然性が要る(Γ(C,piece(U))の代数構造がα.appLE ∘ pieceRingEquiv.symmのため)。

材料一覧(自然性まわり): appLE_preimage_naturality(任意の射)、
piece_appLE_naturality(pullback.fst版)、pieceRingEquiv_appLE(a⊗ₜ1成分)、
pieceRingEquiv_appLE_snd(1⊗ₜr成分)、piecePullbackIso_inv_fst/_snd・
piecePullbackIso_inv_isoSpec_hom(Spec側の特徴づけ)。

次の一手: 2成分からpieceRingEquiv.symmの自然性(環準同型の一致)を組み立て、
α側の自然性と合わせて「Γ(C,piece(U))の代数構造が制限と可換」を出す。そこまで
行けばdescend2_of_mapのhψが検証でき新設計のf i jが構成できる。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き34
pieceRingEquiv_symm_naturality完成(commit 95f44bc1、sorry無し、Instance4
ビルドerror 0件・EXIT:0、Found側sorry 0件):
  restr_ExtX (pieceRingEquiv_U.symm (a ⊗ₜ r)) = pieceRingEquiv_V.symm ((restr_X a) ⊗ₜ r)
証明は a⊗ₜr = (a⊗ₜ1)(1⊗ₜr) と分解し、pieceRingEquiv_symm_tmul_one・
pieceRingEquiv_symm_one_tmul(今回併せて追加、続き30・32のCommRingCat射の等式を
元レベルへ落としたもの)を当て、appLE側の自然性で移すだけ。
Algebra.TensorProduct.extのような重い道具は不要だった。

続き22で「3つ目の不透明な同型」として立ちはだかった障害は、Spec側・環レベル・
両成分・自然性のすべてで解消した。

次の一手: Γ(C,piece(U))のA_U⊗ℝ-代数構造はα.appLE ∘ pieceRingEquiv.symmなので、
今回の自然性とScheme.Hom.appLE_preimage_naturality(α側)を合わせて「代数構造が
制限と可換」を出し、descend2_of_mapのhψの検証に使ってf i jを構成する。
集計は10/24——§4は0/2。

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★2026-09-05夜さらに続き35
pieceAlgebraMap_naturality完成(commit aa08ea53、sorry無し、Instance4ビルド
error 0件・EXIT:0、Found側sorry 0件):
  restr_C (α.appLE_U (pieceRingEquiv_U.symm (a ⊗ₜ r)))
    = α.appLE_V (pieceRingEquiv_V.symm ((restr_X a) ⊗ₜ r))
Γ(C,piece(U))の代数構造は定義からα.appLE ∘ pieceRingEquiv.symmなので、
続き34の自然性と続き33のα側の自然性を合成するだけ(証明はcongrArg 2つ)。

続き22〜35の総括: 「f i jを降ろすには代数構造が制限と可換であることが要る」→
「pieceRingEquivの自然性が要る」→「その値が要る」→「piecePullbackIso(6段の
calc)の値が要る」と遡り、定義の書き換え(▸→pullback.congrHom)から始めて
Spec側の特徴づけ(fst・snd)→環レベル→純テンソル2成分→自然性→代数構造の
可換性、と14段階を積み上げて解消した。

次の一手: この可換性でdescend2_of_mapのhψを実際に検証し、f i j(R'レベルの
環準同型のSpec)を構成する。集計は10/24——§4は0/2。

## 2026-09-05夜さらに続き36 — hψ検証の計算核 eval₂_map_aeval_eq_zero

FieldLimit.lean に純粋計算の補題を追加(error 0件・EXIT:0、Found側sorry 0件):
  eval₂_map_aeval_eq_zero
    (hcomm : restr.comp algU = algV.comp φ')
    (hψval : ∀ i, eval₂ algV valV (ψ i) = restr (valU i))
    (hp : eval₂ algU valU p = 0) :
    eval₂ algV valV (eval₂ C ψ (map φ' p)) = 0
= descend2_of_map の hψ そのもの。hcomm には続き35の
pieceAlgebraMap_naturality がそのまま入る。証明は eval₂_comp_left と
eval₂_map を各2回だけ。

詰まり2つ:
(1) rw [eval₂_comp_left ...] がパターン不一致。ゴールが素の eval₂ f g x、
    補題は (eval₂Hom f g) x の適用形。→ show で頭を eval₂Hom に揃える。
(2) (fun i => restr (valU i)) vs (⇑restr ∘ valU) の食い違い。
    → have の型の方を ⇑restr ∘ valU で書き funext hψval で作る。

次の一手: この補題に pieceAlgebraMap_naturality を hcomm として食わせ、
descend2_of_map を呼んで f i j : V (i,j) ⟶ U i を構成する。
集計は10/24——§4は0/2。

## 2026-09-05夜さらに続き37 — hψ を descend2_of_map にそのまま食わせる形へ

相棒2本を FieldLimit.lean に追加(error 0件・EXIT:0、Found側sorry 0件):
- mem_ideal_of_eval₂_eq_zero: 表示 e : MvPolynomial ι' 𝔹 ⧸ J ≃+* T があれば
  eval₂ algV valV = e ∘ mk(ringHom_ext で C b と X i だけ見る)。よって
  eval₂ = 0 ⇒ mk r = 0 ⇒ r ∈ J。一発通過(0.84秒)。
- aeval_map_mem_ideal_of_relation: 続き36と合成し、hψ の形そのもの
  (aeval ψ (map φ' p) ∈ J)。aeval_def + algebraMap_eq の2書き換えだけ。

次の一手: descendPieceR のデータに halg/hval/hψval を具体的に用意し
(pieceRingEquiv の値は続き33〜34で計算済み)、descend2_of_map を呼んで
f i j : V (i,j) ⟶ U i を得る。集計は10/24——§4は0/2。

## 2026-09-05夜さらに続き38 — hψ が hcomm だけで自動になった

mathlib の Algebra.Generators は切断 σ を持つ(aeval_val_σ)。よって
ψ i := PV.σ (restr (PU.val i)) と置けば hψval は aeval_val_σ そのもの。
FieldLimit.lean に2本追加(error 0件・EXIT:0、Found側sorry 0件、両方一発通過):
- aeval_map_relation_mem_ker (Generators 版、結論は PV.ker)
- aeval_map_relation_mem_span (Presentation 版、結論は span (range relation))
  ← descend2_of_map の hψ の形そのもの。

意味: §4 の f i j 降下に要る hψ は hcomm = pieceAlgebraMap_naturality
(続き35)だけから出る。続き22〜35で積んだ14段階がここで効いた。

次の一手: descendPieceR の R' レベルのデータで Presentation の relation と
descend2_of_map の q/q₂ を突き合わせ、f i j を得る。集計は10/24——§4は0/2。
