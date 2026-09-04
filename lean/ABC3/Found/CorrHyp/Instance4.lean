/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.QcqsSpace
import ABC3.Found.CorrHyp.ExtLimit
import ABC3.Skeleton.CorrHyp.Section4

/-!
# [CorrHyp] `Space := QcqsSpace`(qcqs スキーム)による `HyperbolicCurveData` の具体化

`corrHypInstance3`(`Instance3.lean`、`Space := Over BaseK`、一般のスキームを
許す)は `Lemma 4.1` の証明に要る `CompactSpace`/`QuasiSeparatedSpace` の
前提を満たせなかった(`ExtLimit.lean`/`corrhyp-goal.md` §4 に記録)。
`corrHypInstance4` はこれを `QcqsSpace`(`QcqsSpace.lean`)に絞ることで
解決した——双曲曲線は常に有限型(qcqs)なので原文の意図には忠実。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Interface.CorrHyp CategoryTheory AlgebraicGeometry

/-- `BaseK`(`Spec ℚ`)はアフィン。 -/
theorem BaseK_isAffine : IsAffine BaseK := by
  show IsAffine (Spec (CommRingCat.of ℚ)); infer_instance

/-- `QcqsSpace` の適当な基点(`ModuliStack` 等、`Lemma 4.1` が読まない
フィールドの placeholder 用)——`BaseK` 自身がアフィンなので qcqs。 -/
noncomputable def basePt4 : QcqsSpace :=
  haveI := BaseK_isAffine
  ⟨Over.mk (𝟙 BaseK), by show CompactSpace BaseK; infer_instance,
    by show QuasiSeparatedSpace BaseK; infer_instance⟩

/-- **`corrHypInstance4`**——`Space := QcqsSpace`(`k := ℚ`-スキームのうち
qcqs なもの全体)による `HyperbolicCurveData` の具体化。`FEt`・`idFEt`・
`comp`・`pullback`・`pbFst`・`pbSnd`・`Ext`・`extFEt` は `QcqsSpace.lean`
の本物のデータ、それ以外(`Lemma 4.1` の主張が読まない)は placeholder。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrHypInstance4 : HyperbolicCurveData where
  Space := QcqsSpace
  FEt := QcqsFEt
  idFEt := QcqsFEt.idFEt
  comp := QcqsFEt.comp
  pullback := QcqsFEt.pullback
  pbFst := QcqsFEt.pbFst
  pbSnd := QcqsFEt.pbSnd
  Iso X Y := X = Y
  type _ := (0, 0)
  Fuchsian := Unit
  Gamma _ := ()
  IsDiscrete _ := True
  Gamma_isDiscrete _ := trivial
  Comm := id
  FiniteIndexIn _ _ := True
  Sub _ _ := True
  MargulisArithmetic _ := False
  ShimuraArithmetic _ := False
  IsConnected _ := True
  core := id
  coreMap X := QcqsFEt.idFEt X
  Aut _ := PUnit
  idAut _ := PUnit.unit
  Ext := QcqsExt
  extFEt f := QcqsExtFEt f
  stackType _ := ⟨0, 0, Empty, inferInstance, Empty.elim⟩
  ModuliStack _ _ := basePt4
  IsGenericallyScheme _ := True
  deg _ := 0
  IsOpenDenseIn _ _ _ := True
  IsConstructibleIn _ _ _ := True

open ABC3.Skeleton.CorrHyp in
/-- **`Lemma 4.1` の statement は `corrHypInstance4` に対して意味を持つ**
(型検査が通ることの記録)——本体はまだ `sorry`。`corrHypInstance3` との
違いは、ここでは `X : QcqsSpace` なので `X.2.1`/`X.2.2` から
`CompactSpace`/`QuasiSeparatedSpace X.1.left` が常に取り出せること
——`AffineTransitionLimit.lean` の被覆補題群の前提を満たす。 -/
theorem lemma_4_1_at_instance4_statement :
    (∀ (X ZK : corrHypInstance4.Space) (_c : Corr corrHypInstance4
      (corrHypInstance4.Ext X) ZK), True) := fun _ _ _ => trivial

/-! ## ロードマップ項目(c')——`corrPieceGlueData`を`Corr`の実データへ代入する

`ExtLimit.lean`の`corrPieceGlueData`(`X:Over BaseK`・`C:Scheme`・
`α:C⟶(ExtF.obj X).left`・`[IsFinite α][Etale α]`からGlueDataを直接
与える)を、`corrHypInstance4`の下での`Corr corrHypInstance4 (QcqsExt X) ZK`
の実データ(`c.C`・`c.α`)へ実際に代入する。`QcqsFEt A B := FEtK A.1 B.1`
なので`c.α : QcqsFEt c.C (QcqsExt X)`は`{f : c.C.1 ⟶ (QcqsExt X).1 //
IsFinite f.left ∧ Etale f.left}`の要素——`c.α.1.left`・`c.α.2.1`・
`c.α.2.2`がそのまま`corrPieceGlueData`の`α`・`[IsFinite α]`・`[Etale α]`
に対応する。

★配管の罠(新発見、`lean-idioms.md`#31「instances透明度の壁」の新しい
現れ方): `c.α.1.left`を型注釈なしでそのまま`corrPieceGlueData`へ渡すと、
`[IsFinite α][Etale α]`の型クラス探索が`(QcqsExt X).1.left`と
`(ExtF.obj X.1).left`(定義的には等しいが`QcqsExt`は`@[reducible]`でない
`def`)を`instances`透明度で見分けられず失敗する
(`failed to synthesize instance ... IsFinite (Over.Hom.left ↑c.α)`、
たとえ`haveI := c.α.2.1`が直前に効いていても)。**修正**:
`letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left`という
**明示的な型注釈付きの`letI`**で`hα`の型を構文的に望む形へ固定して
から使う——これで`instances`透明度でも一致が見える。 -/

open scoped TensorProduct in
open ABC3.Skeleton.CorrHyp CategoryTheory.Limits in
/-- **`corrPieceGlueData`を`Corr corrHypInstance4 (QcqsExt X) ZK`の実データ
(`c.α`)へ接続する**: `X.1.left`のアフィン開`U`ごとに、`c.C`側で対応する
`c.α⁻¹(piece)`のGlueDataを直接与える——`Lemma 4.1`の`c'.C`構成に向けた、
`corrHypInstance4`・`QcqsSpace`・`Corr`の実データへの初めての接続。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrPieceGlueDataOfCorr (X ZK : QcqsSpace)
    (c : Corr corrHypInstance4 (QcqsExt X) ZK)
    (U : X.1.left.Opens) (hU : IsAffineOpen U) : Scheme.GlueData :=
  letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left
  letI : IsFinite hα := c.α.2.1
  letI : Etale hα := c.α.2.2
  corrPieceGlueData X.1 U hU c.C.1.left hα

open scoped TensorProduct in
open ABC3.Skeleton.CorrHyp CategoryTheory.Limits in
/-- `corrPieceGlueDataOfCorr`が使う族は、実際に`c.α⁻¹(piece)`(`c.C`側で
`X.1.left`のアフィン開`U`に対応する片)を覆う——`corrPieceGlueData_cover`
をそのまま呼ぶだけ。 -/
theorem corrPieceGlueDataOfCorr_cover (X ZK : QcqsSpace)
    (c : Corr corrHypInstance4 (QcqsExt X) ZK)
    (U : X.1.left.Opens) (hU : IsAffineOpen U) :
    letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left
    letI : IsFinite hα := c.α.2.1
    letI : Etale hα := c.α.2.2
    letI := pieceAlgebra X.1 U hU
    letI : Algebra (Γ(X.1.left, U) ⊗[ℚ] ℝ) Γ(c.C.1.left, hα ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U)) :=
      ((Scheme.Hom.appLE hα (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U)
        (hα ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U)) le_rfl).hom.comp
        (pieceRingEquiv X.1 U hU).symm.toRingHom).toAlgebra
    haveI : Algebra.Etale (Γ(X.1.left, U) ⊗[ℚ] ℝ) Γ(c.C.1.left, hα ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U)) :=
      piece_algebraEtale_tensor X.1 U hU c.C.1.left hα
    letI h := exists_finite_standardEtaleCover (Γ(X.1.left, U) ⊗[ℚ] ℝ)
      Γ(c.C.1.left, hα ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U))
    (⨆ i ∈ h.choose_spec.choose,
      c.C.1.left.basicOpen (h.choose_spec.choose_spec.choose i))
      = hα ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U) := by
  letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left
  letI : IsFinite hα := c.α.2.1
  letI : Etale hα := c.α.2.2
  exact corrPieceGlueData_cover X.1 U hU c.C.1.left hα

open scoped TensorProduct in
open ABC3.Skeleton.CorrHyp CategoryTheory.Limits in
/-- **`corrPieceGlueDataOfCorr(...).glued`が`c.α⁻¹(piece)`(`c.C`側で
`X.1.left`のアフィン開`U`に対応する片)に実際に同型である**——
`corrPieceGlueData_glued_iso`(`ExtLimit.lean`、既に一般に確立済み)を
`corrPieceGlueDataOfCorr`自身の定義(`corrPieceGlueData`をそのまま
呼ぶだけ)からそのまま特殊化する。これで`Corr`の実データ(`c.C`・
`c.α`)から具体的に構成した候補片の族が、`X.1.left`の**任意の1つの
アフィン開**について、実際に`c.C`の対応する部分を再構成することが
証明できたことになる——`Lemma 4.1`の構成的降下のうち「単一アフィン片
レベル」の全工程(項目(a)〜(c')・項目(b)の特殊化)がこれで完結した。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def corrPieceGlueDataOfCorr_glued_iso (X ZK : QcqsSpace)
    (c : Corr corrHypInstance4 (QcqsExt X) ZK)
    (U : X.1.left.Opens) (hU : IsAffineOpen U) :
    (corrPieceGlueDataOfCorr X ZK c U hU).glued ≅
      (c.α.1.left ⁻¹ᵁ (pullback.fst X.1.hom toBaseK ⁻¹ᵁ U) : Scheme) := by
  letI hα : c.C.1.left ⟶ (ExtF.obj X.1).left := c.α.1.left
  letI : IsFinite hα := c.α.2.1
  letI : Etale hα := c.α.2.2
  show (corrPieceGlueData X.1 U hU c.C.1.left hα).glued ≅ _
  exact corrPieceGlueData_glued_iso X.1 U hU c.C.1.left hα

/- ★★次の一手(未着手、2026-09-04に理解を訂正): 「`c.C.1.left`の外側の
貼り合わせ」ではない——`c.α⁻¹(Ext V_k)`は`c.C`自身の中に既に開集合
として存在する(トートロジー、貼り合わせ不要)。真に必要なのは、
`descendPiece`(`corrHypGlueData`の`Z i`)がまだ`ℝ`レベルのまま
(`piece_descends_iso`が与える`Spec(P₀.Ring)`という**`R`レベル**の
候補片が`descendPiece`内部で使い捨てられ、外へ取り出されていない)
ことを踏まえ、この`R`レベルの候補片自体を`R`レベル(`FgSubalgebra
ℚ ℝ`の圏)で貼り合わせる新しい層——異なる添字・異なるアフィン片
ごとに異なりうる`R_i`を共通精密化`R_i⊔R_j`へ持ち上げて比較する新しい
議論が要り、`ℝ`レベルで構築した`transitionElem`/`gdT`/`cocycle`一式
に相当する困難をもう一段繰り返すことを意味する(新しい独立した規模
の数学的内容、既存部品の配線だけでは済まない)。

さらに(e)`α・β`脚と整合性の等式`h▸extCorr D c'=c`の構成——ただし
`h:ZK=D.Ext Z`という文字通りの命題的等号は、`Corr`定義の`Nonempty C`
欠落および`QcqsSpace`が同型類の商でないことと組み合わさると証明不可能
(あるいは偽)になりうるという構造的懸念が判明している(`corrhyp-
goal.md`2026-09-04の該当エントリに詳細記録、拙速な`Corr`修正は既に
完成済みの§1 5/5を後退させかねないため見送り中)。 -/

end ABC3.Found.CorrHyp
