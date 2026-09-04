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

end ABC3.Found.CorrHyp
