/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.ShimuraArithmeticData

/-!
# [CorrHyp] `corrHypInstance` の `ShimuraArithmetic` を本物に差し替えた版

`Instance.lean` の `corrHypInstance` は `ShimuraArithmetic := fun _ ↦ False`
という placeholder(`Proposition 3.2` が読まないので無害、という理由づけ
だった)を使っていた。`ShimuraArithmeticData.lean` で `Definition 2.3` の
非空虚性(`shimuraArithmeticWitness_SL2Z`)が確立できたので、ここで
`ShimuraArithmetic` を **本物の述語**(`ShimuraArithmeticWitness`)に
差し替えた具体的な項 `corrHypInstance2` を作り、`Definition 2.3` の
`.src` を正当に主張できる状態にする。

★`corrHypInstance` を直接書き換えず(すでに `.src` を持つ
`corr_witness`・`isIsogenous_witness`・`prop_3_2_at_instance` 等が依拠する
ため)、別の項として並べる——`Instance.lean` はそのまま、この項は
`ShimuraArithmetic` だけが異なる**新しい具体化**。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Interface.CorrHyp

/-- `corrHypInstance` と(`ShimuraArithmetic` 以外)同じ、`HyperbolicCurveData`
の具体化——`ShimuraArithmetic` だけを `ShimuraArithmeticWitness`(本物の
述語)に差し替えた。 -/
noncomputable def corrHypInstance2 : HyperbolicCurveData where
  Space := FuchsianGroup
  FEt C X := PLift (IsFiniteIndexIn C X)
  idFEt X := PLift.up (isFiniteIndexIn_refl X)
  comp f g := PLift.up (isFiniteIndexIn_comp f.down g.down)
  pullback _f _g := inter _ _
  pbFst f g := PLift.up (isFiniteIndexIn_pbFst f.down g.down)
  pbSnd f g := PLift.up (isFiniteIndexIn_pbSnd f.down g.down)
  Iso X Y := X = Y
  type _ := (0, 0)
  Fuchsian := Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ)
  Gamma X := X.Γ
  IsDiscrete Γ := DiscreteTopology Γ
  Gamma_isDiscrete X := X.discrete
  Comm Γ := Subgroup.Commensurable.commensurator Γ
  FiniteIndexIn Γ1 Γ2 := Γ1 ≤ Γ2 ∧ (Γ1.subgroupOf Γ2).FiniteIndex
  Sub Γ1 Γ2 := Γ1 ≤ Γ2
  MargulisArithmetic _ := False
  ShimuraArithmetic Γ := ShimuraArithmeticWitness Γ
  IsConnected _ := True
  core := id
  coreMap X := PLift.up (isFiniteIndexIn_refl X)
  Aut _ := PUnit
  idAut _ := PUnit.unit
  Ext := id
  extFEt f := f
  stackType _ := ⟨0, 0, Empty, inferInstance, Empty.elim⟩
  ModuliStack _ _ := ⟨⊥, inferInstance⟩
  IsGenericallyScheme _ := True
  deg _ := 0
  IsOpenDenseIn _ _ _ := True
  IsConstructibleIn _ _ _ := True

open ABC3.Skeleton.CorrHyp in
/-- **安全性の確認**——`ShimuraArithmetic` の値を `False` から本物の述語に
差し替えても `Proposition 3.2` は壊れない(証明は arithmeticity の値を
一切使わないので、`corrHypInstance` のときと文字通り同じ証明が通る)。
`.src` は付けない(`Proposition 3.2` 自体はすでに `Instance.lean` の
`prop_3_2_at_instance` で登記済み——二重登記を避ける)。

★**sorry 無し**。標準3公理のみ。 -/
theorem prop_3_2_at_instance2 (X Z : corrHypInstance2.Space)
    (_hX : ¬ Arithmetic corrHypInstance2 (corrHypInstance2.Gamma X))
    (c : Corr corrHypInstance2 X Z) (_hC : corrHypInstance2.IsConnected c.C) :
    corrHypInstance2.Sub (corrHypInstance2.Gamma Z)
      (corrHypInstance2.Comm (corrHypInstance2.Gamma X)) :=
  corrhyp_prop_3_2 c.α.down c.β.down

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 2.3`(Shimura-arithmetic)の `corrHypInstance2`
における実現**——`Skeleton.CorrHyp.ShimuraArithmetic corrHypInstance2 Γ_SL2Z`
は定義から `corrHypInstance2.ShimuraArithmetic Γ_SL2Z = ShimuraArithmeticWitness
Γ_SL2Z` に展開され、`shimuraArithmeticWitness_SL2Z` がそのまま証明になる。

★**sorry 無し**。標準3公理のみ。 -/
theorem shimuraArithmetic_at_instance2 : ShimuraArithmetic corrHypInstance2 Γ_SL2Z :=
  shimuraArithmeticWitness_SL2Z

def shimuraArithmetic_at_instance2.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 5, item := "Definition 2.3", sectionId := "corrhyp-def-2-3" }

end ABC3.Found.CorrHyp
