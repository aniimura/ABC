/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.TensorIndex
import ABC3.Meta.Claim

/-!
# 可逆加群でテンソルしても全射性は落ちない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か——底変換の道具

底変換 `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の**有限側**を出すには

    `μ_L : 𝓞_K ⊗_{𝓞_F} Γ(L) → Γ(f^*L)`（切断の引き戻し）

が同型であることが要る。★台帳 `presheaf-pullback-global-sections` に書いたとおり、
mathlib の前層引き戻しは随伴の左随伴としてしか定義されていないので**直接には計算できない**。

★★そこで**テンソルで回す**:

    `μ_{A ⊗ B} ≅ μ_A ⊗ μ_B`（在庫の `delta_pullSec`、`§9-743`）
    `μ_{Ō} ` は同型（在庫の `pullbackUnitPreIso`、`§9-745`）
    `A ⊗ A⁻¹ ≅ Ō`（`AInv.isInv`）

なので `μ_A ⊗ μ_{A⁻¹}` は同型、したがって**全射**である。
★★★本ファイルはそこから `μ_A` が全射であることを出す道具を用意する
——あとは `Module.Invertible.bijective_of_surjective` で同型に上がる。

## ★★★★★機構

    `map α β = rTensor B' α ∘ lTensor A β`

なので `map α β` が全射なら `rTensor B' α` が全射。
★そのとき `coker α ⊗ B' = 0` であり、`B'` が可逆なら `- ⊗ B'` は忠実なので `coker α = 0`。
★★忠実性は `C ≅ (C ⊗ B') ⊗ Bᵛ` という同型で出る（`Module.Invertible.linearEquiv`）。
-/

namespace ABC3.Found.Arakelov

open scoped TensorProduct

/-- ★★★★★★**可逆加群でテンソルすると自明性が戻る**（`- ⊗ B'` は忠実）。

★`C ≅ C ⊗ R ≅ C ⊗ (Bᵛ ⊗ B') ≅ (C ⊗ B') ⊗ Bᵛ` という同型で出る。 -/
theorem subsingleton_of_tensor_invertible (R C B' : Type) [CommRing R] [AddCommGroup C] [Module R C]
    [AddCommGroup B'] [Module R B'] [Module.Invertible R B'] (h : Subsingleton (C ⊗[R] B')) :
    Subsingleton C := by
  have e : C ≃ₗ[R] (C ⊗[R] B') ⊗[R] Module.Dual R B' :=
    (TensorProduct.rid R C).symm ≪≫ₗ
      (TensorProduct.congr (LinearEquiv.refl R C) (Module.Invertible.linearEquiv R B').symm) ≪≫ₗ
      (TensorProduct.congr (LinearEquiv.refl R C) (TensorProduct.comm R _ _)) ≪≫ₗ
      (TensorProduct.assoc R C B' (Module.Dual R B')).symm
  haveI : Subsingleton ((C ⊗[R] B') ⊗[R] Module.Dual R B') := by infer_instance
  exact e.injective.subsingleton

/-- ★★★★★★★**`α ⊗ id_{B'}` が全射なら `α` が全射**（`B'` は可逆）。

★`coker α ⊗ B' = 0` を出し、`- ⊗ B'` の忠実性で `coker α = 0` にする。 -/
theorem surjective_of_rTensor_surjective (R A A' B' : Type) [CommRing R]
    [AddCommGroup A] [Module R A] [AddCommGroup A'] [Module R A']
    [AddCommGroup B'] [Module R B'] [Module.Invertible R B']
    (α : A →ₗ[R] A') (h : Function.Surjective (LinearMap.rTensor B' α)) :
    Function.Surjective α := by
  have hc : (LinearMap.range α).mkQ ∘ₗ α = 0 := by
    ext x
    exact (Submodule.Quotient.mk_eq_zero _).mpr ⟨x, rfl⟩
  have hzero : LinearMap.rTensor B' (LinearMap.range α).mkQ = 0 := by
    refine LinearMap.ext fun z => ?_
    obtain ⟨w, rfl⟩ := h z
    show (LinearMap.rTensor B' (LinearMap.range α).mkQ) ((LinearMap.rTensor B' α) w) = 0
    rw [← LinearMap.comp_apply, ← LinearMap.rTensor_comp, hc, LinearMap.rTensor_zero]
    rfl
  haveI hsub : Subsingleton ((A' ⧸ LinearMap.range α) ⊗[R] B') := by
    constructor
    intro x y
    have hsurj : Function.Surjective (LinearMap.rTensor B' (LinearMap.range α).mkQ) :=
      LinearMap.rTensor_surjective B' (Submodule.mkQ_surjective _)
    obtain ⟨u, hu⟩ := hsurj x
    obtain ⟨v, hv⟩ := hsurj y
    rw [← hu, ← hv, hzero]
    simp
  haveI := subsingleton_of_tensor_invertible R (A' ⧸ LinearMap.range α) B' hsub
  intro a'
  have hq : (Submodule.Quotient.mk a' : A' ⧸ LinearMap.range α) = 0 := Subsingleton.elim _ _
  exact (Submodule.Quotient.mk_eq_zero _).mp hq

/-- ★★★★★★★★**`α ⊗ β` が全射なら `α` が全射**（`B'` は可逆）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★`map α β = rTensor B' α ∘ lTensor A β` なので、合成が全射なら外側が全射である。
★★これが底変換で `μ_A ⊗ μ_{A⁻¹}` が同型であることから `μ_A` を取り出す道具である。 -/
theorem surjective_of_map_surjective (R A A' B B' : Type) [CommRing R]
    [AddCommGroup A] [Module R A] [AddCommGroup A'] [Module R A']
    [AddCommGroup B] [Module R B] [AddCommGroup B'] [Module R B'] [Module.Invertible R B']
    (α : A →ₗ[R] A') (β : B →ₗ[R] B')
    (h : Function.Surjective (TensorProduct.map α β)) : Function.Surjective α := by
  refine surjective_of_rTensor_surjective R A A' B' α ?_
  have hcomp : TensorProduct.map α β
      = (LinearMap.rTensor B' α) ∘ₗ (LinearMap.lTensor A β) := by
    ext a b
    rfl
  rw [hcomp] at h
  exact Function.Surjective.of_comp (g := (LinearMap.lTensor A β)) h

/-! ### ★出典の紐付け(`.src`) -/

def surjective_of_map_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(可逆加群でテンソルしても全射性は落ちない——底変換の道具)",
    sectionId := "genell-def-1-1-ii" }

def surjective_of_map_surjective.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Module.Invertible.linearEquiv(Bᵛ ⊗ B ≃ₗ R)"
      (.inMathlib "Module.Invertible.linearEquiv") 4,
    .citation "[mathlib]" "LinearMap.rTensor_surjective(全射のテンソルは全射)"
      (.inMathlib "LinearMap.rTensor_surjective") 4,
    .implicitStep
      ("★これは底変換 deg_K(L|_K) = deg_F(L) の**有限側**のための道具である。" ++
       "★★台帳 presheaf-pullback-global-sections に書いたとおり、" ++
       "mathlib の前層引き戻しは随伴の左随伴としてしか定義されていないので " ++
       "Γ(f^*L) を直接計算できない。★★★そこで μ_{A ⊗ B} ≅ μ_A ⊗ μ_B(delta_pullSec、§9-743)と " ++
       "A ⊗ A⁻¹ ≅ Ō で回し、μ_A ⊗ μ_{A⁻¹} が同型であることから μ_A の全射性を取り出す") 4 ]

end ABC3.Found.Arakelov
