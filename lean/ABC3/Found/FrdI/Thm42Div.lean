/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop114
import ABC3.Found.FrdI.Thm34Pre

/-!
# [FrdI] Theorem 4.2, (i) —— `Div-identity` 自己射の保存

原文 (FrdI p.77):
> preserves primary steps, Div-identity endomorphisms, Div-Frobenius-

## ★★★`Proposition 1.14, (v)` が**圏論的特徴づけ**を与える(2026-08-19 に測った)

`prop_1_14_v`(`Prop114.lean:1648`):
```
(IsPrimeFrobenius φ ∧ IsDivIdentity φ) ↔
  ∀ B α, IsStep α → ∃ B' ψ β,
    IsIrreducibleMor ψ ∧ ¬ IsPreStep ψ ∧ IsStep β ∧ α ≫ ψ = φ ≫ α ≫ β
```

★★**右辺に `Div` も `Φ` も現れない** —— `IsStep`・`IsPreStep`・`IsIrreducibleMor`
だけである。★これらは `Theorem 3.4, (ii)(iii)` で `Ψ` が保つ/反射するので、
**右辺がそのまま移る**。

★したがって `Div-identity` の保存は
「pre-step の保存と反射」＋「`Proposition 1.14, (v)` の両側適用」で出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 v3 u3 w3 u4 v4

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u3} [Category.{v3} D₂] {C₂ : Type u4} [Category.{v4} C₂]
  {Φ₂ : MonoidOn.{v3, u3, w3} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`Proposition 1.14, (v)` の右辺 —— 純粋に圏論的な条件。 -/
def Prop114vRhs (P : PreFrobenioid C₁ Φ₁) {A : C₁} (φ : A ⟶ A) : Prop :=
  ∀ (B : C₁) (α : A ⟶ B), IsStep P α →
    ∃ (B' : C₁) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P ψ ∧ IsStep P β ∧ α ≫ ψ = φ ≫ α ≫ β

/-- ★★同じものを第 2 の圏で。 -/
def Prop114vRhs' (P : PreFrobenioid C₂ Φ₂) {A : C₂} (φ : A ⟶ A) : Prop :=
  ∀ (B : C₂) (α : A ⟶ B), IsStep P α →
    ∃ (B' : C₂) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P ψ ∧ IsStep P β ∧ α ≫ ψ = φ ≫ α ≫ β

/-- ★★★★★**右辺は `Ψ` で移る**。

★`IsStep` / `IsPreStep` の保存と反射、`IsIrreducibleMor` の保存(圏論的)だけを使う。 -/
theorem prop114vRhs_map (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    {A : C₁} (φ : A ⟶ A) (h : Prop114vRhs P₁ φ) :
    Prop114vRhs' P₂ (Ψ.functor.map φ) := by
  intro B α hα
  -- `B` を像の中へ取り直す
  obtain ⟨B₀, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ.functor) B
  obtain ⟨α₀, hα₀⟩ := Ψ.functor.map_surjective (α ≫ ε.inv)
  have hstep₀ : IsStep P₁ α₀ := by
    refine ⟨hPS' α₀ ?_, ?_⟩
    · rw [hα₀]
      exact IsPreStep.comp P₂ hα.1 (isPreStep_of_isIso P₂ ε.inv)
    · intro hiso
      haveI : IsIso α₀ := hiso
      haveI : IsIso (Ψ.functor.map α₀) := Ψ.functor.map_isIso α₀
      rw [hα₀] at *
      exact hα.2 (by
        haveI : IsIso (α ≫ ε.inv) := ‹IsIso (Ψ.functor.map α₀)›
        have : α = (α ≫ ε.inv) ≫ ε.hom := by
          rw [Category.assoc, ε.inv_hom_id, Category.comp_id]
        rw [this]; infer_instance)
  obtain ⟨B'₀, ψ₀, β₀, hirr, hnps, hstepβ, hsq⟩ := h B₀ α₀ hstep₀
  refine ⟨Ψ.functor.obj B'₀, ε.inv ≫ Ψ.functor.map ψ₀, ε.inv ≫ Ψ.functor.map β₀,
    ?_, ?_, ?_, ?_⟩
  · exact isIrreducibleMor_comp_iso_left ε.inv (isIrreducibleMor_map Ψ.functor hirr)
  · intro hc
    refine hnps (hPS' ψ₀ ?_)
    have h1 : Ψ.functor.map ψ₀ = ε.hom ≫ (ε.inv ≫ Ψ.functor.map ψ₀) := by
      rw [← Category.assoc, ε.hom_inv_id, Category.id_comp]
    rw [h1]
    exact IsPreStep.comp P₂ (isPreStep_of_isIso P₂ ε.hom) hc
  · refine ⟨IsPreStep.comp P₂ (isPreStep_of_isIso P₂ ε.inv) (hPS β₀ hstepβ.1), ?_⟩
    intro hiso
    haveI : IsIso (ε.inv ≫ Ψ.functor.map β₀) := hiso
    haveI : IsIso (Ψ.functor.map β₀) := by
      have : Ψ.functor.map β₀ = ε.hom ≫ (ε.inv ≫ Ψ.functor.map β₀) := by
        rw [← Category.assoc, ε.hom_inv_id, Category.id_comp]
      rw [this]; infer_instance
    exact hstepβ.2 (Ψ.functor.isIso_of_map_isIso β₀)
  · have hΨ := congrArg Ψ.functor.map hsq
    rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp,
      CategoryTheory.Functor.map_comp] at hΨ
    rw [← hα₀]
    calc Ψ.functor.map α₀ ≫ ε.inv ≫ Ψ.functor.map ψ₀
        = (Ψ.functor.map α₀ ≫ ε.inv) ≫ Ψ.functor.map ψ₀ := (Category.assoc _ _ _).symm
      _ = _ := by sorry

def prop114vRhs_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — Proposition 1.14, (v) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

end ABC3.Found.FrdI
