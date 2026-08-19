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

/-- ★★★★★**右辺は `Ψ` で移る**(像の中の対象について)。

★`IsStep` / `IsPreStep` の保存と反射、`IsIrreducibleMor` の保存(圏論的)だけを使う。 -/
theorem prop114vRhs_map_image (Ψ : C₁ ≌ C₂)
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.functor.map f))
    (hPS' : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P₁ f)
    {A : C₁} (φ : A ⟶ A) (h : Prop114vRhs P₁ φ)
    (B₀ : C₁) (α : Ψ.functor.obj A ⟶ Ψ.functor.obj B₀) (hα : IsStep P₂ α) :
    ∃ (B' : C₂) (ψ : Ψ.functor.obj B₀ ⟶ B') (β : Ψ.functor.obj B₀ ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P₂ ψ ∧ IsStep P₂ β ∧
        α ≫ ψ = Ψ.functor.map φ ≫ α ≫ β := by
  obtain ⟨α₀, rfl⟩ := Ψ.functor.map_surjective α
  have hstep₀ : IsStep P₁ α₀ := by
    refine ⟨hPS' α₀ hα.1, fun hiso => hα.2 ?_⟩
    haveI := hiso
    exact Ψ.functor.map_isIso α₀
  obtain ⟨B'₀, ψ₀, β₀, hirr, hnps, hstepβ, hsq⟩ := h B₀ α₀ hstep₀
  refine ⟨Ψ.functor.obj B'₀, Ψ.functor.map ψ₀, Ψ.functor.map β₀,
    isIrreducibleMor_map Ψ.functor hirr, fun hc => hnps (hPS' ψ₀ hc),
    ⟨hPS β₀ hstepβ.1, fun hiso => hstepβ.2 ?_⟩, ?_⟩
  · haveI := hiso
    exact isIso_of_reflects_iso β₀ Ψ.functor
  · rw [← CategoryTheory.Functor.map_comp, ← CategoryTheory.Functor.map_comp,
      ← CategoryTheory.Functor.map_comp, hsq]

def prop114vRhs_map_image.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (i) — Proposition 1.14, (v) の右辺は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

end ABC3.Found.FrdI
