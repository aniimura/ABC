/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.AppA
import ABC3.Found.FrdI.Prop19

/-!
# [FrdI] Theorem 3.4, (v) —— 底の射の 3 分解

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★★`Ψ_Base : 𝒟₁ → 𝒟₂` の構成には、`𝒟` の任意の射を
**pre-step 2 本と pull-back 1 本**に分解する必要がある。

## ★★★在庫 2 本で出る(2026-08-19 に測った)

| 段 | 在庫 |
|---|---|
| pull-back を取る | `plBkEquiv`(`Definition 1.3, (i), (c)`)の**本質的全射性** |
| 残りの同型を span にする | `preStepSpan`(`Definition 1.3, (i), (b)`) |

★`plBkOverFunctor A' : Over (⟨A'⟩ : PlBk P) ⥤ Over (Base A')` が圏同値なので、
`φ_𝒟 : Base A ⟶ Base A'` を `Over (Base A')` の対象と見て本質的全射性を使えば、
pull-back 射 `ψ : B ⟶ A'` と同型 `θ : Base B ≅ Base A` で
`Base ψ = θ ≫ φ_𝒟` なるものが取れる。
★あとは `preStepSpan` を `θ^{-1}` に当てるだけ。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★1. pull-back で `φ_𝒟` を実現する -/

include F in
/-- ★★★**`𝒟` の任意の射は pull-back の底に同型を継いだもの**。

★`plBkEquiv` の**本質的全射性**そのもの。 -/
theorem exists_pullBack_realizing {A A' : C}
    (φ𝒟 : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ (B : C) (ψ : B ⟶ A') (_ : IsPullBack P ψ)
      (θ : (P.toElem.obj B).base ≅ (P.toElem.obj A).base),
      P.Base ψ = θ.hom ≫ φ𝒟 := by
  haveI := F.plBkEquiv A'
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := plBkOverFunctor P A')
    (Over.mk φ𝒟)
  refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property,
    (Over.forget ((P.toElem.obj A').base)).mapIso e, ?_⟩
  exact (Over.w e.hom).symm

/-! ## ★2. 3 分解 -/

include F in
/-- ★★★★★**[FrdI] Theorem 3.4, (v) の 3 分解**。

原文 (FrdI p.68):
> φD = Base(ψ) ◦ Base(γ) ◦ Base(α)−1

★`α : X ⟶ A`、`γ : X ⟶ B` は pre-step、`ψ : B ⟶ A'` は pull-back。 -/
theorem base_three_factor {A A' : C}
    (φ𝒟 : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ (X B : C) (α : X ⟶ A) (γ : X ⟶ B) (ψ : B ⟶ A')
      (hα : IsPreStep P α) (_ : IsPreStep P γ) (_ : IsPullBack P ψ),
      φ𝒟 = @inv _ _ _ _ (P.Base α) hα.2 ≫ P.Base γ ≫ P.Base ψ := by
  obtain ⟨B, ψ, hψ, θ, hθ⟩ := exists_pullBack_realizing P F φ𝒟
  haveI : IsIso θ.inv := inferInstance
  obtain ⟨X, α, γ, hα, hγ, hspan⟩ := F.preStepSpan A B θ.inv (by infer_instance)
  refine ⟨X, B, α, γ, ψ, hα, hγ, hψ, ?_⟩
  have h1 : θ.inv ≫ P.Base ψ = φ𝒟 := by
    rw [hθ, ← Category.assoc, θ.inv_hom_id, Category.id_comp]
  rw [← h1, hspan, Category.assoc]

def base_three_factor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 底の射の 3 分解",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
