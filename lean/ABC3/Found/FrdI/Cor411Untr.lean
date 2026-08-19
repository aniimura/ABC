/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34EndBs

/-!
# [FrdI] `Ψ^un-tr` を **底の 1-可換図式だけ**に還元する

★★★★★★**本ファイルの主張**: `Ψ^un-tr` を作るのに要るのは
**`𝒪^▷` の保存ではなく `𝒪^×` の保存**であり、後者は
**底の 1-可換図式 `P₁.proj ⋙ Ψ_Base ≅ Ψ ⋙ P₂.proj` だけ**から出る。

## ★測って分かったこと

在庫の `toElem_map_congr_of_congr` は
`hOTri : ∀ X δ, δ ∈ 𝒪^▷(X) → Ψ δ ∈ 𝒪^▷(Ψ X)` を要求していた。
★しかしその証明を読むと、`hOTri` が当てられる `δ` は
`IsUnitEquivalent` の定義に現れるもの、すなわち**つねに `𝒪^×` の元**である。

★★そして `𝒪^×` については **`IsLinear` が自動**である ——
`𝒪^× = {IsBaseIdentity ∧ IsLinear ∧ IsUnit}` だが、単元は同型なので
`degFr = 1` は `degFr_of_isIso` からただちに出る。
★★★したがって `𝒪^×` の保存に要るのは **`IsBaseIdentity` の保存だけ**であり、
それは在庫の `isBaseIdentity_map_of_baseSquare` がまさに与える。

## ★★★★★★これで `Corollary 4.11, (i)` の `psi-untr` の依存が変わる

| | 以前の見立て | 実測後 |
|---|---|---|
| `psi-untr` が要るもの | `𝒪^▷` の一般対象での特徴づけ(`otriEndCond_of_oTri` の仮定落とし) | **`v-psibase`(底の 1-可換図式)だけ** |

★`OTriGenCond` を Frobenius-trivial でない対象へ降ろす作業は
**`Ψ^un-tr` には要らなかった**。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

section OTimes

variable {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
  {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)

/-- ★★★**単元については `IsLinear` が自動**。

★`End A` の単元は圏の同型(`CategoryTheory.isUnit_iff_isIso`)であり、
同型の Frobenius 次数は 1(`degFr_of_isIso`)。 -/
theorem isLinear_of_isUnit {A : Cc} {δ : End A} (h : IsUnit δ) :
    IsLinear P (δ : A ⟶ A) := by
  haveI : IsIso (δ : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso δ).mp h
  exact degFr_of_isIso P (δ : A ⟶ A)

/-- ★★★★**`𝒪^×` の特徴づけ** —— `IsLinear` の条は落ちる。 -/
theorem mem_otimes_iff {A : Cc} (δ : End A) :
    δ ∈ OTimes P A ↔ IsUnit δ ∧ IsBaseIdentity P (δ : A ⟶ A) :=
  ⟨fun h => ⟨h.2, h.1.1⟩, fun h => ⟨⟨h.2, isLinear_of_isUnit P h.1⟩, h.1⟩⟩

def mem_otimes_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21, item := "Definition 1.2, (ii) — 𝒪^×(A)",
    sectionId := "frdi-def-1-2-ii" }

end OTimes

section Transport

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★**底の 1-可換図式だけで `Ψ` は `𝒪^×` を保つ**。 -/
theorem otimes_map_of_baseSquare (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj)
    (X : C₁) (δ : End X) (h : δ ∈ OTimes P₁ X) :
    (Ψ.map (δ : X ⟶ X) : End (Ψ.obj X)) ∈ OTimes P₂ (Ψ.obj X) :=
  (mem_otimes_iff P₂ _).mpr
    ⟨h.2.map (Functor.mapEnd X Ψ),
     isBaseIdentity_map_of_baseSquare Ψ ΨBase sq δ h.1.1⟩

/-- ★★★★★**`isUnitEquivalent_map` の仮定を `𝒪^×` だけに弱めたもの**。

★在庫版は `𝒪^▷` 全体の保存を要求していたが、
`IsUnitEquivalent` に現れる `δ` はつねに `𝒪^×` の元である。 -/
theorem isUnitEquivalent_map_of_otimes (Ψ : C₁ ⥤ C₂)
    (hOTimes : ∀ (X : C₁) (δ : End X), δ ∈ OTimes P₁ X →
      (Ψ.map (δ : X ⟶ X) : End (Ψ.obj X)) ∈ OTimes P₂ (Ψ.obj X))
    {A B : C₁} {α₁ α₂ : A ⟶ B} (h : IsUnitEquivalent P₁ α₁ α₂) :
    IsUnitEquivalent P₂ (Ψ.map α₁) (Ψ.map α₂) := by
  obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ := h
  refine ⟨Ψ.obj Cc, Ψ.map γ, Ψ.map β, Ψ.map (δ : Cc ⟶ Cc), hOTimes Cc δ hδ, ?_, ?_⟩
  · rw [h₁, Ψ.map_comp]
  · rw [h₂, Ψ.map_comp, Ψ.map_comp]

/-- ★★★★★**`Ψ` は `𝒞^un-tr` の同値関係を保つ**(仮定は `𝒪^×` の保存だけ)。 -/
theorem toElem_map_congr_of_otimes (Ψ : C₁ ⥤ C₂) (Fc₁ : FrobenioidCore P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (hOTimes : ∀ (X : C₁) (δ : End X), δ ∈ OTimes P₁ X →
      (Ψ.map (δ : X ⟶ X) : End (Ψ.obj X)) ∈ OTimes P₂ (Ψ.obj X))
    {A B : C₁} (α₁ α₂ : A ⟶ B) (h : P₁.toElem.map α₁ = P₁.toElem.map α₂) :
    P₂.toElem.map (Ψ.map α₁) = P₂.toElem.map (Ψ.map α₂) :=
  prop_3_3_ii_toElem P₂ (isUnitEquivalent_map_of_otimes Ψ hOTimes
    ((prop_3_3_ii P₁ Fc₁ hiso₁ α₁ α₂).mpr
      ⟨congrArg ElemFrobCat.Hom.deg h, congrArg ElemFrobCat.Hom.div h,
        congrArg ElemFrobCat.Hom.base h⟩))

/-- ★★★★★★**`Ψ^un-tr` は底の 1-可換図式だけから作れる**。

★★これが `Corollary 4.11, (i)` の `psi-untr` の本体である ——
`𝒪^▷` の一般対象での特徴づけ(`OTriGenCond` の降下)は**要らない**。 -/
noncomputable def psiUnTrOfBaseSquare (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
    (Fc₁ : FrobenioidCore P₁)
    (h₁ : IsOfQuasiIsotropicType C₁ P₁) (h₂ : IsOfQuasiIsotropicType C₂ P₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X)
    (ΨBase : D₁ ⥤ D₂) (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj) :
    UnTr P₁ ⥤ UnTr P₂ :=
  psiUnTr Ψ h₁ h₂
    (fun α₁ α₂ h => toElem_map_congr_of_otimes Ψ Fc₁ hiso₁
      (otimes_map_of_baseSquare Ψ ΨBase sq) α₁ α₂ h)

def psiUnTrOfBaseSquare.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.11, (i) — Ψ^un-tr（底の 1-可換図式からの構成）",
    sectionId := "frdi-cor-4-11" }

end Transport

end ABC3.Found.FrdI
