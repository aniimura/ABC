/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop25iii

/-!
# [FrdI] Proposition 4.4, (ii) の最後の 1 条 —— 十分条件を測る

★`Prop44Core.lean` で `𝒞^birat` の `FrobenioidCore` 21 条のうち **20 条**が済み、
残るのは `otriBase`(`Definition 1.3, (iii), (c)` の「全単射は `Base(φ)` にしか依らない」)
1 条だけである(`Gap/FrdI/Prop44.lean`)。

★★**測定(2026-08-18)**: `𝒞^birat` では co-angular pre-step が**同型**なので、

  `otriBase` ⟺ `𝒪^▷(A^birat)` が**可換**

であった。そして `Prop25iii.lean` の `otri_mul_comm` は

  **`A` が Frobenius-normalized ⟹ `𝒪^▷(A)` は可換**

を与える。したがって

  ★★★**`𝒞` が birationally Frobenius-normalized 型なら `otriBase` は出る** ——
  すなわち **`𝒞^birat` は `FrobenioidCore` の 21 条をすべて満たす**。

★★**穴が残るのは「birationally Frobenius-normalized でない Frobenioid」だけ**である。
原文 `Example 4.6` はそういう Frobenioid が実在することを示しているので、
**穴は空虚ではない**。★逆に、`Theorem 5.2` の model Frobenioid や
`Example 6.1`・`6.3` の幾何的・数論的な例はすべて model 型(したがって
birationally Frobenius-normalized 型)なので、**そこでは本条は埋まっている**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★1. `𝒪^▷(A^birat)` が可換なら `otriBase` -/

include P G in
/-- ★★**`𝒪^▷` の可換性から `otriBase`** ——
`𝒞^birat` の co-angular pre-step は同型なので、`φ' = φ ≫ w`(`w ∈ 𝒪^▷(B)`)と書けて、
結論は `w` と `β` の可換性に落ちる。 -/
theorem birat_otriBase_of_comm
    (hcomm : ∀ (X : BiratCat P G) (x y : OTri (biratPre P G) X), x * y = y * x)
    {A B : BiratCat P G} (φ φ' : A ⟶ B)
    (hcφ : IsCoAngular (biratPre P G) φ) (hsφ : IsPreStep (biratPre P G) φ)
    (_hcφ' : IsCoAngular (biratPre P G) φ') (hsφ' : IsPreStep (biratPre P G) φ')
    (hbase : (biratPre P G).Base φ = (biratPre P G).Base φ')
    (α : End A) (_hα : α ∈ OTri (biratPre P G) A)
    (β : End B) (hβ : β ∈ OTri (biratPre P G) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' := by
  haveI hiso : IsIso φ := birat_isIso_of_coaPre_birat φ hcφ hsφ
  -- ★`w := φ⁻¹ ≫ φ'` は `𝒪^▷(B)` の元
  have hwb : (biratPre P G).Base (inv φ ≫ φ') = (biratPre P G).Base (𝟙 B) := by
    rw [(biratPre P G).Base_comp, (biratPre P G).Base_id, ← hbase]
    haveI : IsIso ((biratPre P G).Base φ) := hsφ.2
    have hinv : (biratPre P G).Base (inv φ) = inv ((biratPre P G).Base φ) := by
      refine IsIso.eq_inv_of_hom_inv_id ?_
      rw [← (biratPre P G).Base_comp, IsIso.hom_inv_id, (biratPre P G).Base_id]
    rw [hinv, IsIso.inv_hom_id]
  have hwl : (biratPre P G).degFr (inv φ ≫ φ') = 1 := by
    rw [(biratPre P G).degFr_comp, hsφ'.1,
      birat_degFr_of_inv (P := P) (G := G) (f := φ) (f' := inv φ)
        (IsIso.hom_inv_id φ) hsφ.1, mul_one]
  have hw : (inv φ ≫ φ' : End B) ∈ OTri (biratPre P G) B := ⟨hwb, hwl⟩
  -- ★可換性を `w` と `β` に当てる
  have hc := hcomm B ⟨(inv φ ≫ φ' : End B), hw⟩ ⟨β, hβ⟩
  have hc' : (β : B ⟶ B) ≫ (inv φ ≫ φ') = (inv φ ≫ φ') ≫ (β : B ⟶ B) :=
    congrArg (fun t : OTri (biratPre P G) B => ((t : End B) : B ⟶ B)) hc
  have hφ' : φ' = φ ≫ (inv φ ≫ φ') := by
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  calc (φ' ≫ β : A ⟶ B) = (φ ≫ (inv φ ≫ φ')) ≫ β := by rw [← hφ']
    _ = φ ≫ ((inv φ ≫ φ') ≫ β) := Category.assoc _ _ _
    _ = φ ≫ ((β : B ⟶ B) ≫ (inv φ ≫ φ')) := by rw [hc']
    _ = (φ ≫ β) ≫ (inv φ ≫ φ') := (Category.assoc _ _ _).symm
    _ = ((α : A ⟶ A) ≫ φ) ≫ (inv φ ≫ φ') := by rw [h]
    _ = (α : A ⟶ A) ≫ (φ ≫ (inv φ ≫ φ')) := Category.assoc _ _ _
    _ = (α : A ⟶ A) ≫ φ' := by rw [← hφ']

/-! ## ★2. birationally Frobenius-normalized 型なら 21 条がすべて揃う -/

include P G in
/-- ★★★**`𝒞` が birationally Frobenius-normalized 型なら `otriBase` が出る**。 -/
theorem birat_otriBase_of_frobNormalized
    (hbfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {A B : BiratCat P G} (φ φ' : A ⟶ B)
    (hcφ : IsCoAngular (biratPre P G) φ) (hsφ : IsPreStep (biratPre P G) φ)
    (hcφ' : IsCoAngular (biratPre P G) φ') (hsφ' : IsPreStep (biratPre P G) φ')
    (hbase : (biratPre P G).Base φ = (biratPre P G).Base φ')
    (α : End A) (hα : α ∈ OTri (biratPre P G) A)
    (β : End B) (hβ : β ∈ OTri (biratPre P G) B)
    (h : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' :=
  birat_otriBase_of_comm P G
    (fun X x y => otri_mul_comm (biratPre P G) (hbfn X) x y)
    φ φ' hcφ hsφ hcφ' hsφ' hbase α hα β hβ h

include P G in
/-- ★★★★**`𝒞` が birationally Frobenius-normalized 型なら
`𝒞^birat` は `Definition 1.3` の 21 条をすべて満たす**。

★★これが `Proposition 4.4, (ii)` の穴を「そうでない `𝒞`」だけに絞る。 -/
theorem birat_frobenioidCore_of_frobNormalized
    (hbfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    FrobenioidCore (biratPre P G) where
  baseSurj := birat_baseSurj P G
  preStepSpan := birat_preStepSpan P G
  plBkEquiv := birat_plBkEquiv P G
  frobDegSurj := birat_frobDegSurj P G
  frobDegUniq := birat_frobDegUniq P G
  coAngularComp := birat_coAngularComp P G
  coAngularOfPreStep := birat_coAngularOfPreStep P G
  otriFwd := birat_otriFwd P G
  otriBwd := birat_otriBwd P G
  otriBase := fun φ φ' hcφ hsφ hcφ' hsφ' hbase α hα β hβ h =>
    birat_otriBase_of_frobNormalized P G hbfn φ φ' hcφ hsφ hcφ' hsφ' hbase α hα β hβ h
  arbFactor := birat_arbFactor P G
  arbFactorUniq := birat_arbFactorUniq P G
  pullBackLB := birat_pullBackLB P G
  preStepMono := birat_preStepMono P G
  preStepFactor := birat_preStepFactor P G
  preStepFactorUniq := birat_preStepFactorUniq P G
  preStepFactor' := birat_preStepFactor' P G
  preStepFactorUniq' := birat_preStepFactorUniq' P G
  faithfulUpToUnits := birat_faithfulUpToUnits P G
  isotropicHullExists := birat_isotropicHullExists P G
  isotropicClosed := birat_isotropicClosed P G

/-- ★**locator** —— `Proposition 4.4, (ii)` の 21 条(★**条つき**:
`𝒞` が birationally Frobenius-normalized 型のときに限る)。 -/
def birat_frobenioidCore_of_frobNormalized.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85,
    item := "Proposition 4.4, (ii) — Definition 1.3 の 21 条(birat-Frobenius-normalized 型のとき)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
