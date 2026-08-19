/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] Corollary 4.11, (ii) —— base-identity 自己射の輸送(一般形)

原文 (FrdI p.93):
> such that every endomorphism [of an object of Cibirat] induced by

## ★★★★★測って分かった一般形(2026-08-19)

`Cor411Otimes.lean` の `otimes_map_of_divSlim` は `𝒪^×` について書いたが、
★**同じ骨が「base-identity 自己射」一般について通る**。
しかも **Div-slim の試験に使う単系 `Φ₀` は、pre-Frobenioid の単系と別でよい** ——
原文 (ii) が `𝒞^birat`(単系は `0_𝒟`)の base-identity 自己射を
**もとの `Φᵢ`** で試験するのはこの形である。

★手筋(原文どおり):
1. `Base α = 𝟙` を pull-back に沿って一意に持ち上げる(`Proposition 1.11, (iii)`)。
2. その族はスライスの**自己同型**を定める(自然性は pull-back の単射性)。
3. `Φ₀` が恒等へ送る ＋ **Div-slim** ⟹ 恒等 ⟹ 終対象で `Base (Ψ α) = 𝟙`。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BaseIdPull

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-- ★★**base-identity 自己射の pull-back に沿った持ち上げ**。 -/
noncomputable def baseIdPull (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) : B ⟶ B :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose : B ⟶ B)

theorem baseIdPull_base (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    Q.Base (baseIdPull Q z hz α hα) = 𝟙 _ :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose_spec.1).1

theorem baseIdPull_isBaseIdentity (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    IsBaseIdentity Q (baseIdPull Q z hz α hα) := by
  show Q.Base (baseIdPull Q z hz α hα) = Q.Base (𝟙 B)
  rw [Q.Base_id]
  exact baseIdPull_base Q z hz α hα

theorem baseIdPull_spec (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    z ≫ α = baseIdPull Q z hz α hα ≫ z :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose_spec.1).2

/-- ★★★**自然性** —— pull-back の単射性から。 -/
theorem baseIdPull_natural (Q : PreFrobenioid C Φ) {A : C} (α : A ⟶ A)
    (hα : IsBaseIdentity Q α) {Z W : C} (f : Z ⟶ W) (w : W ⟶ A) (hw : IsPullBack Q w)
    (hz : IsPullBack Q (f ≫ w)) :
    baseIdPull Q (f ≫ w) hz α hα ≫ f = f ≫ baseIdPull Q w hw α hα := by
  have hsZ := baseIdPull_spec Q (f ≫ w) hz α hα
  have hsW := baseIdPull_spec Q w hw α hα
  have hcomp : (baseIdPull Q (f ≫ w) hz α hα ≫ f) ≫ w
      = (f ≫ baseIdPull Q w hw α hα) ≫ w := by
    calc (baseIdPull Q (f ≫ w) hz α hα ≫ f) ≫ w
        = baseIdPull Q (f ≫ w) hz α hα ≫ (f ≫ w) := by rw [Category.assoc]
      _ = (f ≫ w) ≫ α := hsZ.symm
      _ = f ≫ (w ≫ α) := by rw [Category.assoc]
      _ = f ≫ (baseIdPull Q w hw α hα ≫ w) := by rw [hsW]
      _ = (f ≫ baseIdPull Q w hw α hα) ≫ w := by rw [Category.assoc]
  have hbase : Q.Base (baseIdPull Q (f ≫ w) hz α hα ≫ f)
      = Q.Base (f ≫ baseIdPull Q w hw α hα) := by
    rw [Q.Base_comp, Q.Base_comp, baseIdPull_base, baseIdPull_base,
      Category.id_comp, Category.comp_id]
  exact (hw Z).1 (Subtype.ext (Prod.ext hcomp hbase))

theorem baseIdPull_congr (Q : PreFrobenioid C Φ) {A B : C} {z z' : B ⟶ A} (h : z = z')
    (hz : IsPullBack Q z) (hz' : IsPullBack Q z') (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    baseIdPull Q z hz α hα = baseIdPull Q z' hz' α hα := by
  subst h; rfl

def baseIdPull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — base-identity 自己射の pull-back に沿った持ち上げ",
    sectionId := "frdi-cor-4-11" }

end BaseIdPull

end ABC3.Found.FrdI
