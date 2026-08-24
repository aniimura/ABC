/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Quasi

/-!
# [FrdI] 同型の四角形を越えて `Definition 1.2` の型が保たれる

★`Thm34Quasi.lean` には `isPreStep_congr_iso` と `isCoAngular_congr_iso` があるが、
`IsLinear` / `IsIsometric` / `IsBaseIsomorphism` の版が無い。
★★これらは `isCoAngular_map_of_equiv`(co-angular の**保存**)を使うのに要る ——
擬逆 `e.inverse` について 4 型の保存を示すには、counit の同型を渡る必要があるからである。

## ★中身

`h : f ≫ β.hom = α.hom ≫ f'` から

* `degFr f = degFr f'`(同型の次数は 1)
* `Base f ≫ Base β.hom = Base α.hom ≫ Base f'`
* `Div f = Φ.map (Base α.hom) (Div f')`(同型の零因子は 0)

が出る。★最後の 1 本が `IsIsometric` の版を与える ——
`Φ.map (Base α.hom)` は `Base α.hom` が同型なので**全単射**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section CongrIso

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★同型の四角形を越えて Frobenius 次数は等しい。 -/
theorem degFr_congr_iso {A B A' B' : C} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    P.degFr f = P.degFr f' := by
  have := congrArg P.degFr h
  rwa [P.degFr_comp, P.degFr_comp, degFr_of_isIso P α.hom, degFr_of_isIso P β.hom,
    mul_one, one_mul] at this

/-- ★★同型の四角形を越えて linear 性は不変。 -/
theorem isLinear_congr_iso {A B A' B' : C} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    IsLinear P f ↔ IsLinear P f' := by
  unfold IsLinear
  rw [degFr_congr_iso P f f' α β h]

/-- ★★同型の四角形を越えて base-isomorphism 性は不変。 -/
theorem isBaseIsomorphism_congr_iso {A B A' B' : C} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    IsBaseIsomorphism P f ↔ IsBaseIsomorphism P f' := by
  have hb : P.Base f ≫ P.Base β.hom = P.Base α.hom ≫ P.Base f' := by
    rw [← P.Base_comp, ← P.Base_comp, h]
  haveI hia : IsIso (P.Base α.hom) := ⟨P.Base α.inv, by
    rw [← P.Base_comp, α.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, α.inv_hom_id, P.Base_id]⟩
  haveI hib : IsIso (P.Base β.hom) := ⟨P.Base β.inv, by
    rw [← P.Base_comp, β.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, β.inv_hom_id, P.Base_id]⟩
  constructor
  · intro hf
    haveI : IsIso (P.Base f) := hf
    haveI : IsIso (P.Base f ≫ P.Base β.hom) := inferInstance
    rw [hb] at this
    exact IsIso.of_isIso_comp_left (P.Base α.hom) (P.Base f')
  · intro hf'
    haveI : IsIso (P.Base f') := hf'
    haveI : IsIso (P.Base α.hom ≫ P.Base f') := inferInstance
    rw [← hb] at this
    exact IsIso.of_isIso_comp_right (P.Base f) (P.Base β.hom)

/-- ★★★同型の四角形を越えて零因子は `Φ.map (Base α)` で移る。 -/
theorem Div_congr_iso {A B A' B' : C} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    P.Div f = Φ.map (P.Base α.hom) (P.Div f') := by
  have hd := congrArg P.Div h
  rw [P.Div_comp, P.Div_comp,
    show P.Div β.hom = 0 from isIsometric_of_isIso P β.hom,
    show P.Div α.hom = 0 from isIsometric_of_isIso P α.hom,
    map_zero, zero_add, smul_zero, add_zero,
    show P.degFr β.hom = 1 from degFr_of_isIso P β.hom] at hd
  simpa using hd

/-- ★★★同型の四角形を越えて isometric 性は不変。

★`Φ.map (Base α.hom)` は `Base α.hom` が同型なので全単射である。 -/
theorem isIsometric_congr_iso {A B A' B' : C} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    IsIsometric P f ↔ IsIsometric P f' := by
  have hd := Div_congr_iso P f f' α β h
  haveI hia : IsIso (P.Base α.hom) := ⟨P.Base α.inv, by
    rw [← P.Base_comp, α.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, α.inv_hom_id, P.Base_id]⟩
  constructor
  · intro hf
    have h0 : Φ.map (P.Base α.hom) (P.Div f') = 0 := by rw [← hd]; exact hf
    have h1 := congrArg (Φ.map (@inv _ _ _ _ (P.Base α.hom) hia)) h0
    rw [map_zero, ← MonoidOn.map_comp,
      show @inv _ _ _ _ (P.Base α.hom) hia ≫ P.Base α.hom = 𝟙 _ from
        @IsIso.inv_hom_id _ _ _ _ (P.Base α.hom) hia, MonoidOn.map_id] at h1
    exact h1
  · intro hf'
    show P.Div f = 0
    rw [hd, show P.Div f' = 0 from hf', map_zero]

end CongrIso

/-! ### ★出典の紐付け -/

/-- ★★★locator —— 同型の四角形を越えた `Definition 1.2` の型の不変性
(`isPreStep_congr_iso` / `isCoAngular_congr_iso` の仲間)。 -/
def isIsometric_congr_iso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 21,
    item := "Definition 1.2, (i) — 同型の四角形を越えて linear / isometric / base-iso は不変",
    sectionId := "frdi-def-1-2-i" }

end ABC3.Found.FrdI
