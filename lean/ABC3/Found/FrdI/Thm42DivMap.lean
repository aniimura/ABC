/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42DivId

/-!
# [FrdI] Theorem 4.2, (ii)(iii) —— `Ψ` が誘導する `Φ(A) ≅ Φ₂(ΨA)`

原文 (FrdI p.78):
> then the last two equivalences of categories of (ii) arise from isomorphisms of

## ★★★★★★測って分かった —— **逸脱 (B) があると直接に単系同型が出る**

原文は `{A(𝒞^coa-pre)}_p`(`A` から出る primary step のなす圏)の同値から
`Order(Φ(A)_p)` の同値を作り、さらに Div-Frobenius-trivial 性を使って
加法性(単系同型)へ上げる。★**primary step どうしは合成できない**(行き先が違う)ので、
加法性を出すのに図式が要る。

★★我々は **`Div : 𝒪^▷(A) ↠ Φ(A)`**(逸脱 (B)、`prop_4_4` で開示済)を持つ。
`𝒪^▷(A)` は**自己射**なので**合成できる**。したがって:
- `Div v = Div v'` ⟹ `v ≅ v'`(コスライスの圏同値、`exists_iso_of_div_eq`)⟹
  `Div (Ψ v) = Div (Ψ v')` —— **well-defined**
- `Div (a * b) = Div a + Div b`(`𝒪^▷` は base-identity かつ linear)⟹ **加法的**
- `Ψ.symm` で同じものを作れば **全単射**

★★★これで `Ψ_Φ` の対象ごとの成分が出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section DivMap

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★`𝒪^▷` の積の `Div` は和。 -/
theorem div_mul_otri (Q : PreFrobenioid C Φ) {A : C} {a b : End A}
    (ha : a ∈ OTri Q A) (hb : b ∈ OTri Q A) :
    Q.Div (((a * b : End A)) : A ⟶ A)
      = Q.Div (((a : End A)) : A ⟶ A) + Q.Div (((b : End A)) : A ⟶ A) := by
  show Q.Div ((((b : End A)) : A ⟶ A) ≫ (((a : End A)) : A ⟶ A)) = _
  rw [Q.Div_comp, show Q.Base (((b : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hb.1, Q.Base_id,
    show Q.degFr (((a : End A)) : A ⟶ A) = 1 from ha.2]
  have h3 : Φ.map (𝟙 ((Q.toElem.obj A).base)) (Q.Div (((a : End A)) : A ⟶ A))
    = Q.Div (((a : End A)) : A ⟶ A) := Φ.map_id _ _
  rw [h3, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]

/-- ★★★**well-defined 性** —— `Div` が等しい `𝒪^▷` の元は `Ψ` の像でも `Div` が等しい。 -/
theorem div_map_eq_of_div_eq (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A : C} {v v' : End A}
    (hv : v ∈ OTri P A) (hv' : v' ∈ OTri P A)
    (h : P.Div (((v : End A)) : A ⟶ A) = P.Div (((v' : End A)) : A ⟶ A)) :
    P₂.Div (Ψ.functor.map (((v : End A)) : A ⟶ A))
      = P₂.Div (Ψ.functor.map (((v' : End A)) : A ⟶ A)) := by
  obtain ⟨θ, hθiso, hθ⟩ := exists_iso_of_div_eq G (((v : End A)) : A ⟶ A)
    (((v' : End A)) : A ⟶ A) (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ hv)
    (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ hv') h
  haveI := hθiso
  haveI : IsIso (Ψ.functor.map θ) := Ψ.functor.map_isIso θ
  rw [← hθ, Ψ.functor.map_comp]
  exact (div_comp_iso (P₂ := P₂) (Ψ.functor.map (((v : End A)) : A ⟶ A))
    (Ψ.functor.map θ)).symm

def div_map_eq_of_div_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (iii) — Ψ が誘導する Div の対応の well-defined 性",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★`Ψ_Φ` の対象ごとの成分 -/

variable (P) in
/-- ★`Div` の値を実現する `𝒪^▷` の自己射(選択)。 -/
noncomputable def realizeDiv
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : C) (x : Φ.val (P.toElem.obj A).base) : End A :=
  (((hdivS A x).choose : End A))

theorem realizeDiv_mem
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : C) (x : Φ.val (P.toElem.obj A).base) : realizeDiv P hdivS A x ∈ OTri P A :=
  (hdivS A x).choose.2

theorem realizeDiv_div
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : C) (x : Φ.val (P.toElem.obj A).base) :
    P.Div (((realizeDiv P hdivS A x) : A ⟶ A)) = x :=
  (hdivS A x).choose_spec

/-- ★★★★★★**`Ψ` が誘導する `Φ(A) → Φ₂(ΨA)`** —— 加法準同型。 -/
noncomputable def divMap (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (A : C) :
    Φ.val (P.toElem.obj A).base →+ Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base where
  toFun x := P₂.Div (Ψ.functor.map (((realizeDiv P hdivS A x) : A ⟶ A)))
  map_zero' := by
    have h0 : P.Div (((realizeDiv P hdivS A 0) : A ⟶ A))
        = P.Div ((((1 : End A)) : A ⟶ A)) := by
      rw [realizeDiv_div]
      exact (P.Div_id A).symm
    rw [div_map_eq_of_div_eq Ψ G hiso (realizeDiv_mem hdivS A 0) (OTri P A).one_mem h0]
    show P₂.Div (Ψ.functor.map (𝟙 A)) = 0
    rw [CategoryTheory.Functor.map_id]
    exact P₂.Div_id _
  map_add' x y := by
    have hmx := realizeDiv_mem hdivS A x
    have hmy := realizeDiv_mem hdivS A y
    have hxy : P.Div (((realizeDiv P hdivS A (x + y)) : A ⟶ A))
        = P.Div ((((realizeDiv P hdivS A x * realizeDiv P hdivS A y : End A)) : A ⟶ A)) := by
      rw [realizeDiv_div, div_mul_otri P hmx hmy, realizeDiv_div, realizeDiv_div]
    rw [div_map_eq_of_div_eq Ψ G hiso (realizeDiv_mem hdivS A (x + y))
      ((OTri P A).mul_mem hmx hmy) hxy]
    have hmul := map_mul (CategoryTheory.Functor.mapEnd A Ψ.functor)
      (realizeDiv P hdivS A x) (realizeDiv P hdivS A y)
    show P₂.Div ((((CategoryTheory.Functor.mapEnd A Ψ.functor)
        (realizeDiv P hdivS A x * realizeDiv P hdivS A y))
        : Ψ.functor.obj A ⟶ Ψ.functor.obj A)) = _
    rw [hmul]
    exact div_mul_otri P₂ (hOTri A _ hmx) (hOTri A _ hmy)

def divMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (iii) — Ψ が誘導する Φ(A) → Φ₂(ΨA)",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★全単射 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`divMap` は単射**。 -/
theorem divMap_injective (Ψ : C ≌ C₂) (G : Frobenioid P) (G₂ : Frobenioid P₂)
    (hiso : ∀ X : C, IsIsotropic P X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (A : C) : Function.Injective (divMap Ψ G hiso hdivS hOTri A) := by
  intro x y hxy
  have hmx := realizeDiv_mem hdivS A x
  have hmy := realizeDiv_mem hdivS A y
  obtain ⟨θ₂, hθ₂iso, hθ₂⟩ := exists_iso_of_div_eq G₂
    (Ψ.functor.map (((realizeDiv P hdivS A x) : A ⟶ A)))
    (Ψ.functor.map (((realizeDiv P hdivS A y) : A ⟶ A)))
    (prop_1_4_i P₂ _ (fun Z _ => hiso₂ Z)) (isPreStep_of_otri _ (hOTri A _ hmx))
    (prop_1_4_i P₂ _ (fun Z _ => hiso₂ Z)) (isPreStep_of_otri _ (hOTri A _ hmy)) hxy
  obtain ⟨θ, hθ⟩ := Ψ.functor.map_surjective θ₂
  haveI := hθ₂iso
  haveI : IsIso θ := by
    haveI : IsIso (Ψ.functor.map θ) := by rw [hθ]; infer_instance
    exact isIso_of_reflects_iso θ Ψ.functor
  have heq : (((realizeDiv P hdivS A x) : A ⟶ A)) ≫ θ
      = (((realizeDiv P hdivS A y) : A ⟶ A)) := by
    refine Ψ.functor.map_injective ?_
    rw [Ψ.functor.map_comp, hθ]
    exact hθ₂
  have h2 : P.Div (((realizeDiv P hdivS A x) : A ⟶ A))
      = P.Div (((realizeDiv P hdivS A y) : A ⟶ A)) := by
    rw [← heq]
    exact (div_comp_iso (P₂ := P) (((realizeDiv P hdivS A x) : A ⟶ A)) θ).symm
  rw [realizeDiv_div, realizeDiv_div] at h2
  exact h2

/-- ★★★★★**`divMap` は全射**。 -/
theorem divMap_surjective (Ψ : C ≌ C₂) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hOTri' : ∀ (Z : C) (δ : End Z),
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z)
    (A : C) : Function.Surjective (divMap Ψ G hiso hdivS hOTri A) := by
  intro y
  obtain ⟨v, hvy⟩ := hdivS₂ (Ψ.functor.obj A) y
  obtain ⟨w, hw⟩ : ∃ t : End A, Ψ.functor.map (((t : End A)) : A ⟶ A)
      = (((v : End (Ψ.functor.obj A))) : Ψ.functor.obj A ⟶ Ψ.functor.obj A) :=
    Ψ.functor.map_surjective _
  have hwm : w ∈ OTri P A := hOTri' A w (by rw [hw]; exact v.2)
  refine ⟨P.Div (((w : End A)) : A ⟶ A), ?_⟩
  show P₂.Div (Ψ.functor.map (((realizeDiv P hdivS A
    (P.Div (((w : End A)) : A ⟶ A))) : A ⟶ A))) = y
  rw [div_map_eq_of_div_eq Ψ G hiso
    (realizeDiv_mem hdivS A (P.Div (((w : End A)) : A ⟶ A))) hwm
    (by rw [realizeDiv_div]), hw]
  exact hvy

/-- ★★★★★★**[FrdI] Theorem 4.2, (iii)** —— `Ψ` は単系同型 `Φ(A) ≅ Φ₂(ΨA)` を誘導する。 -/
noncomputable def divEquiv (Ψ : C ≌ C₂) (G : Frobenioid P) (G₂ : Frobenioid P₂)
    (hiso : ∀ X : C, IsIsotropic P X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C) (δ : End Z), δ ∈ OTri P Z →
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z))
    (hOTri' : ∀ (Z : C) (δ : End Z),
      ((Ψ.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψ.functor.obj Z))
        ∈ OTri P₂ (Ψ.functor.obj Z) → δ ∈ OTri P Z)
    (A : C) :
    Φ.val (P.toElem.obj A).base ≃+ Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base :=
  AddEquiv.ofBijective (divMap Ψ G hiso hdivS hOTri A)
    ⟨divMap_injective Ψ G G₂ hiso hiso₂ hdivS hOTri A,
     divMap_surjective Ψ G hiso hdivS hdivS₂ hOTri hOTri' A⟩

def divEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 78,
    item := "Theorem 4.2, (iii) — Ψ が誘導する単系同型 Φ(A) ≅ Φ₂(ΨA)",
    sectionId := "frdi-thm-4-2" }

end DivMap

end ABC3.Found.FrdI
