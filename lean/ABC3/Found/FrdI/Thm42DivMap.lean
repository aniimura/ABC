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

end DivMap

end ABC3.Found.FrdI
