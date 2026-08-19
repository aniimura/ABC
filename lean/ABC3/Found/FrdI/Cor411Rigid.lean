/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] Corollary 4.11, (i) の rigidity

原文 (FrdI p.93):
> is a Div-identity automorphism. [Indeed, this follows by applying the functoriality

## ★★★★★測って分かった段取り(2026-08-19)

原文は 3 段で書く:
1. `α ∈ Aut(𝒞ᵢ → 𝒞ᵢ^un-tr)` が誘導する自己同型は **Div-identity**
   (co-angular pre-step への関手性 ＋ `Definition 1.3, (iii), (d)` の第 2 圏同値)
2. `𝒟ᵢ` が **Div-slim** なので **base-identity**
3. `𝒞ᵢ^un-tr` は **unit-trivial** 型なので **恒等**

★★1 の中身を測ると、実は**圏同値そのものは要らず**、
「どの `x ∈ Φ(A)` も `A` から出る co-angular pre-step の零因子である」
(＝圏同値の**本質的全射性**)だけでよい。★あとは `Div_comp` の計算である:

- `Div (α ≫ ϵ) = Φ.map (Base α) (Div ϵ) + (degFr ϵ) • Div α = Φ.map (Base α) x`
- `Div (ϵ ≫ β) = Φ.map (Base ϵ) (Div β) + (degFr β) • Div ϵ = x`

(`α`・`β` は同型なので零因子 0、`ϵ` は pre-step なので次数 1。)
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section DivIdOfAut

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**`Definition 1.3, (iii), (d)` の本質的全射性** ——
`Φ(A)` のどの元も `A` から出る co-angular pre-step の零因子である。

★圏同値の像は `Order(Φ(A))` の中で**同型**を与えるだけだが、
`Φ(A)` は integral かつ sharp なので `mle_antisymm` で**等式**に上がる。 -/
theorem exists_coaPreStep_div (G : Frobenioid P) (A : C)
    (x : Φ.val (P.toElem.obj A).base) :
    ∃ (B : C) (ϵ : A ⟶ B), IsCoAngular P ϵ ∧ IsPreStep P ϵ ∧ P.Div ϵ = x := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := coaPreUnderFunctor P A) (toOrderCat x)
  refine ⟨Z.right.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
  have h1 : MLe (P.Div Z.hom.hom) x := leOfHom e.hom
  have h2 : MLe x (P.Div Z.hom.hom) := leOfHom e.inv
  exact mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 h1 h2

/-- ★★★★**自己同型が co-angular pre-step すべてに沿って自然なら Div-identity**。 -/
theorem isDivIdentity_of_preStep_natural (G : Frobenioid P) {A : C} (α : A ⟶ A)
    [IsIso α]
    (hnat : ∀ (B : C) (ϵ : A ⟶ B), IsCoAngular P ϵ → IsPreStep P ϵ →
      ∃ β : B ⟶ B, IsIso β ∧ α ≫ ϵ = ϵ ≫ β) :
    IsDivIdentity P α := by
  show Φ.map (P.Base α) = Φ.map (P.Base (𝟙 A))
  refine AddMonoidHom.ext (fun x => ?_)
  obtain ⟨B, ϵ, hϵc, hϵs, hdiv⟩ := exists_coaPreStep_div G A x
  obtain ⟨β, hβ, hsq⟩ := hnat B ϵ hϵc hϵs
  haveI := hβ
  have hda : P.Div α = 0 := isIsometric_of_isIso P α
  have hdb : P.Div β = 0 := isIsometric_of_isIso P β
  have hL : P.Div (α ≫ ϵ) = Φ.map (P.Base α) x := by
    rw [P.Div_comp, hda, hdiv, smul_zero, add_zero]
  have hR : P.Div (ϵ ≫ β) = x := by
    rw [P.Div_comp, hdb, map_zero, zero_add, degFr_of_isIso P β, hdiv]
    show ((1 : ℕ+) : ℕ) • x = x
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hkey : Φ.map (P.Base α) x = x := by
    rw [← hL, hsq, hR]
  refine hkey.trans ?_
  show x = Φ.map (P.Base (𝟙 A)) x
  rw [P.Base_id, Φ.map_id]

def isDivIdentity_of_preStep_natural.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — 誘導される自己同型は Div-identity",
    sectionId := "frdi-cor-4-11" }

end DivIdOfAut

end ABC3.Found.FrdI
