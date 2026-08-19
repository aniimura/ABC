/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Rigid
import ABC3.Found.FrdI.Prop44Ker

/-!
# [FrdI] Theorem 5.1 の土台 —— `Definition 1.3, (iii), (d)` の切片/余切片を使い切る

`Theorem 5.1` の証明は最初から最後まで「`Definition 1.3, (iii), (d)` の 2 つの圏同値」を
道具として使う。`Cor411Rigid.lean` では**余切片の本質的全射性だけ**を取り出したが、
`Theorem 5.1` では

* 余切片 `_A(𝒞^coa-pre) ≃ Order(Φ(A))` の**充満忠実性**(= `Div` が等しければ余域が同型)
* 切片 `(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp` の**本質的全射性**と**充満忠実性**

の 4 つすべてが要る。ここでそれらを 1 か所にまとめる。

原文 (FrdI p.24):
> co-angular pre-steps. Then the natural functors

原文 (FrdI p.25):
> categories.
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. 広い部分圏の同型を `𝒞` に降ろす -/

/-- 広い部分圏 `𝒞^coa-pre` の同型から `𝒞` の同型を取り出す。

★`WideSubcategory` の射は `hom`(`𝒞` の射)と `property`(性質)の組なので、
合成も恒等射も `hom` 成分では `𝒞` のものと**同じ**である。 -/
def wideIsoDown [MorphismProperty.IsMultiplicative (coaPreProp P)]
    {X Y : WideSubcategory (coaPreProp P)} (e : X ≅ Y) : X.obj ≅ Y.obj where
  hom := e.hom.hom
  inv := e.inv.hom
  hom_inv_id := congrArg InducedWideCategory.Hom.hom e.hom_inv_id
  inv_hom_id := congrArg InducedWideCategory.Hom.hom e.inv_hom_id

@[simp] theorem wideIsoDown_hom [MorphismProperty.IsMultiplicative (coaPreProp P)]
    {X Y : WideSubcategory (coaPreProp P)} (e : X ≅ Y) :
    (wideIsoDown e).hom = e.hom.hom := rfl

@[simp] theorem wideIsoDown_inv [MorphismProperty.IsMultiplicative (coaPreProp P)]
    {X Y : WideSubcategory (coaPreProp P)} (e : X ≅ Y) :
    (wideIsoDown e).inv = e.inv.hom := rfl

/-! ## ★2. 余切片側 —— `A` から出る co-angular pre-step は `Div` で決まる -/

/-- ★★★**余切片の充満忠実性** ——
`A` から出る 2 本の co-angular pre-step の `Div` が等しければ、
余域は `A` の下で同型である。

原文 (FrdI p.25):
> categories.
-/
theorem coaPreStep_iso_of_div_eq (G : Frobenioid P)
    {A B B' : C} (ϵ : A ⟶ B) (ϵ' : A ⟶ B')
    (hp : coaPreProp P ϵ) (hp' : coaPreProp P ϵ')
    (h : P.Div ϵ = P.Div ϵ') :
    ∃ θ : B ≅ B', ϵ ≫ θ.hom = ϵ' := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  let Z : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (⟨ϵ, hp⟩ : (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩)
  let Z' : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Under.mk (⟨ϵ', hp'⟩ : (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B'⟩)
  have hobj : (coaPreUnderFunctor P A).obj Z = (coaPreUnderFunctor P A).obj Z' := by
    show toOrderCat (P.Div ϵ) = toOrderCat (P.Div ϵ')
    rw [h]
  let e : Z ≅ Z' := (coaPreUnderFunctor P A).preimageIso (eqToIso hobj)
  refine ⟨wideIsoDown ((Under.forget _).mapIso e), ?_⟩
  exact congrArg InducedWideCategory.Hom.hom (Under.w e.hom)

/-! ## ★3. 切片側 —— `A` に入る co-angular pre-step は `Φ(Base)⁻¹(Div)` で決まる -/

/-- ★★★**切片の本質的全射性** ——
`Φ(A)` のどの元も、`A` に**入る** co-angular pre-step から
`Φ(Base φ)⁻¹(Div φ)` として得られる。

原文 (FrdI p.25):
> an arrow ψ : B →A the element (ψ∗)−1(Div(ψ)) ∈Φ(A) [since ψ∗: Φ(A) ∼→Φ(B)
-/
theorem exists_coaPreStep_into (G : Frobenioid P) (A : C)
    (x : Φ.val (P.toElem.obj A).base) :
    ∃ (B : C) (φ : B ⟶ A) (hp : coaPreProp P φ),
      Φ.map (@inv _ _ _ _ (P.Base φ) hp.2.2) (P.Div φ) = x := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := coaPreOverFunctor P A) (Opposite.op (toOrderCat x))
  haveI : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
  refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property, ?_⟩
  have e' := e.unop
  have h1 : MLe x (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)) := leOfHom e'.hom
  have h2 : MLe (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)) x := leOfHom e'.inv
  exact mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 h2 h1

/-- ★★★**切片の充満忠実性** ——
`A` に入る 2 本の co-angular pre-step の `Φ(Base)⁻¹(Div)` が等しければ、
域は `A` の上で同型である。 -/
theorem coaPreStep_into_iso_of_div_eq (G : Frobenioid P)
    {A B B' : C} (φ : B ⟶ A) (φ' : B' ⟶ A)
    (hp : coaPreProp P φ) (hp' : coaPreProp P φ')
    (h : Φ.map (@inv _ _ _ _ (P.Base φ) hp.2.2) (P.Div φ)
       = Φ.map (@inv _ _ _ _ (P.Base φ') hp'.2.2) (P.Div φ')) :
    ∃ θ : B ≅ B', θ.hom ≫ φ' = φ := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv A
  let Z : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (⟨φ, hp⟩ : (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩)
  let Z' : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
    Over.mk (⟨φ', hp'⟩ : (⟨B'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩)
  have hobj : (coaPreOverFunctor P A).obj Z = (coaPreOverFunctor P A).obj Z' := by
    haveI : IsIso (P.Base φ) := hp.2.2
    haveI : IsIso (P.Base φ') := hp'.2.2
    show Opposite.op (toOrderCat (Φ.map (inv (P.Base φ)) (P.Div φ)))
      = Opposite.op (toOrderCat (Φ.map (inv (P.Base φ')) (P.Div φ')))
    rw [h]
  let e : Z ≅ Z' := (coaPreOverFunctor P A).preimageIso (eqToIso hobj)
  refine ⟨wideIsoDown ((Over.forget _).mapIso e), ?_⟩
  exact congrArg InducedWideCategory.Hom.hom (Over.w e.hom)

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Definition 1.3, (iii), (d)` の圏同値を使い切る条。 -/
def coaPreStep_iso_of_div_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 25,
    item := "Definition 1.3, (iii), (d) — 圏同値の充満忠実性",
    sectionId := "frdi-def-1-3-iiid" }

end ABC3.Found.FrdI
