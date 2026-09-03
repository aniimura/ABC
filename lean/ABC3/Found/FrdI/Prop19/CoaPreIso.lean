/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop18
import ABC3.Found.FrdI.Prop19.Prep

/-!
# Prop19 —— co-angular pre-isomorphism

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)


section CoaPreIso

variable [MorphismProperty.IsMultiplicative (coaPreProp P)]

/-- ★★**同じ `Div` を持つ2つの co-angular pre-step は同型で移り合う**。

★使う成分は **充満性と忠実性**。★1.8 では (iii)(d) の忠実性を使わなかったが、
**ここでは使う** —— 「1.8 で使わなかっただけ」であって、
`Definition 1.3, (iii), (d)` の忠実性は**不要な条件ではない**。 -/
theorem coaPre_iso_of_div_eq (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence)
    {A : C} (Z W : Under (⟨A⟩ : WideSubcategory (coaPreProp P)))
    (h : P.Div Z.hom.hom = P.Div W.hom.hom) : Nonempty (Z ≅ W) := by
  haveI := (hequiv A).full
  haveI := (hequiv A).faithful
  have hobj : (coaPreUnderFunctor P A).obj Z = (coaPreUnderFunctor P A).obj W :=
    congrArg toOrderCat h
  obtain ⟨f, -⟩ := (coaPreUnderFunctor P A).map_surjective (eqToHom hobj)
  obtain ⟨g, -⟩ := (coaPreUnderFunctor P A).map_surjective (eqToHom hobj.symm)
  haveI := Preorder.subsingleton_hom ((coaPreUnderFunctor P A).obj Z)
    ((coaPreUnderFunctor P A).obj Z)
  haveI := Preorder.subsingleton_hom ((coaPreUnderFunctor P A).obj W)
    ((coaPreUnderFunctor P A).obj W)
  exact ⟨⟨f, g, (coaPreUnderFunctor P A).map_injective (Subsingleton.elim _ _),
    (coaPreUnderFunctor P A).map_injective (Subsingleton.elim _ _)⟩⟩

/-- ★★**同じ `Div` を持つ co-angular pre-step は、`𝒞` の同型 1 本だけ違う**。

★上の `coaPre_iso_of_div_eq` はコスライスの同型を与える。
★**ここではその `𝒞` 成分を取り出す** —— 使うときはこの形が要る。

★★**`Proposition 1.6, (v)` の `⟸` を証明する道の第 2 歩**である ——
`A` の co-angular pre-step 自己射の**底**が、`Base(Aut_𝒞(A))` の
同じ剰余類に入ることが、これで言える。 -/
theorem coaPre_base_diff (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence)
    {A B E : C} (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : coaPreProp P φ) (hψ : coaPreProp P ψ) (h : P.Div φ = P.Div ψ) :
    ∃ f : B ⟶ E, IsIso f ∧ φ ≫ f = ψ := by
  obtain ⟨e⟩ := coaPre_iso_of_div_eq P hequiv
    (Under.mk (Y := (⟨B⟩ : WideSubcategory (coaPreProp P)))
      (⟨φ, hφ⟩ : (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩))
    (Under.mk (Y := (⟨E⟩ : WideSubcategory (coaPreProp P)))
      (⟨ψ, hψ⟩ : (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨E⟩)) h
  refine ⟨(CommaMorphism.right e.hom).hom, ⟨(CommaMorphism.right e.inv).hom, ?_, ?_⟩, ?_⟩
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (CommaMorphism.right t)) e.hom_inv_id
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (CommaMorphism.right t)) e.inv_hom_id
  · exact congrArg InducedWideCategory.Hom.hom (Under.w e.hom)

end CoaPreIso

/-! ## ★`𝒞^imtr-pre` —— isometric pre-step が定める部分圏

原文 (FrdI p.31):
> Write Cimtr-pre ⊆C for the subcategory determined by the isometric pre-steps
-/

instance : (isometricPreStepProp P).ContainsIdentities :=
  ⟨fun A => ⟨isIsometric_id P A, isPreStep_id P A⟩⟩

instance : (isometricPreStepProp P).IsStableUnderComposition :=
  ⟨fun _ _ hf hg => ⟨IsIsometric.comp P hf.1 hg.1, IsPreStep.comp P hf.2 hg.2⟩⟩

instance : (isometricPreStepProp P).IsMultiplicative where

/-- **`𝒞^imtr-pre`** —— isometric pre-step が定める広い部分圏。

★`𝒞^coa-pre`(`Definition 1.3, (iii), (d)`)と同じ形で作れる ——
`ContainsIdentities` と `IsStableUnderComposition` が
`isIsometric_id` / `isPreStep_id` / `IsIsometric.comp` / `IsPreStep.comp` から
**そのまま出る**(co-angular と違って `Definition 1.3` の条項を引かない)。 -/
abbrev ImtrPre : Type u2 := WideSubcategory (isometricPreStepProp P)

/-- ★**`φ_*` の対象への割り当て** —— `Proposition 1.9, (i)` の分解の
「isometric pre-step 側」を取る。

★★**ここで選択が要る。** `prop_1_9_ii_obj` は **`∃` であって `∃!` ではない** ——
終域は `prop_1_9_i_uniq` により**同型を除いてしか決まらない**。
★対照的に**射の割り当てには選択が要らない**(`imtrPre_hom_uniq`)。
`#print axioms` で `Classical.choice` が入るのはこの定義である。 -/
noncomputable def pushObj (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) : C :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose

/-- ★上で選んだ対象への isometric pre-step。 -/
noncomputable def pushHom (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    pushObj P F φ hφ ε hεi hεs ⟶ B :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose

/-- ★`pushHom` は本当に isometric pre-step。 -/
theorem pushHom_spec (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsIsometric P (pushHom P F φ hφ ε hεi hεs) ∧
      IsPreStep P (pushHom P F φ hφ ε hεi hεs) :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.2.2

/-! ## ★第3段 —— (v) の機械: `𝒞^istr` と isotropification

原文 (FrdI p.32):
> (v) Cistr [equipped with the restriction to C of the given functor C →FΦ] is a

## ★★測定: 部分圏の形が第2段と違う

* `𝒞^imtr-pre` は**射**の条件 → `WideSubcategory`(`InducedWideCategory`)
* `𝒞^istr` は**対象**の条件 → `ObjectProperty.FullSubcategory`(`InducedCategory`)

★**第2段で作った部分圏の機械はここでは使えない。** 見込みどおりだった。
★ただし `InducedCategory.homEquiv` があるので、射の扱いは `WideSubcategory` より**軽い**
(性質のフィールドが無い)。

## ★mathlib の実測(S1–S4、2026-08-15)

| 要るもの | mathlib | 判定 |
|---|---|---|
| 対象の条件による充満部分圏 | `ObjectProperty.FullSubcategory` / `.ι` | ★**使う** |
| induced category の射 | `InducedCategory.homEquiv` / `homMk` | ★**使う** |
| hom 同値から左随伴を作る | `Adjunction.leftAdjointOfEquiv` / `adjunctionOfEquivLeft` | ★★**使う** |

★★**`leftAdjointOfEquiv` は関手と随伴を一度に作る。** isotropic hull の普遍性(`∃!`)が
そのまま hom 同値になるので、**関手を手で組み立てる必要がない**。
-/

end ABC3.Found.FrdI
