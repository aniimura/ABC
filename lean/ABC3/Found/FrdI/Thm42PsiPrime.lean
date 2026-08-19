/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Prop41v

/-!
# [FrdI] Theorem 4.2, (ii) —— `Ψ_Prime` の構成

原文 (FrdI p.80):
> class of primary steps to or from Ai [where the correspondence between elements of

★★`Prime(Φ(A))` は **primary 元の商**(`Quotient (precSetoid _)`)であり、
`Div : 𝒪^▷ ↠ Φ` により **どの素点も `A` の自己射で実現される**。
★★★そして「同じ素点を定めるか」は `SamePrimeCat`(圏論的)で判定できる
(`toPrime_preStepVal_eq_iff`)。

★したがって `Ψ_Prime` は
「素点 `p` を実現する自己射 `u` を取り、`Ψ u` の定める素点を返す」で定義でき、
`SamePrimeCat` の移送(`samePrimeCat_map` / `_rev`)が
**well-defined 性と単射性**を、`primaryInIff_all` が**全射性**を与える。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w u2 v2

section PsiPrime

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

variable (P P₂) in
/-- ★★`Ψ_Prime` を作るのに要る仮定を束ねたもの。

★`primaryIff` は `primaryInIff_all`(`Thm42Prop41v.lean`)が供給する。 -/
structure PrimeCtx (G : Frobenioid P) (G₂ : Frobenioid P₂) (Ψ : C ≌ C₂) : Prop where
  /-- `𝒞` は isotropic 型 -/
  iso : ∀ X : C, IsIsotropic P X
  /-- `𝒞₂` は isotropic 型 -/
  iso₂ : ∀ X : C₂, IsIsotropic P₂ X
  /-- `Div : 𝒪^▷ ↠ Φ`(逸脱 (B)) -/
  divS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
    ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a
  /-- 同上(第 2 の圏) -/
  divS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
    ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a
  /-- `Φ` は perf-factorial -/
  pf : ∀ Y : C, ∃ ι, IsPerfFactorialWith (Φ.val (P.toElem.obj Y).base) ι
  /-- 同上(第 2 の圏) -/
  pf₂ : ∀ Y : C₂, ∃ ι, IsPerfFactorialWith (Φ₂.val (P₂.toElem.obj Y).base) ι
  /-- `Φ` の値は perfect -/
  perfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base)
  /-- 同上(第 2 の圏) -/
  perfM₂ : ∀ Y : C₂, IsPerfectMonoid (Φ₂.val (P₂.toElem.obj Y).base)
  /-- `Ψ` は co-angular を保つ -/
  CA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.functor.map f)
  /-- `Ψ` は pre-step を保つ -/
  PS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.functor.map f)
  /-- `Ψ` は pre-step を反射する -/
  PS' : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P₂ (Ψ.functor.map f) → IsPreStep P f
  /-- `Ψ` は step を保つ -/
  ST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.functor.map f)
  /-- `Ψ` は step を反射する -/
  ST' : ∀ {X Y : C} (f : X ⟶ Y), IsStep P₂ (Ψ.functor.map f) → IsStep P f
  /-- ★`Theorem 4.2, (i)` の結論(`primaryInIff_all` が供給する) -/
  primaryIff : ∀ (A E : C) (ϵ : E ⟶ A), IsPreStep P ϵ →
    (IsPrimaryElt (P.Div ϵ) ↔ IsPrimaryElt (P₂.Div (Ψ.functor.map ϵ)))

variable {G : Frobenioid P} {G₂ : Frobenioid P₂} {Ψ : C ≌ C₂}

/-- ★素点 `p` を実現する `A` の自己射(選択)。 -/
noncomputable def realizeIn (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) : End A :=
  (exists_otri_realizing ctx.divS p).choose

theorem realizeIn_mem (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) : realizeIn ctx A p ∈ OTri P A :=
  (exists_otri_realizing ctx.divS p).choose_spec.choose

theorem realizeIn_preStep (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    IsPreStep P (((realizeIn ctx A p) : A ⟶ A)) :=
  isPreStep_of_otri _ (realizeIn_mem ctx A p)

theorem realizeIn_primary (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    IsPrimaryElt (preStepVal P (((realizeIn ctx A p) : A ⟶ A)) (realizeIn_preStep ctx A p)) :=
  (exists_otri_realizing ctx.divS p).choose_spec.choose_spec.choose_spec.choose

theorem realizeIn_toPrime (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    toPrime _ (preStepVal P (((realizeIn ctx A p) : A ⟶ A)) (realizeIn_preStep ctx A p))
      (realizeIn_primary ctx A p) = p :=
  (exists_otri_realizing ctx.divS p).choose_spec.choose_spec.choose_spec.choose_spec

def realizeIn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — 素点を実現する自己射",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★`Ψ_Prime` 本体 -/

/-- ★`Ψ u` も primary。 -/
theorem map_realizeIn_primary (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    IsPrimaryElt (preStepVal P₂ (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
      (ctx.PS _ (realizeIn_preStep ctx A p))) := by
  refine (isPrimaryElt_preStepVal_iff _ _).mpr ?_
  refine (ctx.primaryIff A A _ (realizeIn_preStep ctx A p)).mp ?_
  exact (isPrimaryElt_preStepVal_iff _ (realizeIn_preStep ctx A p)).mp
    (realizeIn_primary ctx A p)

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii)** —— `Ψ_Prime`。 -/
noncomputable def psiPrime (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    Prime (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base) :=
  toPrime _ (preStepVal P₂ (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
    (ctx.PS _ (realizeIn_preStep ctx A p))) (map_realizeIn_primary ctx A p)

/-- ★★★★★★**`Ψ_Prime` の特徴づけ** ——
どの primary pre-step で実現しても同じ答えになる。 -/
theorem psiPrime_spec (ctx : PrimeCtx P P₂ G G₂ Ψ) {A E : C} (ϵ : E ⟶ A)
    (hsϵ : IsPreStep P ϵ) (hpϵ : IsPrimaryElt (preStepVal P ϵ hsϵ)) :
    psiPrime ctx A (toPrime _ (preStepVal P ϵ hsϵ) hpϵ)
      = toPrime _ (preStepVal P₂ (Ψ.functor.map ϵ) (ctx.PS ϵ hsϵ))
          ((isPrimaryElt_preStepVal_iff _ _).mpr
            ((ctx.primaryIff A E ϵ hsϵ).mp
              ((isPrimaryElt_preStepVal_iff ϵ hsϵ).mp hpϵ))) := by
  obtain ⟨ι₁, H₁⟩ := ctx.pf A
  obtain ⟨ι₂, H₂⟩ := ctx.pf₂ (Ψ.functor.obj A)
  set p := toPrime _ (preStepVal P ϵ hsϵ) hpϵ with hp
  have hsame : SamePrimeCat P (((realizeIn ctx A p) : A ⟶ A)) ϵ := by
    refine (toPrime_preStepVal_eq_iff G ctx.iso H₁ (ctx.perfM A) (P.divisorial _)
      (((realizeIn ctx A p) : A ⟶ A)) (prop_1_4_i P _ (fun Y _ => ctx.iso Y))
      (realizeIn_preStep ctx A p) (realizeIn_primary ctx A p)
      ϵ (prop_1_4_i P _ (fun Y _ => ctx.iso Y)) hsϵ hpϵ).mp ?_
    rw [realizeIn_toPrime ctx A p, hp]
  have hsame₂ : SamePrimeCat P₂ (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
      (Ψ.functor.map ϵ) := samePrimeCat_map Ψ.functor ctx.CA ctx.PS ctx.ST hsame
  exact (toPrime_preStepVal_eq_iff G₂ ctx.iso₂ H₂ (ctx.perfM₂ _) (P₂.divisorial _)
    (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
    (prop_1_4_i P₂ _ (fun Y _ => ctx.iso₂ Y)) _ (map_realizeIn_primary ctx A p)
    (Ψ.functor.map ϵ) (prop_1_4_i P₂ _ (fun Y _ => ctx.iso₂ Y)) _ _).mpr hsame₂

def psiPrime.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — Ψ_Prime",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★全単射 -/

/-- ★`toPrime` は値の等式で移る。 -/
theorem toPrime_congr {M : Type w} [AddCommMonoid M] {a b : M} (h : a = b)
    {ha : IsPrimaryElt a} {hb : IsPrimaryElt b} : toPrime M a ha = toPrime M b hb := by
  subst h; rfl

/-- ★`preStepVal` は射の等式で移る。 -/
theorem preStepVal_congr {A B : C} {f g : B ⟶ A} (h : f = g)
    (hf : IsPreStep P f) (hg : IsPreStep P g) :
    preStepVal P f hf = preStepVal P g hg := by
  subst h; rfl

/-- ★★★★★**`Ψ_Prime` は単射**。 -/
theorem psiPrime_injective (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C) :
    Function.Injective (psiPrime ctx A) := by
  intro p q h
  obtain ⟨ι₁, H₁⟩ := ctx.pf A
  obtain ⟨ι₂, H₂⟩ := ctx.pf₂ (Ψ.functor.obj A)
  rw [← realizeIn_toPrime ctx A p, ← realizeIn_toPrime ctx A q]
  refine (toPrime_preStepVal_eq_iff G ctx.iso H₁ (ctx.perfM A) (P.divisorial _)
    (((realizeIn ctx A p) : A ⟶ A)) (prop_1_4_i P _ (fun Y _ => ctx.iso Y))
    (realizeIn_preStep ctx A p) (realizeIn_primary ctx A p)
    (((realizeIn ctx A q) : A ⟶ A)) (prop_1_4_i P _ (fun Y _ => ctx.iso Y))
    (realizeIn_preStep ctx A q) (realizeIn_primary ctx A q)).mpr ?_
  refine samePrimeCat_map_rev Ψ ctx.iso ctx.PS' ctx.ST' ?_
  refine (toPrime_preStepVal_eq_iff G₂ ctx.iso₂ H₂ (ctx.perfM₂ _) (P₂.divisorial _)
    (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
    (prop_1_4_i P₂ _ (fun Y _ => ctx.iso₂ Y)) _ (map_realizeIn_primary ctx A p)
    (Ψ.functor.map (((realizeIn ctx A q) : A ⟶ A)))
    (prop_1_4_i P₂ _ (fun Y _ => ctx.iso₂ Y)) _ (map_realizeIn_primary ctx A q)).mp h

/-- ★★★★★**`Ψ_Prime` は全射**。 -/
theorem psiPrime_surjective (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C) :
    Function.Surjective (psiPrime ctx A) := by
  intro q
  obtain ⟨v, hv, hvs, hvp, hvq⟩ := exists_otri_realizing ctx.divS₂ q
  obtain ⟨w, hw⟩ := Ψ.functor.map_surjective (((v : End (Ψ.functor.obj A))
    : Ψ.functor.obj A ⟶ Ψ.functor.obj A))
  have hsw : IsPreStep P w := by
    refine ctx.PS' w ?_
    rw [hw]; exact hvs
  have hpw : IsPrimaryElt (preStepVal P w hsw) := by
    refine (isPrimaryElt_preStepVal_iff w hsw).mpr ?_
    refine (ctx.primaryIff A A w hsw).mpr ?_
    rw [hw]
    exact (isPrimaryElt_preStepVal_iff _ hvs).mp hvp
  refine ⟨toPrime _ (preStepVal P w hsw) hpw, ?_⟩
  rw [psiPrime_spec ctx w hsw hpw, ← hvq]
  exact toPrime_congr (preStepVal_congr hw (ctx.PS w hsw) hvs)

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii)** —— `Ψ_Prime` は**全単射**。 -/
theorem psiPrime_bijective (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C) :
    Function.Bijective (psiPrime ctx A) :=
  ⟨psiPrime_injective ctx A, psiPrime_surjective ctx A⟩

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii)** —— `Prime` の間の同型。 -/
noncomputable def psiPrimeEquiv (ctx : PrimeCtx P P₂ G G₂ Ψ) (A : C) :
    Prime (Φ.val (P.toElem.obj A).base)
      ≃ Prime (Φ₂.val (P₂.toElem.obj (Ψ.functor.obj A)).base) :=
  Equiv.ofBijective _ (psiPrime_bijective ctx A)

def psiPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (ii) — Ψ_Prime は全単射",
    sectionId := "frdi-thm-4-2" }

end PsiPrime

end ABC3.Found.FrdI
