/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Prop41v
import ABC3.Found.FrdI.MonoidTransport

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

/-! ## ★★★★★★関手性の下ごしらえ -/

/-- ★base-isomorphism が定める `Φ` の同型。 -/
noncomputable def baseEquivOf (Q : PreFrobenioid C Φ) {A B : C} (f : B ⟶ A)
    (hf : IsIso (Q.Base f)) :
    Φ.val (Q.toElem.obj A).base ≃+ Φ.val (Q.toElem.obj B).base :=
  AddEquiv.ofBijective (Φ.map (Q.Base f))
    (Φ.map_bijective_of_iso (@asIso _ _ _ _ (Q.Base f) hf))

@[simp] theorem baseEquivOf_apply (Q : PreFrobenioid C Φ) {A B : C} (f : B ⟶ A)
    (hf : IsIso (Q.Base f)) (x : Φ.val (Q.toElem.obj A).base) :
    baseEquivOf Q f hf x = Φ.map (Q.Base f) x := rfl

/-- ★`primeMap` の `toPrime` での計算。 -/
theorem primeMap_toPrime {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]
    (e : M ≃+ N) {a : M} (ha : IsPrimaryElt a) :
    primeMap e (toPrime M a ha) = toPrime N (e a) (isPrimaryElt_map e ha) := rfl

/-- ★★★**`𝒪^▷` の四角形の `Div`** —— `f` が linear なら `Φ.map (Base f) (Div u) = Div u'`。 -/
theorem map_base_div_otri (Q : PreFrobenioid C Φ) {A B : C} (f : B ⟶ A)
    (hfl : Q.degFr f = 1)
    {u : End A} (hu : u ∈ OTri Q A) {u' : End B} (hu' : u' ∈ OTri Q B)
    (hsq : f ≫ (((u : End A)) : A ⟶ A) = (((u' : End B)) : B ⟶ B) ≫ f) :
    Φ.map (Q.Base f) (Q.Div (((u : End A)) : A ⟶ A)) = Q.Div (((u' : End B)) : B ⟶ B) := by
  haveI := isCancelAdd_of_isIntegralMonoid _ (Q.divisorial (Q.toElem.obj B).base).1.1
  have h := congrArg Q.Div hsq
  rw [Q.Div_comp, Q.Div_comp,
    show Q.degFr (((u : End A)) : A ⟶ A) = 1 from hu.2,
    show Q.Base (((u' : End B)) : B ⟶ B) = Q.Base (𝟙 B) from hu'.1, Q.Base_id, hfl] at h
  have h2 : Φ.map (Q.Base f) (Q.Div (((u : End A)) : A ⟶ A)) + Q.Div f
      = Q.Div (((u' : End B)) : B ⟶ B) + Q.Div f := by
    have h3 : Φ.map (𝟙 ((Q.toElem.obj B).base)) (Q.Div f) = Q.Div f := Φ.map_id _ _
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul, one_smul, h3] at h
    rw [h, add_comm]
  exact add_right_cancel h2

def map_base_div_otri.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — 𝒪^▷ の四角形の Div",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★関手性 —— pre-step の場合 -/

/-- ★★★★★★**`Ψ_Prime` の関手性(co-angular pre-step の場合)**。

★`otriPull`(`Proposition 1.11, (iv)`)が `f ≫ u = u' ≫ f` を与え、
`Div` を当てると `Φ.map (Base f) (Div u) = Div u'`(`f` が linear だから)。
★★同じ計算を `𝒞₂` 側にも当てればよい。 -/
theorem psiPrime_naturality_preStep (ctx : PrimeCtx P P₂ G G₂ Ψ) (F : FrobenioidCore P)
    (hOTri : ∀ (X : C) (δ : End X), δ ∈ OTri P X →
      ((Ψ.functor.map (((δ : End X)) : X ⟶ X)) : End (Ψ.functor.obj X))
        ∈ OTri P₂ (Ψ.functor.obj X))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    {A B : C} (f : B ⟶ A) (hfc : IsCoAngular P f) (hfs : IsPreStep P f)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    psiPrime ctx B (primeMap (baseEquivOf P f hfs.2) p)
      = primeMap (baseEquivOf P₂ (Ψ.functor.map f) (ctx.PS f hfs).2) (psiPrime ctx A p) := by
  have hu : realizeIn ctx A p ∈ OTri P A := realizeIn_mem ctx A p
  set u' := otriPull P F f hfc hfs.1 ⟨realizeIn ctx A p, hu⟩ with hu'def
  have hsq : f ≫ (((realizeIn ctx A p) : A ⟶ A)) = (((u' : End B)) : B ⟶ B) ≫ f :=
    otriPull_spec P F f hfc hfs.1 ⟨realizeIn ctx A p, hu⟩
  have hdiv : Φ.map (P.Base f) (P.Div (((realizeIn ctx A p) : A ⟶ A)))
      = P.Div (((u' : End B)) : B ⟶ B) :=
    map_base_div_otri P f hfs.1 hu u'.2 hsq
  have hsu' : IsPreStep P (((u' : End B)) : B ⟶ B) := isPreStep_of_otri _ u'.2
  have hsu : IsPreStep P (((realizeIn ctx A p) : A ⟶ A)) := realizeIn_preStep ctx A p
  -- ★左辺の素点を `u'` で書き直す
  have hpu' : IsPrimaryElt (preStepVal P (((u' : End B)) : B ⟶ B) hsu') := by
    rw [preStepVal_of_otri _ u'.2 hsu', ← hdiv, ← preStepVal_of_otri _ hu hsu]
    exact isPrimaryElt_map (baseEquivOf P f hfs.2) (realizeIn_primary ctx A p)
  have hLp : primeMap (baseEquivOf P f hfs.2) p
      = toPrime _ (preStepVal P (((u' : End B)) : B ⟶ B) hsu') hpu' := by
    rw [← realizeIn_toPrime ctx A p, primeMap_toPrime]
    refine toPrime_congr ?_
    show Φ.map (P.Base f) (preStepVal P (((realizeIn ctx A p) : A ⟶ A)) hsu)
      = preStepVal P (((u' : End B)) : B ⟶ B) hsu'
    rw [preStepVal_of_otri _ hu hsu, preStepVal_of_otri _ u'.2 hsu']
    exact hdiv
  rw [hLp, psiPrime_spec ctx _ hsu' hpu']
  -- ★右辺
  rw [← realizeIn_toPrime ctx A p, psiPrime_spec ctx _ hsu (realizeIn_primary ctx A p),
    primeMap_toPrime]
  refine (toPrime_congr ?_).symm
  -- ★`𝒞₂` 側の同じ計算
  have hsq₂ : Ψ.functor.map f ≫ (((Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A)))
        : End (Ψ.functor.obj A)) : Ψ.functor.obj A ⟶ Ψ.functor.obj A)
      = (((Ψ.functor.map ((((u' : End B)) : B ⟶ B))) : End (Ψ.functor.obj B))
        : Ψ.functor.obj B ⟶ Ψ.functor.obj B) ≫ Ψ.functor.map f := by
    rw [← Ψ.functor.map_comp, ← Ψ.functor.map_comp, hsq]
  have hdiv₂ : Φ₂.map (P₂.Base (Ψ.functor.map f))
        (P₂.Div (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A))))
      = P₂.Div (Ψ.functor.map ((((u' : End B)) : B ⟶ B))) :=
    map_base_div_otri P₂ (Ψ.functor.map f) (by rw [hdeg]; exact hfs.1)
      (hOTri A _ hu) (hOTri B _ u'.2) hsq₂
  show Φ₂.map (P₂.Base (Ψ.functor.map f))
      (preStepVal P₂ (Ψ.functor.map (((realizeIn ctx A p) : A ⟶ A))) (ctx.PS _ hsu))
    = preStepVal P₂ (Ψ.functor.map ((((u' : End B)) : B ⟶ B))) (ctx.PS _ hsu')
  rw [preStepVal_of_otri _ (hOTri A _ hu) (ctx.PS _ hsu),
    preStepVal_of_otri _ (hOTri B _ u'.2) (ctx.PS _ hsu')]
  exact hdiv₂

def psiPrime_naturality_preStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — Ψ_Prime の関手性（pre-step）",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★関手性 —— Frobenius 型の場合 -/

/-- ★★★**Frobenius 型の四角形の `Div`** —— `n = degFr v` として
`Φ.map (Base v) (Div u) = n • Div w`。 -/
theorem div_square_frob (Q : PreFrobenioid C Φ) {A B : C} (v : B ⟶ A)
    {w : End B} (hw : w ∈ OTri Q B) {u : End A} (hu : u ∈ OTri Q A)
    (hsq : (((w : End B)) : B ⟶ B) ≫ v = v ≫ (((u : End A)) : A ⟶ A)) :
    Φ.map (Q.Base v) (Q.Div (((u : End A)) : A ⟶ A))
      = ((Q.degFr v : ℕ+) : ℕ) • Q.Div (((w : End B)) : B ⟶ B) := by
  haveI := isCancelAdd_of_isIntegralMonoid _ (Q.divisorial (Q.toElem.obj B).base).1.1
  have h := congrArg Q.Div hsq
  rw [Q.Div_comp, Q.Div_comp,
    show Q.Base (((w : End B)) : B ⟶ B) = Q.Base (𝟙 B) from hw.1, Q.Base_id,
    show Q.degFr (((u : End A)) : A ⟶ A) = 1 from hu.2] at h
  have h3 : Φ.map (𝟙 ((Q.toElem.obj B).base)) (Q.Div v) = Q.Div v := Φ.map_id _ _
  rw [h3, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul] at h
  exact (add_left_cancel (a := Q.Div v)
    (by rw [h, add_comm (Φ.map (Q.Base v) (Q.Div (((u : End A)) : A ⟶ A)))])).symm

def div_square_frob.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — Frobenius 型の四角形の Div",
    sectionId := "frdi-thm-4-2" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`Ψ_Prime` の関手性(Frobenius 型の場合)**。

★`Φ` が perfect なので `Φ.map (Base v) (Div u)` を `n = degFr v` で割れる。
その商を `Div` に持つ `w ∈ 𝒪^▷(B)` を取り、
`Proposition 1.10, (i)` で四角形 `w ≫ v = v ≫ u₀` を立てる。
★★`Div u₀ = Div u` が出るので `faithfulUpToUnits` で `u = u₀ ≫ ε`(単元)。
★★★`Prime` は `MPrec` の商なので `n •` は素点を変えない。 -/
theorem psiPrime_naturality_frobType (ctx : PrimeCtx P P₂ G G₂ Ψ) (F : FrobenioidCore P)
    (hOTri : ∀ (X : C) (δ : End X), δ ∈ OTri P X →
      ((Ψ.functor.map (((δ : End X)) : X ⟶ X)) : End (Ψ.functor.obj X))
        ∈ OTri P₂ (Ψ.functor.obj X))
    (hdeg : ∀ {X Y : C} (g : X ⟶ Y), P₂.degFr (Ψ.functor.map g) = P.degFr g)
    (hFT : ∀ {X Y : C} (g : X ⟶ Y), IsFrobeniusType P g →
      IsFrobeniusType P₂ (Ψ.functor.map g))
    {A B : C} (v : B ⟶ A) (hv : IsFrobeniusType P v)
    (p : Prime (Φ.val (P.toElem.obj A).base)) :
    psiPrime ctx B (primeMap (baseEquivOf P v hv.2) p)
      = primeMap (baseEquivOf P₂ (Ψ.functor.map v) (hFT v hv).2) (psiPrime ctx A p) := by
  classical
  have hu : realizeIn ctx A p ∈ OTri P A := realizeIn_mem ctx A p
  have hsu : IsPreStep P (((realizeIn ctx A p) : A ⟶ A)) := realizeIn_preStep ctx A p
  haveI hvi : IsIso (P.Base v) := hv.2
  -- ★perfect 除法
  obtain ⟨d, hd⟩ := (ctx.perfM B (P.degFr v)).2
    (Φ.map (P.Base v) (P.Div (((realizeIn ctx A p) : A ⟶ A))))
  obtain ⟨w, hwd⟩ := ctx.divS B d
  have hwm : (w : End B) ∈ OTri P B := w.2
  have hsw : IsPreStep P (((w : End B) : B ⟶ B)) := isPreStep_of_otri _ hwm
  -- ★四角形
  obtain ⟨u₀, hsq, -⟩ :=
    prop_1_10_i_exists_given P F (((w : End B) : B ⟶ B)) v hv v hv rfl
  -- ★`u₀ ∈ 𝒪^▷(A)`
  have hb : P.Base (((w : End B) : B ⟶ B)) ≫ P.Base v = P.Base v ≫ P.Base u₀ := by
    rw [← P.Base_comp, ← P.Base_comp, hsq]
  have hu₀b : IsBaseIdentity P u₀ := by
    show P.Base u₀ = P.Base (𝟙 A)
    rw [P.Base_id]
    rw [show P.Base (((w : End B) : B ⟶ B)) = P.Base (𝟙 B) from hwm.1, P.Base_id,
      Category.id_comp] at hb
    exact ((cancel_epi (P.Base v)).mp (by rw [Category.comp_id]; exact hb)).symm
  have hu₀l : IsLinear P u₀ := by
    have hdd : P.degFr (((w : End B) : B ⟶ B) ≫ v) = P.degFr (v ≫ u₀) := by rw [hsq]
    rw [P.degFr_comp, P.degFr_comp, show P.degFr (((w : End B) : B ⟶ B)) = 1 from hwm.2,
      mul_one] at hdd
    show P.degFr u₀ = 1
    exact (mul_right_cancel (b := P.degFr v) (by rw [one_mul]; exact hdd)).symm
  have hu₀m : u₀ ∈ OTri P A := ⟨hu₀b, hu₀l⟩
  have hsu₀ : IsPreStep P u₀ := isPreStep_of_otri _ hu₀m
  -- ★`Div u₀ = Div u`
  have hdiv0 : Φ.map (P.Base v) (P.Div u₀)
      = ((P.degFr v : ℕ+) : ℕ) • P.Div (((w : End B) : B ⟶ B)) :=
    div_square_frob P v hwm hu₀m hsq
  have hdivu : P.Div u₀ = P.Div (((realizeIn ctx A p) : A ⟶ A)) := by
    refine (Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base v) hvi)).1 ?_
    show Φ.map (P.Base v) (P.Div u₀)
      = Φ.map (P.Base v) (P.Div (((realizeIn ctx A p) : A ⟶ A)))
    rw [hdiv0, hwd]
    exact hd
  -- ★単元でずれるだけ
  obtain ⟨ε, hεm, hεeq⟩ := F.faithfulUpToUnits (((realizeIn ctx A p) : A ⟶ A)) u₀
    (show P.Base (((realizeIn ctx A p) : A ⟶ A)) = P.Base u₀ by
      rw [show P.Base (((realizeIn ctx A p) : A ⟶ A)) = P.Base (𝟙 A) from hu.1,
        show P.Base u₀ = P.Base (𝟙 A) from hu₀b])
    hdivu.symm (prop_1_4_i P _ (fun Y _ => ctx.iso Y)) hsu
    (prop_1_4_i P _ (fun Y _ => ctx.iso Y)) hsu₀
  sorry

def psiPrime_naturality_frobType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Theorem 4.2, (ii) — Ψ_Prime の関手性（Frobenius 型）",
    sectionId := "frdi-thm-4-2" }

end PsiPrime

end ABC3.Found.FrdI
