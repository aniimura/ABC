/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34VBase
import ABC3.Found.FrdI.Prop16

/-!
# [FrdI] Theorem 3.4, (v) —— 任意射が誘導するスライス関手 `(𝒞^pl-bk)_A ⥤ (𝒞^pl-bk)_A'`

原文 (FrdI p.67):
> by sending an object φ : C → A of (Pi)A to the object C → A of (Pi)A which

★★★原文の言うとおり、`A ⟶ A'` の射 `f` は
「`p ≫ f` の `Definition 1.3, (iv), (a)` 分解の **pull-back 部分**を取る」
ことでスライス関手を誘導する。

## ★設計 —— 「選択」は対象だけ、射は**普遍性で一意に決まる**

`arbFactor` は分解の**存在**しか言わないので対象の像には選択が要る。
★しかし**射の像は選択を含まない** ——
`α_W : Y_W ⟶ A'` が pull-back なので、`Definition 1.2, (ii)` の全単射から
「合成と底が一致する射」は**ちょうど 1 本**である。
★★★したがって `map_id` / `map_comp` は**一意性から自動的に出る**。

## ★底との 1-可換図式

原文 (FrdI p.67):
> that this functor ﬁts into a natural 1-commutative diagram

★分解の `γ ≫ β`(Frobenius 型 ＋ pre-step)は**どちらも base-isomorphism** なので、
`Base (γ ≫ β)` が `Base Z ≅ Base Y_Z` を与える。これが四角形の成分である。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

section SlicePush

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★**`Definition 1.3, (iv), (a)` の 3 分解を束ねた構造**。 -/
structure ArbFac {A B : C} (φ : A ⟶ B) where
  /-- Frobenius 型の行き先 -/
  mid : C
  /-- pre-step の行き先 -/
  top : C
  /-- Frobenius 型の因子 -/
  frb : A ⟶ mid
  /-- pre-step の因子 -/
  pre : mid ⟶ top
  /-- pull-back の因子 -/
  plb : top ⟶ B
  /-- 分解式 -/
  fac : φ = frb ≫ pre ≫ plb
  /-- 第 1 因子は Frobenius 型 -/
  hfrb : IsFrobeniusType P frb
  /-- 第 2 因子は pre-step -/
  hpre : IsPreStep P pre
  /-- 第 3 因子は pull-back -/
  hplb : IsPullBack P plb

/-- ★`frb ≫ pre` は base-isomorphism —— Frobenius 型も pre-step も base-iso だから。 -/
theorem ArbFac.baseIso {A B : C} {φ : A ⟶ B} (t : ArbFac P φ) :
    IsIso (P.Base (t.frb ≫ t.pre)) := by
  haveI : IsIso (P.Base t.frb) := t.hfrb.2
  haveI : IsIso (P.Base t.pre) := t.hpre.2
  rw [P.Base_comp]
  infer_instance

/-- ★`Base Z ≅ Base (top)`。 -/
noncomputable def ArbFac.baseEquiv {A B : C} {φ : A ⟶ B} (t : ArbFac P φ) :
    (P.toElem.obj A).base ≅ (P.toElem.obj t.top).base :=
  haveI := t.baseIso P
  asIso (P.Base (t.frb ≫ t.pre))

/-- ★★**分解の底の式** —— `Base φ = baseEquiv.hom ≫ Base plb`。 -/
theorem ArbFac.base_fac {A B : C} {φ : A ⟶ B} (t : ArbFac P φ) :
    P.Base φ = (t.baseEquiv P).hom ≫ P.Base t.plb := by
  show P.Base φ = P.Base (t.frb ≫ t.pre) ≫ P.Base t.plb
  rw [← P.Base_comp, Category.assoc, ← t.fac]

variable (F : FrobenioidCore P)

include F in
theorem nonempty_arbFac {A B : C} (φ : A ⟶ B) : Nonempty (ArbFac P φ) := by
  obtain ⟨X, Y, γ, β, α, hfac, hγ, hβ, hα⟩ := F.arbFactor φ
  exact ⟨⟨X, Y, γ, β, α, hfac, hγ, hβ, hα⟩⟩

/-- ★選択した 3 分解。 -/
noncomputable def someArbFac {A B : C} (φ : A ⟶ B) : ArbFac P φ :=
  (nonempty_arbFac P F φ).some

/-! ## ★pull-back の普遍性 —— 射の像はここで一意に決まる -/

/-- ★★**pull-back への射は「合成 ＋ 底」で決まる**(`Definition 1.2, (ii)` の単射性)。 -/
theorem pullBack_hom_ext {Y' A' : C} {α : Y' ⟶ A'} (hα : IsPullBack P α) {T : C}
    {k₁ k₂ : T ⟶ Y'} (h1 : k₁ ≫ α = k₂ ≫ α) (h2 : P.Base k₁ = P.Base k₂) : k₁ = k₂ :=
  (hα T).1 (Subtype.ext (Prod.ext h1 h2))

/-- ★★**pull-back への射の存在**(`Definition 1.2, (ii)` の全射性)。 -/
theorem pullBack_lift {Y' A' : C} {α : Y' ⟶ A'} (hα : IsPullBack P α) {T : C}
    (g : T ⟶ A') (b : (P.toElem.obj T).base ⟶ (P.toElem.obj Y').base)
    (hb : P.Base g = b ≫ P.Base α) :
    ∃ k : T ⟶ Y', k ≫ α = g ∧ P.Base k = b := by
  obtain ⟨k, hk⟩ := (hα T).2 ⟨(g, b), hb⟩
  have hp := Subtype.ext_iff.mp hk
  exact ⟨k, congrArg Prod.fst hp, congrArg Prod.snd hp⟩

/-! ## ★★スライス関手の構成 -/

variable {A A' : C}

/-- ★対象 `Z` について選んだ分解。 -/
noncomputable abbrev slPushFac (f : A ⟶ A') (Z : Over (⟨A⟩ : PlBk P)) :
    ArbFac P (Z.hom.hom ≫ f) :=
  someArbFac P F (Z.hom.hom ≫ f)

/-- ★★射 `u` が誘導する**底**の射。 -/
noncomputable def slPushBase (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    (P.toElem.obj (slPushFac P F f Z).top).base ⟶ (P.toElem.obj (slPushFac P F f W).top).base :=
  ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
    ≫ ((slPushFac P F f W).baseEquiv P).hom

/-- ★★底の両立条件 —— これが pull-back の普遍性を起動する。 -/
theorem slPushBase_compat (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    P.Base (slPushFac P F f Z).plb
      = slPushBase P F f u ≫ P.Base (slPushFac P F f W).plb := by
  have hu : u.left.hom ≫ W.hom.hom = Z.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Over.w u)
  have hZ := (slPushFac P F f Z).base_fac P
  have hW := (slPushFac P F f W).base_fac P
  have key : ((slPushFac P F f Z).baseEquiv P).hom ≫ P.Base (slPushFac P F f Z).plb
      = P.Base u.left.hom
        ≫ (((slPushFac P F f W).baseEquiv P).hom ≫ P.Base (slPushFac P F f W).plb) := by
    rw [← hZ, ← hW, ← P.Base_comp, ← Category.assoc, hu]
  show P.Base (slPushFac P F f Z).plb
    = (((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
        ≫ ((slPushFac P F f W).baseEquiv P).hom) ≫ P.Base (slPushFac P F f W).plb
  rw [Category.assoc, Category.assoc]
  exact (Iso.eq_inv_comp _).mpr key

/-- ★★射の像 —— **選択を含まない**(pull-back の普遍性で一意)。 -/
noncomputable def slPushHom (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    (slPushFac P F f Z).top ⟶ (slPushFac P F f W).top :=
  (pullBack_lift P (slPushFac P F f W).hplb (slPushFac P F f Z).plb
    (slPushBase P F f u) (slPushBase_compat P F f u)).choose

theorem slPushHom_comp (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    slPushHom P F f u ≫ (slPushFac P F f W).plb = (slPushFac P F f Z).plb :=
  (pullBack_lift P (slPushFac P F f W).hplb (slPushFac P F f Z).plb
    (slPushBase P F f u) (slPushBase_compat P F f u)).choose_spec.1

theorem slPushHom_base (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    P.Base (slPushHom P F f u) = slPushBase P F f u :=
  (pullBack_lift P (slPushFac P F f W).hplb (slPushFac P F f Z).plb
    (slPushBase P F f u) (slPushBase_compat P F f u)).choose_spec.2

/-- ★像も pull-back —— `isPullBack_of_comp_left`(左簡約)から。 -/
theorem slPushHom_isPullBack (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    IsPullBack P (slPushHom P F f u) :=
  isPullBack_of_comp_left P (slPushHom P F f u) (slPushFac P F f W).plb
    (slPushFac P F f W).hplb
    (by rw [slPushHom_comp]; exact (slPushFac P F f Z).hplb)

/-- ★★★★★**任意射が誘導するスライス関手**(原文 p.67)。

★`map_id` / `map_comp` は **pull-back への射の一意性**から出る。 -/
noncomputable def slicePush (f : A ⟶ A') :
    Over (⟨A⟩ : PlBk P) ⥤ Over (⟨A'⟩ : PlBk P) where
  obj Z := Over.mk (⟨(slPushFac P F f Z).plb, (slPushFac P F f Z).hplb⟩ :
    (⟨(slPushFac P F f Z).top⟩ : PlBk P) ⟶ ⟨A'⟩)
  map {Z W} u := Over.homMk
    (⟨slPushHom P F f u, slPushHom_isPullBack P F f u⟩ :
      (⟨(slPushFac P F f Z).top⟩ : PlBk P) ⟶ ⟨(slPushFac P F f W).top⟩)
    (WideSubcategory.hom_ext _ (slPushHom_comp P F f u))
  map_id Z := by
    refine Over.OverMorphism.ext (WideSubcategory.hom_ext _ ?_)
    show slPushHom P F f (𝟙 Z) = 𝟙 (slPushFac P F f Z).top
    refine pullBack_hom_ext P (slPushFac P F f Z).hplb ?_ ?_
    · rw [slPushHom_comp, Category.id_comp]
    · rw [slPushHom_base, P.Base_id]
      show ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base (𝟙 Z.left.obj)
        ≫ ((slPushFac P F f Z).baseEquiv P).hom = 𝟙 _
      rw [P.Base_id, Category.id_comp, Iso.inv_hom_id]
  map_comp {Z W V} u u' := by
    refine Over.OverMorphism.ext (WideSubcategory.hom_ext _ ?_)
    show slPushHom P F f (u ≫ u') = slPushHom P F f u ≫ slPushHom P F f u'
    refine pullBack_hom_ext P (slPushFac P F f V).hplb ?_ ?_
    · rw [slPushHom_comp, Category.assoc, slPushHom_comp, slPushHom_comp]
    · rw [slPushHom_base, P.Base_comp, slPushHom_base, slPushHom_base]
      show ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base (u ≫ u').left.hom
          ≫ ((slPushFac P F f V).baseEquiv P).hom
        = (((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
            ≫ ((slPushFac P F f W).baseEquiv P).hom)
          ≫ (((slPushFac P F f W).baseEquiv P).inv ≫ P.Base u'.left.hom
            ≫ ((slPushFac P F f V).baseEquiv P).hom)
      show ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base (u.left.hom ≫ u'.left.hom)
          ≫ ((slPushFac P F f V).baseEquiv P).hom = _
      rw [P.Base_comp]
      simp [Category.assoc]

def slicePush.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 67,
    item := "Theorem 3.4, (v) — 任意射が誘導する (𝒞^pl-bk)_A ⥤ (𝒞^pl-bk)_A'",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★底との 1-可換図式

原文 (FrdI p.67):
> that this functor ﬁts into a natural 1-commutative diagram
-/

/-- ★成分 —— `Base Z ≅ Base (top)` の逆。 -/
noncomputable def slicePushSqApp (f : A ⟶ A') (Z : Over (⟨A⟩ : PlBk P)) :
    (slicePush P F f ⋙ plBkOverFunctor P A').obj Z
      ≅ (plBkOverFunctor P A ⋙ Over.map (P.Base f)).obj Z :=
  Over.isoMk (((slPushFac P F f Z).baseEquiv P).symm) (by
    show ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base Z.hom.hom ≫ P.Base f
      = P.Base (slPushFac P F f Z).plb
    have h := (slPushFac P F f Z).base_fac P
    rw [P.Base_comp] at h
    rw [h, Iso.inv_hom_id_assoc])

/-- ★★★★★**[FrdI] Theorem 3.4, (v)** —— スライス関手は底の `Over.map` と 1-可換。 -/
noncomputable def slicePushSquare (f : A ⟶ A') :
    slicePush P F f ⋙ plBkOverFunctor P A'
      ≅ plBkOverFunctor P A ⋙ Over.map (P.Base f) :=
  NatIso.ofComponents (slicePushSqApp P F f) (fun {Z W} u => by
    refine Over.OverMorphism.ext ?_
    show P.Base (slPushHom P F f u) ≫ ((slPushFac P F f W).baseEquiv P).inv
      = ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
    rw [slPushHom_base]
    show (((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
        ≫ ((slPushFac P F f W).baseEquiv P).hom)
        ≫ ((slPushFac P F f W).baseEquiv P).inv
      = ((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
    rw [Category.assoc, Category.assoc, Iso.hom_inv_id, Category.comp_id])

def slicePushSquare.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 67,
    item := "Theorem 3.4, (v) — スライス関手と底の Over.map の 1-可換図式",
    sectionId := "frdi-thm-3-4" }

end SlicePush

end ABC3.Found.FrdI
