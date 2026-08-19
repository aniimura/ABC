/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34VBase
import ABC3.Found.FrdI.Prop16
import ABC3.Found.FrdI.Thm34IvGen

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

/-! ## ★同値による簡約 -/

section Cancel

universe v9 u9 v8 u8 v7 u7

/-- ★★**同値で前合成した同型は、元の同型を与える**。 -/
noncomputable def cancelLeftEquiv {J : Type u9} [Category.{v9} J] {K : Type u8}
    [Category.{v8} K] {L : Type u7} [Category.{v7} L]
    (e : J ≌ K) {X Y : K ⥤ L} (iso : e.functor ⋙ X ≅ e.functor ⋙ Y) : X ≅ Y :=
  (e.invFunIdAssoc X).symm ≪≫ Functor.isoWhiskerLeft e.inverse iso ≪≫ e.invFunIdAssoc Y

end Cancel

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

/-! ## ★★★2 つの 3 分解の canonical な比較

★★`arbFactorUniq`(原文の一意性)は**使わない** ——
pull-back の普遍性だけで比較同型が**canonical に**作れる。
★これで「選択に依存しない」ことが `Definition 1.2, (ii)` だけから出る。 -/

section Compare

variable {A B : C}

/-- ★比較射の底の成分。 -/
noncomputable def cmpBase {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') :
    (P.toElem.obj t.top).base ⟶ (P.toElem.obj t'.top).base :=
  (t.baseEquiv P).inv ≫ (t'.baseEquiv P).hom

theorem cmpBase_compat {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    P.Base t.plb = cmpBase P t t' ≫ P.Base t'.plb := by
  have hZ := t.base_fac P
  have hW := t'.base_fac P
  have key : (t.baseEquiv P).hom ≫ P.Base t.plb = (t'.baseEquiv P).hom ≫ P.Base t'.plb := by
    rw [← hZ, ← hW, h]
  show P.Base t.plb = ((t.baseEquiv P).inv ≫ (t'.baseEquiv P).hom) ≫ P.Base t'.plb
  rw [Category.assoc]
  exact (Iso.eq_inv_comp _).mpr key

/-- ★比較射(片道)。 -/
noncomputable def cmpHom {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    t.top ⟶ t'.top :=
  (pullBack_lift P t'.hplb t.plb (cmpBase P t t') (cmpBase_compat P t t' h)).choose

theorem cmpHom_comp {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    cmpHom P t t' h ≫ t'.plb = t.plb :=
  (pullBack_lift P t'.hplb t.plb (cmpBase P t t') (cmpBase_compat P t t' h)).choose_spec.1

theorem cmpHom_base {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    P.Base (cmpHom P t t' h) = cmpBase P t t' :=
  (pullBack_lift P t'.hplb t.plb (cmpBase P t t') (cmpBase_compat P t t' h)).choose_spec.2

/-- ★★★★**比較同型** —— 逆向きの比較射と合わせると恒等になる(一意性)。 -/
noncomputable def cmpIso {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    t.top ≅ t'.top where
  hom := cmpHom P t t' h
  inv := cmpHom P t' t h.symm
  hom_inv_id := by
    refine pullBack_hom_ext P t.hplb ?_ ?_
    · rw [Category.assoc, cmpHom_comp, cmpHom_comp, Category.id_comp]
    · rw [P.Base_comp, cmpHom_base, cmpHom_base, P.Base_id]
      show ((t.baseEquiv P).inv ≫ (t'.baseEquiv P).hom)
        ≫ ((t'.baseEquiv P).inv ≫ (t.baseEquiv P).hom) = 𝟙 _
      rw [Category.assoc, ← Category.assoc (t'.baseEquiv P).hom, Iso.hom_inv_id,
        Category.id_comp, Iso.inv_hom_id]
  inv_hom_id := by
    refine pullBack_hom_ext P t'.hplb ?_ ?_
    · rw [Category.assoc, cmpHom_comp, cmpHom_comp, Category.id_comp]
    · rw [P.Base_comp, cmpHom_base, cmpHom_base, P.Base_id]
      show ((t'.baseEquiv P).inv ≫ (t.baseEquiv P).hom)
        ≫ ((t.baseEquiv P).inv ≫ (t'.baseEquiv P).hom) = 𝟙 _
      rw [Category.assoc, ← Category.assoc (t.baseEquiv P).hom, Iso.hom_inv_id,
        Category.id_comp, Iso.inv_hom_id]

/-- ★比較同型は pull-back 部分と可換。 -/
theorem cmpIso_comp {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    (cmpIso P t t' h).hom ≫ t'.plb = t.plb :=
  cmpHom_comp P t t' h

/-- ★比較同型は pull-back なので `𝒞^pl-bk` の同型を与える。 -/
theorem cmpHom_isPullBack {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    IsPullBack P (cmpHom P t t' h) :=
  isPullBack_of_comp_left P (cmpHom P t t' h) t'.plb t'.hplb
    (by rw [cmpHom_comp]; exact t.hplb)

/-- ★★★★★**比較同型の自然性(抽象形)** ——
`k₁`・`k₂` が同じ底の射 `w` を持つなら、比較同型と可換。

★★具体的な構成を一切含まないので、`Ψ` の側でそのまま使える。 -/
theorem cmpHom_natural {X Y B' : C} {φ₁ φ₁' : X ⟶ B'} {φ₂ φ₂' : Y ⟶ B'}
    (t₁ : ArbFac P φ₁) (t₁' : ArbFac P φ₁') (t₂ : ArbFac P φ₂) (t₂' : ArbFac P φ₂')
    (h₁ : φ₁ = φ₁') (h₂ : φ₂ = φ₂')
    (w : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base)
    (k₁ : t₁.top ⟶ t₂.top) (hk₁c : k₁ ≫ t₂.plb = t₁.plb)
    (hk₁b : P.Base k₁ = (t₁.baseEquiv P).inv ≫ w ≫ (t₂.baseEquiv P).hom)
    (k₂ : t₁'.top ⟶ t₂'.top) (hk₂c : k₂ ≫ t₂'.plb = t₁'.plb)
    (hk₂b : P.Base k₂ = (t₁'.baseEquiv P).inv ≫ w ≫ (t₂'.baseEquiv P).hom) :
    k₁ ≫ cmpHom P t₂ t₂' h₂ = cmpHom P t₁ t₁' h₁ ≫ k₂ := by
  refine pullBack_hom_ext P t₂'.hplb ?_ ?_
  · rw [Category.assoc, cmpHom_comp, hk₁c, Category.assoc, hk₂c, cmpHom_comp]
  · rw [P.Base_comp, P.Base_comp, hk₁b, hk₂b, cmpHom_base, cmpHom_base, cmpBase, cmpBase]
    simp only [Category.assoc, Iso.hom_inv_id_assoc]

/-- ★比較同型を `𝒞^pl-bk` の同型へ持ち上げる。 -/
noncomputable def cmpIsoPlBk {φ φ' : A ⟶ B} (t : ArbFac P φ) (t' : ArbFac P φ') (h : φ = φ') :
    (⟨t.top⟩ : PlBk P) ≅ ⟨t'.top⟩ where
  hom := ⟨cmpHom P t t' h, cmpHom_isPullBack P t t' h⟩
  inv := ⟨cmpHom P t' t h.symm, cmpHom_isPullBack P t' t h.symm⟩
  hom_inv_id := WideSubcategory.hom_ext _ (cmpIso P t t' h).hom_inv_id
  inv_hom_id := WideSubcategory.hom_ext _ (cmpIso P t t' h).inv_hom_id

def cmpIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 25,
    item := "Definition 1.3, (iv), (a) — 3 分解の pull-back 部分の canonical 比較",
    sectionId := "frdi-def-1-3-iva" }

end Compare

/-! ## ★★★★「はしご」の等式 —— `Ψ` との両立の要 -/

/-- ★★★★★**はしごの等式** ——
`frb_Z ≫ pre_Z ≫ slPushHom u = u ≫ frb_W ≫ pre_W`。

★★これも **pull-back の一意性だけ**から出る。
★★★これが「`Ψ` を当てても関係が保たれる」ことの土台である ——
`P₂.Base (Ψ g)` を `P.Base g` で書く手段は無い(それこそ `Ψ_Base` である)が、
**`𝒞₁` の中の射の等式**なら `Ψ` を当てるだけで移る。 -/
theorem slPush_ladder (f : A ⟶ A') {Z W : Over (⟨A⟩ : PlBk P)} (u : Z ⟶ W) :
    (slPushFac P F f Z).frb ≫ (slPushFac P F f Z).pre ≫ slPushHom P F f u
      = u.left.hom ≫ (slPushFac P F f W).frb ≫ (slPushFac P F f W).pre := by
  have hu : u.left.hom ≫ W.hom.hom = Z.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (Over.w u)
  refine pullBack_hom_ext P (slPushFac P F f W).hplb ?_ ?_
  · have e1 : ((slPushFac P F f Z).frb ≫ (slPushFac P F f Z).pre ≫ slPushHom P F f u)
        ≫ (slPushFac P F f W).plb = Z.hom.hom ≫ f := by
      rw [Category.assoc, Category.assoc, slPushHom_comp]
      exact (slPushFac P F f Z).fac.symm
    have e2 : (u.left.hom ≫ (slPushFac P F f W).frb ≫ (slPushFac P F f W).pre)
        ≫ (slPushFac P F f W).plb = Z.hom.hom ≫ f := by
      rw [Category.assoc, Category.assoc, ← (slPushFac P F f W).fac, ← Category.assoc, hu]
    exact e1.trans e2.symm
  · have hL : P.Base ((slPushFac P F f Z).frb ≫ (slPushFac P F f Z).pre
          ≫ slPushHom P F f u)
        = ((slPushFac P F f Z).baseEquiv P).hom ≫ slPushBase P F f u := by
      rw [← Category.assoc, P.Base_comp, slPushHom_base]
      rfl
    have hR : P.Base (u.left.hom ≫ (slPushFac P F f W).frb ≫ (slPushFac P F f W).pre)
        = P.Base u.left.hom ≫ ((slPushFac P F f W).baseEquiv P).hom := by
      rw [P.Base_comp]
      rfl
    rw [hL, hR]
    show ((slPushFac P F f Z).baseEquiv P).hom
        ≫ (((slPushFac P F f Z).baseEquiv P).inv ≫ P.Base u.left.hom
          ≫ ((slPushFac P F f W).baseEquiv P).hom)
      = P.Base u.left.hom ≫ ((slPushFac P F f W).baseEquiv P).hom
    exact Iso.hom_inv_id_assoc _ _

/-! ## ★★★★★`Ψ` とスライス関手の 1-可換図式

原文 (FrdI p.68):
> back morphisms, and factorizations as in Deﬁnition 1.3, (iv), (a), it follows that Ψ

★★`Ψ` は 3 分解を 3 分解へ移すので、`Ψ` を当てた分解と
`𝒞₂` 側で選んだ分解が **canonical に比較できる**(`cmpIso`)。 -/

section Psi

variable {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂)
  (Ψ : C ⥤ C₂)
  (hPB : ∀ {X Y : C} (g : X ⟶ Y), IsPullBack P g → IsPullBack P₂ (Ψ.map g))
  (hFT : ∀ {X Y : C} (g : X ⟶ Y), IsFrobeniusType P g → IsFrobeniusType P₂ (Ψ.map g))
  (hPS : ∀ {X Y : C} (g : X ⟶ Y), IsPreStep P g → IsPreStep P₂ (Ψ.map g))

/-- ★★`Ψ` は `Definition 1.3, (iv), (a)` の 3 分解を移す。 -/
def ArbFac.mapPsi {X Y : C} {φ : X ⟶ Y} (t : ArbFac P φ) : ArbFac P₂ (Ψ.map φ) where
  mid := Ψ.obj t.mid
  top := Ψ.obj t.top
  frb := Ψ.map t.frb
  pre := Ψ.map t.pre
  plb := Ψ.map t.plb
  fac := by
    rw [← Ψ.map_comp, ← Ψ.map_comp]
    exact congrArg Ψ.map t.fac
  hfrb := hFT _ t.hfrb
  hpre := hPS _ t.hpre
  hplb := hPB _ t.hplb

/-- ★左辺の分解(`Ψ` を当てたもの)。 -/
noncomputable abbrev psiFacL (f : A ⟶ A') (Z : Over (⟨A⟩ : PlBk P)) :
    ArbFac P₂ (Ψ.map (Z.hom.hom ≫ f)) :=
  (slPushFac P F f Z).mapPsi P P₂ Ψ hPB hFT hPS

/-- ★右辺の分解(`𝒞₂` 側で選んだもの)。 -/
noncomputable abbrev psiFacR (f : A ⟶ A') (Z : Over (⟨A⟩ : PlBk P)) :
    ArbFac P₂ (Ψ.map Z.hom.hom ≫ Ψ.map f) :=
  slPushFac P₂ F₂ (Ψ.map f) ((plBkSlicePsi P P₂ Ψ hPB A).obj Z)

/-- ★★★★★**`Ψ` とスライス関手の 1-可換図式**。 -/
noncomputable def slicePushPsi (f : A ⟶ A') :
    slicePush P F f ⋙ plBkSlicePsi P P₂ Ψ hPB A'
      ≅ plBkSlicePsi P P₂ Ψ hPB A ⋙ slicePush P₂ F₂ (Ψ.map f) :=
  NatIso.ofComponents
    (fun Z => Over.isoMk
      (cmpIsoPlBk P₂ (psiFacL P F P₂ Ψ hPB hFT hPS f Z)
        (psiFacR P P₂ F₂ Ψ hPB f Z) (Ψ.map_comp _ _))
      (WideSubcategory.hom_ext _
        (cmpIso_comp P₂ (psiFacL P F P₂ Ψ hPB hFT hPS f Z)
          (psiFacR P P₂ F₂ Ψ hPB f Z) (Ψ.map_comp _ _))))
    (fun {Z W} u => by
      refine Over.OverMorphism.ext (WideSubcategory.hom_ext _ ?_)
      have hk1c : Ψ.map (slPushHom P F f u) ≫ Ψ.map (slPushFac P F f W).plb
          = Ψ.map (slPushFac P F f Z).plb := by
        rw [← Ψ.map_comp, slPushHom_comp]
      have key : ((psiFacL P F P₂ Ψ hPB hFT hPS f Z).baseEquiv P₂).hom
            ≫ P₂.Base (Ψ.map (slPushHom P F f u))
          = P₂.Base (Ψ.map u.left.hom)
            ≫ ((psiFacL P F P₂ Ψ hPB hFT hPS f W).baseEquiv P₂).hom := by
        show P₂.Base (Ψ.map (slPushFac P F f Z).frb ≫ Ψ.map (slPushFac P F f Z).pre)
            ≫ P₂.Base (Ψ.map (slPushHom P F f u))
          = P₂.Base (Ψ.map u.left.hom)
            ≫ P₂.Base (Ψ.map (slPushFac P F f W).frb ≫ Ψ.map (slPushFac P F f W).pre)
        rw [← P₂.Base_comp, ← P₂.Base_comp, ← Ψ.map_comp, ← Ψ.map_comp, ← Ψ.map_comp,
          ← Ψ.map_comp, Category.assoc]
        exact congrArg (fun g => P₂.Base (Ψ.map g)) (slPush_ladder P F f u)
      have hk1b : P₂.Base (Ψ.map (slPushHom P F f u))
          = ((psiFacL P F P₂ Ψ hPB hFT hPS f Z).baseEquiv P₂).inv
            ≫ P₂.Base (Ψ.map u.left.hom)
            ≫ ((psiFacL P F P₂ Ψ hPB hFT hPS f W).baseEquiv P₂).hom :=
        (Iso.eq_inv_comp _).mpr key
      exact cmpHom_natural P₂ (psiFacL P F P₂ Ψ hPB hFT hPS f Z)
        (psiFacR P P₂ F₂ Ψ hPB f Z) (psiFacL P F P₂ Ψ hPB hFT hPS f W)
        (psiFacR P P₂ F₂ Ψ hPB f W) (Ψ.map_comp _ _) (Ψ.map_comp _ _)
        (P₂.Base (Ψ.map u.left.hom))
        (Ψ.map (slPushHom P F f u)) hk1c hk1b
        (slPushHom P₂ F₂ (Ψ.map f) ((plBkSlicePsi P P₂ Ψ hPB A).map u))
        (slPushHom_comp P₂ F₂ (Ψ.map f) ((plBkSlicePsi P P₂ Ψ hPB A).map u))
        (slPushHom_base P₂ F₂ (Ψ.map f) ((plBkSlicePsi P P₂ Ψ hPB A).map u)))

def slicePushPsi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Ψ とスライス関手の 1-可換図式",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★共役関手の自然性 —— `Ψ_Base` の中核 -/

/-- ★★★★★**共役関手 `basePsiSlice` は `Over.map (Base f)` と可換**。

★★これが原典の `R_i` の「新しい関手」の 1-可換性にあたる。
★3 つの四角形(底・`Ψ`・`𝒞₂` 側の底)を貼り合わせ、
`Definition 1.3, (i), (c)` の圏同値で前合成を簡約する。 -/
noncomputable def basePsiSliceSquare (f : A ⟶ A') :
    Over.map (P.Base f) ⋙ basePsiSlice P F P₂ Ψ hPB A'
      ≅ basePsiSlice P F P₂ Ψ hPB A ⋙ Over.map (P₂.Base (Ψ.map f)) :=
  haveI := F.plBkEquiv A
  haveI := F.plBkEquiv A'
  cancelLeftEquiv (plBkOverFunctor P A).asEquivalence
    (Functor.isoWhiskerRight (slicePushSquare P F f).symm (basePsiSlice P F P₂ Ψ hPB A')
      ≪≫ Functor.isoWhiskerLeft (slicePush P F f)
          (Equivalence.funInvIdAssoc (plBkOverFunctor P A').asEquivalence
            (plBkSlicePsi P P₂ Ψ hPB A' ⋙ plBkOverFunctor P₂ (Ψ.obj A')))
      ≪≫ Functor.isoWhiskerRight (slicePushPsi P F P₂ F₂ Ψ hPB hFT hPS f)
          (plBkOverFunctor P₂ (Ψ.obj A'))
      ≪≫ Functor.isoWhiskerLeft (plBkSlicePsi P P₂ Ψ hPB A)
          (slicePushSquare P₂ F₂ (Ψ.map f))
      ≪≫ (Equivalence.funInvIdAssoc (plBkOverFunctor P A).asEquivalence
            (plBkSlicePsi P P₂ Ψ hPB A ⋙ plBkOverFunctor P₂ (Ψ.obj A)
              ⋙ Over.map (P₂.Base (Ψ.map f)))).symm)

def basePsiSliceSquare.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — 共役関手と Over.map の 1-可換図式",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★★`Ψ_Base` の射 —— `Base` の任意の射を実現する -/

/-- ★★**`ξ` が `φD` を実現する** —— 共役関手の四角形が立つこと。 -/
def RealizesBase {X Y : C} (φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base)
    (ξ : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base) : Prop :=
  Nonempty (Over.map φD ⋙ basePsiSlice P F P₂ Ψ hPB Y
    ≅ basePsiSlice P F P₂ Ψ hPB X ⋙ Over.map ξ)

include F₂ hFT hPS in
/-- ★`Base f` は `Base (Ψ f)` で実現される。 -/
theorem realizesBase_of_map {X Y : C} (f : X ⟶ Y) :
    RealizesBase P F P₂ Ψ hPB (P.Base f) (P₂.Base (Ψ.map f)) :=
  ⟨basePsiSliceSquare P F P₂ F₂ Ψ hPB hFT hPS f⟩

/-- ★★実現は合成で閉じる。 -/
theorem realizesBase_comp {X Y Z : C}
    {φ₁ : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base}
    {φ₂ : (P.toElem.obj Y).base ⟶ (P.toElem.obj Z).base}
    {ξ₁ : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base}
    {ξ₂ : (P₂.toElem.obj (Ψ.obj Y)).base ⟶ (P₂.toElem.obj (Ψ.obj Z)).base}
    (h₁ : RealizesBase P F P₂ Ψ hPB φ₁ ξ₁) (h₂ : RealizesBase P F P₂ Ψ hPB φ₂ ξ₂) :
    RealizesBase P F P₂ Ψ hPB (φ₁ ≫ φ₂) (ξ₁ ≫ ξ₂) := by
  obtain ⟨e₁⟩ := h₁
  obtain ⟨e₂⟩ := h₂
  exact ⟨Functor.isoWhiskerRight (Over.mapComp φ₁ φ₂) (basePsiSlice P F P₂ Ψ hPB Z)
    ≪≫ Functor.isoWhiskerLeft (Over.map φ₁) e₂
    ≪≫ Functor.isoWhiskerRight e₁ (Over.map ξ₂)
    ≪≫ Functor.isoWhiskerLeft (basePsiSlice P F P₂ Ψ hPB X) (Over.mapComp ξ₁ ξ₂).symm⟩

/-- ★`Over.map φ ⋙ Over.map (inv φ) ≅ 𝟭`。 -/
noncomputable def overMapInvIso {Y₁ Y₂ : D} (φ : Y₁ ⟶ Y₂) [IsIso φ] :
    Over.map φ ⋙ Over.map (inv φ) ≅ 𝟭 (Over Y₁) :=
  (Over.mapComp φ (inv φ)).symm ≪≫ eqToIso (congrArg Over.map (IsIso.hom_inv_id φ))
    ≪≫ Over.mapId _

/-- ★`Over.map φ ⋙ Over.map (inv φ) ≅ 𝟭`(`𝒟₂` 側)。 -/
noncomputable def overMapInvIso₂ {Y₁ Y₂ : D₂} (φ : Y₁ ⟶ Y₂) [IsIso φ] :
    Over.map φ ⋙ Over.map (inv φ) ≅ 𝟭 (Over Y₁) :=
  (Over.mapComp φ (inv φ)).symm ≪≫ eqToIso (congrArg Over.map (IsIso.hom_inv_id φ))
    ≪≫ Over.mapId _

set_option maxHeartbeats 2000000 in
/-- ★★★実現は逆射でも成り立つ。 -/
theorem realizesBase_inv {X Y : C}
    (φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base) [IsIso φD]
    (ξ : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base) [IsIso ξ]
    (h : RealizesBase P F P₂ Ψ hPB φD ξ) :
    RealizesBase P F P₂ Ψ hPB (inv φD) (inv ξ) := by
  obtain ⟨e⟩ := h
  refine ⟨cancelLeftEquiv (Over.map φD).asEquivalence ?_⟩
  refine (Functor.isoWhiskerRight (overMapInvIso φD) (basePsiSlice P F P₂ Ψ hPB X)
      ≪≫ (basePsiSlice P F P₂ Ψ hPB X).leftUnitor) ≪≫ ?_
  refine ((Functor.isoWhiskerRight e (Over.map (inv ξ))
      ≪≫ Functor.isoWhiskerLeft (basePsiSlice P F P₂ Ψ hPB X) (overMapInvIso₂ ξ)
      ≪≫ (basePsiSlice P F P₂ Ψ hPB X).rightUnitor)).symm

include F₂ hFT hPS in
set_option maxHeartbeats 2000000 in
/-- ★★★★★★**`Base` の任意の射は実現される**(原文 p.68 の 3 分解)。 -/
theorem exists_realizesBase (φD : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ ξ, RealizesBase P F P₂ Ψ hPB φD ξ := by
  obtain ⟨X, B, α, γ, ψ, hα, hγ, hψ, hfac⟩ := base_three_factor P F φD
  haveI : IsIso (P.Base α) := hα.2
  haveI : IsIso (P₂.Base (Ψ.map α)) := (hPS α hα).2
  refine ⟨inv (P₂.Base (Ψ.map α)) ≫ P₂.Base (Ψ.map γ) ≫ P₂.Base (Ψ.map ψ), ?_⟩
  rw [hfac]
  exact realizesBase_comp P F P₂ Ψ hPB
    (realizesBase_inv P F P₂ Ψ hPB (P.Base α) (P₂.Base (Ψ.map α))
      (realizesBase_of_map P F P₂ F₂ Ψ hPB hFT hPS α))
    (realizesBase_comp P F P₂ Ψ hPB
      (realizesBase_of_map P F P₂ F₂ Ψ hPB hFT hPS γ)
      (realizesBase_of_map P F P₂ F₂ Ψ hPB hFT hPS ψ))

def exists_realizesBase.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Base の任意の射は Ψ 側で実現される",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★★一意性と `psiBaseHom` -/

include F₂ in
set_option maxHeartbeats 2000000 in
/-- ★★★★★**実現する `ξ` は一意** —— `𝒟₂` の slim 性(`Proposition A.2`)から。 -/
theorem realizesBase_unique [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f)
    {X Y : C} {φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base}
    {ξ ξ' : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base}
    (h : RealizesBase P F P₂ Ψ hPB φD ξ) (h' : RealizesBase P F P₂ Ψ hPB φD ξ') :
    ξ = ξ' := by
  obtain ⟨e⟩ := h
  obtain ⟨e'⟩ := h'
  haveI := basePsiSlice_isEquivalence P F P₂ Ψ hPB F₂ hPB' X
  exact overMap_injective hslim₂
    (cancelLeftEquiv (basePsiSlice P F P₂ Ψ hPB X).asEquivalence (e.symm ≪≫ e'))

/-- ★恒等射は恒等射で実現される。 -/
theorem realizesBase_id (X : C) :
    RealizesBase P F P₂ Ψ hPB (𝟙 ((P.toElem.obj X).base))
      (𝟙 ((P₂.toElem.obj (Ψ.obj X)).base)) :=
  ⟨Functor.isoWhiskerRight (Over.mapId _) (basePsiSlice P F P₂ Ψ hPB X)
    ≪≫ (basePsiSlice P F P₂ Ψ hPB X).leftUnitor
    ≪≫ (basePsiSlice P F P₂ Ψ hPB X).rightUnitor.symm
    ≪≫ Functor.isoWhiskerLeft (basePsiSlice P F P₂ Ψ hPB X) (Over.mapId _).symm⟩

include F₂ hFT hPS in
/-- ★★★★★★**`Ψ_Base` の射** —— 選択で取り出すが、一意性で決まっている。 -/
noncomputable def psiBaseHom {X Y : C}
    (φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base) :
    (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base :=
  (exists_realizesBase P F P₂ F₂ Ψ hPB hFT hPS φD).choose

include F₂ hFT hPS in
theorem psiBaseHom_realizes {X Y : C}
    (φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base) :
    RealizesBase P F P₂ Ψ hPB φD (psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS φD) :=
  (exists_realizesBase P F P₂ F₂ Ψ hPB hFT hPS φD).choose_spec

include F₂ hFT hPS in
/-- ★実現する射は `psiBaseHom` に一致する。 -/
theorem psiBaseHom_eq [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f)
    {X Y : C} {φD : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base}
    {ξ : (P₂.toElem.obj (Ψ.obj X)).base ⟶ (P₂.toElem.obj (Ψ.obj Y)).base}
    (h : RealizesBase P F P₂ Ψ hPB φD ξ) :
    psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS φD = ξ :=
  realizesBase_unique P F P₂ F₂ Ψ hPB hslim₂ hPB'
    (psiBaseHom_realizes P F P₂ F₂ Ψ hPB hFT hPS φD) h

include F₂ hFT hPS in
/-- ★★`psiBaseHom` は恒等射を保つ。 -/
theorem psiBaseHom_id [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) (X : C) :
    psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS (𝟙 ((P.toElem.obj X).base))
      = 𝟙 ((P₂.toElem.obj (Ψ.obj X)).base) :=
  psiBaseHom_eq P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'
    (realizesBase_id P F P₂ Ψ hPB X)

include F₂ hFT hPS in
/-- ★★`psiBaseHom` は合成を保つ。 -/
theorem psiBaseHom_comp [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f)
    {X Y Z : C} (φ₁ : (P.toElem.obj X).base ⟶ (P.toElem.obj Y).base)
    (φ₂ : (P.toElem.obj Y).base ⟶ (P.toElem.obj Z).base) :
    psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS (φ₁ ≫ φ₂)
      = psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS φ₁
        ≫ psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS φ₂ :=
  psiBaseHom_eq P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'
    (realizesBase_comp P F P₂ Ψ hPB
      (psiBaseHom_realizes P F P₂ F₂ Ψ hPB hFT hPS φ₁)
      (psiBaseHom_realizes P F P₂ F₂ Ψ hPB hFT hPS φ₂))

include F₂ hFT hPS in
/-- ★★★★★**`psiBaseHom (Base f) = Base (Ψ f)`** —— 1-可換図式の中身。 -/
theorem psiBaseHom_base [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f)
    {X Y : C} (f : X ⟶ Y) :
    psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS (P.Base f) = P₂.Base (Ψ.map f) :=
  psiBaseHom_eq P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'
    (realizesBase_of_map P F P₂ F₂ Ψ hPB hFT hPS f)

def psiBaseHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Ψ_Base の射（存在と一意性）",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★★`Ψ_Base : 𝒟₁ ⥤ 𝒟₂` -/

include F₂ hFT hPS in
/-- ★★★★★★**[FrdI] Theorem 3.4, (v)** —— `Ψ_Base : 𝒟₁ ⥤ 𝒟₂` の構成。

★対象は `Definition 1.3, (i), (a)` が選ぶ `𝒞` の対象を経由し、
射は `psiBaseHom`(存在は 3 分解、一意性は `𝒟₂` の slim 性)で送る。 -/
noncomputable def psiBase [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) :
    D ⥤ D₂ where
  obj Y := (P₂.toElem.obj (Ψ.obj (chosenObj P F Y))).base
  map {Y Y'} g := psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
    ((chosenIso P F Y).hom ≫ g ≫ (chosenIso P F Y').inv)
  map_id Y := by
    have h : (chosenIso P F Y).hom ≫ 𝟙 Y ≫ (chosenIso P F Y).inv
        = 𝟙 ((P.toElem.obj (chosenObj P F Y)).base) := by
      rw [Category.id_comp, Iso.hom_inv_id]
    rw [h]
    exact psiBaseHom_id P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB' _
  map_comp {Y Y' Y''} g g' := by
    have h : (chosenIso P F Y).hom ≫ (g ≫ g') ≫ (chosenIso P F Y'').inv
        = ((chosenIso P F Y).hom ≫ g ≫ (chosenIso P F Y').inv)
          ≫ ((chosenIso P F Y').hom ≫ g' ≫ (chosenIso P F Y'').inv) := by
      simp
    rw [h]
    exact psiBaseHom_comp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB' _ _

include F₂ hFT hPS in
/-- ★1-可換図式の成分。 -/
noncomputable def psiBaseSqApp [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) (A : C) :
    (P.proj ⋙ psiBase P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB').obj A
      ≅ (Ψ ⋙ P₂.proj).obj A where
  hom := psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
    ((chosenIso P F ((P.toElem.obj A).base)).hom)
  inv := psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
    ((chosenIso P F ((P.toElem.obj A).base)).inv)
  hom_inv_id := by
    have h := psiBaseHom_comp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'
      ((chosenIso P F ((P.toElem.obj A).base)).hom)
      ((chosenIso P F ((P.toElem.obj A).base)).inv)
    rw [Iso.hom_inv_id, psiBaseHom_id P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'] at h
    exact h.symm
  inv_hom_id := by
    have h := psiBaseHom_comp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'
      ((chosenIso P F ((P.toElem.obj A).base)).inv)
      ((chosenIso P F ((P.toElem.obj A).base)).hom)
    rw [Iso.inv_hom_id, psiBaseHom_id P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB'] at h
    exact h.symm

include F₂ hFT hPS in
set_option maxHeartbeats 2000000 in
/-- ★★★★★★**[FrdI] Theorem 3.4, (v)** —— `Ψ_Base` の **1-可換図式**。 -/
noncomputable def psiBaseSquare [Ψ.IsEquivalence] (hslim₂ : IsSlimCat D₂)
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.map f) → IsPullBack P f) :
    P.proj ⋙ psiBase P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB' ≅ Ψ ⋙ P₂.proj :=
  NatIso.ofComponents (psiBaseSqApp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB')
    (fun {X Y} f => by
      show psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
            ((chosenIso P F ((P.toElem.obj X).base)).hom ≫ P.Base f
              ≫ (chosenIso P F ((P.toElem.obj Y).base)).inv)
          ≫ psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
            ((chosenIso P F ((P.toElem.obj Y).base)).hom)
        = psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS
            ((chosenIso P F ((P.toElem.obj X).base)).hom)
          ≫ P₂.Base (Ψ.map f)
      rw [← psiBaseHom_comp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB',
        ← psiBaseHom_base P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB' f,
        ← psiBaseHom_comp P F P₂ F₂ Ψ hPB hFT hPS hslim₂ hPB']
      refine congrArg (psiBaseHom P F P₂ F₂ Ψ hPB hFT hPS) ?_
      simp)

def psiBase.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 68,
    item := "Theorem 3.4, (v) — Ψ_Base の構成と 1-可換図式",
    sectionId := "frdi-thm-3-4" }

end Psi

end SlicePush

/-! ## ★★★★★★★`Theorem 3.4` の locator

★★★(i)〜(v) がすべて揃った:

| 条 | 中身 | 実装 |
|---|---|---|
| (i) | `Ψ^istr` | `Thm34.lean` |
| (ii) | pre-step・pull-back・base-isomorphism の保存 | `Thm34Pre.lean` |
| (iii) | `Ψ^pf` | `Thm34Pf.lean` |
| (iv) | `𝒪^▷` / `𝒪^×` の保存と `Ψ^un-tr` | `Thm34Quasi.lean` ＋ `Thm34EndBs.lean` ＋ ★`Thm34IvGen.lean` |
| (v) | `Ψ_Base` の構成・1-可換図式・rigidity | ★`Thm34Slice.lean`(本ファイル) ＋ `Thm34VBase.lean` |

原文 (FrdI p.62):
> Theorem 3.4. (Category-theoreticity of the Base and Frobenius De-

原文 (FrdI p.62):
> an equivalence of categories. Then:

原文 (FrdI p.63):
> rows are equivalences of categories]. Finally, each of the composite functors of this
-/

/-- ★★★★★**[FrdI] Theorem 3.4, (v)** —— 条つきの locator。 -/
def thm_3_4_v.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 63, item := "Theorem 3.4, (v)",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★★**[FrdI] Theorem 3.4** —— 条なしの locator。 -/
def thm_3_4.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62, item := "Theorem 3.4",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
