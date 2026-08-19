/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34VBase

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
  rw [← P.Base_comp, ← Category.assoc, ← t.fac]

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

end SlicePush

end ABC3.Found.FrdI
