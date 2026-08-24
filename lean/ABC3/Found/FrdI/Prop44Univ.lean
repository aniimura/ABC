/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 4.4 —— `Hom^birat` の**普遍性**

原文 (FrdI p.80):
> Proposition 4.4. (Birationalizations) Let C be a Frobenioid. Then:

★★在庫の `HomBirat` には `mk` / `exists_rep` / `sound` / `eq_iff`
(帰納極限**への**射と等号判定)しかなく、**極限からの射**(普遍性)が無かった。
本ファイルはそれを置く。

## ★なぜ要るか

`Corollary 4.10` の `Ψ^birat`(`psiBiratCor`)は**手で組んである**
(`biratPsiMap` / `biratPsiMap_id` / `biratPsiMap_comp`)。
★同じことを毎回手で組むのは高くつく ——
`Proposition 5.5, (ii)` の梱包でも同じ形が要る
(`Prop55ScaleRootCoa.lean` の「逃げ道 B」)。

★★普遍性があれば「**co-angular pre-step を同型に送る写像は `Hom^birat` を経由する**」
と 1 行で書ける。`HomBirat = HomColim (homFunctorBirat …)` で
`HomColim` は mathlib の `colimit` なので、`HomColim.desc` がそのまま乗る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BiratUniv

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★遷移射との両立を `HomColim` の形に持ち上げたもの。 -/
theorem HomBirat.descAux {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ)
    {Z W : IdxBirat P G A} (u : Z ⟶ W) (x : (homFunctorBirat P G A B).obj Z) :
    f W ((homFunctorBirat P G A B).map u x).down = f Z x.down := by
  rw [homFunctorBirat_map]
  exact hf u x.down

/-- ★★★★★**`Hom^birat` の普遍性** ——
各添字での写像が遷移射(＝前合成)と両立すれば、`Hom^birat` から降りる。 -/
noncomputable def HomBirat.desc {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ) :
    HomBirat P G A B → T :=
  HomColim.desc (homFunctorBirat P G A B) (fun Z x => f Z x.down)
    (fun {_ _} u x => HomBirat.descAux f hf u x)

@[simp] theorem HomBirat.desc_mk {A B : C} {T : Type (max v u2 v2)}
    (f : ∀ Z : IdxBirat P G A, (Z.unop.left.obj ⟶ B) → T)
    (hf : ∀ {Z W : IdxBirat P G A} (u : Z ⟶ W) (φ : Z.unop.left.obj ⟶ B),
      f W (u.unop.left.hom ≫ φ) = f Z φ)
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    HomBirat.desc f hf (HomBirat.mk Z φ) = f Z φ :=
  HomColim.desc_mk (homFunctorBirat P G A B) (fun Z x => f Z x.down)
    (fun {_ _} u x => HomBirat.descAux f hf u x) Z (ULift.up φ)

/-- ★★**`Hom^birat` からの 2 つの写像は代表元で決まる**。 -/
theorem HomBirat.desc_ext {A B : C} {T : Type (max v u2 v2)}
    (g h : HomBirat P G A B → T)
    (hgh : ∀ (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B),
      g (HomBirat.mk Z φ) = h (HomBirat.mk Z φ)) : g = h :=
  HomColim.desc_ext (homFunctorBirat P G A B) g h (fun Z x => hgh Z x.down)

end BiratUniv

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Hom^birat` の普遍性(帰納極限からの降下)。 -/
def HomBirat.desc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 80,
    item := "Proposition 4.4 — Hom^birat の普遍性(代表元からの降下)",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
