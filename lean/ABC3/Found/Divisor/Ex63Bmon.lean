/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Datum
import ABC3.Found.Divisor.ArithGp
import Mathlib.CategoryTheory.ConcreteCategory.EpiMono

/-!
# `Example 6.3` の有理関数の単系 `B(L) = L^×`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★`Theorem 5.2` の入力のもう半分

`Ex63Datum.lean` が因子の側(`Φ`)を与えたので、残るのは
**有理関数の単系** `B(L) = L^×` である。`MonoidOn 𝒟` として立てる:

| フィールド | 中身 |
|---|---|
| `functor` | `L ↦ Additive (L^×)`、射は `Units.map` |
| `charInj` | 体の射は単射 ＋ ★**値が群なので char 成分は無料**(`charMap_injective_of_addGroup`) |
| `fsmIso` | `𝒟` が of FSM-type(`finSubOp_isOfFSMType`)なので FSM 射は同型 |

★`B(L)` が **group-like** であることは、値が群であることそのものである。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `unitsHom` / `unitsHom_injective` | 射が誘導する単数群の射 |
| `unitsMonoidFunctor` | `L ↦ Additive (L^×)` の反変関手 |
| `bmonGalois` | ★★★**`B = L^×` は `𝒟` 上の monoid** |
| `bmonGalois_isGroupLike` | ★`B(L)` は group-like |
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI

open scoped NumberField

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]

/-! ## ★1. 射が誘導する単数群の射 -/

/-- ★**射 `f : L → M` が誘導する `L^× → M^×`**(加法的に書いたもの)。 -/
def unitsHom {L M : FinSub F Kbar} (f : L ⟶ M) :
    Additive (L.toIF)ˣ →+ Additive (M.toIF)ˣ :=
  AddMonoidHom.mk'
    (fun x => Additive.ofMul
      (Units.map (FinSub.hom f).toRingHom.toMonoidHom (Additive.toMul x)))
    (fun x y => by
      show Additive.ofMul (Units.map _ (Additive.toMul x * Additive.toMul y))
        = Additive.ofMul (Units.map _ (Additive.toMul x))
          + Additive.ofMul (Units.map _ (Additive.toMul y))
      rw [map_mul]
      rfl)

omit [NumberField F] in
@[simp] theorem unitsHom_apply {L M : FinSub F Kbar} (f : L ⟶ M) (x : Additive (L.toIF)ˣ) :
    Additive.toMul (unitsHom f x)
      = Units.map (FinSub.hom f).toRingHom.toMonoidHom (Additive.toMul x) := rfl

omit [NumberField F] in
theorem unitsHom_injective {L M : FinSub F Kbar} (f : L ⟶ M) :
    Function.Injective (unitsHom f) := by
  intro x y h
  have h' := congrArg Additive.toMul h
  rw [unitsHom_apply, unitsHom_apply] at h'
  refine Additive.toMul.injective (Units.ext ?_)
  exact FinSub.hom_injective f (congrArg Units.val h')

omit [NumberField F] in
theorem unitsHom_id (L : FinSub F Kbar) : unitsHom (𝟙 L) = AddMonoidHom.id _ := by
  ext x
  rfl

omit [NumberField F] in
theorem unitsHom_comp {L M N : FinSub F Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    unitsHom (f ≫ g) = (unitsHom g).comp (unitsHom f) := by
  ext x
  rfl

/-! ## ★2. 反変関手 -/

/-- ★★**`L ↦ Additive (L^×)`** —— `𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫`。 -/
def unitsMonoidFunctor : ((FinSub F Kbar)ᵒᵖ)ᵒᵖ ⥤ AddCommMonCat.{0} where
  obj X := AddCommMonCat.of (Additive (X.unop.unop.toIF)ˣ)
  map {_ _} f := AddCommMonCat.ofHom (unitsHom f.unop.unop)
  map_id X := by
    show AddCommMonCat.ofHom (unitsHom (𝟙 X.unop.unop)) = _
    rw [unitsHom_id]
    rfl
  map_comp {_ _ _} f g := by
    show AddCommMonCat.ofHom (unitsHom (f.unop.unop ≫ g.unop.unop)) = _
    rw [unitsHom_comp]
    rfl

/-! ## ★3. `MonoidOn` -/

variable (F Kbar) [IsGalois F Kbar]

/-- ★★★★**有理関数の単系 `B(L) = L^×` は `𝒟` 上の monoid**。

★`charInj` の第 2 成分(char 写像の単射性)は**値が群なので無料**である。
★`fsmIso` は `𝒟` が of FSM-type であることから(FSM 射は同型)。 -/
noncomputable def bmonGalois : MonoidOn.{0, 0, 0} (FinSub F Kbar)ᵒᵖ where
  functor := unitsMonoidFunctor
  charInj {_ _} α := by
    refine ⟨unitsHom_injective _, ?_⟩
    exact charMap_injective_of_addGroup (unitsHom α.unop)
  fsmIso {A B} α hα := by
    haveI : IsIso α := finSubOp_isOfFSMType B A α hα
    haveI : IsIso α.op := inferInstance
    exact ConcreteCategory.bijective_of_isIso _

omit [NumberField F] in
/-- ★**`B(L)` は group-like** —— 値が群だから。 -/
theorem bmonGalois_isGroupLike (A : (FinSub F Kbar)ᵒᵖ) :
    IsGroupLike ((bmonGalois F Kbar).val A) := by
  show IsGroupLike (Additive (A.unop.toIF)ˣ)
  refine (isGroupLike_iff _).mpr (fun a => ?_)
  exact ⟨⟨a, -a, by simp, by simp⟩, rfl⟩

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の有理関数の単系 `B(L) = L^×`。 -/
def bmonGalois.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 有理関数の単系 B(L) = L^× は 𝒟 上の monoid",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
