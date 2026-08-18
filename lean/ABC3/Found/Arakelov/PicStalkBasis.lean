import ABC3.Found.Arakelov.PicAwayIso

/-!
# Arakelov (B1) 第 90 ブロック —— **基底で同型なら茎で同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 89 と第 77 を繋ぐ器具

第 89 ブロックで「比較射は**基本開集合の上で**全単射」が出た。
★第 77 ブロックの器具は「**茎で**同型なら同型」である。
★★両者を繋ぐのが本ブロック——**基底で全単射なら茎で全単射**。

## ★★機構 —— mathlib の在庫 + 1 手

| 向き | 在庫 |
|---|---|
| 単射 | ★mathlib `stalkFunctor_map_injective_of_isBasis` |
| 全射 | ★★**本ブロックで建てる**(`exists_mem_germ_eq_of_isBasis` から一発) |

★★★全射側は「茎の元は基底の開集合の切断の芽で書ける」(`exists_mem_germ_eq_of_isBasis`)
を使うだけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `stalkFunctor_map_surjective_of_isBasis` | ★★基底で全射 ⟹ 茎で全射 |
| `stalkFunctor_map_bijective_of_isBasis` | ★★★★**基底で全単射 ⟹ 茎で全単射** |
| `isIso_stalkFunctor_map_of_isBasis` | ★★★★**基底で全単射 ⟹ 茎で同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory TopologicalSpace Opposite TopCat.Presheaf

variable {X : TopCat.{u}} {B : Set (Opens X)} (hB : Opens.IsBasis B)

include hB

/-- ★★**基底の上で全射なら茎で全射**。 -/
theorem stalkFunctor_map_surjective_of_isBasis {F G : X.Presheaf AddCommGrpCat.{u}} {α : F ⟶ G}
    (hα : ∀ U ∈ B, Function.Surjective (α.app (op U))) (x : X) :
    Function.Surjective ((stalkFunctor _ x).map α) := by
  intro t
  obtain ⟨U, hxU, hU, s, rfl⟩ := exists_mem_germ_eq_of_isBasis hB G x t
  obtain ⟨s', rfl⟩ := hα U hU s
  exact ⟨F.germ U x hxU s', by rw [stalkFunctor_map_germ_apply]⟩

/-- ★★★★**基底の上で全単射なら茎で全単射**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが第 89 ブロック(基本開集合で全単射)と
第 77 ブロック(茎で同型なら同型)を繋ぐ。 -/
theorem stalkFunctor_map_bijective_of_isBasis {F G : X.Presheaf AddCommGrpCat.{u}} {α : F ⟶ G}
    (hα : ∀ U ∈ B, Function.Bijective (α.app (op U))) (x : X) :
    Function.Bijective ((stalkFunctor _ x).map α) :=
  ⟨stalkFunctor_map_injective_of_isBasis hB (fun U hU => (hα U hU).1) x,
    stalkFunctor_map_surjective_of_isBasis hB (fun U hU => (hα U hU).2) x⟩

/-- ★★★★**基底の上で全単射なら茎で同型**。 -/
theorem isIso_stalkFunctor_map_of_isBasis {F G : X.Presheaf AddCommGrpCat.{u}} {α : F ⟶ G}
    (hα : ∀ U ∈ B, Function.Bijective (α.app (op U))) (x : X) :
    IsIso ((stalkFunctor _ x).map α) :=
  (ConcreteCategory.isIso_iff_bijective _).2 (stalkFunctor_map_bijective_of_isBasis hB hα x)

/-! ## ★出典の紐付け(`.src`) -/

def stalkFunctor_map_bijective_of_isBasis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基底で全単射なら茎で全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
