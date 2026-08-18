import ABC3.Found.Arakelov.PicIdealLocal

/-!
# Arakelov (B2) 第 151 ブロック —— ★★★★★★**イデアル層は層加群である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`IdealSheafData → X.Modules` が繋がった

mathlib には `IdealSheafData` と `SheafOfModules` の接続が**一切無かった**
(2026-08-18 実測)。★第 148–151 の 4 ブロックでそれを作った。

| # | 内容 |
|---|---|
| 148 | 切断 `idealSections`(アフィン開では元のデータと一致) |
| 149 | 前層加群 `idealPresheaf` |
| 150 | 局所性 `idealSections_of_local` |
| 151 | ★★層であること、`X.Modules` へ |

## ★★層の条件の筋

`TopCat.Presheaf.isSheaf_iff_isSheafUniqueGluing` に落とし、

| 段 | 内容 |
|---|---|
| 1 | 台の切断は `𝒪_X` で貼れる(`X.sheaf.existsUnique_gluing`) |
| 2 | 貼った切断が条件を満たす(★第 150) |
| 3 | 一意性は `𝒪_X` の一意性から |

★★`rw` ではなく **`refine (… ).2 ?_`** で入る必要があった
——`rw` は `X.ringCatSheaf.obj` の型検査で落ちる。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

theorem isSheaf_idealPresheaf : Presheaf.IsSheaf (Opens.grothendieckTopology X) (idealPresheaf D).presheaf := by
  refine (TopCat.Presheaf.isSheaf_iff_isSheafUniqueGluing _).2 ?_
  intro ι U sf hsf
  have hsf2 : TopCat.Presheaf.IsCompatible X.presheaf U (fun i => (sf i).1) := by
    intro i j
    exact congrArg Subtype.val (hsf i j)
  obtain ⟨s, hs, huniq⟩ := X.sheaf.existsUnique_gluing U (fun i => (sf i).1) hsf2
  have hmem : s ∈ idealSections D (⨆ i, U i) := by
    refine idealSections_of_local D U (fun i => le_iSup U i) le_rfl s (fun i => ?_)
    have h2 : (X.presheaf.map (homOfLE (le_iSup U i)).op).hom s = (sf i).1 := hs i
    rw [h2]
    exact (sf i).2
  refine ⟨⟨s, hmem⟩, fun i => Subtype.ext (hs i), fun t ht => Subtype.ext ?_⟩
  exact huniq t.1 (fun i => congrArg Subtype.val (ht i))


/-- ★★★★★★**イデアル層を層加群として見る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが B2 の土台である。 -/
noncomputable def idealSheaf : X.Modules where
  val := idealPresheaf D
  isSheaf := isSheaf_idealPresheaf D

/-! ## ★出典の紐付け(`.src`) -/

def idealSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——イデアル層は層加群)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
