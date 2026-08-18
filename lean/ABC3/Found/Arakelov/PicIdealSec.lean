import ABC3.Found.Arakelov.PicSpecWitness
import Mathlib.AlgebraicGeometry.IdealSheaf.Basic

/-!
# Arakelov (B2) 第 148 ブロック —— **イデアル層の切断**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★B2 の第 1 段——イデアル層を**層加群**にする

mathlib の `Scheme.IdealSheafData` は**アフィン開ごとのイデアル**のデータであり、
`SheafOfModules` との接続は**一切無い**(2026-08-18 実測)。

★そこで任意の開集合 `U` に対し

    idealSections D U := {s ∈ Γ(X,U) | ∀ アフィン開 V ≤ U, s|_V ∈ D.ideal V}

と置く。★★これは `Γ(X,U)` のイデアルであり、制限で保たれ、
**アフィン開では元のデータと一致する**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `idealSections` | ★`U` での切断(イデアル) |
| `idealSections_res` | ★制限は切断を保つ |
| `idealSections_affine` | ★★**アフィン開では元のイデアルと一致** |

★★★第 3 のものが要点である——これが無いと `idealSections` が
`IdealSheafData` の情報を保っていることが言えない。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)


/-- ★イデアル層の `U` での切断。 -/
def idealSections (U : X.Opens) : Ideal (Γ(X, U) : Type u) where
  carrier := {s | ∀ (V : X.affineOpens) (h : V.1 ≤ U),
    (X.presheaf.map (homOfLE h).op).hom s ∈ D.ideal V}
  zero_mem' := by intro V h; simp
  add_mem' := by
    intro a b ha hb V h
    simpa using (D.ideal V).add_mem (ha V h) (hb V h)
  smul_mem' := by
    intro c a ha V h
    have := (D.ideal V).smul_mem ((X.presheaf.map (homOfLE h).op).hom c) (ha V h)
    simpa [smul_eq_mul, map_mul] using this

/-- ★制限は切断を保つ。 -/
theorem idealSections_res {U U' : X.Opens} (h : U' ≤ U)
    (s : (Γ(X, U) : Type u)) (hs : s ∈ idealSections D U) :
    (X.presheaf.map (homOfLE h).op).hom s ∈ idealSections D U' := by
  intro V hV
  rw [← CommRingCat.comp_apply, ← Functor.map_comp, ← op_comp]
  exact hs V (le_trans hV h)

/-- ★★アフィン開では元のイデアルと一致する。 -/
theorem idealSections_affine (U : X.affineOpens) :
    idealSections D U.1 = D.ideal U := by
  refine le_antisymm (fun s hs => ?_) (fun s hs V hV => ?_)
  · have := hs U le_rfl
    simpa using this
  · exact D.ideal_le_comap_ideal hV hs


/-! ## ★出典の紐付け(`.src`) -/

def idealSections_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——イデアル層の切断)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
