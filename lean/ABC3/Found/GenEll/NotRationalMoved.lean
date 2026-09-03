/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalPointCoord
import ABC3.Meta.Claim

/-!
# 第 1305 ブロック —— **基礎体の点でなければ動かされる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——第 1298 の言い換え（点の言葉で）

第 1298 は `x` 座標が基礎体に無いことを仮定していた。
★本ブロックは**点そのものが基礎体から来ない**という自然な形にする:

> どの `σ` も `P` を固定するなら、両座標は基礎体にあり、
> したがって `P` は基礎体の点の像である。

☆対偶が「基礎体の点でなければ、どれかの `σ` が動かす」である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

/-- ★★★★★★★★★★★★
**基礎体の点でなければどれかの `σ` が動かす**——★**無条件**（第 1305）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`σ` がすべて固定するなら両座標が基礎体にあり、方程式も基礎体で成り立つ。 -/
theorem exists_galPoint_ne_of_not_rational [IsGalois K₀ M]
    (W : WeierstrassCurve K₀) (P : (W.baseChange M).toAffine.Point)
    (h : ¬ ∃ P₀ : W.toAffine.Point, rhPoint (algebraMap K₀ M) W P₀ = P) :
    ∃ σ : M ≃ₐ[K₀] M, galPoint W σ P ≠ P := by
  by_contra hcon
  push_neg at hcon
  refine h ?_
  cases P with
  | zero => exact ⟨0, rfl⟩
  | some x y hxy =>
      have hx : x ∈ Set.range (algebraMap K₀ M) := by
        refine (InfiniteGalois.mem_range_algebraMap_iff_fixed x).2 ?_
        intro σ
        have h1 := congrArg pointCoords (hcon σ)
        rw [pointCoords_galPoint] at h1
        exact congrArg Prod.fst h1
      have hy : y ∈ Set.range (algebraMap K₀ M) := by
        refine (InfiniteGalois.mem_range_algebraMap_iff_fixed y).2 ?_
        intro σ
        have h1 := congrArg pointCoords (hcon σ)
        rw [pointCoords_galPoint] at h1
        exact congrArg Prod.snd h1
      obtain ⟨x₀, hx₀⟩ := hx
      obtain ⟨y₀, hy₀⟩ := hy
      have hns : W.toAffine.Nonsingular x₀ y₀ := by
        refine (W.toAffine.map_nonsingular (f := algebraMap K₀ M)
          (algebraMap K₀ M).injective x₀ y₀).1 ?_
        rw [hx₀, hy₀]
        exact hxy
      refine ⟨Point.some x₀ y₀ hns, ?_⟩
      show Point.some (algebraMap K₀ M x₀) (algebraMap K₀ M y₀) _ = Point.some x y hxy
      simp only [hx₀, hy₀]

/-! ## ★出典の紐付け(`.src`) -/

def exists_galPoint_ne_of_not_rational.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体の点でなければどれかの σ が動かす。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_ne_of_not_rational.needs : List ProofObligation :=
  [ .citation "[mathlib]" "InfiniteGalois.mem_range_algebraMap_iff_fixed"
      (.inMathlib "InfiniteGalois.mem_range_algebraMap_iff_fixed") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1305）**——第 1298 を点の言葉に直した形である。" ++
       "☆`σ` がすべて固定するなら両座標が基礎体にあり、方程式も基礎体で成り立つ。") 2 ]

end ABC3.Found.GenEll
