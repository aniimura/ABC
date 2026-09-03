/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalPointCoord
import ABC3.Meta.Claim

/-!
# 第 1299 ブロック —— **基礎体の点は `σ` に固定される**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——`h2` の入力（配管 P1 の要）

第 1296（`h2` の最短形）は「`σ` が固定する位数 `l` の点」を要求する。
★**基礎体 `K₀` の上の点はどの `σ : M ≃ₐ[K₀] M` にも固定される**ので、
`μ_l ⊂ E(K₀)`（第 1297）がそのまま入力になる。

☆座標が `algebraMap K₀ M` の像にあり、`σ` は `K₀` を固定するからである。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

/-- ★★★★★★★★★★★★
**基礎体の点は `σ` に固定される**——★**無条件**（第 1299）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆座標が `algebraMap K₀ M` の像にあり、`σ` は `K₀` を固定する。 -/
theorem galPoint_rhPoint_eq (W : WeierstrassCurve K₀) (σ : M ≃ₐ[K₀] M)
    (P₀ : W.toAffine.Point) :
    galPoint W σ (rhPoint (algebraMap K₀ M) W P₀) = rhPoint (algebraMap K₀ M) W P₀ := by
  cases P₀ with
  | zero =>
      show galPoint W σ 0 = 0
      exact map_zero _
  | some x y h =>
      show galPoint W σ (Point.some (algebraMap K₀ M x) (algebraMap K₀ M y) _)
        = Point.some (algebraMap K₀ M x) (algebraMap K₀ M y) _
      show Point.some (σ (algebraMap K₀ M x)) (σ (algebraMap K₀ M y)) _
        = Point.some (algebraMap K₀ M x) (algebraMap K₀ M y) _
      simp only [AlgEquiv.commutes]

/-! ## ★出典の紐付け(`.src`) -/

def galPoint_rhPoint_eq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体の点は σ に固定される。★無条件)",
    sectionId := "genell-thm-3-8" }

def galPoint_rhPoint_eq.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1299）**——`h2` の入力（配管 P1 の要）である。" ++
       "☆`μ_l ⊂ E(K₀)`（第 1297）がそのまま「`σ` が固定する位数 `l` の点」になる。") 2 ]

end ABC3.Found.GenEll
