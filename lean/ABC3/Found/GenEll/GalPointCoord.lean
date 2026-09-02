/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38SigmaExists
import ABC3.Found.GenEll.PointCoordNatural
import ABC3.Interface.GaloisRep.Representation
import ABC3.Meta.Claim

/-!
# 第 1298 ブロック —— **基礎体に無い座標をもつ点は動かされる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——`h1` の入力（配管 P2）

第 1270（`exists_galTate_ne_of_galPoint`）は「`σ` が動かす位数 `l` の点」を要求する。
★本ブロックはそれを**座標が基礎体に無いこと**から出す:

* `galPoint` は座標に `σ` を当てるだけ
* 基礎体に無い元はどれかの `σ` が動かす（第 1280）

☆`E[l] ⊄ E(K₀)` から `h1` が出る、という道の一段目である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

/-- ★★★★★★**`galPoint` は座標に `σ` を当てるだけ**——★**無条件**（第 1298）。 -/
theorem pointCoords_galPoint (W : WeierstrassCurve K₀) (σ : M ≃ₐ[K₀] M)
    (P : (W.baseChange M).toAffine.Point) :
    pointCoords (galPoint W σ P) = (σ (pointCoords P).1, σ (pointCoords P).2) := by
  cases P with
  | zero =>
      show ((0 : M), (0 : M)) = (σ (0 : M), σ (0 : M))
      rw [map_zero]
  | some x y h => rfl

/-- ★★★★★★★★★★★★
**座標が基礎体に無ければ、どれかの `σ` がその点を動かす**——★**無条件**（第 1298）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1280（基礎体に無い元は動かされる）を座標に当てるだけである。 -/
theorem exists_galPoint_ne_of_coord_not_mem [IsGalois K₀ M]
    (W : WeierstrassCurve K₀) (P : (W.baseChange M).toAffine.Point)
    (h : (pointCoords P).1 ∉ Set.range (algebraMap K₀ M)) :
    ∃ σ : M ≃ₐ[K₀] M, galPoint W σ P ≠ P := by
  obtain ⟨σ, hσ⟩ := exists_algEquiv_move ((pointCoords P).1) h
  refine ⟨σ, ?_⟩
  intro hcon
  refine hσ ?_
  have h1 := congrArg pointCoords hcon
  rw [pointCoords_galPoint] at h1
  exact congrArg Prod.fst h1

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_galPoint.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galPoint は座標に σ を当てるだけ。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_ne_of_coord_not_mem.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標が基礎体に無ければどれかの σ が点を動かす。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_ne_of_coord_not_mem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_algEquiv_move(第 1280、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_algEquiv_move") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1298）**——`h1` の入力（配管 P2）である。" ++
       "☆`E[l] ⊄ E(K₀)` から `h1` が出る、という道の一段目。") 2 ]

end ABC3.Found.GenEll
