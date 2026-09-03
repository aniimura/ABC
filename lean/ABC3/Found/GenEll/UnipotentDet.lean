/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv
import Mathlib.LinearAlgebra.Matrix.Trace
import Mathlib.Tactic.LinearCombination
import ABC3.Meta.Claim

/-!
# 第 1292 ブロック —— **行列式 1 で固定ベクトルがあれば幂単**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——幂単性の**別の道**

第 1289 までの道は Tate 一意化の同変性を使っていたが、
**幂単性はもっと安く出る**:

> `μ_l ⊂ E[l]` は Tate 曲線では**基礎体の上で有理**（`ζ ∈ K₀` だから）。
> したがって `σ` は `E[l]` の中の直線を**点ごとに**固定する。
> 一方 `det`（Galois 表現の行列式）は円分指標であり、`ζ_l ∈ L` なら **1**。
> `2×2` 行列で「固定ベクトルがあり `det = 1`」なら、跡は `2`、
> したがって Cayley–Hamilton より `(M − 1)² = 0`。

★★★これで**局所体の拡大（第 1291 の節点 (f)）は要らなくなる**
——要るのは `K₀ = p.adicCompletion L` の上の Tate 一意化だけである。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

open scoped Matrix

/-- ★★★★★★★★★★★★★★★★★★★★
**`det = 1` で固定ベクトルがあれば `(M − 1)² = 0`**——★**無条件**（第 1292）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`M v = v`（`v ≠ 0`）から `det (M − 1) = 0`、`det M = 1` と合わせて跡は `2`。
★`2×2` なので Cayley–Hamilton は成分計算で足りる。 -/
theorem sq_sub_one_eq_zero_of_det_one_of_fixed {F : Type*} [Field F] [DecidableEq F]
    (M : Matrix (Fin 2) (Fin 2) F) (hdet : M.det = 1)
    (v : Fin 2 → F) (hv : v ≠ 0) (hMv : M *ᵥ v = v) :
    (M - 1) * (M - 1) = 0 := by
  have hker : (M - 1) *ᵥ v = 0 := by
    rw [Matrix.sub_mulVec, hMv, Matrix.one_mulVec, sub_self]
  have hdet1 : (M - 1).det = 0 := Matrix.exists_mulVec_eq_zero_iff.1 ⟨v, hv, hker⟩
  have h1 : (M - 1).det = (M 0 0 - 1) * (M 1 1 - 1) - M 0 1 * M 1 0 := by
    rw [Matrix.det_fin_two]
    simp [Matrix.sub_apply, Matrix.one_apply]
  have h2 : M.det = M 0 0 * M 1 1 - M 0 1 * M 1 0 := Matrix.det_fin_two M
  rw [h1] at hdet1
  rw [h2] at hdet
  have htr : M 0 0 + M 1 1 = 2 := by linear_combination hdet - hdet1
  have hbc : M 0 1 * M 1 0 = (M 0 0 - 1) * (M 1 1 - 1) := by linear_combination - hdet1
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp only [Matrix.mul_apply, Fin.sum_univ_two, Matrix.sub_apply, Matrix.one_apply_eq,
      Matrix.one_apply_ne, Matrix.zero_apply, Fin.isValue, ne_eq, Fin.zero_eta, Fin.mk_one,
      one_ne_zero, not_false_eq_true, zero_ne_one]
  · linear_combination hbc + (M 0 0 - 1) * htr
  · linear_combination (M 0 1) * htr
  · linear_combination (M 1 0) * htr
  · linear_combination hbc + (M 1 1 - 1) * htr

/-! ## ★出典の紐付け(`.src`) -/

def sq_sub_one_eq_zero_of_det_one_of_fixed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(det = 1 で固定ベクトルがあれば (M − 1)² = 0。★無条件)",
    sectionId := "genell-thm-3-8" }

def sq_sub_one_eq_zero_of_det_one_of_fixed.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1292）**——幂単性の**別の道**である。" ++
       "☆`μ_l ⊂ E[l]` は Tate 曲線では基礎体の上で有理（`ζ ∈ K₀`）なので、" ++
       "`σ` は直線を点ごとに固定する。★`det` は円分指標であり `ζ_l ∈ L` なら `1`。" ++
       "★★★これで**局所体の拡大（第 1291 の節点 (f)）は要らなくなる**" ++
       "——要るのは `K₀ = p.adicCompletion L` の上の Tate 一意化だけである。") 3 ]

end ABC3.Found.GenEll
