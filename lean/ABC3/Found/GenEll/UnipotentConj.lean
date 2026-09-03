/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Notation
import ABC3.Meta.Claim

/-!
# 第 1229 ブロック —— **自明でない冪単行列は `α` に共役**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——`ζ, π` の基底を作らずに済ませる道

`α ∈ (galRep の像).map (glRedPadic l)` を出すには、第 1228 により
「ある基底 `e₀` で `σ` の行列が `α`」が要る。

☆`galRep` は基底の取り替えで**共役**に変わる（`galRep_basisChange`、在庫）ので、
`σ` の `mod l` の行列が `α` に**共役**でありさえすれば、
基底を取り替えて `α` そのものにできる。

★本ブロックはその線型代数の核を取る:

    (M − 1)² = 0 かつ M ≠ 1  ⟹  M は (1 1 / 0 1) に共役

☆証明は `u ≔ (M−1) v`（`v` は `(M−1) v ≠ 0` なる標準基底ベクトル）と `v` が
基底になること——`(M−1) u = 0` と `u ≠ 0` から独立性が出る。

★★★これで **`ζ, π` に合わせた基底 `e₀` を作らなくてよい**。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★★★★★
**自明でない冪単な 2×2 行列は `(1 1 / 0 1)` に共役**——★**無条件**（第 1229）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`u ≔ (M−1) v` と `v` が基底になる（`(M−1) u = 0` と `u ≠ 0` から）。
その基底で `M` は `e₁ ↦ e₁`、`e₂ ↦ e₁ + e₂` として作用する。

★★★これで `α` を出すのに `ζ, π` に合わせた基底を作らなくてよい。 -/
theorem exists_conj_upperOne {F : Type*} [Field F] (M : Matrix (Fin 2) (Fin 2) F)
    (hnil : (M - 1) * (M - 1) = 0) (hne : M ≠ 1) :
    ∃ A : Matrix (Fin 2) (Fin 2) F, A.det ≠ 0 ∧ M * A = A * !![1, 1; 0, 1] := by
  classical
  set N := M - 1 with hNdef
  have hNne : N ≠ 0 := sub_ne_zero.2 hne
  have hMN : M = N + 1 := by simp [hNdef]
  obtain ⟨j, hj⟩ : ∃ j : Fin 2, N.mulVec (Pi.single j 1) ≠ 0 := by
    by_contra hc
    push_neg at hc
    refine hNne (Matrix.ext_of_mulVec_single ?_)
    intro i
    rw [hc i, Matrix.zero_mulVec]
  set v : Fin 2 → F := Pi.single j 1 with hvdef
  set u : Fin 2 → F := N.mulVec v with hudef
  have hNu : N.mulVec u = 0 := by
    rw [hudef, Matrix.mulVec_mulVec, hnil, Matrix.zero_mulVec]
  have hMu : M.mulVec u = u := by
    rw [hMN, Matrix.add_mulVec, hNu, Matrix.one_mulVec, zero_add]
  have hMv : M.mulVec v = u + v := by
    rw [hMN, Matrix.add_mulVec, Matrix.one_mulVec, ← hudef]
  have hindep : ∀ c : F, u ≠ c • v := by
    intro c hc
    have h1 : c • u = 0 := by
      have h2 : N.mulVec (c • v) = c • N.mulVec v := Matrix.mulVec_smul _ _ _
      rw [← hc, hNu] at h2
      exact h2.symm
    have hc0 : c = 0 := by
      rcases smul_eq_zero.mp h1 with h | h
      · exact h
      · exact absurd h hj
    rw [hc0, zero_smul] at hc
    exact hj hc
  have hmu : ∀ i, M i 0 * u 0 + M i 1 * u 1 = u i := by
    intro i
    have h := congrFun hMu i
    simpa [Matrix.mulVec, dotProduct, Fin.sum_univ_two] using h
  have hmv : ∀ i, M i 0 * v 0 + M i 1 * v 1 = u i + v i := by
    intro i
    have h := congrFun hMv i
    simpa [Matrix.mulVec, dotProduct, Fin.sum_univ_two] using h
  refine ⟨!![u 0, v 0; u 1, v 1], ?_, ?_⟩
  · rw [Matrix.det_fin_two_of]
    intro hdet
    fin_cases j
    · have hv0 : v 0 = 1 := by simp [hvdef]
      have hv1 : v 1 = 0 := by simp [hvdef]
      rw [hv0, hv1, mul_zero, one_mul, zero_sub, neg_eq_zero] at hdet
      refine hindep (u 0) (funext fun i => ?_)
      fin_cases i
      · simp [hv0]
      · simp [hv1, hdet]
    · have hv0 : v 0 = 0 := by simp [hvdef]
      have hv1 : v 1 = 1 := by simp [hvdef]
      rw [hv0, hv1, mul_one, zero_mul, sub_zero] at hdet
      refine hindep (u 1) (funext fun i => ?_)
      · fin_cases i
        · simp [hv0, hdet]
        · simp [hv1]
  · refine Matrix.ext fun i k => ?_
    fin_cases i <;> fin_cases k <;>
      simp [Matrix.mul_apply, Fin.sum_univ_two, hmu, hmv]

/-! ## ★出典の紐付け(`.src`) -/

def exists_conj_upperOne.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(自明でない冪単な 2×2 行列は (1 1 / 0 1) に共役。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_conj_upperOne.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1229）**——`galRep` は基底の取り替えで**共役**に" ++
       "変わる（`galRep_basisChange`、在庫）ので、`σ` の `mod l` の行列が `α` に" ++
       "共役でありさえすれば、基底を取り替えて `α` そのものにできる。" ++
       "★★★これで `α` を出すのに **`ζ, π` に合わせた基底 `e₀` を作らなくてよい**。") 3 ]

end ABC3.Found.GenEll
