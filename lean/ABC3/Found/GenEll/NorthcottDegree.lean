import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings

/-!
# [GenEll] Proposition 1.4, (iv) —— **次数有界**の代数的整数の有限性(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★原文は `X(ℚ̄)^{≤d}` を走る —— 体は固定されていない

`Found/GenEll/DenominatorBound.lean` は**固定した数体 `K`** 上で
`{y : K | H(y) ≤ B}` の有限性を取った。
★しかし原文が走らせるのは `X(ℚ̄)^{≤d}` ——
**`[F:ℚ] ≤ d` なる数体すべての和集合**である。
★★**次数 `≤ d` の数体は無限個ある**(`ℚ(√p)` を考えよ)から、
固定体の結果を並べるだけでは届かない。

## ★本ファイルが取るもの

> **`ℂ` の代数的整数で、次数が `d` 以下かつ共役がすべて `B` 以下のものは有限個**

★mathlib は**固定した数体 `K` の中**での版
(`NumberField.Embeddings.finite_of_norm_le`)しか持たない。
本ファイルは**体を固定しない**版である。

★機構は mathlib の `finite_of_norm_le` と同じ:
`Polynomial.coeff_bdd_of_roots_le`(根が有界なら係数が有界)と
`Polynomial.bUnion_roots_finite`(次数と係数が有界な整数係数多項式の根は有限個)。
★**次数を `[K:ℚ]` ではなく引数 `d` にするだけ**で体の固定が外れる。

## ★★まだ届いていないもの

原文は「**高さ**が `C` 以下」であって「共役が `B` 以下」ではない。
★両者を繋ぐには `α` の定義体 `ℚ(α)` を作り、
そこで `v(α) ≤ H(α)`(`NorthcottNF.lean` の `infinitePlace_le_mulHeight₁`)と
**分母の評価**(`DenominatorBound.lean`)を使う必要がある。
★★**本ファイルはその手前までである。** `.src` は条つきである。
-/

namespace ABC3.Found.GenEll

open Polynomial

/-- ★★**次数有界・共役有界の代数的整数は有限個**(体を固定しない版)。

★mathlib の `NumberField.Embeddings.finite_of_norm_le` は
**固定した数体 `K` の中**でしか言っていない(2026-08-17 実測)。
本定理は `d` を引数にすることで体の固定を外す。 -/
theorem finite_isIntegral_natDegree_le_norm_le (d : ℕ) (B : ℝ) :
    {α : ℂ | IsIntegral ℤ α ∧ (minpoly ℤ α).natDegree ≤ d ∧
      ∀ z ∈ ((minpoly ℤ α).map (algebraMap ℤ ℂ)).roots, ‖z‖ ≤ B}.Finite := by
  classical
  set C : ℕ := ⌈max B 1 ^ d * (d.choose (d / 2) : ℝ)⌉₊ with hC
  refine Set.Finite.subset (Polynomial.bUnion_roots_finite (algebraMap ℤ ℂ) d
    (Set.finite_Icc (-(C : ℤ)) (C : ℤ))) ?_
  rintro α ⟨hint, hdeg, hroots⟩
  simp only [Set.mem_iUnion]
  refine ⟨minpoly ℤ α, ⟨hdeg, fun i => ?_⟩, ?_⟩
  · -- ★係数の評価
    have hbdd := Polynomial.coeff_bdd_of_roots_le (B := B) (d := d) (algebraMap ℤ ℂ)
      (minpoly.monic hint) (IsAlgClosed.splits _) hdeg hroots i
    have h2 : |(((minpoly ℤ α).coeff i : ℤ) : ℝ)| ≤ (C : ℝ) := by
      refine le_trans ?_ (Nat.le_ceil _)
      calc |(((minpoly ℤ α).coeff i : ℤ) : ℝ)|
          = ‖((minpoly ℤ α).map (algebraMap ℤ ℂ)).coeff i‖ := by
            rw [coeff_map, eq_intCast]
            simp
        _ ≤ max B 1 ^ d * (d.choose (d / 2) : ℝ) := hbdd
    rw [Set.mem_Icc, ← abs_le, ← @Int.cast_le ℝ, Int.cast_abs]
    push_cast
    exact h2
  · -- ★`α` はその多項式の根である
    simp only [Finset.mem_coe, Multiset.mem_toFinset, mem_roots']
    refine ⟨((minpoly.monic hint).map (algebraMap ℤ ℂ)).ne_zero, ?_⟩
    have h := minpoly.aeval ℤ α
    rw [aeval_def, algebraMap_int_eq] at h
    rw [IsRoot, eval_map, algebraMap_int_eq]
    exact h

/-! ## ★出典の紐付け(`.src`) -/

def finite_isIntegral_natDegree_le_norm_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数有界の代数的整数の有限性のみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
