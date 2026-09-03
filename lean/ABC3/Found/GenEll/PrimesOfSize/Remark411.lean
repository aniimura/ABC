/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.Chebyshev
import ABC3.Found.GenEll.PrimesOfSize.Lemma41

/-!
# PrimesOfSize —— `[GenEll] Remark 4.1.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-! ## ★★[GenEll] Remark 4.1.1 の**第 1 項** —— 条件 (ii) は初等的である

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

★Remark 4.1.1 は 2 つのことを言う:

1. 条件 (ii) を満たす `ϵ, x_ϵ, C_ϵ` を見つけるのは **entirely elementary**
2. 条件 (i) は **素数定理の帰結**である(外部文献 [Edw] p.76 を指すだけ)

★★**本節は 1 を取る。** 2 は `θ(x)/x → 1` そのもので、mathlib に無い
(`Mathlib/NumberTheory/Chebyshev.lean` は Chebyshev 型の評価しか持たない、
2026-08-16 実測)。公開プロジェクト `PrimeNumberTheoremAnd` が持つ。

★★**したがって `Skeleton/GenEll/Section4.lean` の `remark_4_1_1`
(2 項の連言)は本節では閉じない。** `.src` も条つきにしてある。
**半分取れたことを全部取れたと読まない。**
-/

/-- ★★**Remark 4.1.1 の第 1 項** —— 条件 (ii) を満たす `x_ϵ` は必ず取れる。

> 任意の `M` と `ϵ > 0` に対し、`x ≥ x_ϵ` なる全ての `x` で `M·log(x) ≤ ϵ·x`

★原文が "entirely elementary" と書いた中身は、`log` の増大度が
恒等写像より真に小さいこと(`Real.isLittleO_log_id_atTop`)である。
★`M + 1` で割るのは `M = 0` でも通すため。 -/
theorem exists_xeps_cond_ii (M : ℕ) (eps : ℝ) (heps : 0 < eps) :
    ∃ xeps : ℝ, 0 < xeps ∧ ∀ x : ℝ, xeps ≤ x → (M : ℝ) * Real.log x ≤ eps * x := by
  have hc : 0 < eps / (M + 1) := by positivity
  have hbound := Real.isLittleO_log_id_atTop.bound hc
  obtain ⟨x₀, hx₀⟩ := Filter.eventually_atTop.mp hbound
  refine ⟨max x₀ 1, lt_of_lt_of_le zero_lt_one (le_max_right _ _), fun x hx => ?_⟩
  have hx1 : (1 : ℝ) ≤ x := le_trans (le_max_right _ _) hx
  have hxpos : (0 : ℝ) < x := lt_of_lt_of_le zero_lt_one hx1
  have h := hx₀ x (le_trans (le_max_left _ _) hx)
  simp only [Real.norm_eq_abs, id_eq, abs_of_nonneg hxpos.le] at h
  have hlog : Real.log x ≤ eps / (M + 1) * x := le_trans (le_abs_self _) h
  have hM0 : (0 : ℝ) ≤ (M : ℝ) := Nat.cast_nonneg M
  have hMlt : (M : ℝ) < (M : ℝ) + 1 := lt_add_one _
  have hMpos : (0 : ℝ) < (M : ℝ) + 1 := by positivity
  calc (M : ℝ) * Real.log x ≤ (M : ℝ) * (eps / ((M : ℝ) + 1) * x) := by
        exact mul_le_mul_of_nonneg_left hlog hM0
    _ = ((M : ℝ) / ((M : ℝ) + 1)) * (eps * x) := by field_simp
    _ ≤ 1 * (eps * x) := by
        refine mul_le_mul_of_nonneg_right ?_ (by positivity)
        rw [div_le_one hMpos]
        exact le_of_lt hMlt
    _ = eps * x := one_mul _

/-! ## ★出典の紐付け(`.src`) -/

def exists_xeps_cond_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Remark 4.1.1(第 1 項——条件 (ii) の初等性のみ)",
    sectionId := "genell-rem-4-1-1" }

end ABC3.Found.GenEll
