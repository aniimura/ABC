/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Theorem

/-!
# six exponentials を当てるための道具

★[FrdI] `Lemma 6.5, (ii)` に `six_exponentials` を当てるとき要る 3 つ。

* `isAlgebraic_exp_rat_mul_log` —— **`exp(μ·log q)` は代数的**(`μ ∈ ℚ`、`q` は正の自然数)。
  `z^{den} = q^{num}` が有理数なので `X^{den} − C(q^{num})` の根である。
* `linearIndependent_ofReal` —— 実数の族の ℚ 上一次独立性は複素数へ移る。
* `linearIndependent_pair_one_real` —— `t ∉ ℚ` なら `{1, t}` は ℚ 上一次独立。
-/

namespace ABC3.Found.SixExp

open Complex

/-- ★`z^n` が有理数なら `z` は代数的。 -/
theorem isAlgebraic_of_pow_eq_rat {z : ℂ} {n : ℕ} (hn : 0 < n) {r : ℚ}
    (h : z ^ n = (r : ℂ)) : IsAlgebraic ℚ z := by
  refine ⟨Polynomial.X ^ n - Polynomial.C r, Polynomial.X_pow_sub_C_ne_zero hn r, ?_⟩
  simp [h]

/-- ★★**`exp(μ·log q)` は代数的**(`q` は正の自然数、`μ` は有理数)。

★`z^{den} = q^{num}` が有理数なので `X^{den} − C(q^{num})` の根。 -/
theorem isAlgebraic_exp_rat_mul_log {q : ℕ} (hq : 0 < q) (μ : ℚ) :
    IsAlgebraic ℚ (Complex.exp ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ)) := by
  set z : ℂ := Complex.exp ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ) with hzdef
  have hqR : (0:ℝ) < ((q : ℕ) : ℝ) := by exact_mod_cast hq
  have hlogq : Complex.exp (((Real.log ((q : ℕ) : ℝ) : ℝ)) : ℂ) = ((q : ℕ) : ℂ) := by
    rw [← Complex.ofReal_exp, Real.exp_log hqR]
    norm_num
  have hdQ : ((μ.den : ℚ)) * μ = (μ.num : ℚ) := by
    have h := Rat.num_div_den μ
    have hden : ((μ.den : ℚ)) ≠ 0 := Nat.cast_ne_zero.mpr μ.den_nz
    field_simp at h
    linarith
  have hdC : ((μ.den : ℂ)) * ((μ : ℚ) : ℂ) = ((μ.num : ℤ) : ℂ) := by
    exact_mod_cast congrArg (fun t : ℚ => (t : ℂ)) hdQ
  refine isAlgebraic_of_pow_eq_rat μ.pos (r := ((q : ℚ)) ^ (μ.num)) ?_
  have hstep : z ^ (μ.den)
      = Complex.exp ((μ.den : ℂ) * ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ)) := by
    rw [hzdef, ← Complex.exp_nat_mul]
  have hnum : (μ.den : ℂ) * ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ)
      = ((μ.num : ℤ) : ℂ) * (((Real.log ((q : ℕ) : ℝ) : ℝ)) : ℂ) := by
    calc (μ.den : ℂ) * ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ)
        = ((μ.den : ℂ) * ((μ : ℚ) : ℂ)) * (((Real.log ((q : ℕ) : ℝ) : ℝ)) : ℂ) := by
          push_cast; ring
      _ = ((μ.num : ℤ) : ℂ) * (((Real.log ((q : ℕ) : ℝ) : ℝ)) : ℂ) := by rw [hdC]
  rw [hstep, hnum, Complex.exp_int_mul, hlogq]
  push_cast
  ring

/-- ★実数の族が ℚ 上一次独立なら、複素数に送っても一次独立。 -/
theorem linearIndependent_ofReal {ι : Type*} {v : ι → ℝ} (h : LinearIndependent ℚ v) :
    LinearIndependent ℚ (fun i => ((v i : ℝ) : ℂ)) := by
  rw [linearIndependent_iff'] at h ⊢
  intro t g hsum i hi
  refine h t g ?_ i hi
  have hcast : ((∑ j ∈ t, g j • v j : ℝ) : ℂ) = ∑ j ∈ t, g j • ((v j : ℝ) : ℂ) := by
    push_cast [Rat.smul_def]
    rfl
  have hz : ((∑ j ∈ t, g j • v j : ℝ) : ℂ) = 0 := by rw [hcast, hsum]
  exact_mod_cast hz

/-- ★`t` が有理数でなければ `{1, t}` は ℚ 上一次独立(実数版)。 -/
theorem linearIndependent_pair_one_real {t : ℝ} (ht : ∀ r : ℚ, t ≠ (r : ℝ)) :
    LinearIndependent ℚ ![(1 : ℝ), t] := by
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have hg' : ((g 0 : ℚ) : ℝ) + ((g 1 : ℚ) : ℝ) * t = 0 := by
    rw [Fin.sum_univ_two] at hg
    simpa [Rat.smul_def] using hg
  have hg1 : g 1 = 0 := by
    by_contra h1
    have hc : ((g 1 : ℚ) : ℝ) ≠ 0 := by simpa using h1
    have hteq : t = ((-g 0 / g 1 : ℚ) : ℝ) := by
      push_cast
      field_simp
      linear_combination hg'
    exact ht _ hteq
  have hg0 : g 0 = 0 := by
    rw [hg1] at hg'
    simp at hg'
    exact_mod_cast hg'
  intro i
  fin_cases i
  · exact hg0
  · exact hg1

end ABC3.Found.SixExp
