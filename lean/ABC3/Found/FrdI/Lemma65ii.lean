/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Algebraic
import ABC3.Found.FrdI.Lemma65

/-!
# [FrdI] Lemma 6.5, (ii) —— 6 素数の対数の比の関係式は不可能

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116–117。

> (ii) Let `p₁, p₂, …, p₆` be distinct prime numbers. Then there do not exist
> `λ₁, λ₂ ∈ ℚ>0` such that:
> `log(p₁)/log(p₂) = λ₁·log(p₃)/log(p₄) = λ₂·log(p₅)/log(p₆)`.

原文の証明:

> Assertion (ii) is a consequence of a theorem of Lang: Indeed, since the `log(pᵢ)`
> are linearly independent over ℚ [by assertion (i)], it follows that each of the
> following two sets of numbers is also linearly independent over ℚ:
> `{log(p₂), log(p₄), log(p₆)}`; `{1, log(p₃)/log(p₄)}`.
> Moreover, all six products of one element from the first set and one element from
> the second set are of the form `μ·log(pᵢ)`, where `μ ∈ ℚ>0`. Thus, the exponential
> of each of these products is algebraic, in contradiction to Lang's theorem.

★**Lang の定理 = six exponentials theorem** で、`Found/SixExp/Theorem.lean` に
`six_exponentials` として**証明済み**である。本ファイルはそれを当てる。

## ★6 つの積

`x = {1, t}`(`t = log p₃/log p₄`)、`y = {log p₂, log p₄, log p₆}` として

| 積 | `μ` | `pᵢ` |
|---|---|---|
| `1·log p₂` | `1` | `p₂` |
| `1·log p₄` | `1` | `p₄` |
| `1·log p₆` | `1` | `p₆` |
| `t·log p₄` | `1` | `p₃` |
| `t·log p₂` | `1/λ₁` | `p₁` |
| `t·log p₆` | `λ₂/λ₁` | `p₅` |

★後ろ 2 つが仮定の関係式を使う所である。
-/

namespace ABC3.Found.FrdI

open ABC3.Found.SixExp

/-- ★★相異なる素数の対数の比は有理数でない。 -/
theorem log_prime_ne_rat_mul {p q : Nat.Primes} (hpq : p ≠ q) (r : ℚ) :
    Real.log ((p : ℕ) : ℝ) ≠ (r : ℝ) * Real.log ((q : ℕ) : ℝ) := by
  classical
  intro h
  have hli := log_primes_linearIndependent
  rw [linearIndependent_iff'] at hli
  set g : Nat.Primes → ℚ := fun z => if z = p then 1 else if z = q then -r else 0 with hg
  have hgp : g p = 1 := by simp [hg]
  have hgq : g q = -r := by simp [hg, hpq.symm]
  have hsum : ∑ z ∈ ({p, q} : Finset Nat.Primes), g z • Real.log ((z : ℕ) : ℝ) = 0 := by
    rw [Finset.sum_pair hpq, hgp, hgq]
    simp only [Rat.smul_def]
    push_cast
    linarith [h]
  have hz := hli ({p, q} : Finset Nat.Primes) g hsum p (by simp)
  rw [hgp] at hz
  exact one_ne_zero hz

/-- ★★★★★★**[FrdI] Lemma 6.5, (ii)** —— 相異なる 6 素数について、
`log(p₀)/log(p₁) = λ₁·log(p₂)/log(p₃) = λ₂·log(p₄)/log(p₅)` は不可能。

★添字は `0`–`5`(原文の `1`–`6`)。 -/
theorem lemma_6_5_ii (P : Fin 6 → Nat.Primes) (hP : Function.Injective P)
    (l1 l2 : ℚ) (hl1 : 0 < l1) (hl2 : 0 < l2)
    (h1 : Real.log ((P 0 : ℕ) : ℝ) / Real.log ((P 1 : ℕ) : ℝ)
        = (l1 : ℝ) * (Real.log ((P 2 : ℕ) : ℝ) / Real.log ((P 3 : ℕ) : ℝ)))
    (h2 : Real.log ((P 0 : ℕ) : ℝ) / Real.log ((P 1 : ℕ) : ℝ)
        = (l2 : ℝ) * (Real.log ((P 4 : ℕ) : ℝ) / Real.log ((P 5 : ℕ) : ℝ))) :
    False := by
  classical
  set L : Fin 6 → ℝ := fun i => Real.log ((P i : ℕ) : ℝ) with hLdef
  have hLpos : ∀ i, 0 < L i := by
    intro i
    rw [hLdef]
    exact Real.log_pos (by exact_mod_cast (P i).2.one_lt)
  have hl1R : (0:ℝ) < (l1 : ℝ) := by exact_mod_cast hl1
  have hl2R : (0:ℝ) < (l2 : ℝ) := by exact_mod_cast hl2
  have h1' : L 0 / L 1 = (l1 : ℝ) * (L 2 / L 3) := h1
  have h2' : L 0 / L 1 = (l2 : ℝ) * (L 4 / L 5) := h2
  -- ★`t = log p₃ / log p₄`(原文の添字では `p₃/p₄`、ここでは `L 2 / L 3`)
  set t : ℝ := L 2 / L 3 with htdef
  -- ★`t` は有理数でない
  have hne23 : P 2 ≠ P 3 := fun h => absurd (hP h) (by decide)
  have htirr : ∀ r : ℚ, t ≠ (r : ℝ) := by
    intro r hr
    have hL3 : L 3 ≠ 0 := (hLpos 3).ne'
    have hkey : L 2 = (r : ℝ) * L 3 := by
      rw [htdef] at hr
      field_simp at hr
      linarith [hr]
    exact log_prime_ne_rat_mul hne23 r hkey
  -- ★2 つの一次独立性
  set xr : Fin 2 → ℝ := ![1, t] with hxrdef
  have hidx : Function.Injective (![1, 3, 5] : Fin 3 → Fin 6) := by decide
  set Q : Fin 3 → Nat.Primes := fun i => P (![1, 3, 5] i) with hQdef
  have hQinj : Function.Injective Q := hP.comp hidx
  set yr : Fin 3 → ℝ := fun i => Real.log ((Q i : ℕ) : ℝ) with hyrdef
  have hxli : LinearIndependent ℚ (fun j => ((xr j : ℝ) : ℂ)) :=
    linearIndependent_ofReal (linearIndependent_pair_one_real htirr)
  have hyli : LinearIndependent ℚ (fun i => ((yr i : ℝ) : ℂ)) :=
    linearIndependent_ofReal (log_primes_linearIndependent.comp Q hQinj)
  -- ★6 つの積の形
  have key : ∀ (μ : ℚ) (q : Nat.Primes) (j : Fin 2) (i : Fin 3),
      xr j * yr i = (μ : ℝ) * Real.log ((q : ℕ) : ℝ) →
      IsAlgebraic ℚ (Complex.exp (((xr j : ℝ) : ℂ) * ((yr i : ℝ) : ℂ))) := by
    intro μ q j i h
    have hmul : ((xr j : ℝ) : ℂ) * ((yr i : ℝ) : ℂ)
        = ((((μ : ℝ) * Real.log ((q : ℕ) : ℝ)) : ℝ) : ℂ) := by
      rw [← Complex.ofReal_mul, h]
    rw [hmul]
    exact isAlgebraic_exp_rat_mul_log q.2.pos μ
  have hL1 : L 1 ≠ 0 := (hLpos 1).ne'
  have hL3 : L 3 ≠ 0 := (hLpos 3).ne'
  have hL5 : L 5 ≠ 0 := (hLpos 5).ne'
  have e11 : xr 1 * yr 1 = ((1 : ℚ) : ℝ) * Real.log ((P 2 : ℕ) : ℝ) := by
    show t * L 3 = ((1 : ℚ) : ℝ) * L 2
    rw [htdef]
    push_cast
    field_simp
  have e10 : xr 1 * yr 0 = ((1 / l1 : ℚ) : ℝ) * Real.log ((P 0 : ℕ) : ℝ) := by
    show t * L 1 = ((1 / l1 : ℚ) : ℝ) * L 0
    have hl1' : (l1 : ℝ) ≠ 0 := ne_of_gt hl1R
    have hh : L 0 / L 1 = (l1 : ℝ) * (L 2 / L 3) := by rw [h1', htdef]
    rw [htdef]
    push_cast
    field_simp at hh ⊢
    linarith [hh]
  have e12 : xr 1 * yr 2 = ((l2 / l1 : ℚ) : ℝ) * Real.log ((P 4 : ℕ) : ℝ) := by
    show t * L 5 = ((l2 / l1 : ℚ) : ℝ) * L 4
    have hl1' : (l1 : ℝ) ≠ 0 := ne_of_gt hl1R
    have h12 : (l1 : ℝ) * (L 2 / L 3) = (l2 : ℝ) * (L 4 / L 5) := by
      rw [← htdef, ← h1', h2']
    rw [htdef]
    push_cast
    field_simp at h12 ⊢
    linarith [h12]
  have e00 : xr 0 * yr 0 = ((1 : ℚ) : ℝ) * Real.log ((P 1 : ℕ) : ℝ) := by
    show (1 : ℝ) * L 1 = ((1 : ℚ) : ℝ) * L 1
    norm_num
  have e01 : xr 0 * yr 1 = ((1 : ℚ) : ℝ) * Real.log ((P 3 : ℕ) : ℝ) := by
    show (1 : ℝ) * L 3 = ((1 : ℚ) : ℝ) * L 3
    norm_num
  have e02 : xr 0 * yr 2 = ((1 : ℚ) : ℝ) * Real.log ((P 5 : ℕ) : ℝ) := by
    show (1 : ℝ) * L 5 = ((1 : ℚ) : ℝ) * L 5
    norm_num
  -- ★six exponentials を当てる
  obtain ⟨j, i, htr⟩ := six_exponentials hxli hyli
  refine htr ?_
  fin_cases j <;> fin_cases i
  · exact key 1 (P 1) 0 0 e00
  · exact key 1 (P 3) 0 1 e01
  · exact key 1 (P 5) 0 2 e02
  · exact key (1 / l1) (P 0) 1 0 e10
  · exact key 1 (P 2) 1 1 e11
  · exact key (l2 / l1) (P 4) 1 2 e12

/-- ★locator —— `Lemma 6.5, (ii)`。 -/
def lemma_6_5_ii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5, (ii)",
    sectionId := "frdi-lemma-6-5" }

/-- ★★★★★★**[FrdI] Lemma 6.5** —— (i)(ii) が揃った。 -/
def lemma_6_5.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5",
    sectionId := "frdi-lemma-6-5" }

end ABC3.Found.FrdI
