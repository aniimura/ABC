/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Tactic.NormNum.Prime
import ABC3.Found.GenEll.PrimeNumberTheorem

/-!
# PrimeConstants —— `[GenEll] Remark 4.1.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-! ## ★`θ` の単調性 -/

/-- ★**`θ` は非減少**である——被和関数が非負だから。 -/
theorem theta_mono {x y : ℝ} (h : x ≤ y) : Chebyshev.theta x ≤ Chebyshev.theta y := by
  have hf : ⌊x⌋₊ ≤ ⌊y⌋₊ := Nat.floor_le_floor h
  have hsub := theta_sub_theta x y hf
  have hnn : (0:ℝ) ≤ ∑ p ∈ Finset.Ioc ⌊x⌋₊ ⌊y⌋₊, logPrime p :=
    Finset.sum_nonneg (fun p _ => logPrime_nonneg p)
  linarith

theorem theta_nonneg (x : ℝ) : 0 ≤ Chebyshev.theta x := by
  rw [theta_eq_sum_logPrime]
  exact Finset.sum_nonneg (fun p _ => logPrime_nonneg p)

/-! ## ★★★★条件 (i)(ii) を満たす定数の存在 -/

/-- ★★★★**`ϵ = 1/6` に対して `Lemma 4.1` の条件 (i)(ii) を満たす `x_ϵ, C_ϵ` が存在する**。

原文 (GenEll p.21):
> which satisfy condition (ii) of Lemma 4.1 is entirely elementary. On the other hand, with regard to condition (i), the fact that

★(i) の 2 式はここでは**素数でない `x`** という制限なしに取れるので、
`lemma_4_1` の仮説にそのまま渡せる。 -/
theorem exists_cond_i_ii (M : ℕ) :
    ∃ xeps Ceps : ℝ, 0 < xeps ∧ 0 < Ceps ∧ Ceps < (1/6 : ℝ) * xeps
      ∧ (∀ x : ℝ, 0 < x → Chebyshev.theta x < 5 / 4 * x + Ceps)
      ∧ (∀ x : ℝ, 0 < x → xeps ≤ x → (1 - (1/6 : ℝ)) * x < Chebyshev.theta x)
      ∧ (∀ x : ℝ, 0 < x → xeps ≤ x → (M : ℝ) * Real.log x ≤ (1/6 : ℝ) * x) := by
  have htend := theta_div_tendsto_one
  have hev : ∀ᶠ x : ℝ in Filter.atTop, |Chebyshev.theta x / x - 1| < 1/6 := by
    have := Metric.tendsto_nhds.1 htend (1/6) (by norm_num)
    simpa [Real.dist_eq] using this
  obtain ⟨X₁, hX₁⟩ := Filter.eventually_atTop.mp hev
  set X₀ : ℝ := max X₁ 1 with hX₀def
  have hX₀pos : (0:ℝ) < X₀ := lt_of_lt_of_le zero_lt_one (le_max_right _ _)
  have hkey : ∀ x : ℝ, X₀ ≤ x → 0 < x →
      (5/6 : ℝ) * x < Chebyshev.theta x ∧ Chebyshev.theta x < 5/4 * x := by
    intro x hx hxpos
    have h := hX₁ x (le_trans (le_max_left _ _) hx)
    have h2 := abs_lt.1 h
    constructor
    · have hlt : (1 - 1/6 : ℝ) < Chebyshev.theta x / x := by linarith [h2.1]
      rw [lt_div_iff₀ hxpos] at hlt
      linarith
    · have hlt : Chebyshev.theta x / x < 1 + 1/6 := by linarith [h2.2]
      rw [div_lt_iff₀ hxpos] at hlt
      linarith
  set Ceps : ℝ := Chebyshev.theta X₀ + 1 with hCdef
  have hCpos : 0 < Ceps := by
    have := theta_nonneg X₀
    simp only [hCdef]
    linarith
  obtain ⟨y₀, hy₀pos, hy₀⟩ := exists_xeps_cond_ii M (1/6 : ℝ) (by norm_num)
  refine ⟨max (max X₀ y₀) (6 * (Ceps + 1)), Ceps, ?_, hCpos, ?_, ?_, ?_, ?_⟩
  · exact lt_of_lt_of_le hX₀pos (le_trans (le_max_left _ _) (le_max_left _ _))
  · have h6 : 6 * (Ceps + 1) ≤ max (max X₀ y₀) (6 * (Ceps + 1)) := le_max_right _ _
    linarith
  · intro x hxpos
    rcases le_or_gt X₀ x with hx | hx
    · have hup := (hkey x hx hxpos).2
      linarith
    · have hle : Chebyshev.theta x ≤ Chebyshev.theta X₀ := theta_mono hx.le
      have hq : (0:ℝ) < 5 / 4 * x := by linarith
      simp only [hCdef]
      linarith
  · intro x hxpos hx
    have hx0 : X₀ ≤ x := le_trans (le_trans (le_max_left _ _) (le_max_left _ _)) hx
    have hlow := (hkey x hx0 hxpos).1
    linarith
  · intro x _ hx
    exact hy₀ x (le_trans (le_trans (le_max_right _ _) (le_max_left _ _)) hx)

/-- ★素数 `l` が `2, 3, 5` のいずれでもなければ `30` と互いに素である。

★`Corollary 4.4` が `Theorem 3.8` の条件 (b) を使うときに要る(`30 = 2·3·5`)。 -/
theorem coprime_thirty_of_prime {l : ℕ} (hl : l.Prime) (h2 : l ≠ 2) (h3 : l ≠ 3) (h5 : l ≠ 5) :
    Nat.Coprime l 30 := by
  rw [Nat.Prime.coprime_iff_not_dvd hl]
  intro hd
  have h : l ∣ 2 * (3 * 5) := by
    have h30 : (2 * (3 * 5) : ℕ) = 30 := by norm_num
    rw [h30]; exact hd
  rcases (Nat.Prime.dvd_mul hl).1 h with hh | hh
  · exact h2 ((Nat.prime_dvd_prime_iff_eq hl Nat.prime_two).1 hh)
  · rcases (Nat.Prime.dvd_mul hl).1 hh with hh' | hh'
    · exact h3 ((Nat.prime_dvd_prime_iff_eq hl Nat.prime_three).1 hh')
    · exact h5 ((Nat.prime_dvd_prime_iff_eq hl (by norm_num)).1 hh')

/-! ## ★出典の紐付け(`.src`) -/

def exists_cond_i_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Remark 4.1.1(条件 (i)(ii) を満たす定数の存在)",
    sectionId := "genell-rem-4-1-1" }

end ABC3.Found.GenEll
