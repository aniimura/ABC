/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Tactic.NormNum.Prime
import ABC3.Found.GenEll.PrimeNumberTheorem
import ABC3.Found.GenEll.PrimeConstants.Remark411

/-!
# PrimeConstants —— `[GenEll] Corollary 4.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real Finset

/-! ## ★★★対数和の大きい有限素数集合 -/

/-- ★★★**いくらでも大きい対数和を持つ有限素数集合がある**。

原文 (GenEll p.23):
> Corollary 4.3, except that instead of applying condition (a) of Theorem 3.8, we

★原文は『by enlarging `S` [and possibly increasing the "C" of condition (c)],
we may always assume that `x_ϵ ≤ x_S`』と 1 文で済ませている。
★★その中身は `θ(x) → ∞` であり、`θ(x) > (5/6)·x` から直ちに出る。 -/
theorem exists_finset_primes_sum_log_gt (y : ℝ) :
    ∃ A : Finset ℕ, (∀ p ∈ A, p.Prime) ∧ y < ∑ p ∈ A, Real.log p := by
  obtain ⟨xeps, Ceps, hxpos, _, _, _, hlow, _⟩ := exists_cond_i_ii 1
  set x : ℝ := max xeps (max (2 * y + 1) 1) with hxdef
  have hxpos' : (0:ℝ) < x :=
    lt_of_lt_of_le zero_lt_one (le_trans (le_max_right _ _) (le_max_right _ _))
  have hx : xeps ≤ x := le_max_left _ _
  have h := hlow x hxpos' hx
  have hy : 2 * y + 1 ≤ x := le_trans (le_max_left _ _) (le_max_right _ _)
  refine ⟨(Finset.Ioc 0 ⌊x⌋₊).filter Nat.Prime, fun p hp => (Finset.mem_filter.1 hp).2, ?_⟩
  have hsum : ∑ p ∈ (Finset.Ioc 0 ⌊x⌋₊).filter Nat.Prime, Real.log p
      = Chebyshev.theta x := by
    rw [theta_eq_sum_logPrime, Finset.sum_filter]
    exact Finset.sum_congr rfl (fun p _ => by
      by_cases hp : p.Prime
      · simp [hp, logPrime_of_prime hp]
      · simp [hp, logPrime])
  rw [hsum]
  linarith

/-! ## ★★★★★`Corollary 4.3` / `4.4` の (c) を出す数値の段 -/

/-- ★★★★★**原文 p.22-23 の不等式の連鎖を 1 本にしたもの**。

原文 (GenEll p.23):
> Corollary 4.3, except that instead of applying condition (a) of Theorem 3.8, we

引数の意味: `d = [L:ℚ]`、`F = ht^Falt`、`dinf = deg∞`、`x_S`、`x_bad = x_{S∘} − x_S`、
`x_T`(『enlarging S』の分)、`extra`(`l•` のとき `3d·log-diff`)、`H`(`Lemma 4.1` の `h`)、
`l`、`A`(`Proposition 3.4` の定数)、`B`(`ht^Falt` の下界の絶対値)、`P = d^{1+ϵ}`。

★係数の帳尻は原文 p.23 の『`2·3·12 + 8·100 ≤ 100 + 800 = 900`』である:

    2·x_bad ≤ 5·23040·d·deg∞ ≤ 72·23040·d·F + …   (Prop 3.4 を 1+ϵ = 6/5 で)
    8·H     ≤ 800·23040·d·F + …
    72 + 800 = 872 ≤ 900

★★`F < 0` のとき `872 ≤ 900` は逆向きになるが、差 `28·23040·d·(F+B) ≥ 0` がちょうど埋める
——これが定数の `828 = 800 + 28` の出所である。 -/
theorem cor4_numeric (d F dinf xS xbad xT extra H l A B C₈ P : ℝ)
    (hd1 : 1 ≤ d) (hP1 : 1 ≤ P) (hdP : d ≤ P)
    (hA : dinf ≤ 12 * (1 + 1/5) * F + A)
    (hxbad : xbad ≤ 5/2 * (23040 * (d * dinf)))
    (hB : -B ≤ F) (hBnn : 0 ≤ B) (hC₈ : 0 < C₈) (hxT : 0 ≤ xT)
    (hH : H ≤ 23040 * 100 * (d * F) + 23040 * 100 * (C₈ * P) + 23040 * 100 * (d * B))
    (hl : l ≤ 2 * (xS + xbad + xT + extra) + 8 * H) :
    l ≤ 23040 * 900 * (d * F) + 2 * extra + 2 * xS
        + (23040 * (5 * |A| + 800 * C₈ + 828 * B) + 2 * xT + 1) * P := by
  have hdpos : (0:ℝ) < d := by linarith
  have hAabs : A ≤ |A| := le_abs_self A
  have hdA : d * dinf ≤ 12 * (1 + 1/5) * (d * F) + d * A := by
    have hm := mul_le_mul_of_nonneg_left hA hdpos.le
    linarith [hm]
  have h1 : d * A ≤ d * |A| := mul_le_mul_of_nonneg_left hAabs hdpos.le
  have h2 : d * |A| ≤ P * |A| := mul_le_mul_of_nonneg_right hdP (abs_nonneg A)
  have hdBP : d * B ≤ P * B := mul_le_mul_of_nonneg_right hdP hBnn
  have hxTP : xT ≤ xT * P := by nlinarith
  have hFB : (0:ℝ) ≤ d * F + d * B := by nlinarith
  have hPnn : (0:ℝ) ≤ P := by linarith
  linarith [hl, hxbad, hdA, h1, h2, hdBP, hxTP, hFB, hH, hPnn]

def exists_finset_primes_sum_log_gt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.3 の証明(by enlarging S … we may always assume that x_ϵ ≤ x_S)",
    sectionId := "genell-cor-4-3" }

def cor4_numeric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 23,
    item := "Corollary 4.3 の証明(2·3·12 + 8·100 ≤ 100 + 800 = 900)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
