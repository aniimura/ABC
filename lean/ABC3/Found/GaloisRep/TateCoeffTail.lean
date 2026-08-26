import ABC3.Found.GaloisRep.TateUWBridge

/-!
# Galois (G6) 第 233 ブロック —— **★★★★★`a₄`・`a₆` の切り詰めとの差**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★係数の側も揃える

第 231 で `X`・`Y` の差が揃った。`tateDefectTrunc` にはもう二つ
`partialEval tateA4` と `partialEval tateA6` が入っている。
解析側の対応物は第 220 の `−5·s₃` と `−(5s₃ + 7s₅)/12` である。

## ★★★★`σ_k` の一般の評価

第 230 の `σ₁ ≤ N²` を一般化する:

    σ_k(N) = ∑_{d|N} dᵏ ≤ |divisors| · Nᵏ ≤ N·Nᵏ = N^{k+1} ≤ (2^{k+1})ᴺ

★`N < 2ᴺ` から `N^{k+1} ≤ (2ᴺ)^{k+1} = (2^{k+1})ᴺ`。
★★`k = 3` なら底は `16`、`k = 5` なら `64`。よって `‖q‖ ≤ 1/128` を仮定に置けば
どちらも第 227 の一般の尾の評価が使える。

## ★★★`a₆` は 12 で割った係数

`tateA6` の係数は `−(5σ₃(N) + 7σ₅(N))/12`(整数——第 94 ブロックの `twelve_dvd`)。
★**割り算のまま扱わず、`12·tateA6 = −(5s₃ + 7s₅)`(`twelve_mul_tateA6`)を係数に降ろす**。
そうすれば `partialEval` の側も `12·partialEval tateA6 = −(5·PE₃ + 7·PE₅)` になり、
差は `−(5(s₃−PE₃) + 7(s₅−PE₅))/12` と書ける。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sigma_le_pow` | ★★`σ_k(N) ≤ N^{k+1}` |
| `pow_succ_le_two_pow_pow` | ★`N^{k+1} ≤ (2^{k+1})ᴺ` |
| `sigmaSum_eq_tsum_nat` | ★`sigmaSum` を `ℕ` 添字で |
| `norm_sigma_k_tail_le` | ★★★★★`σ_k` 級数の尾の評価 |
| `norm_sigmaSum_sub_partialEval_le` | ★★★★★★`sigmaSum k` と切り詰めの差 |
| `coeff_twelve_tateA6` / `partialEval_tateA6` | ★★★`12·a₆ = −(5s₃+7s₅)` を係数に降ろす |
| `norm_a4_sub_partialEval_le` | ★★★★★**`a₄` の差** |
| `norm_a6_sub_partialEval_le` | ★★★★★**`a₆` の差** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★`σ_k` の評価 -/

/-- ★★**`σ_k(N) ≤ N^{k+1}`**。 -/
theorem sigma_le_pow (k N : ℕ) : ArithmeticFunction.sigma k N ≤ N ^ (k + 1) := by
  rw [ArithmeticFunction.sigma_apply]
  have hsub : N.divisors ⊆ Finset.Ico 1 (N + 1) := by
    intro d hd
    rw [Nat.mem_divisors] at hd
    rw [Finset.mem_Ico]
    exact ⟨Nat.one_le_iff_ne_zero.2 (fun h => hd.2 (by simpa [h] using hd.1)),
      Nat.lt_succ_of_le (Nat.le_of_dvd (Nat.pos_of_ne_zero hd.2) hd.1)⟩
  have hcard : N.divisors.card ≤ N := by
    have hc := Finset.card_le_card hsub
    simpa using hc
  have hle : ∀ d ∈ N.divisors, d ^ k ≤ N ^ k := by
    intro d hd
    exact Nat.pow_le_pow_left (Nat.divisor_le hd) k
  calc ∑ d ∈ N.divisors, d ^ k ≤ ∑ _d ∈ N.divisors, N ^ k := Finset.sum_le_sum hle
    _ = N.divisors.card * N ^ k := by rw [Finset.sum_const, smul_eq_mul]
    _ ≤ N * N ^ k := Nat.mul_le_mul_right _ hcard
    _ = N ^ (k + 1) := by ring

/-- ★`N^{k+1} ≤ (2^{k+1})ᴺ`。 -/
theorem pow_succ_le_two_pow_pow (k N : ℕ) : N ^ (k + 1) ≤ (2 ^ (k + 1)) ^ N := by
  have h : N ≤ 2 ^ N := Nat.lt_two_pow_self.le
  calc N ^ (k + 1) ≤ (2 ^ N) ^ (k + 1) := Nat.pow_le_pow_left h (k + 1)
    _ = (2 ^ (k + 1)) ^ N := by rw [← pow_mul, mul_comm, pow_mul]

/-! ## ★`sigmaSum` を `ℕ` 添字で -/

theorem sigmaSum_eq_tsum_nat (k : ℕ) (τ : UpperHalfPlane) :
    sigmaSum k τ
      = ∑' N : ℕ, ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ)
        * Complex.exp (2 * ↑π * I * τ) ^ (N + 1) := by
  rw [sigmaSum, ← (Equiv.pnatEquivNat.symm).tsum_eq
    (fun n : ℕ+ => ((ArithmeticFunction.sigma k (n : ℕ) : ℕ) : ℂ)
      * Complex.exp (2 * ↑π * I * τ) ^ ((n : ℕ) : ℤ))]
  refine tsum_congr fun N => ?_
  have hcoe : ((Equiv.pnatEquivNat.symm N : ℕ+) : ℕ) = N + 1 := by
    simp [Equiv.pnatEquivNat]
  rw [hcoe, zpow_natCast]

/-! ## ★★★★★`σ_k` 級数の尾 -/

theorem norm_sigma_k_term_le (k : ℕ) (q : ℂ) (N : ℕ) :
    ‖((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1)‖
      ≤ ((2 ^ (k + 1) : ℝ) * ‖q‖) ^ (N + 1) := by
  rw [norm_mul, norm_pow, Complex.norm_natCast]
  have h1 : (ArithmeticFunction.sigma k (N + 1) : ℕ) ≤ (2 ^ (k + 1)) ^ (N + 1) :=
    le_trans (sigma_le_pow k (N + 1)) (pow_succ_le_two_pow_pow k (N + 1))
  have h1' : ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℝ) ≤ ((2 : ℝ) ^ (k + 1)) ^ (N + 1) := by
    exact_mod_cast h1
  calc ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℝ) * ‖q‖ ^ (N + 1)
      ≤ ((2 : ℝ) ^ (k + 1)) ^ (N + 1) * ‖q‖ ^ (N + 1) :=
        mul_le_mul_of_nonneg_right h1' (pow_nonneg (norm_nonneg q) _)
    _ = ((2 ^ (k + 1) : ℝ) * ‖q‖) ^ (N + 1) := by rw [mul_pow]

theorem summable_sigma_k_series (k : ℕ) (q : ℂ) (hq : (2 ^ (k + 1) : ℝ) * ‖q‖ ≤ 1 / 2) :
    Summable fun N : ℕ => ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1) := by
  have h40 : (0 : ℝ) ≤ (2 ^ (k + 1) : ℝ) * ‖q‖ := by positivity
  have h4 : (2 ^ (k + 1) : ℝ) * ‖q‖ < 1 := by linarith
  have hg : Summable fun N : ℕ => ((2 ^ (k + 1) : ℝ) * ‖q‖) ^ (N + 1) :=
    (summable_nat_add_iff 1).2 (summable_geometric_of_lt_one h40 h4)
  exact Summable.of_norm_bounded hg (fun N => norm_sigma_k_term_le k q N)

/-- ★★★★★**`σ_k` 級数の尾の評価**。 -/
theorem norm_sigma_k_tail_le (k : ℕ) (q : ℂ) (hq : (2 ^ (k + 1) : ℝ) * ‖q‖ ≤ 1 / 2) (n : ℕ) :
    ‖(∑' N : ℕ, ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1))
        - partialSum (fun N => ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1)) n‖
      ≤ 2 * ((2 ^ (k + 1) : ℝ) * ‖q‖) ^ (n + 1) := by
  have h40 : (0 : ℝ) ≤ (2 ^ (k + 1) : ℝ) * ‖q‖ := by positivity
  have h := norm_tsum_sub_partialSum_le
    (fun N => ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1))
    ((2 ^ (k + 1) : ℝ) * ‖q‖) 1 h40 hq (by norm_num) (summable_sigma_k_series k q hq)
    (fun N => by simpa using norm_sigma_k_term_le k q N) n
  simpa using h

theorem partialEval_sigmaSeries_k_succ (k : ℕ) (q : ℂ) (n : ℕ) :
    partialEval (sigmaSeries k) q (n + 1)
      = partialSum (fun N => ((ArithmeticFunction.sigma k (N + 1) : ℕ) : ℂ) * q ^ (N + 1)) n := by
  rw [partialEval, partialSum, Finset.sum_range_succ']
  have h0 : (((PowerSeries.coeff 0 (sigmaSeries k) : ℤ)) : ℂ) * q ^ 0 = 0 := by
    rw [coeff_sigmaSeries]
    simp
  rw [h0, add_zero]
  refine Finset.sum_congr rfl fun N _ => ?_
  rw [coeff_sigmaSeries, if_neg (Nat.succ_ne_zero N)]
  push_cast
  ring

/-- ★★★★★★**`sigmaSum k` と `partialEval (sigmaSeries k)` の差**。 -/
theorem norm_sigmaSum_sub_partialEval_le (k : ℕ) (τ : UpperHalfPlane)
    (hq : (2 ^ (k + 1) : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2) (n : ℕ) :
    ‖sigmaSum k τ - partialEval (sigmaSeries k) (Complex.exp (2 * ↑π * I * τ)) (n + 1)‖
      ≤ 2 * ((2 ^ (k + 1) : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1) := by
  rw [sigmaSum_eq_tsum_nat k τ, partialEval_sigmaSeries_k_succ]
  exact norm_sigma_k_tail_le k _ hq n

/-! ## ★★★`a₄`・`a₆` の係数 -/

theorem coeff_twelve_tateA6 (N : ℕ) :
    (12 : ℤ) * PowerSeries.coeff N tateA6
      = -(5 * PowerSeries.coeff N (sigmaSeries 3) + 7 * PowerSeries.coeff N (sigmaSeries 5)) := by
  have h := congrArg (PowerSeries.coeff N) twelve_mul_tateA6
  rw [PowerSeries.coeff_C_mul, map_neg, map_add, PowerSeries.coeff_C_mul,
    PowerSeries.coeff_C_mul] at h
  exact h

theorem partialEval_tateA4 (q : ℂ) (m : ℕ) :
    partialEval tateA4 q m = -5 * partialEval (sigmaSeries 3) q m := by
  rw [partialEval, partialEval, Finset.mul_sum]
  refine Finset.sum_congr rfl fun N _ => ?_
  rw [tateA4, PowerSeries.coeff_C_mul]
  push_cast
  ring

/-- ★★★**`12·partialEval a₆ = −(5·PE₃ + 7·PE₅)`**——割り算のまま扱わない。 -/
theorem partialEval_tateA6 (q : ℂ) (m : ℕ) :
    (12 : ℂ) * partialEval tateA6 q m
      = -(5 * partialEval (sigmaSeries 3) q m + 7 * partialEval (sigmaSeries 5) q m) := by
  simp only [partialEval, Finset.mul_sum, ← Finset.sum_add_distrib, ← Finset.sum_neg_distrib]
  refine Finset.sum_congr rfl fun N _ => ?_
  have hc := coeff_twelve_tateA6 N
  have hc' : (12 : ℂ) * ((PowerSeries.coeff N tateA6 : ℤ) : ℂ)
      = -(5 * ((PowerSeries.coeff N (sigmaSeries 3) : ℤ) : ℂ)
        + 7 * ((PowerSeries.coeff N (sigmaSeries 5) : ℤ) : ℂ)) := by
    exact_mod_cast congrArg (fun x : ℤ => (x : ℂ)) hc
  linear_combination (q ^ N) * hc'

/-! ## ★★★★★`a₄`・`a₆` の差 -/

/-- ★★★★★**`a₄` の解析値と切り詰めの差**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_a4_sub_partialEval_le (τ : UpperHalfPlane)
    (hq : (16 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2) (n : ℕ) :
    ‖(-5 * sigmaSum 3 τ)
        - partialEval tateA4 (Complex.exp (2 * ↑π * I * τ)) (n + 1)‖
      ≤ 10 * ((16 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1) := by
  have hq' : ((2 : ℝ) ^ (3 + 1)) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2 := by
    norm_num
    linarith
  have h := norm_sigmaSum_sub_partialEval_le 3 τ hq' n
  rw [partialEval_tateA4]
  have hkey : (-5 * sigmaSum 3 τ)
      - (-5) * partialEval (sigmaSeries 3) (Complex.exp (2 * ↑π * I * τ)) (n + 1)
      = (-5) * (sigmaSum 3 τ
        - partialEval (sigmaSeries 3) (Complex.exp (2 * ↑π * I * τ)) (n + 1)) := by ring
  rw [hkey, norm_mul]
  have h5 : ‖(-5 : ℂ)‖ = 5 := by simp
  rw [h5]
  have hnum : ((2 : ℝ) ^ (3 + 1)) = 16 := by norm_num
  rw [hnum] at h
  linarith

/-- ★★★★★**`a₆` の解析値と切り詰めの差**。 -/
theorem norm_a6_sub_partialEval_le (τ : UpperHalfPlane)
    (hq : (64 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2) (n : ℕ) :
    ‖(-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12)
        - partialEval tateA6 (Complex.exp (2 * ↑π * I * τ)) (n + 1)‖
      ≤ 2 * ((64 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1) := by
  have hr0 : (0 : ℝ) ≤ ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ := norm_nonneg _
  have hq3' : ((2 : ℝ) ^ (3 + 1)) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2 := by
    norm_num
    linarith
  have hq5' : ((2 : ℝ) ^ (5 + 1)) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ ≤ 1 / 2 := by
    norm_num
    linarith
  have h3 := norm_sigmaSum_sub_partialEval_le 3 τ hq3' n
  have h5 := norm_sigmaSum_sub_partialEval_le 5 τ hq5' n
  have hn3 : ((2 : ℝ) ^ (3 + 1)) = 16 := by norm_num
  have hn5 : ((2 : ℝ) ^ (5 + 1)) = 64 := by norm_num
  rw [hn3] at h3
  rw [hn5] at h5
  have hpe := partialEval_tateA6 (Complex.exp (2 * ↑π * I * (τ : ℂ))) (n + 1)
  have hkey : (-(5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12)
      - partialEval tateA6 (Complex.exp (2 * ↑π * I * τ)) (n + 1)
      = (-5 / 12) * (sigmaSum 3 τ
          - partialEval (sigmaSeries 3) (Complex.exp (2 * ↑π * I * τ)) (n + 1))
        + (-7 / 12) * (sigmaSum 5 τ
          - partialEval (sigmaSeries 5) (Complex.exp (2 * ↑π * I * τ)) (n + 1)) := by
    linear_combination (-1 / 12 : ℂ) * hpe
  rw [hkey]
  have hcomb := norm_add_le
    ((-5 / 12 : ℂ) * (sigmaSum 3 τ
      - partialEval (sigmaSeries 3) (Complex.exp (2 * ↑π * I * τ)) (n + 1)))
    ((-7 / 12 : ℂ) * (sigmaSum 5 τ
      - partialEval (sigmaSeries 5) (Complex.exp (2 * ↑π * I * τ)) (n + 1)))
  rw [norm_mul, norm_mul] at hcomb
  have hc5 : ‖(-5 / 12 : ℂ)‖ = 5 / 12 := by
    rw [norm_div]
    simp
  have hc7 : ‖(-7 / 12 : ℂ)‖ = 7 / 12 := by
    rw [norm_div]
    simp
  rw [hc5, hc7] at hcomb
  have hmono : ((16 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1)
      ≤ ((64 : ℝ) * ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖) ^ (n + 1) :=
    pow_le_pow_left₀ (by positivity) (by linarith) _
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def sigma_le_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——sigma_k の一般の評価)",
    sectionId := "genell-def-3-3" }

def norm_a4_sub_partialEval_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——a4 の解析値と切り詰めの差)",
    sectionId := "genell-def-3-3" }

def norm_a6_sub_partialEval_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——a6 の解析値と切り詰めの差)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
