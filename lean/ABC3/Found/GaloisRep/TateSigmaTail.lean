import ABC3.Found.GaloisRep.TateSigmaComplex

/-!
# Galois (G6) 第 230 ブロック —— **★★★★★★`s₁` の切り詰めとの差**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★最後の部品

第 227 で `X`・`Y` の尾の評価、第 229 で `s₁` の q 展開(ℂ の上)が取れた。
残るのは **`s₁` の解析値と `partialEval (sigmaSeries 1) q n` の差**である。

    ‖∑_{m≥1}f(qᵐ) − partialEval (sigmaSeries 1) q (n+1)‖ ≤ 2·(4‖q‖)^{n+1}

## ★★★★係数は `4ᴺ` で押さえる

`σ₁(N) = ∑_{d|N} d ≤ N·|divisors(N)| ≤ N·N = N²`(約数はすべて `N` 以下、個数も `N` 以下)。
さらに `N < 2ᴺ` なので `N² ≤ 4ᴺ`。よって

    ‖σ₁(N)qᴺ‖ ≤ (4‖q‖)ᴺ

★★**`‖q‖ ≤ 1/8` を仮定に置く**と `4‖q‖ ≤ 1/2` になり、第 227 の一般の尾の評価が
`r := 4‖q‖` でそのまま使える。

★★★**定数が `n` に依存してよい**ことに注意——`(4‖q‖)^{n+1} = 4^{n+1}‖q‖^{n+1}` は
`‖q‖^{n+1}` の定数倍(定数は `4^{n+1}`)である。各 `n` を固定して `q → 0` を見るので、
`n` に依存する定数で十分である。

## ★`partialEval` と部分和の対応

`sigmaSeries 1` の `0` 次係数は `0` なので

    partialEval (sigmaSeries 1) q (n+1) = ∑_{k<n} σ₁(k+1)q^{k+1}

——`Finset.sum_range_succ'` で先頭を切り出せばよい(`partialEval_sigmaSeries_succ`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sigma_one_le_sq` | ★★`σ₁(N) ≤ N²` |
| `sq_le_four_pow` | ★`N² ≤ 4ᴺ` |
| `norm_sigma_term_le` | ★★`‖σ₁(N)qᴺ‖ ≤ (4‖q‖)ᴺ` |
| `norm_sigma_tail_le` | ★★★★★`σ₁` 級数の尾の評価 |
| `partialEval_sigmaSeries_succ` | ★`partialEval` は部分和である |
| `norm_tateXtail_one_sub_partialEval_le` | ★★★★★★**`s₁` の解析値と切り詰めの差** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★係数の評価 -/

/-- ★★**`σ₁(N) ≤ N²`**——約数はすべて `N` 以下、個数も `N` 以下。 -/
theorem sigma_one_le_sq (N : ℕ) : ArithmeticFunction.sigma 1 N ≤ N ^ 2 := by
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
  have hle : ∀ d ∈ N.divisors, d ^ 1 ≤ N := by
    intro d hd
    rw [pow_one]
    exact Nat.divisor_le hd
  calc ∑ d ∈ N.divisors, d ^ 1 ≤ ∑ _d ∈ N.divisors, N := Finset.sum_le_sum hle
    _ = N.divisors.card * N := by rw [Finset.sum_const, smul_eq_mul]
    _ ≤ N * N := Nat.mul_le_mul_right N hcard
    _ = N ^ 2 := by ring

/-- ★`N² ≤ 4ᴺ`(`N < 2ᴺ` から)。 -/
theorem sq_le_four_pow (N : ℕ) : N ^ 2 ≤ 4 ^ N := by
  have h : N ≤ 2 ^ N := Nat.lt_two_pow_self.le
  calc N ^ 2 ≤ (2 ^ N) ^ 2 := Nat.pow_le_pow_left h 2
    _ = 4 ^ N := by rw [← pow_mul, mul_comm, pow_mul]; norm_num

theorem norm_sigma_term_le (q : ℂ) (k : ℕ) :
    ‖((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1)‖
      ≤ (4 * ‖q‖) ^ (k + 1) := by
  rw [norm_mul, norm_pow, Complex.norm_natCast]
  have h1 : (ArithmeticFunction.sigma 1 (k + 1) : ℕ) ≤ 4 ^ (k + 1) :=
    le_trans (sigma_one_le_sq (k + 1)) (sq_le_four_pow (k + 1))
  have h1' : ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℝ) ≤ (4 : ℝ) ^ (k + 1) := by
    exact_mod_cast h1
  calc ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℝ) * ‖q‖ ^ (k + 1)
      ≤ (4 : ℝ) ^ (k + 1) * ‖q‖ ^ (k + 1) :=
        mul_le_mul_of_nonneg_right h1' (pow_nonneg (norm_nonneg q) _)
    _ = (4 * ‖q‖) ^ (k + 1) := by rw [mul_pow]

/-! ## ★★★★★尾の評価 -/

theorem summable_sigma_series (q : ℂ) (hq : ‖q‖ ≤ 1 / 8) :
    Summable fun k : ℕ => ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1) := by
  have h4 : 4 * ‖q‖ < 1 := by linarith
  have h40 : (0 : ℝ) ≤ 4 * ‖q‖ := by positivity
  have hg : Summable fun k : ℕ => (4 * ‖q‖) ^ (k + 1) :=
    (summable_nat_add_iff 1).2 (summable_geometric_of_lt_one h40 h4)
  exact Summable.of_norm_bounded hg (fun k => norm_sigma_term_le q k)

/-- ★★★★★**`σ₁` 級数の尾の評価**——定数は `n` に依存してよい。 -/
theorem norm_sigma_tail_le (q : ℂ) (hq : ‖q‖ ≤ 1 / 8) (n : ℕ) :
    ‖(∑' k : ℕ, ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1))
        - partialSum (fun k => ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1)) n‖
      ≤ 2 * (4 * ‖q‖) ^ (n + 1) := by
  have h40 : (0 : ℝ) ≤ 4 * ‖q‖ := by positivity
  have h4 : 4 * ‖q‖ ≤ 1 / 2 := by linarith
  have h := norm_tsum_sub_partialSum_le
    (fun k => ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1))
    (4 * ‖q‖) 1 h40 h4 (by norm_num) (summable_sigma_series q hq)
    (fun k => by simpa using norm_sigma_term_le q k) n
  simpa using h

/-! ## ★`partialEval` と部分和 -/

theorem partialEval_sigmaSeries_succ (q : ℂ) (n : ℕ) :
    partialEval (sigmaSeries 1) q (n + 1)
      = partialSum (fun k => ((ArithmeticFunction.sigma 1 (k + 1) : ℕ) : ℂ) * q ^ (k + 1)) n := by
  rw [partialEval, partialSum, Finset.sum_range_succ']
  have h0 : (((PowerSeries.coeff 0 (sigmaSeries 1) : ℤ)) : ℂ) * q ^ 0 = 0 := by
    rw [coeff_sigmaSeries]
    simp
  rw [h0, add_zero]
  refine Finset.sum_congr rfl fun k _ => ?_
  rw [coeff_sigmaSeries, if_neg (Nat.succ_ne_zero k)]
  push_cast
  ring

/-- ★★★★★★**`s₁` の解析値と切り詰めの差**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem norm_tateXtail_one_sub_partialEval_le (q : ℂ) (hq : ‖q‖ ≤ 1 / 8) (n : ℕ) :
    ‖(∑' m : ℕ, tateXterm (q ^ (m + 1))) - partialEval (sigmaSeries 1) q (n + 1)‖
      ≤ 2 * (4 * ‖q‖) ^ (n + 1) := by
  have hq1 : ‖q‖ < 1 := by linarith
  rw [tsum_nat_tateXterm_eq_sigma q hq1, partialEval_sigmaSeries_succ]
  exact norm_sigma_tail_le q hq n

/-! ## ★出典の紐付け(`.src`) -/

def sigma_one_le_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——sigma_1 の評価)",
    sectionId := "genell-def-3-3" }

def norm_tateXtail_one_sub_partialEval_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——s1 の解析値と切り詰めの差)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
