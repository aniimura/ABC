/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.IdealPowBound
import Mathlib.Analysis.SpecialFunctions.Log.Summable
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.Analysis.PSeries

/-!
# `log ζ_L(s) = Σ_p a_K(p)·p^{-s} + O(1)`(鎖 `cheb` の `cheb-log-zeta`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何をするか

`Theorem 6.4, (iv)` が使う Tchebotarev 密度定理を、**類体論を使わずに**
「完全分解する素点の密度は `1/[L:K]`」の形で出すための二段目である。

★★★**要点**: 原文の `Σ_{f(𝔭)=1} N𝔭^{-s}` は、**有理素数ごとの計数**で書き直せる ——
`f(𝔭) = 1` の素イデアルはちょうど「ノルムが有理素数 `p` のイデアル」であり、
その個数が `a_K(p)` だからである。したがって

  `log ζ_L(s) = Σ_p a_K(p)·p^{-s} + O(1)`

を示せばよく、**素イデアルへの局所分解は要らない**。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `zetaR` / `zetaR_eulerProduct_hasProd` | ★**実の** Euler 積(`ℝ` は完備ノルム環) |
| `tsum_log_localFactor` | ★★★`Σ_p log L_p(s) = log ζ_L(s)` |
| `zetaSummandR_pow_le` | `e ≥ 2` の項の一様評価(`a(p^e) ≤ (ed+1)^d` を使う) |
| `abs_log_one_add_sub_le` | `\|log(1+A+R) − A\| ≤ R + (A+R)²` |
| `abs_log_localFactor_sub_le` | 素数ごとの一様評価 `≤ C/p²` |
| `abs_log_zetaR_sub_tsum_le` | ★★★★★★**`\|log ζ_L(s) − Σ_p a(p)p^{-s}\| ≤ C`** |

★★`log ζ_L(s) = Σ_p log L_p` は「正項の無限積の対数」であり、
`Σ_{p ∈ t} log L_p = log ∏_{p∈t} L_p ≤ log ζ_L(s)` で**部分和が一様有界**だから
和が収束し、`exp` を取って一意性で一致する。
-/

namespace ABC3.Found.NF

open NumberField Ideal Filter Topology Complex

/-! ## ★0. 一般の補題 -/

open Filter Topology in
/-- ★`1` 以上の項の無限積では、有限部分積は極限以下。 -/
theorem prod_le_of_hasProd_of_one_le {ι : Type*} [DecidableEq ι] {f : ι → ℝ} {a : ℝ}
    (h : HasProd f a) (h1 : ∀ i, 1 ≤ f i) (s : Finset ι) : ∏ i ∈ s, f i ≤ a := by
  rw [HasProd, SummationFilter.unconditional] at h
  refine ge_of_tendsto h ?_
  filter_upwards [Filter.eventually_ge_atTop s] with t ht
  have hs1 : (1 : ℝ) ≤ ∏ i ∈ s, f i := Finset.one_le_prod (fun i _ => h1 i)
  have hd1 : (1 : ℝ) ≤ ∏ i ∈ t \ s, f i := Finset.one_le_prod (fun i _ => h1 i)
  calc ∏ i ∈ s, f i ≤ (∏ i ∈ s, f i) * ∏ i ∈ t \ s, f i :=
        le_mul_of_one_le_right (by linarith) hd1
    _ = ∏ i ∈ t, f i := by
        rw [← Finset.prod_union Finset.disjoint_sdiff, Finset.union_sdiff_of_subset ht]

/-- ★対数の和が収束すれば、もとの族の無限積は `exp` である。 -/
theorem hasProd_of_hasSum_log' {ι : Type*} {f : ι → ℝ} {a : ℝ} (hfn : ∀ i, 0 < f i)
    (hf : HasSum (fun i ↦ Real.log (f i)) a) : HasProd f (Real.exp a) :=
  hf.rexp.congr (by simp [Real.exp_log, hfn])

/-- ★★**`|log(1 + A + R) − A| ≤ R + (A+R)²`**(`A, R ≥ 0`)。

★上からは `log(1+T) ≤ T`、下からは `log(1+T) ≥ 1 − (1+T)⁻¹ ≥ T − T²`。 -/
theorem abs_log_one_add_sub_le {A R : ℝ} (hA : 0 ≤ A) (hR : 0 ≤ R) :
    |Real.log (1 + (A + R)) - A| ≤ R + (A + R) ^ 2 := by
  have hT0 : (0 : ℝ) ≤ A + R := by linarith
  have hpos : (0 : ℝ) < 1 + (A + R) := by linarith
  have hup : Real.log (1 + (A + R)) ≤ A + R := by
    have := Real.log_le_sub_one_of_pos hpos
    linarith
  have hlow : 1 - (1 + (A + R))⁻¹ ≤ Real.log (1 + (A + R)) :=
    Real.one_sub_inv_le_log_of_pos hpos
  have hinv : (1 + (A + R))⁻¹ ≤ 1 - (A + R) + (A + R) ^ 2 := by
    rw [inv_eq_one_div, div_le_iff₀ hpos]
    nlinarith
  rw [abs_le]
  constructor <;> nlinarith

/-! ## ★1. 実の Dedekind ζ -/

variable (K : Type*) [Field K] [NumberField K] {s : ℝ}

/-- ★実の加数 `a_K(n)·n^{-s}`。 -/
noncomputable def zetaSummandR (K : Type*) [Field K] [NumberField K] (s : ℝ) (n : ℕ) : ℝ :=
  (idealCount K n : ℝ) * (n : ℝ) ^ (-s)

theorem zetaSummandR_nonneg (n : ℕ) : 0 ≤ zetaSummandR K s n :=
  mul_nonneg (Nat.cast_nonneg _) (Real.rpow_nonneg (Nat.cast_nonneg _) _)

theorem ofReal_zetaSummandR (n : ℕ) :
    ((zetaSummandR K s n : ℝ) : ℂ) = zetaSummand K (s : ℂ) n := by
  rw [zetaSummandR, zetaSummand, Complex.ofReal_mul, Complex.ofReal_natCast]
  congr 1
  rw [Complex.ofReal_cpow (Nat.cast_nonneg n) (-s)]
  push_cast
  ring

theorem norm_zetaSummand (n : ℕ) :
    ‖zetaSummand K (s : ℂ) n‖ = zetaSummandR K s n := by
  rw [← ofReal_zetaSummandR, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg (zetaSummandR_nonneg K n)]

theorem summable_zetaSummandR (hs : 1 < s) : Summable (zetaSummandR K s) := by
  have hre : (1 : ℝ) < ((s : ℂ)).re := by simpa using hs
  exact (summable_zetaSummand K hre).congr fun n ↦ norm_zetaSummand K n

/-- ★実の Dedekind ζ。 -/
noncomputable def zetaR (K : Type*) [Field K] [NumberField K] (s : ℝ) : ℝ :=
  ∑' n : ℕ, zetaSummandR K s n

theorem ofReal_zetaR (hs : 1 < s) : ((zetaR K s : ℝ) : ℂ) = dedekindZeta K s := by
  rw [zetaR, Complex.ofReal_tsum,
    ← tsum_zetaSummand K (show (1:ℝ) < ((s:ℂ)).re by simpa using hs)]
  exact tsum_congr fun n ↦ ofReal_zetaSummandR K n

@[simp] theorem zetaSummandR_one : zetaSummandR K s 1 = 1 := by
  rw [zetaSummandR, idealCount_one]
  simp

@[simp] theorem zetaSummandR_zero (hs : s ≠ 0) : zetaSummandR K s 0 = 0 := by
  rw [zetaSummandR, Nat.cast_zero, Real.zero_rpow (neg_ne_zero.mpr hs), mul_zero]

theorem zetaSummandR_mul (hs : s ≠ 0) {m n : ℕ} (h : Nat.Coprime m n) :
    zetaSummandR K s (m * n) = zetaSummandR K s m * zetaSummandR K s n := by
  rcases eq_or_ne m 0 with rfl | hm
  · rw [Nat.coprime_zero_left] at h
    subst h
    rw [zero_mul, zetaSummandR_zero K hs, zero_mul]
  rcases eq_or_ne n 0 with rfl | hn
  · rw [Nat.coprime_zero_right] at h
    subst h
    rw [one_mul, zetaSummandR_zero K hs, mul_zero]
  · rw [zetaSummandR, zetaSummandR, zetaSummandR, idealCount_mul K h hm hn,
      Nat.cast_mul, Nat.cast_mul, Real.mul_rpow (Nat.cast_nonneg m) (Nat.cast_nonneg n)]
    ring

theorem summable_norm_zetaSummandR (hs : 1 < s) :
    Summable (fun n ↦ ‖zetaSummandR K s n‖) := by
  refine (summable_zetaSummandR K hs).congr fun n ↦ ?_
  rw [Real.norm_eq_abs, abs_of_nonneg (zetaSummandR_nonneg K n)]

/-- ★★★★★**実の Euler 積**(`ℝ` は完備ノルム可換環なので同じ機構が使える)。 -/
theorem zetaR_eulerProduct_hasProd (hs : 1 < s) :
    HasProd (fun p : Nat.Primes ↦ ∑' e : ℕ, zetaSummandR K s (p ^ e)) (zetaR K s) :=
  EulerProduct.eulerProduct_hasProd (zetaSummandR_one K)
    (fun {_ _} h ↦ zetaSummandR_mul K (by linarith) h)
    (summable_norm_zetaSummandR K hs) (zetaSummandR_zero K (by linarith))

/-! ## ★2. `e ≥ 2` の項の一様評価 -/

theorem zetaSummandR_le_div (hs : 1 ≤ s) {n : ℕ} (hn : 1 ≤ n) :
    zetaSummandR K s n ≤ (idealCount K n : ℝ) / n := by
  have hx : (1 : ℝ) ≤ (n : ℝ) := by exact_mod_cast hn
  have h1 : ((n : ℝ)) ^ (-s) ≤ ((n : ℝ)) ^ (-1 : ℝ) :=
    Real.rpow_le_rpow_of_exponent_le hx (by linarith)
  rw [Real.rpow_neg_one] at h1
  rw [zetaSummandR, div_eq_mul_inv]
  exact mul_le_mul_of_nonneg_left h1 (Nat.cast_nonneg _)

/-- ★`e ≥ 2` の一様定数の被加数。 -/
noncomputable def tailTerm (K : Type*) [Field K] [NumberField K] (e : ℕ) : ℝ :=
  ((((e + 2) * Module.finrank ℤ (𝓞 K) + 1 : ℕ) : ℝ)) ^ Module.finrank ℤ (𝓞 K) * (1/2 : ℝ) ^ e

theorem tailTerm_nonneg (e : ℕ) : 0 ≤ tailTerm K e :=
  mul_nonneg (pow_nonneg (Nat.cast_nonneg _) _) (pow_nonneg (by norm_num) _)

/-- ★多項式 × 等比なので総和は収束する。 -/
theorem summable_tailTerm : Summable (tailTerm K) := by
  set d := Module.finrank ℤ (𝓞 K) with hd
  have hbase : Summable (fun e : ℕ => ((e : ℝ) + 1) ^ d * (1/2 : ℝ) ^ e) := by
    have h := summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) d
      (r := (1/2 : ℝ)) (by rw [Real.norm_eq_abs]; norm_num)
    have h2 := h.comp_injective (add_left_injective 1)
    have h3 : Summable (fun e : ℕ => (2:ℝ) * (((e+1 : ℕ) : ℝ) ^ d * (1/2 : ℝ) ^ (e+1))) :=
      h2.mul_left 2
    refine h3.congr fun e => ?_
    push_cast
    ring
  refine Summable.of_nonneg_of_le (tailTerm_nonneg K)
    (fun e => ?_) (hbase.mul_left (((3 * d + 1 : ℕ) : ℝ) ^ d))
  have hle : (((e + 2) * d + 1 : ℕ) : ℝ) ≤ ((3 * d + 1 : ℕ) : ℝ) * ((e : ℝ) + 1) := by
    push_cast
    nlinarith [Nat.cast_nonneg (α := ℝ) e, Nat.cast_nonneg (α := ℝ) d]
  calc tailTerm K e
      = (((e + 2) * d + 1 : ℕ) : ℝ) ^ d * (1/2 : ℝ) ^ e := rfl
    _ ≤ (((3 * d + 1 : ℕ) : ℝ) * ((e : ℝ) + 1)) ^ d * (1/2 : ℝ) ^ e := by
        refine mul_le_mul_of_nonneg_right (pow_le_pow_left₀ (Nat.cast_nonneg _) hle d) ?_
        positivity
    _ = ((3 * d + 1 : ℕ) : ℝ) ^ d * (((e : ℝ) + 1) ^ d * (1/2 : ℝ) ^ e) := by
        rw [mul_pow]; ring

/-- ★`e ≥ 2` の項の一様定数。 -/
noncomputable def tailConst (K : Type*) [Field K] [NumberField K] : ℝ :=
  ∑' e : ℕ, tailTerm K e

theorem tailConst_nonneg : 0 ≤ tailConst K :=
  tsum_nonneg (tailTerm_nonneg K)

/-- ★★`a(p^{e+2})·p^{-(e+2)s} ≤ tailTerm(e)/p²`。 -/
theorem zetaSummandR_pow_le (hs : 1 ≤ s) {p : ℕ} (hp : p.Prime) (e : ℕ) :
    zetaSummandR K s (p ^ (e + 2)) ≤ tailTerm K e * ((p : ℝ) ^ 2)⁻¹ := by
  have hp2 : 2 ≤ p := hp.two_le
  have hppos : (0 : ℝ) < (p : ℝ) := by
    have : (0 : ℕ) < p := hp.pos
    exact_mod_cast this
  have h1 : zetaSummandR K s (p ^ (e + 2))
      ≤ (idealCount K (p ^ (e + 2)) : ℝ) / ((p ^ (e + 2) : ℕ) : ℝ) :=
    zetaSummandR_le_div K hs (Nat.one_le_pow _ _ hp.pos)
  rw [Nat.cast_pow] at h1
  have h2 : (idealCount K (p ^ (e + 2)) : ℝ)
      ≤ (((e + 2) * Module.finrank ℤ (𝓞 K) + 1 : ℕ) : ℝ) ^ Module.finrank ℤ (𝓞 K) := by
    have := idealCount_prime_pow_le K hp (e + 2)
    exact_mod_cast this
  have h3 : ((2 : ℝ) ^ e * (p : ℝ) ^ 2) ≤ ((p : ℝ) ^ (e + 2)) := by
    rw [pow_add]
    refine mul_le_mul_of_nonneg_right ?_ (by positivity)
    exact pow_le_pow_left₀ (by norm_num) (by exact_mod_cast hp2) e
  have hden : (0 : ℝ) < (2 : ℝ) ^ e * (p : ℝ) ^ 2 := by positivity
  calc zetaSummandR K s (p ^ (e + 2))
      ≤ (idealCount K (p ^ (e + 2)) : ℝ) / ((p : ℝ) ^ (e + 2)) := h1
    _ ≤ ((((e + 2) * Module.finrank ℤ (𝓞 K) + 1 : ℕ) : ℝ) ^ Module.finrank ℤ (𝓞 K))
          / ((2 : ℝ) ^ e * (p : ℝ) ^ 2) := by
        rw [div_le_div_iff₀ (lt_of_lt_of_le hden h3) hden]
        nlinarith [pow_nonneg (le_of_lt hppos) (e + 2)]
    _ = tailTerm K e * ((p : ℝ) ^ 2)⁻¹ := by
        rw [tailTerm, div_pow, one_pow]
        field_simp

/-! ## ★3. 局所因子 -/

/-- ★局所因子 `L_p(s) = Σ_e a(p^e) p^{-es}`。 -/
noncomputable def localFactor (K : Type*) [Field K] [NumberField K] (s : ℝ) (p : ℕ) : ℝ :=
  ∑' e : ℕ, zetaSummandR K s (p ^ e)

theorem summable_localFactor (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    Summable (fun e : ℕ ↦ zetaSummandR K s (p ^ e)) :=
  (summable_zetaSummandR K hs).comp_injective (Nat.pow_right_injective hp)

theorem localFactor_eq_one_add (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    localFactor K s p = 1 + ∑' e : ℕ, zetaSummandR K s (p ^ (e + 1)) := by
  rw [localFactor, (summable_localFactor K hs hp).tsum_eq_zero_add]
  simp

theorem summable_shift1 (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    Summable (fun e : ℕ ↦ zetaSummandR K s (p ^ (e + 1))) :=
  (summable_localFactor K hs hp).comp_injective (add_left_injective 1)

theorem summable_shift2 (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    Summable (fun e : ℕ ↦ zetaSummandR K s (p ^ (e + 2))) :=
  (summable_localFactor K hs hp).comp_injective (add_left_injective 2)

theorem tsum_shift_eq (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    ∑' e : ℕ, zetaSummandR K s (p ^ (e + 1))
      = zetaSummandR K s p + ∑' e : ℕ, zetaSummandR K s (p ^ (e + 2)) := by
  rw [(summable_shift1 K hs hp).tsum_eq_zero_add]
  simp [pow_one]

theorem localFactor_eq (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    localFactor K s p
      = 1 + (zetaSummandR K s p + ∑' e : ℕ, zetaSummandR K s (p ^ (e + 2))) := by
  rw [localFactor_eq_one_add K hs hp, tsum_shift_eq K hs hp]

theorem one_le_localFactor (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) : 1 ≤ localFactor K s p := by
  rw [localFactor]
  have h := (summable_localFactor K hs hp).le_tsum 0
    (fun e _ => zetaSummandR_nonneg K (p ^ e))
  simpa using h

theorem hasProd_localFactor (hs : 1 < s) :
    HasProd (fun p : Nat.Primes ↦ localFactor K s p) (zetaR K s) :=
  zetaR_eulerProduct_hasProd K hs

theorem localFactor_pos (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) : 0 < localFactor K s p :=
  lt_of_lt_of_le zero_lt_one (one_le_localFactor K hs hp)

theorem one_le_zetaR (hs : 1 < s) : 1 ≤ zetaR K s := by
  have h := prod_le_of_hasProd_of_one_le (hasProd_localFactor K hs)
    (fun p : Nat.Primes => one_le_localFactor K hs p.2.two_le) ∅
  simpa using h

/-! ## ★4. `Σ_p log L_p(s) = log ζ_L(s)` -/

/-- ★★部分和が `log ζ_L(s)` で一様に押さえられるので、対数の和は収束する。 -/
theorem summable_log_localFactor (hs : 1 < s) :
    Summable (fun p : Nat.Primes ↦ Real.log (localFactor K s p)) := by
  refine summable_of_sum_le (c := Real.log (zetaR K s))
    (fun p => Real.log_nonneg (one_le_localFactor K hs p.2.two_le)) (fun t => ?_)
  have hprodpos : (0 : ℝ) < ∏ p ∈ t, localFactor K s p :=
    Finset.prod_pos (fun p _ => localFactor_pos K hs p.2.two_le)
  have hlog : Real.log (∏ p ∈ t, localFactor K s p)
      = ∑ p ∈ t, Real.log (localFactor K s p) :=
    Real.log_prod (fun p _ => ne_of_gt (localFactor_pos K hs p.2.two_le))
  rw [← hlog]
  exact Real.log_le_log hprodpos (prod_le_of_hasProd_of_one_le (hasProd_localFactor K hs)
    (fun p : Nat.Primes => one_le_localFactor K hs p.2.two_le) t)

/-- ★★★★**`Σ_p log L_p(s) = log ζ_L(s)`**。 -/
theorem tsum_log_localFactor (hs : 1 < s) :
    ∑' p : Nat.Primes, Real.log (localFactor K s p) = Real.log (zetaR K s) := by
  have hsum := summable_log_localFactor K hs
  have h1 : HasProd (fun p : Nat.Primes ↦ localFactor K s p)
      (Real.exp (∑' p : Nat.Primes, Real.log (localFactor K s p))) :=
    hasProd_of_hasSum_log' (fun p => localFactor_pos K hs p.2.two_le) hsum.hasSum
  have h2 := (hasProd_localFactor K hs).unique h1
  rw [h2, Real.log_exp]

/-! ## ★5. 素数ごとの一様評価 -/

theorem tail2_le (hs : 1 < s) {p : ℕ} (hp : p.Prime) :
    ∑' e : ℕ, zetaSummandR K s (p ^ (e + 2)) ≤ tailConst K * ((p : ℝ) ^ 2)⁻¹ := by
  have h := Summable.tsum_mono (summable_shift2 K hs hp.two_le)
    ((summable_tailTerm K).mul_right ((p : ℝ) ^ 2)⁻¹)
    (fun e => zetaSummandR_pow_le K (le_of_lt hs) hp e)
  rwa [tsum_mul_right, ← tailConst] at h

/-- ★`a(p)` の一様上界に使う定数。 -/
noncomputable def primeConst (K : Type*) [Field K] [NumberField K] : ℝ :=
  (((1 * Module.finrank ℤ (𝓞 K) + 1 : ℕ)) : ℝ) ^ Module.finrank ℤ (𝓞 K)

theorem primeConst_nonneg : 0 ≤ primeConst K := by
  rw [primeConst]; positivity

theorem zetaSummandR_prime_le (hs : 1 ≤ s) {p : ℕ} (hp : p.Prime) :
    zetaSummandR K s p ≤ primeConst K * ((p : ℝ))⁻¹ := by
  have h1 := zetaSummandR_le_div K hs (le_of_lt hp.one_lt)
  have h2 : (idealCount K p : ℝ) ≤ primeConst K := by
    have := idealCount_prime_pow_le K hp 1
    rw [pow_one] at this
    rw [primeConst]
    exact_mod_cast this
  have hppos : (0 : ℝ) < (p : ℝ) := by exact_mod_cast hp.pos
  calc zetaSummandR K s p ≤ (idealCount K p : ℝ) / (p : ℝ) := h1
    _ ≤ primeConst K / (p : ℝ) := div_le_div_of_nonneg_right h2 (le_of_lt hppos)
    _ = primeConst K * ((p : ℝ))⁻¹ := div_eq_mul_inv _ _

/-- ★`L_p(s) − 1` の一様上界に使う定数。 -/
noncomputable def mainConst (K : Type*) [Field K] [NumberField K] : ℝ :=
  primeConst K + tailConst K

theorem mainConst_nonneg : 0 ≤ mainConst K :=
  add_nonneg (primeConst_nonneg K) (tailConst_nonneg K)

theorem inv_sq_le_inv {p : ℕ} (hp : 1 ≤ p) : ((p : ℝ) ^ 2)⁻¹ ≤ ((p : ℝ))⁻¹ := by
  have hp1 : (1 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  rw [inv_le_inv₀ (by positivity) (by linarith)]
  nlinarith

theorem localFactor_sub_one_eq (hs : 1 < s) {p : ℕ} (hp : 2 ≤ p) :
    localFactor K s p - 1
      = zetaSummandR K s p + ∑' e : ℕ, zetaSummandR K s (p ^ (e + 2)) := by
  rw [localFactor_eq K hs hp]; ring

theorem localFactor_sub_one_le (hs : 1 < s) {p : ℕ} (hp : p.Prime) :
    localFactor K s p - 1 ≤ mainConst K * ((p : ℝ))⁻¹ := by
  have hA := zetaSummandR_prime_le K (le_of_lt hs) hp
  have hR := tail2_le K hs hp
  have hle : tailConst K * ((p : ℝ) ^ 2)⁻¹ ≤ tailConst K * ((p : ℝ))⁻¹ :=
    mul_le_mul_of_nonneg_left (inv_sq_le_inv hp.one_lt.le) (tailConst_nonneg K)
  rw [localFactor_sub_one_eq K hs hp.two_le, mainConst, add_mul]
  linarith

/-- ★★★**素数ごとの一様評価** `|log L_p(s) − a(p)p^{-s}| ≤ C/p²`。 -/
theorem abs_log_localFactor_sub_le (hs : 1 < s) {p : ℕ} (hp : p.Prime) :
    |Real.log (localFactor K s p) - zetaSummandR K s p|
      ≤ (tailConst K + mainConst K ^ 2) * ((p : ℝ) ^ 2)⁻¹ := by
  have hA0 : 0 ≤ zetaSummandR K s p := zetaSummandR_nonneg K p
  have hR0 : 0 ≤ ∑' e : ℕ, zetaSummandR K s (p ^ (e + 2)) :=
    tsum_nonneg (fun e => zetaSummandR_nonneg K _)
  have hbase := abs_log_one_add_sub_le hA0 hR0
  rw [← localFactor_eq K hs hp.two_le, ← localFactor_sub_one_eq K hs hp.two_le] at hbase
  have hR := tail2_le K hs hp
  have hT := localFactor_sub_one_le K hs hp
  have hT0 : 0 ≤ localFactor K s p - 1 := by linarith [one_le_localFactor K hs hp.two_le]
  have hsq : (localFactor K s p - 1) ^ 2 ≤ mainConst K ^ 2 * ((p : ℝ) ^ 2)⁻¹ := by
    calc (localFactor K s p - 1) ^ 2 ≤ (mainConst K * ((p : ℝ))⁻¹) ^ 2 :=
          pow_le_pow_left₀ hT0 hT 2
      _ = mainConst K ^ 2 * ((p : ℝ) ^ 2)⁻¹ := by rw [mul_pow, ← inv_pow]
  calc |Real.log (localFactor K s p) - zetaSummandR K s p|
      ≤ (∑' e : ℕ, zetaSummandR K s (p ^ (e + 2))) + (localFactor K s p - 1) ^ 2 := hbase
    _ ≤ tailConst K * ((p : ℝ) ^ 2)⁻¹ + mainConst K ^ 2 * ((p : ℝ) ^ 2)⁻¹ := by linarith
    _ = (tailConst K + mainConst K ^ 2) * ((p : ℝ) ^ 2)⁻¹ := by ring

/-! ## ★6. `log ζ_L(s) = Σ_p a(p)p^{-s} + O(1)` -/

theorem summable_primes_inv_sq : Summable (fun p : Nat.Primes ↦ ((p : ℝ) ^ 2)⁻¹) :=
  (Real.summable_nat_pow_inv.mpr (by norm_num)).comp_injective
    (fun _ _ h => Subtype.ext h)

theorem summable_zetaSummandR_primes (hs : 1 < s) :
    Summable (fun p : Nat.Primes ↦ zetaSummandR K s p) :=
  (summable_zetaSummandR K hs).comp_injective (fun _ _ h => Subtype.ext h)

/-- ★★`log ζ_L(s) − Σ_p a(p)p^{-s}` の一様上界。 -/
noncomputable def logConst (K : Type*) [Field K] [NumberField K] : ℝ :=
  (tailConst K + mainConst K ^ 2) * ∑' p : Nat.Primes, ((p : ℝ) ^ 2)⁻¹

/-- ★★★★★★**[cheb-log-zeta] `log ζ_L(s) = Σ_p a_K(p)·p^{-s} + O(1)`**
—— `O(1)` は `s > 1` について**一様**。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★`a_K(p) = #{𝔭 : N𝔭 = p} = #{𝔭 ∣ p : f(𝔭) = 1}` なので、これは原文の
`log ζ_L(s) = Σ_{f(𝔭)=1} N𝔭^{-s} + O(1)` そのものである。 -/
theorem abs_log_zetaR_sub_tsum_le (hs : 1 < s) :
    |Real.log (zetaR K s) - ∑' p : Nat.Primes, zetaSummandR K s p| ≤ logConst K := by
  have hsum1 := summable_log_localFactor K hs
  have hsum2 := summable_zetaSummandR_primes K hs
  have habs : Summable (fun p : Nat.Primes ↦
      ‖Real.log (localFactor K s p) - zetaSummandR K s p‖) := by
    refine Summable.of_nonneg_of_le (fun p => norm_nonneg _) (fun p => ?_)
      ((summable_primes_inv_sq).mul_left (tailConst K + mainConst K ^ 2))
    rw [Real.norm_eq_abs]
    exact abs_log_localFactor_sub_le K hs p.2
  rw [← tsum_log_localFactor K hs, ← hsum1.tsum_sub hsum2]
  calc ‖∑' p : Nat.Primes, (Real.log (localFactor K s p) - zetaSummandR K s p)‖
      ≤ ∑' p : Nat.Primes, ‖Real.log (localFactor K s p) - zetaSummandR K s p‖ :=
        norm_tsum_le_tsum_norm habs
    _ ≤ ∑' p : Nat.Primes, (tailConst K + mainConst K ^ 2) * ((p : ℝ) ^ 2)⁻¹ := by
        refine Summable.tsum_mono habs ((summable_primes_inv_sq).mul_left _) (fun p => ?_)
        rw [Real.norm_eq_abs]
        exact abs_log_localFactor_sub_le K hs p.2
    _ = logConst K := by rw [logConst, tsum_mul_left]

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Theorem 6.4, (iv)` が使う Tchebotarev の二段目。 -/
def abs_log_zetaR_sub_tsum_le.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — log ζ_L(s) = Σ_p a(p)p^{-s} + O(1)",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
