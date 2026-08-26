import ABC3.Found.GaloisRep.WeierstrassQExp

/-!
# Galois (G6) 第 219 ブロック —— **★★★★★★★`℘'` の q 展開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★段 5 のもう半分

第 218 で `℘` の q 展開が取れた。Tate の方程式には `℘'` も要る。
同じ手順(段に切る → 各段に Lipschitz 公式 → 段の和の可和性)を**重み 3** で繰り返す。

★★`℘` と違って**補正項が要らない**——`∑_l 1/(z−l)³` は絶対収束するので、
`℘` のときの `−1/l²` に当たるものが最初から無い。定数項も出ない。

## ★★★★★重み 3 の Lipschitz 公式は `h = f + 2g` を生む

mathlib の `qExpansion_identity` を `k = 2` で使うと右辺は `∑_{n≥0} n² tⁿ` になる。
これを閉じた形に直すと:

    ∑_{n≥0} n² tⁿ = t(1+t)/(1−t)³ = f(t) + 2g(t)     (`f = tateXterm`, `g = tateYterm`)

★`n² = 2·C(n+2,2) − 3·C(n+1,1) + 1` と分けて `hasSum_choose_mul_geometric_of_norm_lt_one`
を 2 回使えば出る。★★この `h := f + 2g` が **`℘'` の側に出る項**である(`tateDterm`)。

## ★★★★★★符号は `h` の反転で吸収される

奇数べきなので折り返しで符号が変わる(`∑_m 1/(w+m)³ = −∑_m 1/(−w+m)³`)。
そのため上向きの段と下向きの段は**符号が逆**に出る。一方

    h(1/t) = −h(t)                                    (`tateDterm_inv`)

なので、**両者はちょうど 1 つの形にまとまる**(`rowDP`):

    ∑_m −2/(z − (nτ+m))³ = (2πi)³ · h(q^{−n} u)

★`f` は反転で不変(`tateXterm_inv`)、`h` は反転で符号が変わる——
偶数べきと奇数べきの違いがそのまま出ている。

## ★★到達点

    ℘'(z) = (2πi)³ · ∑_{n∈ℤ} h(q^{−n} u)

| 定理 | 内容 |
|---|---|
| `tateDterm` | `h(t) = f(t) + 2g(t)` |
| `tateDterm_inv` | ★★★★★**`h(1/t) = −h(t)`** |
| `hasSum_sq_mul_geometric` | ★★★★★`∑_{n≥0} n² tⁿ = h(t)` |
| `lipschitz_three` | ★★★★★**Lipschitz 公式(重み 3)** |
| `tsum_inv_cube_neg` | ★★★★奇数べきの折り返しは符号を変える |
| `rowDP` | ★★★★★★**段の値の一様形** |
| `summable_int_tateDterm` | ★★★★★★段の和は可和 |
| `derivWeierstrassP_qExpansion` | ★★★★★★★**`℘'` の q 展開** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real Filter PeriodPair

/-! ## ★`℘'` の側に出る項 `h = f + 2g` -/

/-- ★ℂ の上では `tateYterm t = t²/(1−t)³`。 -/
theorem tateYterm_field (t : ℂ) : tateYterm t = t ^ 2 / (1 - t) ^ 3 := by
  simp [tateYterm, Ring.inverse_eq_inv', div_eq_mul_inv, inv_pow]

/-- ★`℘'` の側に出る項 `h(t) = f(t) + 2g(t) = t(1+t)/(1−t)³`。 -/
noncomputable def tateDterm {R : Type} [CommRing R] (t : R) : R :=
  tateXterm t + 2 * tateYterm t

theorem tateDterm_field (t : ℂ) : tateDterm t = (t + t ^ 2) / (1 - t) ^ 3 := by
  rcases eq_or_ne t 1 with rfl | ht
  · simp [tateDterm, tateXterm, tateYterm]
  have h1 : (1 : ℂ) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht)
  simp only [tateDterm, tateXterm_field, tateYterm_field]
  field_simp
  ring

/-- ★★★★★**`h` は反転で符号が変わる**——`h(1/t) = −h(t)`。

★`f` は反転で不変(`tateXterm_inv`)だが `h` は符号を変える。
これで `℘'` の上向きの段と下向きの段(符号が逆)が 1 つの形にまとまる。 -/
theorem tateDterm_inv {t : ℂ} (ht0 : t ≠ 0) (ht1 : t ≠ 1) :
    tateDterm t⁻¹ = -tateDterm t := by
  have h1 : (1 : ℂ) - t ≠ 0 := sub_ne_zero.2 (Ne.symm ht1)
  simp only [tateDterm]
  rw [tateXterm_inv ht0 ht1, tateYterm_inv ht0 ht1, tateXterm_field, tateYterm_field]
  field_simp
  ring

theorem tateDterm_one_complex : tateDterm (1 : ℂ) = 0 := by
  simp [tateDterm, tateXterm, tateYterm]

/-- ★★★★★**`∑_{n≥0} n² tⁿ = h(t)`**——重み 3 の Lipschitz 公式の右辺がこの形になる。

★`n² = 2·C(n+2,2) − 3·C(n+1,1) + 1` と分けて `choose` の幾何級数を 2 回使う。 -/
theorem hasSum_sq_mul_geometric {r : ℂ} (hr : ‖r‖ < 1) :
    HasSum (fun n : ℕ => (n : ℂ) ^ 2 * r ^ n) (tateDterm r) := by
  have hr1 : (1 : ℂ) - r ≠ 0 := by
    intro h
    have hh : r = 1 := by linear_combination -h
    rw [hh] at hr
    simp at hr
  have h2 := (hasSum_choose_mul_geometric_of_norm_lt_one 2 hr).mul_left 2
  have h1 := (hasSum_choose_mul_geometric_of_norm_lt_one 1 hr).mul_left 3
  have h0 := hasSum_geometric_of_norm_lt_one hr
  have hsum := (h2.sub h1).add h0
  have hval : (2 : ℂ) * (1 / (1 - r) ^ (2 + 1)) - 3 * (1 / (1 - r) ^ (1 + 1)) + (1 - r)⁻¹
      = tateDterm r := by
    rw [tateDterm_field]
    field_simp
    ring
  rw [hval] at hsum
  refine hsum.congr_fun fun n => ?_
  have hc2 : ((((n + 2).choose 2 : ℕ)) : ℂ) = ((n : ℂ) + 2) * ((n : ℂ) + 1) / 2 := by
    rw [Nat.cast_choose_two]
    push_cast
    ring
  have hc1 : ((((n + 1).choose 1 : ℕ)) : ℂ) = (n : ℂ) + 1 := by
    rw [Nat.choose_one_right]
    push_cast
    ring
  simp only [hc2, hc1]
  ring

/-! ## ★★★★★重み 3 の Lipschitz 公式と段への当てはめ -/

/-- ★★★★★**Lipschitz 公式(重み 3)**——右辺が `h = f + 2g` の形になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem lipschitz_three (z : UpperHalfPlane) :
    (∑' m : ℤ, 1 / ((z : ℂ) + (m : ℂ)) ^ 3)
      = -((2 * ↑π * I) ^ 3 / 2) * tateDterm (Complex.exp (2 * ↑π * I * z)) := by
  have hnorm : ‖Complex.exp (2 * ↑π * I * (z : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one z
  have h := EisensteinSeries.qExpansion_identity (k := 2) (by norm_num) z
  rw [h, (hasSum_sq_mul_geometric hnorm).tsum_eq]
  congr 1
  norm_num [Nat.factorial]
  ring

/-- ★★★★**奇数べきは折り返しで符号が変わる**——`∑_m 1/(w+m)³ = −∑_m 1/(−w+m)³`。 -/
theorem tsum_inv_cube_neg (w : ℂ) :
    (∑' m : ℤ, 1 / (w + (m : ℂ)) ^ 3) = -∑' m : ℤ, 1 / (-w + (m : ℂ)) ^ 3 := by
  rw [← (Equiv.neg ℤ).tsum_eq (fun m : ℤ => 1 / (-w + (m : ℂ)) ^ 3), ← tsum_neg]
  refine tsum_congr fun m => ?_
  simp only [Equiv.neg_apply, Int.cast_neg]
  rw [show (-w + -(m : ℂ)) ^ 3 = -((w + m) ^ 3) by ring]
  rw [div_neg, neg_neg]

/-- ★`∑_m 1/(w−m)³ = ∑_m 1/(w+m)³`(`m` だけを反転するので符号は出ない)。 -/
theorem tsum_one_div_cube_sub_eq_add (w : ℂ) :
    (∑' m : ℤ, 1 / (w - (m : ℂ)) ^ 3) = ∑' m : ℤ, 1 / (w + (m : ℂ)) ^ 3 := by
  rw [← (Equiv.neg ℤ).tsum_eq (fun m : ℤ => 1 / (w + (m : ℂ)) ^ 3)]
  refine tsum_congr fun m => ?_
  simp only [Equiv.neg_apply, Int.cast_neg]
  rw [show w + -(m : ℂ) = w - m by ring]

/-- ★★★★★**上向きの段の Lipschitz 公式(重み 3)**。 -/
theorem lipschitz_shift_three (z τ : UpperHalfPlane) (n : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) + (n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 3)
      = -((2 * ↑π * I) ^ 3 / 2) *
        tateDterm (Complex.exp (2 * ↑π * I * τ) ^ n * Complex.exp (2 * ↑π * I * z)) := by
  have h := lipschitz_three (shiftUp z τ n)
  rw [shiftUp_coe] at h
  rw [h, exp_add_nat_mul]

/-- ★★★★★★**下向きの段の Lipschitz 公式(重み 3)**——奇数べきなので符号が反転する。 -/
theorem lipschitz_shift_neg_three (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im)
    (n : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - ((n : ℕ) + 1 : ℕ) * (τ : ℂ) + (m : ℂ)) ^ 3)
      = ((2 * ↑π * I) ^ 3 / 2) *
        tateDterm (Complex.exp (2 * ↑π * I * τ) ^ (n + 1)
          * (Complex.exp (2 * ↑π * I * z))⁻¹) := by
  have hrefl := tsum_inv_cube_neg ((z : ℂ) - ((n : ℕ) + 1 : ℕ) * (τ : ℂ))
  have h := lipschitz_three (shiftDown z τ him n)
  rw [shiftDown_coe] at h
  rw [hrefl, show -((z : ℂ) - ((n : ℕ) + 1 : ℕ) * (τ : ℂ))
      = ((n : ℕ) + 1 : ℕ) * (τ : ℂ) - (z : ℂ) by ring, h, exp_sub_nat_mul]
  ring

/-! ## ★★★★★`℘'` を段に切る -/

/-- ★★★★**`℘'` の和を `ℤ × ℤ` の上に移す**。 -/
theorem hasSum_derivWeierstrassP_prod (τ : UpperHalfPlane) (z : ℂ) :
    HasSum (fun x : ℤ × ℤ =>
      -2 / (z - ((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ))) ^ 3)
      (PeriodPair.derivWeierstrassP (tauPair τ) z) := by
  have h := (tauPair τ).hasSum_derivWeierstrassP z
  rw [← Equiv.hasSum_iff (tauPair τ).latticeEquivProd.toEquiv.symm] at h
  refine h.congr_fun fun x => ?_
  simp only [Function.comp_apply]
  have heq : (tauPair τ).latticeEquivProd.toEquiv.symm x
      = (tauPair τ).latticeEquivProd.symm x := rfl
  rw [heq, tauPair_symm_apply]

/-- ★★★★★**`℘'` を τ の段に切る**。 -/
theorem derivWeierstrassP_eq_tsum_rows (τ : UpperHalfPlane) (z : ℂ) :
    PeriodPair.derivWeierstrassP (tauPair τ) z
      = ∑' n : ℤ, ∑' m : ℤ, -2 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 3 := by
  have h := hasSum_derivWeierstrassP_prod τ z
  rw [← h.tsum_eq, h.summable.tsum_prod]

/-- ★段の和から定数 `−2` を外に出す。 -/
theorem rowDP_eq_smul (τ : UpperHalfPlane) (z : ℂ) (n : ℤ) :
    (∑' m : ℤ, -2 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 3)
      = -2 * ∑' m : ℤ, 1 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 3 := by
  rw [← tsum_mul_left]
  refine tsum_congr fun m => ?_
  rw [mul_one_div]

/-- ★★★★★**`n ≤ 0` の段の `℘'` 側**。 -/
theorem rowDP_nonpos (z τ : UpperHalfPlane) (k : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - (((-(k : ℤ)) : ℤ) * (τ : ℂ) + (m : ℂ))) ^ 3)
      = -((2 * ↑π * I) ^ 3 / 2) * tateDterm
        (Complex.exp (2 * ↑π * I * τ) ^ k * Complex.exp (2 * ↑π * I * z)) := by
  have hcong : ∀ m : ℤ,
      1 / ((z : ℂ) - (((-(k : ℤ)) : ℤ) * (τ : ℂ) + (m : ℂ))) ^ 3
        = 1 / (((z : ℂ) + (k : ℂ) * (τ : ℂ)) - (m : ℂ)) ^ 3 := by
    intro m
    push_cast
    ring_nf
  rw [tsum_congr hcong, tsum_one_div_cube_sub_eq_add]
  exact lipschitz_shift_three z τ k

/-- ★★★★★★**`n ≥ 1` の段の `℘'` 側**——符号が逆になる。 -/
theorem rowDP_pos (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (j : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - (((j : ℤ) + 1) * (τ : ℂ) + (m : ℂ))) ^ 3)
      = ((2 * ↑π * I) ^ 3 / 2) * tateDterm
        (Complex.exp (2 * ↑π * I * τ) ^ (j + 1) * (Complex.exp (2 * ↑π * I * z))⁻¹) := by
  have hcong : ∀ m : ℤ,
      1 / ((z : ℂ) - (((j : ℤ) + 1) * (τ : ℂ) + (m : ℂ))) ^ 3
        = 1 / (((z : ℂ) - ((j : ℕ) + 1 : ℕ) * (τ : ℂ)) - (m : ℂ)) ^ 3 := by
    intro m
    push_cast
    ring_nf
  rw [tsum_congr hcong, tsum_one_div_cube_sub_eq_add]
  exact lipschitz_shift_neg_three z τ him j

/-- ★★★★★★**`℘'` 側の段の値の一様形**——`n` 段目は `−((2πi)³/2)·h(q^{−n}u)`。

★上向きと下向きは符号が逆だが、`h` が反転で符号を変える(`tateDterm_inv`)ので
1 つの形にまとまる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem rowDP (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (n : ℤ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 3)
      = -((2 * ↑π * I) ^ 3 / 2) * tateDterm
        (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z)) := by
  have hq : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hu : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  rcases le_or_gt n 0 with hn0 | hn0
  · obtain ⟨k, rfl⟩ : ∃ k : ℕ, n = -(k : ℤ) := ⟨n.natAbs, by omega⟩
    rw [rowDP_nonpos z τ k]
    congr 2
    rw [neg_neg, zpow_natCast]
  · obtain ⟨j, rfl⟩ : ∃ j : ℕ, n = (j : ℤ) + 1 := ⟨(n - 1).toNat, by omega⟩
    have h := rowDP_pos z τ him j
    push_cast at h ⊢
    have hn : ‖Complex.exp (2 * ↑π * I * τ) ^ (j + 1)
        * (Complex.exp (2 * ↑π * I * z))⁻¹‖ < 1 := norm_down_lt_one z τ him j
    have hne1 : Complex.exp (2 * ↑π * I * τ) ^ (j + 1)
        * (Complex.exp (2 * ↑π * I * z))⁻¹ ≠ 1 := by
      intro hc
      rw [hc] at hn
      simp at hn
    have hne0 : Complex.exp (2 * ↑π * I * τ) ^ (j + 1)
        * (Complex.exp (2 * ↑π * I * z))⁻¹ ≠ 0 :=
      mul_ne_zero (pow_ne_zero _ hq) (inv_ne_zero hu)
    have hinv : (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (j + 1)
        * (Complex.exp (2 * ↑π * I * (z : ℂ)))⁻¹)⁻¹
          = Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-((j : ℤ) + 1))
            * Complex.exp (2 * ↑π * I * (z : ℂ)) := by
      rw [mul_inv, inv_inv, ← zpow_natCast (Complex.exp (2 * ↑π * I * (τ : ℂ))) (j + 1),
        ← zpow_neg]
      push_cast
      ring_nf
    have hd := tateDterm_inv hne0 hne1
    rw [hinv] at hd
    rw [h, hd]
    ring

/-! ## ★★★★★段の和の可和性 -/

theorem norm_tateDterm (t : ℂ) : ‖tateDterm t‖ = ‖t + t ^ 2‖ / ‖1 - t‖ ^ 3 := by
  rw [tateDterm_field, norm_div, norm_pow]

/-- ★★★★`‖t‖ ≤ 1/2` なら `‖h(t)‖ ≤ 12‖t‖`。 -/
theorem norm_tateDterm_le_of_small {t : ℂ} (ht : ‖t‖ ≤ 1 / 2) :
    ‖tateDterm t‖ ≤ 12 * ‖t‖ := by
  have h1 : (1 : ℝ) / 2 ≤ ‖1 - t‖ := by
    have := norm_sub_norm_le (1 : ℂ) t
    simp only [norm_one] at this
    linarith
  have h2 : (0 : ℝ) < ‖1 - t‖ := by linarith
  have hnum : ‖t + t ^ 2‖ ≤ ‖t‖ + ‖t‖ ^ 2 := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_pow]
  have hcube : (1 / 8 : ℝ) ≤ ‖1 - t‖ ^ 3 := by
    have hp := pow_le_pow_left₀ (by norm_num : (0:ℝ) ≤ 1 / 2) h1 3
    norm_num at hp
    linarith
  have hsq : ‖t‖ ^ 2 ≤ ‖t‖ / 2 := by nlinarith [norm_nonneg t]
  rw [norm_tateDterm, div_le_iff₀ (by positivity)]
  nlinarith [mul_nonneg (norm_nonneg t) (by linarith : (0:ℝ) ≤ ‖1 - t‖ ^ 3 - 1 / 8)]

/-- ★★★★`2 ≤ ‖t‖` なら `‖h(t)‖ ≤ 12/‖t‖`。 -/
theorem norm_tateDterm_le_of_large {t : ℂ} (ht : 2 ≤ ‖t‖) :
    ‖tateDterm t‖ ≤ 12 / ‖t‖ := by
  have h1 : ‖t‖ / 2 ≤ ‖1 - t‖ := by
    have := norm_sub_norm_le t (1 : ℂ)
    simp only [norm_one] at this
    rw [show (1 : ℂ) - t = -(t - 1) by ring, norm_neg]
    linarith
  have h2 : (0 : ℝ) < ‖t‖ := by linarith
  have h3 : (0 : ℝ) < ‖1 - t‖ := by linarith
  have hnum : ‖t + t ^ 2‖ ≤ ‖t‖ + ‖t‖ ^ 2 := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_pow]
  have hnum2 : ‖t + t ^ 2‖ ≤ 3 / 2 * ‖t‖ ^ 2 := by nlinarith
  rw [norm_tateDterm, div_le_div_iff₀ (by positivity) h2]
  nlinarith [pow_le_pow_left₀ (by positivity : (0:ℝ) ≤ ‖t‖ / 2) h1 3, norm_nonneg (t + t ^ 2)]

theorem summable_tateDterm_small (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun k : ℕ => tateDterm (q ^ ((k : ℤ)) * c) := by
  have hten : Tendsto (fun k : ℕ => ‖q‖ ^ k * ‖c‖) atTop (nhds 0) := by
    have := tendsto_pow_atTop_nhds_zero_of_lt_one (norm_nonneg q) hq
    simpa using this.mul_const ‖c‖
  have hev : ∀ᶠ k : ℕ in atTop, ‖q‖ ^ k * ‖c‖ ≤ 1 / 2 :=
    hten.eventually_le_const (by norm_num)
  refine Summable.of_norm_bounded_eventually
    (g := fun k : ℕ => 12 * ‖c‖ * ‖q‖ ^ k) ?_ ?_
  · exact ((summable_geometric_of_lt_one (norm_nonneg q) hq).mul_left (12 * ‖c‖))
  · rw [Nat.cofinite_eq_atTop]
    filter_upwards [hev] with k hk
    have hn : ‖q ^ ((k : ℤ)) * c‖ = ‖q‖ ^ k * ‖c‖ := by
      rw [norm_mul, zpow_natCast, norm_pow]
    have h := norm_tateDterm_le_of_small (t := q ^ ((k : ℤ)) * c) (by rw [hn]; exact hk)
    rw [hn] at h
    linarith [h]

theorem summable_tateDterm_large (q c : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) (hc0 : c ≠ 0) :
    Summable fun k : ℕ => tateDterm (q ^ (-(k : ℤ)) * c) := by
  have hqpos : (0 : ℝ) < ‖q‖ := norm_pos_iff.2 hq0
  have hcpos : (0 : ℝ) < ‖c‖ := norm_pos_iff.2 hc0
  have hnorm : ∀ k : ℕ, ‖q ^ (-(k : ℤ)) * c‖ = ‖c‖ / ‖q‖ ^ k := by
    intro k
    rw [norm_mul, norm_zpow, zpow_neg, zpow_natCast, div_eq_mul_inv, mul_comm]
  have hten : Tendsto (fun k : ℕ => ‖q‖ ^ k) atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (norm_nonneg q) hq
  have hev : ∀ᶠ k : ℕ in atTop, ‖q‖ ^ k ≤ ‖c‖ / 2 :=
    hten.eventually_le_const (by positivity)
  refine Summable.of_norm_bounded_eventually
    (g := fun k : ℕ => 12 / ‖c‖ * ‖q‖ ^ k) ?_ ?_
  · exact ((summable_geometric_of_lt_one (norm_nonneg q) hq).mul_left (12 / ‖c‖))
  · rw [Nat.cofinite_eq_atTop]
    filter_upwards [hev] with k hk
    have hkpos : (0 : ℝ) < ‖q‖ ^ k := by positivity
    have hbig : (2 : ℝ) ≤ ‖q ^ (-(k : ℤ)) * c‖ := by
      rw [hnorm k, le_div_iff₀ hkpos]
      linarith
    have h := norm_tateDterm_le_of_large hbig
    rw [hnorm k] at h
    refine h.trans (le_of_eq ?_)
    field_simp

/-- ★★★★★★**`∑_{n ∈ ℤ} h(q^{−n} u)` は可和**。 -/
theorem summable_int_tateDterm (q u : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) (hu0 : u ≠ 0) :
    Summable fun n : ℤ => tateDterm (q ^ (-n) * u) := by
  refine Summable.of_nat_of_neg ?_ ?_
  · exact summable_tateDterm_large q u hq hq0 hu0
  · refine (summable_tateDterm_small q u hq).congr fun k => ?_
    rw [neg_neg]

/-! ## ★★★★★★★`℘'` の q 展開 -/

/-- ★★★★★★★**`℘'` の q 展開**——`℘'(z) = (2πi)³·∑_{n∈ℤ} h(q^{−n}u)`。

★`℘` と違って補正項も定数項も出ない(重み 3 の和は絶対収束するので
`−1/l²` に当たるものが最初から無い)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem derivWeierstrassP_qExpansion (z τ : UpperHalfPlane)
    (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.derivWeierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 3 * ∑' n : ℤ,
        tateDterm (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z)) := by
  rw [derivWeierstrassP_eq_tsum_rows τ (z : ℂ)]
  have hpt : ∀ n : ℤ,
      (∑' m : ℤ, -2 / ((z : ℂ) - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 3)
        = (2 * ↑π * I) ^ 3 * tateDterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
            * Complex.exp (2 * ↑π * I * (z : ℂ))) := by
    intro n
    rw [rowDP_eq_smul τ (z : ℂ) n, rowDP z τ him n]
    ring
  rw [tsum_congr hpt, tsum_mul_left]

/-- ★★★★★★**`℘'` の q 展開(正べきの形)**——形式側と同じ向きに揃えた。 -/
theorem derivWeierstrassP_qExpansion_pos (z τ : UpperHalfPlane)
    (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.derivWeierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 3 * ∑' n : ℤ,
        tateDterm (Complex.exp (2 * ↑π * I * τ) ^ n * Complex.exp (2 * ↑π * I * z)) := by
  rw [derivWeierstrassP_qExpansion z τ him]
  congr 1
  rw [← (Equiv.neg ℤ).tsum_eq (fun n : ℤ =>
    tateDterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ n * Complex.exp (2 * ↑π * I * (z : ℂ))))]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def tateDterm_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘' 側の項の反転則)",
    sectionId := "genell-def-3-3" }

def lipschitz_three.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——重み 3 の Lipschitz 公式)",
    sectionId := "genell-def-3-3" }

def rowDP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘' 側の段の値の一様形)",
    sectionId := "genell-def-3-3" }

def derivWeierstrassP_qExpansion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘' の q 展開)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
