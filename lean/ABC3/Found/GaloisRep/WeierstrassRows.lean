import ABC3.Found.GaloisRep.LipschitzRows

/-!
# Galois (G6) 第 217 ブロック —— **★★★★★★℘ を段に切って各段の値を出す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★道 β の段 5(後半)

第 216 で各段の Lipschitz 公式を取った。本ブロックは **`℘` の二重和をその段に切る**。

| 段 | 内容 | 状態 |
|---|---|---|
| 4 | 格子 `ℤ + τℤ` の `g₂, g₃` を Eisenstein の q 展開に繋ぐ | 済(第 215) |
| 5 | `℘` の q 展開 | ★**本ブロックで各段の値まで** |
| 6 | 「`ℤ` 係数の形式級数が関数として 0 なら形式的に 0」 | 未着手 |

## ★★★★★段に切る

mathlib の `℘` は**格子の元にわたる**無条件可和な和である:

    ℘(z) = ∑' l : L.lattice, (1/(z−l)² − 1/l²)

`latticeEquivProd` で `ℤ × ℤ` に移し、`Summable.tsum_prod` で `n`(τ の係数)を外側にすると
**横一列ごとの和**になる。★段ごとには `1/(z−l)²` と `1/l²` が**別々に可和**なので
差に分けられる(`row_eq_sub`)——これは二重和全体では成り立たない
(`∑_l 1/l²` は絶対収束しない)。**段に切ってから分けるのが要点**である。

## ★★★★★★段の値は 1 つの形にまとまる

| 段 | ℘ 側 | 格子側 |
|---|---|---|
| `n ≤ 0` | `(2πi)²·f(q^{−n}u)`(第 216 の上向き) | `(2πi)²·f(q^{−n})` |
| `n ≥ 1` | `(2πi)²·f(q^{n}u⁻¹)`(第 216 の下向き) | `(2πi)²·f(qⁿ)` |
| `n = 0`(格子側) | —— | `∑_m 1/m² = 2ζ(2) = π²/3` |

★★★`tateXterm` は**反転で不変**(`tateXterm_inv`、既存)なので
`f(qⁿu⁻¹) = f(q^{−n}u)`——**上向きと下向きが 1 つの形 `f(q^{−n}u)` にまとまる**
(`rowP`)。`t ≠ 1` は `‖t‖ < 1`(`norm_down_lt_one`)から出る。
格子側も同様に `f(q^{−n})` にまとまる(`rowG`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `summable_one_div_sq_add` / `_sub` | 段ごとの可和性 |
| `hasSum_weierstrassP_prod` | ★★★★`℘` の和を `ℤ × ℤ` に移す |
| `weierstrassP_eq_tsum_rows` | ★★★★★**`℘` を段に切る** |
| `row_eq_sub` | ★★★★段ごとには差に分けられる |
| `rowP_nonpos` / `rowP_pos` / `rowP` | ★★★★★★**℘ 側の段の値** |
| `rowG_pos` / `rowG_neg` / `rowG_zero` / `rowG` | ★★★★★★**格子側の段の値** |
| `tsum_int_one_div_sq` | ★★★★★`∑_{m∈ℤ} 1/m² = π²/3` |
| `row_value` / `row_value_zero` | ★★★★★★**段の値(差の形)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-! ## ★段ごとの可和性 -/

/-- ★`∑_m 1/(w+m)²` は可和(重み 2 は 1 次元の段には十分)。 -/
theorem summable_one_div_sq_add (w : ℂ) : Summable fun m : ℤ => 1 / (w + (m : ℂ)) ^ 2 := by
  have h := EisensteinSeries.linear_right_summable w 1 (k := 2) (by norm_num)
  refine h.congr fun m => ?_
  push_cast
  rw [one_mul, one_div]
  rfl

/-- ★`∑_m 1/(w−m)²` は可和。 -/
theorem summable_one_div_sq_sub (w : ℂ) : Summable fun m : ℤ => 1 / (w - (m : ℂ)) ^ 2 := by
  rw [← (Equiv.neg ℤ).summable_iff]
  refine (summable_one_div_sq_add w).congr fun m => ?_
  simp only [Function.comp_apply, Equiv.neg_apply, Int.cast_neg]
  rw [show w - -(m : ℂ) = w + m by ring]

/-- ★`∑_m 1/(w−m)² = ∑_m 1/(w+m)²`(`m ↦ −m` の折り返し)。 -/
theorem tsum_one_div_sq_sub_eq_add (w : ℂ) :
    (∑' m : ℤ, 1 / (w - (m : ℂ)) ^ 2) = ∑' m : ℤ, 1 / (w + (m : ℂ)) ^ 2 := by
  rw [← (Equiv.neg ℤ).tsum_eq (fun m : ℤ => 1 / (w + (m : ℂ)) ^ 2)]
  refine tsum_congr fun m => ?_
  simp only [Equiv.neg_apply, Int.cast_neg]
  rw [show w + -(m : ℂ) = w - m by ring]

/-! ## ★★★★★℘ を段に切る -/

/-- ★周期対 `⟨τ, 1⟩` の格子の元は `n·τ + m` である。 -/
theorem tauPair_symm_apply (τ : UpperHalfPlane) (x : ℤ × ℤ) :
    (((tauPair τ).latticeEquivProd.symm x : (tauPair τ).lattice) : ℂ)
      = (x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ) := by
  simpa [tauPair] using (tauPair τ).latticeEquiv_symm_apply x

/-- ★★★★**`℘` の和を `ℤ × ℤ` の上に移す**。 -/
theorem hasSum_weierstrassP_prod (τ : UpperHalfPlane) (z : ℂ) :
    HasSum (fun x : ℤ × ℤ =>
      1 / (z - ((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ))) ^ 2
        - 1 / ((x.1 : ℂ) * (τ : ℂ) + (x.2 : ℂ)) ^ 2)
      (PeriodPair.weierstrassP (tauPair τ) z) := by
  have h := (tauPair τ).hasSum_weierstrassP z
  rw [← Equiv.hasSum_iff (tauPair τ).latticeEquivProd.toEquiv.symm] at h
  refine h.congr_fun fun x => ?_
  simp only [Function.comp_apply]
  have heq : (tauPair τ).latticeEquivProd.toEquiv.symm x
      = (tauPair τ).latticeEquivProd.symm x := rfl
  rw [heq, tauPair_symm_apply]

/-- ★★★★★**℘ を τ の段に切る**——二重和の内側が横一列である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem weierstrassP_eq_tsum_rows (τ : UpperHalfPlane) (z : ℂ) :
    PeriodPair.weierstrassP (tauPair τ) z
      = ∑' n : ℤ, ∑' m : ℤ,
          (1 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2
            - 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2) := by
  have h := hasSum_weierstrassP_prod τ z
  rw [← h.tsum_eq, h.summable.tsum_prod]

/-- ★★★★**各段は 2 つの和の差に分かれる**。

★二重和全体では成り立たない(`∑_l 1/l²` は絶対収束しない)。段に切ってから分けるのが要点。 -/
theorem row_eq_sub (τ : UpperHalfPlane) (z : ℂ) (n : ℤ) :
    (∑' m : ℤ, (1 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2
        - 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2))
      = (∑' m : ℤ, 1 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2)
        - ∑' m : ℤ, 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2 := by
  have h1 : Summable fun m : ℤ => 1 / (z - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2 := by
    refine (summable_one_div_sq_sub (z - (n : ℂ) * (τ : ℂ))).congr fun m => ?_
    rw [show z - (n : ℂ) * (τ : ℂ) - (m : ℂ) = z - ((n : ℂ) * (τ : ℂ) + (m : ℂ)) by ring]
  have h2 : Summable fun m : ℤ => 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2 :=
    summable_one_div_sq_add ((n : ℂ) * (τ : ℂ))
  exact h1.tsum_sub h2

/-! ## ★★★★★★℘ 側の段の値 -/

/-- ★★★★★**`n ≤ 0` の段の ℘ 側**——`n = −k` の段は `tateXterm (qᵏ u)` になる。 -/
theorem rowP_nonpos (z τ : UpperHalfPlane) (k : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - (((-(k : ℤ)) : ℤ) * (τ : ℂ) + (m : ℂ))) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm
        (Complex.exp (2 * ↑π * I * τ) ^ k * Complex.exp (2 * ↑π * I * z)) := by
  have hcong : ∀ m : ℤ,
      1 / ((z : ℂ) - (((-(k : ℤ)) : ℤ) * (τ : ℂ) + (m : ℂ))) ^ 2
        = 1 / (((z : ℂ) + (k : ℂ) * (τ : ℂ)) - (m : ℂ)) ^ 2 := by
    intro m
    push_cast
    ring_nf
  rw [tsum_congr hcong, tsum_one_div_sq_sub_eq_add]
  exact lipschitz_shift z τ k

/-- ★★★★★★**`n ≥ 1` の段の ℘ 側**——`n = j+1` の段は `tateXterm (q^{j+1} u⁻¹)` になる。 -/
theorem rowP_pos (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (j : ℕ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - (((j : ℤ) + 1) * (τ : ℂ) + (m : ℂ))) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm
        (Complex.exp (2 * ↑π * I * τ) ^ (j + 1) * (Complex.exp (2 * ↑π * I * z))⁻¹) := by
  have hcong : ∀ m : ℤ,
      1 / ((z : ℂ) - (((j : ℤ) + 1) * (τ : ℂ) + (m : ℂ))) ^ 2
        = 1 / (((z : ℂ) - ((j : ℕ) + 1 : ℕ) * (τ : ℂ)) - (m : ℂ)) ^ 2 := by
    intro m
    push_cast
    ring_nf
  rw [tsum_congr hcong, tsum_one_div_sq_sub_eq_add]
  exact lipschitz_shift_neg z τ him j

/-- ★下向きの段の `t = q^{j+1} u⁻¹` は `‖t‖ < 1`——特に `t ≠ 1`。 -/
theorem norm_down_lt_one (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (j : ℕ) :
    ‖Complex.exp (2 * ↑π * I * τ) ^ (j + 1) * (Complex.exp (2 * ↑π * I * z))⁻¹‖ < 1 := by
  have h := UpperHalfPlane.norm_exp_two_pi_I_lt_one (shiftDown z τ him j)
  rw [shiftDown_coe, exp_sub_nat_mul] at h
  exact h

/-- ★★★★★★**段の値の一様形**——`n` 段目は `(2πi)² · tateXterm (q^{−n} u)` である。

★上向き(`n ≤ 0`)は素直に、下向き(`n ≥ 1`)は `tateXterm` の反転不変性
(`tateXterm_inv`)で同じ形になる。`t ≠ 1` は `‖t‖ < 1` から出る。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem rowP (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) (n : ℤ) :
    (∑' m : ℤ, 1 / ((z : ℂ) - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm
        (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z)) := by
  have hq : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hu : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  rcases le_or_gt n 0 with hn0 | hn0
  · obtain ⟨k, rfl⟩ : ∃ k : ℕ, n = -(k : ℤ) := ⟨n.natAbs, by omega⟩
    rw [rowP_nonpos z τ k]
    congr 2
    rw [neg_neg, zpow_natCast]
  · obtain ⟨j, rfl⟩ : ∃ j : ℕ, n = (j : ℤ) + 1 := ⟨(n - 1).toNat, by omega⟩
    have h := rowP_pos z τ him j
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
    rw [h, ← tateXterm_inv hne0 hne1]
    congr 2
    rw [mul_inv, inv_inv, ← zpow_natCast (Complex.exp (2 * ↑π * I * τ)) (j + 1), ← zpow_neg]
    push_cast
    ring_nf

/-! ## ★★★★★★格子側の段の値 -/

/-- ★★★★★**`n ≥ 1` の段の格子側**——`∑_m 1/(nτ+m)² = (2πi)² tateXterm (qⁿ)`。 -/
theorem rowG_pos (τ : UpperHalfPlane) (j : ℕ) :
    (∑' m : ℤ, 1 / ((((j : ℤ) + 1) : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (j + 1)) := by
  have h := lipschitz_shift τ τ j
  have hcong : ∀ m : ℤ,
      1 / ((((j : ℤ) + 1) : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2
        = 1 / ((τ : ℂ) + (j : ℕ) * (τ : ℂ) + (m : ℂ)) ^ 2 := by
    intro m
    push_cast
    ring_nf
  rw [tsum_congr hcong, h]
  congr 2

/-- ★★★★**`n ≤ −1` の段の格子側は `n ≥ 1` と同じ**(折り返し)。 -/
theorem rowG_neg (τ : UpperHalfPlane) (j : ℕ) :
    (∑' m : ℤ, 1 / ((((-((j : ℤ) + 1))) : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (j + 1)) := by
  rw [← rowG_pos τ j, tsum_inv_sq_neg ((((j : ℤ) + 1) : ℂ) * (τ : ℂ))]
  refine tsum_congr fun m => ?_
  push_cast
  ring_nf

/-- ★`∑_{n:ℕ} 1/(n+1)² = ζ(2) = π²/6`。 -/
theorem tsum_nat_one_div_sq : (∑' n : ℕ, 1 / ((n : ℂ) + 1) ^ (2 : ℕ)) = (π : ℂ) ^ 2 / 6 := by
  have h := zeta_eq_tsum_one_div_nat_add_one_cpow (s := 2) (by norm_num)
  rw [riemannZeta_two] at h
  have h2 : (∑' n : ℕ, 1 / ((n : ℂ) + 1) ^ (2 : ℂ))
      = ∑' n : ℕ, 1 / ((n : ℂ) + 1) ^ (2 : ℕ) := by
    refine tsum_congr fun n => ?_
    rw [show ((2 : ℂ)) = ((2 : ℕ) : ℂ) by norm_num, Complex.cpow_natCast]
  rw [h2] at h
  exact h.symm

/-- ★`∑_{n:ℕ} 1/(n+1)²` は可和。 -/
theorem summable_nat_one_div_sq_shift :
    Summable fun n : ℕ => 1 / ((n : ℂ) + 1) ^ (2 : ℕ) := by
  have hinj : Function.Injective (fun n : ℕ => (n : ℤ) + 1) := by
    intro a b hab
    simpa using hab
  have h := (summable_one_div_sq_add (0 : ℂ)).comp_injective hinj
  refine h.congr fun n => ?_
  simp only [Function.comp_apply]
  push_cast
  ring_nf

/-- ★★★★★**`∑_{m ∈ ℤ} 1/m² = 2ζ(2) = π²/3`**(`m = 0` の項は `1/0 = 0`)。 -/
theorem tsum_int_one_div_sq : (∑' m : ℤ, 1 / ((m : ℂ)) ^ 2) = (π : ℂ) ^ 2 / 3 := by
  have h1 : Summable fun n : ℕ => 1 / ((((n : ℤ) + 1) : ℤ) : ℂ) ^ 2 := by
    refine summable_nat_one_div_sq_shift.congr fun n => ?_
    push_cast
    ring_nf
  have h2 : Summable fun n : ℕ => 1 / (((-((n : ℤ) + 1)) : ℤ) : ℂ) ^ 2 := by
    refine summable_nat_one_div_sq_shift.congr fun n => ?_
    push_cast
    rw [neg_pow]
    norm_num
  have key := tsum_of_add_one_of_neg_add_one (f := fun m : ℤ => 1 / ((m : ℂ)) ^ 2) h1 h2
  rw [key]
  have e1 : (∑' n : ℕ, 1 / ((((n : ℤ) + 1) : ℤ) : ℂ) ^ 2) = (π : ℂ) ^ 2 / 6 := by
    rw [← tsum_nat_one_div_sq]
    refine tsum_congr fun n => ?_
    push_cast
    ring_nf
  have e2 : (∑' n : ℕ, 1 / (((-((n : ℤ) + 1)) : ℤ) : ℂ) ^ 2) = (π : ℂ) ^ 2 / 6 := by
    rw [← tsum_nat_one_div_sq]
    refine tsum_congr fun n => ?_
    push_cast
    rw [neg_pow]
    norm_num
  simp only [e1, e2]
  norm_num
  ring

/-- ★★★★★**`n = 0` の段の格子側**——`∑_m 1/m² = π²/3`。 -/
theorem rowG_zero (τ : UpperHalfPlane) :
    (∑' m : ℤ, 1 / ((((0 : ℤ)) : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2) = (π : ℂ) ^ 2 / 3 := by
  rw [← tsum_int_one_div_sq]
  refine tsum_congr fun m => ?_
  push_cast
  ring_nf

/-- ★★★★★★**格子側の段の値の一様形**——`n ≠ 0` の段は `(2πi)² · tateXterm (q^{−n})`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem rowG (τ : UpperHalfPlane) {n : ℤ} (hn : n ≠ 0) :
    (∑' m : ℤ, 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2)
      = (2 * ↑π * I) ^ 2 * tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n)) := by
  have hq : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hqn : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  rcases lt_or_gt_of_ne hn with hn0 | hn0
  · obtain ⟨j, rfl⟩ : ∃ j : ℕ, n = -((j : ℤ) + 1) := ⟨(-n - 1).toNat, by omega⟩
    have h := rowG_neg τ j
    push_cast at h ⊢
    rw [h]
    congr 2
  · obtain ⟨j, rfl⟩ : ∃ j : ℕ, n = (j : ℤ) + 1 := ⟨(n - 1).toNat, by omega⟩
    have h := rowG_pos τ j
    push_cast at h ⊢
    have hpow1 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (j + 1) ≠ 1 := by
      intro hc
      have hlt : ‖Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (j + 1)‖ < 1 := by
        rw [norm_pow]
        exact pow_lt_one₀ (norm_nonneg _) hqn (Nat.succ_ne_zero j)
      rw [hc] at hlt
      simp at hlt
    have hpow0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (j + 1) ≠ 0 := pow_ne_zero _ hq
    rw [h, ← tateXterm_inv hpow0 hpow1]
    congr 2

/-! ## ★★★★★★段の値(差の形) -/

/-- ★★★★★★**`n ≠ 0` の段の値**——`(2πi)²·(f(q^{−n}u) − f(q^{−n}))`。 -/
theorem row_value (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) {n : ℤ} (hn : n ≠ 0) :
    (∑' m : ℤ, (1 / ((z : ℂ) - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2
        - 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2))
      = (2 * ↑π * I) ^ 2 *
        (tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z))
          - tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n))) := by
  rw [row_eq_sub τ (z : ℂ) n, rowP z τ him n, rowG τ hn]
  ring

/-- ★★★★★★**`n = 0` の段の値**——`(2πi)²·f(u) − π²/3`。 -/
theorem row_value_zero (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    (∑' m : ℤ, (1 / ((z : ℂ) - (((0 : ℤ) : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2
        - 1 / (((0 : ℤ) : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2))
      = (2 * ↑π * I) ^ 2 * tateXterm (Complex.exp (2 * ↑π * I * z)) - (π : ℂ) ^ 2 / 3 := by
  rw [row_eq_sub τ (z : ℂ) 0, rowP z τ him 0, rowG_zero τ]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def weierstrassP_eq_tsum_rows.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘ を τ の段に切る)",
    sectionId := "genell-def-3-3" }

def rowP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘ 側の段の値の一様形)",
    sectionId := "genell-def-3-3" }

def rowG.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子側の段の値の一様形)",
    sectionId := "genell-def-3-3" }

def row_value.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——段の値)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
