import ABC3.Found.GaloisRep.WeierstrassRows

/-!
# Galois (G6) 第 218 ブロック —— **★★★★★★★`℘` の q 展開(道 β の段 5 が閉じた)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★道 β の段 5 が閉じた

| 段 | 内容 | 状態 |
|---|---|---|
| 4 | 格子 `ℤ + τℤ` の `g₂, g₃` を Eisenstein の q 展開に繋ぐ | 済(第 215) |
| 5 | `℘` の q 展開 | ★★**本ブロックで完了** |
| 6 | 「`ℤ` 係数の形式級数が関数として 0 なら形式的に 0」 | 未着手 |

第 217 で各段の値が出た。残るのは**段の和が収束すること**である。

## ★★★★★両側とも幾何級数で押さえられる

`f(t) = t/(1−t)²` は `t → 0` でも `t → ∞` でも 0 に落ちる:

| 範囲 | 押さえ |
|---|---|
| `‖t‖ ≤ 1/2` | `‖1−t‖ ≥ 1/2` なので `‖f(t)‖ ≤ 4‖t‖` |
| `‖t‖ ≥ 2` | `‖1−t‖ ≥ ‖t‖/2` なので `‖f(t)‖ ≤ 4/‖t‖` |

★★`n → +∞` では `‖q^{−n}u‖ → ∞`、`n → −∞` では `‖q^{−n}u‖ → 0` なので、
**どちらの端も `4·‖q‖^{|n|}` 型の幾何級数**で押さえられる。
★★★大きい側で `tateXterm_inv`(反転不変性)を使う必要は無かった——
`‖f(t)‖ ≤ 4/‖t‖` を直接出せば `t ≠ 1` の仮定が要らない。

## ★★★★★★定数項 `1/12` の出どころ

`n ≠ 0` の段は `(2πi)²(f(q^{−n}u) − f(q^{−n}))` だが、`n = 0` の段だけ格子側が
`∑_m 1/m² = 2ζ(2) = π²/3` になる。一方 `f(q^0) = f(1) = 0`(Lean の
`Ring.inverse 0 = 0` の規約)なので、**一様形との差はちょうど `−π²/3` の 1 項だけ**である。

    ℘(z) = (2πi)²·∑_{n∈ℤ}(f(q^{−n}u) − f(q^{−n})) − π²/3

★`−π²/3 / (2πi)² = 1/12`——これが古典的な `℘/(2πi)² = X(u,q) + 1/12` の `1/12` である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXterm_field` / `norm_tateXterm` | ℂ の上では `f(t) = t/(1−t)²` |
| `norm_tateXterm_le_of_small` / `_large` | ★★★★両端の押さえ |
| `summable_tateXterm_small` / `_large` | ★★★★★片側ずつの可和性 |
| `summable_int_tateXterm` | ★★★★★★**`∑_{n∈ℤ} f(q^{−n}u)` は可和** |
| `weierstrassP_qExpansion` | ★★★★★★★**`℘` の q 展開** |
| `weierstrassP_qExpansion_pos` | ★★★★★★同じものを正べきで(形式側と同じ向き) |
-/

namespace ABC3.Found.GaloisRep

open Complex Real Filter PeriodPair

/-! ## ★`f(t) = t/(1−t)²` の両端の押さえ -/

/-- ★ℂ の上では `tateXterm t = t/(1−t)²`。 -/
theorem tateXterm_field (t : ℂ) : tateXterm t = t / (1 - t) ^ 2 := by
  simp [tateXterm, Ring.inverse_eq_inv', div_eq_mul_inv, inv_pow]

/-- ★`‖f(t)‖ = ‖t‖/‖1−t‖²`。 -/
theorem norm_tateXterm (t : ℂ) : ‖tateXterm t‖ = ‖t‖ / ‖1 - t‖ ^ 2 := by
  rw [tateXterm_field, norm_div, norm_pow]

/-- ★★★★`‖t‖ ≤ 1/2` なら `‖f(t)‖ ≤ 4‖t‖`。 -/
theorem norm_tateXterm_le_of_small {t : ℂ} (ht : ‖t‖ ≤ 1 / 2) :
    ‖tateXterm t‖ ≤ 4 * ‖t‖ := by
  have h1 : (1 : ℝ) / 2 ≤ ‖1 - t‖ := by
    have := norm_sub_norm_le (1 : ℂ) t
    simp only [norm_one] at this
    linarith
  have h2 : (0 : ℝ) < ‖1 - t‖ := by linarith
  rw [norm_tateXterm, div_le_iff₀ (by positivity)]
  nlinarith [norm_nonneg t, sq_nonneg (‖1 - t‖ - 1/2)]

/-- ★★★★`2 ≤ ‖t‖` なら `‖f(t)‖ ≤ 4/‖t‖`。

★これを直接出せば、大きい側で反転不変性(`tateXterm_inv`)を使わずに済み、
`t ≠ 1` の仮定が要らない。 -/
theorem norm_tateXterm_le_of_large {t : ℂ} (ht : 2 ≤ ‖t‖) :
    ‖tateXterm t‖ ≤ 4 / ‖t‖ := by
  have h1 : ‖t‖ / 2 ≤ ‖1 - t‖ := by
    have := norm_sub_norm_le t (1 : ℂ)
    simp only [norm_one] at this
    rw [show (1 : ℂ) - t = -(t - 1) by ring, norm_neg]
    linarith
  have h2 : (0 : ℝ) < ‖t‖ := by linarith
  have h3 : (0 : ℝ) < ‖1 - t‖ := by linarith
  rw [norm_tateXterm, div_le_div_iff₀ (by positivity) h2]
  nlinarith [sq_nonneg (‖1 - t‖ - ‖t‖ / 2)]

/-! ## ★★★★★段の和の可和性 -/

/-- ★★★★★`‖q‖ < 1` なら `∑_k f(qᵏ c)` は可和(小さい側)。 -/
theorem summable_tateXterm_small (q c : ℂ) (hq : ‖q‖ < 1) :
    Summable fun k : ℕ => tateXterm (q ^ ((k : ℤ)) * c) := by
  have hten : Tendsto (fun k : ℕ => ‖q‖ ^ k * ‖c‖) atTop (nhds 0) := by
    have := tendsto_pow_atTop_nhds_zero_of_lt_one (norm_nonneg q) hq
    simpa using this.mul_const ‖c‖
  have hev : ∀ᶠ k : ℕ in atTop, ‖q‖ ^ k * ‖c‖ ≤ 1 / 2 :=
    hten.eventually_le_const (by norm_num)
  refine Summable.of_norm_bounded_eventually
    (g := fun k : ℕ => 4 * ‖c‖ * ‖q‖ ^ k) ?_ ?_
  · exact ((summable_geometric_of_lt_one (norm_nonneg q) hq).mul_left (4 * ‖c‖))
  · rw [Nat.cofinite_eq_atTop]
    filter_upwards [hev] with k hk
    have hn : ‖q ^ ((k : ℤ)) * c‖ = ‖q‖ ^ k * ‖c‖ := by
      rw [norm_mul, zpow_natCast, norm_pow]
    have := norm_tateXterm_le_of_small (t := q ^ ((k : ℤ)) * c) (by rw [hn]; exact hk)
    rw [hn] at this
    linarith [this]

/-- ★★★★★`‖q‖ < 1`、`q ≠ 0`、`c ≠ 0` なら `∑_k f(q^{−k} c)` は可和(大きい側)。 -/
theorem summable_tateXterm_large (q c : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) (hc0 : c ≠ 0) :
    Summable fun k : ℕ => tateXterm (q ^ (-(k : ℤ)) * c) := by
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
    (g := fun k : ℕ => 4 / ‖c‖ * ‖q‖ ^ k) ?_ ?_
  · exact ((summable_geometric_of_lt_one (norm_nonneg q) hq).mul_left (4 / ‖c‖))
  · rw [Nat.cofinite_eq_atTop]
    filter_upwards [hev] with k hk
    have hkpos : (0 : ℝ) < ‖q‖ ^ k := by positivity
    have hbig : (2 : ℝ) ≤ ‖q ^ (-(k : ℤ)) * c‖ := by
      rw [hnorm k, le_div_iff₀ hkpos]
      linarith
    have h := norm_tateXterm_le_of_large hbig
    rw [hnorm k] at h
    refine h.trans (le_of_eq ?_)
    field_simp

/-- ★★★★★★**`∑_{n ∈ ℤ} f(q^{−n} u)` は可和**——両側とも幾何級数で押さえられる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem summable_int_tateXterm (q u : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) (hu0 : u ≠ 0) :
    Summable fun n : ℤ => tateXterm (q ^ (-n) * u) := by
  refine Summable.of_nat_of_neg ?_ ?_
  · exact summable_tateXterm_large q u hq hq0 hu0
  · refine (summable_tateXterm_small q u hq).congr fun k => ?_
    rw [neg_neg]

/-! ## ★★★★★★★`℘` の q 展開 -/

/-- ★`f(1) = 0`(Lean の `Ring.inverse 0 = 0` の規約)。 -/
theorem tateXterm_one_complex : tateXterm (1 : ℂ) = 0 := by simp [tateXterm]

/-- ★★★★★★★**`℘` の q 展開**。

    ℘(z) = (2πi)²·∑_{n∈ℤ}(f(q^{−n}u) − f(q^{−n})) − π²/3

★`n = 0` の段だけ格子側が `π²/3 = 2ζ(2)` になり、それが定数項を生む。
★★`f(q^0) = f(1) = 0` なので一様形との差はちょうど 1 項だけである。
★★★`−π²/3 / (2πi)² = 1/12`——古典的な `℘/(2πi)² = X(u,q) + 1/12` の `1/12`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem weierstrassP_qExpansion (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.weierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 2 * (∑' n : ℤ,
          (tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z))
            - tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n))))
        - (π : ℂ) ^ 2 / 3 := by
  have hq0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hu0 : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hqn : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  have hF := summable_int_tateXterm _ _ hqn hq0 hu0
  have hG : Summable fun n : ℤ =>
      tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)) := by
    refine (summable_int_tateXterm _ 1 hqn hq0 one_ne_zero).congr fun n => ?_
    rw [mul_one]
  have hA : Summable fun n : ℤ =>
      (2 * ↑π * I) ^ 2 * (tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
        * Complex.exp (2 * ↑π * I * (z : ℂ)))
        - tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n))) :=
    (hF.sub hG).mul_left _
  have hB : Summable fun n : ℤ => (if n = 0 then -((π : ℂ) ^ 2 / 3) else 0) :=
    (hasSum_ite_eq 0 (-((π : ℂ) ^ 2 / 3))).summable
  have hpt : ∀ n : ℤ,
      (∑' m : ℤ, (1 / ((z : ℂ) - ((n : ℂ) * (τ : ℂ) + (m : ℂ))) ^ 2
          - 1 / ((n : ℂ) * (τ : ℂ) + (m : ℂ)) ^ 2))
        = (2 * ↑π * I) ^ 2 * (tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
            * Complex.exp (2 * ↑π * I * (z : ℂ)))
            - tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)))
          + (if n = 0 then -((π : ℂ) ^ 2 / 3) else 0) := by
    intro n
    rcases eq_or_ne n 0 with rfl | hn
    · rw [row_value_zero z τ him, if_pos rfl, neg_zero, zpow_zero, one_mul,
        tateXterm_one_complex]
      ring
    · rw [row_value z τ him hn, if_neg hn, add_zero]
  rw [weierstrassP_eq_tsum_rows τ (z : ℂ), tsum_congr hpt, hA.tsum_add hB,
    (hasSum_ite_eq 0 (-((π : ℂ) ^ 2 / 3))).tsum_eq, tsum_mul_left]
  ring

/-- ★★★★★★**`℘` の q 展開(正べきの形)**——形式側の `tateXpair` と同じ向きに揃えた。

`∑_{n∈ℤ}` は `n ↦ −n` で不変なので、`q^{−n}` を `qⁿ` に書き換えられる。 -/
theorem weierstrassP_qExpansion_pos (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.weierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 2 * (∑' n : ℤ,
          (tateXterm (Complex.exp (2 * ↑π * I * τ) ^ n * Complex.exp (2 * ↑π * I * z))
            - tateXterm (Complex.exp (2 * ↑π * I * τ) ^ n)))
        - (π : ℂ) ^ 2 / 3 := by
  rw [weierstrassP_qExpansion z τ him]
  congr 2
  rw [← (Equiv.neg ℤ).tsum_eq (fun n : ℤ =>
    tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ n * Complex.exp (2 * ↑π * I * (z : ℂ)))
      - tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ n))]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def summable_int_tateXterm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——段の和の可和性)",
    sectionId := "genell-def-3-3" }

def weierstrassP_qExpansion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘ の q 展開)",
    sectionId := "genell-def-3-3" }

def weierstrassP_qExpansion_pos.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——℘ の q 展開、正べきの形)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
