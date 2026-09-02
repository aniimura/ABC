/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TorsionIntegralGood
import ABC3.Meta.Claim

/-!
# 第 1394 ブロック —— **核のノルムの因子の付値は 3 の倍数**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`v_p(N)` の帳簿を機械に載せる

`veluKernelNorm`（第 1386）の因子は `t_P = 2 y_P + a₁ x_P + a₃ = y_P − negY(x_P, y_P)`
である。★良い素点（`v_p(Δ) = 0`）で `v_p(2) = 0` のとき、**どの点でも `3 ∣ v_p(t_P)`**
が言える——これが第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` が要求する入力である。

☆場合は 2 つしかない:

| 場合 | `v_p(t_P)` | 道 |
|---|---|---|
| `v_p(x_P) < 0`（深い点） | `= v_p(y_P) = −3m` | 第 1070 `exists_depth_of_valAdd_x_neg` |
| `x_P` が整 | `= 0` | ★本ブロックの新しい段 |

## ★★★★★★★★整な点で `v_p(t_P) = 0` になる理由——**Bézout**

`t_P² = Ψ₂Sq(x_P)`（第 942 `psi2Sq_eval`）であり、倍化公式は

    x(2P) · Ψ₂Sq(x_P) = Φ₂(x_P)   （第 1049 `mulOK_two`）

である。★ここで **`Ψ₂Sq` と `Φ₂` の Bézout 係数が判別式を出す**:

    (12x³ − b₂x² − 10b₄x + b₂b₄ − 27b₆)·Ψ₂Sq(x)
      + (−48x² − 8b₂x + b₂² − 32b₄)·Φ₂(x) = Δ

☆（`Res(Ψ₂Sq, Φ₂) = Δ²` だが、`Δ` 自身が既にイデアルに入る。
係数は `tools` の sympy 計算で求め、`linear_combination … * W.b_relation` で検証した。）

★★したがって `v_p(Δ) = 0` なら **`Ψ₂Sq(x_P)` と `Φ₂(x_P)` が同時に深くなれない**。
`v_p(t_P) > 0` なら `Ψ₂Sq(x_P)` が深いので `Φ₂(x_P)` は単元、よって

    v_p(x(2P)) = 0 − 2 v_p(t_P) < 0

——`2P` が深い点になる。☆`2P` の座標が整と分かっていれば矛盾である。

★★★これで「`x_P` も `x(2P)` も整 ⟹ `v_p(t_P) = 0`」が出る。
`p ∤ l` の良い素点では第 1073（`mem_primeSubring_x_of_addOrderOf_prime`）が
核の全点の `x` を整にするので、`2P` の側も自動で満たされる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial IsDedekindDomain NumberField ABC3.Meta
open WeierstrassCurve.Affine

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★`valAdd` の積と冪（`Lˣ` を経由しない形） -/

/-- ☆**積の `valAdd`**——`Units.mk0` の中で掛け算する形。 -/
theorem valAdd_mulL (p : HeightOneSpectrum (𝓞 L)) {a b : L} (ha : a ≠ 0) (hb : b ≠ 0)
    (hab : a * b ≠ 0) :
    valAdd p (Units.mk0 (a * b) hab)
      = valAdd p (Units.mk0 a ha) + valAdd p (Units.mk0 b hb) := by
  rw [← valAdd_mul p (Units.mk0 a ha) (Units.mk0 b hb)]
  exact valAdd_eq_of_valuation_eq p _ _ (by simp)

/-- ☆**冪の `valAdd`**——`Units.mk0` の中で冪をとる形。 -/
theorem valAdd_powL (p : HeightOneSpectrum (𝓞 L)) {a : L} (ha : a ≠ 0) (n : ℕ)
    (han : a ^ n ≠ 0) :
    valAdd p (Units.mk0 (a ^ n) han) = (n : ℤ) * valAdd p (Units.mk0 a ha) := by
  rw [← valAdd_pow p (Units.mk0 a ha) n]
  exact valAdd_eq_of_valuation_eq p _ _ (by simp)

/-! ## ★★★★★★★★**Bézout**——`Ψ₂Sq` と `Φ₂` は判別式を生成する -/

/-- ☆`Ψ₂Sq` の値を展開した形。 -/
theorem psi2Sq_eval_expand {R : Type} [CommRing R] (W : WeierstrassCurve R) (x : R) :
    W.Ψ₂Sq.eval x = 4 * x ^ 3 + W.b₂ * x ^ 2 + 2 * W.b₄ * x + W.b₆ := by
  simp [WeierstrassCurve.Ψ₂Sq]

/-- ☆`Φ₂` の値を展開した形。 -/
theorem Phi2_eval_expand {R : Type} [CommRing R] (W : WeierstrassCurve R) (x : R) :
    (W.Φ 2).eval x = x ^ 4 - W.b₄ * x ^ 2 - 2 * W.b₆ * x - W.b₈ := by
  simp [WeierstrassCurve.Φ_two]

/-- ★★★★★★★★★★★★★★★★
**`Ψ₂Sq` と `Φ₂` の Bézout 恒等式**——★**無条件**（第 1394）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    (12x³ − b₂x² − 10b₄x + b₂b₄ − 27b₆)·Ψ₂Sq(x)
      + (−48x² − 8b₂x + b₂² − 32b₄)·Φ₂(x) = Δ

★★★これが「良い素点では `Ψ₂Sq(x)` と `Φ₂(x)` が同時に深くなれない」の根である。 -/
theorem bezout_psi2Sq_Phi2 {R : Type} [CommRing R] (W : WeierstrassCurve R) (x : R) :
    (12 * x ^ 3 - W.b₂ * x ^ 2 - 10 * W.b₄ * x + W.b₂ * W.b₄ - 27 * W.b₆)
        * (W.Ψ₂Sq.eval x)
      + (-(48 * x ^ 2) - 8 * W.b₂ * x + W.b₂ ^ 2 - 32 * W.b₄) * ((W.Φ 2).eval x)
      = W.Δ := by
  rw [psi2Sq_eval_expand, Phi2_eval_expand, WeierstrassCurve.Δ]
  linear_combination (12 * x ^ 2 + 2 * W.b₂ * x + 8 * W.b₄) * W.b_relation

/-! ## ★★★★整モデルでは `b` 不変量も整 -/

/-- ☆整モデルの `b` 不変量は付値環の元。 -/
theorem mem_primeSubring_b (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] :
    W.b₂ ∈ primeSubring p ∧ W.b₄ ∈ primeSubring p ∧ W.b₆ ∈ primeSubring p
      ∧ W.b₈ ∈ primeSubring p := by
  obtain ⟨h1, h2, h3, h4, h6⟩ := mem_primeSubring_of_isIntegral p W
  have h4' : (4 : L) ∈ primeSubring p := by simp
  have h2' : (2 : L) ∈ primeSubring p := by simp
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    simp only [WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
      WeierstrassCurve.b₈]
  · exact add_mem (pow_mem h1 2) (mul_mem h4' h2)
  · exact add_mem (mul_mem h2' h4) (mul_mem h1 h3)
  · exact add_mem (pow_mem h3 2) (mul_mem h4' h6)
  · exact sub_mem (add_mem (sub_mem (add_mem (mul_mem (pow_mem h1 2) h6)
      (mul_mem (mul_mem h4' h2) h6)) (mul_mem (mul_mem h1 h3) h4))
      (mul_mem h2 (pow_mem h3 2))) (pow_mem h4 2)

/-- ☆整な点では `t = y − negY(x,y) = 2y + a₁x + a₃` も整。 -/
theorem mem_primeSubring_negYdiff (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {x y : L} (hx : x ∈ primeSubring p) (hy : y ∈ primeSubring p) :
    (y - W.toAffine.negY x y) ∈ primeSubring p := by
  obtain ⟨h1, h2, h3, h4, h6⟩ := mem_primeSubring_of_isIntegral p W
  rw [WeierstrassCurve.Affine.negY]
  exact sub_mem hy (sub_mem (sub_mem (neg_mem hy) (mul_mem h1 hx)) h3)

/-! ## ★★★★★★★★★★★★良い素点では `Φ₂(x)` は単元 -/

/-- ★★★★★★★★★★★★
**`Ψ₂Sq(x)` が深いなら `Φ₂(x)` は単元**——★（第 1394）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆Bézout 恒等式で両方が深いと `Δ` も深くなり、`v_p(Δ) = 0` に反する。 -/
theorem valAdd_Phi2_eq_zero_of_psi_pos (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    {x : L} (hx : x ∈ primeSubring p)
    (hpsi : ValAtLeast p 1 (W.Ψ₂Sq.eval x)) :
    ∃ hne : (W.Φ 2).eval x ≠ 0, valAdd p (Units.mk0 ((W.Φ 2).eval x) hne) = 0 := by
  obtain ⟨hb2, hb4, hb6, hb8⟩ := mem_primeSubring_b p W
  have hc10 : (10 : L) ∈ primeSubring p := by simp
  have hc12 : (12 : L) ∈ primeSubring p := by simp
  have hc27 : (27 : L) ∈ primeSubring p := by simp
  have hc48 : (48 : L) ∈ primeSubring p := by simp
  have hc8 : (8 : L) ∈ primeSubring p := by simp
  have hc32 : (32 : L) ∈ primeSubring p := by simp
  have hc2 : (2 : L) ∈ primeSubring p := by simp
  have hAmem : (12 * x ^ 3 - W.b₂ * x ^ 2 - 10 * W.b₄ * x + W.b₂ * W.b₄ - 27 * W.b₆)
      ∈ primeSubring p :=
    sub_mem (add_mem (sub_mem (sub_mem (mul_mem hc12 (pow_mem hx 3))
      (mul_mem hb2 (pow_mem hx 2))) (mul_mem (mul_mem hc10 hb4) hx))
      (mul_mem hb2 hb4)) (mul_mem hc27 hb6)
  have hBmem : (-(48 * x ^ 2) - 8 * W.b₂ * x + W.b₂ ^ 2 - 32 * W.b₄) ∈ primeSubring p :=
    sub_mem (add_mem (sub_mem (neg_mem (mul_mem hc48 (pow_mem hx 2)))
      (mul_mem (mul_mem hc8 hb2) hx)) (pow_mem hb2 2)) (mul_mem hc32 hb4)
  have hPhiMem : (W.Φ 2).eval x ∈ primeSubring p := by
    rw [Phi2_eval_expand]
    exact sub_mem (sub_mem (sub_mem (pow_mem hx 4) (mul_mem hb4 (pow_mem hx 2)))
      (mul_mem (mul_mem hc2 hb6) hx)) hb8
  have hnot : ¬ ValAtLeast p 1 ((W.Φ 2).eval x) := by
    intro hcon
    have hsum : ValAtLeast p 1 W.Δ := by
      rw [← bezout_psi2Sq_Phi2 W x]
      refine valAtLeast_add ?_ ?_
      · exact valAtLeast_mono (by omega) (valAtLeast_mul (valAtLeast_of_mem hAmem) hpsi)
      · exact valAtLeast_mono (by omega) (valAtLeast_mul (valAtLeast_of_mem hBmem) hcon)
    have := hsum hΔ
    omega
  simp only [ValAtLeast, not_forall, not_le] at hnot
  obtain ⟨hne, hlt⟩ := hnot
  refine ⟨hne, ?_⟩
  have hge := valAtLeast_of_mem hPhiMem hne
  omega

/-! ## ★★★★★★★★★★★★★★★★倍化で深い点が出る -/

/-- ★★★★★★★★★★★★★★★★
**`v_p(t_P) > 0` なら `2P` は深い点**——★（第 1394）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`x(2P)·Ψ₂Sq(x) = Φ₂(x)`（倍化公式）で、右辺は単元・`Ψ₂Sq(x) = t²` なので

    v_p(x(2P)) = −2 v_p(t_P) < 0 -/
theorem valAdd_x_two_smul_eq (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) (hx : x ∈ primeSubring p)
    (hty : y ≠ W.toAffine.negY x y)
    (hpos : 0 < valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty))) :
    ∃ (x₂ y₂ : L) (h₂ : W.toAffine.Nonsingular x₂ y₂) (hx₂ : x₂ ≠ 0),
      (2 : ℕ) • (Point.some x y h) = Point.some x₂ y₂ h₂ ∧
      valAdd p (Units.mk0 x₂ hx₂)
        = -2 * valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty)) := by
  have htne : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
  have hsq : W.Ψ₂Sq.eval x = (y - W.toAffine.negY x y) ^ 2 := psi2Sq_eval W h.left
  have hpsine : W.Ψ₂Sq.eval x ≠ 0 := by rw [hsq]; exact pow_ne_zero 2 htne
  have hsqne : (y - W.toAffine.negY x y) ^ 2 ≠ 0 := pow_ne_zero 2 htne
  have hpsival : valAdd p (Units.mk0 (W.Ψ₂Sq.eval x) hpsine)
      = 2 * valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) := by
    rw [valAdd_congr p hpsine hsqne hsq, valAdd_powL p htne 2 hsqne]
    norm_num
  have hpsiVAL : ValAtLeast p 1 (W.Ψ₂Sq.eval x) := by
    intro hz
    have : valAdd p (Units.mk0 (W.Ψ₂Sq.eval x) hz)
        = valAdd p (Units.mk0 (W.Ψ₂Sq.eval x) hpsine) := valAdd_congr p hz hpsine rfl
    omega
  obtain ⟨hΦne, hΦ0⟩ := valAdd_Phi2_eq_zero_of_psi_pos p W hΔ hΔ0 hx hpsiVAL
  have h2ne : (W.ΨSq ((2 : ℕ) : ℤ)).eval x ≠ 0 := by
    rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two]; exact hpsine
  obtain ⟨x₂, y₂, h₂, hP2, hx2eq, -⟩ := mulOK_two W h h2ne
  rw [show ((2 : ℕ) : ℤ) = 2 from rfl, WeierstrassCurve.ΨSq_two] at hx2eq
  have hx₂ne : x₂ ≠ 0 := by
    intro hc
    apply hΦne
    rw [← hx2eq, hc, zero_mul]
  refine ⟨x₂, y₂, h₂, hx₂ne, hP2, ?_⟩
  have hprod : x₂ * W.Ψ₂Sq.eval x ≠ 0 := by rw [hx2eq]; exact hΦne
  have hval := valAdd_mulL p hx₂ne hpsine hprod
  rw [valAdd_congr p hprod hΦne hx2eq, hΦ0, hpsival] at hval
  omega

/-! ## ★★★★★★★★★★★★★★★★★★★★結論——`3 ∣ v_p(t_P)` -/

/-- ★★★★★★★★★★★★★★★★
**`P` も `2P` も整なら `v_p(t_P) = 0`**——★（第 1394）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v_p(t_P) > 0` なら `2P` が深い点になり、整性に反する。 -/
theorem valAdd_negYdiff_eq_zero_of_two_smul_integral
    (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y)
    (hx : x ∈ primeSubring p) (hy : y ∈ primeSubring p)
    (hty : y ≠ W.toAffine.negY x y)
    (hdbl : ∀ (x₂ y₂ : L) (h₂ : W.toAffine.Nonsingular x₂ y₂),
      (2 : ℕ) • (Point.some x y h) = Point.some x₂ y₂ h₂ → x₂ ∈ primeSubring p) :
    valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty)) = 0 := by
  have htne : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
  have hge : 0 ≤ valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) :=
    valAtLeast_of_mem (mem_primeSubring_negYdiff p W hx hy) htne
  by_contra hne
  have hpos : 0 < valAdd p (Units.mk0 (y - W.toAffine.negY x y) htne) := by omega
  obtain ⟨x₂, y₂, h₂, hx₂ne, hP2, hval⟩ :=
    valAdd_x_two_smul_eq p W hΔ hΔ0 h hx hty hpos
  have hmem := hdbl x₂ y₂ h₂ hP2
  have := valAtLeast_of_mem hmem hx₂ne
  omega

/-- ★★★★★★★★★★★★
**深い点では `v_p(t_P) = v_p(y_P) = −3m`**——★（第 1394）。

☆`v_p(2) = 0` なら `2y` が他の項より厳密に深いので和の付値を決める。 -/
theorem valAdd_negYdiff_eq_of_deep (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (hx : x ≠ 0) (hy : y ≠ 0) (heq : W.toAffine.Equation x y)
    (hneg : valAdd p (Units.mk0 x hx) < 0)
    (hty : y ≠ W.toAffine.negY x y) :
    valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty))
      = valAdd p (Units.mk0 y hy) := by
  obtain ⟨h1m, h2m, h3m, h4m, h6m⟩ := mem_primeSubring_of_isIntegral p W
  obtain ⟨m, hm0, hxm, hym⟩ := exists_depth_of_valAdd_x_neg p W hx hy heq hneg
  have htne : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
  have hsplit : y - W.toAffine.negY x y = 2 * y + (W.a₁ * x + W.a₃) := by
    rw [WeierstrassCurve.Affine.negY]; ring
  have h2yne : (2 : L) * y ≠ 0 := mul_ne_zero two_ne_zero hy
  have h2yval : valAdd p (Units.mk0 ((2 : L) * y) h2yne) = valAdd p (Units.mk0 y hy) := by
    rw [valAdd_mulL p two_ne_zero hy h2yne, h2]; ring
  have hrest : ValAtLeast p (-2 * (m : ℤ)) (W.a₁ * x + W.a₃) := by
    refine valAtLeast_add ?_ ?_
    · have hxA : ValAtLeast p (-2 * (m : ℤ)) x := by
        have := valAtLeast_unit p (Units.mk0 x hx)
        rw [hxm] at this; exact this
      exact valAtLeast_mono (by omega) (valAtLeast_mul (valAtLeast_of_mem h1m) hxA)
    · exact valAtLeast_mono (by omega) (valAtLeast_of_mem h3m)
  have hlt : valAdd p (Units.mk0 ((2 : L) * y) h2yne) < -2 * (m : ℤ) := by
    rw [h2yval, hym]; omega
  have hsumne : (2 : L) * y + (W.a₁ * x + W.a₃) ≠ 0 := by rw [← hsplit]; exact htne
  rw [valAdd_congr p htne hsumne hsplit, valAdd_add_eq_of_lt h2yne hsumne hrest hlt, h2yval]

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] 核のノルムの因子の付値は 3 の倍数**——★（第 1394）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆整な点では `0`、深い点では `−3m`——どちらでも `3 ∣ v_p(t_P)` である。

★★★これが第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` が要求する `3 ∣ v_p(N)` の
**点ごとの段**である。 -/
theorem three_dvd_valAdd_negYdiff (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [WeierstrassCurve.IsIntegral (primeSubring p) W] (hΔ : W.Δ ≠ 0)
    (hΔ0 : valAdd p (Units.mk0 W.Δ hΔ) = 0)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {x y : L} (h : W.toAffine.Nonsingular x y) (hty : y ≠ W.toAffine.negY x y)
    (hdbl : ∀ (x₂ y₂ : L) (h₂ : W.toAffine.Nonsingular x₂ y₂),
      (2 : ℕ) • (Point.some x y h) = Point.some x₂ y₂ h₂ → x₂ ∈ primeSubring p) :
    (3 : ℤ) ∣ valAdd p (Units.mk0 (y - W.toAffine.negY x y) (sub_ne_zero.mpr hty)) := by
  by_cases hx : x ∈ primeSubring p
  · have hy := mem_primeSubring_y_of_mem_x p W h.left hx
    rw [valAdd_negYdiff_eq_zero_of_two_smul_integral p W hΔ hΔ0 h hx hy hty hdbl]
    exact dvd_zero 3
  · have hx0 : x ≠ 0 := fun h0 => hx (by rw [h0]; exact zero_mem _)
    have hneg : valAdd p (Units.mk0 x hx0) < 0 := by
      by_contra hge
      rw [not_lt] at hge
      exact hx ((mem_primeSubring_iff p x).2 ((valAdd_nonneg_iff p _).1 hge))
    have hy0 : y ≠ 0 := y_ne_zero_of_valAdd_x_neg p W hx0 h.left hneg
    rw [valAdd_negYdiff_eq_of_deep p W h2 hx0 hy0 h.left hneg hty]
    obtain ⟨m, hm0, hxm, hym⟩ := exists_depth_of_valAdd_x_neg p W hx0 hy0 h.left hneg
    exact ⟨-(m : ℤ), by omega⟩

/-! ## ★出典の紐付け(`.src`) -/

def bezout_psi2Sq_Phi2.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ψ₂Sq と Φ₂ の Bézout 恒等式——係数は Δ を出す。★無条件)",
    sectionId := "genell-lemma-3-5" }

def bezout_psi2Sq_Phi2.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1394）**——係数は sympy で求めた" ++
       "（`Res(Ψ₂Sq, Φ₂) = Δ²` だが `Δ` 自身がイデアルに入る）。" ++
       "☆Lean 側は `linear_combination … * W.b_relation` で検証している。") 17 ]

def valAdd_Phi2_eq_zero_of_psi_pos.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ψ₂Sq(x) が深いなら Φ₂(x) は単元。★良い素点)",
    sectionId := "genell-lemma-3-5" }

def valAdd_x_two_smul_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(t_P) > 0 なら 2P は深い点。★良い素点)",
    sectionId := "genell-lemma-3-5" }

def valAdd_negYdiff_eq_zero_of_two_smul_integral.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(P も 2P も整なら v_p(t_P) = 0。★良い素点)",
    sectionId := "genell-lemma-3-5" }

def valAdd_negYdiff_eq_of_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い点では v_p(t_P) = v_p(y_P)。★v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_negYdiff.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核のノルムの因子の付値は 3 の倍数。★良い素点・v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_negYdiff.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_depth_of_valAdd_x_neg(第 1070、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_depth_of_valAdd_x_neg") 1,
    .citation "[ABC3]" "mulOK_two(第 1049、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.mulOK_two") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1394）**——`p ∣ l`（良い素点）の道の入力である。" ++
       "☆第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` は `3 ∣ v_p(N)` を要求し、" ++
       "`N = ∏ t_P` なので点ごとに `3 ∣ v_p(t_P)` を言えばよい。" ++
       "★整な点では `0`（Bézout ＋ 倍化）、深い点では `−3m`（第 1070）。") 17 ]

end ABC3.Found.GaloisRep
