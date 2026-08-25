import ABC3.Found.GaloisRep.DerivWeierstrassQExp

/-!
# Galois (G6) 第 220 ブロック —— **★★★★★★★★解析側の Tate 方程式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★道 β の到達点

第 218・219 で `℘` と `℘'` の q 展開が取れた。本ブロックで**解析側の Tate 方程式**を出す。

    Y² + XY = X³ − 5s₃·X − (5s₃ + 7s₅)/12

★★これは形式側の `tateA4 = −5·s₃`、`12·tateA6 = −(5s₃ + 7s₅)` と**同じ形**である。
残るのは段 6(`ℤ` 係数の形式級数が関数として 0 なら形式的に 0)だけになった。

## ★★★★★変数変換

`℘'² = 4℘³ − g₂℘ − g₃` に

| 置き換え | 出どころ |
|---|---|
| `℘ = (2πi)²(X + 1/12)` | 第 218(定数項 `−π²/3 = (2πi)²/12`) |
| `℘' = (2πi)³(2Y + X)` | 第 219(`h = f + 2g` の分解) |
| `g₂ = (2πi)⁴(1 + 240s₃)/12` | 第 215 + `ζ(4) = π⁴/90`、`B₄ = −1/30` |
| `g₃ = −(2πi)⁶(1 − 504s₅)/216` | 第 215 + `ζ(6) = π⁶/945`、`B₆ = 1/42` |

を入れて `(2πi)⁶` で割ると、`X²` の項がちょうど相殺して Tate 方程式が残る。

★★★**`1/12` の平行移動と `2Y + X` の捻りが両方要る**——`℘` を平行移動するだけでは
`X²` の項が消えず、`℘'` を `2Y + X` と読み替えて初めて `Y² + XY` の形になる。

## ★★★mathlib に無かったもの

| 要るもの | 状況 |
|---|---|
| `ζ(4) = π⁴/90` | ★`riemannZeta_four` ✓ |
| `ζ(6) = π⁶/945` | ★無い——`riemannZeta_two_mul_nat` から作った |
| `B₄ = −1/30`、`B₆ = 1/42` | ★無い——`bernoulli'_def` を 1 段ずつ展開して作った |

★`decide` も `norm_num [bernoulli]` も通らない(`bernoulli'` の再帰が展開されない)。
`bernoulli'_def` を 1 回だけ `rw` して `Finset.sum_range_succ` と `Nat.choose` で潰す、を
`n = 3, 4, 5, 6` と積み上げるのが通る形だった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `bernoulli'_three` … `bernoulli_six_val` | ★Bernoulli 数の値 |
| `riemannZeta_six` | ★★`ζ(6) = π⁶/945` |
| `sigmaSum` | 解析側の `s_k(τ) = ∑_{n≥1} σ_k(n)qⁿ` |
| `g2_normalized` / `g3_normalized` | ★★★★★★`g₂, g₃` の正規化 |
| `tateXfun` / `tateYfun` | ★★★★★★解析側の `X, Y` |
| `weierstrassP_eq_tateXfun` | ★★★★★★`℘ = (2πi)²(X + 1/12)` |
| `derivWeierstrassP_eq_tateYfun` | ★★★★★★`℘' = (2πi)³(2Y + X)` |
| `tate_equation_analytic` | ★★★★★★★★**解析側の Tate 方程式** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-! ## ★Bernoulli 数と `ζ(6)` -/

theorem bernoulli'_three : (bernoulli' 3 : ℚ) = 0 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two, Nat.choose]

theorem bernoulli'_four : (bernoulli' 4 : ℚ) = -1/30 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two,
    bernoulli'_three, Nat.choose]

theorem bernoulli'_five : (bernoulli' 5 : ℚ) = 0 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two,
    bernoulli'_three, bernoulli'_four, Nat.choose]

theorem bernoulli'_six : (bernoulli' 6 : ℚ) = 1/42 := by
  rw [bernoulli'_def]
  norm_num [Finset.sum_range_succ, bernoulli'_zero, bernoulli'_one, bernoulli'_two,
    bernoulli'_three, bernoulli'_four, bernoulli'_five, Nat.choose]

theorem bernoulli_four_val : (bernoulli 4 : ℚ) = -1/30 := by
  rw [bernoulli_eq_bernoulli'_of_ne_one (by norm_num), bernoulli'_four]

theorem bernoulli_six_val : (bernoulli 6 : ℚ) = 1/42 := by
  rw [bernoulli_eq_bernoulli'_of_ne_one (by norm_num), bernoulli'_six]

/-- ★★`ζ(6) = π⁶/945`(mathlib には `riemannZeta_four` はあるが 6 は無かった)。 -/
theorem riemannZeta_six : riemannZeta 6 = (π : ℂ) ^ 6 / 945 := by
  have h := riemannZeta_two_mul_nat (k := 3) (by norm_num)
  norm_num [bernoulli_six_val] at h
  rw [show ((6 : ℂ)) = ((2 : ℕ) * (3 : ℕ) : ℕ) by norm_num]
  push_cast
  rw [h]
  ring

/-! ## ★★★★★★`g₂, g₃` の正規化 -/

/-- ★`s_k(τ) = ∑_{n≥1} σ_k(n) qⁿ`(解析側)。 -/
noncomputable def sigmaSum (k : ℕ) (τ : UpperHalfPlane) : ℂ :=
  ∑' n : ℕ+, ((ArithmeticFunction.sigma k (n : ℕ) : ℕ) : ℂ)
    * Complex.exp (2 * ↑π * I * τ) ^ ((n : ℕ) : ℤ)

theorem two_pi_I_pow_four : (2 * (π : ℂ) * I) ^ 4 = 16 * (π : ℂ) ^ 4 := by
  have hI : (I : ℂ) ^ 4 = 1 := by
    rw [show (4 : ℕ) = 2 * 2 by norm_num, pow_mul, Complex.I_sq]
    norm_num
  rw [mul_pow, mul_pow, hI]
  ring

theorem two_pi_I_pow_six : (2 * (π : ℂ) * I) ^ 6 = -(64 * (π : ℂ) ^ 6) := by
  have hI : (I : ℂ) ^ 6 = -1 := by
    rw [show (6 : ℕ) = 2 * 3 by norm_num, pow_mul, Complex.I_sq]
    norm_num
  rw [mul_pow, mul_pow, hI]
  ring

/-- ★★★★★★**`g₂` の正規化**——`g₂ = (2πi)⁴·(1 + 240 s₃)/12`。 -/
theorem g2_normalized (τ : UpperHalfPlane) :
    (tauPair τ).g₂ = (2 * ↑π * I) ^ 4 * ((1 + 240 * sigmaSum 3 τ) / 12) := by
  rw [g2_qExpansion τ, riemannZeta_four, two_pi_I_pow_four, sigmaSum]
  rw [show (((bernoulli 4 : ℚ)) : ℂ) = -1/30 by rw [bernoulli_four_val]; push_cast; ring]
  field_simp
  ring

/-- ★★★★★★**`g₃` の正規化**——`g₃ = −(2πi)⁶·(1 − 504 s₅)/216`。 -/
theorem g3_normalized (τ : UpperHalfPlane) :
    (tauPair τ).g₃ = (2 * ↑π * I) ^ 6 * (-((1 - 504 * sigmaSum 5 τ) / 216)) := by
  rw [g3_qExpansion τ, riemannZeta_six, two_pi_I_pow_six, sigmaSum]
  rw [show (((bernoulli 6 : ℚ)) : ℂ) = 1/42 by rw [bernoulli_six_val]; push_cast; ring]
  field_simp
  ring

/-! ## ★★★★★★解析側の `X, Y` -/

theorem tateYterm_eq_half (t : ℂ) : tateYterm t = (tateDterm t - tateXterm t) / 2 := by
  simp only [tateDterm]
  ring

/-- ★★★★★★**解析側の `X(u,q)`**。 -/
noncomputable def tateXfun (z τ : UpperHalfPlane) : ℂ :=
  (∑' n : ℤ, tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z)))
    - ∑' n : ℤ, tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n))

/-- ★★★★★★**解析側の `Y(u,q)`**——定数項は `s₁(q) = ∑_{n≥1} f(qⁿ)` にあたる。 -/
noncomputable def tateYfun (z τ : UpperHalfPlane) : ℂ :=
  (∑' n : ℤ, tateYterm (Complex.exp (2 * ↑π * I * τ) ^ (-n) * Complex.exp (2 * ↑π * I * z)))
    + (1 / 2) * ∑' n : ℤ, tateXterm (Complex.exp (2 * ↑π * I * τ) ^ (-n))

/-- ★`∑_n g(q^{−n}u)` は可和——`2g = h − f` から出る(新しい評価は要らない)。 -/
theorem summable_int_tateYterm (q u : ℂ) (hq : ‖q‖ < 1) (hq0 : q ≠ 0) (hu0 : u ≠ 0) :
    Summable fun n : ℤ => tateYterm (q ^ (-n) * u) := by
  have h := ((summable_int_tateDterm q u hq hq0 hu0).sub
    (summable_int_tateXterm q u hq hq0 hu0)).mul_right (1 / 2 : ℂ)
  refine h.congr fun n => ?_
  rw [tateYterm_eq_half]
  ring

/-- ★★★★★★**`℘ = (2πi)²(X + 1/12)`**。 -/
theorem weierstrassP_eq_tateXfun (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.weierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 2 * (tateXfun z τ + 1 / 12) := by
  have hq0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hu0 : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hqn : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  have hF := summable_int_tateXterm _ _ hqn hq0 hu0
  have hG : Summable fun n : ℤ =>
      tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)) := by
    refine (summable_int_tateXterm _ 1 hqn hq0 one_ne_zero).congr fun n => ?_
    rw [mul_one]
  rw [weierstrassP_qExpansion z τ him, tateXfun, ← hF.tsum_sub hG]
  have hpi : (2 * (π : ℂ) * I) ^ 2 = -(4 * (π : ℂ) ^ 2) := by
    rw [mul_pow, mul_pow, Complex.I_sq]
    ring
  rw [hpi]
  ring

/-- ★★★★★★**`℘' = (2πi)³(2Y + X)`**——`h = f + 2g` の分解がそのまま効く。 -/
theorem derivWeierstrassP_eq_tateYfun (z τ : UpperHalfPlane)
    (him : (z : ℂ).im < (τ : ℂ).im) :
    PeriodPair.derivWeierstrassP (tauPair τ) (z : ℂ)
      = (2 * ↑π * I) ^ 3 * (2 * tateYfun z τ + tateXfun z τ) := by
  have hq0 : Complex.exp (2 * ↑π * I * (τ : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hu0 : Complex.exp (2 * ↑π * I * (z : ℂ)) ≠ 0 := Complex.exp_ne_zero _
  have hqn : ‖Complex.exp (2 * ↑π * I * (τ : ℂ))‖ < 1 :=
    UpperHalfPlane.norm_exp_two_pi_I_lt_one τ
  have hF := summable_int_tateXterm _ _ hqn hq0 hu0
  have hY := summable_int_tateYterm _ _ hqn hq0 hu0
  rw [derivWeierstrassP_qExpansion z τ him, tateXfun, tateYfun]
  congr 1
  have hsplit : (∑' n : ℤ,
      tateDterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
        * Complex.exp (2 * ↑π * I * (z : ℂ))))
      = (∑' n : ℤ, tateXterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
          * Complex.exp (2 * ↑π * I * (z : ℂ))))
        + 2 * ∑' n : ℤ, tateYterm (Complex.exp (2 * ↑π * I * (τ : ℂ)) ^ (-n)
          * Complex.exp (2 * ↑π * I * (z : ℂ))) := by
    rw [← hY.tsum_mul_left, ← hF.tsum_add (hY.mul_left 2)]
    rfl
  rw [hsplit]
  ring

/-! ## ★★★★★★★★解析側の Tate 方程式 -/

theorem two_pi_I_ne_zero : (2 * (π : ℂ) * I) ≠ 0 := by
  have h1 : (π : ℂ) ≠ 0 := by
    exact_mod_cast Real.pi_ne_zero
  exact mul_ne_zero (mul_ne_zero two_ne_zero h1) Complex.I_ne_zero

/-- ★★★★★★★★**解析側の Tate 方程式**。

`℘'² = 4℘³ − g₂℘ − g₃` に `℘ = (2πi)²(X+1/12)`、`℘' = (2πi)³(2Y+X)`、
`g₂ = (2πi)⁴(1+240s₃)/12`、`g₃ = −(2πi)⁶(1−504s₅)/216` を入れて `(2πi)⁶` で割ると
`X²` の項が相殺して:

    Y² + XY = X³ − 5s₃·X − (5s₃ + 7s₅)/12

★★これは形式側の `tateA4 = −5s₃`、`12·tateA6 = −(5s₃ + 7s₅)` と同じ形である。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_equation_analytic (z τ : UpperHalfPlane) (him : (z : ℂ).im < (τ : ℂ).im)
    (hz : (z : ℂ) ∉ (tauPair τ).lattice) :
    tateYfun z τ ^ 2 + tateXfun z τ * tateYfun z τ
      = tateXfun z τ ^ 3 + (-5 * sigmaSum 3 τ) * tateXfun z τ
        - (5 * sigmaSum 3 τ + 7 * sigmaSum 5 τ) / 12 := by
  have h := (tauPair τ).derivWeierstrassP_sq (z : ℂ) hz
  rw [derivWeierstrassP_eq_tateYfun z τ him, weierstrassP_eq_tateXfun z τ him,
    g2_normalized τ, g3_normalized τ] at h
  have hne : (4 * (2 * (π : ℂ) * I) ^ 6 : ℂ) ≠ 0 :=
    mul_ne_zero (by norm_num) (pow_ne_zero _ two_pi_I_ne_zero)
  refine mul_left_cancel₀ hne ?_
  linear_combination h

/-! ## ★出典の紐付け(`.src`) -/

def g2_normalized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——g2 の正規化)",
    sectionId := "genell-def-3-3" }

def weierstrassP_eq_tateXfun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——(℘, ℘') から (X, Y) への変数変換)",
    sectionId := "genell-def-3-3" }

def tate_equation_analytic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——解析側の Tate 方程式)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
