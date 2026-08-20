import ABC3.Found.GaloisRep.Infinity

/-!
# Galois (G5) 第 162 ブロック —— **★★★★★★★還元の核が部分群になる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★打ち消しの起きない証明が見つかった

還元 `red : E(F(W)) → E(F)` が群準同型であることの残りは
**「無限遠に落ちる点の全体 `E₁` が加法で閉じている」**であった。

素朴に `addX = ℓ² + a₁ℓ − a₂ − x₁ − x₂` の付値を評価しようとすると、
`w(ℓ²)` と `w(x₁)` が**ちょうど一致して打ち消す**ので進めない
(これは形式群を持ち出す通常の証明が避けている難所である)。

★**Vieta の第 2 基本対称式を使うと打ち消しが起きない**:

    x₁x₂ + x₁x₃ + x₂x₃ = −(2ℓm + a₁m + a₃ℓ − a₄),   m := y₁ − ℓx₁

* **左辺**: `w(x₁x₂) = w(x₁)·w(x₂)` が単独で最大(`w(x₃) ≤ 1` だから)。
* **右辺**: `m = y₃ − ℓx₃` と書き直せるので `w(m) ≤ max(1, w(ℓ))`、
  したがって `w(右辺) ≤ max(1, w(ℓ), w(ℓ)²) ≤ max(w x₁, w x₂)`。

★★`w(x₁)·w(x₂) ≤ max(w x₁, w x₂)` は `w(x₁), w(x₂) > 1` に矛盾する。

★★★**値群が ℤ であることも、平方根も、形式群も使っていない**——
乗法的な付値のまま完結する。

## ★★★★★★混合の場合は第 160 ブロックから**無料で出る**

`S₁ ∈ E₁`、`S₂ ∉ E₁` のときは、本ブロックの `val_addX_le_of_mixed` で
`S₃ := −(S₁+S₂)` も付値環に入ることが分かる。★そこで

    S₂ + S₃ = −S₁ ∈ E₁

に**第 160 の整な加法性**(`redPoint_add_some`)を当てると
`red S₂ + red S₃ = red(−S₁) = 0`、すなわち `red(S₁+S₂) = −red S₃ = red S₂`。
★★**混合の場合を直接評価する必要はない**(打ち消しを完全に回避できる)。

## ★★★Vieta は mathlib にあった

`addPolynomial_slope`(3 次式の**完全な因数分解**)と
`addPolynomial_eq`・`Cubic.prod_X_sub_C_eq`・`Cubic.toPoly_injective` を繋ぐと、
第 2 基本対称式が 4 行で出る。

## ★★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `val_y_le_of_val_x_le` | ★★★曲線上の点で `w x ≤ 1 ⟹ w y ≤ 1` |
| `one_lt_val_x_iff` | ★★★`w x > 1 ⟺ w y > 1` |
| `vieta_e2` | ★★★★★**第 2 基本対称式**(mathlib の因数分解から) |
| `val_slope_sq_le` | ★★★★`w(ℓ)² ≤ max(w x₁, w x₂)` |
| `one_lt_val_addX_of_infinity` | ★★★★★★★**`E₁` が加法で閉じている** |
| `val_slope_sq_eq` | ★★★★片方だけ無限遠なら `w(ℓ)² = w(x₁)` |
| `val_addX_le_of_mixed` | ★★★★★★**混合なら和は無限遠に落ちない** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★★Vieta の第 2 基本対称式(mathlib の因数分解から) -/

/-- ★★★★★**曲線上の 2 点と第 3 点の第 2 基本対称式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `addPolynomial_slope`(3 次式の完全な因数分解)と `addPolynomial_eq` を
`Cubic.toPoly_injective` で突き合わせるだけで出る。★★弦・接線を問わない。 -/
theorem vieta_e2 {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)
    {x₁ x₂ y₁ y₂ : F} (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂)
    (hxy : ¬(x₁ = x₂ ∧ y₁ = W.negY x₂ y₂)) :
    x₁ * x₂ + x₁ * (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂))
      + x₂ * (W.addX x₁ x₂ (W.slope x₁ x₂ y₁ y₂))
      = 2 * x₁ * (W.slope x₁ x₂ y₁ y₂) ^ 2
        + (W.a₁ * x₁ - 2 * y₁ - W.a₃) * (W.slope x₁ x₂ y₁ y₂) + (-W.a₁ * y₁ + W.a₄) := by
  have h := addPolynomial_slope h₁ h₂ hxy
  rw [addPolynomial_eq, Cubic.prod_X_sub_C_eq, neg_inj, Cubic.toPoly_injective] at h
  have hc := congrArg Cubic.c h
  simpa using hc.symm

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-! ## ★★★座標の付値は連動する -/

/-- ★★★曲線上の点では `w x ≤ 1 ⟹ w y ≤ 1`。

★`w y > 1` なら左辺の付値は `w(y)²` になり、右辺(`x` の多項式)の `≤ 1` に矛盾する。 -/
theorem val_y_le_of_val_x_le {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y)
    (hx : v.valuation W.FunctionField x ≤ 1) :
    v.valuation W.FunctionField y ≤ 1 := by
  set w := v.valuation W.FunctionField with hw
  by_contra hcon
  rw [not_le] at hcon
  rw [equation_iff] at heq
  simp only [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₃,
    WeierstrassCurve.map_a₄, WeierstrassCurve.map_a₆] at heq
  have hy2 : w y < w y ^ 2 := by
    have hp := pow_lt_pow_right₀ hcon (by norm_num : 1 < 2); rwa [pow_one] at hp
  have hlhs : w (y ^ 2 + algebraMap F W.FunctionField W.a₁ * x * y
      + algebraMap F W.FunctionField W.a₃ * y) = w y ^ 2 := by
    have hgroup : y ^ 2 + algebraMap F W.FunctionField W.a₁ * x * y
        + algebraMap F W.FunctionField W.a₃ * y
        = y ^ 2 + (algebraMap F W.FunctionField W.a₁ * x * y
          + algebraMap F W.FunctionField W.a₃ * y) := by ring
    have hrest : w (algebraMap F W.FunctionField W.a₁ * x * y
        + algebraMap F W.FunctionField W.a₃ * y) < w y ^ 2 := by
      refine lt_of_le_of_lt (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) hy2
      · rw [Valuation.map_mul, Valuation.map_mul]
        calc w (algebraMap F W.FunctionField W.a₁) * w x * w y ≤ 1 * 1 * w y :=
              mul_le_mul' (mul_le_mul' (valuation_algebraMap_field W v _) hx) le_rfl
          _ = w y := by rw [one_mul, one_mul]
      · rw [Valuation.map_mul]
        calc w (algebraMap F W.FunctionField W.a₃) * w y ≤ 1 * w y :=
              mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
          _ = w y := one_mul _
    rw [hgroup, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_pow]; exact ne_of_gt hrest), Valuation.map_pow,
      max_eq_left (le_of_lt hrest)]
  have hrhs : w (x ^ 3 + algebraMap F W.FunctionField W.a₂ * x ^ 2
      + algebraMap F W.FunctionField W.a₄ * x + algebraMap F W.FunctionField W.a₆) ≤ 1 :=
    val_add_le W v (val_add_le W v (val_add_le W v (val_pow_le W v hx 3)
      (val_mul_le W v (valuation_algebraMap_field W v _) (val_pow_le W v hx 2)))
      (val_mul_le W v (valuation_algebraMap_field W v _) hx))
      (valuation_algebraMap_field W v _)
  rw [heq] at hlhs
  rw [hlhs] at hrhs
  exact absurd hrhs (not_le.2 (lt_trans hcon hy2))

/-- ★★★曲線上の点では `w x > 1 ⟺ w y > 1`。

★(⟸) は `val_y_le_of_val_x_le` の対偶、(⟹) は第 161 の `w(y)² = w(x)³`。 -/
theorem one_lt_val_x_iff (h2 : IsUnit (2 : F)) {x y : W.FunctionField}
    (heq : (W.map (algebraMap F W.FunctionField)).Equation x y) :
    1 < v.valuation W.FunctionField x ↔ 1 < v.valuation W.FunctionField y := by
  constructor
  · intro hx
    have hkey := val_y_sq_eq_val_x_cube W v h2 heq hx
    by_contra hcon
    rw [not_lt] at hcon
    have h1 : v.valuation W.FunctionField y ^ 2 ≤ 1 := by
      calc v.valuation W.FunctionField y ^ 2 ≤ 1 ^ 2 := pow_le_pow_left' hcon 2
        _ = 1 := one_pow 2
    rw [hkey] at h1
    exact absurd h1 (not_le.2 (by
      have hp := pow_lt_pow_right₀ hx (by norm_num : 0 < 3)
      rwa [pow_zero] at hp))
  · intro hy
    by_contra hcon
    rw [not_lt] at hcon
    exact absurd (val_y_le_of_val_x_le W v heq hcon) (not_le.2 hy)

/-! ## ★★★★両方が無限遠に落ちる場合 -/

/-- ★★★★**傾きの 2 乗は `max(w x₁, w x₂)` を超えない**。

★`ℓ² = x₃ + x₁ + x₂ − a₁ℓ + a₂` の付値を読むだけ。
★★`w(ℓ) > max` と仮定すると `w(ℓ)² ≤ w(ℓ)` になって潰れる。 -/
theorem val_slope_sq_le {x₁ x₂ x₃ l : W.FunctionField}
    (hx₃ : x₃ = (W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l)
    (h1 : 1 < v.valuation W.FunctionField x₁)
    (h3 : v.valuation W.FunctionField x₃ ≤ 1) :
    v.valuation W.FunctionField l ^ 2
      ≤ max (v.valuation W.FunctionField x₁) (v.valuation W.FunctionField x₂) := by
  set w := v.valuation W.FunctionField with hw
  set X := max (w x₁) (w x₂) with hX
  have hX1 : 1 < X := lt_of_lt_of_le h1 (le_max_left _ _)
  have hdec : l ^ 2 = x₃ + x₁ + x₂
      - algebraMap F W.FunctionField W.a₁ * l + algebraMap F W.FunctionField W.a₂ := by
    rw [hx₃, WeierstrassCurve.Affine.addX, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂]; ring
  have hb : w (l ^ 2) ≤ max X (w l) := by
    rw [hdec]
    refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
    · refine le_trans (Valuation.map_sub w _ _) (max_le ?_ ?_)
      · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
        · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
          · exact le_trans (le_trans h3 (le_of_lt hX1)) (le_max_left _ _)
          · exact le_trans (le_max_left _ _) (le_max_left _ _)
        · exact le_trans (le_max_right _ _) (le_max_left _ _)
      · rw [Valuation.map_mul]
        exact le_trans (le_trans (mul_le_mul' (valuation_algebraMap_field W v _) le_rfl)
          (le_of_eq (one_mul _))) (le_max_right _ _)
    · exact le_trans (le_trans (valuation_algebraMap_field W v _) (le_of_lt hX1)) (le_max_left _ _)
  rw [Valuation.map_pow] at hb
  by_cases hl : w l ≤ X
  · rwa [max_eq_left hl] at hb
  · exfalso
    rw [not_le] at hl
    rw [max_eq_right (le_of_lt hl)] at hb
    have hlt : w l < w l ^ 2 := by
      have hp := pow_lt_pow_right₀ (lt_trans hX1 hl) (by norm_num : 1 < 2); rwa [pow_one] at hp
    exact absurd hb (not_le.2 hlt)

variable [DecidableEq F]

/-- ★★★★★★★**無限遠に落ちる 2 点の和も無限遠に落ちる**——`E₁` は加法で閉じている。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★Vieta の第 2 基本対称式で読むと、左辺は `w(x₁)·w(x₂)`(打ち消しなし)、
右辺は `max(w x₁, w x₂)` 以下になり、`w(x₁), w(x₂) > 1` に矛盾する。
★★形式群も値群の構造も使わない。 -/
theorem one_lt_val_addX_of_infinity {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Equation x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Equation x₂ y₂)
    (hxy : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂))
    (hv₁ : 1 < v.valuation W.FunctionField x₁) (hv₂ : 1 < v.valuation W.FunctionField x₂) :
    1 < v.valuation W.FunctionField ((W.map (algebraMap F W.FunctionField)).addX x₁ x₂
      ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂)) := by
  set w := v.valuation W.FunctionField with hw
  set l := (W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂ with hl
  set x₃ := (W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l with hx₃
  set X := max (w x₁) (w x₂) with hX
  have hX1 : 1 < X := lt_of_lt_of_le hv₁ (le_max_left _ _)
  by_contra hcon
  rw [not_lt] at hcon
  set y₃ := (W.map (algebraMap F W.FunctionField)).negAddY x₁ x₂ y₁ l with hy₃def
  have heq₃ : (W.map (algebraMap F W.FunctionField)).Equation x₃ y₃ := equation_negAdd h₁ h₂ hxy
  have hy₃ : w y₃ ≤ 1 := val_y_le_of_val_x_le W v heq₃ hcon
  set m := y₁ - l * x₁ with hm
  have hm₃ : m = y₃ - l * x₃ := by rw [hm, hy₃def, WeierstrassCurve.Affine.negAddY, ← hx₃]; ring
  have hmb : w m ≤ max 1 (w l) := by
    rw [hm₃]
    refine le_trans (Valuation.map_sub w _ _) (max_le (le_trans hy₃ (le_max_left _ _)) ?_)
    rw [Valuation.map_mul]
    exact le_trans (le_trans (mul_le_mul' le_rfl hcon) (le_of_eq (mul_one _))) (le_max_right _ _)
  have hl2 : w l ^ 2 ≤ X := val_slope_sq_le W v hx₃ hv₁ hcon
  have hlX : w l ≤ X := by
    by_cases hlo : w l ≤ 1
    · exact le_trans hlo (le_of_lt hX1)
    · rw [not_le] at hlo
      refine le_trans (le_of_lt ?_) hl2
      have hp := pow_lt_pow_right₀ hlo (by norm_num : 1 < 2); rwa [pow_one] at hp
  have hmX : w m ≤ X := le_trans hmb (max_le (le_of_lt hX1) hlX)
  have hlm : w l * w m ≤ X := by
    by_cases hlo : w l ≤ 1
    · calc w l * w m ≤ 1 * 1 := mul_le_mul' hlo (le_trans hmb (max_le le_rfl hlo))
        _ = 1 := one_mul 1
        _ ≤ X := le_of_lt hX1
    · rw [not_le] at hlo
      calc w l * w m ≤ w l * w l :=
            mul_le_mul' le_rfl (le_trans hmb (max_le (le_of_lt hlo) le_rfl))
        _ = w l ^ 2 := (sq (w l)).symm
        _ ≤ X := hl2
  have h2le : w (2 : W.FunctionField) ≤ 1 := by
    have he : (2 : W.FunctionField) = algebraMap F W.FunctionField (2 : F) := by rw [map_ofNat]
    rw [he]; exact valuation_algebraMap_field W v _
  have hvi := vieta_e2 (W.map (algebraMap F W.FunctionField)) h₁ h₂ hxy
  rw [← hl, ← hx₃] at hvi
  have hvi' : x₁ * x₂ + x₁ * x₃ + x₂ * x₃
      = -(2 * l * m + (W.map (algebraMap F W.FunctionField)).a₁ * m
        + (W.map (algebraMap F W.FunctionField)).a₃ * l
        - (W.map (algebraMap F W.FunctionField)).a₄) := by rw [hvi, hm]; ring
  have hRHS : w (x₁ * x₂ + x₁ * x₃ + x₂ * x₃) ≤ X := by
    rw [hvi', Valuation.map_neg]
    refine le_trans (Valuation.map_sub w _ _) (max_le ?_ ?_)
    · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
      · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
        · rw [Valuation.map_mul, Valuation.map_mul]
          calc w 2 * w l * w m ≤ 1 * (w l * w m) := by
                rw [mul_assoc]; exact mul_le_mul' h2le le_rfl
            _ = w l * w m := one_mul _
            _ ≤ X := hlm
        · rw [Valuation.map_mul]
          exact le_trans (mul_le_mul' (by
            rw [WeierstrassCurve.map_a₁]; exact valuation_algebraMap_field W v _) hmX)
            (le_of_eq (one_mul _))
      · rw [Valuation.map_mul]
        exact le_trans (mul_le_mul' (by
          rw [WeierstrassCurve.map_a₃]; exact valuation_algebraMap_field W v _) hlX)
          (le_of_eq (one_mul _))
    · rw [WeierstrassCurve.map_a₄]
      exact le_trans (valuation_algebraMap_field W v _) (le_of_lt hX1)
  have hgt : X < w x₁ * w x₂ :=
    max_lt (lt_mul_right (lt_trans zero_lt_one hv₁) hv₂)
      (by rw [mul_comm]; exact lt_mul_right (lt_trans zero_lt_one hv₂) hv₁)
  have hsmall : w (x₁ * x₃ + x₂ * x₃) < w (x₁ * x₂) := by
    rw [Valuation.map_mul]
    refine lt_of_le_of_lt (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) hgt
    · rw [Valuation.map_mul]
      exact le_trans (le_trans (mul_le_mul' le_rfl hcon) (le_of_eq (mul_one _)))
        (le_max_left _ _)
    · rw [Valuation.map_mul]
      exact le_trans (le_trans (mul_le_mul' le_rfl hcon) (le_of_eq (mul_one _)))
        (le_max_right _ _)
  rw [show x₁ * x₂ + x₁ * x₃ + x₂ * x₃ = x₁ * x₂ + (x₁ * x₃ + x₂ * x₃) from by ring,
    Valuation.map_add_of_distinct_val w (ne_of_gt hsmall),
    max_eq_left (le_of_lt hsmall), Valuation.map_mul] at hRHS
  exact absurd hRHS (not_le.2 hgt)

/-! ## ★★★★片方だけ無限遠に落ちる場合 -/

/-- ★★★★**片方だけ無限遠なら、傾きの 2 乗はちょうど `w(x₁)`**。

★`w(ℓ)·w(x₁) = w(y₁)` と第 161 の `w(y₁)² = w(x₁)³` から `w(ℓ)² = w(x₁)`。 -/
theorem val_slope_sq_eq (h2 : IsUnit (2 : F)) {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Equation x₁ y₁)
    (hv₁ : 1 < v.valuation W.FunctionField x₁)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    v.valuation W.FunctionField ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂) ^ 2
      = v.valuation W.FunctionField x₁ := by
  set w := v.valuation W.FunctionField with hw
  have hy₁ : 1 < w y₁ := (one_lt_val_x_iff W v h2 h₁).1 hv₁
  have hxne : x₁ ≠ x₂ := fun hh => absurd (hh ▸ hv₁) (not_lt.2 hx₂)
  have hsub : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hxne
  have hwx : w (x₁ - x₂) = w x₁ := by
    rw [sub_eq_add_neg, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg]; exact ne_of_gt (lt_of_le_of_lt hx₂ hv₁)),
      Valuation.map_neg, max_eq_left (le_of_lt (lt_of_le_of_lt hx₂ hv₁))]
  have hwy : w (y₁ - y₂) = w y₁ := by
    rw [sub_eq_add_neg, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg]; exact ne_of_gt (lt_of_le_of_lt hy₂ hy₁)),
      Valuation.map_neg, max_eq_left (le_of_lt (lt_of_le_of_lt hy₂ hy₁))]
  have hkey : w ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂) * w x₁ = w y₁ := by
    rw [slope_of_X_ne hxne, ← hwx, ← Valuation.map_mul, div_mul_cancel₀ _ hsub, hwy]
  have hcube := val_y_sq_eq_val_x_cube W v h2 h₁ hv₁
  have hx0 : w x₁ ≠ 0 := ne_of_gt (lt_trans zero_lt_one hv₁)
  have hsq : (w ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂)) ^ 2 * (w x₁) ^ 2
      = w x₁ * (w x₁) ^ 2 := by
    rw [← mul_pow, hkey, hcube, pow_succ, mul_comm]
  exact mul_right_cancel₀ (pow_ne_zero 2 hx0) hsq

/-- ★★★★★★**片方だけ無限遠なら、和は無限遠に落ちない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`w(ℓ)² = w(x₁)` と `m = y₂ − ℓx₂`(切片は `S₂` 側からも書ける)から
Vieta の右辺は `w(x₁)` 以下。★★`w(x₃) > 1` と仮定すると左辺は `w(x₁)·w(x₃)` になり矛盾。

★★★これと第 160 の整な加法性を合わせると、混合の場合の還元の加法性が
**直接評価なしで**従う(`S₂ + S₃ = −S₁`)。 -/
theorem val_addX_le_of_mixed (h2 : IsUnit (2 : F)) {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Equation x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Equation x₂ y₂)
    (hv₁ : 1 < v.valuation W.FunctionField x₁)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hxy : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂)) :
    v.valuation W.FunctionField ((W.map (algebraMap F W.FunctionField)).addX x₁ x₂
      ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂)) ≤ 1 := by
  set w := v.valuation W.FunctionField with hw
  set l := (W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂ with hl
  set x₃ := (W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l with hx₃
  have hxne : x₁ ≠ x₂ := fun hh => absurd (hh ▸ hv₁) (not_lt.2 hx₂)
  have hsub : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hxne
  have hls : w l ^ 2 = w x₁ := val_slope_sq_eq W v h2 h₁ hv₁ hx₂ hy₂
  have hl1 : 1 < w l := by
    by_contra hc
    rw [not_lt] at hc
    have hle : w l ^ 2 ≤ 1 := by
      calc w l ^ 2 ≤ 1 ^ 2 := pow_le_pow_left' hc 2
        _ = 1 := one_pow 2
    rw [hls] at hle
    exact absurd hle (not_le.2 hv₁)
  set m := y₁ - l * x₁ with hm
  have hmeq : m = y₂ - l * x₂ := by
    have hlx : l * (x₁ - x₂) = y₁ - y₂ := by
      rw [hl, slope_of_X_ne hxne, div_mul_cancel₀ _ hsub]
    rw [hm]; linear_combination -hlx
  have hmb : w m ≤ w l := by
    rw [hmeq]
    refine le_trans (Valuation.map_sub w _ _) (max_le (le_trans hy₂ (le_of_lt hl1)) ?_)
    rw [Valuation.map_mul]
    exact le_trans (mul_le_mul' le_rfl hx₂) (le_of_eq (mul_one _))
  have h2le : w (2 : W.FunctionField) ≤ 1 := by
    have he : (2 : W.FunctionField) = algebraMap F W.FunctionField (2 : F) := by rw [map_ofNat]
    rw [he]; exact valuation_algebraMap_field W v _
  have hvi := vieta_e2 (W.map (algebraMap F W.FunctionField)) h₁ h₂ hxy
  rw [← hl, ← hx₃] at hvi
  have hvi' : x₁ * x₂ + x₁ * x₃ + x₂ * x₃
      = -(2 * l * m + (W.map (algebraMap F W.FunctionField)).a₁ * m
        + (W.map (algebraMap F W.FunctionField)).a₃ * l
        - (W.map (algebraMap F W.FunctionField)).a₄) := by rw [hvi, hm]; ring
  have hll : w l ≤ w l ^ 2 := by
    have hp := pow_lt_pow_right₀ hl1 (by norm_num : 1 < 2)
    rw [pow_one] at hp; exact le_of_lt hp
  have hRHS : w (x₁ * x₂ + x₁ * x₃ + x₂ * x₃) ≤ w x₁ := by
    rw [hvi', Valuation.map_neg, ← hls]
    refine le_trans (Valuation.map_sub w _ _) (max_le ?_ ?_)
    · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
      · refine le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)
        · rw [Valuation.map_mul, Valuation.map_mul]
          calc w 2 * w l * w m ≤ 1 * (w l * w l) := by
                rw [mul_assoc]; exact mul_le_mul' h2le (mul_le_mul' le_rfl hmb)
            _ = w l ^ 2 := by rw [one_mul, ← sq]
        · rw [Valuation.map_mul]
          exact le_trans (mul_le_mul' (by
            rw [WeierstrassCurve.map_a₁]; exact valuation_algebraMap_field W v _)
            (le_trans hmb hll)) (le_of_eq (one_mul _))
      · rw [Valuation.map_mul]
        exact le_trans (mul_le_mul' (by
          rw [WeierstrassCurve.map_a₃]; exact valuation_algebraMap_field W v _) hll)
          (le_of_eq (one_mul _))
    · rw [WeierstrassCurve.map_a₄]
      exact le_trans (valuation_algebraMap_field W v _) (le_trans (le_of_lt hl1) hll)
  by_contra hcon
  rw [not_le] at hcon
  have hx₃0 : (0 : WithZero (Multiplicative ℤ)) < w x₃ := lt_trans zero_lt_one hcon
  have hx₁0 : (0 : WithZero (Multiplicative ℤ)) < w x₁ := lt_trans zero_lt_one hv₁
  have hA : w x₃ < w x₁ * w x₃ := by rw [mul_comm]; exact lt_mul_right hx₃0 hv₁
  have hB : w x₁ < w x₁ * w x₃ := lt_mul_right hx₁0 hcon
  have hne : w (x₂ * x₃) < w (x₁ * x₃) := by
    rw [Valuation.map_mul, Valuation.map_mul]
    exact lt_of_le_of_lt (mul_le_mul' hx₂ le_rfl) (by rw [one_mul]; exact hA)
  have hbig : w (x₁ * x₃ + x₂ * x₃) = w x₁ * w x₃ := by
    rw [Valuation.map_add_of_distinct_val w (ne_of_gt hne), max_eq_left (le_of_lt hne),
      Valuation.map_mul]
  have hsmall : w (x₁ * x₂) < w (x₁ * x₃ + x₂ * x₃) := by
    rw [hbig, Valuation.map_mul]
    exact lt_of_le_of_lt (le_trans (mul_le_mul' le_rfl hx₂) (le_of_eq (mul_one _))) hB
  rw [show x₁ * x₂ + x₁ * x₃ + x₂ * x₃ = x₁ * x₂ + (x₁ * x₃ + x₂ * x₃) from by ring,
    Valuation.map_add_of_distinct_val w (ne_of_lt hsmall), max_eq_right (le_of_lt hsmall),
    hbig] at hRHS
  exact absurd hRHS (not_le.2 hB)

/-! ## ★出典の紐付け(`.src`) -/

def one_lt_val_addX_of_infinity.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——還元の核が加法で閉じていること)",
    sectionId := "genell-thm-3-8" }

def val_addX_le_of_mixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——片方だけ無限遠なら和は無限遠に落ちないこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
