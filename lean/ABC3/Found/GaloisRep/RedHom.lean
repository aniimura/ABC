import ABC3.Found.GaloisRep.RedZero

/-!
# Galois (G5) 第 160 ブロック —— **★★★★★★★還元は加法を保つ(全場合)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★点の還元が加法を保つ

    redPoint (S₁ + S₂) = redPoint S₁ + redPoint S₂

★場合分けは 4 つで、いずれも既存のブロックで閉じる:

| 場合 | 出所 |
|---|---|
| `S₁ + S₂ = 0`(`x₁ = x₂` かつ `y₁ = negY x₂ y₂`) | 第 155 `redConst_negY` |
| `red D ≠ 0` | 第 158(弦・接線を問わない) |
| `red D = 0` かつ `red x₁ ≠ red x₂` | 第 156 |
| `red D = 0` かつ `red x₁ = red x₂` | **本ブロック**(両辺 0) |

## ★★★★★最後の場合——両辺が 0 になる

`red D = 0` かつ `red x₁ = red x₂` のとき、**傾きは必ず極を持つ**:

* `red y₁ ≠ red y₂` なら `w(y₁−y₂) = 1 > w(x₁−x₂)` で `w(slope) > 1`。
* `red y₁ = red y₂` なら還元先の点が **2 等分点**になるので、
  **非特異性**から分子 `N` の還元が 0 でない(`polynomialY = 0` なら `polynomialX ≠ 0`)。
  したがって `w(N) = 1 > w(D)` で `w(slope) > 1`。

★したがって `w(addX) = w(slope)² > 1` となり `redPoint = 0`(第 159)。
★★右辺も `red y₁ = negY(red x₂)(red y₂)` から 0(第 159)。

## ★★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `slopeNum_ne_zero_of_slopeDen_eq_zero` | ★★★★★★**非特異性から分子が消えない** |
| `val_lt_one_of_redConst_eq_zero` | ★★`red t = 0` なら `w t < 1` |
| `one_lt_val_slope_of_redD_zero` | ★★★★★★**`red D = 0` なら傾きは極を持つ** |
| `redPoint_add_some` | ★★★★★★★**還元は加法を保つ(全場合)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)

/-- ★★★★★★**非特異点で分母が消えるなら分子は消えない**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`polynomialY = 2y + a₁x + a₃` が 0 なら、非特異性から `polynomialX ≠ 0`。 -/
theorem slopeNum_ne_zero_of_slopeDen_eq_zero {x y : F} (hns : W.Nonsingular x y)
    (hD : y + y + W.a₁ * x + W.a₃ = 0) :
    x ^ 2 + x * x + x ^ 2 + W.a₂ * (x + x) + W.a₄ - W.a₁ * y ≠ 0 := by
  rcases hns.2 with hX | hY
  · rw [evalEval_polynomialX] at hX
    intro h0
    exact hX (by linear_combination -h0)
  · exfalso
    rw [evalEval_polynomialY] at hY
    exact hY (by linear_combination hD)

variable [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

theorem val_lt_one_of_redConst_eq_zero {t : W.FunctionField}
    (ht : v.valuation W.FunctionField t ≤ 1) (h0 : redConst W v h hv t = 0) :
    v.valuation W.FunctionField t < 1 :=
  lt_of_le_of_ne ht (fun heq => (redConst_ne_zero_of_eq_one W v h hv heq) h0)

variable [DecidableEq F] [W.IsElliptic]

/-- ★★★★★★**`red D = 0` のとき傾きは極を持つ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`red y₁ = red y₂` の側では**非特異性**が効く。 -/
theorem one_lt_val_slope_of_redD_zero {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂))
    (hxe : redConst W v h hv x₁ = redConst W v h hv x₂)
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) = 0) :
    1 < v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂) := by
  set w := v.valuation W.FunctionField with hw
  by_cases hye : redConst W v h hv y₁ = redConst W v h hv y₂
  · have hR2 : W.Nonsingular (redConst W v h hv x₂) (redConst W v h hv y₂) :=
      equation_iff_nonsingular.mp (equation_redConst W v h hv h₂.1 hx₂ hy₂)
    have hD' := hD
    rw [redConst_slopeDen W v h hv hx₁ hy₁ hy₂, hye, hxe] at hD'
    have hNred := slopeNum_ne_zero_of_slopeDen_eq_zero W hR2 hD'
    have hNv : w (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2
        + (W.map (algebraMap F W.FunctionField)).a₂ * (x₁ + x₂)
        + (W.map (algebraMap F W.FunctionField)).a₄
        - (W.map (algebraMap F W.FunctionField)).a₁ * y₂) = 1 := by
      refine val_eq_one_of_redConst_ne_zero W v h hv (val_slopeNum_le W v hx₁ hx₂ hy₂) ?_
      rw [redConst_slopeNum W v h hv hx₁ hx₂ hy₂, hxe]
      exact hNred
    have hNne : x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2
        + (W.map (algebraMap F W.FunctionField)).a₂ * (x₁ + x₂)
        + (W.map (algebraMap F W.FunctionField)).a₄
        - (W.map (algebraMap F W.FunctionField)).a₁ * y₂ ≠ 0 := by
      intro h0
      rw [h0, Valuation.map_zero] at hNv
      exact absurd hNv (by simp)
    have hDne : y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃ ≠ 0 := by
      intro h0
      have hid := slope_num_identity (W.map (algebraMap F W.FunctionField)) h₁.1 h₂.1
      rw [h0, mul_zero] at hid
      rcases mul_eq_zero.1 hid.symm with hz | hz
      · have hxeq : x₁ = x₂ := sub_eq_zero.1 hz
        have hyy : y₁ = y₂ := Y_eq_of_Y_ne h₁.1 h₂.1 hxeq (fun hh => hne ⟨hxeq, hh⟩)
        refine hne ⟨hxeq, ?_⟩
        rw [WeierstrassCurve.Affine.negY]
        subst hxeq; subst hyy
        linear_combination h0
      · exact hNne hz
    have hDlt : w (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) < 1 :=
      val_lt_one_of_redConst_eq_zero W v h hv (val_slopeDen_le W v hx₁ hy₁ hy₂) hD
    have hDpos : (0 : WithZero (Multiplicative ℤ)) <
        w (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
          + (W.map (algebraMap F W.FunctionField)).a₃) :=
      lt_of_le_of_ne zero_le (Ne.symm (by rw [Ne, Valuation.zero_iff]; exact hDne))
    rw [slope_eq_div _ h₁.1 h₂.1 hne hDne, div_eq_mul_inv, Valuation.map_mul, map_inv₀, hNv,
      one_mul, one_lt_inv₀ hDpos]
    exact hDlt
  · have hxne : x₁ ≠ x₂ := by
      intro hxeq
      exact hye (by rw [Y_eq_of_Y_ne h₁.1 h₂.1 hxeq (fun hh => hne ⟨hxeq, hh⟩)])
    have hsubx : w (x₁ - x₂) < 1 := by
      refine val_lt_one_of_redConst_eq_zero W v h hv (val_sub_le W v hx₁ hx₂) ?_
      rw [redConst_sub W v h hv hx₁ hx₂, hxe, sub_self]
    have hsuby : w (y₁ - y₂) = 1 := by
      refine val_eq_one_of_redConst_ne_zero W v h hv (val_sub_le W v hy₁ hy₂) ?_
      rw [redConst_sub W v h hv hy₁ hy₂]
      exact sub_ne_zero.2 hye
    have hxpos : (0 : WithZero (Multiplicative ℤ)) < w (x₁ - x₂) :=
      lt_of_le_of_ne zero_le (Ne.symm (by
        rw [Ne, Valuation.zero_iff]; exact sub_ne_zero.2 hxne))
    rw [slope_of_X_ne hxne, div_eq_mul_inv, Valuation.map_mul, map_inv₀, hsuby, one_mul,
      one_lt_inv₀ hxpos]
    exact hsubx

/-- ★★★★★★★**還元は加法を保つ**(全場合)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★4 つの場合に分かれ、いずれも第 155-159 で閉じる。 -/
theorem redPoint_add_some {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    redPoint W v h hv (Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂)
      = redPoint W v h hv (Point.some x₁ y₁ h₁) + redPoint W v h hv (Point.some x₂ y₂ h₂) := by
  by_cases hzero : x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂
  · rw [Point.add_of_Y_eq hzero.1 hzero.2, redPoint_zero,
      redPoint_some W v h hv h₁ hx₁ hy₁, redPoint_some W v h hv h₂ hx₂ hy₂]
    refine (Point.add_of_Y_eq (by rw [hzero.1]) ?_).symm
    rw [hzero.2, redConst_negY W v h hv hx₂ hy₂]
  · by_cases hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) = 0
    · by_cases hxe : redConst W v h hv x₁ = redConst W v h hv x₂
      · have hpole := one_lt_val_slope_of_redD_zero W v h hv h₁ h₂ hx₁ hy₁ hx₂ hy₂ hzero hxe hD
        rw [Point.add_some hzero, redPoint_eq_zero_of_not_le W v h hv _ (by
            intro hcon
            exact absurd hcon.1 (not_le.2 (one_lt_val_addX W v hx₁ hx₂ hpole))),
          redPoint_add_eq_zero_of_redD W v h hv h₁ h₂ hx₁ hy₁ hx₂ hy₂ hxe hD]
      · exact redPoint_add_of_redX_ne W v h hv h₁ h₂ hx₁ hy₁ hx₂ hy₂ hxe
    · exact redPoint_add_of_redD_ne_zero W v h hv h₁ h₂ hx₁ hy₁ hx₂ hy₂ hzero hD

/-! ## ★出典の紐付け(`.src`) -/

def redPoint_add_some.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——点の還元が加法を保つこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
