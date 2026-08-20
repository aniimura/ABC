import ABC3.Found.GaloisRep.SlopeIdent

/-!
# Galois (G5) 第 158 ブロック —— **★★★★★★★還元は加法を保つ(分母が消えない場合)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★場合分けが 1 つに潰れた

第 157 の**統一表示**

    slope = N / D,    N = x₁²+x₁x₂+x₂²+a₂(x₁+x₂)+a₄−a₁y₂,   D = y₁+y₂+a₁x₁+a₃

により、**弦と接線を区別せずに**還元が扱える:

    red(D) ≠ 0   ⟹   redPoint (S₁ + S₂) = redPoint S₁ + redPoint S₂

★これは第 156(`red x₁ ≠ red x₂`)を**含む**——`red x₁ ≠ red x₂` でも
`red x₁ = red x₂` でも、`red D ≠ 0` でありさえすればよい。

## ★★★★機構

1. `red D ≠ 0` ⟹ `w(D) = 1`(第 155)。
2. `N` は多項式だから `w(N) ≤ 1`、したがって `w(slope) ≤ 1` かつ
   `red(slope) = red N / red D`(第 153 の `redConst_div`)。
3. 還元先でも `¬(red x₁ = red x₂ ∧ red y₁ = negY (red x₂) (red y₂))` が成り立つ
   (成り立たなければ `red D = 0` になる)ので、統一表示が使える。
4. `addX`・`addY` は第 155 で還元と可換。

## ★★★残るのは `red D = 0` の場合だけ

そこでは `red y₁ = negY (red x₂) (red y₂)` が直ちに出るので
**還元先で `red S₁ + red S₂ = 0`**。上でも傾きが極を持つので `redPoint = 0` になる。

## ★★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `redConst_slopeDen` / `val_slopeDen_le` | ★★分母の還元と付値 |
| `redConst_slopeNum` | ★★分子の還元 |
| `redD_ne_zero_imp` | ★★★★**`red D ≠ 0` なら還元先も退化しない** |
| `redConst_slope_of_redD_ne_zero` | ★★★★★★**傾きが還元と可換**(弦・接線を問わない) |
| `val_slope_le_of_redD_ne_zero` | ★★★傾きの付値は 1 以下 |
| `redPoint_add_of_redD_ne_zero` | ★★★★★★★**還元は加法を保つ** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

/-! ## ★★分母と分子の還元 -/

theorem redConst_slopeDen {x₁ y₁ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃)
      = redConst W v h hv y₁ + redConst W v h hv y₂ + W.a₁ * redConst W v h hv x₁ + W.a₃ := by
  have hA1 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁) ≤ 1 :=
    valuation_algebraMap_field W v _
  have hA3 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₃) ≤ 1 :=
    valuation_algebraMap_field W v _
  rw [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃,
    redConst_add W v h hv (val_add_le W v (val_add_le W v hy₁ hy₂)
      (val_mul_le W v hA1 hx₁)) hA3,
    redConst_add W v h hv (val_add_le W v hy₁ hy₂) (val_mul_le W v hA1 hx₁),
    redConst_add W v h hv hy₁ hy₂, redConst_mul W v h hv hA1 hx₁,
    redConst_algebraMap, redConst_algebraMap]

theorem val_slopeDen_le {x₁ y₁ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    v.valuation W.FunctionField (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
      + (W.map (algebraMap F W.FunctionField)).a₃) ≤ 1 := by
  rw [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃]
  exact val_add_le W v (val_add_le W v (val_add_le W v hy₁ hy₂)
    (val_mul_le W v (valuation_algebraMap_field W v _) hx₁))
    (valuation_algebraMap_field W v _)

theorem val_slopeNum_le {x₁ x₂ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    v.valuation W.FunctionField (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2
      + (W.map (algebraMap F W.FunctionField)).a₂ * (x₁ + x₂)
      + (W.map (algebraMap F W.FunctionField)).a₄
      - (W.map (algebraMap F W.FunctionField)).a₁ * y₂) ≤ 1 := by
  rw [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₄]
  exact val_sub_le W v (val_add_le W v (val_add_le W v (val_add_le W v
    (val_add_le W v (val_pow_le W v hx₁ 2) (val_mul_le W v hx₁ hx₂))
    (val_pow_le W v hx₂ 2))
    (val_mul_le W v (valuation_algebraMap_field W v _) (val_add_le W v hx₁ hx₂)))
    (valuation_algebraMap_field W v _))
    (val_mul_le W v (valuation_algebraMap_field W v _) hy₂)

theorem redConst_slopeNum {x₁ x₂ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1) :
    redConst W v h hv (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2
        + (W.map (algebraMap F W.FunctionField)).a₂ * (x₁ + x₂)
        + (W.map (algebraMap F W.FunctionField)).a₄
        - (W.map (algebraMap F W.FunctionField)).a₁ * y₂)
      = (redConst W v h hv x₁) ^ 2 + redConst W v h hv x₁ * redConst W v h hv x₂
        + (redConst W v h hv x₂) ^ 2
        + W.a₂ * (redConst W v h hv x₁ + redConst W v h hv x₂) + W.a₄
        - W.a₁ * redConst W v h hv y₂ := by
  have hA1 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁) ≤ 1 :=
    valuation_algebraMap_field W v _
  have hA2 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₂) ≤ 1 :=
    valuation_algebraMap_field W v _
  have hA4 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₄) ≤ 1 :=
    valuation_algebraMap_field W v _
  rw [WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂, WeierstrassCurve.map_a₄]
  have e1 : v.valuation W.FunctionField (x₁ ^ 2) ≤ 1 := val_pow_le W v hx₁ 2
  have e2 : v.valuation W.FunctionField (x₁ * x₂) ≤ 1 := val_mul_le W v hx₁ hx₂
  have e3 : v.valuation W.FunctionField (x₂ ^ 2) ≤ 1 := val_pow_le W v hx₂ 2
  have e4 : v.valuation W.FunctionField
      (algebraMap F W.FunctionField W.a₂ * (x₁ + x₂)) ≤ 1 :=
    val_mul_le W v hA2 (val_add_le W v hx₁ hx₂)
  have e5 : v.valuation W.FunctionField
      (algebraMap F W.FunctionField W.a₁ * y₂) ≤ 1 := val_mul_le W v hA1 hy₂
  rw [redConst_sub W v h hv (val_add_le W v (val_add_le W v (val_add_le W v
        (val_add_le W v e1 e2) e3) e4) hA4) e5,
    redConst_add W v h hv (val_add_le W v (val_add_le W v (val_add_le W v e1 e2) e3) e4) hA4,
    redConst_add W v h hv (val_add_le W v (val_add_le W v e1 e2) e3) e4,
    redConst_add W v h hv (val_add_le W v e1 e2) e3,
    redConst_add W v h hv e1 e2,
    redConst_pow W v h hv hx₁ 2, redConst_mul W v h hv hx₁ hx₂,
    redConst_pow W v h hv hx₂ 2,
    redConst_mul W v h hv hA2 (val_add_le W v hx₁ hx₂),
    redConst_add W v h hv hx₁ hx₂,
    redConst_mul W v h hv hA1 hy₂,
    redConst_algebraMap, redConst_algebraMap, redConst_algebraMap]

/-! ## ★★★★還元先も退化しない -/

/-- ★★★★**`red D ≠ 0` なら還元先の 2 点は互いに負元でない**。 -/
theorem redD_ne_zero_imp {x₁ x₂ y₁ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) ≠ 0) :
    ¬(redConst W v h hv x₁ = redConst W v h hv x₂ ∧
      redConst W v h hv y₁ = W.negY (redConst W v h hv x₂) (redConst W v h hv y₂)) := by
  rintro ⟨hxe, hye⟩
  refine hD ?_
  rw [redConst_slopeDen W v h hv hx₁ hy₁ hy₂, hye, hxe, WeierstrassCurve.Affine.negY]
  ring

/-! ## ★★★★★★傾きが還元と可換 -/

variable [DecidableEq F]

/-- ★★★★★★**`red D ≠ 0` なら傾きは還元と可換である**(弦・接線を問わない)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 157 の統一表示により、場合分けが要らなくなった。 -/
theorem redConst_slope_of_redD_ne_zero {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Equation x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Equation x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂))
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) ≠ 0) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂)
      = W.slope (redConst W v h hv x₁) (redConst W v h hv x₂)
        (redConst W v h hv y₁) (redConst W v h hv y₂) := by
  set D := y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
    + (W.map (algebraMap F W.FunctionField)).a₃ with hDdef
  have hDv : v.valuation W.FunctionField D ≤ 1 := val_slopeDen_le W v hx₁ hy₁ hy₂
  have hDone : v.valuation W.FunctionField D = 1 :=
    val_eq_one_of_redConst_ne_zero W v h hv hDv hD
  have hDne : D ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hDone
    exact absurd hDone (by simp)
  rw [slope_eq_div _ h₁ h₂ hne hDne,
    redConst_div W v h hv (val_slopeNum_le W v hx₁ hx₂ hy₂) hDone,
    redConst_slopeNum W v h hv hx₁ hx₂ hy₂, redConst_slopeDen W v h hv hx₁ hy₁ hy₂,
    slope_eq_div W (equation_redConst W v h hv h₁ hx₁ hy₁)
      (equation_redConst W v h hv h₂ hx₂ hy₂)
      (redD_ne_zero_imp W v h hv hx₁ hy₁ hy₂ hD)
      (by rw [← redConst_slopeDen W v h hv hx₁ hy₁ hy₂]; exact hD)]

theorem val_slope_le_of_redD_ne_zero {x₁ x₂ y₁ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Equation x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Equation x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂))
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) ≠ 0) :
    v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂) ≤ 1 := by
  set D := y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
    + (W.map (algebraMap F W.FunctionField)).a₃ with hDdef
  have hDv : v.valuation W.FunctionField D ≤ 1 := val_slopeDen_le W v hx₁ hy₁ hy₂
  have hDone : v.valuation W.FunctionField D = 1 :=
    val_eq_one_of_redConst_ne_zero W v h hv hDv hD
  have hDne : D ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hDone
    exact absurd hDone (by simp)
  rw [slope_eq_div _ h₁ h₂ hne hDne, div_eq_mul_inv, Valuation.map_mul, map_inv₀,
    hDone, inv_one, mul_one]
  exact val_slopeNum_le W v hx₁ hx₂ hy₂

/-- ★★★★★★★**還元は加法を保つ**(`red D ≠ 0` の場合)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★弦か接線かを問わない。第 156 の場合を**含む**。 -/
theorem redPoint_add_of_redD_ne_zero [W.IsElliptic] {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂))
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) ≠ 0) :
    redPoint W v h hv (Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂)
      = redPoint W v h hv (Point.some x₁ y₁ h₁) + redPoint W v h hv (Point.some x₂ y₂ h₂) := by
  set l := (W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂ with hl
  have hlv : v.valuation W.FunctionField l ≤ 1 :=
    val_slope_le_of_redD_ne_zero W v h hv h₁.1 h₂.1 hx₁ hx₂ hy₁ hy₂ hne hD
  rw [Point.add_some hne, redPoint_some W v h hv _ (val_addX_le W v hx₁ hx₂ hlv)
      (val_addY_le W v hx₁ hx₂ hy₁ hlv),
    redPoint_some W v h hv h₁ hx₁ hy₁, redPoint_some W v h hv h₂ hx₂ hy₂,
    Point.add_some (redD_ne_zero_imp W v h hv hx₁ hy₁ hy₂ hD)]
  refine point_some_congr ?_ ?_ _ _
  · rw [redConst_addX W v h hv hx₁ hx₂ hlv, hl,
      redConst_slope_of_redD_ne_zero W v h hv h₁.1 h₂.1 hx₁ hx₂ hy₁ hy₂ hne hD]
  · rw [redConst_addY W v h hv hx₁ hx₂ hy₁ hlv, hl,
      redConst_slope_of_redD_ne_zero W v h hv h₁.1 h₂.1 hx₁ hx₂ hy₁ hy₂ hne hD]

/-! ## ★出典の紐付け(`.src`) -/

def redPoint_add_of_redD_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——還元が加法を保つこと(分母が消えない場合))",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
