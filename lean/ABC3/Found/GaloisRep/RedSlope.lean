import ABC3.Found.GaloisRep.RedFormula

/-!
# Galois (G5) 第 156 ブロック —— **★★★★★★還元は加法を保つ(弦の場合)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★`red x₁ ≠ red x₂` の場合が閉じた

    redPoint (S₁ + S₂) = redPoint S₁ + redPoint S₂        (red x₁ ≠ red x₂)

★機構は 3 段である:

1. `red x₁ ≠ red x₂` ⟹ `red(x₁ − x₂) ≠ 0` ⟹ **`w(x₁ − x₂) = 1`**(第 155)。
2. したがって `slope = (y₁ − y₂)/(x₁ − x₂)` は付値が 1 以下で、
   第 153 の `redConst_div` から **`red(slope) = slope(red …)`**。
3. `addX`・`addY` は第 155 で還元と可換だから、両辺が同じ点になる。

★★**還元先でも弦の公式が使える**(`red x₁ ≠ red x₂` だから)のが要点である。

## ★★★残るのは `red x₁ = red x₂` の場合

そこでは還元先で **2 倍公式**に切り替わる。★`x₁ ≠ x₂` なのに `red x₁ = red x₂` の
場合(`w(x₁ − x₂) < 1`)が本質的に難しい——傾きが `O_v` に入らないことがある。

## ★★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `redConst_slope_of_redX_ne` | ★★★★★**傾きが還元と可換** |
| `val_slope_le_of_redX_ne` | ★★★傾きの付値は 1 以下 |
| `redPoint_add_of_redX_ne` | ★★★★★★**還元は加法を保つ(弦の場合)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)
  {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

/-- ★★★★★**`red x₁ ≠ red x₂` のとき傾きは還元と可換である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`red(x₁ − x₂) ≠ 0` から `w(x₁ − x₂) = 1` が出て、`redConst_div` が効く。 -/
theorem redConst_slope_of_redX_ne {x₁ x₂ y₁ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : redConst W v h hv x₁ ≠ redConst W v h hv x₂) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂)
      = W.slope (redConst W v h hv x₁) (redConst W v h hv x₂)
        (redConst W v h hv y₁) (redConst W v h hv y₂) := by
  have hxne : x₁ ≠ x₂ := fun hxeq => hne (by rw [hxeq])
  have hsub : v.valuation W.FunctionField (x₁ - x₂) ≤ 1 := val_sub_le W v hx₁ hx₂
  have hredsub : redConst W v h hv (x₁ - x₂) ≠ 0 := by
    rw [redConst_sub W v h hv hx₁ hx₂]
    exact sub_ne_zero.2 hne
  have hone : v.valuation W.FunctionField (x₁ - x₂) = 1 :=
    val_eq_one_of_redConst_ne_zero W v h hv hsub hredsub
  rw [slope_of_X_ne hxne, slope_of_X_ne hne,
    redConst_div W v h hv (val_sub_le W v hy₁ hy₂) hone,
    redConst_sub W v h hv hy₁ hy₂, redConst_sub W v h hv hx₁ hx₂]

/-- ★★★`red x₁ ≠ red x₂` のとき傾きの付値は 1 以下。 -/
theorem val_slope_le_of_redX_ne {x₁ x₂ y₁ y₂ : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : redConst W v h hv x₁ ≠ redConst W v h hv x₂) :
    v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂) ≤ 1 := by
  have hxne : x₁ ≠ x₂ := fun hxeq => hne (by rw [hxeq])
  have hsub : v.valuation W.FunctionField (x₁ - x₂) ≤ 1 := val_sub_le W v hx₁ hx₂
  have hredsub : redConst W v h hv (x₁ - x₂) ≠ 0 := by
    rw [redConst_sub W v h hv hx₁ hx₂]
    exact sub_ne_zero.2 hne
  have hone : v.valuation W.FunctionField (x₁ - x₂) = 1 :=
    val_eq_one_of_redConst_ne_zero W v h hv hsub hredsub
  rw [slope_of_X_ne hxne, div_eq_mul_inv, Valuation.map_mul, map_inv₀, hone, inv_one, mul_one]
  exact val_sub_le W v hy₁ hy₂

/-- ★★★★★★**還元は加法を保つ**(`red x₁ ≠ red x₂` の場合)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★還元先でも弦の公式が使えるので、`addX`・`addY` の可換性(第 155)で閉じる。 -/
theorem redPoint_add_of_redX_ne [W.IsElliptic] {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hne : redConst W v h hv x₁ ≠ redConst W v h hv x₂) :
    redPoint W v h hv (Point.some x₁ y₁ h₁ + Point.some x₂ y₂ h₂)
      = redPoint W v h hv (Point.some x₁ y₁ h₁) + redPoint W v h hv (Point.some x₂ y₂ h₂) := by
  have hxne : x₁ ≠ x₂ := fun hxeq => hne (by rw [hxeq])
  have hxy : ¬(x₁ = x₂ ∧ y₁ = (W.map (algebraMap F W.FunctionField)).negY x₂ y₂) :=
    fun hh => hxne hh.1
  have hxyF : ¬(redConst W v h hv x₁ = redConst W v h hv x₂ ∧
      redConst W v h hv y₁ = W.negY (redConst W v h hv x₂) (redConst W v h hv y₂)) :=
    fun hh => hne hh.1
  set l := (W.map (algebraMap F W.FunctionField)).slope x₁ x₂ y₁ y₂ with hl
  have hlv : v.valuation W.FunctionField l ≤ 1 :=
    val_slope_le_of_redX_ne W v h hv hx₁ hx₂ hy₁ hy₂ hne
  rw [Point.add_some hxy, redPoint_some W v h hv _ (val_addX_le W v hx₁ hx₂ hlv)
      (val_addY_le W v hx₁ hx₂ hy₁ hlv),
    redPoint_some W v h hv h₁ hx₁ hy₁, redPoint_some W v h hv h₂ hx₂ hy₂,
    Point.add_some hxyF]
  refine point_some_congr ?_ ?_ _ _
  · rw [redConst_addX W v h hv hx₁ hx₂ hlv, hl,
      redConst_slope_of_redX_ne W v h hv hx₁ hx₂ hy₁ hy₂ hne]
  · rw [redConst_addY W v h hv hx₁ hx₂ hy₁ hlv, hl,
      redConst_slope_of_redX_ne W v h hv hx₁ hx₂ hy₁ hy₂ hne]

/-! ## ★出典の紐付け(`.src`) -/

def redPoint_add_of_redX_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——還元が加法を保つこと(弦の場合))",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
