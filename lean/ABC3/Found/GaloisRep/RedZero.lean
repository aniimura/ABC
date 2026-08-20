import ABC3.Found.GaloisRep.RedAdd

/-!
# Galois (G5) 第 159 ブロック —— **★★★★★`red D = 0` の場合の部品**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★残る場合の両辺はどちらも 0 になる

第 158 で `red D ≠ 0` の場合が閉じた。★残るのは

    red x₁ = red x₂   かつ   red D = 0

の場合で、そこでは**両辺がどちらも 0** になる。

### ★★★右辺(還元先の和)

`red D = 0` は `red y₁ + red y₂ + a₁·red x₁ + a₃ = 0`、すなわち
`red y₁ = negY (red x₁) (red y₂)`。★`red x₁ = red x₂` と合わせて
mathlib の `Point.add_of_Y_eq` がそのまま効く。

### ★★★左辺(和の還元)

`w(D) < 1` なので傾き `slope = N/D` は**極を持つ**(`w(N) = 1` のとき)。
★そのとき `addX = ℓ² + a₁ℓ − a₂ − x₁ − x₂` の付値は `w(ℓ)² > 1` になる
——`ℓ²` の項が**単独で最大**だからである。
★★座標が付値環に入らないので `redPoint = 0`(定義の `else` 枝)。

## ★★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `one_lt_val_addX` | ★★★★★**傾きが極を持てば `addX` も極を持つ** |
| `redPoint_eq_zero_of_not_le` | ★★座標が付値環に入らなければ `redPoint = 0` |
| `redPoint_add_eq_zero_of_redD` | ★★★★**`red D = 0` なら還元先の和は 0** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-- ★★★★★**傾きが極を持てば `addX` も極を持つ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ℓ²` の項が単独で最大になるので、`w(addX) = w(ℓ)² > 1`。 -/
theorem one_lt_val_addX {x₁ x₂ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hl : 1 < v.valuation W.FunctionField l) :
    1 < v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l) := by
  set w := v.valuation W.FunctionField with hw
  have hb : ∀ {t : W.FunctionField}, w t ≤ 1 → w t ≤ w l := fun ht => le_trans ht (le_of_lt hl)
  have hsq : w l < w l ^ 2 := by
    have hp := pow_lt_pow_right₀ hl (by norm_num : 1 < 2)
    rwa [pow_one] at hp
  have hrest : w (algebraMap F W.FunctionField W.a₁ * l
      - algebraMap F W.FunctionField W.a₂ - x₁ - x₂) ≤ w l := by
    refine le_trans (Valuation.map_sub w _ _) (max_le ?_ (hb hx₂))
    refine le_trans (Valuation.map_sub w _ _) (max_le ?_ (hb hx₁))
    refine le_trans (Valuation.map_sub w _ _) (max_le ?_ (hb (valuation_algebraMap_field W v _)))
    rw [Valuation.map_mul]
    calc w (algebraMap F W.FunctionField W.a₁) * w l ≤ 1 * w l :=
          mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
      _ = w l := one_mul _
  have hdec : (W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l
      = l ^ 2 + (algebraMap F W.FunctionField W.a₁ * l
        - algebraMap F W.FunctionField W.a₂ - x₁ - x₂) := by
    rw [WeierstrassCurve.Affine.addX, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂]
    ring
  rw [hdec, Valuation.map_add_of_distinct_val w (by
    rw [Valuation.map_pow]
    exact ne_of_gt (lt_of_le_of_lt hrest hsq)), Valuation.map_pow,
    max_eq_left (le_of_lt (lt_of_le_of_lt hrest hsq))]
  exact lt_trans hl hsq

variable {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

/-- ★★座標が付値環に入らなければ `redPoint = 0`。 -/
theorem redPoint_eq_zero_of_not_le [W.IsElliptic] {x y : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular x y)
    (hxy : ¬(v.valuation W.FunctionField x ≤ 1 ∧ v.valuation W.FunctionField y ≤ 1)) :
    redPoint W v h hv (Point.some x y hns) = 0 := by
  rw [redPoint, dif_neg hxy]

variable [DecidableEq F]

/-- ★★★★**`red D = 0` かつ `red x₁ = red x₂` なら還元先の和は 0**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`red D = 0` から `red y₁ = negY (red x₂) (red y₂)` が直ちに出る。 -/
theorem redPoint_add_eq_zero_of_redD [W.IsElliptic] {x₁ y₁ x₂ y₂ : W.FunctionField}
    (h₁ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₁ y₁)
    (h₂ : (W.map (algebraMap F W.FunctionField)).Nonsingular x₂ y₂)
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hy₁ : v.valuation W.FunctionField y₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hy₂ : v.valuation W.FunctionField y₂ ≤ 1)
    (hxe : redConst W v h hv x₁ = redConst W v h hv x₂)
    (hD : redConst W v h hv (y₁ + y₂ + (W.map (algebraMap F W.FunctionField)).a₁ * x₁
        + (W.map (algebraMap F W.FunctionField)).a₃) = 0) :
    redPoint W v h hv (Point.some x₁ y₁ h₁) + redPoint W v h hv (Point.some x₂ y₂ h₂) = 0 := by
  rw [redPoint_some W v h hv h₁ hx₁ hy₁, redPoint_some W v h hv h₂ hx₂ hy₂]
  refine Point.add_of_Y_eq hxe ?_
  rw [redConst_slopeDen W v h hv hx₁ hy₁ hy₂] at hD
  rw [WeierstrassCurve.Affine.negY, ← hxe]
  linear_combination hD

/-! ## ★出典の紐付け(`.src`) -/

def one_lt_val_addX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——傾きが極を持てば addX も極を持つこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
