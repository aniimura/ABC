import ABC3.Found.GaloisRep.RedPoint

/-!
# Galois (G5) 第 155 ブロック —— **★★★★★★加法公式の各部品が還元と可換**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★加法公式は多項式だから還元と可換である

mathlib の加法公式は次の 4 つで書かれている:

    negY x y      = −y − a₁x − a₃
    addX x₁ x₂ ℓ  = ℓ² + a₁ℓ − a₂ − x₁ − x₂
    negAddY ...   = ℓ(addX − x₁) + y₁
    addY ...      = negY (addX) (negAddY)

★**どれも係数が `F` にある多項式**なので、第 153 の `red` が環準同型であることから
そのまま還元と可換になる。★★付値が 1 以下であることも保たれる。

## ★★★★これで残るのは傾き `ℓ` だけ

    red(slope x₁ x₂ y₁ y₂) = slope (red x₁) (red x₂) (red y₁) (red y₂)

★`red x₁ ≠ red x₂` の場合は `slope = (y₁−y₂)/(x₁−x₂)` で分母の付値が 1 だから
第 153 の `redConst_div` がそのまま効く。
★★残るのは `red x₁ = red x₂` の場合(還元先で 2 倍公式に切り替わる)である。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `val_neg_le` / `val_sub_le` | ★付値が 1 以下であることの保存 |
| `redConst_eq_zero_of_lt_one` | ★★`w t < 1 ⟹ red t = 0` |
| `val_eq_one_of_redConst_ne_zero` | ★★★**`red t ≠ 0 ⟹ w t = 1`**(傾きの分母に効く) |
| `redConst_negY` / `redConst_addX` | ★★★★還元との可換 |
| `redConst_negAddY` / `redConst_addY` | ★★★★★同上 |
| `val_negY_le` / `val_addX_le` / `val_negAddY_le` / `val_addY_le` | ★★付値の保存 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing] (v : HeightOneSpectrum W.CoordinateRing)

/-! ## ★付値が 1 以下であることの保存 -/

theorem val_neg_le {t : W.FunctionField} (ht : v.valuation W.FunctionField t ≤ 1) :
    v.valuation W.FunctionField (-t) ≤ 1 := by rwa [Valuation.map_neg]

theorem val_sub_le {s t : W.FunctionField} (hs : v.valuation W.FunctionField s ≤ 1)
    (ht : v.valuation W.FunctionField t ≤ 1) :
    v.valuation W.FunctionField (s - t) ≤ 1 := by
  rw [sub_eq_add_neg]; exact val_add_le W v hs (val_neg_le W v ht)

theorem val_negY_le {x y : W.FunctionField} (hx : v.valuation W.FunctionField x ≤ 1)
    (hy : v.valuation W.FunctionField y ≤ 1) :
    v.valuation W.FunctionField ((W.map (algebraMap F W.FunctionField)).negY x y) ≤ 1 := by
  rw [WeierstrassCurve.Affine.negY, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃]
  exact val_sub_le W v (val_sub_le W v (val_neg_le W v hy)
    (val_mul_le W v (valuation_algebraMap_field W v _) hx)) (valuation_algebraMap_field W v _)

theorem val_addX_le {x₁ x₂ l : W.FunctionField} (hx₁ : v.valuation W.FunctionField x₁ ≤ 1)
    (hx₂ : v.valuation W.FunctionField x₂ ≤ 1) (hl : v.valuation W.FunctionField l ≤ 1) :
    v.valuation W.FunctionField ((W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l) ≤ 1 := by
  rw [WeierstrassCurve.Affine.addX, WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂]
  exact val_sub_le W v (val_sub_le W v (val_sub_le W v
    (val_add_le W v (val_pow_le W v hl 2)
      (val_mul_le W v (valuation_algebraMap_field W v _) hl))
    (valuation_algebraMap_field W v _)) hx₁) hx₂

theorem val_negAddY_le {x₁ x₂ y₁ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hl : v.valuation W.FunctionField l ≤ 1) :
    v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).negAddY x₁ x₂ y₁ l) ≤ 1 := by
  rw [WeierstrassCurve.Affine.negAddY]
  exact val_add_le W v (val_mul_le W v hl
    (val_sub_le W v (val_addX_le W v hx₁ hx₂ hl) hx₁)) hy₁

theorem val_addY_le {x₁ x₂ y₁ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hl : v.valuation W.FunctionField l ≤ 1) :
    v.valuation W.FunctionField
      ((W.map (algebraMap F W.FunctionField)).addY x₁ x₂ y₁ l) ≤ 1 := by
  rw [WeierstrassCurve.Affine.addY]
  exact val_negY_le W v (val_addX_le W v hx₁ hx₂ hl) (val_negAddY_le W v hx₁ hx₂ hy₁ hl)

/-! ## ★★★還元が 0 になる条件 -/

variable {c y₀ : F} (h : W.Equation c y₀)
  (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))

theorem redConst_eq_zero_of_lt_one {t : W.FunctionField}
    (ht : v.valuation W.FunctionField t < 1) : redConst W v h hv t = 0 := by
  refine redConst_eq W v h hv (le_of_lt ht) ?_
  rwa [map_zero, sub_zero]

/-- ★★★**`red t ≠ 0` なら `w t = 1`**——傾きの分母に効く。 -/
theorem val_eq_one_of_redConst_ne_zero {t : W.FunctionField}
    (ht : v.valuation W.FunctionField t ≤ 1) (hne : redConst W v h hv t ≠ 0) :
    v.valuation W.FunctionField t = 1 := by
  refine le_antisymm ht ?_
  by_contra hlt
  rw [not_le] at hlt
  exact hne (redConst_eq_zero_of_lt_one W v h hv hlt)

/-! ## ★★★★加法公式の各部品が還元と可換 -/

/-- ★★★★`red(negY x y) = negY (red x) (red y)`。 -/
theorem redConst_negY {x y : W.FunctionField}
    (hx : v.valuation W.FunctionField x ≤ 1) (hy : v.valuation W.FunctionField y ≤ 1) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).negY x y)
      = W.negY (redConst W v h hv x) (redConst W v h hv y) := by
  have hA1 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁) ≤ 1 :=
    valuation_algebraMap_field W v _
  have hA3 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₃) ≤ 1 :=
    valuation_algebraMap_field W v _
  have h1 : v.valuation W.FunctionField (-y) ≤ 1 := val_neg_le W v hy
  have h2 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁ * x) ≤ 1 :=
    val_mul_le W v hA1 hx
  have h3 : v.valuation W.FunctionField
      (-y - algebraMap F W.FunctionField W.a₁ * x) ≤ 1 := val_sub_le W v h1 h2
  rw [WeierstrassCurve.Affine.negY, WeierstrassCurve.Affine.negY,
    WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₃,
    redConst_sub W v h hv h3 hA3, redConst_sub W v h hv h1 h2,
    redConst_neg W v h hv hy, redConst_mul W v h hv hA1 hx,
    redConst_algebraMap, redConst_algebraMap]

/-- ★★★★`red(addX x₁ x₂ ℓ) = addX (red x₁) (red x₂) (red ℓ)`。 -/
theorem redConst_addX {x₁ x₂ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hl : v.valuation W.FunctionField l ≤ 1) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).addX x₁ x₂ l)
      = W.addX (redConst W v h hv x₁) (redConst W v h hv x₂) (redConst W v h hv l) := by
  have hA1 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁) ≤ 1 :=
    valuation_algebraMap_field W v _
  have hA2 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₂) ≤ 1 :=
    valuation_algebraMap_field W v _
  have h1 : v.valuation W.FunctionField (l ^ 2) ≤ 1 := val_pow_le W v hl 2
  have h2 : v.valuation W.FunctionField (algebraMap F W.FunctionField W.a₁ * l) ≤ 1 :=
    val_mul_le W v hA1 hl
  have h3 : v.valuation W.FunctionField (l ^ 2 + algebraMap F W.FunctionField W.a₁ * l) ≤ 1 :=
    val_add_le W v h1 h2
  have h4 : v.valuation W.FunctionField
      (l ^ 2 + algebraMap F W.FunctionField W.a₁ * l
        - algebraMap F W.FunctionField W.a₂) ≤ 1 := val_sub_le W v h3 hA2
  have h5 : v.valuation W.FunctionField
      (l ^ 2 + algebraMap F W.FunctionField W.a₁ * l
        - algebraMap F W.FunctionField W.a₂ - x₁) ≤ 1 := val_sub_le W v h4 hx₁
  rw [WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.addX,
    WeierstrassCurve.map_a₁, WeierstrassCurve.map_a₂,
    redConst_sub W v h hv h5 hx₂, redConst_sub W v h hv h4 hx₁,
    redConst_sub W v h hv h3 hA2, redConst_add W v h hv h1 h2,
    redConst_pow W v h hv hl 2, redConst_mul W v h hv hA1 hl,
    redConst_algebraMap, redConst_algebraMap]

/-- ★★★★★`red(negAddY ...) = negAddY (red ...)`。 -/
theorem redConst_negAddY {x₁ x₂ y₁ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hl : v.valuation W.FunctionField l ≤ 1) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).negAddY x₁ x₂ y₁ l)
      = W.negAddY (redConst W v h hv x₁) (redConst W v h hv x₂) (redConst W v h hv y₁)
        (redConst W v h hv l) := by
  rw [WeierstrassCurve.Affine.negAddY, WeierstrassCurve.Affine.negAddY,
    redConst_add W v h hv (val_mul_le W v hl
      (val_sub_le W v (val_addX_le W v hx₁ hx₂ hl) hx₁)) hy₁,
    redConst_mul W v h hv hl (val_sub_le W v (val_addX_le W v hx₁ hx₂ hl) hx₁),
    redConst_sub W v h hv (val_addX_le W v hx₁ hx₂ hl) hx₁,
    redConst_addX W v h hv hx₁ hx₂ hl]

/-- ★★★★★`red(addY ...) = addY (red ...)`。 -/
theorem redConst_addY {x₁ x₂ y₁ l : W.FunctionField}
    (hx₁ : v.valuation W.FunctionField x₁ ≤ 1) (hx₂ : v.valuation W.FunctionField x₂ ≤ 1)
    (hy₁ : v.valuation W.FunctionField y₁ ≤ 1) (hl : v.valuation W.FunctionField l ≤ 1) :
    redConst W v h hv ((W.map (algebraMap F W.FunctionField)).addY x₁ x₂ y₁ l)
      = W.addY (redConst W v h hv x₁) (redConst W v h hv x₂) (redConst W v h hv y₁)
        (redConst W v h hv l) := by
  rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.addY,
    redConst_negY W v h hv (val_addX_le W v hx₁ hx₂ hl)
      (val_negAddY_le W v hx₁ hx₂ hy₁ hl),
    redConst_addX W v h hv hx₁ hx₂ hl,
    redConst_negAddY W v h hv hx₁ hx₂ hy₁ hl]

/-! ## ★出典の紐付け(`.src`) -/

def redConst_addY.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——加法公式の各部品が還元と可換であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
