import ABC3.Found.GaloisRep.TateCollinear

/-!
# Galois (G6) 第 257 ブロック —— **★★★★★★★★共線性から群法則へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★行列式 0 から `P₁ + P₂ + P₃ = 0` へ

第 256 で Tate 級数の 3 点の**共線性の行列式が 0** であることが出た。
本ブロックは、それを mathlib の群法則 `WeierstrassCurve.Affine.Point.add` に翻訳する。

    x₁(y₂−y₃) + x₂(y₃−y₁) + x₃(y₁−y₂) = 0
      かつ 3 点が曲線上、`x` 座標が相異なる
      ⟹ P₁ + P₂ + P₃ = 0

## ★★★★★★Vieta は「差の差」で出る

`yᵢ = ℓxᵢ + m` を曲線の式に代入すると、各 `xᵢ` は

    g(x) := x³ − B x² − C x − D,   B = ℓ² + a₁ℓ − a₂

の根になる。★★★**3 次式の根と係数の関係は、根の差で 2 回割るだけで出る**:

    (x₁−x₂)·(x₁²+x₁x₂+x₂² − B(x₁+x₂) − C) = 0     (= g(x₁) − g(x₂))
    (x₂−x₃)·(x₂²+x₂x₃+x₃² − B(x₂+x₃) − C) = 0
    ⟹ (x₁−x₃)·(x₁+x₂+x₃ − B) = 0

★因数分解も多項式環も要らない——`linear_combination` 3 回である。

## ★`x` 座標が相異なることは落とせない

`x₃ = x₁`(かつ `y₃ = y₁`、つまり `P₃ = P₁`)のとき、行列式は**自動的に 0** になり
何の情報も持たない。★このとき本当の第 3 根は別の点なので、結論は成り立たない。
**3 つの `x` 座標が相異なることは本質的な仮定**である
(Tate 側では葉 (d)(単射性)から来る)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_line_of_collinear` | ★★★★行列式 0 なら同一直線上 |
| `addX_eq_of_line` | ★★★★★★**Vieta**——`addX x₁ x₂ ℓ = x₃` |
| `add_eq_neg_of_line` | ★★★★★★★`P₁ + P₂ = −P₃` |
| `add_add_eq_zero_of_collDet` | ★★★★★★★★**行列式 0 ⟹ 和が 0** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

variable {F : Type} [Field F]

/-! ## ★★★★同一直線上にあること -/

/-- ★★★★**行列式が 0 で `x₁ ≠ x₂` なら 3 点は同じ直線の上にある**。 -/
theorem exists_line_of_collinear {x₁ y₁ x₂ y₂ x₃ y₃ : F} (h12 : x₁ ≠ x₂)
    (hd : x₁ * (y₂ - y₃) + x₂ * (y₃ - y₁) + x₃ * (y₁ - y₂) = 0) :
    ∃ l m : F, y₁ = l * x₁ + m ∧ y₂ = l * x₂ + m ∧ y₃ = l * x₃ + m := by
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 h12
  refine ⟨(y₁ - y₂) / (x₁ - x₂), y₁ - (y₁ - y₂) / (x₁ - x₂) * x₁, by ring, ?_, ?_⟩
  · field_simp
    ring
  · field_simp
    linear_combination -hd

/-! ## ★★★★★★Vieta -/

/-- ★★★★★★**共線な 3 点は Vieta を満たす**——`addX x₁ x₂ ℓ = x₃`。

★3 次式の根と係数の関係を、根の差で 2 回割るだけで出す。 -/
theorem addX_eq_of_line (W : WeierstrassCurve.Affine F) {x₁ y₁ x₂ y₂ x₃ y₃ l m : F}
    (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂) (h₃ : W.Equation x₃ y₃)
    (hy₁ : y₁ = l * x₁ + m) (hy₂ : y₂ = l * x₂ + m) (hy₃ : y₃ = l * x₃ + m)
    (h12 : x₁ ≠ x₂) (h13 : x₁ ≠ x₃) (h23 : x₂ ≠ x₃) :
    W.addX x₁ x₂ l = x₃ := by
  rw [equation_iff] at h₁ h₂ h₃
  subst hy₁ hy₂ hy₃
  set B : F := l ^ 2 + W.a₁ * l - W.a₂ with hB
  set C : F := 2 * l * m + W.a₁ * m + W.a₃ * l - W.a₄ with hC
  have g1 : x₁ ^ 3 - B * x₁ ^ 2 - C * x₁ - (m ^ 2 + W.a₃ * m - W.a₆) = 0 := by
    simp only [hB, hC]; linear_combination -h₁
  have g2 : x₂ ^ 3 - B * x₂ ^ 2 - C * x₂ - (m ^ 2 + W.a₃ * m - W.a₆) = 0 := by
    simp only [hB, hC]; linear_combination -h₂
  have g3 : x₃ ^ 3 - B * x₃ ^ 2 - C * x₃ - (m ^ 2 + W.a₃ * m - W.a₆) = 0 := by
    simp only [hB, hC]; linear_combination -h₃
  have e12 : (x₁ - x₂) * (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 - B * (x₁ + x₂) - C) = 0 := by
    linear_combination g1 - g2
  have e23 : (x₂ - x₃) * (x₂ ^ 2 + x₂ * x₃ + x₃ ^ 2 - B * (x₂ + x₃) - C) = 0 := by
    linear_combination g2 - g3
  have f12 : x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 - B * (x₁ + x₂) - C = 0 :=
    (mul_eq_zero.1 e12).resolve_left (sub_ne_zero.2 h12)
  have f23 : x₂ ^ 2 + x₂ * x₃ + x₃ ^ 2 - B * (x₂ + x₃) - C = 0 :=
    (mul_eq_zero.1 e23).resolve_left (sub_ne_zero.2 h23)
  have key : (x₁ - x₃) * (x₁ + x₂ + x₃ - B) = 0 := by linear_combination f12 - f23
  have hsum : x₁ + x₂ + x₃ - B = 0 :=
    (mul_eq_zero.1 key).resolve_left (sub_ne_zero.2 h13)
  rw [WeierstrassCurve.Affine.addX]
  simp only [hB] at hsum
  linear_combination -hsum

/-! ## ★★★★★★★★群法則 -/

/-- ★★★★★★★**共線な 3 点**——`P₁ + P₂ = −P₃`。 -/
theorem add_eq_neg_of_line [DecidableEq F] {W : WeierstrassCurve.Affine F}
    {x₁ y₁ x₂ y₂ x₃ y₃ l m : F}
    (n₁ : W.Nonsingular x₁ y₁) (n₂ : W.Nonsingular x₂ y₂) (n₃ : W.Nonsingular x₃ y₃)
    (hy₁ : y₁ = l * x₁ + m) (hy₂ : y₂ = l * x₂ + m) (hy₃ : y₃ = l * x₃ + m)
    (h12 : x₁ ≠ x₂) (h13 : x₁ ≠ x₃) (h23 : x₂ ≠ x₃) :
    Point.some x₁ y₁ n₁ + Point.some x₂ y₂ n₂ = -Point.some x₃ y₃ n₃ := by
  have e₁ : W.Equation x₁ y₁ := n₁.left
  have e₂ : W.Equation x₂ y₂ := n₂.left
  have e₃ : W.Equation x₃ y₃ := n₃.left
  have hne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 h12
  have hX := addX_eq_of_line W e₁ e₂ e₃ hy₁ hy₂ hy₃ h12 h13 h23
  have hsl : W.slope x₁ x₂ y₁ y₂ = l := by
    rw [slope_of_X_ne h12, hy₁, hy₂, div_eq_iff hne]
    ring
  rw [Point.add_of_X_ne h12, Point.neg_some]
  simp only [Point.some.injEq]
  refine ⟨by rw [hsl]; exact hX, ?_⟩
  rw [hsl, WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY, hX]
  congr 1
  rw [hy₁, hy₃]
  ring

/-- ★★★★★★★★**行列式が 0 なら 3 点の和は 0**——第 256 の共線性を群法則に翻訳する形。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem add_add_eq_zero_of_collDet [DecidableEq F] {W : WeierstrassCurve.Affine F}
    {x₁ y₁ x₂ y₂ x₃ y₃ : F}
    (n₁ : W.Nonsingular x₁ y₁) (n₂ : W.Nonsingular x₂ y₂) (n₃ : W.Nonsingular x₃ y₃)
    (hd : x₁ * (y₂ - y₃) + x₂ * (y₃ - y₁) + x₃ * (y₁ - y₂) = 0)
    (h12 : x₁ ≠ x₂) (h13 : x₁ ≠ x₃) (h23 : x₂ ≠ x₃) :
    Point.some x₁ y₁ n₁ + Point.some x₂ y₂ n₂ + Point.some x₃ y₃ n₃ = 0 := by
  obtain ⟨l, m, hy₁, hy₂, hy₃⟩ := exists_line_of_collinear h12 hd
  rw [add_eq_neg_of_line n₁ n₂ n₃ hy₁ hy₂ hy₃ h12 h13 h23]
  exact neg_add_cancel _

/-! ## ★出典の紐付け(`.src`) -/

def addX_eq_of_line.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線な 3 点の Vieta)",
    sectionId := "genell-def-3-3" }

def add_add_eq_zero_of_collDet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——行列式 0 から群法則へ)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
