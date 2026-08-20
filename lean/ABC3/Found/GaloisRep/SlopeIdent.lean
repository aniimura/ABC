import ABC3.Found.GaloisRep.RedSlope

/-!
# Galois (G5) 第 157 ブロック —— **★★★★★★傾きの統一表示**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★弦と接線が**同じ有理式**で書ける

mathlib の `slope` は場合分けで定義されている:

* `x₁ ≠ x₂`(弦): `(y₁ − y₂)/(x₁ − x₂)`
* `x₁ = x₂`、`y₁ ≠ negY x₂ y₂`(接線): `(3x₁² + 2a₂x₁ + a₄ − a₁y₁)/(2y₁ + a₁x₁ + a₃)`

★**曲線上の 2 点に対しては、どちらも次の 1 つの式に等しい**:

    slope = (x₁² + x₁x₂ + x₂² + a₂(x₁+x₂) + a₄ − a₁y₂) / (y₁ + y₂ + a₁x₁ + a₃)

★★根拠は 2 つの方程式の差から出る恒等式

    (y₁ − y₂)(y₁ + y₂ + a₁x₁ + a₃)
      = (x₁ − x₂)(x₁² + x₁x₂ + x₂² + a₂(x₁+x₂) + a₄ − a₁y₂)

である(`linear_combination h₁ − h₂` の 1 行)。

## ★★★★★これが還元の場合分けを潰す

分子・分母が**多項式**なので、第 153–155 の道具がそのまま効く:

    red(D) ≠ 0  ⟹  red(slope) = slope(red …)

★弦か接線かを問わない。★★`red x₁ = red x₂` でも `red(D) ≠ 0` なら
還元先の 2 倍公式がそのまま出る——**これが D2 の残りの山を崩す**。

★★★`red(D) = 0` の場合は、還元先で `red S₁ + red S₂ = 0` になる
(`red y₁ = negY (red x₂) (red y₂)` が `red D = 0` から直ちに出る)。

## ★★★これが (G7) でも効く

(G7) 半安定モデルも点の還元を要求する。★ここで積んだものはそのまま流用できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `slope_num_identity` | ★★★★★★**恒等式**(任意の可換環で成立) |
| `slope_eq_of_denom_ne_zero` | ★★★★★弦の場合の第 2 表示 |
| `slope_self_eq` | ★★★★接線の場合も同じ形 |
| `slope_eq_div` | ★★★★★★**統一表示** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine

/-- ★★★★★★**曲線上の 2 点に対する傾きの恒等式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★2 つの Weierstrass 方程式の差から `linear_combination` の 1 行で出る。
★★任意の可換環で成立する。 -/
theorem slope_num_identity {R : Type} [CommRing R] (W : WeierstrassCurve.Affine R)
    {x₁ x₂ y₁ y₂ : R} (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂) :
    (y₁ - y₂) * (y₁ + y₂ + W.a₁ * x₁ + W.a₃)
      = (x₁ - x₂) * (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 + W.a₂ * (x₁ + x₂) + W.a₄ - W.a₁ * y₂) := by
  rw [equation_iff] at h₁ h₂
  linear_combination h₁ - h₂

/-- ★★★★★弦の場合の第 2 表示。 -/
theorem slope_eq_of_denom_ne_zero {F : Type} [Field F] [DecidableEq F]
    (W : WeierstrassCurve.Affine F) {x₁ x₂ y₁ y₂ : F}
    (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂) (hx : x₁ ≠ x₂)
    (hd : y₁ + y₂ + W.a₁ * x₁ + W.a₃ ≠ 0) :
    W.slope x₁ x₂ y₁ y₂
      = (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 + W.a₂ * (x₁ + x₂) + W.a₄ - W.a₁ * y₂)
        / (y₁ + y₂ + W.a₁ * x₁ + W.a₃) := by
  have hid := slope_num_identity W h₁ h₂
  have hxne : x₁ - x₂ ≠ 0 := sub_ne_zero.2 hx
  rw [slope_of_X_ne hx, div_eq_div_iff hxne hd]
  linear_combination hid

/-- ★★★★接線の場合も同じ形になる。 -/
theorem slope_self_eq {F : Type} [Field F] [DecidableEq F]
    (W : WeierstrassCurve.Affine F) {x y : F} (hy : y ≠ W.negY x y) :
    W.slope x x y y
      = (x ^ 2 + x * x + x ^ 2 + W.a₂ * (x + x) + W.a₄ - W.a₁ * y)
        / (y + y + W.a₁ * x + W.a₃) := by
  rw [slope_of_Y_ne rfl hy, WeierstrassCurve.Affine.negY]
  congr 1 <;> ring

/-- ★★★★★★**傾きの統一表示**——弦も接線も同じ有理式で書ける。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★分子・分母が多項式なので、還元との可換性が場合分けなしに従う。 -/
theorem slope_eq_div {F : Type} [Field F] [DecidableEq F]
    (W : WeierstrassCurve.Affine F) {x₁ x₂ y₁ y₂ : F}
    (h₁ : W.Equation x₁ y₁) (h₂ : W.Equation x₂ y₂)
    (hne : ¬(x₁ = x₂ ∧ y₁ = W.negY x₂ y₂))
    (hd : y₁ + y₂ + W.a₁ * x₁ + W.a₃ ≠ 0) :
    W.slope x₁ x₂ y₁ y₂
      = (x₁ ^ 2 + x₁ * x₂ + x₂ ^ 2 + W.a₂ * (x₁ + x₂) + W.a₄ - W.a₁ * y₂)
        / (y₁ + y₂ + W.a₁ * x₁ + W.a₃) := by
  by_cases hx : x₁ = x₂
  · have hy : y₁ ≠ W.negY x₂ y₂ := fun hh => hne ⟨hx, hh⟩
    have hyy : y₁ = y₂ := Y_eq_of_Y_ne h₁ h₂ hx hy
    subst hx; subst hyy
    rw [slope_of_Y_ne rfl hy, WeierstrassCurve.Affine.negY]
    congr 1 <;> ring
  · exact slope_eq_of_denom_ne_zero W h₁ h₂ hx hd

/-! ## ★出典の紐付け(`.src`) -/

def slope_eq_div.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——弦と接線の傾きの統一表示)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
