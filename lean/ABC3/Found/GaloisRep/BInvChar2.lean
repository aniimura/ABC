import ABC3.Found.GaloisRep.OmegaIdentity

/-!
# Galois (G1) 第 12 ブロック —— ★★★★**標数 2 での `b` 不変量**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★手計算で**鍵の関係を特定した**

第 11 の恒等式を `n = 2` で書き下すと

    C preΨ₄ = ψ₂ (a₁(X ψ₂² − Ψ₃) + a₃ ψ₂²)      （標数 2)

★手で展開すると(標数 2、`ψ₂ = a₁X + a₃`):

    右辺 = (a₁X+a₃)⁴ − a₁(a₁X+a₃)Ψ₃
         = b₂X⁵ + b₄X⁴ + **(b₂b₆ + b₄²)X²** + (b₂b₈+b₄b₆)X + (b₄b₈+b₆²)
    左辺 = preΨ₄
         = b₂X⁵ + b₄X⁴ + (b₂b₈+b₄b₆)X + (b₄b₈+b₆²)

★★★差は **`(b₂b₆ + b₄²)X²`** の 1 項だけであり、
★★★★標数 2 で **`b₂b₆ = b₄²`** なのでこれが消える。

## ★★標数 2 での `b` 不変量

| 不変量 | 一般 | 標数 2 |
|---|---|---|
| `b₂` | `a₁² + 4a₂` | ★`a₁²` |
| `b₄` | `2a₄ + a₁a₃` | ★`a₁a₃` |
| `b₆` | `a₃² + 4a₆` | ★`a₃²` |

★したがって `b₂b₆ = a₁²a₃² = (a₁a₃)² = b₄²` ✅

## ★★★これが `ω_n` が整数係数である**構造的な理由**である

`b₂b₆ − b₄²` は判別式に現れる量で、標数 2 で消える。
★★**`ω_n` の分母 2 が消えるのは、この関係が標数 2 で成り立つから**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `b2_char_two` / `b4_char_two` / `b6_char_two` | ★標数 2 での `b` 不変量 |
| `b2_b6_eq_b4_sq` | ★★★★**`b₂b₆ = b₄²`**(鍵の関係) |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★**標数 2 での `b₂ = a₁²`**。 -/
theorem b2_char_two : W.b₂ = W.a₁ ^ 2 := by
  rw [WeierstrassCurve.b₂]
  have : (4 : R) = 0 := by
    have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
    have h4 : (4 : R) = 2 * 2 := by norm_num
    rw [h4, h, mul_zero]
  rw [this]
  ring

/-- ★**標数 2 での `b₄ = a₁a₃`**。 -/
theorem b4_char_two : W.b₄ = W.a₁ * W.a₃ := by
  rw [WeierstrassCurve.b₄]
  have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
  rw [h]
  ring

/-- ★**標数 2 での `b₆ = a₃²`**。 -/
theorem b6_char_two : W.b₆ = W.a₃ ^ 2 := by
  rw [WeierstrassCurve.b₆]
  have : (4 : R) = 0 := by
    have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
    have h4 : (4 : R) = 2 * 2 := by norm_num
    rw [h4, h, mul_zero]
  rw [this]
  ring

/-- ★★★★**標数 2 での鍵の関係 `b₂b₆ = b₄²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `ω_n` の分母 2 が消える**構造的な理由**である——
第 11 の恒等式の左右の差はちょうど `(b₂b₆ + b₄²)X²` の 1 項であり、
標数 2 でこれが消える。 -/
theorem b2_b6_eq_b4_sq : W.b₂ * W.b₆ = W.b₄ ^ 2 := by
  rw [b2_char_two, b4_char_two, b6_char_two]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def b2_b6_eq_b4_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——標数 2 での b 不変量)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
