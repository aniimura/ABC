import ABC3.Found.GaloisRep.Psi2Eval

/-!
# Galois (G1) 第 30 ブロック —— **★★★★★倍化公式 `x(2P) = Φ₂/Ψ₂Sq`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★乗法公式の基底

`[n](x,y) = (Φₙ/ΨSqₙ, …)` の帰納の基底が `n = 2` である。

    x(2P) · Ψ₂Sq(x) = Φ₂(x)

★分母を払った形で述べるので、`Ψ₂Sq(x) = 0`(= 2-捩れ)の場合も**式として正しい**。

## ★★通し方——**除算を先に消す**

`field_simp` を最後にかけると式の並びが変わって `linear_combination` が刺さらない
(§9-368 で実測)。★★**`d := y − negY x y` を先に払ってから**曲線の方程式を使う:

| 段 | 内容 |
|---|---|
| 1 | `addX x x (N/d) · d² = N² + a₁Nd − (a₂+2x)d²`(ここだけ `field_simp`) |
| 2 | `Ψ₂Sq(x) = d²`(第 29) |
| 3 | 残りは多項式恒等式——証明証明書 `c = −(a₁² + 4a₂ + 8x)` |

★★★係数 `c` は Python で曲線の方程式による剰余計算から得た(余り 0)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `doubling_x` | ★★★★★**`x(2P) · Ψ₂Sq(x) = Φ₂(x)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★**倍化公式**(分母を払った形)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが乗法公式の帰納の基底である。 -/
theorem doubling_x {x y : F} (h : W.toAffine.Equation x y) (hy : y ≠ W.toAffine.negY x y) :
    W.toAffine.addX x x (W.toAffine.slope x x y y) * W.Ψ₂Sq.eval x = (W.Φ 2).eval x := by
  have hd : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hy
  have haddX : W.toAffine.addX x x
        ((3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y) / (y - W.toAffine.negY x y))
      * (y - W.toAffine.negY x y) ^ 2
      = (3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y) ^ 2
        + W.a₁ * (3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y) * (y - W.toAffine.negY x y)
        - (W.a₂ + 2 * x) * (y - W.toAffine.negY x y) ^ 2 := by
    rw [WeierstrassCurve.Affine.addX]
    field_simp
    ring
  rw [W.toAffine.slope_of_Y_ne rfl hy, psi2Sq_eval W h, haddX, WeierstrassCurve.Φ_two,
    WeierstrassCurve.Affine.negY]
  rw [WeierstrassCurve.Affine.Equation] at h
  simp only [WeierstrassCurve.Affine.polynomial, evalEval, eval_add, eval_sub, eval_mul,
    eval_pow, eval_C, eval_X] at h
  simp only [eval_sub, eval_mul, eval_pow, eval_C, eval_X]
  rw [WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  linear_combination (-(W.a₁ ^ 2 + 4 * W.a₂ + 8 * x)) * h


/-! ## ★出典の紐付け(`.src`) -/

def doubling_x.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の基底——倍化公式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
