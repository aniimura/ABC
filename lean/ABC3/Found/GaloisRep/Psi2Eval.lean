import ABC3.Found.GaloisRep.OmegaCurve
import Mathlib.AlgebraicGeometry.EllipticCurve.DivisionPolynomial.Basic

/-!
# Galois (G1) 第 29 ブロック —— **曲線の上では `Ψ₂Sq(x) = (y − negY)²`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★乗法公式の第 1 歩

`[n](x,y) = (Φₙ/ΨSqₙ, …)` を示すには、まず **`ΨSq₂ = Ψ₂Sq` が
点の上で `(2y + a₁x + a₃)²` になる**ことが要る。

★mathlib は `ψ₂_sq : ψ₂² = C Ψ₂Sq + 4·polynomial` を多項式として持つが、
**点で評価した形**は無い。

★★曲線の方程式 `E := y² + a₁xy + a₃y − x³ − a₂x² − a₄x − a₆ = 0` を使って
`linear_combination (-4) * E` で出る。

## ★★これが効く場所

`y ≠ negY x y`(= `2y + a₁x + a₃ ≠ 0`)は倍化が定義できる条件そのものなので、
★**`Ψ₂Sq(x) ≠ 0` と同値**になる——`E[2]` の記述に直結する。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★★★★**曲線の上では `Ψ₂Sq(x) = (y − negY x y)²`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★`2y + a₁x + a₃` が `Ψ₂Sq` の平方根になる——2-捩れの判定はこれである。 -/
theorem psi2Sq_eval {x y : F} (h : W.toAffine.Equation x y) :
    W.Ψ₂Sq.eval x = (y - W.toAffine.negY x y) ^ 2 := by
  rw [WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.Affine.negY]
  rw [WeierstrassCurve.Affine.Equation] at h
  simp only [WeierstrassCurve.Affine.polynomial, evalEval, eval_add, eval_sub, eval_mul,
    eval_pow, eval_C, eval_X] at h
  simp only [eval_add, eval_mul, eval_pow, eval_C, eval_X]
  rw [WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆]
  linear_combination (-4 : F) * h

/-- ★**2-捩れの判定**——`Ψ₂Sq(x) = 0` は `y = negY x y` と同値。 -/
theorem psi2Sq_eval_eq_zero_iff {x y : F} (h : W.toAffine.Equation x y) :
    W.Ψ₂Sq.eval x = 0 ↔ y = W.toAffine.negY x y := by
  rw [psi2Sq_eval W h, pow_eq_zero_iff two_ne_zero, sub_eq_zero]

/-! ## ★出典の紐付け(`.src`) -/

def psi2Sq_eval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式——Ψ₂Sq の点での評価)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
