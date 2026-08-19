import ABC3.Found.GaloisRep.FiniteTwoTorsion

/-!
# Galois (G1) 第 34 ブロック —— **1 つの `x` の上の `y` は有限**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★一般の `n` の有限性に要る 3 つのうち 2 つ目

`torsion_finite` を一般の `n` で示すには 3 つ要る(§9-372 で測った):

| 要るもの | 状態 |
|---|---|
| `n • P = 0 ⟺ ΨSqₙ(x_P) = 0` | ★乗法公式が要る(残り) |
| **1 つの `x` の上の `y` は有限** | ★★**本ブロック** |
| `ΨSqₙ` の根の有限性 | ★mathlib の在庫 |

## ★★Weierstrass の方程式は `y` の 2 次式である

    y² + (a₁x + a₃)y − (x³ + a₂x² + a₄x + a₆) = 0

★モニックなので零多項式でない——根は有限個。

★★`n = 2` では `y` が**一意**に決まった(第 32)が、一般の `n` では
**高々 2 つ**で十分である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- 与えられた `x` に対する `y` の方程式(`y` の 2 次式)。 -/
noncomputable def yPoly (x : F) : F[X] :=
  X ^ 2 + C (W.toAffine.a₁ * x + W.toAffine.a₃) * X
    - C (x ^ 3 + W.toAffine.a₂ * x ^ 2 + W.toAffine.a₄ * x + W.toAffine.a₆)

theorem yPoly_ne_zero (x : F) : yPoly W x ≠ 0 := by
  intro h
  have hd : (yPoly W x).coeff 2 = 1 := by
    simp [yPoly, coeff_sub, coeff_add, coeff_X_pow, coeff_C_mul, coeff_X, coeff_C,
      -map_add, -map_mul, -map_pow, -map_ofNat]
  rw [h] at hd
  simp at hd

theorem equation_iff_isRoot (x y : F) :
    W.toAffine.Equation x y ↔ (yPoly W x).IsRoot y := by
  rw [WeierstrassCurve.Affine.Equation, yPoly, Polynomial.IsRoot]
  simp only [WeierstrassCurve.Affine.polynomial, evalEval, eval_add, eval_sub, eval_mul,
    eval_pow, eval_C, eval_X]

theorem finite_y_of_x (x : F) : {y : F | W.toAffine.Equation x y}.Finite := by
  have : {y : F | W.toAffine.Equation x y} = {y : F | (yPoly W x).IsRoot y} := by
    ext y; exact equation_iff_isRoot W x y
  rw [this]
  exact Polynomial.finite_setOf_isRoot (yPoly_ne_zero W x)


/-! ## ★出典の紐付け(`.src`) -/

def finite_y_of_x.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(1 つの x の上の y は有限)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
