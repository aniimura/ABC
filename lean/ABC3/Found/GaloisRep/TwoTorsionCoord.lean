import ABC3.Found.GaloisRep.TwoTorsion
import Mathlib.AlgebraicGeometry.EllipticCurve.DivisionPolynomial.Degree

/-!
# Galois (G1) 第 32 ブロック —— **2-捩れは `x` で決まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★`E[2]` の有限性の骨

第 31 で `2 • P = 0 ⟺ Ψ₂Sq(x_P) = 0` を得た。★`Ψ₂Sq` は `X` だけの多項式なので
根は有限個。★★あとは **`x` から `y` が決まる**ことが言えれば `E[2]` は有限である。

    y = negY x y  ⟹  2y = −(a₁x + a₃)

★★★標数 ≠ 2 なら `y` は一意に決まる。

## ★★除算を書かないのが要点

`y = −(a₁x + a₃)/2` と書くと `field_simp` が要り、`W.a₁` と `W.toAffine.a₁` の
**綴りの違い**で `ring` が刺さらなかった(実測)。
★**`2 * y = …` の形にすれば除算が消えて `linear_combination` 一発**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `two_torsion_two_mul_y` | ★`2y = −(a₁x + a₃)` |
| `two_torsion_y_unique` | ★★★**2-捩れの `y` は `x` で決まる** |
| `psi2Sq_ne_zero_of_two` | ★`Ψ₂Sq ≠ 0`(標数 ≠ 2) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★`y = negY x y` なら `2y = −(a₁x + a₃)`。 -/
theorem two_torsion_two_mul_y {x y : F} (hy : y = W.toAffine.negY x y) :
    2 * y = -(W.toAffine.a₁ * x + W.toAffine.a₃) := by
  rw [WeierstrassCurve.Affine.negY] at hy
  linear_combination hy

/-- ★★★**2-捩れの `y` は `x` で決まる**(標数 ≠ 2)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これで `E[2]` の有限性が `Ψ₂Sq` の根の有限性に帰着する。 -/
theorem two_torsion_y_unique (h2 : (2 : F) ≠ 0) {x y y' : F}
    (hy : y = W.toAffine.negY x y) (hy' : y' = W.toAffine.negY x y') : y = y' :=
  mul_left_cancel₀ h2 <| by
    rw [two_torsion_two_mul_y W hy, two_torsion_two_mul_y W hy']

/-- ★標数 ≠ 2 なら `Ψ₂Sq ≠ 0`。 -/
theorem psi2Sq_ne_zero_of_two (h2 : (2 : F) ≠ 0) : W.Ψ₂Sq ≠ 0 :=
  W.Ψ₂Sq_ne_zero (by
    intro h
    have hmul : (2 : F) * 2 = 0 := by linear_combination h
    rcases mul_eq_zero.1 hmul with h' | h' <;> exact h2 h')

/-- ★`Ψ₂Sq` の根は有限個。 -/
theorem finite_psi2Sq_roots (h2 : (2 : F) ≠ 0) : {x : F | W.Ψ₂Sq.IsRoot x}.Finite :=
  Polynomial.finite_setOf_isRoot (psi2Sq_ne_zero_of_two W h2)

/-! ## ★出典の紐付け(`.src`) -/

def two_torsion_y_unique.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(2-捩れの y は x で決まること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
