import ABC3.Found.GaloisRep.Doubling

/-!
# Galois (G1) 第 31 ブロック —— **★★★★★2-捩れの判定が群法則に繋がった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`n = 2` で `E[n]` と分点多項式が結ばれた

    2 • P = 0  ⟺  Ψ₂Sq(x_P) = 0

★これが `torsion_finite`(G1 の第 2 欄)の**型そのもの**である
——一般の `n` でも `n • P = 0 ⟺ ΨSqₙ(x_P) = 0` を示せばよい。

## ★★機構

| 段 | 内容 |
|---|---|
| 1 | `2 • P = 0 ⟺ P = −P`(群論) |
| 2 | `−(x,y) = (x, negY x y)`(mathlib `Point.neg_some`) |
| 3 | `Point.some.injEq` で座標の等式に落とす |
| 4 | `y = negY x y ⟺ Ψ₂Sq(x) = 0`(第 29) |

★★★**`Ψ₂Sq` が `X` だけの多項式**であることが効く——根が有限個なので
`E[2]` の有限性が直ちに出る(§9-366 で測ったとおり)。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★**2-捩れの判定**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `E[n]` と分点多項式を結ぶ形の `n = 2` の場合である。 -/
theorem two_smul_eq_zero_iff {x y : F} (h : W.toAffine.Nonsingular x y) :
    (2 : ℕ) • (Point.some _ _ h) = 0 ↔ W.Ψ₂Sq.eval x = 0 := by
  rw [psi2Sq_eval_eq_zero_iff W h.left]
  rw [two_smul, add_eq_zero_iff_eq_neg, Point.neg_some h, Point.some.injEq]
  exact ⟨fun hP => hP.2, fun hP => ⟨rfl, hP⟩⟩


/-! ## ★出典の紐付け(`.src`) -/

def two_smul_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(2-捩れの判定——E[n] と分点多項式の接続)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
