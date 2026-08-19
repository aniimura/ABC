import ABC3.Found.GaloisRep.OmegaCurve
import Mathlib.AlgebraicGeometry.EllipticCurve.DivisionPolynomial.Degree

/-!
# Galois (G1) 第 29 ブロック —— **`preΨₙ ≠ 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★`E[n]` の有限性に要る土台

`torsion_finite`(G1 の第 2 欄)は「`n`-捩れ点の `x` 座標は `ψₙ` の根」から出る。
★そのためには **`ψₙ ≠ 0`** が要る——さもないと根が無限個になりうる。

## ★★mathlib に在庫があった

`DivisionPolynomial/Degree.lean` の

    coeff_preΨ_ne_zero {n : ℤ} (h : (n : R) ≠ 0) :
      (W.preΨ n).coeff ((n.natAbs ^ 2 - if Even n then 4 else 1) / 2) ≠ 0

★**仮定 `(n : R) ≠ 0` は G1 の主張の仮定そのもの**である
(`structure_eq` / `torsion_finite` がどちらも `(n : K) ≠ 0` を課している)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preΨ_ne_zero` | ★★★**`preΨₙ ≠ 0`**(`(n : R) ≠ 0` のとき) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★★★**`preΨₙ ≠ 0`**——係数の 1 つが非零だから。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `E[n]` の有限性の土台である。 -/
theorem preΨ_ne_zero {n : ℤ} (h : (n : R) ≠ 0) : W.preΨ n ≠ 0 := by
  intro hz
  exact W.coeff_preΨ_ne_zero h (by rw [hz]; simp)

/-! ## ★出典の紐付け(`.src`) -/

def preΨ_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(n-捩れの有限性——preΨₙ が非零であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
