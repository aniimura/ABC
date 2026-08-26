import ABC3.Found.GenEll.LatticeCurve
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point

/-!
# GenEll 第 336 ブロック —— **★★★★★★`(℘, ℘'/2)` は曲線の点である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★一意化の (ii) の前半

一意化の連鎖の (ii)「`ℂ/Λ ≅ E(ℂ)` の群同型」のうち、**写像を作る側**を取る:

    z ∉ Λ  ⟹  (℘(z), ℘'(z)/2) は latticeCurve P の点

★これは `℘'² = 4℘³ − g₂℘ − g₃`(mathlib の `derivWeierstrassP_sq`)を
`y² = x³ − (g₂/4)x − (g₃/4)` に読み替えるだけで、**`linear_combination h/4` で出る**。
★★非特異性は楕円曲線なので `equation_iff_nonsingular` から自動。

★★★さらに **`℘` が偶・`℘'` が奇**であることから、点の対応は**符号を保つ**:

    latticePoint P (-z) = -(latticePoint P z)

★`latticeCurve P` は `a₁ = a₃ = 0` なので `negY x y = -y` であり、そのまま合う。

## ★★残り((ii) の後半)

★加法定理——`latticePoint P (z + w) = latticePoint P z + latticePoint P w`。
★★mathlib に `℘` の加法定理は無い(2026-08-26 実測、`weierstrassP_add` で 0 件)。
★★★古典的には「`℘, ℘'` の生成する体は `Λ` 周期の楕円関数体である」ことか、
あるいは 3 点共線条件を直接示す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `equation_weierstrassP` | ★★★★★**`(℘, ℘'/2)` は方程式を満たす** |
| `nonsingular_weierstrassP` | ★★★非特異 |
| `latticePoint` | ★★★★**曲線の点としての `z`** |
| `latticePoint_neg` | ★★★★**`z ↦ -z` は点の符号に対応** |
-/

namespace ABC3.Found.GenEll

open PeriodPair WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★方程式を満たすこと -/

/-- ★★★★★**`(℘(z), ℘'(z)/2)` は `latticeCurve P` の方程式を満たす**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`℘'² = 4℘³ − g₂℘ − g₃` を `y² = x³ − (g₂/4)x − (g₃/4)` に読み替えるだけ。 -/
theorem equation_weierstrassP (P : PeriodPair) (z : ℂ) (hz : z ∉ P.lattice) :
    (latticeCurve P).toAffine.Equation (P.weierstrassP z) (P.derivWeierstrassP z / 2) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  have h := P.derivWeierstrassP_sq z hz
  simp only [latticeCurve]
  linear_combination h / 4

/-- ★★★非特異である(楕円曲線だから)。 -/
theorem nonsingular_weierstrassP (P : PeriodPair) [(latticeCurve P).IsElliptic]
    (z : ℂ) (hz : z ∉ P.lattice) :
    (latticeCurve P).toAffine.Nonsingular (P.weierstrassP z) (P.derivWeierstrassP z / 2) :=
  (WeierstrassCurve.Affine.equation_iff_nonsingular).1 (equation_weierstrassP P z hz)

/-! ## ★★★★曲線の点としての `z` -/

/-- ★★★★**`z ∉ Λ` に対応する曲線の点**。 -/
noncomputable def latticePoint (P : PeriodPair) [(latticeCurve P).IsElliptic]
    (z : ℂ) (hz : z ∉ P.lattice) : (latticeCurve P).toAffine.Point :=
  Point.some (P.weierstrassP z) (P.derivWeierstrassP z / 2) (nonsingular_weierstrassP P z hz)

/-- ★★★★**`z ↦ -z` は点の符号に対応する**——`℘` が偶、`℘'` が奇だから。 -/
theorem latticePoint_neg (P : PeriodPair) [(latticeCurve P).IsElliptic]
    (z : ℂ) (hz : z ∉ P.lattice) (hz' : -z ∉ P.lattice) :
    latticePoint P (-z) hz' = -latticePoint P z hz := by
  rw [latticePoint, latticePoint, Point.neg_some]
  congr 1
  · exact P.weierstrassP_neg z
  · rw [P.derivWeierstrassP_neg]
    show -P.derivWeierstrassP z / 2
      = (latticeCurve P).toAffine.negY (P.weierstrassP z) (P.derivWeierstrassP z / 2)
    simp only [WeierstrassCurve.Affine.negY, latticeCurve]
    ring

/-! ## ★出典の紐付け(`.src`) -/

def equation_weierstrassP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def latticePoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
