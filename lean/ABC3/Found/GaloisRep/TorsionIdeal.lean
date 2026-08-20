import ABC3.Found.GaloisRep.GalRepWitness
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point

/-!
# Galois (G5) 第 113 ブロック —— **★★★★★捩れ点のイデアルは単項**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★Weil 対の第一段は mathlib から出る

Weil 対 `e_n` の古典的な構成は、`P ∈ E[n]` に対して

    div(f_P) = n(P) − n(O)

なる関数 `f_P` を取ることから始まる。★これは「因子 `n(P) − n(O)` が主因子である」
という主張であり、**イデアル類群の言葉に直せば mathlib にある**:

| mathlib | 内容 |
|---|---|
| `Point.toClass` | `W.Point →+ Additive (ClassGroup W.CoordinateRing)`、**単射** |
| `Point.toClass_some` | `= ClassGroup.mk (XYIdeal' h)` |
| `ClassGroup.mk_eq_one_iff` | `mk I = 1 ↔ I が単項` |

★★したがって `n • P = 0` から `(XYIdeal' P)^n` が単項であることが**直ちに出る**
——その生成元が `f_P` である。

## ★★これは Skeleton ではなく Found である

2026-08-20 の在庫調査で判明した。★当初は「因子の層から積む」と見積もっていたが、
**mathlib は楕円曲線の群法則を類群経由で証明している**ので、
`toClass` がそのまま使える。★★見積もりの下方修正を `.src` とともに記録する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `classGroup_pow_eq_one` | ★★★`n • P = 0` ⟹ 類の `n` 乗が自明 |
| `xyIdeal_pow_isPrincipal` | ★★★★★**`(XYIdeal' P)^n` は単項**(= `f_P` の存在) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve.Affine

variable {F : Type} [Field F] [DecidableEq F] {W : WeierstrassCurve.Affine F}

/-- ★★★**`n • P = 0` なら `P` のイデアル類の `n` 乗は自明**。

★`toClass` が加法準同型であることの直接の帰結である。 -/
theorem classGroup_pow_eq_one {x y : F} (h : W.Nonsingular x y) (n : ℕ)
    (hP : n • (Point.some x y h) = 0) :
    (ClassGroup.mk W.FunctionField (CoordinateRing.XYIdeal' (W := W) h)) ^ n = 1 := by
  have h1 : Point.toClass (Point.some x y h)
      = Additive.ofMul (ClassGroup.mk W.FunctionField (CoordinateRing.XYIdeal' (W := W) h)) :=
    Point.toClass_some h
  have h2 : (n : ℕ) • Point.toClass (Point.some x y h) = 0 := by
    rw [← map_nsmul, hP, map_zero]
  rw [h1, ← ofMul_pow] at h2
  exact ofMul_eq_zero.1 h2

/-- ★★★★★**捩れ点のイデアルの冪は単項である**——これが `f_P` の存在である。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★古典的には `div(f_P) = n(P) − n(O)` と書かれるものを、
イデアル類群の言葉で述べたものである。★★生成元が `f_P` にあたる。 -/
theorem xyIdeal_pow_isPrincipal {x y : F} (h : W.Nonsingular x y) (n : ℕ)
    (hP : n • (Point.some x y h) = 0) :
    ((((CoordinateRing.XYIdeal' h) ^ n : (FractionalIdeal (nonZeroDivisors W.CoordinateRing)
        W.FunctionField)ˣ) : FractionalIdeal (nonZeroDivisors W.CoordinateRing)
        W.FunctionField) : Submodule W.CoordinateRing W.FunctionField).IsPrincipal := by
  rw [← ClassGroup.mk_eq_one_iff (K := W.FunctionField), map_pow]
  exact classGroup_pow_eq_one h n hP

/-! ## ★出典の紐付け(`.src`) -/

def xyIdeal_pow_isPrincipal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成の第一段——f_P の存在)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
