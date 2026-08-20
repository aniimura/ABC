import ABC3.Found.GaloisRep.GenericNotTorsion

/-!
# Galois (G5) 第 126 ブロック —— **★★★★★`f_P` は座標環の元である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★層 3 の第 1 段

第 113 ブロックで「`(XYIdeal' P)^n` は**分数イデアルとして**単項」を得た。
★Weil 対の構成では `f_P` を `[n]^*`(第 118・125 の `μ : F[W] →+* F(W)`)に
食わせるので、**`f_P` が座標環の元であること**が要る。

★★これは mathlib の 2 つの事実で出る:

| mathlib | 内容 |
|---|---|
| `CoordinateRing.XYIdeal'_eq` | `XYIdeal'` は**整イデアル** `XYIdeal W x (C y)` に等しい |
| `ClassGroup.mk_eq_one_of_coe_ideal` | 類が自明 ⟺ その整イデアルが単項 |

★★★したがって `f_P ∈ F[W]` が**直ちに出る**——古典的には
「`div(f_P) = n(P) − n(O)` だから `f_P` はアフィン曲線上正則」に対応する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `xyIdeal_pow_isPrincipal_integral` | ★★★★★**`(XYIdeal)^n = (f_P)` in `F[W]`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F] [DecidableEq F] {W : WeierstrassCurve.Affine F}

/-- ★★★★★**捩れ点のイデアルの冪は座標環の中で単項である**——`f_P ∈ F[W]`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 113 ブロックの分数イデアル版を、mathlib の `XYIdeal'_eq`(整イデアルとの一致)で
**整イデアルの単項性**に落としたものである。
★★これにより `f_P` を `[n]` の引き戻し `μ : F[W] →+* F(W)` に食わせられる。 -/
theorem xyIdeal_pow_isPrincipal_integral {x y : F} (h : W.Nonsingular x y) (n : ℕ)
    (hP : n • (Point.some x y h) = 0) :
    ∃ f : W.CoordinateRing, f ≠ 0 ∧
      (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {f} := by
  refine (ClassGroup.mk_eq_one_of_coe_ideal (I := (CoordinateRing.XYIdeal' h) ^ n)
    (I' := (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n) ?_).1 ?_
  · rw [FractionalIdeal.coeIdeal_pow, ← CoordinateRing.XYIdeal'_eq h, Units.val_pow_eq_pow_val]
  · rw [map_pow]
    exact classGroup_pow_eq_one h n hP

/-! ## ★出典の紐付け(`.src`) -/

def xyIdeal_pow_isPrincipal_integral.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——f_P が座標環の元であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
