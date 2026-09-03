/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim

/-!
# 第 1284 ブロック —— **座標は写像と可換**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——節点 8 の下ごしらえ

第 1283 で Tate 一意化の同変性を**座標の言葉**で取った。
★同じ言葉で `rhPoint`（体の準同型）と `vcPoint`（変数変換）を書いておけば、
「Tate モデルの点」と「`E` の点」の間で `σ` の作用が対応することが言える。

☆本ブロックはその 3 本である:

| 定理 | 内容 |
|---|---|
| `pointCoords_rhPoint'` | `rhPoint f` は座標に `f` を当てるだけ(原点も含む形。在庫の `pointCoords_rhPoint` は `some` の場合だけ) |
| `map_vcX`・`map_vcY` | 変数変換の座標式は環準同型と可換 |
| `map_vcX_fixed`・`map_vcY_fixed` | `σ` が `C` を固定するなら `σ` は `vcX`・`vcY` と可換 |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Meta

/-- ★★★★★★**`rhPoint f` は座標に `f` を当てるだけ**——★**無条件**（第 1284）。

☆在庫の `pointCoords_rhPoint` は `.some` の場合だけなので、原点も含む形を取る。 -/
theorem pointCoords_rhPoint' {F K : Type*} [Field F] [Field K] (f : F →+* K)
    (W : WeierstrassCurve F) (P : W.toAffine.Point) :
    pointCoords (rhPoint f W P) = (f (pointCoords P).1, f (pointCoords P).2) := by
  cases P with
  | zero =>
      show ((0 : K), (0 : K)) = (f (0 : F), f (0 : F))
      rw [map_zero]
  | some x y h => rfl

/-- ★★★★★★**`vcX` は環準同型と可換**——★**無条件**（第 1284）。 -/
theorem map_vcX {R S : Type*} [CommRing R] [CommRing S] (f : R →+* S)
    (C : VariableChange R) (x : R) :
    f (vcX C x) = vcX (C.map f) (f x) := by
  have hu : f (((C.u⁻¹ : Rˣ) : R)) = (((C.map f).u⁻¹ : Sˣ) : S) := by
    show f (((C.u⁻¹ : Rˣ) : R)) = (((Units.map (f : R →* S) C.u)⁻¹ : Sˣ) : S)
    rw [← map_inv]
    rfl
  show f (((C.u⁻¹ : Rˣ) : R) ^ 2 * (x - C.r))
    = (((C.map f).u⁻¹ : Sˣ) : S) ^ 2 * (f x - (C.map f).r)
  rw [map_mul, map_pow, map_sub, hu]
  rfl

/-- ★★★★★★**`vcY` は環準同型と可換**——★**無条件**（第 1284）。 -/
theorem map_vcY {R S : Type*} [CommRing R] [CommRing S] (f : R →+* S)
    (C : VariableChange R) (x y : R) :
    f (vcY C x y) = vcY (C.map f) (f x) (f y) := by
  have hu : f (((C.u⁻¹ : Rˣ) : R)) = (((C.map f).u⁻¹ : Sˣ) : S) := by
    show f (((C.u⁻¹ : Rˣ) : R)) = (((Units.map (f : R →* S) C.u)⁻¹ : Sˣ) : S)
    rw [← map_inv]
    rfl
  show f (((C.u⁻¹ : Rˣ) : R) ^ 3 * (y - C.s * (x - C.r) - C.t))
    = (((C.map f).u⁻¹ : Sˣ) : S) ^ 3 * (f y - (C.map f).s * (f x - (C.map f).r) - (C.map f).t)
  rw [map_mul, map_pow, map_sub, map_sub, map_mul, map_sub, hu]
  rfl

/-- ★★★★★★★★**`σ` が `C` を固定するなら `σ` は `vcX` と可換**——★**無条件**（第 1284）。 -/
theorem map_vcX_fixed {R : Type*} [CommRing R] (f : R →+* R)
    (C : VariableChange R) (hC : C.map f = C) (x : R) :
    f (vcX C x) = vcX C (f x) := by
  rw [map_vcX, hC]

/-- ★★★★★★★★**`σ` が `C` を固定するなら `σ` は `vcY` と可換**——★**無条件**（第 1284）。 -/
theorem map_vcY_fixed {R : Type*} [CommRing R] (f : R →+* R)
    (C : VariableChange R) (hC : C.map f = C) (x y : R) :
    f (vcY C x y) = vcY C (f x) (f y) := by
  rw [map_vcY, hC]

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_rhPoint'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(rhPoint は座標に f を当てるだけ。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_vcX.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(vcX は環準同型と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_vcY.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(vcY は環準同型と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_vcX_fixed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が C を固定するなら σ は vcX と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_vcY_fixed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が C を固定するなら σ は vcY と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
