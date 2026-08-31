/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.OverlapAwayHom
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★重なりの環準同型は `i`・`j` について対称（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★これは何か —— 段 E1d の心臓

`§9-909` は `A⁰_{x_i x_j} →+* Γ(X, W)` を **`i` を主に取って**作った。
★同じものを **`j` を主に取って**も作れる。★★**その 2 つが等しい**というのが、
「2 つのチャート射が重なりで一致する」ことの代数的な中身である。

## ★★★機構 —— 局所化の一意性 ＋ 生成元での確認 ＋ コサイクル則

★`A⁰_{x_i x_j}` は `A⁰_{x_i}` の局所化なので、
そこからの 2 本の環準同型は **`A⁰_{x_i}` への制限が一致すれば等しい**
（`IsLocalization.ringHom_ext`）。

★★制限の一致は `ext_of_projCoord`（`A⁰_{x_i}` からの環準同型は
定数と正規化座標で決まる、`§9-AwayGenerated`）に落ちる:

| 元 | 要る等式 |
|---|---|
| 定数 `c` | `awayMap_j(定数) = awayMap_i(定数)`（本ファイル） |
| `x_k/x_i` | `awayMap_j(x_k/x_i)·awayMap_i(x_i/x_j) = awayMap_i(x_k/x_j)`（本ファイル） |

★★★後者に `j` 側の特徴づけを当てると

    `χ(x_k/x_i) · (s_i/s_j) = s_k/s_j`

になり、**コサイクル則**（`§9-908`）`s_k/s_j = (s_k/s_i)·(s_i/s_j)` と
`s_i/s_j` が単元であること（`§9-907`）から `χ(x_k/x_i) = s_k/s_i` が出る。
★これは `i` 側の特徴づけそのものである。

## ★`x` を等式で受ける形に一般化した

`§9-909` の `overlapAwayHom` は始域が `A⁰_{x_i·x_j}` に固定されていた。
★`j` を主に取った版の始域は `A⁰_{x_j·x_i}` になってしまい**型が違う**ので、
`hx : x = x_i·x_j` を引数で受ける `overlapAwayHomOf` に一般化した。
★★これで両方の始域を `A⁰_{x_i·x_j}` に揃えられる（`hx` を `mul_comm` に取る）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★始域を等式で受ける形 -/

/-- ★★★★★★★★★★★★**重なりの環準同型**（`x = x_i·x_j` を等式で受ける形）。

★`§9-909` の `overlapAwayHom` は `x` が `x_i·x_j` に固定されていた。
`j` を主に取った版と型を揃えるため、`hx` を引数にした。 -/
noncomputable def overlapAwayHomOf (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    {x : MvPolynomial (Fin (N + 1)) R} (hx : x = MvPolynomial.X i * MvPolynomial.X j)
    {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) x
      →+* (Γ(X, W) : Type) :=
  letI := (HomogeneousLocalization.awayMap
    (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (MvPolynomial.isHomogeneous_X R j) hx).toAlgebra
  haveI := HomogeneousLocalization.Away.isLocalization_mul
    (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (f := MvPolynomial.X i) (g := MvPolynomial.X j) (d := 1) (e := 1)
    (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j) hx
    one_ne_zero
  IsLocalization.Away.lift _ (isUnit_globalAwayHom_isLocalizationElem M hM φ s i j hi hj)

/-- ★★★**`§9-909` の形は `hx := rfl` の場合である**。 -/
theorem overlapAwayHom_eq_of (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    overlapAwayHom M hM φ s i j hi hj = overlapAwayHomOf M hM φ s i j rfl hi hj := rfl

/-- ★★★★★★**主に取った側の特徴づけ**。 -/
theorem overlapAwayHomOf_awayMap (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1))
    {x : MvPolynomial (Fin (N + 1)) R} (hx : x = MvPolynomial.X i * MvPolynomial.X j)
    {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j))
    (z : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) :
    overlapAwayHomOf M hM φ s i j hx hi hj (HomogeneousLocalization.awayMap
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R j) hx z)
      = X.presheaf.map (homOfLE hi).op (globalAwayHom M hM φ s i z) := by
  letI := (HomogeneousLocalization.awayMap
    (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (MvPolynomial.isHomogeneous_X R j) hx).toAlgebra
  haveI := HomogeneousLocalization.Away.isLocalization_mul
    (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (f := MvPolynomial.X i) (g := MvPolynomial.X j) (d := 1) (e := 1)
    (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j) hx
    one_ne_zero
  exact IsLocalization.Away.lift_eq _
    (isUnit_globalAwayHom_isLocalizationElem M hM φ s i j hi hj) z

/-! ## ★★★★★2 つの `awayMap` の間の代数的な関係 -/

/-- ★★★★★**`awayMap_j(x_k/x_i)·awayMap_i(x_i/x_j) = awayMap_i(x_k/x_j)`**。

★`A⁰_{x_i x_j}` の中での分数の計算そのもの:
`(x_k x_j)/(x_i x_j) · (x_i x_i)/(x_i x_j) = (x_k x_i)/(x_i x_j)`。 -/
theorem awayMap_projCoord_mul (N : ℕ) (R : Type) [CommRing R] (i j k : Fin (N + 1)) :
    HomogeneousLocalization.awayMap (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R j)
        (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
          = MvPolynomial.X i * MvPolynomial.X j) (projCoord N R i k)
      * HomogeneousLocalization.awayMap (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R i)
        (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j))
        (projCoord N R j i)
      = HomogeneousLocalization.awayMap (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R i)
        (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j))
        (projCoord N R j k) := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [HomogeneousLocalization.val_mul]
  rw [projCoord, projCoord, projCoord, awayMap_mk, awayMap_mk, awayMap_mk,
    Away.val_mk, Away.val_mk, Away.val_mk]
  rw [Localization.mk_mul, Localization.mk_eq_mk_iff, Localization.r_iff_exists]
  refine ⟨1, ?_⟩
  simp only [OneMemClass.coe_one, one_mul, Submonoid.coe_mul, pow_one]
  ring

/-- ★★★**定数はどちらの `awayMap` でも同じ所に行く**。 -/
theorem awayMap_awayConst_comm (N : ℕ) (R : Type) [CommRing R] (i j : Fin (N + 1)) (c : R) :
    HomogeneousLocalization.awayMap (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R j)
        (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
          = MvPolynomial.X i * MvPolynomial.X j) (awayConst N R i c)
      = HomogeneousLocalization.awayMap (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R i)
        (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j))
        (awayConst N R j c) := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [awayConst, awayConst, awayMap_mk, awayMap_mk, Away.val_mk, Away.val_mk]
  simp

/-! ## ★★★★★★★★★★★★★★対称性 -/

/-- ★★★★★★★★★★★★★★**重なりの環準同型は `i`・`j` について対称である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これが段 E1d の**代数的な中身**である——`i` を主に取って作ったものと
`j` を主に取って作ったものが**同じ写像**である。

★機構は 3 段:
1. `A⁰_{x_i x_j}` は `A⁰_{x_i}` の局所化なので `IsLocalization.ringHom_ext` で
   `A⁰_{x_i}` への制限の一致に落とす。
2. `ext_of_projCoord` で定数と正規化座標での一致に落とす。
3. 正規化座標では `awayMap_projCoord_mul` と**コサイクル則**（`§9-908`）を使い、
   `s_i/s_j` が単元であること（`§9-907`）で割る。 -/
theorem overlapAwayHomOf_symm (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    overlapAwayHomOf M hM φ s i j rfl hi hj
      = overlapAwayHomOf M hM φ s j i
          (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R)
            (MvPolynomial.X j)) hj hi := by
  letI := (HomogeneousLocalization.awayMap
    (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (MvPolynomial.isHomogeneous_X R j)
    (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
      = MvPolynomial.X i * MvPolynomial.X j)).toAlgebra
  haveI := HomogeneousLocalization.Away.isLocalization_mul
    (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (f := MvPolynomial.X i) (g := MvPolynomial.X j) (d := 1) (e := 1)
    (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j)
    (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
      = MvPolynomial.X i * MvPolynomial.X j)
    one_ne_zero
  refine IsLocalization.ringHom_ext
    (Submonoid.powers (HomogeneousLocalization.Away.isLocalizationElem
      (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j))) ?_
  refine ext_of_projCoord i _ _ ?_ ?_
  · intro c
    show overlapAwayHomOf M hM φ s i j rfl hi hj
        (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R j) rfl
          (awayConst N R i c))
      = overlapAwayHomOf M hM φ s j i (mul_comm _ _) hj hi
        (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R j) rfl
          (awayConst N R i c))
    rw [overlapAwayHomOf_awayMap, awayMap_awayConst_comm, overlapAwayHomOf_awayMap,
      globalAwayHom_awayConst, globalAwayHom_awayConst, res_trans, res_trans]
  · intro k
    show overlapAwayHomOf M hM φ s i j rfl hi hj
        (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R j) rfl
          (projCoord N R i k))
      = overlapAwayHomOf M hM φ s j i (mul_comm _ _) hj hi
        (HomogeneousLocalization.awayMap _ (MvPolynomial.isHomogeneous_X R j) rfl
          (projCoord N R i k))
    rw [overlapAwayHomOf_awayMap, globalAwayHom_projCoord]
    have hkey := congrArg (overlapAwayHomOf M hM φ s j i
      (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j)) hj hi)
      (awayMap_projCoord_mul N R i j k)
    rw [map_mul, overlapAwayHomOf_awayMap M hM φ s j i _ hj hi,
      overlapAwayHomOf_awayMap M hM φ s j i _ hj hi,
      globalAwayHom_projCoord, globalAwayHom_projCoord] at hkey
    have hco := globalRatio_cocycle_of_le M hM (s k) (s i) (s j) hi hj
    have hu := isUnit_globalRatio_of_le M hM (s i) (s j) hj hi
    refine (hu.mul_right_cancel ?_).symm
    rw [hkey, ← hco]

/-- ★★★★★★★★★★★★★★★**もう一方の側の特徴づけ** ——
`overlapAwayHom ∘ awayMap(x_i 側) = res ∘ globalAwayHom j`。

★`overlapAwayHomOf_symm` と `overlapAwayHomOf_awayMap` を繋いだだけである。
★★これで `Proj.SpecMap_awayMap_awayι` に**両側から**渡せる。 -/
theorem overlapAwayHom_awayMap_symm (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j))
    (z : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X j)) :
    overlapAwayHom M hM φ s i j hi hj (HomogeneousLocalization.awayMap
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R i)
        (mul_comm (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) (MvPolynomial.X j)) z)
      = X.presheaf.map (homOfLE hj).op (globalAwayHom M hM φ s j z) := by
  rw [overlapAwayHom_eq_of, overlapAwayHomOf_symm]
  exact overlapAwayHomOf_awayMap M hM φ s j i _ hj hi z

/-! ## ★出典の紐付け(`.src`) -/

def overlapAwayHomOf_symm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型は i・j について対称である)",
    sectionId := "genell-prop-1-4" }

def overlapAwayHom_awayMap_symm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型——もう一方の側の特徴づけ)",
    sectionId := "genell-prop-1-4" }

def overlapAwayHomOf_symm.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsLocalization.ringHom_ext(局所化からの環準同型は制限で決まる)"
      (.inMathlib "IsLocalization.ringHom_ext") 1,
    .citation "[ABC3]" "ext_of_projCoord(A⁰_{x_i} からの環準同型は定数と正規化座標で決まる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_projCoord") 2,
    .citation "[ABC3]" "globalRatio_cocycle_of_le(コサイクル則、§9-908)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_cocycle_of_le") 2,
    .citation "[ABC3]" "isUnit_globalRatio_of_le(比は W の上で単元、§9-907)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isUnit_globalRatio_of_le") 2,
    .implicitStep
      ("★★これが段 E1d の**代数的な中身**である" ++
       "——i を主に取って作った重なりの環準同型と j を主に取って作ったものが" ++
       "**同じ写像**である") 4,
    .implicitStep
      ("★残るのは幾何の側の配管: Proj.SpecMap_awayMap_awayι で" ++
       "2 つのチャート射を同じ awayι 𝒜 (x_i·x_j) へ落とし、glueOpens(§9-836)に渡す") 4 ]

def overlapAwayHomOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型——始域を等式で受ける形)",
    sectionId := "genell-prop-1-4" }

def overlapAwayHomOf_awayMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型——主に取った側の特徴づけ)",
    sectionId := "genell-prop-1-4" }

def awayMap_projCoord_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(2 つの awayMap の間の関係——正規化座標)",
    sectionId := "genell-prop-1-4" }

def awayMap_awayConst_comm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(2 つの awayMap の間の関係——定数)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
