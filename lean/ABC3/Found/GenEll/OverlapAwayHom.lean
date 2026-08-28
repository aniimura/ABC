/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalRatioCocycle
import ABC3.Found.GenEll.AwayGenerated
import ABC3.Found.GenEll.GlobalChartSurjective
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★重なりの環準同型 `A⁰_{x_i x_j} →+* Γ(X, W)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★これは何か —— 段 E1d の配管 (1)

段 E1d は「2 つのチャート射 `X_{s_i} ⟶ ℙᴺ_R` が重なりで一致する」ことである。
★★mathlib の `Proj.SpecMap_awayMap_awayι`

    `Spec.map (awayMap 𝒜 hg hx) ≫ Proj.awayι 𝒜 f = Proj.awayι 𝒜 (f·g)`

を使うと、**両方のチャート射を同じ `awayι 𝒜 (x_i·x_j)` へ落とせる**。
★そのために要るのが **`A⁰_{x_i x_j} →+* Γ(X, W)`**（`W ⊆ X_{s_i} ∩ X_{s_j}`）である。

## ★★★機構 —— 局所化の普遍性

mathlib の `HomogeneousLocalization.Away.isLocalization_mul`:

    `A⁰_{f g}` は `A⁰_f` を `g^{deg f}/f^{deg g}` で局所化したもの

★`f = x_i`、`g = x_j`（どちらも 1 次）なら、その元は `x_j/x_i = projCoord i j` である
（`isLocalizationElem_eq_projCoord'`）。
★★`globalAwayHom … (x_j/x_i) = s_j/s_i`（`§9-842`）なので、
普遍性に渡す条件は **`W` の上で `s_j/s_i` が単元**——`§9-907` が取ってある。

## ★★一般の `W` に直した

`§9-907`・`§9-908` は `W = X_t ⊓ X_s` の形で取っていたが、
★**`W ≤ X_t` と `W ≤ X_s` だけがあればよい**ので一般形に直した
（`globalRatio_mul_of_le`・`isUnit_globalRatio_of_le`・`globalRatio_cocycle_of_le`）。
★★これで `i` と `j` を入れ替えた形も**同じ `W` の上で**言える
——`inf_comm` の書き換えが要らなくなる。

## ★残っている段（明示）

★★★**対称な等式** `overlapAwayHom ∘ awayMap(x_i 側) = res ∘ globalAwayHom j` が残る。
★これは `ext_of_projCoord`（`A⁰_{x_i}` からの環準同型は定数と正規化座標で決まる）と
`globalRatio_cocycle`（`s_k/s_j = (s_k/s_i)·(s_i/s_j)`）で出る見込みである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★一般の `W` に直した比の恒等式 -/

/-- ★★★**`(s/t)·(t/s) = 1`**（`W ≤ X_t`・`W ≤ X_s` だけを仮定した形）。 -/
theorem globalRatio_mul_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) {W : X.Opens}
    (h1 : W ≤ nonVanishing M t) (h2 : W ≤ nonVanishing M s) :
    X.presheaf.map (homOfLE h1).op (globalRatio M hM s t)
      * X.presheaf.map (homOfLE h2).op (globalRatio M hM t s) = 1 := by
  have hle : W ≤ nonVanishing M t ⊓ nonVanishing M s := le_inf h1 h2
  have key := congrArg (X.presheaf.map (homOfLE hle).op) (globalRatio_mul_globalRatio M hM s t)
  rw [map_mul, map_one, res_trans, res_trans] at key
  exact key

/-- ★★★★**比は `W` の上で単元**（一般形）。 -/
theorem isUnit_globalRatio_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t : (M.obj (op ⊤) : Type)) {W : X.Opens}
    (h1 : W ≤ nonVanishing M t) (h2 : W ≤ nonVanishing M s) :
    IsUnit (X.presheaf.map (homOfLE h1).op (globalRatio M hM s t)) :=
  IsUnit.of_mul_eq_one _ (globalRatio_mul_of_le M hM s t h1 h2)

/-- ★★★★**コサイクル則 `s/u = (s/t)·(t/u)`**（一般形）。 -/
theorem globalRatio_cocycle_of_le (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (s t u : (M.obj (op ⊤) : Type)) {W : X.Opens}
    (h1 : W ≤ nonVanishing M t) (h2 : W ≤ nonVanishing M u) :
    X.presheaf.map (homOfLE h2).op (globalRatio M hM s u)
      = X.presheaf.map (homOfLE h1).op (globalRatio M hM s t)
        * X.presheaf.map (homOfLE h2).op (globalRatio M hM t u) := by
  have hle : W ≤ nonVanishing M t ⊓ nonVanishing M u := le_inf h1 h2
  have key := congrArg (X.presheaf.map (homOfLE hle).op) (globalRatio_cocycle M hM s t u)
  rw [map_mul, res_trans, res_trans, res_trans] at key
  exact key

/-! ## ★局所化の元は正規化座標である -/

/-- ★**`isLocalizationElem` は正規化座標である**（一般の係数環 `R`）。

★`§9-871` の同名の補題は `R = ℤ` に固定されていた。本補題はそれを一般化したものである。 -/
theorem isLocalizationElem_eq_projCoord' (N : ℕ) (R : Type) [CommRing R] (i j : Fin (N + 1)) :
    HomogeneousLocalization.Away.isLocalizationElem
        (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j)
      = projCoord N R i j := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [HomogeneousLocalization.Away.isLocalizationElem, projCoord,
    Away.val_mk, Away.val_mk, pow_one]

/-! ## ★★★★★局所化の普遍性に渡す条件 -/

/-- ★★★★★**`x_j/x_i` の像は `W` の上で単元である**。

★`isLocalizationElem_eq_projCoord'` と `globalAwayHom_projCoord`（`§9-842`）で
`s_j/s_i` に直し、`isUnit_globalRatio_of_le`（`§9-907`）を当てるだけである。 -/
theorem isUnit_globalAwayHom_isLocalizationElem (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    IsUnit (((X.presheaf.map (homOfLE hi).op).hom.comp (globalAwayHom M hM φ s i))
      (HomogeneousLocalization.Away.isLocalizationElem
        (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j))) := by
  show IsUnit (X.presheaf.map (homOfLE hi).op
    (globalAwayHom M hM φ s i (HomogeneousLocalization.Away.isLocalizationElem
      (MvPolynomial.isHomogeneous_X R i) (MvPolynomial.isHomogeneous_X R j))))
  rw [isLocalizationElem_eq_projCoord', globalAwayHom_projCoord]
  exact isUnit_globalRatio_of_le M hM (s j) (s i) hi hj

/-! ## ★★★★★★★★★★★★★重なりの環準同型 -/

/-- ★★★★★★★★★★★★★**重なりの環準同型** `A⁰_{x_i x_j} →+* Γ(X, W)`。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★`A⁰_{x_i x_j}` は `A⁰_{x_i}` を `x_j/x_i` で局所化したもの
（mathlib の `HomogeneousLocalization.Away.isLocalization_mul`）なので、
`x_j/x_i` の像が単元であれば **`A⁰_{x_i} → Γ(X, W)` が一意に延びる**。 -/
noncomputable def overlapAwayHom (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j)) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i * MvPolynomial.X j) →+* (Γ(X, W) : Type) :=
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
  IsLocalization.Away.lift _ (isUnit_globalAwayHom_isLocalizationElem M hM φ s i j hi hj)

/-- ★★★★★★★★★★★**`i` の側の特徴づけ** ——
`overlapAwayHom ∘ awayMap = res ∘ globalAwayHom i`。

★これが `Proj.SpecMap_awayMap_awayι` に渡す等式である。 -/
theorem overlapAwayHom_awayMap (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) {W : X.Opens}
    (hi : W ≤ nonVanishing M (s i)) (hj : W ≤ nonVanishing M (s j))
    (z : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i)) :
    overlapAwayHom M hM φ s i j hi hj (HomogeneousLocalization.awayMap
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.isHomogeneous_X R j)
        (rfl : (MvPolynomial.X i * MvPolynomial.X j : MvPolynomial (Fin (N + 1)) R)
          = MvPolynomial.X i * MvPolynomial.X j) z)
      = X.presheaf.map (homOfLE hi).op (globalAwayHom M hM φ s i z) := by
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
  exact IsLocalization.Away.lift_eq _
    (isUnit_globalAwayHom_isLocalizationElem M hM φ s i j hi hj) z

/-! ## ★出典の紐付け(`.src`) -/

def globalRatio_mul_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)((s/t)·(t/s) = 1——一般の W)",
    sectionId := "genell-prop-1-4" }

def isUnit_globalRatio_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(比は W の上で単元——一般の W)",
    sectionId := "genell-prop-1-4" }

def globalRatio_cocycle_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(コサイクル則——一般の W)",
    sectionId := "genell-prop-1-4" }

def isLocalizationElem_eq_projCoord'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(isLocalizationElem は正規化座標である——一般の係数環)",
    sectionId := "genell-prop-1-4" }

def isUnit_globalAwayHom_isLocalizationElem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(x_j/x_i の像は W の上で単元である)",
    sectionId := "genell-prop-1-4" }

def overlapAwayHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型 A⁰_{x_i x_j} →+* Γ(X, W))",
    sectionId := "genell-prop-1-4" }

def overlapAwayHom_awayMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの環準同型の i 側の特徴づけ)",
    sectionId := "genell-prop-1-4" }

def overlapAwayHom.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "HomogeneousLocalization.Away.isLocalization_mul(A⁰_{fg} は A⁰_f の局所化)"
      (.inMathlib "HomogeneousLocalization.Away.isLocalization_mul") 1,
    .citation "[ABC3]" "isUnit_globalRatio_of_le(比は W の上で単元、§9-907)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isUnit_globalRatio_of_le") 2,
    .citation "[ABC3]" "globalAwayHom_projCoord(x_j/x_i ↦ s_j/s_i、§9-842)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalAwayHom_projCoord") 1,
    .implicitStep
      ("★★一般の W に直した(W ≤ X_{s_i} と W ≤ X_{s_j} だけを仮定)。" ++
       "これで i と j を入れ替えた形も**同じ W の上で**言える" ++
       "——inf_comm の書き換えが要らなくなる") 2,
    .implicitStep
      ("★★★**対称な等式** overlapAwayHom ∘ awayMap(x_i 側) = res ∘ globalAwayHom j が残る。" ++
       "ext_of_projCoord(A⁰_{x_i} からの環準同型は定数と正規化座標で決まる)と " ++
       "globalRatio_cocycle(s_k/s_j = (s_k/s_i)·(s_i/s_j))で出る見込みである") 5 ]

end ABC3.Found.GenEll
