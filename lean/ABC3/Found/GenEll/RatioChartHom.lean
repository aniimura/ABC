/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SpanLocalization
import ABC3.Found.GenEll.ProjPointOfCoords
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★比の組から `ℙᴺ` の点を作る —— 任意の可換環で（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★これは何か —— `§9-873` の一般化

`§9-873`（`awayHomOfCoords`）は座標の組 `a : Fin (N+1) → ℂ` と `a_i ≠ 0` から
`A⁰_{x_i} →+* ℂ` を作った。★そこでは `a_k/a_i` を作るのに**割り算**を使っていた。

★★★★本ファイルは**割り算を使わない形**に一般化する:

    **比の組** `r : Fin (N+1) → S`（`S` 任意の可換環）で `r_i = 1` なら
    `A⁰_{x_i} →+* S` が定まり、`x_k/x_i ↦ r_k`

★★これで `S = (𝓞_F)_Q`（局所化）が使える——`§9-943` が
`x_k/x_j ∈ (𝓞_F)_Q` を与えているからである。

## ★★★機構 —— `A⁰_{x_i} ≅ ℤ[x]/(x_i − 1)`

`§9-859`（段 C2a-2）の `awayQuotHom` に `Ideal.Quotient.lift` を繋ぐだけ。
★`x_i − 1` が核に入るのは **`r_i = 1`** だからで、
`§9-873` の `a_i/a_i = 1` が**仮定そのもの**に置き換わっている。

## ★これで何が言えるか

★★★`§9-943` の `j`（座標の最小割り切り成分）に対し `r_k ≔ x_k/x_j ∈ (𝓞_F)_Q` と置けば

    `Spec (𝓞_F)_Q ⟶ D₊(x_j) ⊆ ℙᴺ`

が**構成できる**。★残るのはこれが「局所化した点」と一致することであり、
道筋は生成点での一致（`§9-941`）＋ `ℙᴺ` の分離性である。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★比の組による評価 -/

/-- ★**比の組から `ℤ[x] →+* S`**（`x_k ↦ r_k`）——割り算を使わない形。 -/
noncomputable def evalRatios (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S) :
    MvPolynomial (Fin (N+1)) ℤ →+* S :=
  MvPolynomial.eval₂Hom (Int.castRingHom S) r

theorem evalRatios_X (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S) (k : Fin (N+1)) :
    evalRatios N S r (MvPolynomial.X k) = r k := by
  rw [evalRatios, MvPolynomial.eval₂Hom_X']

/-- ★**`x_i − 1` は核に入る** —— ここで `r_i = 1` を使う
（`§9-873` の `a_i/a_i = 1` が仮定そのものになった）。 -/
theorem span_le_ker_evalRatios (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S)
    (i : Fin (N+1)) (hi : r i = 1) :
    Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) ℤ)}
      ≤ RingHom.ker (evalRatios N S r) := by
  rw [Ideal.span_le]
  rintro z rfl
  rw [SetLike.mem_coe, RingHom.mem_ker, map_sub, map_one, evalRatios_X, hi, sub_self]

/-! ## ★★★★★★★★★★★★★任意の可換環へのチャート射 -/

/-- ★★★★★★★★★★★★★**比の組が定める `A⁰_{x_i} →+* S`**（`S` は任意の可換環）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-873` は `S = ℂ` で `a_k/a_i` を割り算で作っていたが、
比を**直接与える**ことにすれば任意の可換環で動く。
★★これで `S = (𝓞_F)_Q` が使える（`§9-943` が `x_k/x_j ∈ (𝓞_F)_Q` を与える）。 -/
noncomputable def awayHomOfRatios (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S)
    (i : Fin (N+1)) (hi : r i = 1) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
      (MvPolynomial.X i) →+* S :=
  (Ideal.Quotient.lift _ (evalRatios N S r) (fun _ hz =>
    RingHom.mem_ker.1 (span_le_ker_evalRatios N S r i hi hz))).comp (awayQuotHom N ℤ i)

/-- ★★★**正規化座標は `r_k` へ行く**。 -/
theorem awayHomOfRatios_projCoord (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S)
    (i : Fin (N+1)) (hi : r i = 1) (k : Fin (N+1)) :
    awayHomOfRatios N S r i hi (projCoord N ℤ i k) = r k := by
  have hq : awayQuotHom N ℤ i (projCoord N ℤ i k)
      = Ideal.Quotient.mk _ (MvPolynomial.X k) := by
    have h := congrArg (fun f : MvPolynomial (Fin (N+1)) ℤ →+*
        (MvPolynomial (Fin (N+1)) ℤ ⧸ Ideal.span
          {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) ℤ)}) => f (MvPolynomial.X k))
      (awayQuotHom_comp_awayEval N ℤ i)
    simp only [RingHom.comp_apply] at h
    rw [show awayEval N ℤ i (MvPolynomial.X k) = projCoord N ℤ i k from by
      rw [awayEval, MvPolynomial.eval₂Hom_X']] at h
    exact h
  rw [awayHomOfRatios, RingHom.comp_apply, hq, Ideal.Quotient.lift_mk, evalRatios_X]

/-! ## ★★★★★★★★比の組が定める `ℙᴺ` の点 -/

/-- ★★★★★★★★**比の組が定める `Spec S ⟶ ℙᴺ`**（チャート `D₊(x_i)` を通る）。 -/
noncomputable def projPointOfRatios (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S)
    (i : Fin (N+1)) (hi : r i = 1) :
    Spec (CommRingCat.of S) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) :=
  Spec.map (CommRingCat.ofHom (awayHomOfRatios N S r i hi)) ≫ chartA N ℤ i

/-- ★★**構成した点は定義からチャートを通る**。 -/
theorem range_projPointOfRatios (N : ℕ) (S : Type) [CommRing S] (r : Fin (N+1) → S)
    (i : Fin (N+1)) (hi : r i = 1) :
    Set.range (projPointOfRatios N S r i hi).base ⊆ Set.range (chartA N ℤ i).base := by
  rintro _ ⟨y, rfl⟩
  exact ⟨(Spec.map (CommRingCat.ofHom (awayHomOfRatios N S r i hi))).base y, rfl⟩

/-! ## ★出典の紐付け(`.src`) -/

def evalRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(比の組から ℤ[x] →+* S——割り算を使わない形)",
    sectionId := "genell-prop-1-4" }

def span_le_ker_evalRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(x_i − 1 は核に入る——r_i = 1 から)",
    sectionId := "genell-prop-1-4" }

def awayHomOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(比の組が定める A⁰_{x_i} →+* S——任意の可換環で)",
    sectionId := "genell-prop-1-4" }

def awayHomOfRatios_projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(正規化座標は r_k へ行く)",
    sectionId := "genell-prop-1-4" }

def projPointOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(比の組が定める Spec S ⟶ ℙᴺ)",
    sectionId := "genell-prop-1-4" }

def range_projPointOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(構成した点は定義からチャートを通る)",
    sectionId := "genell-prop-1-4" }

def awayHomOfRatios.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayQuotHom(A⁰_{x_i} ≅ ℤ[x]/(x_i − 1)、段 C2a-2、§9-859)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayQuotHom") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-873(awayHomOfCoords)が ℂ 専用だったのは" ++
       "**a_k/a_i を割り算で作っていたから**であり、比を直接与えることにすれば" ++
       "任意の可換環で動く。x_i − 1 が核に入るのは r_i = 1 という**仮定そのもの**である") 4,
    .implicitStep
      ("★これで S = (𝓞_F)_Q が使える(§9-943 が x_k/x_j ∈ (𝓞_F)_Q を与える)。" ++
       "★★残るのは構成した点が「局所化した点」と一致することであり、" ++
       "道筋は生成点での一致(§9-941)＋ ℙᴺ の分離性である") 4 ]

end ABC3.Found.GenEll
