/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.FiniteCompat
import ABC3.Found.GenEll.HomogeneousCoordsArch
import ABC3.Found.GenEll.NorthcottPacked
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★同次座標だけで Northcott（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 素点ごとの整合が消えた

`§9-940`（束ねた形）は点ごとに**存在命題 2 本**（有限素点の整合・無限素点の整合）を
仮定として受けていた。★★★★どちらも定理になった:

* 無限素点——`§9-942`（`harch_of_homogeneous_coords`）
* 有限素点——`§9-951`（`hfin_of_homogeneous_coords`）

★★★★★★**本ファイルはそれを差し込む**。残る点の側のデータは

    点 `xF p`、その**同次座標** `x p`（`§9-941` で必ず取れる）、
    そして `hidx`・`hemb`・`hpt`（点の族が相異なる複素点を与えること、`§9-936`）

だけである。

## ★★★これが `Proposition 1.4, (iv)` の現在地

| 側 | 仮定 | 状態 |
|---|---|---|
| 層 | `IsLocallyTrivial`・`IsAmple` | ★原文の仮定そのもの |
| 幾何 | 自明化アフィン開被覆 | ★`§9-938` で消えた |
| 点 | 大域チャート | ★`§9-930`〜`935`＋`§9-939` で消えた |
| 点 | 素点ごとの整合 | ★★**本ファイルで消えた** |
| 点 | 同次座標 `x` | ★`§9-941` で構成物（存在は定理） |
| 点 | `hpt`（相異なる複素点） | ☆点の族の定義に属する条件 |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial HomogeneousLocalization NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★到達 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott まで
—— 点の側は同次座標だけ**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★点ごとに要るのは**同次座標**（`§9-941` で必ず取れる）と、
点の族が相異なる複素点を与えること（`§9-936`）だけである。
★★★★素点ごとの整合（`§9-940` の存在命題 2 本）は
`§9-942`・`§9-951` で**定理になった**。 -/
theorem northcott_of_isAmple_coords (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (hample : IsAmple M) (x₀ : (X : Type))
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type)) (d : ℕ) :
    ∃ (L : ℕ) (_ : 0 < L) (N : ℕ)
      (s : Fin (N + 1) → ((presheafTensorPow M L).obj (op ⊤) : Type))
      (hcov : (⨆ k, nonVanishing (presheafTensorPow M L) (s k)) = ⊤),
      (∀ i, IsAffineOpen (nonVanishing (presheafTensorPow M L) (s i))) ∧
      ∀ (D : ArithCartier X)
        (_ : D.divisor = divisorOfSection (presheafTensorPow M L) (s 0))
        (_ : D.green = fun p => greenFS N (p ≫ globalToProj (presheafTensorPow M L)
          (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        {P : Type}
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
        (_ : ∀ p, Module.finrank ℚ (fld p) ≤ d)
        (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
        (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
        (_ : ∀ p, x p ≠ 0) (_ : ∀ p, x p 0 ≠ 0)
        -- ★同次座標であること（`§9-941` で必ず取れる）
        (idx0 : P → Fin (N + 1))
        (hx : ∀ p, haveI := hnf p; Set.range
          (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
            ≫ globalToProj (presheafTensorPow M L)
              (isLocallyTrivial_presheafTensorPow hM L) φ s hcov).base
          ⊆ Set.range (chartA N ℤ (idx0 p)).base)
        (_ : ∀ p, x p (idx0 p) ≠ 0)
        (_ : ∀ p, haveI := hnf p; ∀ k, ((x p k : fld p))
          = projPointCoord N ℤ (fld p)
              (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
                ≫ globalToProj (presheafTensorPow M L)
                  (isLocallyTrivial_presheafTensorPow hM L) φ s hcov)
              (idx0 p) (hx p) k * ((x p (idx0 p) : fld p)))
        -- ★点の族が相異なる複素点を与えること
        (v₀ : ∀ p, haveI := hnf p; InfinitePlace (fld p))
        (idx : Fin (N + 1)) (_ : ∀ p, x p idx ≠ 0)
        (_ : ∀ p, haveI := hnf p; ∀ k,
          archRingHom (fld p) (v₀ p) (x p k) = ((x p k : fld p) : ℂ))
        (_ : Function.Injective (fun p => haveI := hnf p;
          archPoint (xF p) (v₀ p) ≫ globalToProj (presheafTensorPow M L)
            (isLocallyTrivial_presheafTensorPow hM L) φ s hcov))
        (C : ℝ),
        {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  obtain ⟨L, hL, N, s, hcov, haffs, key⟩ := northcott_of_isAmple_packed M hM hample x₀ φ d
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro D hDdiv hDgreen P fld hnf hdeg xF x hx0 h0 idx0 hx hxi hrel v₀ idx hidx hemb hpt C
  refine key D hDdiv hDgreen fld hnf hdeg xF x hx0 h0 ?_ ?_ v₀ idx hidx hemb hpt C
  · intro p Q hQ
    haveI := hnf p
    exact hfin_of_homogeneous_coords (presheafTensorPow M L)
      (isLocallyTrivial_presheafTensorPow hM L) φ s hcov (fld p) (xF p) (idx0 p) (hx p)
      (x p) (hxi p) (hrel p) Q hQ
  · intro p v
    haveI := hnf p
    exact harch_of_homogeneous_coords N (fld p) (xF p)
      (globalToProj (presheafTensorPow M L)
        (isLocallyTrivial_presheafTensorPow hM L) φ s hcov)
      (idx0 p) (hx p) (x p) (hxi p) (hrel p) v

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_isAmple_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott まで——点の側は同次座標だけ)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_isAmple_packed(束ねた形、§9-940)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_isAmple_packed") 4,
    .citation "[ABC3]" "hfin_of_homogeneous_coords(有限素点の整合、§9-951)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hfin_of_homogeneous_coords") 4,
    .citation "[ABC3]" "harch_of_homogeneous_coords(無限素点の整合、§9-942)"
      (.inProject "ABC3" "ABC3.Found.GenEll.harch_of_homogeneous_coords") 4,
    .citation "[ABC3]" "exists_homogeneous_coords(同次座標は必ず取れる、§9-941)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_homogeneous_coords") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): §9-940 が点ごとに受けていた存在命題 2 本" ++
       "(有限素点・無限素点の整合)は**どちらも定理になった**。" ++
       "★残る点の側のデータは、点とその同次座標(§9-941 で必ず取れる)と、" ++
       "点の族が相異なる複素点を与えること(§9-936)だけである") 6,
    .implicitStep
      ("★★Proposition 1.4, (iv) の現在地: 層の仮定(IsLocallyTrivial・IsAmple)は原文そのもの、" ++
       "幾何の仮定は §9-938 で消え、大域チャートは §9-930〜935＋§9-939 で消え、" ++
       "素点ごとの整合は本ファイルで消えた") 5 ]

end ABC3.Found.GenEll
