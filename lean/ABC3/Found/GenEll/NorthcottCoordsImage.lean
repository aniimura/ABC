/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottImage
import ABC3.Found.GenEll.NorthcottCoords
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★単射性なしの `Proposition 1.4, (iv)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 点の側の仮定が座標だけになった

`§9-952`（`northcott_of_isAmple_coords`）は点の側に

* 同次座標 `x`（`§9-941` で必ず取れる）
* `hidx`・`hemb`・`hpt`（点の族が相異なる複素点を与えること）

を要求していた。★★`§9-960` で測ったとおり、**`hpt` は結論を添字集合 `P` に移すためだけ**の
条件であり、Northcott の内容ではない。

★★★★★★**本ファイルは結論を「座標の像」にすることで `hidx`・`hemb`・`hpt` を全部落とす**:

    `{p | ht(p) ≤ C}` の**正規化座標の像**は有限

★★★これで `Proposition 1.4, (iv)` の点の側に残るのは

    点 `xF p` と、その**同次座標** `x p`（`§9-941` で必ず取れる）

**だけ**になった。

## ★★★機構 —— 3 本を繋ぐだけ

| 段 | 内容 |
|---|---|
| `§9-951` | 有限素点の整合（座標から） |
| `§9-942` | 無限素点の整合（座標から） |
| `§9-933` | その 2 つから `ht = log H/[F:ℚ]` |
| `§9-960` | 高さが素朴高さなら**座標の像**は有限（`hinj` なし） |

★`choose` で `§9-951`・`§9-942` の存在命題をほどき、`§9-933` に渡し、`§9-960` で閉じる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov
open MvPolynomial NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★座標の像の有限性 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**貼った射に沿った Northcott（像の側）**
—— 単射性の仮定なし。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★点の側に要るのは**点とその同次座標**だけである
（素点ごとの整合は `§9-942`・`§9-951` で定理、`hpt` は `§9-960` で不要になった）。 -/
theorem northcott_globalToProj_of_local_image (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (d : ℕ)
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
    (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
    (hx0 : ∀ p, x p ≠ 0) (h0 : ∀ p, x p 0 ≠ 0)
    (idx0 : P → Fin (N + 1))
    (hx : ∀ p, haveI := hnf p; Set.range
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
        ≫ globalToProj M hM φ s hcov).base
      ⊆ Set.range (chartA N ℤ (idx0 p)).base)
    (hxi : ∀ p, x p (idx0 p) ≠ 0)
    (hrel : ∀ p, haveI := hnf p; ∀ k, ((x p k : fld p))
      = projPointCoord N ℤ (fld p)
          (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 (fld p)) (fld p))) ≫ xF p
            ≫ globalToProj M hM φ s hcov)
          (idx0 p) (hx p) k * ((x p (idx0 p) : fld p)))
    (idx : Fin (N + 1)) (C : ℝ) :
    ((fun (p : P) (k : Fin (N + 1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)) ''
      {p | haveI := hnf p;
        htArith (fld p) ((hyperplaneArith N).comap (globalToProj M hM φ s hcov)) (xF p)
          ≤ C}).Finite := by
  choose chartOf y hy hspan hw using fun p => haveI := hnf p;
    hfin_of_homogeneous_coords M hM φ s hcov (fld p) (xF p) (idx0 p) (hx p) (x p) (hxi p) (hrel p)
  choose archChart ρ hfac hcv hiv using fun p (v : haveI := hnf p; InfinitePlace (fld p)) =>
    haveI := hnf p;
    harch_of_homogeneous_coords N (fld p) (xF p) (globalToProj M hM φ s hcov)
      (idx0 p) (hx p) (x p) (hxi p) (hrel p) v
  refine northcott_of_log_mulHeight_image d fld hnf hdeg x _ (fun p => ?_) idx C
  haveI := hnf p
  exact htArith_globalToProj_eq_log_mulHeight_of_local M hM φ s hcov (fld p) (xF p)
    (x p) (hx0 p) (h0 p) (chartOf p) (y p) (hy p) (hspan p) (hw p)
    (archChart p) (ρ p) (hfac p) (hcv p) (hiv p)

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★`IsAmple` からの到達 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**`IsAmple` から Northcott
—— 点の側は同次座標だけ、単射性なし**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★点の側のデータは**点とその同次座標**（`§9-941` で必ず取れる）だけである
——`hidx`・`hemb`・`hpt` は結論を像にすることで**全部落ちた**。 -/
theorem northcott_of_isAmple_coords_image (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M)
    (hample : IsAmple M) (x₀ : (X : Type))
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type)) (d : ℕ) :
    ∃ (L : ℕ) (_ : 0 < L) (N : ℕ)
      (s : Fin (N + 1) → ((presheafTensorPow M L).obj (op ⊤) : Type))
      (hcov : (⨆ k, nonVanishing (presheafTensorPow M L) (s k)) = ⊤),
      (∀ i, IsAffineOpen (nonVanishing (presheafTensorPow M L) (s i))) ∧
      ∀ {P : Type}
        (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
        (_ : ∀ p, Module.finrank ℚ (fld p) ≤ d)
        (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
        (x : ∀ p, Fin (N + 1) → NumberField.RingOfIntegers (fld p))
        (_ : ∀ p, x p ≠ 0) (_ : ∀ p, x p 0 ≠ 0)
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
        (idx : Fin (N + 1)) (C : ℝ),
        ((fun (p : P) (k : Fin (N + 1)) =>
          ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)) ''
          {p | haveI := hnf p;
            htArith (fld p) ((hyperplaneArith N).comap
              (globalToProj (presheafTensorPow M L)
                (isLocallyTrivial_presheafTensorPow hM L) φ s hcov)) (xF p) ≤ C}).Finite := by
  obtain ⟨L, hL, N, s, hcov, haffs⟩ := exists_fin_cover_of_isAmple M hM hample x₀
  refine ⟨L, hL, N, s, hcov, haffs, ?_⟩
  intro P fld hnf hdeg xF x hx0 h0 idx0 hx hxi hrel idx C
  exact northcott_globalToProj_of_local_image (presheafTensorPow M L)
    (isLocallyTrivial_presheafTensorPow hM L) φ s hcov d fld hnf hdeg xF x hx0 h0
    idx0 hx hxi hrel idx C

/-! ## ★出典の紐付け(`.src`) -/

def northcott_globalToProj_of_local_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(貼った射に沿った Northcott——像の側、単射性なし)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_coords_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(IsAmple から Northcott——点の側は同次座標だけ、単射性なし)",
    sectionId := "genell-prop-1-4" }

def northcott_of_isAmple_coords_image.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_log_mulHeight_image(hinj の要らない Northcott、§9-960)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_log_mulHeight_image") 4,
    .citation "[ABC3]" "htArith_globalToProj_eq_log_mulHeight_of_local(局所形の高さ、§9-933)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_globalToProj_eq_log_mulHeight_of_local") 4,
    .citation "[ABC3]" "hfin_of_homogeneous_coords(有限素点の整合、§9-951)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hfin_of_homogeneous_coords") 3,
    .citation "[ABC3]" "harch_of_homogeneous_coords(無限素点の整合、§9-942)"
      (.inProject "ABC3" "ABC3.Found.GenEll.harch_of_homogeneous_coords") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 結論を「座標の像」にすることで" ++
       "hidx・hemb・hpt が**全部落ちた**" ++
       "——§9-960 で測ったとおり、それらは結論を添字集合 P に移すためだけの条件であり、" ++
       "Northcott の内容ではない") 6,
    .implicitStep
      ("★★★これで Proposition 1.4, (iv) の点の側に残るのは" ++
       "**点とその同次座標**(§9-941 で必ず取れる)だけである") 5 ]

end ABC3.Found.GenEll
