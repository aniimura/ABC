/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.RatiosLocalization
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★有限素点の整合は定理である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★★これは何か —— `§9-940` の 1 本目が消える

`§9-940`（束ねた形）の**有限素点の整合**

    `∀ Q, ∃ i y, (局所化した点が X_{s_i} を通る) ∧ (𝔞_Q = (x_i)) ∧ ((s_0/s_i)(y)·x_i = x_0)`

は、`§9-941` の同次座標を与えれば**定理として出る**。
★★`§9-942` で無限素点の側は既に定理になっているので、
これで `§9-940` の**存在命題 2 本がどちらも消えた**。

## ★★★機構 —— 本日の道具立てを 1 本に繋ぐ

| 段 | 使う場所 |
|---|---|
| `§9-950` | 座標から素点ごとの比の組 `r`（`r_j = 1`）と `𝔞_Q = (x_j)`、`r_0·x_j = x_0` |
| `§9-947` | 局所化した点 ＝ `projPointOfRatios N (𝓞_F)_Q r j` |
| `§9-948` | したがって `X_{s_j}` を通る（`y` が取れる） |
| `§9-949` | `(s_0/s_j)(y)·x_j = x_0` |

★配線の要は `IsScalarTower`——`𝓞_F → (𝓞_F)_Q → F` が `𝓞_F → F` と一致することで、
局所化した点を `F` へ落とすと**生成点になる**。

## ★これで `Proposition 1.4, (iv)` に残った点の側の仮定

★★★同次座標 `x`（`§9-941` で構成物）と、`hidx`・`hemb`・`hpt`
（点の族が相異なる複素点を与えること、`§9-936`）だけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial NumberField ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★有限素点の整合 -/

set_option maxHeartbeats 1000000 in
set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★★★★**`§9-940` の有限素点の整合は
同次座標から出る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★本日の道具立て（`§9-947`・`§9-948`・`§9-949`・`§9-950`）を 1 本に繋いだものである。
★★配線の要は `IsScalarTower`——局所化した点を `F` へ落とすと**生成点になる**。 -/
theorem hfin_of_homogeneous_coords (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤)
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ X)
    (i : Fin (N + 1))
    (hx : Set.range (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF
        ≫ globalToProj M hM φ s hcov).base ⊆ Set.range (chartA N ℤ i).base)
    (x : Fin (N + 1) → 𝓞 F) (hxi : x i ≠ 0)
    (hrel : ∀ k, ((x k : F)) = projPointCoord N ℤ F
      (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ globalToProj M hM φ s hcov)
      i hx k * ((x i : F)))
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal) :
    haveI := hQ.isPrime
    ∃ (j : Fin (N + 1)) (y : Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶
        (nonVanishing M (s j)).toScheme),
      Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
          = y ≫ (nonVanishing M (s j)).ι ∧
      Ideal.map (algebraMap (𝓞 F) (Localization.AtPrime Q)) (Ideal.span (Set.range x))
          = Ideal.span {algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)} ∧
      (Spec.preimage (y ≫ (nonVanishing M (s j)).toScheme.toSpecΓ)).hom
          ((nonVanishing M (s j)).topIso.inv.hom (globalRatio M hM (s 0) (s j)))
        * algebraMap (𝓞 F) (Localization.AtPrime Q) (x j)
        = algebraMap (𝓞 F) (Localization.AtPrime Q) (x 0) := by
  haveI := hQ.isPrime
  obtain ⟨j, r, hrj, hcj, hspan, hr0, hrF⟩ :=
    exists_ratios_localization F x i hxi
      (fun k => projPointCoord N ℤ F
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F)) ≫ xF ≫ globalToProj M hM φ s hcov)
        i hx k) hrel Q hQ
  have hgen : Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F))
      ≫ (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q)))
        ≫ xF ≫ globalToProj M hM φ s hcov)
      = Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) F))
          ≫ xF ≫ globalToProj M hM φ s hcov := by
    rw [← Category.assoc, ← Spec.map_comp, ← CommRingCat.ofHom_comp,
      ← IsScalarTower.algebraMap_eq]
  have hid : (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q)))
        ≫ xF ≫ globalToProj M hM φ s hcov)
      = projPointOfRatios N (Localization.AtPrime Q) r j hrj := by
    have key : ∀ (p : Spec (CommRingCat.of F) ⟶
          Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
        (hp : Set.range p.base ⊆ Set.range (chartA N ℤ i).base),
        Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F))
            ≫ (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q)))
              ≫ xF ≫ globalToProj M hM φ s hcov) = p →
        projPointCoord N ℤ F p i hp j ≠ 0 →
        (∀ k, algebraMap (Localization.AtPrime Q) F (r k) * projPointCoord N ℤ F p i hp j
          = projPointCoord N ℤ F p i hp k) →
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (Localization.AtPrime Q)))
            ≫ xF ≫ globalToProj M hM φ s hcov)
          = projPointOfRatios N (Localization.AtPrime Q) r j hrj := by
      rintro p hp rfl hcj' hr'
      exact localized_eq_projPointOfRatios N F Q hQ _ i hp j hcj' r hrj hr'
    exact key _ hx hgen hcj hrF
  have hrange : Set.range (Spec.map (CommRingCat.ofHom
        (algebraMap (𝓞 F) (Localization.AtPrime Q))) ≫ xF
      ≫ globalToProj M hM φ s hcov).base ⊆ Set.range (chartA N ℤ j).base := by
    rw [hid]
    exact range_projPointOfRatios N (Localization.AtPrime Q) r j hrj
  obtain ⟨y, hy⟩ := exists_localChart_at_index M hM φ s hcov F Q hQ xF j hrange
  have hβ' : y ≫ ((nonVanishing M (s j)).ι ≫ globalToProj M hM φ s hcov)
      = projPointOfRatios N (Localization.AtPrime Q) r j hrj := by
    rw [← Category.assoc, ← hy, Category.assoc]
    exact hid
  exact ⟨j, y, hy, hspan,
    hw_of_localized M hM φ s hcov j (Localization.AtPrime Q) y r hrj hβ' _ _ hr0⟩

/-! ## ★出典の紐付け(`.src`) -/

def hfin_of_homogeneous_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(有限素点の整合は同次座標から出る)",
    sectionId := "genell-prop-1-4" }

def hfin_of_homogeneous_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_ratios_localization(座標から素点ごとの比の組、§9-950)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_ratios_localization") 3,
    .citation "[ABC3]" "localized_eq_projPointOfRatios(局所化した点の同定、§9-947)"
      (.inProject "ABC3" "ABC3.Found.GenEll.localized_eq_projPointOfRatios") 3,
    .citation "[ABC3]" "exists_localChart_at_index(添字を指定した局所チャート、§9-948)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_localChart_at_index") 3,
    .citation "[ABC3]" "hw_of_localized((s_0/s_j)(y_Q)·x_j = x_0、§9-949)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hw_of_localized") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-940 の**存在命題 2 本がどちらも消えた**" ++
       "——無限素点の側は §9-942、有限素点の側は本ファイルである。" ++
       "★配線の要は IsScalarTower(𝓞_F → (𝓞_F)_Q → F が 𝓞_F → F と一致すること)で、" ++
       "局所化した点を F へ落とすと**生成点になる**") 6,
    .implicitStep
      ("★★これで Proposition 1.4, (iv) の点の側に残るのは" ++
       "同次座標 x(§9-941 で構成物)と hidx・hemb・hpt" ++
       "(点の族が相異なる複素点を与えること、§9-936)だけである") 5 ]

end ABC3.Found.GenEll
