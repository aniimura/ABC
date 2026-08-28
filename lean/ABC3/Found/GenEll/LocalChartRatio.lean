/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LocalChartAtIndex
import ABC3.Found.GenEll.PullbackChartLocal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★チャート射の値は座標の比である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— `§9-940` の `hw`

`§9-940` の有限素点の整合に残っていた最後の一本は

    **`(s_0/s_j)(y_Q) · x_j = x_0`**

である。★`§9-943` は `g·x_j = x_0` なる `g` を与えているので、
言うべきことは **`(s_0/s_j)(y_Q) = g`** だけである。

## ★★★機構 —— チャート射は `ℙᴺ` の点で決まる

`§9-947` で「局所化した点 `β` は `projPointOfRatios N R r j`」が分かっている。
一方 `β = y ≫ ι ≫ ψ = Spec(チャート射) ≫ chartA j` でもある。
★★`chartA j` は**開埋め込み＝モノ**だから `Spec.map` の引数が一致し、
`Spec.map` は忠実なので

    **チャート射 ＝ `awayHomOfRatios N R r j`**

★したがって `projCoord j 0` の値は `r_0 = x_0/x_j` である。
★★あとは `§9-886`（チャート射の値は比の切断の値）でそれを `(s_0/s_j)(y_Q)` に読み替える。

## ★これで `Proposition 1.4, (iv)` の有限素点の整合は 3 本とも揃った

| `§9-940` の要求 | 出どころ |
|---|---|
| 局所化した点が `X_{s_j}` を通る | `§9-947`＋`§9-948` |
| `𝔞_Q = (x_j)` | `§9-943` |
| `(s_0/s_j)(y_Q)·x_j = x_0` | ★**本ファイル** |
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial NumberField ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★★★★★★チャート射の同定 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★**点のチャート射は比の組が定める射である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`chartA j` は開埋め込み＝**モノ**だから `Spec.map` の引数が一致し、
`Spec.map` は忠実なので環準同型の水準で一致する。 -/
theorem chartHom_eq_awayHomOfRatios (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (j : Fin (N + 1))
    (R : Type) [CommRing R]
    (y : Spec (CommRingCat.of R) ⟶ (nonVanishing M (s j)).toScheme)
    (r : Fin (N + 1) → R) (hrj : r j = 1)
    (hβ : y ≫ ((nonVanishing M (s j)).ι ≫ globalToProj M hM φ s hcov)
      = projPointOfRatios N R r j hrj) :
    Spec.preimage (y ≫ globalChartMorphism M hM φ s j)
      = CommRingCat.ofHom (awayHomOfRatios N R r j hrj) := by
  apply Spec.map_injective
  rw [Spec.map_preimage]
  have h1 : y ≫ globalChartMorphism M hM φ s j ≫ chartA N ℤ j
      = y ≫ ((nonVanishing M (s j)).ι ≫ globalToProj M hM φ s hcov) := by
    rw [ι_globalToProj, globalChartToProj]
  rw [← cancel_mono (chartA N ℤ j), Category.assoc, h1, hβ, projPointOfRatios]

/-! ## ★★★★★★★★★★★★★★★★★★★★比の切断の値 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★★★★★**`(s_0/s_j)` の点での値は `r_0` である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-886`（チャート射の値は比の切断の値）と `chartHom_eq_awayHomOfRatios` を繋ぐだけである。 -/
theorem globalRatio_value_eq (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (j : Fin (N + 1))
    (R : Type) [CommRing R]
    (y : Spec (CommRingCat.of R) ⟶ (nonVanishing M (s j)).toScheme)
    (r : Fin (N + 1) → R) (hrj : r j = 1)
    (hβ : y ≫ ((nonVanishing M (s j)).ι ≫ globalToProj M hM φ s hcov)
      = projPointOfRatios N R r j hrj) :
    (Spec.preimage (y ≫ (nonVanishing M (s j)).toScheme.toSpecΓ)).hom
      ((nonVanishing M (s j)).topIso.inv.hom (globalRatio M hM (s 0) (s j))) = r 0 := by
  rw [← preimage_globalChartMorphism_projCoord' M hM φ s j 0 (CommRingCat.of R) y,
    chartHom_eq_awayHomOfRatios M hM φ s hcov j R y r hrj hβ, CommRingCat.hom_ofHom,
    awayHomOfRatios_projCoord]

/-- ★★★★★★★★★★★★★★★★★★★★★★**`§9-940` の `hw` はこれで出る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-943` が `r_0·x_j = x_0` を与え、本ファイルが `(s_0/s_j)(y_Q) = r_0` を与える。 -/
theorem hw_of_localized (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (hcov : (⨆ k, nonVanishing M (s k)) = ⊤) (j : Fin (N + 1))
    (R : Type) [CommRing R]
    (y : Spec (CommRingCat.of R) ⟶ (nonVanishing M (s j)).toScheme)
    (r : Fin (N + 1) → R) (hrj : r j = 1)
    (hβ : y ≫ ((nonVanishing M (s j)).ι ≫ globalToProj M hM φ s hcov)
      = projPointOfRatios N R r j hrj)
    (a b : R) (hg : r 0 * a = b) :
    (Spec.preimage (y ≫ (nonVanishing M (s j)).toScheme.toSpecΓ)).hom
        ((nonVanishing M (s j)).topIso.inv.hom (globalRatio M hM (s 0) (s j))) * a = b := by
  rw [globalRatio_value_eq M hM φ s hcov j R y r hrj hβ]
  exact hg

/-! ## ★出典の紐付け(`.src`) -/

def chartHom_eq_awayHomOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(点のチャート射は比の組が定める射である)",
    sectionId := "genell-prop-1-4" }

def globalRatio_value_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)((s_0/s_j) の点での値は r_0 である)",
    sectionId := "genell-prop-1-4" }

def hw_of_localized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)((s_0/s_j)(y_Q)·x_j = x_0)",
    sectionId := "genell-prop-1-4" }

def hw_of_localized.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "localized_eq_projPointOfRatios(局所化した点の同定、§9-947)"
      (.inProject "ABC3" "ABC3.Found.GenEll.localized_eq_projPointOfRatios") 3,
    .citation "[ABC3]" "preimage_globalChartMorphism_projCoord'(チャート射の値は比の切断の値、§9-930)"
      (.inProject "ABC3" "ABC3.Found.GenEll.preimage_globalChartMorphism_projCoord'") 2,
    .citation "[ABC3]" "exists_span_and_ratio_localization(r_0·x_j = x_0、§9-943)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_span_and_ratio_localization") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-940 の有限素点の整合に残っていた最後の一本" ++
       "((s_0/s_j)(y_Q)·x_j = x_0)は、**チャート射が ℙᴺ の点で決まる**ことから出る。" ++
       "chartA j は開埋め込み＝モノで Spec.map は忠実だから、" ++
       "点のチャート射は awayHomOfRatios と一致する") 5,
    .implicitStep
      ("★★これで Proposition 1.4, (iv) の有限素点の整合は 3 本とも揃った: " ++
       "局所化した点が X_{s_j} を通る(§9-947＋948)、𝔞_Q = (x_j)(§9-943)、" ++
       "(s_0/s_j)(y_Q)·x_j = x_0(本ファイル)") 5 ]

end ABC3.Found.GenEll
