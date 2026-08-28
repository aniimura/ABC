/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartChangeField
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★局所化した点はそのチャートを通る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★★これは何か —— 有限素点の紐が結ばれた

`§9-943`（座標の 1 つが他を割る）で得た `j` に対し、比 `r_k ≔ x_k/x_j ∈ (𝓞_F)_Q` が取れる。
`§9-944` はそこから `Spec (𝓞_F)_Q ⟶ D₊(x_j) ⊆ ℙᴺ` を**構成**した。

★★★★本ファイルはそれが**局所化した点そのものである**ことを示す:

    `β = projPointOfRatios N (𝓞_F)_Q r j`

——したがって★★**局所化した点は `D₊(x_j)` を通る**。

## ★★★機構 —— 生成点に降ろして比べる

1. `projPointOfRatios` は環準同型に沿って自然（`specMap_projPointOfRatios`）
   ——`(𝓞_F)_Q → F` で送れば `F` 係数の比の組の点になる
2. `F` の水準では `§9-946`（比の組が点を決める）が使える
3. ★`§9-945`（分離的なら生成点で決まる）で `(𝓞_F)_Q` の水準に戻す

★★どの段も `ℙᴺ` の分離性と `(𝓞_F)_Q` が付値環であることしか使っていない。

## ★これで `Proposition 1.4, (iv)` に残った段

★★★残るのは**紐の付け替え**だけである——`ℙᴺ` の側で得たこの事実を
`X` の側（`X_{s_j}` を通ること）に移すこと。機構は `§9-913`
（`ψ⁻¹(D₊(x_j)) = X_{s_j}`）である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MvPolynomial NumberField

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★比の組の点は環準同型に沿って自然である -/

/-- ★**`awayHomOfRatios` は後合成と可換**。 -/
theorem awayHomOfRatios_comp (N : ℕ) (R S : Type) [CommRing R] [CommRing S] (g : R →+* S)
    (r : Fin (N+1) → R) (j : Fin (N+1)) (hrj : r j = 1) :
    g.comp (awayHomOfRatios N R r j hrj)
      = awayHomOfRatios N S (fun k => g (r k)) j (by simp [hrj]) := by
  refine ext_of_projCoord j _ _ (fun c => ?_) (fun k => ?_)
  · simp [awayConst_eq_intCast]
  · rw [RingHom.comp_apply, awayHomOfRatios_projCoord, awayHomOfRatios_projCoord]

/-- ★★**比の組が定める点は環準同型に沿って自然である**。 -/
theorem specMap_projPointOfRatios (N : ℕ) (R S : Type) [CommRing R] [CommRing S] (g : R →+* S)
    (r : Fin (N+1) → R) (j : Fin (N+1)) (hrj : r j = 1) :
    Spec.map (CommRingCat.ofHom g) ≫ projPointOfRatios N R r j hrj
      = projPointOfRatios N S (fun k => g (r k)) j (by simp [hrj]) := by
  rw [projPointOfRatios, projPointOfRatios, ← Category.assoc, ← Spec.map_comp,
    ← CommRingCat.ofHom_comp, awayHomOfRatios_comp]

/-! ## ★★★★★★★★★★★★★★★★★★★★局所化した点の同定 -/

/-- ★★★★★★★★★★★★★★★★★★★★**局所化した点は比の組が定める点である**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は 3 段——`projPointOfRatios` の自然性、`§9-946`（比の組が点を決める、体の上で）、
`§9-945`（分離的なら生成点で決まる）。 -/
theorem localized_eq_projPointOfRatios (N : ℕ) (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal)
    (β : haveI := hQ.isPrime; Spec (CommRingCat.of (Localization.AtPrime Q))
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i : Fin (N + 1))
    (hx : haveI := hQ.isPrime; Set.range
      (Spec.map (CommRingCat.ofHom
        (algebraMap (Localization.AtPrime Q) F)) ≫ β).base
      ⊆ Set.range (chartA N ℤ i).base)
    (j : Fin (N + 1))
    (hcj : haveI := hQ.isPrime; projPointCoord N ℤ F
      (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx j ≠ 0)
    (rR : haveI := hQ.isPrime; Fin (N + 1) → Localization.AtPrime Q)
    (hrj : haveI := hQ.isPrime; rR j = 1)
    (hr : haveI := hQ.isPrime; ∀ k,
      algebraMap (Localization.AtPrime Q) F (rR k)
        * projPointCoord N ℤ F
            (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx j
      = projPointCoord N ℤ F
          (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx k) :
    haveI := hQ.isPrime
    β = projPointOfRatios N (Localization.AtPrime Q) rR j hrj := by
  haveI := hQ.isPrime
  refine eq_of_generic_eq_atPrime F Q hQ β _ ?_
  rw [specMap_projPointOfRatios]
  exact projPointOfRatios_eq_of_coords N F _ i hx j hcj
    (fun k => algebraMap (Localization.AtPrime Q) F (rR k)) (by simp [hrj]) hr

/-- ★★★★★★★★★★★★★★★★★★★★★**局所化した点は `D₊(x_j)` を通る**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★これが `§9-940` の有限素点の整合の**幾何の紐**である
——`§9-943` が与える `j`（座標の最小割り切り成分）で、点は実際にチャートを通る。 -/
theorem range_localized_subset_chart (N : ℕ) (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal)
    (β : haveI := hQ.isPrime; Spec (CommRingCat.of (Localization.AtPrime Q))
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i : Fin (N + 1))
    (hx : haveI := hQ.isPrime; Set.range
      (Spec.map (CommRingCat.ofHom
        (algebraMap (Localization.AtPrime Q) F)) ≫ β).base
      ⊆ Set.range (chartA N ℤ i).base)
    (j : Fin (N + 1))
    (hcj : haveI := hQ.isPrime; projPointCoord N ℤ F
      (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx j ≠ 0)
    (rR : haveI := hQ.isPrime; Fin (N + 1) → Localization.AtPrime Q)
    (hrj : haveI := hQ.isPrime; rR j = 1)
    (hr : haveI := hQ.isPrime; ∀ k,
      algebraMap (Localization.AtPrime Q) F (rR k)
        * projPointCoord N ℤ F
            (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx j
      = projPointCoord N ℤ F
          (Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) i hx k) :
    haveI := hQ.isPrime
    Set.range β.base ⊆ Set.range (chartA N ℤ j).base := by
  haveI := hQ.isPrime
  rw [localized_eq_projPointOfRatios N F Q hQ β i hx j hcj rR hrj hr]
  exact range_projPointOfRatios N (Localization.AtPrime Q) rR j hrj

/-! ## ★出典の紐付け(`.src`) -/

def awayHomOfRatios_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(awayHomOfRatios は後合成と可換)",
    sectionId := "genell-prop-1-4" }

def specMap_projPointOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(比の組が定める点は環準同型に沿って自然である)",
    sectionId := "genell-prop-1-4" }

def localized_eq_projPointOfRatios.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所化した点は比の組が定める点である)",
    sectionId := "genell-prop-1-4" }

def range_localized_subset_chart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(局所化した点は D₊(x_j) を通る)",
    sectionId := "genell-prop-1-4" }

def range_localized_subset_chart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "eq_of_generic_eq_atPrime(分離的なら生成点で決まる、§9-945)"
      (.inProject "ABC3" "ABC3.Found.GenEll.eq_of_generic_eq_atPrime") 3,
    .citation "[ABC3]" "projPointOfRatios_eq_of_coords(比の組が点を決める、§9-946)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projPointOfRatios_eq_of_coords") 3,
    .citation "[ABC3]" "exists_span_and_ratio_localization(座標の 1 つが他を割る、§9-943)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_span_and_ratio_localization") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-940 の有限素点の整合の**幾何の紐**が結ばれた" ++
       "——§9-943 が与える j(座標の最小割り切り成分)で、" ++
       "局所化した点は実際に D₊(x_j) を通る。" ++
       "機構は projPointOfRatios の自然性 ＋ 体の水準での同定(§9-946) ＋ " ++
       "分離性の付値判定法(§9-945)の 3 段である") 5,
    .implicitStep
      ("★残るのは紐の付け替えだけである——ℙᴺ の側で得たこの事実を " ++
       "X の側(X_{s_j} を通ること)に移すこと。機構は §9-913(ψ⁻¹(D₊(x_j)) = X_{s_j})である") 4 ]

end ABC3.Found.GenEll
