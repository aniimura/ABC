/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.SeparatedUnique
import ABC3.Found.GenEll.ProjArcModel
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★チャートの取り替えは任意の体で動く（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★★これは何か —— `ℂ` 専用だった 3 本を一般化する

`§9-871`〜`§9-874` の

* `specMap_mem_basicOpen_iff`（`Spec` の点が基本開集合に入る条件）
* `range_of_projCoord_ne_zero`（座標が `0` でなければチャートに入る）
* `projPointCoord_proportional`・`projPointCoord_congr`（チャートを変えると定数倍）

は `ℂ` に特殊化して書かれていたが、★**証明はどれも体（あるいは可換環）で通る**。
本ファイルはそれを一般化し、その上で

    ★★★★**比の組 `r`（`r_j = 1`）が点の正規化座標と合っていれば、
    その点は `projPointOfRatios N F r j` そのものである**

を取る。

## ★★★機構 —— チャート `j` に移してから `ext_of_projCoord`

1. `c_j ≠ 0` なら点は `D₊(x_j)` にも入る（`range_of_projCoord_ne_zero'`）
2. そこでの正規化座標は `c_k/c_j`（`projPointCoord_congr'`）——仮定よりこれが `r_k`
3. `A⁰_{x_j}` からの環準同型は `projCoord` の像で決まる（`§9-849` の `ext_of_projCoord`）
4. ★したがって点のチャート射は `awayHomOfRatios`（`§9-944`）と一致する

## ★これで何が閉じるか

★★★`§9-943`（有限素点では座標の 1 つが他を割る）＋`§9-944`（比の組から点を作る）
＋`§9-945`（分離的なら生成点で決まる）＋**本ファイル**で、

    `Spec (𝓞_F)_Q ⟶ ℙᴺ` が `D₊(x_j)` を通る

が言える道が全部揃った。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★体の点と基本開集合 -/

/-- ★**体の点が基本開集合に入るのは値が `0` でないとき**——`§9-871` の一般化。 -/
theorem specMap_mem_basicOpen_iff' {A : CommRingCat.{0}} (F : Type) [Field F]
    (φ : A ⟶ CommRingCat.of F) (a : A) :
    (Spec.map φ).base (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of F))
        ∈ PrimeSpectrum.basicOpen a ↔ φ.hom a ≠ 0 := by
  show PrimeSpectrum.comap φ.hom (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of F))
      ∈ PrimeSpectrum.basicOpen a ↔ _
  rw [PrimeSpectrum.mem_basicOpen]
  show a ∉ Ideal.comap φ.hom (⊥ : Ideal (CommRingCat.of F)) ↔ _
  rw [Ideal.mem_comap, Ideal.mem_bot]

/-- ★★**座標が `0` でなければチャートに入る**——`§9-874` の一般化（任意の体で）。 -/
theorem range_of_projCoord_ne_zero' (N : ℕ) (F : Type) [Field F] (i i' : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of F)
    (hne : φ.hom (projCoord N ℤ i i') ≠ 0) :
    Set.range (Spec.map φ ≫ chartA N ℤ i).base ⊆ Set.range (chartA N ℤ i').base := by
  haveI : Subsingleton (Spec (CommRingCat.of F)) :=
    inferInstanceAs (Subsingleton (PrimeSpectrum F))
  have hpre : (Spec.map φ).base (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of F))
      ∈ PrimeSpectrum.basicOpen (projCoord N ℤ i i') :=
    (specMap_mem_basicOpen_iff' F φ _).2 hne
  rw [← isLocalizationElem_eq_projCoord,
    ← Proj.awayι_preimage_basicOpen _ (MvPolynomial.isHomogeneous_X ℤ i) one_pos
      (MvPolynomial.isHomogeneous_X ℤ i') one_pos] at hpre
  have hmem : (Spec.map φ ≫ chartA N ℤ i).base
      (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of F))
      ∈ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i') := hpre
  rintro _ ⟨y, rfl⟩
  have hy : y = (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of F)) :=
    Subsingleton.elim _ _
  subst hy
  rw [← Proj.opensRange_awayι (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
    (MvPolynomial.X i') (MvPolynomial.isHomogeneous_X ℤ i') one_pos] at hmem
  exact hmem

/-! ## ★★チャートを変えると定数倍（任意の可換環で） -/

/-- ★★**正規化座標はチャートを変えると定数倍**——`§9-874` の一般化（任意の可換環で）。 -/
theorem projPointCoord_proportional' (N : ℕ) (i i' : Fin (N+1)) {C : CommRingCat.{0}}
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i)) ⟶ C)
    (φ' : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i')) ⟶ C)
    (hp : Spec.map φ ≫ chartA N ℤ i = Spec.map φ' ≫ chartA N ℤ i') (k : Fin (N+1)) :
    φ.hom (projCoord N ℤ i k)
      = φ'.hom (projCoord N ℤ i' k) * φ.hom (projCoord N ℤ i i') := by
  obtain ⟨χ, h1, h2⟩ := exists_overlap_factor N i i' φ φ' hp
  have hc : ∀ (a : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i)),
      φ.hom a = χ.hom (toOverlap N i i' a) := by
    intro a; rw [← h1, CommRingCat.hom_comp, CommRingCat.hom_ofHom]; rfl
  have hc' : ∀ (a : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i')),
      φ'.hom a = χ.hom (toOverlap' N i i' a) := by
    intro a; rw [← h2, CommRingCat.hom_comp, CommRingCat.hom_ofHom]; rfl
  rw [hc, hc, hc', toOverlap_eq_mul, map_mul, toOverlap_projCoord_self]

/-- ★★**点の正規化座標はチャートを変えると定数倍**（任意の体で）。 -/
theorem projPointCoord_congr' (N : ℕ) (F : Type) [Field F]
    (p : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i i' : Fin (N+1))
    (hx : Set.range p.base ⊆ Set.range (chartA N ℤ i).base)
    (hx' : Set.range p.base ⊆ Set.range (chartA N ℤ i').base) (k : Fin (N+1)) :
    projPointCoord N ℤ F p i hx k
      = projPointCoord N ℤ F p i' hx' k * projPointCoord N ℤ F p i hx i' := by
  have hfac : Spec.map (projChartHom N ℤ F p i hx) ≫ chartA N ℤ i = p :=
    specMap_projChartHom N ℤ F p i hx
  have hfac' : Spec.map (projChartHom N ℤ F p i' hx') ≫ chartA N ℤ i' = p :=
    specMap_projChartHom N ℤ F p i' hx'
  exact projPointCoord_proportional' N i i' _ _ (by rw [hfac, hfac']) k

/-! ## ★★★★★★★★★★★★★★★★★★★比の組が点を決める -/

/-- ★★★★★★★★★★★★★★★★★★★**比の組が正規化座標と合っていれば、
点は `projPointOfRatios` そのものである**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は 3 段:
* `c_j ≠ 0` なら点は `D₊(x_j)` にも入る（`range_of_projCoord_ne_zero'`）
* そこでの正規化座標は `c_k/c_j`（`projPointCoord_congr'`）——仮定よりこれが `r_k`
* `A⁰_{x_j}` からの環準同型は `projCoord` の像で決まる（`§9-849` の `ext_of_projCoord`）

★★これが `§9-944` で構成した点と実際の点を結ぶ最後の等式である。 -/
theorem projPointOfRatios_eq_of_coords (N : ℕ) (F : Type) [Field F]
    (q : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (i : Fin (N + 1)) (hx : Set.range q.base ⊆ Set.range (chartA N ℤ i).base)
    (j : Fin (N + 1)) (hcj : projPointCoord N ℤ F q i hx j ≠ 0)
    (r : Fin (N + 1) → F) (hrj : r j = 1)
    (hr : ∀ k, r k * projPointCoord N ℤ F q i hx j = projPointCoord N ℤ F q i hx k) :
    q = projPointOfRatios N F r j hrj := by
  have hfac : Spec.map (projChartHom N ℤ F q i hx) ≫ chartA N ℤ i = q :=
    specMap_projChartHom N ℤ F q i hx
  have hxj : Set.range q.base ⊆ Set.range (chartA N ℤ j).base := by
    rw [← hfac]
    exact range_of_projCoord_ne_zero' N F i j (projChartHom N ℤ F q i hx) hcj
  have hcoord : ∀ k, projPointCoord N ℤ F q j hxj k = r k := by
    intro k
    have hcg := projPointCoord_congr' N F q i j hx hxj k
    have h1 : projPointCoord N ℤ F q j hxj k * projPointCoord N ℤ F q i hx j
        = r k * projPointCoord N ℤ F q i hx j := by
      rw [← hcg]; exact (hr k).symm
    exact mul_right_cancel₀ hcj h1
  have hring : projChartHom N ℤ F q j hxj = CommRingCat.ofHom (awayHomOfRatios N F r j hrj) := by
    refine CommRingCat.hom_ext ?_
    refine ext_of_projCoord j _ _ (fun c => ?_) (fun k => ?_)
    · rw [CommRingCat.hom_ofHom, awayConst_eq_intCast, map_intCast, map_intCast]
    · rw [CommRingCat.hom_ofHom, awayHomOfRatios_projCoord]
      exact hcoord k
  rw [projPointOfRatios, ← hring, specMap_projChartHom N ℤ F q j hxj]

/-! ## ★出典の紐付け(`.src`) -/

def specMap_mem_basicOpen_iff'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(体の点が基本開集合に入るのは値が 0 でないとき)",
    sectionId := "genell-prop-1-4" }

def range_of_projCoord_ne_zero'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(座標が 0 でなければチャートに入る——任意の体で)",
    sectionId := "genell-prop-1-4" }

def projPointCoord_proportional'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(正規化座標はチャートを変えると定数倍——任意の可換環で)",
    sectionId := "genell-prop-1-4" }

def projPointCoord_congr'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(点の正規化座標はチャートを変えると定数倍——任意の体で)",
    sectionId := "genell-prop-1-4" }

def projPointOfRatios_eq_of_coords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(比の組が正規化座標と合っていれば点は projPointOfRatios である)",
    sectionId := "genell-prop-1-4" }

def projPointOfRatios_eq_of_coords.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ext_of_projCoord(A⁰_{x_j} からの環準同型は projCoord で決まる、§9-849)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_projCoord") 3,
    .citation "[ABC3]" "awayHomOfRatios(比の組が定めるチャート射、§9-944)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayHomOfRatios") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-871〜874 が ℂ に特殊化して書かれていたのは" ++
       "**書き方の都合**であって、証明はどれも体(あるいは可換環)で通る。" ++
       "一般化は書き直しではなく被せである") 4,
    .implicitStep
      ("★★★これで §9-943(座標の 1 つが他を割る)＋§9-944(比の組から点を作る)" ++
       "＋§9-945(分離的なら生成点で決まる)＋本ファイルで、" ++
       "Spec (𝓞_F)_Q ⟶ ℙᴺ が D₊(x_j) を通ることを言う道が全部揃った") 5 ]

end ABC3.Found.GenEll
