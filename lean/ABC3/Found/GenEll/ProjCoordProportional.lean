/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjPointOfCoords
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★同次座標はチャートを変えると定数倍になる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★これは何か —— `ArcModel` の中核

`ArcModel (ℙᴺ_ℤ) (ℂ^{N+1})` を作るには、複素点に**射影空間の点**を対応させる必要がある。
★正規化座標 `projPointCoord` はチャート `i` の選び方に依るが、
**射影的には依らない**——それが本ファイルの主張である:

    `coord_i(k) = coord_{i'}(k) · coord_i(i')`   （`coord_i(i') ≠ 0`）

★★右の因子は **`k` に依らない**ので、2 つの座標ベクトルは**定数倍**である。

## ★★★機構 —— `§9-869`・`§9-870` をそのまま点に当てる

★`§9-870` の `exists_overlap_factor`（2 つのチャートを通る点は重なりを通る）で
`φ = toOverlap ≫ χ`、`φ' = toOverlap' ≫ χ` と書け、
★★`§9-869` の `toOverlap_eq_mul`（`x_k/x_i = (x_k/x_{i'})·(x_{i'}/x_i)`）を当てるだけである。
★★★遷移単元は `toOverlap (x_{i'}/x_i)` そのものなので、
定数は **`φ(x_{i'}/x_i)` = `coord_i(i')`** と点の座標で書ける。

## ★★★★これで `ArcModel` の部品が揃う

| `ArcModel` の欄 | 部品 |
|---|---|
| `emb` | `projPointCoord`（本ファイルで射影的に well-defined と判る） |
| `emb_range` の全射性 | `projPointOfCoords`（`§9-873`） |
| `cone` | `univ`（閉錐） |
| `emb_injective` | `ext_of_projCoord`（`§9-850`）＋ 本ファイル |
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★★★★★★★★★★★環準同型の水準 -/

/-- ★★★★★★★★★★★★**2 つのチャートの正規化座標は定数倍である**（環準同型の水準）。

    `φ(x_k/x_i) = φ'(x_k/x_{i'}) · φ(x_{i'}/x_i)`

★右の因子は **`k` に依らない**。 -/
theorem projPointCoord_proportional (N : ℕ) (i i' : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (φ' : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i'))
      ⟶ CommRingCat.of ℂ)
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

/-! ## ★★★★★★★★★★★★点の水準 -/

/-- ★★★★★★★★★★★★**複素点の正規化座標はチャートを変えると定数倍になる**。 -/
theorem projPointCoord_congr (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (i i' : Fin (N+1))
    (hx : Set.range p.base ⊆ Set.range (chartA N ℤ i).base)
    (hx' : Set.range p.base ⊆ Set.range (chartA N ℤ i').base) (k : Fin (N+1)) :
    projPointCoord N ℤ ℂ p i hx k
      = projPointCoord N ℤ ℂ p i' hx' k * projPointCoord N ℤ ℂ p i hx i' := by
  have hfac : Spec.map (projChartHom N ℤ ℂ p i hx) ≫ chartA N ℤ i = p :=
    specMap_projChartHom N ℤ ℂ p i hx
  have hfac' : Spec.map (projChartHom N ℤ ℂ p i' hx') ≫ chartA N ℤ i' = p :=
    specMap_projChartHom N ℤ ℂ p i' hx'
  exact projPointCoord_proportional N i i' _ _ (by rw [hfac, hfac']) k

/-- ★★**その定数は `0` でない**。

★「点が `D₊(x_{i'})` にも入っている」ことの座標での言い方（`§9-871`）である。 -/
theorem projPointCoord_ne_zero (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (i i' : Fin (N+1))
    (hx : Set.range p.base ⊆ Set.range (chartA N ℤ i).base)
    (hx' : Set.range p.base ⊆ Set.range (chartA N ℤ i').base) :
    projPointCoord N ℤ ℂ p i hx i' ≠ 0 := by
  have hfac : Spec.map (projChartHom N ℤ ℂ p i hx) ≫ chartA N ℤ i = p :=
    specMap_projChartHom N ℤ ℂ p i hx
  refine projCoord_ne_zero_of_range N i i' (projChartHom N ℤ ℂ p i hx) ?_
  rw [hfac]
  exact hx'

/-- ★★★**座標ベクトルは定数倍で移り合う**（`ArcModel` の `emb` が well-defined な形）。 -/
theorem exists_smul_projPointCoord (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (i i' : Fin (N+1))
    (hx : Set.range p.base ⊆ Set.range (chartA N ℤ i).base)
    (hx' : Set.range p.base ⊆ Set.range (chartA N ℤ i').base) :
    ∃ t : ℂ, t ≠ 0 ∧ ∀ k, projPointCoord N ℤ ℂ p i hx k
      = t * projPointCoord N ℤ ℂ p i' hx' k :=
  ⟨projPointCoord N ℤ ℂ p i hx i', projPointCoord_ne_zero N p i i' hx hx',
    fun k => by rw [projPointCoord_congr N p i i' hx hx' k, mul_comm]⟩

/-! ## ★出典の紐付け(`.src`) -/

def projPointCoord_proportional.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(2 つのチャートの正規化座標は定数倍——環準同型の水準)",
    sectionId := "genell-prop-1-4" }

def projPointCoord_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(複素点の正規化座標はチャートを変えると定数倍)",
    sectionId := "genell-prop-1-4" }

def exists_smul_projPointCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(座標ベクトルは定数倍で移り合う)",
    sectionId := "genell-prop-1-4" }

def projPointCoord_proportional.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_overlap_factor(2 つのチャートを通る点は重なりを通る、§9-870)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_overlap_factor") 2,
    .citation "[ABC3]" "toOverlap_eq_mul(重なりの上での遷移、§9-869)"
      (.inProject "ABC3" "ABC3.Found.GenEll.toOverlap_eq_mul") 2,
    .citation "[ABC3]" "projCoord_ne_zero_of_range(座標の非零性、§9-871)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projCoord_ne_zero_of_range") 2,
    .implicitStep
      ("★遷移単元は toOverlap (x_{i'}/x_i) そのものなので、" ++
       "定数は **φ(x_{i'}/x_i) = coord_i(i')** と点の座標で書ける" ++
       "——存在量化子を使わずに済む") 3,
    .implicitStep
      ("★★これで ArcModel (ℙᴺ_ℤ) (ℂ^{N+1}) の部品が揃った: " ++
       "emb は projPointCoord(本ファイルで射影的に well-defined と判る)、" ++
       "emb_range の全射性は projPointOfCoords(§9-873)、cone は univ(閉錐)、" ++
       "emb_injective は ext_of_projCoord(§9-850)＋本ファイル") 4 ]

end ABC3.Found.GenEll
