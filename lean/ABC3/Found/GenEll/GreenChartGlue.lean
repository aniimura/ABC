/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GreenChartIndep
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★2 つのチャートに入る複素点は重なりを通る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★これは何か —— 大域的な `green` の**スキームの側**

`§9-869` は「Fubini–Study はチャートに依らない」を**代数の側**（重なりのチャートを
通る点について）で取った。★本ファイルはその**スキームの側**を埋める:

> 同じ点が 2 つのチャートを通るなら、**重なりのチャートも通る**

★★これで `greenChartOf` の値がチャートに依らないことが、
「重なりを通る」という仮定なしに言える（`greenChartOf_congr`）。

## ★★★機構 —— mathlib の 3 本

| 道具 | 役割 |
|---|---|
| `Proj.basicOpen_mul` | `D₊(x_i x_{i'}) = D₊(x_i) ⊓ D₊(x_{i'})` |
| `IsOpenImmersion.lift` | 像が入っていれば持ち上がる |
| `Proj.SpecMap_awayMap_awayι` | 重なりのチャートは各チャートを経由する |

★あとは `chartA` が**モノ**（開埋め込み）であることで因子を消すだけである。

## ★★★★★測定の記録 —— 非退化条件は `φ` だけで書ける

★`§9-869` の `greenChartOf_overlap` は「遷移単元 `χ(t) ≠ 0`」を要求するが、
`t = x_{i'}/x_i` は `toOverlap (x_{i'}/x_i)` そのものなので
（`toOverlap_projCoord_self`）、条件は **`φ(x_{i'}/x_i) ≠ 0`** と、
`χ` を持ち出さずに書ける（2026-08-28 実測）。
★★これは「点が `D₊(x_{i'})` にも入っている」ことの座標での言い方である。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★重なりのチャート -/

/-- ★**重なりのチャート** `Spec A⁰_{x_i x_{i'}} ⟶ Proj 𝒜`。 -/
noncomputable abbrev chartOverlap (N : ℕ) (i i' : Fin (N+1)) :
    Spec (.of <| HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
      (MvPolynomial.X i * MvPolynomial.X i'))
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) :=
  Proj.awayι _ (MvPolynomial.X i * MvPolynomial.X i') (X_mul_X_mem N i i') two_pos

/-- ★**重なりのチャートの像は 2 つのチャートの像の共通部分である**。 -/
theorem opensRange_chartOverlap (N : ℕ) (i i' : Fin (N+1)) :
    (chartOverlap N i i').opensRange
      = (chartA N ℤ i).opensRange ⊓ (chartA N ℤ i').opensRange := by
  rw [Proj.opensRange_awayι, Proj.opensRange_awayι, Proj.opensRange_awayι,
    Proj.basicOpen_mul]

/-- ★★**2 つのチャートを通る点の像は重なりの像に入る**。 -/
theorem range_subset_chartOverlap (N : ℕ) (i i' : Fin (N+1)) {C : CommRingCat.{0}}
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i)) ⟶ C)
    (φ' : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i')) ⟶ C)
    (hp : Spec.map φ ≫ chartA N ℤ i = Spec.map φ' ≫ chartA N ℤ i') :
    Set.range (Spec.map φ ≫ chartA N ℤ i).base
      ⊆ Set.range (chartOverlap N i i').base := by
  rintro _ ⟨y, rfl⟩
  have hy1 : (Spec.map φ ≫ chartA N ℤ i).base y ∈ (chartA N ℤ i).opensRange :=
    ⟨(Spec.map φ).base y, rfl⟩
  have hy2 : (Spec.map φ ≫ chartA N ℤ i).base y ∈ (chartA N ℤ i').opensRange := by
    rw [hp]; exact ⟨(Spec.map φ').base y, rfl⟩
  have h : (Spec.map φ ≫ chartA N ℤ i).base y ∈ (chartOverlap N i i').opensRange := by
    rw [opensRange_chartOverlap]; exact ⟨hy1, hy2⟩
  exact h

/-! ## ★★重なりのチャートは各チャートを経由する -/

theorem specMap_toOverlap_comp (N : ℕ) (i i' : Fin (N+1)) :
    Spec.map (CommRingCat.ofHom (toOverlap N i i')) ≫ chartA N ℤ i = chartOverlap N i i' :=
  Proj.SpecMap_awayMap_awayι _ (MvPolynomial.isHomogeneous_X ℤ i) one_pos
    (MvPolynomial.isHomogeneous_X ℤ i') rfl

theorem specMap_toOverlap'_comp (N : ℕ) (i i' : Fin (N+1)) :
    Spec.map (CommRingCat.ofHom (toOverlap' N i i')) ≫ chartA N ℤ i' = chartOverlap N i i' :=
  Proj.SpecMap_awayMap_awayι _ (MvPolynomial.isHomogeneous_X ℤ i') one_pos
    (MvPolynomial.isHomogeneous_X ℤ i) (mul_comm _ _)

/-! ## ★★★★★★★★★★重なりを経由する分解 -/

/-- ★★★★★★★★★★**2 つのチャートを通る点は重なりのチャートを通る**。

★機構は `IsOpenImmersion.lift`（像が入っていれば持ち上がる）と
`chartA` が**モノ**であること（開埋め込み）だけである。 -/
theorem exists_overlap_factor (N : ℕ) (i i' : Fin (N+1)) {C : CommRingCat.{0}}
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i)) ⟶ C)
    (φ' : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i')) ⟶ C)
    (hp : Spec.map φ ≫ chartA N ℤ i = Spec.map φ' ≫ chartA N ℤ i') :
    ∃ χ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i * MvPolynomial.X i')) ⟶ C,
      CommRingCat.ofHom (toOverlap N i i') ≫ χ = φ ∧
      CommRingCat.ofHom (toOverlap' N i i') ≫ χ = φ' := by
  have hrange := range_subset_chartOverlap N i i' φ φ' hp
  have hfac := IsOpenImmersion.lift_fac (chartOverlap N i i')
    (Spec.map φ ≫ chartA N ℤ i) hrange
  refine ⟨Spec.preimage (IsOpenImmersion.lift (chartOverlap N i i')
    (Spec.map φ ≫ chartA N ℤ i) hrange), ?_, ?_⟩
  · apply Spec.map_injective
    rw [Spec.map_comp, Spec.map_preimage, ← cancel_mono (chartA N ℤ i), Category.assoc,
      specMap_toOverlap_comp]
    exact hfac
  · apply Spec.map_injective
    rw [Spec.map_comp, Spec.map_preimage, ← cancel_mono (chartA N ℤ i'), Category.assoc,
      specMap_toOverlap'_comp, hfac, hp]

/-! ## ★★★★★★★★★★★★Fubini–Study はチャートに依らない（仮定なしの形） -/

/-- ★**遷移単元は `x_{i'}/x_i` の像である**。 -/
theorem toOverlap_projCoord_self (N : ℕ) (i i' : Fin (N+1)) :
    toOverlap N i i' (projCoord N ℤ i i') = overlapUnit N i i' := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [projCoord, toOverlap, awayMap_mk, overlapUnit, Away.val_mk, Away.val_mk, pow_one]

/-- ★★★★★★★★★★★★**Fubini–Study はチャートに依らない**（重なりの仮定なし）。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★同じ複素点を 2 つのチャートで書いたとき、Fubini–Study の値は等しい。
★★非退化条件は **`φ(x_{i'}/x_i) ≠ 0`**——「点が `D₊(x_{i'})` にも入っている」ことの
座標での言い方である。 -/
theorem greenChartOf_congr (N : ℕ) (i i' : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (φ' : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i'))
      ⟶ CommRingCat.of ℂ)
    (hp : Spec.map φ ≫ chartA N ℤ i = Spec.map φ' ≫ chartA N ℤ i')
    (hne : φ.hom (projCoord N ℤ i i') ≠ 0) :
    greenChartOf N i φ = greenChartOf N i' φ' := by
  obtain ⟨χ, h1, h2⟩ := exists_overlap_factor N i i' φ φ' hp
  have hχ : χ.hom (overlapUnit N i i') ≠ 0 := by
    rw [← toOverlap_projCoord_self]
    intro hc
    refine hne ?_
    rw [← h1, CommRingCat.hom_comp, CommRingCat.hom_ofHom]
    exact hc
  rw [← h1, ← h2]
  exact greenChartOf_overlap N i i' χ hχ

/-! ## ★出典の紐付け(`.src`) -/

def opensRange_chartOverlap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(重なりのチャートの像は 2 つのチャートの像の共通部分)",
    sectionId := "genell-prop-1-4" }

def exists_overlap_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(2 つのチャートを通る点は重なりのチャートを通る)",
    sectionId := "genell-prop-1-4" }

def greenChartOf_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(Fubini–Study はチャートに依らない——重なりの仮定なし)",
    sectionId := "genell-prop-1-4" }

def greenChartOf_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "greenChartOf_overlap(代数の側、§9-869)"
      (.inProject "ABC3" "ABC3.Found.GenEll.greenChartOf_overlap") 2,
    .citation "[mathlib]" "Proj.SpecMap_awayMap_awayι(重なりのチャートは各チャートを経由)"
      (.inMathlib "AlgebraicGeometry.Proj.SpecMap_awayMap_awayι") 2,
    .citation "[mathlib]" "Proj.basicOpen_mul / IsOpenImmersion.lift"
      (.inMathlib "AlgebraicGeometry.Proj.basicOpen_mul") 2,
    .implicitStep
      ("★★★★★測定: 非退化条件は χ を持ち出さずに **φ(x_{i'}/x_i) ≠ 0** と書ける" ++
       "——遷移単元 t = x_{i'}/x_i は toOverlap (x_{i'}/x_i) そのものだからである" ++
       "(toOverlap_projCoord_self、2026-08-28 実測)。" ++
       "★これは「点が D₊(x_{i'}) にも入っている」ことの座標での言い方である") 3,
    .implicitStep
      ("★★大域的な green を書き下すのに残るのは、" ++
       "(1) 各複素点にチャートを選ぶ配管(exists_chart_range は体値の点について既にある)と " ++
       "(2) 得られた green の**連続性**(htArith_sub_abs_le が要求する ArcModel の位相)である") 4 ]

end ABC3.Found.GenEll
