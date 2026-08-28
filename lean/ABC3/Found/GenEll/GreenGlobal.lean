/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GreenChartGlue
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★大域的な Fubini–Study の Green 関数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5–6。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★★★★★★★★★★★これは何か —— `ArithCartier` の欄に入る `green`

`§9-869`（代数）と `§9-870`（スキーム）でチャート独立性が閉じたので、
★**チャートを 1 つ選んで定義すれば、選び方に依らない**。本ファイルはそれを行う:

    `greenFS N : complexPoints (ℙᴺ_ℤ) → ℝ`

★★これが `ArithCartier` の欄（`green : GreenFn X`）にそのまま入る形である。

## ★★★機構 —— 選択と、選び方に依らないこと

| 補題 | 内容 |
|---|---|
| `exists_chart_range`（`§9-C2b`） | 体値の点はどれかのチャートに入る |
| `specMap_projChartHom` | 選んだチャートで実際に分解する |
| `projCoord_ne_zero_of_range` | ★**別のチャートにも入る ⟹ `φ(x_{i'}/x_i) ≠ 0`** |
| `greenChartOf_congr`（`§9-870`） | ★★チャートに依らない |

★★★3 番目が本ファイルの新しい段である——「点が `D₊(x_{i'})` に入る」を
**座標の非零性**に翻訳する。機構は
`Proj.awayι_preimage_basicOpen`（チャートの上での基本開の逆像）と、
「体の点の基本開への所属は像が `0` でないこと」だけである。

## ★残っている段（明示）

★★★★`greenFS` の**連続性**（`htArith_sub_abs_le` が要求する `ArcModel` の位相）が残る。
★これは `sup` と `log` の連続性であり、`ProjTopology.lean` の測定が要る。
-/

namespace ABC3.Found.GenEll

open MvPolynomial HomogeneousLocalization AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★基本開への所属を座標の非零性に翻訳する -/

/-- ★**`isLocalizationElem` は正規化座標である**（1 次同士の場合）。 -/
theorem isLocalizationElem_eq_projCoord (N : ℕ) (i i' : Fin (N+1)) :
    HomogeneousLocalization.Away.isLocalizationElem
        (𝒜 := MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.isHomogeneous_X ℤ i) (MvPolynomial.isHomogeneous_X ℤ i')
      = projCoord N ℤ i i' := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [HomogeneousLocalization.Away.isLocalizationElem, projCoord,
    Away.val_mk, Away.val_mk, pow_one]

/-- ★**体の点が基本開に入るのは像が `0` でないことである**。 -/
theorem specMap_mem_basicOpen_iff {A : CommRingCat.{0}} (φ : A ⟶ CommRingCat.of ℂ) (a : A) :
    (Spec.map φ).base (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
        ∈ PrimeSpectrum.basicOpen a ↔ φ.hom a ≠ 0 := by
  show PrimeSpectrum.comap φ.hom (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ PrimeSpectrum.basicOpen a ↔ _
  rw [PrimeSpectrum.mem_basicOpen]
  show a ∉ Ideal.comap φ.hom (⊥ : Ideal (CommRingCat.of ℂ)) ↔ _
  rw [Ideal.mem_comap, Ideal.mem_bot]

/-- ★★★★★**別のチャートにも入るなら、その座標は `0` でない**。

★これが「点が `D₊(x_{i'})` に入る」ことの**座標での言い方**である。 -/
theorem projCoord_ne_zero_of_range (N : ℕ) (i i' : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ)
    (hx : Set.range (Spec.map φ ≫ chartA N ℤ i).base
      ⊆ Set.range (chartA N ℤ i').base) :
    φ.hom (projCoord N ℤ i i') ≠ 0 := by
  have hmem : (Spec.map φ ≫ chartA N ℤ i).base
      (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ Set.range (chartA N ℤ i').base := hx ⟨_, rfl⟩
  have hopen : (Spec.map φ ≫ chartA N ℤ i).base
      (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i') := by
    rw [← Proj.opensRange_awayι (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
      (MvPolynomial.X i') (MvPolynomial.isHomogeneous_X ℤ i') one_pos]
    exact hmem
  have hpre : (Spec.map φ).base (⟨⊥, Ideal.isPrime_bot⟩ : PrimeSpectrum (CommRingCat.of ℂ))
      ∈ (chartA N ℤ i) ⁻¹ᵁ (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)
        (MvPolynomial.X i')) := hopen
  rw [Proj.awayι_preimage_basicOpen _ (MvPolynomial.isHomogeneous_X ℤ i) one_pos
    (MvPolynomial.isHomogeneous_X ℤ i') one_pos, isLocalizationElem_eq_projCoord] at hpre
  exact (specMap_mem_basicOpen_iff φ _).1 hpre

/-! ## ★★選んだチャートでの分解 -/

/-- ★★**`projChartHom` は実際にその点を分解する**。 -/
theorem specMap_projChartHom (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (i : Fin (N + 1))
    (hx : Set.range x.base ⊆ Set.range (chartA N R i).base) :
    Spec.map (projChartHom N R F x i hx) ≫ chartA N R i = x := by
  rw [projChartHom, Spec.map_preimage]
  exact IsOpenImmersion.lift_fac _ _ hx

/-! ## ★★★★★★★★★★★★★大域的な Green 関数 -/

/-- ★**複素点に選ぶチャート**。 -/
noncomputable def chartIndexOf (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    Fin (N + 1) :=
  Classical.choose (exists_chart_range N ℤ ℂ p)

theorem chartIndexOf_spec (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) :
    Set.range p.base ⊆ Set.range (chartA N ℤ (chartIndexOf N p)).base :=
  Classical.choose_spec (exists_chart_range N ℤ ℂ p)

/-- ★★★★★★★★★★★★★**大域的な Fubini–Study の Green 関数**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★`ArithCartier` の欄（`green : GreenFn X`）にそのまま入る形である。 -/
noncomputable def greenFS (N : ℕ)
    (p : complexPoints (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))) : ℝ :=
  greenChartOf N (chartIndexOf N p)
    (projChartHom N ℤ ℂ p (chartIndexOf N p) (chartIndexOf_spec N p))

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★**`greenFS` はどのチャートで測っても同じ値である**。

★`§9-870` の `greenChartOf_congr` に、
`projCoord_ne_zero_of_range`（別のチャートにも入る ⟹ 座標が `0` でない）を渡すだけである。 -/
theorem greenFS_eq_greenChartOf (N : ℕ) (i : Fin (N+1))
    (φ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i))
      ⟶ CommRingCat.of ℂ) :
    greenFS N (Spec.map φ ≫ chartA N ℤ i) = greenChartOf N i φ := by
  set p := Spec.map φ ≫ chartA N ℤ i with hp
  set i₀ := chartIndexOf N p with hi₀
  set φ₀ := projChartHom N ℤ ℂ p i₀ (chartIndexOf_spec N p) with hφ₀
  have hfac : Spec.map φ₀ ≫ chartA N ℤ i₀ = p :=
    specMap_projChartHom N ℤ ℂ p i₀ (chartIndexOf_spec N p)
  have hsub : Set.range (Spec.map φ₀ ≫ chartA N ℤ i₀).base
      ⊆ Set.range (chartA N ℤ i).base := by
    rw [hfac, hp]
    rintro _ ⟨y, rfl⟩
    exact ⟨(Spec.map φ).base y, rfl⟩
  have hne : φ₀.hom (projCoord N ℤ i₀ i) ≠ 0 :=
    projCoord_ne_zero_of_range N i₀ i φ₀ hsub
  show greenChartOf N i₀ φ₀ = greenChartOf N i φ
  exact greenChartOf_congr N i₀ i φ₀ φ (by rw [hfac, hp]) hne

/-! ## ★★★★★★★★★★★★★★超平面の算術因子とその高さ -/

/-- ★★★★★★**超平面の算術因子** `(超平面, Fubini–Study)`。 -/
noncomputable def hyperplaneArith (N : ℕ) :
    ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)) :=
  { divisor := hyperplaneIdeal N ℤ
    green := greenFS N }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★**超平面の算術因子の高さは素朴高さである**（仮定なし）。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `htArith F (hyperplaneArith N) (Spec ψ ≫ chartA i₀) = log H(x) / [F : ℚ]`

★`§9-868` の `htArith_hyperplane_eq_log_mulHeight` から仮定 `hdiv`・`hgreenChart` が
**どちらも消えた**形である。★★`hgreenChart` は
`archPoint_chart_factor`（複素点は同じチャートを通る）と
`greenFS_eq_greenChartOf`（`greenFS` はチャートに依らない）で潰れる。 -/
theorem htArith_hyperplaneArith (F : Type) [Field F] [NumberField F]
    (N : ℕ) (i₀ : Fin (N+1))
    (ψ : CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
      ⟶ CommRingCat.of (NumberField.RingOfIntegers F))
    (x : Fin (N+1) → NumberField.RingOfIntegers F) (hx : x ≠ 0)
    (hw : ∀ k, x k = ψ.hom (projCoord N ℤ i₀ k) * x i₀)
    (h0 : x 0 ≠ 0) :
    htArith F (hyperplaneArith N) (Spec.map ψ ≫ chartA N ℤ i₀)
      = Real.log (Height.mulHeight (fun k => ((x k : F)))) / (Module.finrank ℚ F : ℝ) := by
  refine htArith_hyperplane_eq_log_mulHeight F N i₀ ψ (hyperplaneArith N) rfl x hx hw h0 ?_
  intro v
  show greenFS N (archPoint (Spec.map ψ ≫ chartA N ℤ i₀) v) = _
  rw [archPoint_chart_factor F N i₀ ψ v, greenFS_eq_greenChartOf]

/-! ## ★出典の紐付け(`.src`) -/

def hyperplaneArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(超平面の算術因子(超平面, Fubini–Study))",
    sectionId := "genell-prop-1-4" }

def htArith_hyperplaneArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(超平面の算術因子の高さは素朴高さである——仮定なし)",
    sectionId := "genell-prop-1-4" }

def projCoord_ne_zero_of_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(別のチャートにも入るなら、その座標は 0 でない)",
    sectionId := "genell-prop-1-4" }

def specMap_projChartHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(projChartHom は実際にその点を分解する)",
    sectionId := "genell-prop-1-4" }

def greenFS.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(大域的な Fubini–Study の Green 関数)",
    sectionId := "genell-prop-1-4" }

def greenFS_eq_greenChartOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Proposition 1.4(greenFS はどのチャートで測っても同じ値である)",
    sectionId := "genell-prop-1-4" }

def greenFS_eq_greenChartOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "greenChartOf_congr(チャート独立性、§9-870)"
      (.inProject "ABC3" "ABC3.Found.GenEll.greenChartOf_congr") 2,
    .citation "[ABC3]" "exists_chart_range(体値の点はどれかのチャートに入る、§9-C2b)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_chart_range") 2,
    .citation "[mathlib]" "Proj.awayι_preimage_basicOpen(チャートの上での基本開の逆像)"
      (.inMathlib "AlgebraicGeometry.Proj.awayι_preimage_basicOpen") 2,
    .implicitStep
      ("★本ファイルの新しい段は「点が D₊(x_{i'}) に入る」を**座標の非零性**に" ++
       "翻訳することである(projCoord_ne_zero_of_range)。機構は " ++
       "Proj.awayι_preimage_basicOpen と「体の点の基本開への所属は像が 0 でないこと」だけ") 3,
    .implicitStep
      ("★★残るのは greenFS の**連続性**(htArith_sub_abs_le が要求する ArcModel の位相)である。" ++
       "sup と log の連続性であり、ProjTopology.lean の測定が要る") 5 ]

end ABC3.Found.GenEll
