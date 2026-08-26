/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeWeilOrd
import Mathlib.RingTheory.Valuation.ValuationRing

/-!
# `ord_x = 0` は「茎の単元」と同値(鎖 `normalize` の `B-functor` の最後の 1 段の前段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★何のためか

`B(L)` の関手性(`L → M` で `B(L) → B(M)`)を出すには、
「`ord_x(u) = 0`」を**茎の言葉**に翻訳しておく必要がある ——
そうすれば支配的な射 `g : V[M] → V[L]` に沿って
`𝒪_{V[L],g(x)} → 𝒪_{V[M],x}` が単元を単元へ送ることから、直ちに従う。

## ★★★中身は 2 本

| 向き | 根拠 |
|---|---|
| 単元 ⟹ `ord = 0` | `HeightOneSpectrum.valuation_eq_one_iff_notMem` ＋ `IsLocalRing.notMem_maximalIdeal` |
| `ord = 0` ⟹ 単元 | ★DVR は付値環(`ValuationRing.isInteger_or_isInteger`) |

★後者が本ファイルの実質である —— `v(x) ≤ 1` なら `x` は茎に入る、を
**DVR が付値環であること**から出す(`x` か `x⁻¹` のどちらかが整であり、
`x⁻¹` の側なら `v(x) = 1` から分母が単元になる)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `exists_algebraMap_of_valuation_le_one` | ★★`v(x) ≤ 1` ⟹ `x` は整 |
| `ordPt_eq_zero_iff_val` | `ord = 0 ⟺ v = 1` |
| `ordPt_eq_zero_of_isUnit` | ★**茎の単元は `ord = 0`** |
| `exists_unit_of_ordPt_eq_zero` | ★★**`ord = 0` なら茎の単元** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory IsDedekindDomain

universe u

/-! ## ★1. DVR: 付値 `≤ 1` なら整 -/

section DVR

variable {R : Type*} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type*} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★**`v(x) ≤ 1` なら `x` は `R` から来る** —— DVR は付値環である。 -/
theorem exists_algebraMap_of_valuation_le_one {x : K}
    (hx : (dvrSpectrum R).valuation K x ≤ 1) : ∃ r : R, algebraMap R K r = x := by
  rcases ValuationRing.isInteger_or_isInteger R x with ⟨r, hr⟩ | ⟨r, hr⟩
  · exact ⟨r, hr⟩
  rcases eq_or_ne x 0 with rfl | hx0
  · exact ⟨0, by simp⟩
  have hinv : (dvrSpectrum R).valuation K x⁻¹ = ((dvrSpectrum R).valuation K x)⁻¹ :=
    map_inv₀ _ _
  have hne : (dvrSpectrum R).valuation K x ≠ 0 := by
    simpa [Valuation.zero_iff] using hx0
  have hle' : ((dvrSpectrum R).valuation K x)⁻¹ ≤ 1 := by
    rw [← hinv, ← hr]
    exact HeightOneSpectrum.valuation_le_one _ r
  have h1 : (dvrSpectrum R).valuation K x = 1 := by
    refine le_antisymm hx ?_
    calc (1 : WithZero (Multiplicative ℤ))
        = (dvrSpectrum R).valuation K x * ((dvrSpectrum R).valuation K x)⁻¹ :=
          (mul_inv_cancel₀ hne).symm
      _ ≤ (dvrSpectrum R).valuation K x * 1 := mul_le_mul_right hle' _
      _ = (dvrSpectrum R).valuation K x := mul_one _
  have hru : (dvrSpectrum R).valuation K (algebraMap R K r) = 1 := by
    rw [hr, hinv, h1, inv_one]
  have hnot : r ∉ (dvrSpectrum R).asIdeal := by
    rw [← HeightOneSpectrum.valuation_eq_one_iff_notMem (K := K)]
    exact hru
  have hu : IsUnit r := IsLocalRing.notMem_maximalIdeal.mp hnot
  refine ⟨(hu.unit⁻¹ : Rˣ), ?_⟩
  have hmul : algebraMap R K ((hu.unit⁻¹ : Rˣ) : R) * algebraMap R K r = 1 := by
    rw [← map_mul]
    have : ((hu.unit⁻¹ : Rˣ) : R) * r = 1 := by
      simp
    rw [this, map_one]
  rw [hr] at hmul
  field_simp at hmul
  exact hmul

end DVR

/-! ## ★2. `ord = 0 ⟺ 茎の単元` -/

variable {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]

/-- ★`ord_x(u) = 0` は付値が `1` であることと同値。 -/
theorem ordPt_eq_zero_iff_val (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {u : X.functionField} (hu : u ≠ 0) :
    ordPt X hnorm x u = 0 ↔ ordPtVal X hnorm x u = 1 := by
  have hv : ordPtVal X hnorm x u ≠ 0 := ordPtVal_ne_zero hnorm x hu
  rw [ordPt, dif_neg hv, neg_eq_zero]
  have key : (WithZero.unzero hv).toAdd = 0 ↔ WithZero.unzero hv = 1 := by
    constructor
    · intro h
      have := congrArg Multiplicative.ofAdd h
      simpa using this
    · intro h
      rw [h]
      rfl
  rw [key, ← WithZero.coe_inj, WithZero.coe_unzero]
  exact Iff.rfl

/-- ★★**茎の単元は `ord_x = 0`**。 -/
theorem ordPt_eq_zero_of_isUnit (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (r : (X.presheaf.stalk x.1)ˣ) :
    ordPt X hnorm x
        (algebraMap (X.presheaf.stalk x.1) X.functionField (r : X.presheaf.stalk x.1)) = 0 := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  have hne : (algebraMap (X.presheaf.stalk x.1) X.functionField
      (r : X.presheaf.stalk x.1)) ≠ 0 := by
    simp only [ne_eq, map_eq_zero_iff _ (IsFractionRing.injective
      (X.presheaf.stalk x.1) X.functionField)]
    exact r.ne_zero
  rw [ordPt_eq_zero_iff_val hnorm x hne]
  show (dvrSpectrum (X.presheaf.stalk x.1)).valuation X.functionField
      (algebraMap _ _ (r : X.presheaf.stalk x.1)) = 1
  rw [HeightOneSpectrum.valuation_eq_one_iff_notMem (K := X.functionField)]
  exact fun hmem => (IsLocalRing.mem_maximalIdeal _).mp hmem r.isUnit

/-- ★★★**`ord_x = 0` なら茎の単元**。 -/
theorem exists_unit_of_ordPt_eq_zero (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {u : X.functionField} (hu : u ≠ 0) (h : ordPt X hnorm x u = 0) :
    ∃ r : (X.presheaf.stalk x.1)ˣ,
      algebraMap (X.presheaf.stalk x.1) X.functionField (r : X.presheaf.stalk x.1) = u := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  have hval : (dvrSpectrum (X.presheaf.stalk x.1)).valuation X.functionField u = 1 :=
    (ordPt_eq_zero_iff_val hnorm x hu).mp h
  obtain ⟨r, hr⟩ := exists_algebraMap_of_valuation_le_one (le_of_eq hval)
  have hru : (dvrSpectrum (X.presheaf.stalk x.1)).valuation X.functionField
      (algebraMap _ _ r) = 1 := by rw [hr]; exact hval
  have hnot : r ∉ (dvrSpectrum (X.presheaf.stalk x.1)).asIdeal := by
    rw [← HeightOneSpectrum.valuation_eq_one_iff_notMem (K := X.functionField)]
    exact hru
  have hunit : IsUnit r := IsLocalRing.notMem_maximalIdeal.mp hnot
  exact ⟨hunit.unit, by rw [hunit.unit_spec]; exact hr⟩

/-! ## ★3. `ord ≥ 0` は「茎に入る」と同値 -/

/-- ★`ord_x(f) ≥ 0` は付値が `≤ 1` であることと同値。 -/
theorem ordPt_nonneg_iff_val (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {u : X.functionField} (hu : u ≠ 0) :
    0 ≤ ordPt X hnorm x u ↔ ordPtVal X hnorm x u ≤ 1 := by
  have hv : ordPtVal X hnorm x u ≠ 0 := ordPtVal_ne_zero hnorm x hu
  rw [ordPt, dif_neg hv, withZero_le_one_iff_toAdd_nonpos _ hv, neg_nonneg]

/-- ★★**茎の元は `ord_x ≥ 0`**。 -/
theorem ordPt_nonneg_of_mem_stalk (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (r : X.presheaf.stalk x.1)
    (hr : algebraMap (X.presheaf.stalk x.1) X.functionField r ≠ 0) :
    0 ≤ ordPt X hnorm x (algebraMap (X.presheaf.stalk x.1) X.functionField r) := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  rw [ordPt_nonneg_iff_val hnorm x hr]
  show (dvrSpectrum (X.presheaf.stalk x.1)).valuation X.functionField (algebraMap _ _ r) ≤ 1
  exact HeightOneSpectrum.valuation_le_one _ r

/-- ★★★**`ord_x ≥ 0` なら茎から来る**。 -/
theorem exists_mem_stalk_of_ordPt_nonneg (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {u : X.functionField} (hu : u ≠ 0) (h : 0 ≤ ordPt X hnorm x u) :
    ∃ r : X.presheaf.stalk x.1,
      algebraMap (X.presheaf.stalk x.1) X.functionField r = u := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  exact exists_algebraMap_of_valuation_le_one ((ordPt_nonneg_iff_val hnorm x hu).mp h)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `B(L)` の条件を茎の言葉に翻訳したもの。 -/
def exists_unit_of_ordPt_eq_zero.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — ord_x(u) = 0 は茎の単元であることと同値",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
