/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeWeil
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.Spectrum.Prime.Topology
import Mathlib.AlgebraicGeometry.Stalk

/-!
# `ord_x` と アフィン開との両立(鎖 `weil` の `ord-pt` / `affine-compat`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> normal [geometrically integral] variety over a field k; K the function field of V ;

## ★★2 つのこと

1. **`ord_x : K(X)^× → ℤ`** —— 余次元 1 の点 `x` の茎は DVR
   (`isDiscreteValuationRing_stalk_of_codimOne`)で、関数体はその分数体
   (`isFractionRing_stalk`)だから、`Definition 3.1` の付値がそのまま取れる。
   ★DVR の高さ 1 スペクトルは**極大イデアルただ 1 つ**である。

2. **アフィン開との両立** —— 「余次元 1 の点」と「高さ 1 の素イデアル」の対応。
   ★`Spec R` の茎は局所化(`Spec.stalkIso`)、局所化の Krull 次元は高さ
   (`IsLocalization.AtPrime.ringKrullDim_eq_height`)。
   ★開埋め込みでは茎が変わらない(`f.stalkMap` が同型)ので、
   **どのアフィン開で見ても同じ**である。

★★**`SchemeWeil.lean` で訂正したとおり、1 は 2 に依存しない** ——
茎はそれ自身が局所 Noether 整域なので、アフィン開へ降りる必要が無い。
★2 が要るのは、**素因子を環の言葉で数える**とき(`div(f)` の台の有限性)である。

## ★本ファイルで閉じること

| 節点 | 定理 |
|---|---|
| `weil:ord-pt` | `ordPtVal` / `ordPt` / `ordPt_mul` / `ordPt_one` |
| `weil:affine-compat` | `isCodimOnePt_spec_iff` / `isCodimOnePt_iff_of_isOpenImmersion` |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory IsDedekindDomain

universe u

/-! ## ★1. DVR の高さ 1 スペクトル -/

/-- ★★**DVR の高さ 1 スペクトル** —— 極大イデアルただ 1 つ。 -/
def dvrSpectrum (R : Type*) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R] :
    HeightOneSpectrum R :=
  ⟨IsLocalRing.maximalIdeal R, inferInstance, IsDiscreteValuationRing.not_a_field R⟩

/-! ## ★2. `ord_x` -/

variable (X : Scheme.{u}) [IsIntegral X] [IsLocallyNoetherian X]

/-- ★★**余次元 1 の点が定める関数体の付値**。 -/
noncomputable def ordPtVal (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    Valuation X.functionField (WithZero (Multiplicative ℤ)) :=
  letI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  (dvrSpectrum (X.presheaf.stalk x.1)).valuation X.functionField

/-- ★★★**`ord_x : K(X) → ℤ`**(`0` では `0` と置く)。 -/
noncomputable def ordPt (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (f : X.functionField) : ℤ :=
  if h : ordPtVal X hnorm x f = 0 then 0 else -(WithZero.unzero h).toAdd

variable {X}

theorem ordPtVal_ne_zero (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {f : X.functionField} (hf : f ≠ 0) : ordPtVal X hnorm x f ≠ 0 := by
  simpa [Valuation.zero_iff] using hf

@[simp] theorem ordPt_zero (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ordPt X hnorm x 0 = 0 := by
  simp [ordPt]

@[simp] theorem ordPt_one (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ordPt X hnorm x 1 = 0 := by
  have h : ordPtVal X hnorm x 1 ≠ 0 := by simp
  rw [ordPt, dif_neg h, neg_eq_zero]
  have : WithZero.unzero h = 1 := by
    rw [← WithZero.coe_inj, WithZero.coe_unzero]
    simp
  rw [this]
  rfl

/-- ★★**`ord_x` は乗法を加法へ移す**。 -/
theorem ordPt_mul (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    {f g : X.functionField} (hf : f ≠ 0) (hg : g ≠ 0) :
    ordPt X hnorm x (f * g) = ordPt X hnorm x f + ordPt X hnorm x g := by
  have hvf := ordPtVal_ne_zero hnorm x hf
  have hvg := ordPtVal_ne_zero hnorm x hg
  have hvfg := ordPtVal_ne_zero hnorm x (mul_ne_zero hf hg)
  rw [ordPt, ordPt, ordPt, dif_neg hvf, dif_neg hvg, dif_neg hvfg, ← neg_add]
  refine congrArg Neg.neg ?_
  have hmul : WithZero.unzero hvfg = WithZero.unzero hvf * WithZero.unzero hvg := by
    rw [← WithZero.coe_inj, WithZero.coe_unzero, WithZero.coe_mul,
      WithZero.coe_unzero, WithZero.coe_unzero]
    exact map_mul _ f g
  rw [hmul]
  rfl

/-! ## ★3. アフィン開との両立(鎖 `weil` の `affine-compat`) -/

/-- ★★**`Spec R` の余次元 1 の点は、ちょうど高さ 1 の素イデアル**。

★茎は局所化(`Spec.stalkIso`)、局所化の Krull 次元は高さ
(`IsLocalization.AtPrime.ringKrullDim_eq_height`)。 -/
theorem isCodimOnePt_spec_iff (R : CommRingCat.{u}) (x : Spec R) :
    IsCodimOnePt (Spec R) x ↔ x.asIdeal.height = 1 := by
  have h1 : ringKrullDim ((Spec R).presheaf.stalk x)
      = ringKrullDim (Localization.AtPrime x.asIdeal) :=
    ((Spec.stalkIso R x).commRingCatIsoToRingEquiv).ringKrullDim
  have h2 : ringKrullDim (Localization.AtPrime x.asIdeal) = x.asIdeal.height :=
    IsLocalization.AtPrime.ringKrullDim_eq_height x.asIdeal (Localization.AtPrime x.asIdeal)
  rw [IsCodimOnePt, h1, h2]
  norm_cast

/-- ★★**開埋め込みでは余次元 1 性は変わらない** —— 茎が同型だから。

★★これが「どのアフィン開で見ても同じ」ということである。 -/
theorem isCodimOnePt_iff_of_isOpenImmersion {X Y : Scheme.{u}} (f : X ⟶ Y)
    [IsOpenImmersion f] (x : X) :
    IsCodimOnePt X x ↔ IsCodimOnePt Y (f.base x) := by
  have h : ringKrullDim (Y.presheaf.stalk (f.base x)) = ringKrullDim (X.presheaf.stalk x) :=
    ((asIso (f.stalkMap x)).commRingCatIsoToRingEquiv).ringKrullDim
  rw [IsCodimOnePt, IsCodimOnePt, ← h]

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `ord_x`。 -/
def ordPt.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 余次元 1 の点が定める ord_x : K(X)^× → ℤ",
    sectionId := "frdi-example-6-1" }

/-- ★★locator —— `Example 6.1` のアフィン開との両立(余次元 1 の点 ↔ 高さ 1 の素イデアル)。 -/
def isCodimOnePt_spec_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 余次元 1 の点と高さ 1 の素イデアルの対応",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
