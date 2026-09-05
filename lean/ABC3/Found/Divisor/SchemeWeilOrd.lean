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
| `weil:ord-nondeg` | ★`exists_ordPt_eq_one` / `exists_ordPt_eq_one'` / `exists_ordPt_eq` / `ordPt_surjective` / `not_forall_ordPt_eq_zero` |
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

/-! ## ★4. 錨 —— `ord_x` は零写像ではない(鎖 `weil` の `ord-nondeg`)

★★**なぜこの節が要るか** —— `ord_x` の性質として書かれるもの
(`ordPt_mul` / `ordPt_one` / `div(f)` の台の有限性)は、
**`ord_x ≡ 0` と置けばすべて自明に成り立ってしまう**。
`IsNormalScheme` は「正直な定義が可能になる」条件であって、
「零写像を排除する」条件ではない。排除するには**錨**が要る。

★錨は一様化元(uniformizer)である —— 余次元 1 の茎は DVR なので、
`HeightOneSpectrum.valuation_exists_uniformizer` が付値 `exp (-1)` の元を与える。
符号規約 `ordPt := -(unzero v f).toAdd` により `-(-1) = 1`。
-/

/-- ★**`ordPt` は `-log ∘ ordPtVal`** —— `0` の場合も含めて場合分け無しで書ける。

★`WithZero.log 0 = 0` なので、`ordPt` の `if` と一致する。 -/
theorem ordPt_eq_neg_log (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
    (f : X.functionField) :
    ordPt X hnorm x f = -WithZero.log (ordPtVal X hnorm x f) := by
  rw [ordPt]
  split
  · next h => rw [h]; simp
  · next h =>
      refine congrArg Neg.neg ?_
      conv_rhs => rw [← WithZero.coe_unzero h]
      rfl

/-- ★★**一様化元** —— 付値がちょうど `exp (-1)` になる関数体の元。

★`IsDedekindDomain.HeightOneSpectrum.valuation_exists_uniformizer` をそのまま乗せた。 -/
theorem exists_uniformizer_ordPtVal (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ∃ π : X.functionField, π ≠ 0 ∧ ordPtVal X hnorm x π = WithZero.exp (-1 : ℤ) := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  obtain ⟨π, hπ⟩ :=
    (dvrSpectrum (X.presheaf.stalk x.1)).valuation_exists_uniformizer X.functionField
  have hval : ordPtVal X hnorm x π = WithZero.exp (-1 : ℤ) := hπ
  refine ⟨π, ?_, hval⟩
  rintro rfl
  rw [map_zero] at hval
  exact WithZero.exp_ne_zero hval.symm

/-- ★★★**錨(一般形)** —— 任意の `n : ℤ` に対し `ord_x(f) = n` となる `f ∈ K(X)^×` がある。 -/
theorem exists_ordPt_eq (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) (n : ℤ) :
    ∃ f : (X.functionField)ˣ, ordPt X hnorm x (f : X.functionField) = n := by
  obtain ⟨π, hπ0, hval⟩ := exists_uniformizer_ordPtVal hnorm x
  refine ⟨Units.mk0 π hπ0 ^ n, ?_⟩
  have hcoe : ((Units.mk0 π hπ0 ^ n : (X.functionField)ˣ) : X.functionField) = π ^ n := by
    simp
  rw [hcoe, ordPt_eq_neg_log]
  have hz : ordPtVal X hnorm x (π ^ n) = WithZero.exp (-1 : ℤ) ^ n := by
    rw [map_zpow₀, hval]
  rw [hz, WithZero.log_zpow, WithZero.log_exp]
  simp

/-- ★★★★★**錨** —— `ord_x(f) = 1` となる `f ∈ K(X)^×` が存在する。

★★これ 1 本で `ordPt ≡ 0` は死ぬ(`not_forall_ordPt_eq_zero`)。 -/
theorem exists_ordPt_eq_one (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ∃ f : (X.functionField)ˣ, ordPt X hnorm x (f : X.functionField) = 1 :=
  exists_ordPt_eq hnorm x 1

/-- ★★★★★**錨(単元を剥いた形)** —— `Check/FrdI/Ex61OrdDegenerate.lean` の
`HasOrdAnchor` がそのまま消費できる形。 -/
theorem exists_ordPt_eq_one' (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ∃ f : X.functionField, f ≠ 0 ∧ ordPt X hnorm x f = 1 := by
  obtain ⟨f, hf⟩ := exists_ordPt_eq_one hnorm x
  exact ⟨(f : X.functionField), f.ne_zero, hf⟩

/-- ★★**錨(茎の側)** —— 一様化元は茎(正則関数)から取れる。

★これは `ordPt ≥ 0` の側に錨を打つ —— `div(f)` が実際に `x` を係数 `1` で含む。 -/
theorem exists_stalk_ordPt_eq_one (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ∃ r : X.presheaf.stalk x.1,
      algebraMap (X.presheaf.stalk x.1) X.functionField r ≠ 0 ∧
        ordPt X hnorm x (algebraMap (X.presheaf.stalk x.1) X.functionField r) = 1 := by
  haveI := isDiscreteValuationRing_stalk_of_codimOne X hnorm x
  obtain ⟨r, hr⟩ :=
    (dvrSpectrum (X.presheaf.stalk x.1)).valuation_exists_uniformizer' X.functionField
  have hval : ordPtVal X hnorm x (algebraMap _ _ r) = WithZero.exp (-1 : ℤ) := hr
  refine ⟨r, ?_, ?_⟩
  · intro h
    rw [h, map_zero] at hval
    exact WithZero.exp_ne_zero hval.symm
  · rw [ordPt_eq_neg_log, hval, WithZero.log_exp]
    simp

/-- ★★★**`ord_x : K(X)^× → ℤ` は全射**。 -/
theorem ordPt_surjective (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    Function.Surjective fun f : (X.functionField)ˣ => ordPt X hnorm x (f : X.functionField) :=
  fun n => exists_ordPt_eq hnorm x n

/-- ★★★★★**退化の否定** —— `ord_x` は零写像ではない。

★★`Check/**/*Degenerate.lean` の 7 例と同じ罠(条件を落とすと主張が自明になる)を、
ここで**塞いだ**ことの記録である。 -/
theorem not_forall_ordPt_eq_zero (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X) :
    ¬ ∀ f : X.functionField, ordPt X hnorm x f = 0 := by
  intro h
  obtain ⟨f, hf⟩ := exists_ordPt_eq_one hnorm x
  rw [h (f : X.functionField)] at hf
  exact zero_ne_one hf

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `ord_x`。 -/
def ordPt.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 余次元 1 の点が定める ord_x : K(X)^× → ℤ",
    sectionId := "frdi-example-6-1" }

/-- ★★★★★locator —— `ord_x` が零写像でないこと(錨)。

★★★逸脱の記録: **原典はこの主張を明示的に書いていない**。
FrdI p.109 の `Example 6.1` は `ord_x` を「余次元 1 の点の付値」として導入するだけで、
「零写像でない」ことは自明として畳んでいる。
本ファイルではこれを**明示の定理として立てた** —— 理由は、
`ord_x` の性質(乗法性・台の有限性)だけでは `ord_x ≡ 0` が反例として残り、
後続の `div(f)` / `B(L)` の議論が空になるからである。
主張の追加であって前提の追加ではないので、後続の証明への影響は無い。 -/
def exists_ordPt_eq_one.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — ord_x は零写像でない(一様化元による錨)",
    sectionId := "frdi-example-6-1" }

/-- ★★locator —— `Example 6.1` のアフィン開との両立(余次元 1 の点 ↔ 高さ 1 の素イデアル)。 -/
def isCodimOnePt_spec_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 余次元 1 の点と高さ 1 の素イデアルの対応",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
