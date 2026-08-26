/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeOrdUnit
import ABC3.Found.Divisor.HeightOneBridge
import Mathlib.AlgebraicGeometry.Noetherian

/-!
# `div(f)` の台は有限(鎖 `weil` の `div-finite` —— スキームの層)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> that map into [possibly subvarieties of codimension ≥1 of] prime divisors of DK.

## ★★★環の層の結果をスキームへ運ぶ 3 段

| 段 | 根拠 |
|---|---|
| アフィン開の素因子 ⟷ 高さ 1 の素イデアル | `isLocalization_stalk` ＋ `ringKrullDim_eq_height` |
| `f = a/b` で `a, b ∉ 𝔭_v` なら `ord_v f = 0` | ★茎で単元(`ordPt_eq_zero_of_isUnit`) |
| `a` を含む高さ 1 素イデアルは有限個 | `finite_heightOne_primes_containing`(在庫) |

★★**茎の付値と環の付値を同一視する必要は無い** ——
「`ord = 0` ⟺ 茎の単元」(`SchemeOrdUnit.lean`)を経由すれば、
`ordPt` の定義を開かずに済む。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `primeIdealOf_height_eq_one` | 素因子の素イデアルは高さ 1 |
| `ordPt_eq_zero_of_notMem` | ★★`a, b ∉ 𝔭_v` なら `ord_v (a/b) = 0` |
| `finite_ordPt_ne_zero_on_affine` | ★★★**アフィン開の中では台は有限** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory

universe u

variable {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]

/-! ## ★1. 素因子の素イデアルは高さ 1 -/

omit [AlgebraicGeometry.IsIntegral X] [IsLocallyNoetherian X] in
/-- ★**アフィン開の中の素因子は、切断環の高さ 1 の素イデアル**。 -/
theorem primeIdealOf_height_eq_one {U : X.Opens} (hU : IsAffineOpen U) [Nonempty U]
    (v : PrimeDivisorPt X) (hv : (v : X) ∈ U) :
    (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal.height = 1 := by
  letI : Algebra Γ(X, U) (X.presheaf.stalk ((⟨v.1, hv⟩ : U) : X)) :=
    TopCat.Presheaf.algebra_section_stalk X.presheaf ⟨v.1, hv⟩
  haveI := hU.isLocalization_stalk ⟨v.1, hv⟩
  have h := IsLocalization.AtPrime.ringKrullDim_eq_height
    (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal (X.presheaf.stalk ((⟨v.1, hv⟩ : U) : X))
  have hv1 : ringKrullDim (X.presheaf.stalk v.1) = 1 := v.2
  rw [hv1] at h
  exact_mod_cast h.symm

/-! ## ★2. `a, b ∉ 𝔭_v` なら `ord_v (a/b) = 0` -/

/-- ★★**分子も分母も `𝔭_v` に入らなければ `ord_v = 0`** —— 茎で単元になるから。 -/
theorem ordPt_eq_zero_of_notMem (hnorm : IsNormalScheme X) {U : X.Opens} (hU : IsAffineOpen U)
    [Nonempty U] (v : PrimeDivisorPt X) (hv : (v : X) ∈ U) {a b : Γ(X, U)}
    (ha' : a ∉ (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal)
    (hb' : b ∉ (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal)
    {f : X.functionField}
    (hf : f * algebraMap Γ(X, U) X.functionField b
      = algebraMap Γ(X, U) X.functionField a) :
    ordPt X hnorm v f = 0 := by
  letI : Algebra Γ(X, U) (X.presheaf.stalk ((⟨v.1, hv⟩ : U) : X)) :=
    TopCat.Presheaf.algebra_section_stalk X.presheaf ⟨v.1, hv⟩
  haveI := hU.isLocalization_stalk ⟨v.1, hv⟩
  haveI := AlgebraicGeometry.functionField_isScalarTower (X := X) U ⟨v.1, hv⟩
  obtain ⟨ua, hua⟩ : IsUnit (algebraMap Γ(X, U) (X.presheaf.stalk v.1) a) :=
    IsLocalization.map_units (M := (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal.primeCompl) _ ⟨a, ha'⟩
  obtain ⟨ub, hub⟩ : IsUnit (algebraMap Γ(X, U) (X.presheaf.stalk v.1) b) :=
    IsLocalization.map_units (M := (hU.primeIdealOf ⟨v.1, hv⟩).asIdeal.primeCompl) _ ⟨b, hb'⟩
  have htower : ∀ y : Γ(X, U),
      algebraMap (X.presheaf.stalk v.1) X.functionField
        (algebraMap Γ(X, U) (X.presheaf.stalk v.1) y)
      = algebraMap Γ(X, U) X.functionField y := fun y =>
    (IsScalarTower.algebraMap_apply Γ(X, U) (X.presheaf.stalk v.1) X.functionField y).symm
  have hbne : algebraMap Γ(X, U) X.functionField b ≠ 0 := by
    rw [← htower b, map_ne_zero_iff _
      (IsFractionRing.injective (X.presheaf.stalk v.1) X.functionField), ← hub]
    exact ub.ne_zero
  have hmul : ((ua * ub⁻¹ : (X.presheaf.stalk v.1)ˣ) : X.presheaf.stalk v.1)
      * algebraMap Γ(X, U) (X.presheaf.stalk v.1) b
      = algebraMap Γ(X, U) (X.presheaf.stalk v.1) a := by
    rw [← hua, ← hub, Units.val_mul, mul_assoc, ub.inv_mul, mul_one]
  have hkey : algebraMap (X.presheaf.stalk v.1) X.functionField
      ((ua * ub⁻¹ : (X.presheaf.stalk v.1)ˣ) : X.presheaf.stalk v.1) = f := by
    have hc := congrArg (algebraMap (X.presheaf.stalk v.1) X.functionField) hmul
    rw [map_mul, htower a, htower b, ← hf] at hc
    exact mul_right_cancel₀ hbne hc
  rw [← hkey]
  exact ordPt_eq_zero_of_isUnit hnorm v (ua * ub⁻¹)

/-! ## ★3. アフィン開の中では台は有限 -/

open scoped Classical in
/-- ★★★**アフィン開の中では `ord_v f ≠ 0` となる素因子は有限個**。 -/
theorem finite_ordPt_ne_zero_on_affine (hnorm : IsNormalScheme X) {U : X.Opens}
    (hU : IsAffineOpen U) [Nonempty U] {f : X.functionField} (hf : f ≠ 0) :
    {v : PrimeDivisorPt X | (v : X) ∈ U ∧ ordPt X hnorm v f ≠ 0}.Finite := by
  haveI : IsNoetherianRing Γ(X, U) := IsLocallyNoetherian.component_noetherian ⟨U, hU⟩
  haveI := functionField_isFractionRing_of_isAffineOpen (X := X) U hU
  obtain ⟨⟨a, b⟩, hab⟩ := IsLocalization.surj (nonZeroDivisors Γ(X, U)) f
  have hbne : algebraMap Γ(X, U) X.functionField (b : Γ(X, U)) ≠ 0 := by
    rw [map_ne_zero_iff _ (IsFractionRing.injective Γ(X, U) X.functionField)]
    exact nonZeroDivisors.coe_ne_zero b
  have hane : a ≠ 0 := by
    intro h
    rw [h, map_zero] at hab
    exact (mul_ne_zero hf hbne) hab
  have hbne0 : (b : Γ(X, U)) ≠ 0 := nonZeroDivisors.coe_ne_zero b
  -- 素イデアルの側の有限集合
  set P : Set (Ideal Γ(X, U)) :=
    {p | p.IsPrime ∧ (∀ q : Ideal Γ(X, U), q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) ∧ a ∈ p} ∪
    {p | p.IsPrime ∧ (∀ q : Ideal Γ(X, U), q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) ∧ (b : Γ(X, U)) ∈ p}
    with hP
  have hPfin : P.Finite :=
    (finite_heightOne_primes_containing hane).union
      (finite_heightOne_primes_containing hbne0)
  have hP'fin : {p : PrimeSpectrum Γ(X, U) | p.asIdeal ∈ P}.Finite :=
    Set.Finite.preimage (fun _ _ _ _ h => PrimeSpectrum.ext h) hPfin
  refine Set.Finite.of_finite_image ?_ (Subtype.val_injective.injOn)
  refine Set.Finite.subset (hP'fin.image (fun p => (hU.fromSpec).base p)) ?_
  rintro _ ⟨v, ⟨hvU, hord⟩, rfl⟩
  refine ⟨hU.primeIdealOf ⟨v.1, hvU⟩, ?_, hU.fromSpec_primeIdealOf ⟨v.1, hvU⟩⟩
  have hh : (hU.primeIdealOf ⟨v.1, hvU⟩).asIdeal.height = 1 :=
    primeIdealOf_height_eq_one hU v hvU
  have hmin := (height_eq_one_iff (hU.primeIdealOf ⟨v.1, hvU⟩).asIdeal).mp hh
  by_cases hin : a ∈ (hU.primeIdealOf ⟨v.1, hvU⟩).asIdeal
  · exact Or.inl ⟨(hU.primeIdealOf ⟨v.1, hvU⟩).isPrime, hmin.2, hin⟩
  · by_cases hin2 : (b : Γ(X, U)) ∈ (hU.primeIdealOf ⟨v.1, hvU⟩).asIdeal
    · exact Or.inr ⟨(hU.primeIdealOf ⟨v.1, hvU⟩).isPrime, hmin.2, hin2⟩
    · exact absurd (ordPt_eq_zero_of_notMem hnorm hU v hvU hin hin2 hab) hord

/-! ## ★4. 大域(準コンパクト)—— 台は有限 -/

/-- ★★★★**`div(f)` の台は有限**(`X` が準コンパクトなとき)。 -/
theorem finite_ordPt_ne_zero (hnorm : IsNormalScheme X) [CompactSpace X]
    {f : X.functionField} (hf : f ≠ 0) :
    {v : PrimeDivisorPt X | ordPt X hnorm v f ≠ 0}.Finite := by
  obtain ⟨s, hs⟩ := isCompact_univ.elim_finite_subcover
    (fun U : X.affineOpens => (U.1 : Set X)) (fun U => U.1.2) (by
      intro x _
      obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ := X.isBasis_affineOpens.exists_subset_of_mem_open
        (Set.mem_univ x) isOpen_univ
      exact Set.mem_iUnion.mpr ⟨⟨U, hU⟩, hxU⟩)
  refine Set.Finite.subset (Set.Finite.biUnion s.finite_toSet
    (fun U _ => show {v : PrimeDivisorPt X | (v : X) ∈ U.1 ∧ ordPt X hnorm v f ≠ 0}.Finite from ?_))
    ?_
  · by_cases hne : Nonempty U.1
    · exact finite_ordPt_ne_zero_on_affine hnorm U.2 hf
    · refine Set.Finite.subset (Set.finite_empty) ?_
      rintro v ⟨hvU, -⟩
      exact absurd ⟨⟨v.1, hvU⟩⟩ hne
  · intro v hv
    have : (v : X) ∈ (Set.univ : Set X) := Set.mem_univ _
    obtain ⟨U, hU⟩ := Set.mem_iUnion.mp (hs this)
    obtain ⟨hUs, hvU⟩ := Set.mem_iUnion.mp hU
    exact Set.mem_biUnion hUs ⟨hvU, hv⟩

open scoped Classical in
/-- ★★★★★**`div(f) : WeilDiv X`**(一般のスキーム、準コンパクト)。

★これで鎖 `weil` の `weil-group` がアフィンの場合を越える。 -/
noncomputable def weilDivOfFn (hnorm : IsNormalScheme X) [CompactSpace X]
    {f : X.functionField} (hf : f ≠ 0) : WeilDiv X :=
  Finsupp.onFinset (finite_ordPt_ne_zero hnorm hf).toFinset
    (fun v => ordPt X hnorm v f)
    (fun v hv => by simpa using hv)

@[simp] theorem weilDivOfFn_apply (hnorm : IsNormalScheme X) [CompactSpace X]
    {f : X.functionField} (hf : f ≠ 0) (v : PrimeDivisorPt X) :
    weilDivOfFn hnorm hf v = ordPt X hnorm v f := rfl

/-- ★★★**`div` は乗法を加法へ移す**。 -/
theorem weilDivOfFn_mul (hnorm : IsNormalScheme X) [CompactSpace X]
    {f g : X.functionField} (hf : f ≠ 0) (hg : g ≠ 0) :
    weilDivOfFn hnorm (mul_ne_zero hf hg) = weilDivOfFn hnorm hf + weilDivOfFn hnorm hg := by
  refine Finsupp.ext fun v => ?_
  rw [weilDivOfFn_apply, Finsupp.add_apply, weilDivOfFn_apply, weilDivOfFn_apply,
    ordPt_mul hnorm v hf hg]

/-- ★**`div` は群準同型**(単元群から)。 -/
noncomputable def weilDivHom (hnorm : IsNormalScheme X) [CompactSpace X] :
    Additive ((X.functionField)ˣ) →+ WeilDiv X where
  toFun u := weilDivOfFn hnorm (Additive.toMul u).ne_zero
  map_zero' := by
    refine Finsupp.ext fun v => ?_
    show ordPt X hnorm v ((1 : (X.functionField)ˣ) : X.functionField) = 0
    exact ordPt_one hnorm v
  map_add' u w := by
    refine Finsupp.ext fun v => ?_
    show ordPt X hnorm v
        (((Additive.toMul u * Additive.toMul w : (X.functionField)ˣ)) : X.functionField)
      = ordPt X hnorm v ((Additive.toMul u : (X.functionField)ˣ) : X.functionField)
        + ordPt X hnorm v ((Additive.toMul w : (X.functionField)ˣ) : X.functionField)
    exact ordPt_mul hnorm v (Additive.toMul u).ne_zero (Additive.toMul w).ne_zero

/-- ★★★locator —— `Example 6.1` の `div : K(X)^× → WeilDiv X`。 -/
def weilDivHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — div : K(X)^× → WeilDiv X(一般のスキーム)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
