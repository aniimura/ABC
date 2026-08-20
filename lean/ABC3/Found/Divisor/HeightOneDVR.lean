/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.DedekindDomain.Dvr
import Mathlib.RingTheory.DiscreteValuationRing.TFAE
import Mathlib.RingTheory.LocalProperties.IntegrallyClosed
import Mathlib.RingTheory.DedekindDomain.AdicValuation

/-!
# 正規 Noether 整域の余次元 1 の点は DVR

★★★**正規多様体の因子論の礎石**である。これがあって初めて余次元 1 の点での
**位数 `ord`** が定義でき、有理関数の因子 `div(f)` が書ける。

## ★なぜ自分で作るのか(2026-08-20 の在庫測定)

mathlib には **Dedekind 環(Krull 次元 ≤ 1)の場合しか無い**
(`IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain`)。
`Krull 整域`(`IsKrullDomain`)も無い。
★[FrdI] `Example 6.1` は**任意次元の正規多様体**を扱うので、
高さ 1 の素イデアルでの局所化が DVR であることを一般に要する。

## ★中身は 3 行

`IsDiscreteValuationRing.TFAE` の第 4 項は
**「整閉 ∧ 非零素イデアルがちょうど 1 つ」**である。したがって:

* 整閉性は局所化で保たれる(`isIntegrallyClosed_of_isLocalization`)
* 非零素イデアルの一意性は「`p` が**極小非零素**」から素イデアルの対応
  (`IsLocalization.AtPrime.orderIsoOfPrime`)で出る

★「高さ 1」を `Ideal.height` ではなく**極小非零素**という形で受け取る
(整域では同値で、`height` の API 摩擦を避けられる)。
-/

namespace ABC3.Found.Divisor

/-- ★★★★**正規 Noether 整域の極小非零素イデアルでの局所化は DVR**。

★`hmin` は「`p` の下にある非零素イデアルは `p` だけ」、すなわち **`p` の高さが 1**。 -/
theorem isDiscreteValuationRing_atPrime_of_minimal
    {R : Type*} [CommRing R] [IsDomain R] [IsNoetherianRing R] [IsIntegrallyClosed R]
    (p : Ideal R) [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hmin : ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p)
    (S : Type*) [CommRing S] [IsDomain S] [Algebra R S] [IsLocalization.AtPrime S p] :
    IsDiscreteValuationRing S := by
  haveI : IsNoetherianRing S := IsLocalization.isNoetherianRing p.primeCompl S inferInstance
  haveI : IsLocalRing S := IsLocalization.AtPrime.isLocalRing S p
  haveI : IsIntegrallyClosed S :=
    isIntegrallyClosed_of_isLocalization S p.primeCompl p.primeCompl_le_nonZeroDivisors
  have hnf : ¬ IsField S := IsLocalization.AtPrime.not_isField R hp0 S
  have htfae := (IsDiscreteValuationRing.TFAE S hnf).out 3 0
  have huniq : ∃! P : Ideal S, P ≠ ⊥ ∧ P.IsPrime := by
    refine ⟨IsLocalRing.maximalIdeal S,
      ⟨IsLocalRing.isField_iff_maximalIdeal_eq.not.mp hnf, inferInstance⟩, ?_⟩
    rintro P ⟨hP0, hPp⟩
    have hcomap : P.comap (algebraMap R S) ≠ ⊥ := by
      intro hc
      apply hP0
      have h1 : ((IsLocalization.AtPrime.orderIsoOfPrime S p) ⟨P, hPp⟩ : Ideal R)
          = ((IsLocalization.AtPrime.orderIsoOfPrime S p) ⟨⊥, Ideal.isPrime_bot⟩ : Ideal R) := by
        show P.comap (algebraMap R S) = (⊥ : Ideal S).comap (algebraMap R S)
        rw [hc, Ideal.comap_bot_of_injective _
          (IsLocalization.injective S p.primeCompl_le_nonZeroDivisors)]
      have h2 : (⟨P, hPp⟩ : {q : Ideal S // q.IsPrime}) = ⟨⊥, Ideal.isPrime_bot⟩ :=
        (IsLocalization.AtPrime.orderIsoOfPrime S p).injective (Subtype.ext h1)
      exact congrArg Subtype.val h2
    have hle : P.comap (algebraMap R S) ≤ p :=
      ((IsLocalization.AtPrime.orderIsoOfPrime S p) ⟨P, hPp⟩).2.2
    have heq : P.comap (algebraMap R S) = p :=
      hmin _ (Ideal.IsPrime.comap _) hcomap hle
    have hmax : (IsLocalRing.maximalIdeal S).comap (algebraMap R S) = p :=
      IsLocalization.AtPrime.under_maximalIdeal S p
    have h3 : (⟨P, hPp⟩ : {q : Ideal S // q.IsPrime})
        = ⟨IsLocalRing.maximalIdeal S, inferInstance⟩ :=
      (IsLocalization.AtPrime.orderIsoOfPrime S p).injective (Subtype.ext (by
        show P.comap (algebraMap R S) = (IsLocalRing.maximalIdeal S).comap (algebraMap R S)
        rw [heq, hmax]))
    exact congrArg Subtype.val h3
  exact htfae.mp ⟨inferInstance, huniq⟩

/-- ★`Localization.AtPrime` の形での言い換え。 -/
theorem isDiscreteValuationRing_localization_atPrime_of_minimal
    {R : Type*} [CommRing R] [IsDomain R] [IsNoetherianRing R] [IsIntegrallyClosed R]
    (p : Ideal R) [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hmin : ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) :
    IsDiscreteValuationRing (Localization.AtPrime p) :=
  isDiscreteValuationRing_atPrime_of_minimal p hp0 hmin _

/-! ## ★余次元 1 の点での位数 `ord`

★★DVR になったので、**Dedekind 環の adic 付値の在庫がそのまま使える**。
`Localization.AtPrime p` は DVR ゆえ Dedekind 環であり、
その唯一の高さ 1 素点(極大イデアル)の adic 付値を分数体へ延ばせばよい。 -/

open IsDedekindDomain in
/-- ★★★**余次元 1 の点での位数** `ord_p`。 -/
noncomputable def ordAt {R : Type*} [CommRing R] [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] (p : Ideal R) [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hmin : ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) :
    Valuation (FractionRing R) (WithZero (Multiplicative ℤ)) := by
  haveI hdvr : IsDiscreteValuationRing (Localization.AtPrime p) :=
    isDiscreteValuationRing_localization_atPrime_of_minimal p hp0 hmin
  haveI : IsFractionRing (Localization.AtPrime p) (FractionRing R) :=
    IsFractionRing.isFractionRing_of_isLocalization p.primeCompl _ _
      p.primeCompl_le_nonZeroDivisors
  have hnf : ¬ IsField (Localization.AtPrime p) :=
    IsLocalization.AtPrime.not_isField R hp0 _
  exact HeightOneSpectrum.valuation (FractionRing R)
    ⟨IsLocalRing.maximalIdeal (Localization.AtPrime p), inferInstance,
      IsLocalRing.isField_iff_maximalIdeal_eq.not.mp hnf⟩

/-- ★環の元の位数は非負(付値は 1 以下)。 -/
theorem ordAt_le_one {R : Type*} [CommRing R] [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] (p : Ideal R) [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hmin : ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) (r : R) :
    ordAt p hp0 hmin (algebraMap R (FractionRing R) r) ≤ 1 := by
  haveI hdvr : IsDiscreteValuationRing (Localization.AtPrime p) :=
    isDiscreteValuationRing_localization_atPrime_of_minimal p hp0 hmin
  haveI : IsFractionRing (Localization.AtPrime p) (FractionRing R) :=
    IsFractionRing.isFractionRing_of_isLocalization p.primeCompl _ _
      p.primeCompl_le_nonZeroDivisors
  show IsDedekindDomain.HeightOneSpectrum.valuation (FractionRing R) _ _ ≤ 1
  rw [show (algebraMap R (FractionRing R) r)
      = algebraMap (Localization.AtPrime p) (FractionRing R)
        (algebraMap R (Localization.AtPrime p) r) from
    (IsScalarTower.algebraMap_apply R (Localization.AtPrime p) (FractionRing R) r)]
  exact IsDedekindDomain.HeightOneSpectrum.valuation_le_one _ _

/-- ★非零元の位数は非零。 -/
theorem ordAt_ne_zero {R : Type*} [CommRing R] [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] (p : Ideal R) [hp : p.IsPrime] (hp0 : p ≠ ⊥)
    (hmin : ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p) {x : FractionRing R}
    (hx : x ≠ 0) : ordAt p hp0 hmin x ≠ 0 :=
  (ordAt p hp0 hmin).ne_zero_iff.mpr hx

end ABC3.Found.Divisor
