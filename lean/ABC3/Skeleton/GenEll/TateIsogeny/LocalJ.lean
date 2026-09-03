/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.Lemma35Concrete
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.TateDXNeZero
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateSetupDvr
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.DegInfTateParam
import ABC3.Found.GenEll.MuPrimitiveRoot
import ABC3.Found.GenEll.CyclotomicUnits
import ABC3.Found.GaloisRep.TateModelPoint
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.CompletionValuationBridge
import ABC3.Found.GaloisRep.UnramQuad
import ABC3.Found.GaloisRep.TateMuInvolution
import ABC3.Found.GaloisRep.VeluTateDelta
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Found.GaloisRep.VeluMuSum
import ABC3.Found.GenEll.JScale
import ABC3.Meta.Claim
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Skeleton.GenEll.TateIsogeny.Target

/-!
# TateIsogeny —— `ζ` を消す・`j` で受ける捉れ点版・局所の終点・大域の Vélu の商で受ける形

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★★★★★第 947 と 927 を繋ぐ——`ζ` を消す -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 有理な `l`-捉れ点だけで
`q_{E′} = q_E^l`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 948）**——`tateParam_quot_velu_dvr`（第 927）は
`ζ`・`uζ`・`hζ`・`hζu`・`hζl`・`hord` の 6 つを受けていた。
☆本定理はそれを `exists_primitiveRoot_of_torsion_point`（第 947）で埋め、

    **`P` が位数 `l` の点で、`l ∤ v(q)`**

だけに置き換える。★`ζ` はもはや引数に現れない——
Vélu の帳簿（`hu`・`hv`・`hw`）だけが `ζ` について全称で残る。

☆これが `isMuAtBadPrimes_of_veluQuotient` の局所の段そのものである。 -/
theorem tateParam_quot_velu_of_torsion {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0)
    (hvw : ∀ ζ : R, IsPrimitiveRoot ζ l → ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (P : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    (hW' : W' = veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))) :
    tateParamR W' hsplit = q ^ l := by
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point q hq hq0 hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  haveI := hell
  refine tateParam_quot_velu_dvr q hq hq0 hΔ hl hζ hlu
    (ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu)
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw W' hsplit ?_
  rw [hW', hPz]
  rfl

def tateParam_quot_velu_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有理な l-捉れ点だけで q_{E′} = q_E^l)",
    sectionId := "genell-lemma-3-5" }

def tateParam_quot_velu_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_primitiveRoot_of_torsion_point(ζ を作る、第 947、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point") 1,
    .citation "[ABC3]" "tateParam_quot_velu_dvr(Vélu の商の Tate 母数、第 927、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_dvr") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★`j` で受ける形の捉れ点版 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 有理な `l`-捉れ点と `j` の一致だけで
`q_{E′} = q_E^l`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 955）**——第 948 は `W′` が Vélu の商に**等しい**ことを
要求していたが、実際に得られるのは変数変換を除いた `j` の一致である。
☆`j_veluQuotientFull_nsmul_variableChange`（第 950）がその `j` を与えるので、
本定理が**そのままの形で受けられる**。

★`ζ` は引数に現れない——`exists_primitiveRoot_of_torsion_point`（第 947）が作る。 -/
theorem tateParam_quot_velu_j_of_torsion {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) q hq hq0).v
      (mkTateSetup (K := K) q hq hq0).Q))
    (hlu : IsUnit ((l : R)))
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0)
    (hvw : ∀ ζ : R, IsPrimitiveRoot ζ l → ∃ v w : R,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt q hq)
            (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt q hq)
              (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        ∧ ((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    (P : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (hellQ : (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic)
    (W' : WeierstrassCurve K) [W'.IsElliptic] [W'.IsMinimal R]
    (hsplit : W'.HasSplitMultiplicativeReduction R)
    (hW'j : W'.j = (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j) :
    tateParamR W' hsplit = q ^ l := by
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point q hq hq0 hΔ hl hcop P hP hP0
  subst hPz
  haveI hQell : (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic := hellQ
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  haveI := hell
  exact tateParam_quot_velu_j_dvr q hq hq0 hΔ hl hζ hlu
    (ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu)
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw W' hsplit hW'j

def tateParam_quot_velu_j_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有理な l-捉れ点と j の一致だけで q_{E′} = q_E^l)",
    sectionId := "genell-lemma-3-5" }

def tateParam_quot_velu_j_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_primitiveRoot_of_torsion_point(ζ を作る、第 947、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point") 1,
    .citation "[ABC3]" "tateParam_quot_velu_j_dvr(j で受ける形、第 914、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_j_dvr") 1,
    .citation "[ABC3]" "isUnit_one_sub_pow_of_isUnit_natCast(hu は hlu から、第 951、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★捉れ点で受ける形の局所の終点 -/

open ABC3.Skeleton.GenEll ABC3.Found.GaloisRep ABC3.Found.GenEll IsDedekindDomain
  NumberField WeierstrassCurve Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 数体の素点での `Δ_min` の関係——
捉れ点と `j` で受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 965）**——第 904 は `ζ`・`uζ`・`hu`・`v`・`w` を受け、
さらに `W′` が Vélu の商に**等しい**ことを要求していた。
☆本定理はそれを

* 位数 `l` の点 `P`（`P ≠ 0`）と `l ∤ v(q)`（第 947・946）
* `j` の一致（第 950 が与える）
* Vélu の帳簿 `hvw`（第 961・962 が与える）

に置き換える。★これが `isMuAtBadPrimes_of_veluQuotient` の
**各悪い素点での終点**である。 -/
theorem minDeltaExp_eq_mul_of_torsion {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    (P : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Point)
    (hP : l • P = 0) (hP0 : P ≠ 0)
    (hellQ : (veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
          (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic)
    (hW'j : (E'.baseChange (p.adicCompletion L)).j
      = (veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
          (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).j) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hqpow := tateParam_quot_velu_j_of_torsion
    (tateParamR (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)
    hΔ hl hodd hcop hlu hql h2 h2K hvw P hP hP0 hellQ
    (E'.baseChange (p.adicCompletion L)) h' hW'j
  exact minDeltaExp_eq_mul_of_tateParamR (R := p.adicCompletionIntegers L) E E' l h h' p hp
    C hC hc4ne hc4 C' hC' hc4ne' hc4' hqpow

def minDeltaExp_eq_mul_of_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(数体の素点での Δ_min の関係——捉れ点と j で受ける形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_torsion.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateParam_quot_velu_j_of_torsion(q_{E′} = q_E^l、第 955、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.tateParam_quot_velu_j_of_torsion") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_tateParamR(Δ_min へ、第 892、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_tateParamR") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★大域の Vélu の商で受ける形 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 各悪い素点での `Δ_min` の関係——
大域の Vélu の商で受ける形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 972）**——第 965 は Tate モデルの上の点 `P` と
`j` の一致 `hW′j`、それに商の楕円性 `hellQ` を受けていた。
☆本定理はその 4 つを**大域のデータ `Q`・`hQ`・`hE′` から作る**:

* `P`・`hP`・`hP0`・`hW′j` ← `exists_point_j_tateModel`（第 970）
* `ζ`・`uζ` ← `exists_primitiveRoot_of_torsion_point`（第 947）
* `hellQ` ← `isElliptic_veluQuotient_tate_mu`（第 971）

★残るのは**各悪い素点で局所データを供給する段**だけである。 -/
theorem minDeltaExp_eq_mul_of_globalVelu {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (h' : (E'.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletionIntegers L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateXpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYpair (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateXpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYpair (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXpair (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ ((veluCurve (tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)) v w).map
            (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel p E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  have hellQ := isElliptic_veluQuotient_tate_mu
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hlu h2K
      ζ uζ hζ hζu hζl hord v w hv hw hell
  rw [← hPz] at hellQ
  exact minDeltaExp_eq_mul_of_torsion p E E' h h' hp C hC hc4ne hc4 C' hC' hc4ne' hc4'
    hl hodd hΔ hcop hlu hql h2 h2K hvw P hP hP0 hellQ (hj hellQ)



end ABC3.Skeleton.GenEll
