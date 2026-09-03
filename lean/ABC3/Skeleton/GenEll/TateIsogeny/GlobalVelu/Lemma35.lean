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
import ABC3.Skeleton.GenEll.TateIsogeny.LocalJ
import ABC3.Skeleton.GenEll.TateIsogeny.GlobalVelu.Lemma32

/-!
# GlobalVelu —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★★★★★第 999 —— 大域の Vélu の商で受ける形（軽い版）

★第 972 は `E′` について
分裂乗法還元 `h′`・極小モデル `[IsMinimal (E′ ⊗ Lv)]`・`C′`・`hC′`・`hc4ne′`・`hc4′`、
さらに `E` について `C`・`hC`・`hc4ne`・`hc4` を要求していた。

☆第 996（`j` の段で止める）＋第 997（`E′` の仮説は 2 本）＋第 998（`E′.j ≠ 0` は自動）で、
それらは**すべて半安定性 2 本に置き換わる**。★残るのは `E` の側の局所データだけである。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 各悪い素点での `Δ_min` の関係——
大域の Vélu の商で受ける形（`E′` の局所データを要らない版）**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 999）**——第 972 の軽量版。
☆`E′` について要るのは `SemistableAt p E′` だけになった。 -/
theorem minDeltaExp_eq_mul_of_globalVelu' {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
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
  -- ★`μ_l` の形での楕円性（第 971）
  haveI hellMu := isElliptic_veluQuotient_tate_mu
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hlu h2K
      ζ uζ hζ hζu hζl hord v w hv hw hell
  -- ☆`P` の形での楕円性
  haveI hellP : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic := by
    rw [hPz]; exact hellMu
  -- ★曲線の水準で `P` と `μ_l` を繋ぐ（`j` を跨がないので motive が壊れない）
  have hcurveEq : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := p.adicCompletion L)
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ
          (QuotientGroup.mk uζ)))) := by
    rw [hPz]
    rfl
  -- ☆第 996 で `E_{q^l}` に繋ぐ
  have hjtate := j_veluQuot_eq_j_tate_pow
    (tateParamR (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
    (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hζ hlu
    uζ hζu hζl hord hql h2 h2K hodd v w hv hw
  exact minDeltaExp_eq_mul_of_j_tate_pow p hp E E' h hssE hssE' hjneg hql
    ((hj hellP).trans ((ABC3.Found.GenEll.j_congr_curve hcurveEq).trans hjtate))

def minDeltaExp_eq_mul_of_globalVelu'.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各悪い素点での Δ_min の関係——E′ の局所データを要らない版)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_globalVelu'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_j_tateModel(P と hW′j、第 970、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_point_j_tateModel") 1,
    .citation "[ABC3]" "j_veluQuot_eq_j_tate_pow(j の段で止める、第 996、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.j_veluQuot_eq_j_tate_pow") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_j_tate_pow(終点、第 997、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_j_tate_pow") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1000 —— 悪い素点で局所データを供給する

★第 999 が受ける局所データのうち、**自前で作れるものはすべて作る**:

| 入力 | 出どころ |
|---|---|
| `hp`（付値の橋） | 第 964 `valuation_algebraMap_adicCompletion` |
| `hΔ` | 第 977 `tateModel_map_Delta_ne_zero` |
| `hcop`（付値の言葉） | 第 978 `not_dvd_vAdd_tateParam_of_not_dvd_jExp` |
| `hql` | 第 977 `pow_mem_of_mem_ideal` |
| `h2`・`h2K` | 第 977 |
| `E_{q^l} ⊗ Lv` の楕円性 | `tateCurveAt_map_isElliptic` ＋ `tateCurveAt_c4_isUnit` |

☆残るのは **3 本だけ**である:

* `hmin`（完備化で極小）——第 973
* `h`（完備化で分裂乗法還元）——第 976＋993（`p ∣ 2` の非分裂が残件）
* `hlu`（`l` が単元）・`hvw`（Vélu の係数が整）
-/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
局所データを自前で作る形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1000）**——第 999 の 10 本の局所入力を 3 本に絞る。 -/
theorem minDeltaExp_eq_mul_at_bad_prime {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      (E.baseChange (p.adicCompletion L)))
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hlu : IsUnit ((l : (p.adicCompletionIntegers L))))
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
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  have hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性——`c₄` は定数項 `1`、`1/j` は `q^l · 単元`
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  exact minDeltaExp_eq_mul_of_globalVelu' p E E' h hp hssE hssE' hjneg hl hodd
    (tateModel_map_Delta_ne_zero (E.baseChange (p.adicCompletion L)) h)
    (not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hjneg hcop)
    hlu hql (two_ne_zero_adicCompletionIntegers L p) (two_ne_zero_adicCompletion L p)
    hvw hQ hE'

def minDeltaExp_eq_mul_at_bad_prime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——局所データを自前で作る形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_globalVelu′(第 999、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_globalVelu'") 1,
    .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp(hcop の言い換え、第 978、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp") 1 ]


def minDeltaExp_eq_mul_of_globalVelu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各悪い素点での Δ_min の関係——大域の Vélu の商で受ける形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_globalVelu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_j_tateModel(P と hW′j、第 970、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_point_j_tateModel") 1,
    .citation "[ABC3]" "isElliptic_veluQuotient_tate_mu(hellQ、第 971、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluQuotient_tate_mu") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_torsion(終点、第 965、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_torsion") 1 ]

end ABC3.Skeleton.GenEll
