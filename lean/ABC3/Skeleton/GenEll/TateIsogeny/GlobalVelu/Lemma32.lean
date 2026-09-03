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

/-!
# GlobalVelu —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★★★★★第 996 —— `j` の段で止める

★第 995 により `E′` の側に要るのは **`j` の一致だけ**になった。
☆本ブロックは第 914（`tateParam_quot_velu_j_dvr`）の**前半**を取り出す:

    `j(μ_l による Vélu の商) = j(E_{q^l})`

★第 914 はこの後 `tateParamR_eq_of_j_tateCurveAt` に流して `q_{E′} = q_E^l` を出すが、
そこで `E′` の分裂性が要る。☆`j` の段で止めれば要らない。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] `μ_l` による Vélu の商の `j` は
`E_{q^l}` の `j`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 996）**——第 914 の前半である。
☆`veluQuotientFull_tate_mu`（第 890）で商を `veluCurve` の形にし、
`j_velu_tate_mu_map`（第 886）で `E_{q^l}` に繋ぐ。
★`E′` の分裂性も極小モデルも要らない。 -/
theorem j_veluQuot_eq_j_tate_pow {R : Type} [CommRing R] [IsDomain R] [CharZero R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R)))
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : q ^ l ∈ IsLocalRing.maximalIdeal R)
    (h2 : (2 : R) ≠ 0) (h2K : (2 : K) ≠ 0) (hodd : l ≠ 2)
    (v w : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq)
        (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq)
          (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq)
                (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
    [((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic]
    [((tateCurveAt (q ^ l) hql).map (algebraMap R K)).IsElliptic]
    [(veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
              (QuotientGroup.mk uζ))))).IsElliptic] :
    (veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))).j
      = ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).j := by
  haveI i1 : ((veluCurve (tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq) v w).map (algebraMap R K)).IsElliptic :=
    inferInstanceAs (((veluCurve (tateCurveAt q hq) v w).map (algebraMap R K)).IsElliptic)
  haveI i2 : (veluQuotientFull ((tateCurveAt (mkTateSetup (K := K) q hq hq0).q
      (mkTateSetup (K := K) q hq hq0).hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic :=
    inferInstanceAs ((veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic)
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  have hquot := veluQuotientFull_tate_mu (mkTateSetup q hq hq0) hΔ
    (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  -- ☆`(mkTateSetup q hq hq0).q` と `q` は defeq だが構文的には違うので、
  -- `have` でゴールの形に言い直してから `rw` する（第 971 と同じ穴）。
  have hquot' : veluQuotientFull ((tateCurveAt q hq).map (algebraMap R K))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := K) q hq hq0) hΔ (QuotientGroup.mk uζ))))
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R K) := hquot
  rw [ABC3.Found.GenEll.j_congr_curve hquot']
  exact j_velu_tate_mu_map hl hζ hlu hu q hq hql h2
    (fun i hi => tateDXpair_ne_zero_of_mu (mkTateSetup q hq hq0) hΔ
      (dvrTatePhiAddEquiv q hq hq0 hΔ) (fun _ => rfl) hl hodd ζ uζ hζu hζl hord i
      (Finset.mem_erase.1 hi).1 (Finset.mem_range.1 (Finset.mem_erase.1 hi).2) (hu i hi))
    v w hv hw

def j_veluQuot_eq_j_tate_pow.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l による Vélu の商の j は E_{q^l} の j)",
    sectionId := "genell-lemma-3-2" }

def j_veluQuot_eq_j_tate_pow.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu(商を veluCurve の形に、第 890、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu") 1,
    .citation "[ABC3]" "j_velu_tate_mu_map(veluCurve の j は E_{q^l} の j、第 886、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.j_velu_tate_mu_map") 1 ]

end ABC3.Skeleton.GenEll
