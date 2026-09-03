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
import ABC3.Skeleton.GenEll.TateIsogeny.Extension

/-!
# TateIsogeny —— 局所の `Δ` の計算・Vélu の商の整性（悪い素点）

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★★★★★第 1061 —— 局所の `Δ` の計算（組み立て）

★第 1055-1060 の部品を第 1027 と同じ形で組む。
☆結果は **`v(Δ(E′ ⊗ Lv)) = 12·v(l) + l·v(q)`** である。 -/

open Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★**Vélu の商の `Δ` の付値**（第 1061）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`E` が完備化で極小・分裂乗法還元をもつとき、
`v(Δ(E′ ⊗ Lv)) = 12·v(l) + l·v(q_E)` である。
★これが第 1051 の見通し「`neronExp p E′ = v_p(l)`」の中身である。 -/
theorem vAdd_Delta_veluQuotient_tate {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q))
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hΔ' : (E'.baseChange Lv).Δ ≠ 0)
    (hlne : algebraMap R Lv (l : R) ≠ 0)
    (hqne : algebraMap R Lv (tateParamR (E.baseChange Lv) h) ≠ 0) :
    vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).Δ) hΔ')
      = 12 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap R Lv (l : R)) hlne)
        + l * vAdd (tateDvrVal R Lv)
          (Units.mk0 (algebraMap R Lv (tateParamR (E.baseChange Lv) h)) hqne) := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    ABC3.Found.GaloisRep.pow_mem_of_mem_ideal hq hl.pos
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop P hP hP0
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((ABC3.Found.GaloisRep.tateCurveAt_c4_isUnit
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu'u, hueq'⟩ := ABC3.Found.GaloisRep.evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    ABC3.Found.GaloisRep.tateCurveAt_map_isElliptic _ hql hev' hc4T'
  obtain ⟨v, w, hv, hw, hell, h4, h6⟩ :=
    exists_vw_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hodd hlu hql h2 ζ hζ
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  have hquot := ABC3.Found.GaloisRep.veluQuotientFull_tate_mu
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  have hquot' : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ))))
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := hquot
  have hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := by
    rw [← hcurve, hPz]
    exact hquot'
  have hu0 := ABC3.Found.GaloisRep.vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  have hTne : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).Δ ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2
      (ABC3.Found.GaloisRep.tateCurveAt_Delta_ne_zero hql (pow_ne_zero l hq0))
  exact ABC3.Found.GaloisRep.vAdd_Delta_of_veluCurve_eq
    (tateParamR (E.baseChange Lv) h) hq hql v w h4 h6 C₀ (E'.baseChange Lv) hΔ' hEq hu0
    hTne hlne hqne

def vAdd_Delta_veluQuotient_tate.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の Δ の付値は 12·v(l) + l·v(q_E))",
    sectionId := "genell-lemma-3-5" }

def vAdd_Delta_veluQuotient_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "vAdd_Delta_of_veluCurve_eq(帳簿、第 1059、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_Delta_of_veluCurve_eq") 1,
    .citation "[ABC3]" "exists_variableChange_veluQuotient_tateModel(輸送、第 1058、証明済み)"
      (.inProject "ABC3"
        "ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel") 1,
    .citation "[ABC3]" "vAdd_tateModel_u_eq_zero(極小性、第 1056、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_tateModel_u_eq_zero") 1 ]

open Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**極小な素点では Vélu の商の Néron 指数は
`v_p(l)`**（第 1062）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1062）**——第 1051 の見通しの中身である。
☆`minDeltaExp = v_p(Δ) − 12·neronExp` に
第 1061（`v(Δ(E′)) = 12·v(l) + l·v(q)`）と
第 1044（`minDeltaExp p E′ = l·minDeltaExp p E`）を代入する。 -/
theorem neronExp_veluQuotient_eq_of_minimal {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q))
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hminE : neronExp p E = 0)
    (hmu : minDeltaExp p E' = l * minDeltaExp p E)
    (hlL : ((l : ℕ) : L) ≠ 0)
    (hlR : algebraMap R Lv (l : R) ≠ 0)
    (hlLv : algebraMap L Lv ((l : ℕ) : L) ≠ 0)
    (hcast : algebraMap R Lv (l : R) = algebraMap L Lv ((l : ℕ) : L))
    (hqE : vAdd (tateDvrVal R Lv)
        (Units.mk0 (algebraMap R Lv (tateParamR (E.baseChange Lv) h))
          ((map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2
            (tateParamR_ne_zero (E.baseChange Lv) h)))
      = minDeltaExp p E) :
    neronExp p E' = valAdd p (Units.mk0 ((l : ℕ) : L) hlL) := by
  have hΔE : E.Δ ≠ 0 := (inferInstance : E.IsElliptic).isUnit.ne_zero
  have hΔE' : E'.Δ ≠ 0 := (inferInstance : E'.IsElliptic).isUnit.ne_zero
  have hΔ'Lv : (E'.baseChange Lv).Δ ≠ 0 := by
    show (E'.map (algebraMap L Lv)).Δ ≠ 0
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (algebraMap L Lv).injective).2 hΔE'
  have hkey := vAdd_Delta_veluQuotient_tate E E' h hl hodd hcop hlu h2 h2K hQ hE'
    hΔ'Lv hlR ((map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2
      (tateParamR_ne_zero (E.baseChange Lv) h))
  -- ★左辺を `valAdd p (Δ E′)` に直す
  have hlhs : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).Δ) hΔ'Lv)
      = valAdd p (Units.mk0 E'.Δ hΔE') := by
    have he : (Units.mk0 ((E'.baseChange Lv).Δ) hΔ'Lv)
        = Units.mk0 (algebraMap L Lv E'.Δ)
          ((map_ne_zero_iff _ (algebraMap L Lv).injective).2 hΔE') := by
      refine Units.ext ?_
      exact WeierstrassCurve.map_Δ _ _
    rw [he]
    exact ABC3.Found.GaloisRep.vAdd_algebraMap_eq_valAdd p hp E'.Δ hΔE' _
  -- ☆`v(l)` を `valAdd p (l)` に直す
  have hlv : vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap R Lv (l : R)) hlR)
      = valAdd p (Units.mk0 ((l : ℕ) : L) hlL) := by
    have he : (Units.mk0 (algebraMap R Lv (l : R)) hlR)
        = Units.mk0 (algebraMap L Lv ((l : ℕ) : L)) hlLv := by
      refine Units.ext ?_
      exact hcast
    rw [he]
    exact ABC3.Found.GaloisRep.vAdd_algebraMap_eq_valAdd p hp ((l : ℕ) : L) hlL _
  rw [hlhs, hlv, hqE] at hkey
  -- ★`minDeltaExp` の定義に流す
  have hdE : minDeltaExp p E = valAdd p (Units.mk0 E.Δ hΔE) := by
    rw [minDeltaExp, dif_neg hΔE, hminE]
    ring
  have hdE' : minDeltaExp p E' = valAdd p (Units.mk0 E'.Δ hΔE') - 12 * neronExp p E' := by
    rw [minDeltaExp, dif_neg hΔE']
  omega

def neronExp_veluQuotient_eq_of_minimal.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小な素点では Vélu の商の Néron 指数は v_p(l))",
    sectionId := "genell-lemma-3-5" }

def neronExp_veluQuotient_eq_of_minimal.needs : List ProofObligation :=
  [ .citation "[ABC3]" "vAdd_Delta_veluQuotient_tate(第 1061、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.vAdd_Delta_veluQuotient_tate") 1,
    .citation "[ABC3]" "vAdd_algebraMap_eq_valAdd(付値の橋、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_valAdd") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1066 —— Vélu の商の整性（悪い素点）

★第 1003 が Vélu の係数 `v`・`w` を **`R` の元として**与えるので、
`veluCurve (tate q) v w` は `R`-曲線である。
☆第 1058 の輸送で `C₀ • (E′ ⊗ Lv)` がその底変換に等しく、
第 1065 で整性が戻る。 -/

open Finset in
open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**悪い素点では Vélu の商は整**（第 1066）。 -/
theorem isIntegral_veluQuotient_baseChange_of_split {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q))
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    : WeierstrassCurve.IsIntegral R (E'.baseChange Lv) := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    ABC3.Found.GaloisRep.pow_mem_of_mem_ideal hq hl.pos
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop P hP hP0
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((ABC3.Found.GaloisRep.tateCurveAt_c4_isUnit
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu'u, hueq'⟩ := ABC3.Found.GaloisRep.evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    ABC3.Found.GaloisRep.tateCurveAt_map_isElliptic _ hql hev' hc4T'
  obtain ⟨v, w, hv, hw, hell, h4, h6⟩ :=
    exists_vw_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hodd hlu hql h2 ζ hζ
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  have hquot := ABC3.Found.GaloisRep.veluQuotientFull_tate_mu
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  have hquot' : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ))))
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := hquot
  have hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := by
    rw [← hcurve, hPz]
    exact hquot'
  have hu0 := ABC3.Found.GaloisRep.vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  have hTne : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).Δ ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2
      (ABC3.Found.GaloisRep.tateCurveAt_Delta_ne_zero hql (pow_ne_zero l hq0))
  exact ABC3.Found.GaloisRep.isIntegral_of_variableChange_map (E'.baseChange Lv)
    (ABC3.Found.GenEll.veluCurve (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w)
    C₀ hEq

def isIntegral_veluQuotient_baseChange_of_split.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は整)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_veluQuotient_baseChange_of_split.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isIntegral_of_variableChange_map(第 1065、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isIntegral_of_variableChange_map") 1,
    .citation "[ABC3]" "exists_vw_tate_mu(v・w が R の元、第 1003、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_vw_tate_mu") 1 ]

end ABC3.Skeleton.GenEll
