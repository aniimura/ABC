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
import ABC3.Skeleton.GenEll.TateIsogeny.LocalData

/-!
# TateIsogeny —— 一般の局所体に開く・拡大に載せ替える・非分裂の枝を閉じる

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep WeierstrassCurve IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-! ## ★★★★★★★★★★★★★★★★第 1027 —— 連鎖を一般の局所体に開く

★第 1004 は `p.adicCompletion L` に固定されていた。
☆第 1026 で `exists_point_j_tateModel` を一般の `(Lv, R)` に開いたので、
連鎖全体も開ける——他の部品（996・997・971・1003）はもともと一般である。

★★これで**不分岐 2 次拡大 `Lv′ = Frac (R[X]/(f))`（第 1012-1025）を
そのまま通せる**ようになった。 -/

open Finset in
theorem minDeltaExp_eq_mul_at_bad_prime_gen {L : Type} [Field L] [NumberField L]
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
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    pow_mem_of_mem_ideal hq hl.pos
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'.map (algebraMap R Lv)).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  have hcop' := not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hjneg hcop
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel_gen E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop'
      P hP hP0
  obtain ⟨v, w, hv, hw, hell, _, _⟩ :=
    exists_vw_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hodd hlu hql h2 ζ hζ
  haveI hellMu := isElliptic_veluQuotient_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ
    hl hlu h2K ζ uζ hζ hζu hζl hord v w hv hw hell
  haveI hellP : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic := by
    rw [hPz]; exact hellMu
  have hcurveEq : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ)))) := by
    rw [hPz]
    rfl
  have hjtate := j_veluQuot_eq_j_tate_pow (tateParamR (E.baseChange Lv) h) hq hq0 hΔ
    hl hζ hlu uζ hζu hζl hord hql h2 h2K hodd v w hv hw
  exact minDeltaExp_eq_mul_of_j_tate_pow p hp E E' h hssE hssE' hjneg hql
    ((hj hellP).trans ((ABC3.Found.GenEll.j_congr_curve hcurveEq).trans hjtate))

def minDeltaExp_eq_mul_at_bad_prime_gen.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——一般の局所体で)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_gen.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_point_j_tateModel_gen(第 1026、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_point_j_tateModel_gen") 1,
    .citation "[ABC3]" "exists_vw_tate_mu(hvw の中身、第 1003、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_vw_tate_mu") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_j_tate_pow(終点、第 997、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq_mul_of_j_tate_pow") 1 ]

/-! ## ★★★★★★★★★★★★★★★★第 1028 —— 拡大に載せ替える

★第 1027 は一般の `(Lv, R)` で書けている。
☆本ブロックは **`(Lv′, R′) = (Frac (R[X]/(f)), R[X]/(f))`** を実際に流し込み、
付値・完備性・標数・`l` の単元性を**すべて自前で作る**。

☆残るのは `E ⊗ Lv′` の側の 3 本だけである:
極小性・楕円性・**分裂乗法還元**。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
不分岐 2 次拡大に載せ替えた形**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1028）**——第 1027 に第 1012-1025 を流し込む段である。
☆`hp′` は第 1024（付値の橋）、完備性は第 1012＋1013＋1022 で作る。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_ext {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    {f : Polynomial R} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (FractionRing (AdjoinRoot f))).IsElliptic]
    [(E.baseChange (FractionRing (AdjoinRoot f))).IsMinimal (AdjoinRoot f)]
    [(E'.baseChange (FractionRing (AdjoinRoot f))).IsElliptic]
    (h' : (E.baseChange (FractionRing (AdjoinRoot f))).HasSplitMultiplicativeReduction
      (AdjoinRoot f))
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : IsUnit ((l : R)))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hloc : IsLocalRing (AdjoinRoot f) := inferInstance
  -- ★完備性——第 1012（`R`-加群）→ 第 1013（イデアルの側）→ 第 1022（極大の同定）
  haveI hcomp0 : IsAdicComplete (IsLocalRing.maximalIdeal R) (AdjoinRoot f) :=
    isAdicComplete_adjoinRoot (IsLocalRing.maximalIdeal R) hf hdeg
  haveI hcomp : IsAdicComplete (IsLocalRing.maximalIdeal (AdjoinRoot f)) (AdjoinRoot f) := by
    rw [maximalIdeal_adjoinRoot_eq_map hf hdeg hirr]
    exact isAdicComplete_map_algebraMap _
  -- ☆標数
  haveI hchar : CharZero (AdjoinRoot f) :=
    charZero_of_injective_algebraMap (algebraMap_adjoinRoot_injective hf hdeg)
  haveI hcharK : CharZero (FractionRing (AdjoinRoot f)) :=
    IsFractionRing.charZero (AdjoinRoot f)
  -- ★付値の橋
  have hp' : ∀ x : L, (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
      (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f)))
      (algebraMap L (FractionRing (AdjoinRoot f)) x)
      = (HeightOneSpectrum.valuation L p) x := by
    intro x
    rw [halg, valuation_quadFieldHom hf hdeg hirr, hp]
  -- ☆`l` は `R′` でも単元
  have hlu' : IsUnit ((l : AdjoinRoot f)) := by
    have := hlu.map (algebraMap R (AdjoinRoot f))
    rwa [map_natCast] at this
  exact minDeltaExp_eq_mul_at_bad_prime_gen p hp' E E' h' hssE hssE' hjneg hl hodd hcop
    hlu' (NeZero.ne' (2 : AdjoinRoot f)).symm (NeZero.ne' (2 : FractionRing (AdjoinRoot f))).symm
    hQ hE'

def minDeltaExp_eq_mul_at_bad_prime_ext.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——不分岐 2 次拡大に載せ替えた形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_ext.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_gen(第 1027、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_gen") 1,
    .citation "[ABC3]" "valuation_quadFieldHom(付値の橋、第 1024、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.valuation_quadFieldHom") 1,
    .citation "[ABC3]" "isAdicComplete_adjoinRoot(完備性、第 1012、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isAdicComplete_adjoinRoot") 1 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★第 1034 —— 非分裂の枝を閉じる

★第 1033 までで不分岐 2 次拡大の部品はすべて建った。
☆本ブロックは**呼び出し側**である——`f` を作り、instance を並べ、
第 1029（極小・乗法還元）と第 1033（分裂）を第 1028 に流す。

★`E` 側の極小化（`C`）と Vélu の商の輸送（第 969）は第 1002 と同じ形である。 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] 非分裂でも `Δ_min` は `l` 倍**——
不分岐 2 次拡大を通す形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1034）**——`hirr`（拡大が本当に 2 次であること、
すなわち曲線が非分裂であること）だけを受ければ閉じる。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_nonsplit {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    [((C • E).baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv)]
    (hA : IsUnit (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)).c₄)
    (hirr : Irreducible
      ((splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).map
        (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hlu : IsUnit ((l : R)))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  have hfm : (splitQuadPoly
      (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).Monic :=
    monic_splitQuadPoly _ hA
  have hfd : (splitQuadPoly
      (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).natDegree = 2 :=
    natDegree_splitQuadPoly _ hA
  haveI hdom := isDomain_adjoinRoot hfm hirr
  haveI hlocR := isLocalRing_adjoinRoot hfm hfd hirr
  haveI hdvr := isDiscreteValuationRing_adjoinRoot hfm hfd hirr
  letI : Algebra L (FractionRing (AdjoinRoot
      (splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA))) :=
    ((quadFieldHom (K := Lv) hfm hfd).comp (algebraMap L Lv)).toAlgebra
  have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot
      (splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA))) x
      = quadFieldHom hfm hfd (algebraMap L Lv x) := fun _ => rfl
  haveI hmin' := isMinimal_baseChange_ext p hp hfm hfd hirr halg E C hC hc4ne hc4
  haveI hmult' := hasMultiplicativeReduction_ext p hp hfm hfd hirr halg E C hC hc4ne hc4 hjneg
  have h' := hasSplitMultiplicativeReduction_ext (C • E) hA rfl halg
  -- ★点と Vélu の商を `C • E` に運ぶ（第 1002 と同じ）
  have h2L : (2 : L) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq := ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq C E E' hQ h2L hE'
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hjneg
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [jExp_variableChange p E C]; exact hcop
  have hkey := minDeltaExp_eq_mul_at_bad_prime_ext p hp hfm hfd hirr halg (C • E) (C • E') h'
    (semistableAt_variableChange p E C hssE) (semistableAt_variableChange p E' C hssE')
    hjC hl hodd hcopC hlu hQ' hEq
  rwa [minDeltaExp_variableChange p E' C, minDeltaExp_variableChange p E C] at hkey

def minDeltaExp_eq_mul_at_bad_prime_nonsplit.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非分裂でも Δ_min は l 倍——不分岐 2 次拡大を通す形)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_nonsplit.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_ext(第 1028、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_ext") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_ext(分裂性、第 1033、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasSplitMultiplicativeReduction_ext") 1,
    .citation "[ABC3]" "isMinimal_baseChange_ext・hasMultiplicativeReduction_ext(第 1029)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_baseChange_ext") 1 ]

end ABC3.Skeleton.GenEll
