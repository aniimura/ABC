/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateLocalModel
import ABC3.Skeleton.GenEll.TateIsogenyK
import ABC3.Meta.Claim

/-!
# 第 1140 ブロック —— **`IsMuAtBadPrimes` の `p ∤ l` なし版**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

第 1137-1139 で `minDeltaExp` の連鎖が `hlu`（＝`p ∤ l`）なしで組み直された。
★本ブロックはそれを `IsMuAtBadPrimes` まで上げる。

☆これで `isMuAtBadPrimes_of_veluQuotient_of_large`（第 1044）が使っていた
`hd : [L:ℚ] + 1 < l`（★原典に無い仮説）が不要になる。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Found.GenEll WeierstrassCurve IsDedekindDomain
  NumberField Finset
open scoped Classical

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**悪い素点での `Δ_min` の関係（分裂性を仮定しない形、
第 1140、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_any_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletion L p
  haveI := ABC3.Found.GaloisRep.charZero_adicCompletionIntegers L p
  obtain ⟨C, hC, hc4ne, hc4⟩ :=
    ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg p E hssE hjneg
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hmin :=
    ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime p hp E C hC hc4ne hc4
  haveI hmult :=
    ABC3.Found.GaloisRep.hasMultiplicativeReduction_at_bad_prime p hp E C hC hc4ne hc4 hjneg
  have hA := ABC3.Found.GaloisRep.integralModel_c4_isUnit
    (R := p.adicCompletionIntegers L) ((C • E).baseChange (p.adicCompletion L))
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L))
  · exact minDeltaExp_eq_mul_at_bad_prime_full_K p E E' hssE hssE' hjneg hl hcop
      C hmin hs hQ hE'
  · -- ★非分裂——不分岐 2 次拡大を通す（第 1034）
    have hirr := ABC3.Found.GaloisRep.irreducible_map_residue_of_not_splits
      (ABC3.Found.GaloisRep.monic_splitQuadPoly _ hA)
      (ABC3.Found.GaloisRep.natDegree_splitQuadPoly _ hA)
      (fun hsp => hs (ABC3.Found.GaloisRep.hasSplit_of_splits_splitQuadPoly _ hA hsp))
    exact minDeltaExp_eq_mul_at_bad_prime_nonsplit_K p hp E E' C hC hc4ne hc4 hA hirr
      hssE hssE' hjneg hl hcop hQ hE'

def minDeltaExp_eq_mul_at_bad_prime_any_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——分裂性を仮定しない形。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_any_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_full_K(分裂の枝、第 1138、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_full_K") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_nonsplit_K(非分裂の枝、第 1139、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_nonsplit_K") 1 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] `IsMuAtBadPrimes`——仮説は原文自身のものだけ**
（第 1140）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`p ∤ l` も `[L:ℚ] + 1 < l` も取らない**。 -/
theorem isMuAtBadPrimes_of_veluQuotient_of_coprime_K {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) :
    IsMuAtBadPrimes E E' l := fun p hbad =>
  minDeltaExp_eq_mul_at_bad_prime_any_K p E E' (hssE p) (hssE' p) hbad hl
    (hcop p hbad) hQ hE'

def isMuAtBadPrimes_of_veluQuotient_of_coprime_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsMuAtBadPrimes——仮説は原文自身のものだけ。★p ∤ l を取らない)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient_of_coprime_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_any_K(第 1140、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_any_K") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1140）**——`isMuAtBadPrimes_of_veluQuotient_of_large`（第 1044）が" ++
       "使っていた `hd : [L:ℚ] + 1 < l`（★原典に無い仮説）は、" ++
       "`p ∤ l` を出すためだけに必要だった。" ++
       "☆本ブロックで `p ∤ l` が不要になったので、`hd` も不要になる。") 6 ]

open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] `IsMuAtBadPrimes`——原文に無い仮説は無い**（第 1141）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`hd : [L:ℚ] + 1 < l` を取らない**——第 1140 で `p ∤ l` が不要になったので、
それを出すためだけにあった `hd` も落ちた。 -/
theorem isMuAtBadPrimes_of_veluQuotient_nodeg {L : Type} [Field L] [NumberField L]
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hssE : ∀ p, SemistableAt p E) (hssE' : ∀ p, SemistableAt p E')
    (hcop : ∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) :
    IsMuAtBadPrimes E E' l :=
  isMuAtBadPrimes_of_veluQuotient_of_coprime_K E E' hl Q hQ hE' hssE hssE' hcop

def isMuAtBadPrimes_of_veluQuotient_nodeg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(IsMuAtBadPrimes——原文に無い仮説は無い。★[L:ℚ]+1<l を取らない)",
    sectionId := "genell-lemma-3-5" }

def isMuAtBadPrimes_of_veluQuotient_nodeg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_of_coprime_K(第 1140、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_of_coprime_K") 1 ]

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**不等式（第 1141、★`[L:ℚ]+1<l` なし）**。 -/
theorem lemma_3_5_velu_large_K (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ∀ (P : (L →+* ℂ) → PeriodPair) (Cv : (L →+* ℂ) → VariableChange ℂ),
      (∀ σ, latticeDisc (P σ) ≠ 0) →
      (∀ σ, Cv σ • (E.map σ) = latticeCurve (P σ)) →
      (∀ σ : L →+* ℂ, (E.map σ).IsElliptic) →
      (∀ σ : L →+* ℂ, (Cv σ • (E.map σ)).IsElliptic) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := ABC3.Found.GaloisRep.lemma_3_5_velu_bad_delta eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint
    hssE hssE' hcop => ?_⟩
  exact hC L E E' l hl.pos Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint hssE hssE'
    (isMuAtBadPrimes_of_veluQuotient_nodeg E E' hl Q hQ hE' hssE hssE' hcop)

def lemma_3_5_velu_large_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(不等式——★[L:ℚ]+1<l を取らない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_large_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1 ]

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**アルキメデスの周期対を自前で作る形（第 1141）**。 -/
theorem lemma_3_5_velu_arch_K (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p : HeightOneSpectrum (𝓞 L), neronExp p E = 0) →
      (∀ p : HeightOneSpectrum (𝓞 L), E'.IsIntegral (primeSubring p)) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_large_K eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' hmin0 hint hssE hssE' hcop => ?_⟩
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_periodPair_family E
  exact hC L E E' l hl Q hQ hE' P Cv hΔ hPC hell1 hell2 hmin0 hint hssE hssE' hcop

def lemma_3_5_velu_arch_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(アルキメデスの周期対を自前で作る形。★[L:ℚ]+1<l を取らない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_arch_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1 ]

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**`hfin` で受ける形（第 1141）**。 -/
theorem lemma_3_5_velu_defect_K (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      ((∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
          - (∑ᶠ p : HeightOneSpectrum (𝓞 L),
            (neronExp p E' : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
        ≤ (3 / 2) * (Module.finrank ℚ L : ℝ) * Real.log l) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := ABC3.Found.GaloisRep.lemma_3_5_of_isogeny_estimate_le eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' hfin hssE hssE' hcop => ?_⟩
  obtain ⟨P, Cv, hΔ, hPC, hell1, hell2⟩ := exists_periodPair_family E
  refine hC L E E' l hssE' ?_
    (ABC3.Found.GaloisRep.htFalt_veluQuotientFull_le_of_defect E E' l hl.pos Q hQ hE'
      P Cv hΔ hPC hell1 hell2 hfin)
  refine ABC3.Found.GaloisRep.degInfOf_ge_of_local E E' l (fun p => ?_)
  exact ABC3.Found.GaloisRep.minDeltaExp_le_of_bad_delta p E E' (hssE p) l
    (fun hb => isMuAtBadPrimes_of_veluQuotient_nodeg E E' hl Q hQ hE'
      hssE hssE' hcop p hb)

def lemma_3_5_velu_defect_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(hfin で受ける形。★[L:ℚ]+1<l を取らない)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_defect_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1 ]

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★★★★★**[GenEll] `Lemma 3.5`——曲線の水準（第 1141、★`[L:ℚ]+1<l` なし）**。 -/
theorem lemma_3_5_velu_K (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E E' : WeierstrassCurve L)
      [E.IsElliptic] [E'.IsElliptic] (l : ℕ), l.Prime →
      ∀ Q : E.toAffine.Point, addOrderOf Q = l →
      E' = veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q))) →
      (∀ p, SemistableAt p E) →
      (∀ p, SemistableAt p E') →
      (∀ p : HeightOneSpectrum (𝓞 L), jExp p E < 0 → ¬ ((l : ℤ) ∣ jExp p E)) →
      (1 / (12 * (1 + eps))) * (l : ℝ) * degInfOf L E
        ≤ htFaltOf L E + 2 * Real.log l + C := by
  obtain ⟨C, hC⟩ := lemma_3_5_velu_defect_K eps heps
  refine ⟨C, fun L _ _ E E' _ _ l hl Q hQ hE' hssE hssE' hcop => ?_⟩
  exact hC L E E' l hl Q hQ hE'
    (ABC3.Found.GaloisRep.hfin_of_veluQuotientFull E E' hl Q hQ hE') hssE hssE' hcop

def lemma_3_5_velu_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の水準——★[L:ℚ]+1<l を取らない。残る逸脱は l ≠ 2 だけ)",
    sectionId := "genell-lemma-3-5" }

def lemma_3_5_velu_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMuAtBadPrimes_of_veluQuotient_nodeg(第 1141、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.isMuAtBadPrimes_of_veluQuotient_nodeg") 1 ]

end ABC3.Skeleton.GenEll
