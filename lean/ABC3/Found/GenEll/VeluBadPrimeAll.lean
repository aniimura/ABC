/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluSemistableBadMu
import ABC3.Found.GenEll.MuPrimitiveRootOrDeep
import ABC3.Found.GenEll.VeluNotDvdLFree
import ABC3.Meta.Claim

/-!
# 第 1428 ブロック —— **★★★★★★★★★★★★★★★★★★★★★★★★
悪い素点は `p ∣ l` でも無条件**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

第 1424 の二者択一の両側が閉じた:

| 核 | 閉じた場所 |
|---|---|
| 深い核 | 第 1412・1415・1420-1422（`hlu` を一度も使わない） |
| `μ_l` 型 | ★第 1425-1427（短 Weierstrass 形で `p ∣ l` を通す） |

☆したがって**悪い素点では `hlu`（`p ∤ l`）がまったく要らない**。
★残る条件は **`p ∤ 6`**（`h48`・`h864`）だけで、`p ∣ l` の場合は `l ≥ 5` を意味する。

| 定理 | 内容 |
|---|---|
| `semistableAt_veluQuotient_bad_ram_all` | ★★★★★★★★分裂の場合（分岐版）——`hlu` なし |
| `semistableAt_veluQuot_multRed_local_all` | ★★★★★★★★★★★★分裂性も仮定しない形 |

☆残るのは **良い素点で `p ∣ l`** の場合だけである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定**——★**`p ∣ l` を許す**（第 1428）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1424 の二者択一の両側を閉じた形である:

* 深い核——`v ∈ 𝔪` から `c₄` が単元（第 1412）、整性は第 1420-1422
* `μ_l` 型——`c₄ = l⁴·c₄(E_{q^l})` から短 Weierstrass 形で（第 1425-1427）

★残る条件は `p ∤ 6`（`h48`・`h864`）だけである。 -/
theorem semistableAt_veluQuotient_bad_ram_all {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (h2K : (2 : Lv) ≠ 0)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hΔL : E'.Δ ≠ 0) :
    SemistableAt p E' := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE, hcoords⟩ :=
    exists_variableChange_veluQuotient_tateModel_coords E E' h hl hQ h2K hE'
  have hu0 := vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  rcases ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl P hP hP0 with
    ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ | ⟨y, hPz, hdeep⟩
  · -- ★★`μ_l` 型（第 1427）
    exact semistableAt_veluQuotient_bad_mu he p hpe E E' h hl h2K h48 h864
      C₀ P hcurve hu0 ζ uζ hζ hζu hζl hord hPz
  · -- ★★★深い核（第 1424 と同じ）
    obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
    have hTmem := pointCoords_tatePhi_mem_of_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      y hdeep
    rw [← hPz] at hTmem
    have hmem := pointCoords_mem_primeSubring_of_image_mem he p hpe C₀ E Q _ hcoords hTmem
    haveI hint : (veluQuotientFull E (((range (2 * m + 1)).erase 0).image
        (fun k : ℕ => pointCoords (k • Q)))).IsIntegral (primeSubring p) :=
      isIntegral_veluQuotientFull_of_pointCoords_mem p E hl Q hQ hmem
    haveI hintE' : WeierstrassCurve.IsIntegral (primeSubring p) E' := by rw [hE']; exact hint
    have hyl : (2 * m + 1) • tatePhi
        (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
        (QuotientGroup.mk y) = 0 := by rw [← hPz]; exact hP
    obtain ⟨w, hw⟩ := exists_veluW_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      m y hyl hdeep
    have hquot := veluQuotientFull_tate_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      y hdeep _ w h2K rfl hw
    rw [hPz] at hcurve
    exact semistableAt_velu_of_veluCurve_eq_ram he p hpe E' hΔL
      (tateParamR (E.baseChange Lv) h) hq _ w
      (isUnit_c4_add_240_deep (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0)
        y hdeep _ rfl)
      C₀ (hcurve.symm.trans hquot) hu0

def semistableAt_veluQuotient_bad_ram_all.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——分岐版。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_ram_all.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuotient_bad_mu(第 1427、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_veluQuotient_bad_mu") 1,
    .citation "[ABC3]" "isUnit_c4_add_240_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_add_240_deep") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1428）**——第 1424 の二者択一の**両側が閉じた**。" ++
       "☆深い核は第 1412・1415・1420-1422、`μ_l` 型は第 1425-1427 である。" ++
       "★残る条件は `p ∤ 6` だけ——`p ∣ l` の場合は `l ≥ 5` を意味する。") 17 ]

end ABC3.Found.GaloisRep

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset Polynomial
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定**——★**分裂性も `p ∤ l` も仮定しない**（第 1428）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1417 の `semistableAt_veluQuot_multRed_local_free` から
**`hlu`（`p ∤ l`）が落ちた**形である。
★中身は第 1405 と同じ——非分裂なら不分岐 2 次拡大へ上げ、
第 1428 の `semistableAt_veluQuotient_bad_ram_all` を当てるだけ。

☆残る条件は `p ∤ 6`（`h48`・`h864`）だけである。 -/
theorem semistableAt_veluQuot_multRed_local_all
    (L : Type) [Field L] [NumberField L] [inst : DecidableEq L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [W.IsElliptic] (C : WeierstrassCurve.VariableChange L)
    [(C • W).IsElliptic]
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • W))
    (hc4ne : (C • W).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • W).c₄) hc4ne) = 0)
    [((C • W).baseChange Lv).IsElliptic]
    [WeierstrassCurve.HasMultiplicativeReduction R ((C • W).baseChange Lv)]
    (hj : jExp p W < 0) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (h2K : (2 : Lv) ≠ 0)
    (Q : (C • W).toAffine.Point) (hQ : addOrderOf Q = l)
    (hΔL : (veluQuotientFull (C • W)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).Δ ≠ 0) :
    SemistableAt p (veluQuotientFull (C • W)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  haveI hCint : WeierstrassCurve.IsIntegral (primeSubring p) (C • W) := by
    haveI := hC; infer_instance
  haveI hVell := isElliptic_veluQuotientFull_nsmul_nf' L (C • W) hQ
  have hA := integralModel_c4_isUnit (R := R) ((C • W).baseChange Lv)
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv)
  · exact semistableAt_veluQuotient_bad_ram_all he p hpe (C • W) _ hs hl hodd h2K h48 h864
      hQ rfl hΔL
  · set F := splitQuadPoly (WeierstrassCurve.integralModel R ((C • W).baseChange Lv)) hA
      with hFdef
    have hfm : F.Monic := monic_splitQuadPoly _ hA
    have hfd : F.natDegree = 2 := natDegree_splitQuadPoly _ hA
    have hirr := irreducible_map_residue_of_not_splits hfm hfd
      (fun hsp => hs (hasSplit_of_splits_splitQuadPoly _ hA hsp))
    haveI hdom := isDomain_adjoinRoot hfm hirr
    haveI hlocR := isLocalRing_adjoinRoot hfm hfd hirr
    haveI hdvr := isDiscreteValuationRing_adjoinRoot hfm hfd hirr
    haveI hac := isAdicComplete_maximalIdeal_adjoinRoot hfm hfd hirr
    haveI hcz : CharZero (AdjoinRoot F) :=
      charZero_of_injective_algebraMap (algebraMap_adjoinRoot_injective (R := R) hfm hfd)
    haveI hcz2 : CharZero (FractionRing (AdjoinRoot F)) := by
      refine charZero_of_injective_algebraMap (R := AdjoinRoot F) ?_
      exact IsFractionRing.injective (AdjoinRoot F) (FractionRing (AdjoinRoot F))
    letI : Algebra L (FractionRing (AdjoinRoot F)) :=
      ((quadFieldHom (K := Lv) hfm hfd).comp (algebraMap L Lv)).toAlgebra
    have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot F)) x
        = quadFieldHom hfm hfd (algebraMap L Lv x) := fun _ => rfl
    haveI hmin' := isMinimal_baseChange_ext_ram p hpe hfm hfd hirr halg W C hC hc4ne hc4
    haveI hmult' :=
      hasMultiplicativeReduction_ext_ram he p hpe hfm hfd hirr halg W C hC hc4ne hc4 hj
    have h' := hasSplitMultiplicativeReduction_ext (C • W) hA hFdef halg
    haveI hell' : ((C • W).baseChange (FractionRing (AdjoinRoot F))).IsElliptic :=
      isElliptic_baseChange (C • W)
    have h2K' : (2 : FractionRing (AdjoinRoot F)) ≠ 0 := two_ne_zero
    exact semistableAt_veluQuotient_bad_ram_all
      (Lv := FractionRing (AdjoinRoot F)) (R := AdjoinRoot F) he p
      (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) (C • W) _ h' hl hodd h2K'
      h48 h864 hQ rfl hΔL

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_multRed_local_all.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——分裂性も p ∤ l も不要)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_multRed_local_all.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuotient_bad_ram_all(第 1428、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_veluQuotient_bad_ram_all") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1428）**——第 1417 から `hlu`（`p ∤ l`）が落ちた。" ++
       "☆残る条件は `p ∤ 6` だけである。" ++
       "★★★これで**悪い素点は `p ∣ l` でも閉じた**" ++
       "——残るのは良い素点で `p ∣ l` の場合だけである。") 17 ]

end ABC3.Found.GenEll
