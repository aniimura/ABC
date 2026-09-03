/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluSemistableBadRam
import ABC3.Found.GenEll.AlphaSplitDichotomy
import ABC3.Found.GaloisRep.RamifiedBadPrime
import ABC3.Found.GenEll.VeluGoodPrime
import ABC3.Meta.Claim

/-!
# 第 1405 ブロック —— **悪い素点：分裂性を仮定しない形**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——**分裂の二者択一**

第 1404 で悪い素点の半安定性は分岐版になったが、**分裂乗法還元**を仮定していた。
★実際に手に入るのは `HasMultiplicativeReduction` だけである。

☆本ブロックは α の葉の第 1383（`exists_h2_h1_of_multRed_local`）と**同型の議論**で
分裂性を落とす:

* 分裂しているならそのまま第 1404
* 分裂していないなら `R[X]/(f)`（`f = splitQuadPoly`）へ上げる
  ——**不分岐**な 2 次拡大なので `hpe` の `e` は変わらない（第 1382）

★★★これで悪い素点に残るのは **`l ∣ v_p(q)` の場合**（核が `μ_l` 型でない）だけになる。
☆その場合は商が `q^{1/l}` になり、やはり乗法還元だが別の Tate の計算が要る。
★消費側（`hdag_of_stableLine`）は `PrimeToLocalHeights l`（＝`l ∤ jExp p`）を持っているので、
`VeluQuotOK` の側にそれを通せば避けられる可能性がある。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset Polynomial
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定**——★**分裂性を仮定しない**（第 1405）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮定に残るのは `l ∤ jExp p W`（＝`PrimeToLocalHeights`）と `l ∤ e` だけである。

★★★α の葉の第 1383 と**同型の議論**——非分裂なら不分岐 2 次拡大へ上げる。 -/
theorem semistableAt_veluQuot_multRed_local
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
    (hj : jExp p W < 0) {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p W))
    (hlu : IsUnit ((l : primeSubring p)))
    (hluR : IsUnit ((l : R))) (h2R : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
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
  haveI hint := isIntegral_veluQuotientFull_of_addOrderOf_prime p (C • W) hl hlu Q hQ
  haveI hVell := isElliptic_veluQuotientFull_nsmul_nf' L (C • W) hQ
  have hjC : jExp p (C • W) < 0 := by rw [jExp_variableChange]; exact hj
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • W)) := by rw [jExp_variableChange]; exact hcop
  have hA := integralModel_c4_isUnit (R := R) ((C • W).baseChange Lv)
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv)
  · exact semistableAt_veluQuotient_bad_ram he p hpe (C • W) _ hs hl hodd
      (not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram p hpe (C • W) hs hjC hl hle hcopC)
      hluR h2R h2K hQ rfl hΔL
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
    have hluR' : IsUnit ((l : AdjoinRoot F)) := by
      have h := hluR.map (algebraMap R (AdjoinRoot F))
      simpa using h
    have h2R' : (2 : AdjoinRoot F) ≠ 0 := two_ne_zero
    have h2K' : (2 : FractionRing (AdjoinRoot F)) ≠ 0 := two_ne_zero
    exact semistableAt_veluQuotient_bad_ram
      (Lv := FractionRing (AdjoinRoot F)) (R := AdjoinRoot F) he p
      (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) (C • W) _ h' hl hodd
      (not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram p
        (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) (C • W) h' hjC hl hle hcopC)
      hluR' h2R' h2K' hQ rfl hΔL

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_multRed_local.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——分裂性を仮定しない。★l ∤ jExp p)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_multRed_local.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuotient_bad_ram(第 1404、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_veluQuotient_bad_ram") 1,
    .citation "[ABC3]" "hasSplitMultiplicativeReduction_ext(第 1033、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasSplitMultiplicativeReduction_ext") 1,
    .citation "[ABC3]" "not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram(第 1371、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1405）**——α の葉の第 1383 と同型の議論で分裂性を落とした。" ++
       "☆非分裂なら `R[X]/(f)` へ上げる——不分岐なので `hpe` の `e` は変わらない。" ++
       "★★★これで悪い素点に残るのは **`l ∣ v_p(q)` の場合**だけになる" ++
       "（核が `μ_l` 型でない。商は `q^{1/l}` でやはり乗法還元だが別の Tate の計算が要る）。" ++
       "☆消費側の `hdag_of_stableLine` は `PrimeToLocalHeights l` を持っているので、" ++
       "`VeluQuotOK` の側にそれを通せば避けられる可能性がある。") 17 ]

end ABC3.Found.GenEll
