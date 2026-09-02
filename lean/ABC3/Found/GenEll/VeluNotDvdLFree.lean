/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluNotDvdL
import ABC3.Found.GaloisRep.VeluSemistableBadFree
import ABC3.Meta.Claim

/-!
# 第 1417 ブロック —— **★★★★★★★★★★★★★★★★★★★★★★★★
`p ∤ l` の側が完全に無条件で閉じた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

第 1407 の `semistableAt_veluQuot_of_not_dvd` は
`hcop`（`l ∤ jExp p E`）を仮定していた。
★第 1408 で分かったとおり `hcop` は**底変換で保たれない**ので、
`VeluQuotOK`（任意の有限拡大の素点で半安定性を要求する）には使えなかった。

☆第 1416 で局所の側から `hcop` が落ちたので、その伝播を 3 段行う:

| # | 定理 | 落ちた仮定 |
|---|---|---|
| 1 | `semistableAt_veluQuot_multRed_local_free` | `hle`（`l ∤ e`）・`hcop` |
| 2 | `semistableAt_veluQuot_badPrime_free` | `hcop` |
| 3 | `semistableAt_veluQuot_of_not_dvd_free` | `hcop` |

★★★★これで **`p ∤ l` の側は完全に無条件**になった——残る仮定は

* `SemistableAt p E`（もとの曲線が半安定）
* `l` が奇素数、`p ∤ l`（`hlu`）
* `addOrderOf Q = l`

だけである。☆`VeluQuotOK` の側に通すべき仮定は**もう無い**。

## ☆残るのは `p ∣ l` だけ

`blocked-leaves.json` の `pDivLRevised2026_09_02`（形式群）と
`pDivLDualIdea2026_09_02`（双対同種）を参照。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset Polynomial
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-! ## ★★★★段 1——局所（分裂性を仮定しない形） -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定**——★**`hcop` も `l ∤ e` もなし**（第 1417）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1405 の `semistableAt_veluQuot_multRed_local` から
`hle`（`l ∤ e`）と `hcop`（`l ∤ jExp p W`）が落ちた形である。

★中身は第 1405 と同じ——非分裂なら不分岐 2 次拡大へ上げ、
第 1416 の `semistableAt_veluQuotient_bad_ram_free` を当てるだけ。 -/
theorem semistableAt_veluQuot_multRed_local_free
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
  have hA := integralModel_c4_isUnit (R := R) ((C • W).baseChange Lv)
  by_cases hs : WeierstrassCurve.HasSplitMultiplicativeReduction R ((C • W).baseChange Lv)
  · exact semistableAt_veluQuotient_bad_ram_free he p hpe (C • W) _ hs hl hodd
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
    exact semistableAt_veluQuotient_bad_ram_free
      (Lv := FractionRing (AdjoinRoot F)) (R := AdjoinRoot F) he p
      (valuation_algebraMap_ext_ram p hpe hfm hfd hirr halg) (C • W) _ h' hl hodd
      hluR' h2R' h2K' hQ rfl hΔL

/-! ## ★★★★★★★★段 2——悪い素点（大域） -/

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点（`p ∤ l`）では Vélu の商は半安定**——★**無条件**（第 1417）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1406 の `semistableAt_veluQuot_badPrime` から `hcop` が落ちた形である。 -/
theorem semistableAt_veluQuot_badPrime_free [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  obtain ⟨C, hC, hc4ne, hc4⟩ := exists_minimal_c4_unit_of_jExp_neg p E hss hj
  obtain ⟨e, he1, hle, hIe⟩ := exists_locCyc_package p hl
  have hpe := valuation_algebraMap_locCycField p hl he1 hIe
  haveI hCell : (C • E).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero E hΔ C)⟩
  haveI hbcell : ((C • E).baseChange (locCycField p hl)).IsElliptic :=
    isElliptic_baseChange (C • E)
  haveI hmin := isMinimal_baseChange_at_bad_prime_ram (Lv := locCycField p hl)
    (R := locCycRing p hl) p hpe E C hC hc4ne hc4
  haveI hmult := hasMultiplicativeReduction_at_bad_prime_ram (Lv := locCycField p hl)
    (R := locCycRing p hl) he1 p hpe E C hC hc4ne hc4 hj
  have hlL : ((l : L)) ≠ 0 := Nat.cast_ne_zero.2 hl.ne_zero
  have hval := valAdd_natCast_eq_zero_of_isUnit p hlL hlu
  have hlv : ((l : locCycField p hl)) ≠ 0 := Nat.cast_ne_zero.2 hl.ne_zero
  have hluR : IsUnit ((l : locCycRing p hl)) :=
    isUnit_natCast_of_valAdd_eq_zero_ram p hpe l hlL hlv hval
  have h2R : (2 : locCycRing p hl) ≠ 0 := two_ne_zero
  have h2K : (2 : locCycField p hl) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hVell := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  have hΔL' : (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).Δ ≠ 0 :=
    (veluQuotientFull (C • E)
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).isUnit_Δ.ne_zero
  have hres := semistableAt_veluQuot_multRed_local_free L he1 p hpe E C hC hc4ne hc4 hj hl hodd
    hlu hluR h2R h2K (vcPoint C E Q) hQ' hΔL'
  rw [← heq] at hres
  have hfin := semistableAt_variableChange p _ C⁻¹ hres
  rwa [inv_smul_smul] at hfin

/-! ## ★★★★★★★★★★★★★★★★段 3——`p ∤ l` の側が閉じた -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] `p ∤ l` では Vélu の商は半安定**——★★**無条件**（第 1417）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★**2026-09-02（第 1417）**——第 1407 から `hcop`（`l ∤ jExp p E`）が落ちた。
☆残る仮定は `SemistableAt p E`・`l` が奇素数・`p ∤ l`・`addOrderOf Q = l` だけである。

★★★これは `VeluQuotOK`（**任意の**有限拡大の素点で半安定性を要求する）に
そのまま通せる形である——底変換で壊れる仮定を一つも使っていない。 -/
theorem semistableAt_veluQuot_of_not_dvd_free [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  by_cases hj : 0 ≤ jExp p E
  · exact semistableAt_veluQuot_goodPrime p E hss hj hl hodd hlu Q hQ
  · exact semistableAt_veluQuot_badPrime_free p E hss (by omega) hl hodd hlu Q hQ

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_multRed_local_free.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——局所・★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_badPrime_free.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点 p ∤ l では Vélu の商は半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_of_not_dvd_free.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p ∤ l では Vélu の商は半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_of_not_dvd_free.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_goodPrime(第 1403、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_goodPrime") 1,
    .citation "[ABC3]" "semistableAt_veluQuotient_bad_ram_free(第 1416、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_veluQuotient_bad_ram_free") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1417）**——第 1407 から `hcop`（`l ∤ jExp p E`）が落ち、" ++
       "**`p ∤ l` の側が完全に無条件で閉じた**。" ++
       "☆残る仮定は `SemistableAt p E`・`l` が奇素数・`p ∤ l`・`addOrderOf Q = l` だけ。" ++
       "★★★これは `VeluQuotOK`（任意の有限拡大の素点で半安定性を要求する）に" ++
       "そのまま通せる形である——底変換で壊れる仮定を一つも使っていない" ++
       "（`hcopDoesNotDescend2026_09_02` の問題が解消した）。" ++
       "☆残るのは `p ∣ l` だけである。") 17 ]

end ABC3.Found.GenEll
