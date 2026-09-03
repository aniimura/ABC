/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluBadPrimeAll
import ABC3.Found.GenEll.VeluGoodPrimeMem
import ABC3.Found.GenEll.VeluNotDvdLFree
import ABC3.Found.GenEll.PrimeOverL
import ABC3.Meta.Claim

/-!
# 第 1436 ブロック —— **★★★★★★★★★★★★★★★★★★★★★★★★
Vélu の商の半安定性が仮定 1 本になった**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か

第 1410-1435 で積み上げたものを **1 本の定理**に組み立てる:

| 場合 | 使うもの | 状態 |
|---|---|---|
| `p ∤ l` | 第 1417 `semistableAt_veluQuot_of_not_dvd_free` | ★**無条件** |
| `p ∣ l`・悪い素点 | 第 1428 `semistableAt_veluQuot_multRed_local_all` | ★**無条件**（`p ∤ 6` は第 1435） |
| `p ∣ l`・良い素点 | 第 1432 `semistableAt_veluQuot_good_of_mem` | ☆**核の座標の整性**が要る |

☆`p ∣ l` かつ `l ≥ 5` なら `p ∤ 6` は自動である（第 1435）。

★★★したがって `semistableAt_veluQuot_all_of_goodMem` の仮定は

* `SemistableAt p E`
* `l` が素数で `5 ≤ l`（界面の側から取れる——第 1434 の測定）
* `addOrderOf Q = l`
* ★☆**`p ∣ l` かつ良い素点のときだけ**、核の座標が `p` で整であること

の 4 つだけになった。☆最後の 1 つが形式群（`Ê(𝔪)[l] ∩ ⟨Q⟩ = 0`）であり、
第 1431 で **`jExp p E′ ≥ 0` 1 本**にも帰着させてある。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Skeleton.GenEll ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★★★悪い素点（`p ∣ l` も許す） -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**悪い素点では Vélu の商は半安定**——★**`p ∣ l` を許す**（第 1436）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1417 の `semistableAt_veluQuot_badPrime_free` から `hlu`（`p ∤ l`）が落ち、
代わりに `p ∤ 6`（`h48`・`h864`）だけが要る形である。 -/
theorem semistableAt_veluQuot_badPrime_all [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
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
  have h2K : (2 : locCycField p hl) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hVell := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  have hΔL' : (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).Δ ≠ 0 :=
    (veluQuotientFull (C • E)
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).isUnit_Δ.ne_zero
  have hres := semistableAt_veluQuot_multRed_local_all L he1 p hpe E C hC hc4ne hc4 hj hl hodd
    h48 h864 h2K (vcPoint C E Q) hQ' hΔL'
  rw [← heq] at hres
  have hfin := semistableAt_variableChange p _ C⁻¹ hres
  rwa [inv_smul_smul] at hfin

/-! ## ★★★★★★★★良い素点（核の座標が整なら `p ∣ l` も許す） -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★
**良い素点では核の座標が整なら Vélu の商は半安定**——★**`p ∣ l` を許す**（第 1436）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1403 の `semistableAt_veluQuot_goodPrime` の `hlu` を、
極小モデルの上の核の座標の整性に置き換えた形である。 -/
theorem semistableAt_veluQuot_goodPrime_of_mem [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : 0 ≤ jExp p E)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hmem : ∀ (C : WeierstrassCurve.VariableChange L),
      WeierstrassCurve.IsMinimal (primeSubring p) (C • E) →
      ∀ k ∈ (range l).erase 0,
        (pointCoords (k • vcPoint C E Q)).1 ∈ primeSubring p ∧
          (pointCoords (k • vcPoint C E Q)).2 ∈ primeSubring p) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hinst : inst = fun a b => Classical.propDecidable (a = b) := by
    funext a b
    exact Subsingleton.elim _ _
  subst hinst
  have hΔ : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  haveI := hC
  have hmd : minDeltaExp p E = 0 := by
    rw [minDeltaExp_eq_maxJ_of_semistable p E hss]
    exact max_eq_left (by omega)
  haveI hCell : (C • E).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (variableChange_Delta_ne_zero E hΔ C)⟩
  have hΔC : valAdd p (Units.mk0 ((C • E).Δ) (C • E).isUnit_Δ.ne_zero) = 0 := by
    have h := minDeltaExp_eq p E hΔ C hC
    rw [hmd] at h
    rw [← h]
  have hQ' : addOrderOf (vcPoint C E Q) = l := by rw [addOrderOf_vcPoint, hQ]
  haveI hell' := isElliptic_veluQuotientFull_nsmul_nf' L (C • E) hQ'
  have hΔ' : (veluQuotientFull (C • E)
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).Δ ≠ 0 :=
    (veluQuotientFull (C • E)
      (((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).isUnit_Δ.ne_zero
  have hgood := semistableAt_veluQuot_good_of_mem p (C • E) hΔC hl hodd
    (vcPoint C E Q) hQ' (hmem C hC) hΔ'
  have heq := veluQuotientFull_vcPoint_eq C E _ hQ two_ne_zero rfl
  rw [← heq] at hgood
  have hres := semistableAt_variableChange p _ C⁻¹ hgood
  rwa [inv_smul_smul] at hres

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★組み立て -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Vélu の商は半安定**——★**残る仮定は 1 本**（第 1436）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮定は

* `SemistableAt p E`
* `l` が素数で `5 ≤ l`
* `addOrderOf Q = l`
* ★☆**`p ∣ l` かつ良い素点のときだけ**、極小モデルの上で核の座標が `p` で整であること

の 4 つだけである。

★★★最後の 1 つが形式群（`Ê(𝔪)[l] ∩ ⟨Q⟩ = 0`、局所体の分岐が `e ≥ l−1` のときにしか
破れない）であり、第 1431 で **`jExp p E′ ≥ 0` 1 本**にも帰着させてある。 -/
theorem semistableAt_veluQuot_all_of_goodMem [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E)
    {l : ℕ} (hl : l.Prime) (hl5 : 5 ≤ l)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hmem : ¬ IsUnit ((l : primeSubring p)) → 0 ≤ jExp p E →
      ∀ (C : WeierstrassCurve.VariableChange L),
        WeierstrassCurve.IsMinimal (primeSubring p) (C • E) →
        ∀ k ∈ (range l).erase 0,
          (pointCoords (k • vcPoint C E Q)).1 ∈ primeSubring p ∧
            (pointCoords (k • vcPoint C E Q)).2 ∈ primeSubring p) :
    SemistableAt p (veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) := by
  have hodd : l ≠ 2 := by omega
  by_cases hlu : IsUnit ((l : primeSubring p))
  · -- ★`p ∤ l` は第 1417 で無条件
    exact semistableAt_veluQuot_of_not_dvd_free p E hss hl hodd hlu Q hQ
  · -- ★`p ∣ l` なら `l ≥ 5` から `p ∤ 6`（第 1435）
    have h48 := valAdd_48_eq_zero p hl hl5 hlu
    have h864 := valAdd_864_eq_zero p hl hl5 hlu
    by_cases hj : 0 ≤ jExp p E
    · exact semistableAt_veluQuot_goodPrime_of_mem p E hss hj hl hodd Q hQ (hmem hlu hj)
    · exact semistableAt_veluQuot_badPrime_all p E hss (by omega) hl hodd h48 h864 Q hQ

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuot_badPrime_all.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_goodPrime_of_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では核の座標が整なら Vélu の商は半安定。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_all_of_goodMem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は半安定。★残る仮定は良い素点 p ∣ l の核の整性 1 本)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_all_of_goodMem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_of_not_dvd_free(第 1417、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_of_not_dvd_free") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_multRed_local_all(第 1428、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_multRed_local_all") 1,
    .citation "[ABC3]" "semistableAt_veluQuot_good_of_mem(第 1432、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_good_of_mem") 1,
    .citation "[ABC3]" "valAdd_48_eq_zero・valAdd_864_eq_zero(第 1435、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.valAdd_48_eq_zero") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1436）**——第 1410-1435 を組み立てて、" ++
       "`Lemma 3.5` の半安定性に残る仮定を**1 本**にした。" ++
       "☆残るのは「`p ∣ l` かつ良い素点のとき核の座標が `p` で整」だけであり、" ++
       "それは形式群の `l`-捩れ（`Ê(𝔪)[l] ∩ ⟨Q⟩ = 0`）である" ++
       "——局所体の分岐が `e ≥ l−1` のときにしか破れない。" ++
       "★第 1431 で同じ穴を **`jExp p E′ ≥ 0`**（`j` の整性、" ++
       "古典的にはモジュラー多項式）にも帰着させてある。" ++
       "☆`5 ≤ l` は界面の側から取れる（第 1434 の測定）。") 17 ]

end ABC3.Found.GenEll
