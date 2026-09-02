/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluBadPrimeSplit
import ABC3.Found.GenEll.LocalCycPackage
import ABC3.Found.GaloisRep.UnitAtPrime
import ABC3.Meta.Claim

/-!
# 第 1406 ブロック —— **悪い素点の組み立て**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——局所パッケージを繋ぐ

第 1405 は局所体 `Lv`・`R` と `hpe` を仮定として受けていた。
★本ブロックは**与えられた悪い素点 `p` に対して**それを実際に作る:

| 部品 | 出どころ |
|---|---|
| 極小モデル `C` と `v_p(c₄) = 0` | `exists_minimal_c4_unit_of_jExp_neg`（第 954） |
| `L_p(ζ_l)` と `e`（`l ∤ e`） | `exists_locCyc_package`（第 1377） |
| 付値の比較 `v_{Lv}(x) = v_p(x)^e` | `valuation_algebraMap_locCycField`（第 1377） |
| `IsMinimal`・`HasMultiplicativeReduction` | 第 1375 |
| `IsUnit (l : R)` | ★本ブロック（`v_p(l) = 0` から、分岐版） |
| 極小モデルからの戻し | `veluQuotientFull_vcPoint_eq` ＋ `semistableAt_variableChange` |

★★★これで**悪い素点の半安定性が `p` ごとに閉じた**
——残る仮定は `l ∤ jExp p E`（＝`PrimeToLocalHeights`、消費側は持っている）だけである。

☆`c₄(E/C) ≠ 0` は第 1408 で仮定から落とした——Tate モデルとの変数変換
`(C₀ ⊗ Lv) • (E′ ⊗ Lv) = veluCurve(E_q, v, w) ⊗ Lv` と
`c₄(veluCurve) = c₄(E_q) + 240v` の単元性から出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★`v_p(n) = 0` なら `n` は局所環で単元（分岐版） -/

/-- ★★★★★★**`v_p(n) = 0` なら `n` は局所環で単元**——★分岐版（第 1406）。 -/
theorem isUnit_natCast_of_valAdd_eq_zero_ram
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] {e : ℕ}
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (n : ℕ) (hn : (n : L) ≠ 0) (hnv : (n : Lv) ≠ 0)
    (hval : valAdd p (Units.mk0 ((n : L)) hn) = 0) :
    IsUnit ((n : R)) := by
  have hmapL : algebraMap L Lv ((n : L)) = (n : Lv) := by rw [map_natCast]
  have hne' : algebraMap L Lv ((n : L)) ≠ 0 := by rw [hmapL]; exact hnv
  have hbridge := vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe ((n : L)) hn hne'
  rw [hval, mul_zero] at hbridge
  refine isUnit_of_vAdd_eq_zero (R := R) (Units.mk0 (algebraMap L Lv ((n : L))) hne')
    hbridge ((n : R)) ?_
  show algebraMap R Lv ((n : R)) = algebraMap L Lv ((n : L))
  rw [map_natCast, hmapL]

/-- ☆付値環で単元なら `valAdd` は `0`。 -/
theorem valAdd_natCast_eq_zero_of_isUnit (p : HeightOneSpectrum (𝓞 L)) {n : ℕ}
    (hn : (n : L) ≠ 0) (hu : IsUnit ((n : primeSubring p))) :
    valAdd p (Units.mk0 ((n : L)) hn) = 0 := by
  have hspec : ((hu.unit : primeSubring p) : L) = (n : L) := by
    rw [hu.unit_spec]; push_cast; ring
  have hmul : ((hu.unit : primeSubring p) : L) * ((↑hu.unit⁻¹ : primeSubring p) : L) = 1 := by
    have h : ((hu.unit : primeSubring p)) * ((↑hu.unit⁻¹ : primeSubring p)) = 1 := hu.unit.mul_inv
    exact_mod_cast congrArg (fun z : (primeSubring p) => (z : L)) h
  rw [hspec] at hmul
  have hinv : ((n : L))⁻¹ = ((↑hu.unit⁻¹ : primeSubring p) : L) :=
    (eq_inv_of_mul_eq_one_right hmul).symm
  have h1 : ((n : L)) ∈ primeSubring p := by
    rw [← hspec]; exact ((hu.unit : primeSubring p)).2
  have h2 : ((n : L))⁻¹ ∈ primeSubring p := by
    rw [hinv]; exact ((↑hu.unit⁻¹ : primeSubring p)).2
  have hge : 0 ≤ valAdd p (Units.mk0 ((n : L)) hn) :=
    (valAdd_nonneg_iff p _).2 ((mem_primeSubring_iff p _).1 h1)
  have hinvne : ((n : L))⁻¹ ≠ 0 := inv_ne_zero hn
  have hge2 : 0 ≤ valAdd p (Units.mk0 (((n : L))⁻¹) hinvne) :=
    (valAdd_nonneg_iff p _).2 ((mem_primeSubring_iff p _).1 h2)
  have heq : valAdd p (Units.mk0 (((n : L))⁻¹) hinvne)
      = - valAdd p (Units.mk0 ((n : L)) hn) := by
    rw [← valAdd_inv]
    exact valAdd_eq_of_valuation_eq p _ _ (by simp)
  omega

/-! ## ★★★★★★★★★★★★★★★★★★★★悪い素点の組み立て -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点（`p ∤ l`）では Vélu の商は半安定**——★（第 1406）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆仮定は `SemistableAt p E`・`jExp p E < 0`（悪い素点）・`l` 奇素数・`p ∤ l`・
`l ∤ jExp p E`（＝`PrimeToLocalHeights`）である。

★★★これで良い素点（第 1403）と合わせて **`p ∤ l` の側が閉じた**。 -/
theorem semistableAt_veluQuot_badPrime [inst : DecidableEq L]
    (p : HeightOneSpectrum (𝓞 L)) (E : WeierstrassCurve L) [E.IsElliptic]
    (hss : SemistableAt p E) (hj : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
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
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • vcPoint C E Q)))).isUnit_Δ.ne_zero
  have hres := semistableAt_veluQuot_multRed_local L he1 p hpe E C hC hc4ne hc4 hj hl hodd
    hle hcop hlu hluR h2R h2K (vcPoint C E Q) hQ' hΔL'
  rw [← heq] at hres
  have hfin := semistableAt_variableChange p _ C⁻¹ hres
  rwa [inv_smul_smul] at hfin

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_natCast_of_valAdd_eq_zero_ram.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(n) = 0 なら n は局所環で単元——分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valAdd_natCast_eq_zero_of_isUnit.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(付値環で単元なら valAdd は 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_badPrime.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点 p ∤ l では Vélu の商は半安定。★l ∤ jExp p のみ)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuot_badPrime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuot_multRed_local(第 1405、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_veluQuot_multRed_local") 1,
    .citation "[ABC3]" "exists_locCyc_package(第 1377、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_locCyc_package") 1,
    .citation "[ABC3]" "exists_minimal_c4_unit_of_jExp_neg(第 954、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1406）**——第 1405 の局所パッケージを `p` ごとに作った。" ++
       "★★★これで良い素点（第 1403）と合わせて **`p ∤ l` の側が閉じた**。" ++
       "☆残る仮定は `l ∤ jExp p E`（＝`PrimeToLocalHeights`、消費側は持っている）だけである" ++
       "——`c₄(E/C) ≠ 0` は第 1408 で落とした。") 17 ]

end ABC3.Found.GenEll
