/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim
import ABC3.Found.GaloisRep.RamifiedBadPrime.Definition33

/-!
# RamifiedBadPrime —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical
variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★★★★★★★★★★★★★★★★★
**`hcop` の出どころ——分岐した拡大でも `l ∤ e` なら通る**（第 1371）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `local_inputs_of_split`（第 1317）に渡す `hnd` である。 -/
theorem not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram [CharZero Lv]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) {l : ℕ} (hl : l.Prime) (hle : ¬ (l ∣ e))
    (hcop : ¬ ((l : ℤ) ∣ jExp p E)) :
    ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h) (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q) := by
  rw [vAdd_tateParam_eq_neg_mul_jExp p hpe E h hj]
  intro hdvd
  have hdvd' : (l : ℤ) ∣ (e : ℤ) * jExp p E := (dvd_neg).1 hdvd
  have hlp : Prime ((l : ℤ)) := Nat.prime_iff_prime_int.1 hl
  rcases hlp.dvd_mul.1 hdvd' with h1 | h2
  · exact hle (by exact_mod_cast h1)
  · exact hcop h2

def not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(hcop の出どころ。分岐した拡大でも l ∤ e なら通る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def not_dvd_vAdd_tateParam_of_not_dvd_jExp_ram.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "vAdd_algebraMap_eq_mul_valAdd(第 1370、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_mul_valAdd") 1,
    .citation "[ABC3]" "le_finrank_of_pow_eq_map_maximalIdeal(第 1369、証明済み。l ∤ e の出どころ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.le_finrank_of_pow_eq_map_maximalIdeal") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1371）**——第 932／第 978 の `e` 倍版である。" ++
       "☆`l ∤ e` は第 1369 の `e ≤ [L_v′ : L_v]` から出る——" ++
       "`L_p(ζ_l)/L_p` は `≤ l−1` 次、分裂用の 2 次拡大は `≤ 2 < 5 ≤ l`。") 19 ]

end ABC3.Found.GaloisRep
