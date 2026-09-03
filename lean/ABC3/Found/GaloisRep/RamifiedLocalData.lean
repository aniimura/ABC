/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.RamifiedBadPrime
import ABC3.Meta.Claim

/-!
# 第 1375 ブロック —— **拡大の上の局所データ（分岐版）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1372 が要るインスタンス引数

第 1372（`exists_h2_h1_of_bad_prime_ram`）は

* `[WeierstrassCurve.IsMinimal R (E.baseChange Lv)]`
* `HasSplitMultiplicativeReduction R (E.baseChange Lv)`

を**インスタンス引数・仮説**として受ける。
☆不分岐版（`hp`）は第 908／第 973／第 976 にあるが、
分岐版（`hpe`）が要る——★本ブロックはそれを作る。

☆どれも `hp` の使い所は 2 種類しかない:

| 使い所 | 不分岐 | 分岐（`e` 倍） |
|---|---|---|
| `v_p(x) ≤ 1 ⟹ v_R(φ x) ≤ 1` | 等号 | `t ≤ 1 ⟹ t^e ≤ 1` |
| 付値指数の一致 | `vAdd = valAdd` | `vAdd = e·valAdd` |

★`v_p(c₄) = 0` は `e·0 = 0` で保たれ、`0 < v_p(Δ)` は `e ≥ 1` で保たれる。
☆したがって**極小性も乗法還元も分岐で壊れない**。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField

open scoped Classical

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★★★**`p` の整数環の元は `R` から来る**——分岐版（第 1375）。 -/
theorem exists_preimage_of_mem_primeSubring_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (x : L) (hx : x ∈ primeSubring p) :
    ∃ r : R, algebraMap R Lv r = algebraMap L Lv x := by
  refine exists_preimage_of_valuation_le_one (R := R) _ ?_
  rw [hpe x]
  exact pow_le_one₀ zero_le ((mem_primeSubring_iff p x).1 hx)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`p` 上整な曲線は `R` 上整**——分岐版（第 1375）。 -/
theorem isIntegral_baseChange_of_isIntegral_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (V : WeierstrassCurve L) [hV : WeierstrassCurve.IsIntegral (primeSubring p) V] :
    WeierstrassCurve.IsIntegral R (V.baseChange Lv) := by
  obtain ⟨r1, hr1⟩ := exists_preimage_of_mem_primeSubring_ram (Lv := Lv) (R := R) p hpe V.a₁
    (by rw [← WeierstrassCurve.integralModel_a₁_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₁).2)
  obtain ⟨r2, hr2⟩ := exists_preimage_of_mem_primeSubring_ram (Lv := Lv) (R := R) p hpe V.a₂
    (by rw [← WeierstrassCurve.integralModel_a₂_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₂).2)
  obtain ⟨r3, hr3⟩ := exists_preimage_of_mem_primeSubring_ram (Lv := Lv) (R := R) p hpe V.a₃
    (by rw [← WeierstrassCurve.integralModel_a₃_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₃).2)
  obtain ⟨r4, hr4⟩ := exists_preimage_of_mem_primeSubring_ram (Lv := Lv) (R := R) p hpe V.a₄
    (by rw [← WeierstrassCurve.integralModel_a₄_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₄).2)
  obtain ⟨r6, hr6⟩ := exists_preimage_of_mem_primeSubring_ram (Lv := Lv) (R := R) p hpe V.a₆
    (by rw [← WeierstrassCurve.integralModel_a₆_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₆).2)
  refine ⟨⟨⟨r1, r2, r3, r4, r6⟩, ?_⟩⟩
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · show algebraMap L Lv V.a₁ = algebraMap R Lv r1
    exact hr1.symm
  · show algebraMap L Lv V.a₂ = algebraMap R Lv r2
    exact hr2.symm
  · show algebraMap L Lv V.a₃ = algebraMap R Lv r3
    exact hr3.symm
  · show algebraMap L Lv V.a₄ = algebraMap R Lv r4
    exact hr4.symm
  · show algebraMap L Lv V.a₆ = algebraMap R Lv r6
    exact hr6.symm

/-- ★★★★★★★★★★**極小モデルは拡大でも極小**——分岐版（第 1375）。 -/
theorem isMinimal_baseChange_of_c4_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) :
    (W.baseChange Lv).IsMinimal R := by
  haveI hint : WeierstrassCurve.IsIntegral R (W.baseChange Lv) :=
    isIntegral_baseChange_of_isIntegral_ram p hpe W
  have hinj : Function.Injective (algebraMap L Lv) := (algebraMap L Lv).injective
  have hne2 : algebraMap L Lv W.c₄ ≠ 0 := (map_ne_zero_iff _ hinj).2 hc4ne
  have hne3 : algebraMap L Lv W.Δ ≠ 0 := (map_ne_zero_iff _ hinj).2 hΔ
  have hc4eq : (W.baseChange Lv).c₄ = algebraMap L Lv W.c₄ := WeierstrassCurve.map_c₄ _ _
  have hΔeq : (W.baseChange Lv).Δ = algebraMap L Lv W.Δ := WeierstrassCurve.map_Δ _ _
  have hc4ne' : (W.baseChange Lv).c₄ ≠ 0 := by rw [hc4eq]; exact hne2
  have hΔne' : (W.baseChange Lv).Δ ≠ 0 := by rw [hΔeq]; exact hne3
  refine isMinimal_of_c4_vAdd_eq_zero (W.baseChange Lv) hΔne' hc4ne' ?_
  have hu : Units.mk0 ((W.baseChange Lv).c₄) hc4ne'
      = Units.mk0 (algebraMap L Lv W.c₄) hne2 := Units.ext hc4eq
  rw [hu, vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe W.c₄ hc4ne hne2, hc4, mul_zero]

/-- ★★★★★★★★★★★★**悪い素点では拡大でも乗法還元をもつ**——分岐版（第 1375）。 -/
theorem hasMultiplicativeReduction_baseChange_ram {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0)
    (hΔpos : 0 < valAdd p (Units.mk0 W.Δ hΔ))
    [(W.baseChange Lv).IsMinimal R] :
    (W.baseChange Lv).HasMultiplicativeReduction R := by
  have he0 : (0 : ℤ) < (e : ℤ) := by exact_mod_cast he
  have hinj : Function.Injective (algebraMap L Lv) := (algebraMap L Lv).injective
  have hne2 : algebraMap L Lv W.c₄ ≠ 0 := (map_ne_zero_iff _ hinj).2 hc4ne
  have hne3 : algebraMap L Lv W.Δ ≠ 0 := (map_ne_zero_iff _ hinj).2 hΔ
  have hc4eq : (W.baseChange Lv).c₄ = algebraMap L Lv W.c₄ := WeierstrassCurve.map_c₄ _ _
  have hΔeq : (W.baseChange Lv).Δ = algebraMap L Lv W.Δ := WeierstrassCurve.map_Δ _ _
  have hc4ne' : (W.baseChange Lv).c₄ ≠ 0 := by rw [hc4eq]; exact hne2
  have hΔne' : (W.baseChange Lv).Δ ≠ 0 := by rw [hΔeq]; exact hne3
  have huc : Units.mk0 ((W.baseChange Lv).c₄) hc4ne'
      = Units.mk0 (algebraMap L Lv W.c₄) hne2 := Units.ext hc4eq
  have hud : Units.mk0 ((W.baseChange Lv).Δ) hΔne'
      = Units.mk0 (algebraMap L Lv W.Δ) hne3 := Units.ext hΔeq
  refine ⟨?_, ?_⟩
  · have hv : 0 < vAdd (tateDvrVal R Lv) (Units.mk0 ((W.baseChange Lv).Δ) hΔne') := by
      rw [hud, vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe W.Δ hΔ hne3]
      exact mul_pos he0 hΔpos
    have := (tateDvrVal_pos_iff (R := R) (K := Lv)
      (Units.mk0 ((W.baseChange Lv).Δ) hΔne')).1 hv
    simpa using this
  · have hv : vAdd (tateDvrVal R Lv) (Units.mk0 ((W.baseChange Lv).c₄) hc4ne') = 0 := by
      rw [huc, vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe W.c₄ hc4ne hne2, hc4, mul_zero]
    have := valuation_eq_one_of_vAdd_eq_zero (R := R)
      (Units.mk0 ((W.baseChange Lv).c₄) hc4ne') hv
    simpa using this

/-- ★★★★★★★★**極小モデルは拡大でも極小**——`C • E` の形（第 1375）。 -/
theorem isMinimal_baseChange_at_bad_prime_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0) :
    WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv) := by
  haveI := hC
  exact isMinimal_baseChange_of_c4_ram p hpe (C • E)
    (variableChange_Delta_ne_zero E (E.isUnit_Δ.ne_zero) C) hc4ne hc4

/-- ★★★★★★★★★★★★**悪い素点では拡大でも乗法還元**——`C • E` の形（第 1375）。 -/
theorem hasMultiplicativeReduction_at_bad_prime_ram {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (hj : jExp p E < 0)
    [WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv)] :
    WeierstrassCurve.HasMultiplicativeReduction R ((C • E).baseChange Lv) := by
  haveI := hC
  have hΔ : (C • E).Δ ≠ 0 := variableChange_Delta_ne_zero E (E.isUnit_Δ.ne_zero) C
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hj
  exact hasMultiplicativeReduction_baseChange_ram he p hpe (C • E) hΔ hc4ne hc4
    (valAdd_Delta_pos_of_jExp_neg p (C • E) hΔ hc4ne hc4 hjC)

/-! ## ★出典の紐付け(`.src`) -/

def isMinimal_baseChange_of_c4_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小モデルは拡大でも極小。分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def hasMultiplicativeReduction_baseChange_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では拡大でも乗法還元。分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_baseChange_of_isIntegral_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(p 上整な曲線は R 上整。分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def hasMultiplicativeReduction_baseChange_ram.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "vAdd_algebraMap_eq_mul_valAdd(第 1370、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_mul_valAdd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1375）**——`v_p(c₄) = 0` は `e·0 = 0` で保たれ、" ++
       "`0 < v_p(Δ)` は `e ≥ 1` で保たれる。" ++
       "☆したがって**極小性も乗法還元も分岐で壊れない**。") 17 ]

end ABC3.Found.GaloisRep
