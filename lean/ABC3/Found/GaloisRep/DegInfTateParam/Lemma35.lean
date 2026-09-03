/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.HtFinJ
import ABC3.Found.GaloisRep.DegInfTateParam.Lemma32

/-!
# DegInfTateParam —— `[GenEll] Lemma 3.5` の分

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

/-! ## ★★★★★★★★乗法還元なら極小性は完備化に移る -/

/-- ★★★★★★★★**`v_p(c₄) = 0` なら極小性は完備化に移る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-08-31（第 908）**——一般の極小性の移行は**弱近似を要する**が、
`Lemma 3.5` が実際に使うのは**乗法還元の場合だけ**である。
☆そこでは `v_p(c₄) = 0` なので、`isMinimal_of_c4_vAdd_eq_zero`（第 320）が
**直接**極小性を与える——弱近似は要らない。

★機構は 3 行である:

1. 整性は移る（`isIntegral_baseChange_of_isIntegral`、証明済み）
2. `v_p(c₄) = 0` は `Lv` の付値に移る（`vAdd_algebraMap_eq_valAdd`、証明済み）
3. `isMinimal_of_c4_vAdd_eq_zero`（証明済み） -/
theorem isMinimal_baseChange_of_c4 (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0) :
    (W.baseChange Lv).IsMinimal R := by
  haveI hint : WeierstrassCurve.IsIntegral R (W.baseChange Lv) :=
    isIntegral_baseChange_of_isIntegral p hp W
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
  rw [hu, vAdd_algebraMap_eq_valAdd (R := R) p hp W.c₄ hc4ne hne2]
  exact hc4

/-! ## ★★★★★★★★★★乗法還元も完備化に移る -/

/-- ★付値が `0` なら乗法的な付値は `1`。 -/
theorem valuation_eq_one_of_vAdd_eq_zero (x : Lvˣ)
    (h : vAdd (tateDvrVal R Lv) x = 0) :
    (IsDiscreteValuationRing.maximalIdeal R).valuation Lv (x : Lv) = 1 := by
  rw [valuation_eq_ofAdd_neg_vAdd (R := R), h]
  simp

/-- ★★★★★★★★★★**`v_p(c₄) = 0` と `v_p(Δ) > 0` なら、完備化で
乗法還元をもつ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが `SemistableAt`（悪い側）の内容そのものである——
極小モデルで `v_p(c₄) = 0`、そして `jExp p W < 0` なら `v_p(Δ_min) > 0`。

☆残るのは**分裂性**（剰余体で 2 次式が分裂すること）だけである。 -/
theorem hasMultiplicativeReduction_baseChange (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4ne : W.c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 W.c₄ hc4ne) = 0)
    (hΔpos : 0 < valAdd p (Units.mk0 W.Δ hΔ))
    [(W.baseChange Lv).IsMinimal R] :
    (W.baseChange Lv).HasMultiplicativeReduction R := by
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
  · -- ★`v(Δ) < 1`
    have hv : 0 < vAdd (tateDvrVal R Lv) (Units.mk0 ((W.baseChange Lv).Δ) hΔne') := by
      rw [hud, vAdd_algebraMap_eq_valAdd (R := R) p hp W.Δ hΔ hne3]
      exact hΔpos
    have := (tateDvrVal_pos_iff (R := R) (K := Lv)
      (Units.mk0 ((W.baseChange Lv).Δ) hΔne')).1 hv
    simpa using this
  · -- ★`v(c₄) = 1`
    have hv : vAdd (tateDvrVal R Lv) (Units.mk0 ((W.baseChange Lv).c₄) hc4ne') = 0 := by
      rw [huc, vAdd_algebraMap_eq_valAdd (R := R) p hp W.c₄ hc4ne hne2]
      exact hc4
    have := valuation_eq_one_of_vAdd_eq_zero (R := R)
      (Units.mk0 ((W.baseChange Lv).c₄) hc4ne') hv
    simpa using this

def hasMultiplicativeReduction_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(c₄) = 0 と v_p(Δ) > 0 なら完備化で乗法還元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_baseChange_of_c4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(c₄) = 0 なら極小性は完備化に移る。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
