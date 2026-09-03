/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim

/-!
# RamifiedBadPrime —— `[GenEll] Definition 3.3` の分

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

/-- ★★★★★★★★★★★★**`e · v_p(j(E)) = −v(q)`**——分岐版（第 1371）。

☆`e = 1` に取れば第 932（`jExp_eq_neg_vAdd_of_j_tateCurveAt`）に戻る。 -/
theorem mul_jExp_eq_neg_vAdd_of_j_tateCurveAt [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic] [(E.baseChange Lv).IsElliptic]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic]
    (hc4 : algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0)
    (hev : algebraMap R Lv (evalAdic tateJinvSeries q hq) ≠ 0)
    (hqne : algebraMap R Lv q ≠ 0)
    (hjne : E.j ≠ 0)
    (hj : (E.baseChange Lv).j = ((tateCurveAt q hq).map (algebraMap R Lv)).j) :
    (e : ℤ) * jExp p E = - vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap R Lv q) hqne) := by
  have hΔbc : (E.baseChange Lv).Δ = algebraMap L Lv E.Δ := WeierstrassCurve.map_Δ _ _
  have hcbc : (E.baseChange Lv).c₄ = algebraMap L Lv E.c₄ := WeierstrassCurve.map_c₄ _ _
  have hjmap : (E.baseChange Lv).j = algebraMap L Lv E.j := by
    rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, ABC3.Found.GenEll.j_eq_inv_Delta_mul,
      hΔbc, hcbc, map_mul, map_inv₀, map_pow]
  have hjT := j_tateCurveAt_inv (R := R) (Lv := Lv) q hq hc4
  have hjE : algebraMap L Lv E.j
      = (algebraMap R Lv (evalAdic tateJinvSeries q hq))⁻¹ := by
    rw [← hjmap, hj, hjT]
  have hmapne : algebraMap L Lv E.j ≠ 0 := by
    rw [hjE]
    exact inv_ne_zero hev
  have hval := vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe E.j hjne hmapne
  rw [jExp, dif_neg hjne, ← hval]
  have hu : Units.mk0 (algebraMap L Lv E.j) hmapne
      = (Units.mk0 (algebraMap R Lv (evalAdic tateJinvSeries q hq)) hev)⁻¹ := by
    refine Units.ext ?_
    rw [Units.val_inv_eq_inv_val, Units.val_mk0, Units.val_mk0]
    exact hjE
  rw [hu, vAdd_inv,
    vAdd_evalAdic_tateJinvSeries (R := R) (K := Lv) q hq (by
      intro h0
      exact hqne (by rw [h0, map_zero])) hev hqne]

/-- ★★★★★★★★★★★★★★★★**Tate 母数の付値は `−e·jExp`**——分岐版（第 1371）。 -/
theorem vAdd_tateParam_eq_neg_mul_jExp [CharZero Lv]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E : WeierstrassCurve L) [E.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [WeierstrassCurve.IsMinimal R (E.baseChange Lv)]
    (h : WeierstrassCurve.HasSplitMultiplicativeReduction R (E.baseChange Lv))
    (hj : jExp p E < 0) :
    vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h) (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q = - ((e : ℤ) * jExp p E) := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  obtain ⟨hq', C₀, hne, hCE⟩ := tateParamR_spec (E.baseChange Lv) h
  have hbase := tateModel_baseChange (E.baseChange Lv) h hCE
  haveI : ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)).IsElliptic := by rw [hbase]; infer_instance
  have hjne : E.j ≠ 0 := by
    intro hc; rw [jExp, dif_pos hc] at hj; omega
  have hc4E : E.c₄ ≠ 0 := by
    intro hc; apply hjne; rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, hc]; ring
  have hqne : algebraMap R Lv (tateParamR (E.baseChange Lv) h) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 hq0
  obtain ⟨u, hu, hueq⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) (tateParamR (E.baseChange Lv) h) hq
  have hev : algebraMap R Lv
      (evalAdic tateJinvSeries (tateParamR (E.baseChange Lv) h) hq) ≠ 0 := by
    rw [hueq, map_mul]
    exact mul_ne_zero hqne ((hu.map (algebraMap R Lv)).ne_zero)
  have hc4Lv : (E.baseChange Lv).c₄ ≠ 0 := by
    show (E.map (algebraMap L Lv)).c₄ ≠ 0
    rw [WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ (algebraMap L Lv).injective).2 hc4E
  have hc4T : algebraMap R Lv
      (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ ≠ 0 := by
    have hmap : algebraMap R Lv (tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄
        = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).c₄ :=
      (WeierstrassCurve.map_c₄ _ _).symm
    rw [hmap, hbase, WeierstrassCurve.variableChange_c₄]
    exact mul_ne_zero (by simp) hc4Lv
  have hjeq : (E.baseChange Lv).j
      = ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv)).j := by
    rw [ABC3.Found.GenEll.j_congr_curve hbase, WeierstrassCurve.variableChange_j]
  have hkey := mul_jExp_eq_neg_vAdd_of_j_tateCurveAt p hpe E
    (tateParamR (E.baseChange Lv) h) hq hc4T hev hqne hjne hjeq
  rw [hkey, neg_neg]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def mul_jExp_eq_neg_vAdd_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(e·v_p(j) = −v(q)。分岐版。★無条件)",
    sectionId := "genell-def-3-3" }

def vAdd_tateParam_eq_neg_mul_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数の付値は −e·jExp。分岐版。★無条件)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
