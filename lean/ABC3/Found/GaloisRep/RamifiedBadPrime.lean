/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1371 ブロック —— **`hcop` の出どころを分岐した拡大に伸ばす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★これは何か——第 932／第 978 の `e` 倍版

第 1370 で `hp` の分岐版（`hpe`）を作った。★本ブロックはそれを

* 第 932 `jExp_eq_neg_vAdd_of_j_tateCurveAt`（`v_p(j) = −v(q)`）
* 第 978 `vAdd_tateParam_eq_neg_jExp`（Tate 母数の付値は `−jExp`）
* `not_dvd_vAdd_tateParam_of_not_dvd_jExp`（`hcop` の出どころ）

に流し込む。☆結論は `e` 倍になる:

    e · jExp p E = − v(q)

★★★そして **`l ∤ e`** なら `¬ l ∣ jExp p E` から `¬ l ∣ v(q)` が出る。
☆`l ∤ e` は第 1369（`e ≤ [L_v′ : L_v]`）から出る——
`L_p(ζ_l)/L_p` は `≤ l−1` 次、分裂用の 2 次拡大は `≤ 2 < 5 ≤ l`。
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

/-! ## ★出典の紐付け(`.src`) -/

def mul_jExp_eq_neg_vAdd_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(e·v_p(j) = −v(q)。分岐版。★無条件)",
    sectionId := "genell-def-3-3" }

def vAdd_tateParam_eq_neg_mul_jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 母数の付値は −e·jExp。分岐版。★無条件)",
    sectionId := "genell-def-3-3" }

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
