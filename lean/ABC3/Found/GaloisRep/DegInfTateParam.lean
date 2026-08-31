/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.HtFinJ

/-!
# Galois (G6) 第 892 ブロック —— **★★★★★★★★`R` の側の母数で述べた `Δ_min` の関係**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_tateParam`（第 719）は仮説を **`Kˣ` の側**の母数
`tateParamK` で述べているが、葉 1 の連鎖（第 883・891）が与えるのは
**`R` の側**の `tateParamR` である。

★本ブロックはその隔たりを埋める——`tateParamK` は `tateParamR` の像なので
`Units.ext` と `map_pow` だけで済む。

| 定理 | 内容 |
|---|---|
| `tateParamK_pow` | ★`q_{E′} = q_E^l`（`R`）⇒ 同じ（`K`） |
| `minDeltaExp_eq_mul_of_tateParamR` | ★★★★★★★★**`R` の側で述べた `Δ_min` の関係** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★**`R` の側の関係は `Kˣ` の側の関係を与える**。 -/
theorem tateParamK_pow [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (W W' : WeierstrassCurve Lv) [W.IsElliptic] [W.IsMinimal R]
    [W'.IsElliptic] [W'.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (h' : W'.HasSplitMultiplicativeReduction R)
    (l : ℕ) (hq : tateParamR W' h' = (tateParamR W h) ^ l) :
    tateParamK W' h' = (tateParamK W h) ^ l := by
  refine Units.ext ?_
  rw [Units.val_pow_eq_pow_val, tateParamK, tateParamK, Units.val_mk0, Units.val_mk0,
    ← map_pow, hq]

/-- ★★★★★★★★**`R` の側の母数で述べた `Δ_min` の関係**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★これが葉 1（`tateParam_quot_velu`、第 891）の結論を
`Lemma 3.5` の入力 `hloc` に直接繋ぐ形である。 -/
theorem minDeltaExp_eq_mul_of_tateParamR [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E E' : WeierstrassCurve L) (l : ℕ)
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic] [(E'.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (h' : (E'.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (C' : VariableChange L) (hC' : IsMinimal (primeSubring p) (C' • E'))
    (hc4ne' : (C' • E').c₄ ≠ 0) (hc4' : valAdd p (Units.mk0 ((C' • E').c₄) hc4ne') = 0)
    (hq : tateParamR (E'.baseChange Lv) h' = (tateParamR (E.baseChange Lv) h) ^ l) :
    minDeltaExp p E' = l * minDeltaExp p E :=
  minDeltaExp_eq_mul_of_tateParam E E' l h h' p hp C hC hc4ne hc4 C' hC' hc4ne' hc4'
    (tateParamK_pow (R := R) _ _ h h' l hq)

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

/-! ## ★★★★★★★★★★`v_p(j) = −v(q)`——分裂性を要求しない橋 -/

/-- ★★**Tate 曲線の `j` の逆は `1/j` の評価そのもの**。 -/
theorem j_tateCurveAt_inv [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic]
    (hc4 : algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R Lv)).j
      = (algebraMap R Lv (evalAdic tateJinvSeries q hq))⁻¹ := by
  have hkey := evalAdic_tateJinvSeries_mul_c4 (I := IsLocalRing.maximalIdeal R) q hq
  have hΔ : ((tateCurveAt q hq).map (algebraMap R Lv)).Δ
      = algebraMap R Lv (evalAdic tateJinvSeries q hq)
        * (algebraMap R Lv (tateCurveAt q hq).c₄) ^ 3 := by
    rw [WeierstrassCurve.map_Δ, ← hkey, map_mul, map_pow]
  rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, hΔ, WeierstrassCurve.map_c₄]
  field_simp

/-- ★★★★★★★★★★**`v_p(j(E)) = −v(q)`**——`E` の `j` が Tate 曲線 `E_q` の `j` なら。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 932）**——この橋は**分裂性を要求しない**。
`j` の一致だけを受けて `jExp` を出すので、非分裂の場合にもそのまま使える。

☆道は 3 段:

1. `j(E ⊗ Lv) = φ(j(E))`（`j = Δ⁻¹c₄³` なので体準同型と可換）
2. `j(E_q ⊗ Lv) = (φ(1/j の評価))⁻¹`（第 932 前半）
3. `v(1/j の評価) = v(q)`（第 931） -/
theorem jExp_eq_neg_vAdd_of_j_tateCurveAt [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (E : WeierstrassCurve L) [E.IsElliptic] [(E.baseChange Lv).IsElliptic]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic]
    (hc4 : algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0)
    (hev : algebraMap R Lv (evalAdic tateJinvSeries q hq) ≠ 0)
    (hqne : algebraMap R Lv q ≠ 0)
    (hjne : E.j ≠ 0)
    (hj : (E.baseChange Lv).j = ((tateCurveAt q hq).map (algebraMap R Lv)).j) :
    jExp p E = - vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap R Lv q) hqne) := by
  -- ★段 1
  have hΔbc : (E.baseChange Lv).Δ = algebraMap L Lv E.Δ := WeierstrassCurve.map_Δ _ _
  have hcbc : (E.baseChange Lv).c₄ = algebraMap L Lv E.c₄ := WeierstrassCurve.map_c₄ _ _
  have hjmap : (E.baseChange Lv).j = algebraMap L Lv E.j := by
    rw [ABC3.Found.GenEll.j_eq_inv_Delta_mul, ABC3.Found.GenEll.j_eq_inv_Delta_mul,
      hΔbc, hcbc, map_mul, map_inv₀, map_pow]
  -- ★段 2
  have hjT := j_tateCurveAt_inv (R := R) (Lv := Lv) q hq hc4
  have hjE : algebraMap L Lv E.j
      = (algebraMap R Lv (evalAdic tateJinvSeries q hq))⁻¹ := by
    rw [← hjmap, hj, hjT]
  -- ★段 3
  have hmapne : algebraMap L Lv E.j ≠ 0 := by
    rw [hjE]
    exact inv_ne_zero hev
  have hval := vAdd_algebraMap_eq_valAdd (R := R) p hp E.j hjne hmapne
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

/-! ## ★出典の紐付け(`.src`) -/

def jExp_eq_neg_vAdd_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v_p(j) = −v(q)——分裂性を要求しない橋。★無条件)",
    sectionId := "genell-lemma-3-2" }

def hasMultiplicativeReduction_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(c₄) = 0 と v_p(Δ) > 0 なら完備化で乗法還元。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_baseChange_of_c4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(c₄) = 0 なら極小性は完備化に移る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def tateParamK_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の q_{E′} = q_E^l から Kˣ の側へ。★無条件)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_tateParamR.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の母数で述べた Δ_min の関係。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
