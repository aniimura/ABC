/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.TateParamJ
import ABC3.Found.GaloisRep.HtFinJ

/-!
# DegInfTateParam —— `[GenEll] Lemma 3.2` の分

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

/-! ## ★★★★★★★★★★★★`v_p(j) < 0` なら Tate 母数が作れる -/

/-- ★★Tate 曲線の底変換が楼円である十分条件。 -/
theorem tateCurveAt_map_isElliptic [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    (hne : algebraMap R Lv (evalAdic tateJinvSeries q hq) ≠ 0)
    (hc4 : algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic := by
  refine ⟨?_⟩
  have hkey := evalAdic_tateJinvSeries_mul_c4 (I := IsLocalRing.maximalIdeal R) q hq
  have hΔ : ((tateCurveAt q hq).map (algebraMap R Lv)).Δ
      = algebraMap R Lv (evalAdic tateJinvSeries q hq)
        * (algebraMap R Lv (tateCurveAt q hq).c₄) ^ 3 := by
    rw [WeierstrassCurve.map_Δ, ← hkey, map_mul, map_pow]
  rw [hΔ, isUnit_iff_ne_zero]
  exact mul_ne_zero hne (pow_ne_zero 3 hc4)

/-- ★★★★★★★★★★★★**`1/j` が `𝔪` の元なら Tate 母数が作れる**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★★★**2026-09-01（第 934）**——**分裂性を問わない**。
`v_p(j) < 0` なら `1/j ∈ 𝔪` なので、形式逆関数定理（第 933）が
`j(E_q) = j(W)` なる `q` を与える。

☆これと第 932（`jExp = −v(q)`、分裂性不要）を並べると、
非分裂乗法還元の場合も**捧りを経由せずに**扱える。 -/
theorem exists_tateParam_of_inv_j [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (W : WeierstrassCurve Lv) [W.IsElliptic] (hjne : W.j ≠ 0)
    (t : R) (ht : t ∈ IsLocalRing.maximalIdeal R)
    (hjt : algebraMap R Lv t = (W.j)⁻¹)
    (hc4 : ∀ (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R),
      algebraMap R Lv (tateCurveAt q hq).c₄ ≠ 0) :
    ∃ (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
      (_ : ((tateCurveAt q hq).map (algebraMap R Lv)).IsElliptic),
      W.j = ((tateCurveAt q hq).map (algebraMap R Lv)).j := by
  have hc0 : PowerSeries.coeff 0 tateJinvSeries = 0 := by
    rw [PowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact constantCoeff_tateJinvSeries
  obtain ⟨q, hq, hqt⟩ := evalAdic_surjective_of_coeff_one
    (I := IsLocalRing.maximalIdeal R) tateJinvSeries hc0
    coeff_one_tateJinvSeries ht
  have hne : algebraMap R Lv (evalAdic tateJinvSeries q hq) ≠ 0 := by
    rw [hqt, hjt]
    exact inv_ne_zero hjne
  haveI hell := tateCurveAt_map_isElliptic (R := R) (Lv := Lv) q hq hne (hc4 q hq)
  refine ⟨q, hq, hell, ?_⟩
  rw [j_tateCurveAt_inv (R := R) (Lv := Lv) q hq (hc4 q hq), hqt, hjt, inv_inv]

/-! ## ★出典の紐付け(`.src`) -/

def exists_tateParam_of_inv_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(1/j が 𝔪 の元なら Tate 母数が作れる。★無条件)",
    sectionId := "genell-lemma-3-2" }

def jExp_eq_neg_vAdd_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v_p(j) = −v(q)——分裂性を要求しない橋。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateParamK_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の q_{E′} = q_E^l から Kˣ の側へ。★無条件)",
    sectionId := "genell-lemma-3-2" }

def minDeltaExp_eq_mul_of_tateParamR.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(R の側の母数で述べた Δ_min の関係。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
