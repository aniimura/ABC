/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Found.GaloisRep.JVeluTateMuK
import ABC3.Meta.Claim

/-!
# 第 1137 ブロック —— **悪い素点での `Δ_min` の関係の `hlu` なし版**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か

`minDeltaExp_eq_mul_of_globalVelu'`（第 999）は `hlu : IsUnit ((l : R))`（＝`p ∤ l`）と
「Vélu の係数 `v`・`w` が `R` に取れる」を要求する。
★本ブロックは**どちらも要求しない版**である。

☆`v`・`w` は商体 `Lv` に取り、座標は `tateXK`・`tateYK`（分母を払った形）で受ける。

| 段 | 使うもの | 第 |
|---|---|---|
| 商の楕円性 | `isElliptic_veluQuotientFull_tate_mu_K` | 1135 |
| `j` の一致 | `j_veluQuot_eq_j_tate_pow_K` | 1136 |
| `φ(1 − ζ^i) ≠ 0` | `zeta_pow_sub_ne_zero_K` | 1134 |

★`p ∣ l` では Vélu のモデルは極小から `l¹²` 離れるが、`j` はその差を見ない。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Found.GenEll WeierstrassCurve IsDedekindDomain
  NumberField Finset
open scoped Classical

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 悪い素点での `Δ_min` の関係——
`hlu` なし版**（第 1137）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★**`IsUnit ((l : R))` を仮説に置いていない**——`p ∣ l` の悪い素点でも成り立つ。 -/
theorem minDeltaExp_eq_mul_of_globalVelu'_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsMinimal (p.adicCompletionIntegers L)]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x) = (HeightOneSpectrum.valuation L p) x)
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hΔ : ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).toAffine.Δ ≠ 0)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).v
      (mkTateSetup (K := p.adicCompletion L)
        (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).Q))
    (hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    (h2 : (2 : (p.adicCompletionIntegers L)) ≠ 0)
    (h2K : (2 : (p.adicCompletion L)) ≠ 0)
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletion L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (tateXK (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYK (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
                (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
              (tateXK (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYK (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
                      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
                    (tateXK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ (veluCurve ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))) v w).IsElliptic)
    [((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic]
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hinj : Function.Injective (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)) :=
    IsFractionRing.injective _ _
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel p E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ hl hcop P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ := hvw ζ hζ
  haveI := hell
  have hne : ∀ i ∈ (range l).erase 0,
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)) (1 - ζ ^ i) ≠ 0 :=
    fun i hi => ABC3.Found.GaloisRep.zeta_pow_sub_ne_zero_K hinj hζ hi
  -- ★`μ_l` の形での楕円性（第 1135、`hlu` なし）
  haveI hellMu := ABC3.Found.GaloisRep.isElliptic_veluQuotientFull_tate_mu_K
      (mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ)
      (fun _ => rfl) hl.pos ζ uζ hζu hζl hord hne v w h2K hv hw hell
  -- ☆`P` の形での楕円性
  haveI hellP : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic := by
    rw [hPz]; exact hellMu
  -- ★曲線の水準で `P` と `μ_l` を繋ぐ
  have hcurveEq : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
        (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
        (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := p.adicCompletion L)
          (tateParamR (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)
          (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ
          (QuotientGroup.mk uζ)))) := by
    rw [hPz]
    rfl
  -- ☆第 1136 で `E_{q^l}` に繋ぐ（`hlu` なし）
  haveI hE1 : (veluCurve ((tateCurveAt (mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).q
      (mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).hq).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))) v w).IsElliptic := hell
  haveI hE3 : (veluQuotientFull ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ (QuotientGroup.mk uζ))))).IsElliptic := hellMu
  haveI hE2 : ((tateCurveAt ((mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)).q ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic := by assumption
  have hjtate := ABC3.Found.GaloisRep.j_veluQuot_eq_j_tate_pow_K
    (mkTateSetup (K := p.adicCompletion L) (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h)) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange (p.adicCompletion L)) h) (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h) hΔ)
    (fun _ => rfl) hl hodd hζ uζ hζu hζl hord hql h2K v w hv hw
  exact minDeltaExp_eq_mul_of_j_tate_pow p hp E E' h hssE hssE' hjneg hql
    ((hj hellP).trans ((ABC3.Found.GenEll.j_congr_curve hcurveEq).trans hjtate))



set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**各悪い素点での `Δ_min` の関係（第 1138、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (p.adicCompletion L)).IsElliptic]
    [(E'.baseChange (p.adicCompletion L)).IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      (E.baseChange (p.adicCompletion L)))
    (h : (E.baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletion L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
            (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (tateXK (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            (tateYK (ζ ^ i)
              ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
                (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
              (tateXK (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
              (tateYK (ζ ^ i)
                ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR (E.baseChange (p.adicCompletion L)) h)
                (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map
                      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
                    (tateXK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                    (tateYK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h))
                  * tateXK (ζ ^ i)
                      ((tateParamR (E.baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR (E.baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)))
        ∧ (veluCurve ((tateCurveAt (tateParamR (E.baseChange (p.adicCompletion L)) h)
              (tateParamR_mem (E.baseChange (p.adicCompletion L)) h)).map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))) v w).IsElliptic)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hp := ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion L p
  have hql : (tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal (tateParamR_mem (E.baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性——`c₄` は定数項 `1`、`1/j` は `q^l · 単元`
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero (E.baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  exact minDeltaExp_eq_mul_of_globalVelu'_K p E E' h hp hssE hssE' hjneg hl hodd
    (tateModel_map_Delta_ne_zero (E.baseChange (p.adicCompletion L)) h)
    (not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hjneg hcop)
    hql (two_ne_zero_adicCompletionIntegers L p) (two_ne_zero_adicCompletion L p)
    hvw hQ hE'

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**極小化した対に当てる形（第 1138、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_vc_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    (hvw : ∀ ζ : (p.adicCompletionIntegers L), IsPrimitiveRoot ζ l →
      ∃ v w : (p.adicCompletion L),
      v = ∑ i ∈ (range l).erase 0,
          veluV2 ((tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
            (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)).map
              (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
            (tateXK (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            (tateYK (ζ ^ i)
              ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
              (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU ((tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)).map
                (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
              (tateXK (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
              (tateYK (ζ ^ i)
                ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
            + 2 * (veluV2 ((tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                    (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)).map
                      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)))
                    (tateXK (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                    (tateYK (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h))
                  * tateXK (ζ ^ i)
                      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) * (ζ ^ i) ^ (l - 1))
                      (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
                      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)))
        ∧ (veluCurve ((tateCurveAt (tateParamR ((C • E).baseChange (p.adicCompletion L)) h)
              (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h)).map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))) v w).IsElliptic)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  have h2L : (2 : L) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq := ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq C E E' hQ h2L hE'
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hjneg
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [jExp_variableChange p E C]; exact hcop
  have hkey := minDeltaExp_eq_mul_at_bad_prime_K p (C • E) (C • E')
    (semistableAt_variableChange p E C hssE) (semistableAt_variableChange p E' C hssE')
    hjC hl hodd hcopC hmin h hvw hQ' hEq
  rwa [minDeltaExp_variableChange p E' C, minDeltaExp_variableChange p E C] at hkey

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★★★★★**局所データを自前で作る形（第 1138、★`hlu` なし）**。

★★★★**これが節点 5 の本体である**——`p ∣ l` の悪い素点でも `Δ_min` は `l` 倍になる。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_full_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (E E' : WeierstrassCurve L)
    [E.IsElliptic] [E'.IsElliptic]
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (C : WeierstrassCurve.VariableChange L)
    [((C • E).baseChange (p.adicCompletion L)).IsElliptic]
    [((C • E').baseChange (p.adicCompletion L)).IsElliptic]
    (hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers L)
      ((C • E).baseChange (p.adicCompletion L)))
    (h : ((C • E).baseChange (p.adicCompletion L)).HasSplitMultiplicativeReduction
      (p.adicCompletionIntegers L))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI := hmin
  have hql : (tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l
      ∈ IsLocalRing.maximalIdeal (p.adicCompletionIntegers L) :=
    pow_mem_of_mem_ideal
      (tateParamR_mem ((C • E).baseChange (p.adicCompletion L)) h) hl.pos
  -- ★`E_{q^l} ⊗ Lv` の楕円性（第 1000 と同じ作り方）
  have hqlne : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2
      (pow_ne_zero l (tateParamR_ne_zero ((C • E).baseChange (p.adicCompletion L)) h))
  have hc4T' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
        hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal (p.adicCompletionIntegers L))
    ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql
  have hev' : algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L)
      (evalAdic tateJinvSeries
        ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne
      ((hu'.map (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).ne_zero)
  haveI hEllQL : ((tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l)
      hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  have hΔl : ((tateCurveAt ((tateParamR ((C • E).baseChange (p.adicCompletion L)) h) ^ l) hql).map
      (algebraMap (p.adicCompletionIntegers L) (p.adicCompletion L))).Δ ≠ 0 :=
    ((WeierstrassCurve.isElliptic_iff _).1 hEllQL).ne_zero
  exact minDeltaExp_eq_mul_at_bad_prime_vc_K p E E' hssE hssE' hjneg hl hodd hcop C hmin h
    (fun ζ hζ => ABC3.Found.GaloisRep.exists_vw_tate_mu_K hl hodd hζ _ _ hql hΔl)
    hQ hE'


set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**一般の局所体で（第 1139、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_gen_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  have hinj : Function.Injective (algebraMap R Lv) := IsFractionRing.injective R Lv
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    pow_mem_of_mem_ideal hq hl.pos
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective _ _)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu', hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'.map (algebraMap R Lv)).ne_zero)
  haveI hEllQL : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  have hΔl : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).Δ ≠ 0 := ((WeierstrassCurve.isElliptic_iff _).1 hEllQL).ne_zero
  have hcop' := not_dvd_vAdd_tateParam_of_not_dvd_jExp p hp E h hjneg hcop
  obtain ⟨P, hP, hP0, hj⟩ := exists_point_j_tateModel_gen E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    exists_primitiveRoot_of_torsion_point (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop'
      P hP hP0
  obtain ⟨v, w, hv, hw, hell⟩ :=
    ABC3.Found.GaloisRep.exists_vw_tate_mu_K (K := Lv) hl hodd hζ
      (tateParamR (E.baseChange Lv) h) hq hql hΔl
  have hne : ∀ i ∈ (range l).erase 0, algebraMap R Lv (1 - ζ ^ i) ≠ 0 :=
    fun i hi => ABC3.Found.GaloisRep.zeta_pow_sub_ne_zero_K hinj hζ hi
  haveI hellMu := ABC3.Found.GaloisRep.isElliptic_veluQuotientFull_tate_mu_K
    (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ)
    (fun _ => rfl) hl.pos ζ uζ hζu hζl hord hne v w h2K hv hw hell
  haveI hellP : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))).IsElliptic := by
    rw [hPz]; exact hellMu
  have hcurveEq : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ)))) := by
    rw [hPz]
    rfl
  haveI hE1 : (veluCurve ((tateCurveAt
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0).q
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0).hq).map
      (algebraMap R Lv)) v w).IsElliptic := hell
  haveI hE2 : ((tateCurveAt
      ((mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0).q ^ l) hql).map
      (algebraMap R Lv)).IsElliptic := hEllQL
  haveI hE3 : (veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ))))).IsElliptic := hellMu
  have hjtate := ABC3.Found.GaloisRep.j_veluQuot_eq_j_tate_pow_K
    (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ)
    (fun _ => rfl) hl hodd hζ uζ hζu hζl hord hql h2K v w hv hw
  exact minDeltaExp_eq_mul_of_j_tate_pow p hp E E' h hssE hssE' hjneg hql
    ((hj hellP).trans ((ABC3.Found.GenEll.j_congr_curve hcurveEq).trans hjtate))

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**不分岐 2 次拡大に載せ替えた形（第 1139、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_ext_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    {f : Polynomial R} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange (FractionRing (AdjoinRoot f))).IsElliptic]
    [(E.baseChange (FractionRing (AdjoinRoot f))).IsMinimal (AdjoinRoot f)]
    [(E'.baseChange (FractionRing (AdjoinRoot f))).IsElliptic]
    (h' : (E.baseChange (FractionRing (AdjoinRoot f))).HasSplitMultiplicativeReduction
      (AdjoinRoot f))
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hloc : IsLocalRing (AdjoinRoot f) := inferInstance
  -- ★完備性——第 1012（`R`-加群）→ 第 1013（イデアルの側）→ 第 1022（極大の同定）
  haveI hcomp0 : IsAdicComplete (IsLocalRing.maximalIdeal R) (AdjoinRoot f) :=
    isAdicComplete_adjoinRoot (IsLocalRing.maximalIdeal R) hf hdeg
  haveI hcomp : IsAdicComplete (IsLocalRing.maximalIdeal (AdjoinRoot f)) (AdjoinRoot f) := by
    rw [maximalIdeal_adjoinRoot_eq_map hf hdeg hirr]
    exact isAdicComplete_map_algebraMap _
  -- ☆標数
  haveI hchar : CharZero (AdjoinRoot f) :=
    charZero_of_injective_algebraMap (algebraMap_adjoinRoot_injective hf hdeg)
  haveI hcharK : CharZero (FractionRing (AdjoinRoot f)) :=
    IsFractionRing.charZero (AdjoinRoot f)
  -- ★付値の橋
  have hp' : ∀ x : L, (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
      (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f)))
      (algebraMap L (FractionRing (AdjoinRoot f)) x)
      = (HeightOneSpectrum.valuation L p) x := by
    intro x
    rw [halg, valuation_quadFieldHom hf hdeg hirr, hp]
  exact minDeltaExp_eq_mul_at_bad_prime_gen_K p hp' E E' h' hssE hssE' hjneg hl hodd hcop
    (NeZero.ne' (2 : FractionRing (AdjoinRoot f))).symm hQ hE'

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★★★**非分裂の枝（第 1139、★`hlu` なし）**。 -/
theorem minDeltaExp_eq_mul_at_bad_prime_nonsplit_K {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L))
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    [((C • E).baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsMinimal R ((C • E).baseChange Lv)]
    (hA : IsUnit (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)).c₄)
    (hirr : Irreducible
      ((splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).map
        (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    (hssE : SemistableAt p E) (hssE' : SemistableAt p E') (hjneg : jExp p E < 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hcop : ¬ ((l : ℤ) ∣ jExp p E))
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) :
    minDeltaExp p E' = l * minDeltaExp p E := by
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  have hfm : (splitQuadPoly
      (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).Monic :=
    monic_splitQuadPoly _ hA
  have hfd : (splitQuadPoly
      (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA).natDegree = 2 :=
    natDegree_splitQuadPoly _ hA
  haveI hdom := isDomain_adjoinRoot hfm hirr
  haveI hlocR := isLocalRing_adjoinRoot hfm hfd hirr
  haveI hdvr := isDiscreteValuationRing_adjoinRoot hfm hfd hirr
  letI : Algebra L (FractionRing (AdjoinRoot
      (splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA))) :=
    ((quadFieldHom (K := Lv) hfm hfd).comp (algebraMap L Lv)).toAlgebra
  have halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot
      (splitQuadPoly (WeierstrassCurve.integralModel R ((C • E).baseChange Lv)) hA))) x
      = quadFieldHom hfm hfd (algebraMap L Lv x) := fun _ => rfl
  haveI hmin' := isMinimal_baseChange_ext p hp hfm hfd hirr halg E C hC hc4ne hc4
  haveI hmult' := hasMultiplicativeReduction_ext p hp hfm hfd hirr halg E C hC hc4ne hc4 hjneg
  have h' := hasSplitMultiplicativeReduction_ext (C • E) hA rfl halg
  -- ★点と Vélu の商を `C • E` に運ぶ（第 1002 と同じ）
  have h2L : (2 : L) ≠ 0 := two_ne_zero
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq := ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq C E E' hQ h2L hE'
  have hjC : jExp p (C • E) < 0 := by rw [jExp_variableChange p E C]; exact hjneg
  have hcopC : ¬ ((l : ℤ) ∣ jExp p (C • E)) := by
    rw [jExp_variableChange p E C]; exact hcop
  have hkey := minDeltaExp_eq_mul_at_bad_prime_ext_K p hp hfm hfd hirr halg (C • E) (C • E') h'
    (semistableAt_variableChange p E C hssE) (semistableAt_variableChange p E' C hssE')
    hjC hl hodd hcopC hQ' hEq
  rwa [minDeltaExp_variableChange p E' C, minDeltaExp_variableChange p E C] at hkey

def minDeltaExp_eq_mul_of_globalVelu'_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——大域の Vélu の商で受ける形。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_of_globalVelu'_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "j_veluQuot_eq_j_tate_pow_K(j の一致、第 1136、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.j_veluQuot_eq_j_tate_pow_K") 1,
    .citation "[ABC3]" "isElliptic_veluQuotientFull_tate_mu_K(商の楕円性、第 1135、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isElliptic_veluQuotientFull_tate_mu_K") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_j_tate_pow(j から Δ_min へ、第 995、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_j_tate_pow") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1137）**——第 999 の `hlu` と「`v`・`w` が `R` に取れる」の" ++
       "両方が落ちた。☆`v`・`w` は商体に取り、座標は `tateXK`・`tateYK` で受ける。" ++
       "★これで節点 5 の本体が閉じた。") 8 ]

def minDeltaExp_eq_mul_at_bad_prime_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各悪い素点での Δ_min の関係。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_of_globalVelu'_K(第 1137、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_globalVelu'_K") 1,
    .citation "[ABC3]" "evalAdic_tateJinvSeries_eq_mul_unit(j の級数の単元分解、在庫)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.evalAdic_tateJinvSeries_eq_mul_unit") 1 ]

def minDeltaExp_eq_mul_at_bad_prime_vc_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小化した対に当てる形。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_vc_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_K(第 1138、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_K") 1,
    .citation "[ABC3]" "veluQuotientFull_vcPoint_eq(Vélu の商を変数変換で移す、第 969)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1 ]

def minDeltaExp_eq_mul_at_bad_prime_full_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(局所データを自前で作る形。★p ∣ l の悪い素点でも Δ_min は l 倍)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_full_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_vw_tate_mu_K(hvw を無条件に作る、第 1138、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_vw_tate_mu_K") 1,
    .citation "[ABC3]" "minDeltaExp_eq_mul_of_globalVelu'_K(第 1137、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_of_globalVelu'_K") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1138）**——`p ∤ l` を要求していた連鎖が" ++
       "すべて `hlu` なしで組み直せた。☆`hvw` の楕円性は " ++
       "`Δ(veluCurve) = l¹²·Δ(E_{q^l})`（第 1131）から出る。" ++
       "★これで §3 の枠の節点 5 が閉じた。") 6 ]

def minDeltaExp_eq_mul_at_bad_prime_gen_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点での Δ_min の関係——一般の局所体で。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_gen_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_vw_tate_mu_K(hvw を無条件に作る、第 1138、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_vw_tate_mu_K") 1,
    .citation "[ABC3]" "j_veluQuot_eq_j_tate_pow_K(j の一致、第 1136、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.j_veluQuot_eq_j_tate_pow_K") 1 ]

def minDeltaExp_eq_mul_at_bad_prime_ext_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(不分岐 2 次拡大に載せ替えた形。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_ext_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_gen_K(第 1139、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_gen_K") 1,
    .citation "[ABC3]" "isAdicComplete_adjoinRoot(完備性、第 1012、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isAdicComplete_adjoinRoot") 1 ]

def minDeltaExp_eq_mul_at_bad_prime_nonsplit_K.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(非分裂の枝。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-5" }

def minDeltaExp_eq_mul_at_bad_prime_nonsplit_K.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq_mul_at_bad_prime_ext_K(第 1139、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.minDeltaExp_eq_mul_at_bad_prime_ext_K") 1,
    .citation "[ABC3]" "isMinimal_baseChange_ext(拡大での極小性、第 1029、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_baseChange_ext") 1 ]

end ABC3.Skeleton.GenEll
