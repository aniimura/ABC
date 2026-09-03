/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluMuK
import ABC3.Meta.Claim

/-!
# 第 1136 ブロック —— **`j(E_q/μ_l) = j(E_{q^l})` の `hlu` なし版**（`Found`、節点 5）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★これは何か

`j_veluQuot_eq_j_tate_pow`（`Skeleton/GenEll/TateIsogeny.lean:1305`）は
`hlu : IsUnit ((l : R))` を要求する。★本ブロックはそれを**取らない版**を組む。

☆道具はすべて揃っている:

| 段 | 使うもの | 第 |
|---|---|---|
| 商 `= veluCurve (E_q ⊗ K) v w` | `veluQuotientFull_tate_mu_K` | 1135 |
| `v`・`w` は `SV`・`SW` の像を `l⁶`・`2l⁸` で割ったもの | `sum_natCast_pow_mul_veluV2_K` 等 | 1134 |
| `c₄ = l⁴c₄′`・`c₆ = l⁶c₆′` | `c4_velu_tateDF`・`c6_velu_tateDF` | 1129-1130 |
| `c₄`・`c₆` が `l⁴`・`l⁶` 倍なら `j` は等しい | `j_eq_of_c4_c6` | 1132 |

★`p ∣ l` では Vélu のモデルは極小から `l¹²` 離れるが、**`j` はその差を見ない**。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll Finset QuotientGroup
open scoped Classical

section Sums

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]

/-- ★`SV`——`c₄` の分母なし恒等式に現れる和。 -/
noncomputable def veluSV (l : ℕ) (ζ q : R) (hq : q ∈ I) : R :=
  ∑ i ∈ (range l).erase 0,
    veluV2DF l (tateCurveAt q hq)
      (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)

/-- ★`SW`——`c₆` の分母なし恒等式に現れる和。 -/
noncomputable def veluSW (l : ℕ) (ζ q : R) (hq : q ∈ I) : R :=
  ∑ i ∈ (range l).erase 0,
    ((l : R) ^ 2 * veluUDF l (tateCurveAt q hq)
          (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      + 2 * (veluV2DF l (tateCurveAt q hq)
              (tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            * tateXpairDF l (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))

end Sums

/-! ## ★★★★★★★★★★`c₄`・`c₆` の関係（`v`・`w` を与えられた形） -/

section C4C6

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★
**`v`・`w` が `SV`・`SW` の像を割ったものなら `c₄ = l⁴c₄′`・`c₆ = l⁶c₆′`**（第 1136）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★**`IsUnit ((l : R))` を仮説に置いていない**。 -/
theorem c4_c6_veluCurve_tate_field {l : ℕ} (hl : l.Prime) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (v w : K)
    (hv : (l : K) ^ 6 * v = algebraMap R K (veluSV l ζ q hq))
    (hw : (l : K) ^ 8 * (2 * w) = algebraMap R K (veluSW l ζ q hq)) :
    (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).c₄
        = (l : K) ^ 4 * ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).c₄
      ∧ (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).c₆
        = (l : K) ^ 6 * ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).c₆ := by
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hlK : (l : K) ≠ 0 := by
    intro h
    have hR : ((l : ℕ) : R) = 0 := hinj (by rw [map_natCast, h, map_zero])
    exact (Nat.cast_ne_zero.2 hl.pos.ne' : ((l : ℕ) : R) ≠ 0) hR
  have h4 : (l : R) ^ 6 * (tateCurveAt q hq).c₄ + 240 * veluSV l ζ q hq
      = (l : R) ^ 10 * (tateCurveAt (q ^ l) hql).c₄ := c4_velu_tateDF hl hζ q hq hql
  have h6 : (l : R) ^ 8 * (tateCurveAt q hq).c₆ + 504 * (l : R) ^ 2 * veluSV l ζ q hq
        + 3024 * veluSW l ζ q hq
      = (l : R) ^ 14 * (tateCurveAt (q ^ l) hql).c₆ := by
    have h := c6_velu_tateDF hl hζ q hq hql
    rw [c6DFlhs] at h
    exact h
  have h4K := congrArg (algebraMap R K) h4
  have h6K := congrArg (algebraMap R K) h6
  simp only [map_add, map_mul, map_pow, map_natCast, map_ofNat] at h4K h6K
  constructor
  · rw [veluCurve_c₄, WeierstrassCurve.map_c₄, WeierstrassCurve.map_c₄]
    refine mul_left_cancel₀ (pow_ne_zero 6 hlK) ?_
    linear_combination h4K + 240 * hv
  · rw [veluCurve_c₆, WeierstrassCurve.map_c₆, WeierstrassCurve.map_c₆,
      WeierstrassCurve.map_b₂, tateCurveAt_b₂, map_one]
    refine mul_left_cancel₀ (pow_ne_zero 8 hlK) ?_
    linear_combination h6K + 504 * (l : K) ^ 2 * hv + 3024 * hw

end C4C6

/-! ## ★★★★★★★★★★★★★★★★`j` の一致（`hlu` なし） -/

section J

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `j(E_q/μ_l) = j(E_{q^l})`**（第 1136、★`IsUnit ((l : R))` なし）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★これが節点 5 の `j_veluQuot_eq_j_tate_pow` の `hlu` なし版である。
☆`p ∣ l` では Vélu のモデルは極小から `l¹²` 離れるが、**`j` はその差を見ない**。 -/
theorem j_veluQuot_eq_j_tate_pow_K (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (uζ : Kˣ) (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hql : S.q ^ l ∈ I) (h2K : (2 : K) ≠ 0)
    (v w : K)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU ((tateCurveAt S.q S.hq).map (algebraMap R K))
          (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
          (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        + 2 * (veluV2 ((tateCurveAt S.q S.hq).map (algebraMap R K))
                (tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
                (tateYK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
              * tateXK (I := I) (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)))
    [(veluCurve ((tateCurveAt S.q S.hq).map (algebraMap R K)) v w).IsElliptic]
    [((tateCurveAt (S.q ^ l) hql).map (algebraMap R K)).IsElliptic]
    [(veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords
            (k • tatePhi S hΔ (QuotientGroup.mk uζ))))).IsElliptic] :
    (veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk uζ))))).j
      = ((tateCurveAt (S.q ^ l) hql).map (algebraMap R K)).j := by
  haveI : CharZero K := IsFractionRing.charZero_of_isFractionRing R
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hlK : (l : K) ≠ 0 := by
    intro h
    have hR : ((l : ℕ) : R) = 0 := hinj (by rw [map_natCast, h, map_zero])
    exact (Nat.cast_ne_zero.2 hl.pos.ne' : ((l : ℕ) : R) ≠ 0) hR
  have hne : ∀ i ∈ (range l).erase 0, algebraMap R K (1 - ζ ^ i) ≠ 0 :=
    fun i hi => zeta_pow_sub_ne_zero_K hinj hζ hi
  -- ★商を `veluCurve` の形に（第 1135、`IsUnit` なし）
  have hquot := veluQuotientFull_tate_mu_K S hΔ Φ hΦ hl.pos ζ uζ hζu hζl hord hne v w h2K hv hw
  -- ★`v`・`w` を `SV`・`SW` の像で書く（第 1134）
  have hvS : (l : K) ^ 6 * v = algebraMap R K (veluSV l ζ S.q S.hq) := by
    rw [hv]
    exact sum_natCast_pow_mul_veluV2_K hinj hζ (tateCurveAt S.q S.hq) S.q S.hq
  have hwS : (l : K) ^ 8 * (2 * w) = algebraMap R K (veluSW l ζ S.q S.hq) := by
    rw [hw]
    exact sum_natCast_pow_mul_veluW_K hinj hζ (tateCurveAt S.q S.hq) S.q S.hq
  obtain ⟨h4, h6⟩ := c4_c6_veluCurve_tate_field hl hζ S.q S.hq hql v w hvS hwS
  rw [ABC3.Found.GenEll.j_congr_curve hquot]
  exact j_eq_of_c4_c6 _ _ (l : K) hlK h4 h6

end J

/-! ## ★★★★★★★★★★★★★★`hvw` を無条件に作る -/

section ExistsVW

variable {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `μ_l` に対する Vélu の係数は商体に取れて、その商は楕円である**（第 1138）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★★**`IsUnit ((l : R))` を仮説に置いていない**——`exists_vw_tate_mu`（第 1003）の
`hlu` なし版である。
☆楕円性は `Δ(veluCurve) = l¹²·Δ(E_{q^l})`（第 1131 の `Delta_of_c4_c6`）から出る。 -/
theorem exists_vw_tate_mu_K {l : ℕ} (hl : l.Prime) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I)
    (hΔl : ((tateCurveAt (q ^ l) hql).map (algebraMap R K)).Δ ≠ 0) :
    ∃ v w : K,
      v = ∑ i ∈ (range l).erase 0,
          veluV2 ((tateCurveAt q hq).map (algebraMap R K))
            (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        ∧ 2 * w = ∑ i ∈ (range l).erase 0,
          (veluU ((tateCurveAt q hq).map (algebraMap R K))
              (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 ((tateCurveAt q hq).map (algebraMap R K))
                    (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        ∧ (veluCurve ((tateCurveAt q hq).map (algebraMap R K)) v w).IsElliptic := by
  haveI : CharZero K := IsFractionRing.charZero_of_isFractionRing R
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  have hlK : (l : K) ≠ 0 := by
    intro h
    have hR : ((l : ℕ) : R) = 0 := hinj (by rw [map_natCast, h, map_zero])
    exact (Nat.cast_ne_zero.2 hl.pos.ne' : ((l : ℕ) : R) ≠ 0) hR
  have h2K : (2 : K) ≠ 0 := two_ne_zero
  set SVK : K := ∑ i ∈ (range l).erase 0,
      veluV2 ((tateCurveAt q hq).map (algebraMap R K))
        (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) with hSVK
  set SWK : K := ∑ i ∈ (range l).erase 0,
      (veluU ((tateCurveAt q hq).map (algebraMap R K))
          (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 ((tateCurveAt q hq).map (algebraMap R K))
                (tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXK (I := I) (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)) with hSWK
  refine ⟨SVK, SWK / 2, rfl, by field_simp, ?_⟩
  have hvS : (l : K) ^ 6 * SVK = algebraMap R K (veluSV l ζ q hq) :=
    sum_natCast_pow_mul_veluV2_K hinj hζ (tateCurveAt q hq) q hq
  have hwS : (l : K) ^ 8 * (2 * (SWK / 2)) = algebraMap R K (veluSW l ζ q hq) := by
    rw [show (2 : K) * (SWK / 2) = SWK from by field_simp]
    exact sum_natCast_pow_mul_veluW_K hinj hζ (tateCurveAt q hq) q hq
  obtain ⟨h4, h6⟩ := c4_c6_veluCurve_tate_field hl hζ q hq hql SVK (SWK / 2) hvS hwS
  have hΔeq := Delta_of_c4_c6 _ _ ((l : K)) h4 h6
  rw [WeierstrassCurve.isElliptic_iff, isUnit_iff_ne_zero, hΔeq]
  exact mul_ne_zero (pow_ne_zero _ hlK) hΔl

end ExistsVW

/-! ## ★出典の紐付け(`.src`) -/

def j_veluQuot_eq_j_tate_pow_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l による Vélu の商の j は E_{q^l} の j。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-2" }

def c4_c6_veluCurve_tate_field.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Vélu の商の c₄・c₆ は l⁴・l⁶ 倍。★IsUnit (l) を取らない)",
    sectionId := "genell-lemma-3-2" }

def j_veluQuot_eq_j_tate_pow_K.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu_K(商の形、第 1135、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu_K") 1,
    .citation "[ABC3]" "sum_natCast_pow_mul_veluV2_K(v は SV の像を l⁶ で割ったもの、第 1134)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.sum_natCast_pow_mul_veluV2_K") 1,
    .citation "[ABC3]" "j_eq_of_c4_c6(c₄・c₆ が l⁴・l⁶ 倍なら j は等しい、第 1132)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.j_eq_of_c4_c6") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1136）**——`j_veluQuot_eq_j_tate_pow`（第 996）の " ++
       "`hlu` なし版が組めた。☆`p ∣ l` では Vélu のモデルは極小から `l¹²` 離れるが、" ++
       "**`j` はその差を見ない**（`j = Δ⁻¹·c₄³` で `l¹²` が約分される）。" ++
       "★残るのは `minDeltaExp_eq_mul_of_globalVelu′`（第 999）を本定理で組み直す段である。") 8 ]

end ABC3.Found.GaloisRep
