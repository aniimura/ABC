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

end ABC3.Skeleton.GenEll
