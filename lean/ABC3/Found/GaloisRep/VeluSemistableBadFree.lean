/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluSemistableBadRam
import ABC3.Found.GaloisRep.TateDeepInvolution
import ABC3.Found.GenEll.MuPrimitiveRootOrDeep
import ABC3.Meta.Claim

/-!
# 第 1416 ブロック —— **★★★★★★★★★★★★★★★★★★★★悪い素点の半安定性から
`hcop` が落ちた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か

第 1404 の `semistableAt_veluQuotient_bad_ram` は
`hcop`（`l ∤ v(q_E)`）を仮定していた——核が `μ_l` 型であることを強制するためである。

★第 1408 で分かったとおり `hcop` は**底変換で保たれない**
（`jExp P (E⊗M) = e(P∣p)·jExp p E` なので `l ∣ e(P∣p)` で壊れる）。
☆したがって `VeluQuotOK`（任意の有限拡大の素点で半安定性を要求する）には使えない。

★★★本ブロックは `hcop` を**落とす**。道は第 1410-1415 の二者択一である:

| 核 | `c₄(veluCurve) = c₄(E_q) + 240v` が単元である理由 |
|---|---|
| `μ_l` 型 | ☆`c₄ + 240v = l⁴·c₄(E_{q^l})`（第 1388 の `h4`）——`l` が単元だから |
| 深い核 | ★**`v ∈ 𝔪`**（第 1412 の `veluV_deep_mem`）——核の座標がすべて深いから |

☆どちらの場合も `semistableAt_velu_of_veluCurve_eq_ram`（第 1404）に流し込める。

★★★これで悪い素点の半安定性は**局所的には無条件**になった
（残るのは `p ∤ l`——`hlu`——だけである）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定（分岐版）**——★**`hcop` なし**（第 1416）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1404 の `semistableAt_veluQuotient_bad_ram` から `hcop`（`l ∤ v(q_E)`）を落とした形。

★★★道は二者択一（第 1414 の `primitiveRoot_or_deep_of_torsion_point`）:

* `μ_l` 型なら第 1404 と同じ（`c₄ + 240v = l⁴·c₄(E_{q^l})`）
* 深い核なら `v ∈ 𝔪` から直接（第 1412 の `isUnit_c4_add_240_deep`）

☆残る局所仮定は `l` が `R` で単元（`p ∤ l`）と `l` が奇素数だけである。 -/
theorem semistableAt_veluQuotient_bad_ram_free {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E']
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hΔL : E'.Δ ≠ 0) :
    SemistableAt p E' := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  have hu0 := vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  rcases ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl P hP hP0 with
    ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ | ⟨y, hPz, hdeep⟩
  · -- ★★`μ_l` 型——第 1404 と同じ道
    have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
      pow_mem_of_mem_ideal hq hl.pos
    have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
      (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 (pow_ne_zero l hq0)
    have hc4T' : algebraMap R Lv
        (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
      ((tateCurveAt_c4_isUnit
        ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).ne_zero
    obtain ⟨u', hu'u, hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
      (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
    have hev' : algebraMap R Lv
        (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
      rw [hueq', map_mul]
      exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
    haveI : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
        (algebraMap R Lv)).IsElliptic :=
      tateCurveAt_map_isElliptic _ hql hev' hc4T'
    obtain ⟨v, w, hv, hw, hell, h4, h6⟩ :=
      ABC3.Skeleton.GenEll.exists_vw_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hodd
        hlu hql h2 ζ hζ
    have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
    have hquot := veluQuotientFull_tate_mu
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
    rw [hPz] at hcurve
    have hunit : IsUnit ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ + 240 * v) := by
      rw [h4]
      exact (hlu.pow 4).mul
        (tateCurveAt_c4_isUnit ((tateParamR (E.baseChange Lv) h) ^ l) hql)
    exact semistableAt_velu_of_veluCurve_eq_ram he p hpe E' hΔL
      (tateParamR (E.baseChange Lv) h) hq v w hunit C₀ (hcurve.symm.trans hquot) hu0
  · -- ★★★深い核——`v ∈ 𝔪` から `c₄` が単元
    obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
    have hyl : (2 * m + 1) • tatePhi
        (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
        (QuotientGroup.mk y) = 0 := by rw [← hPz]; exact hP
    obtain ⟨w, hw⟩ := exists_veluW_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      m y hyl hdeep
    have hquot := veluQuotientFull_tate_deep
      (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
      (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
      y hdeep _ w h2K rfl hw
    rw [hPz] at hcurve
    exact semistableAt_velu_of_veluCurve_eq_ram he p hpe E' hΔL
      (tateParamR (E.baseChange Lv) h) hq _ w
      (isUnit_c4_add_240_deep (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0)
        y hdeep _ rfl)
      C₀ (hcurve.symm.trans hquot) hu0

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuotient_bad_ram_free.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——分岐版・★hcop なし)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_ram_free.needs : List ProofObligation :=
  [ .citation "[ABC3]" "primitiveRoot_or_deep_of_torsion_point(第 1414、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.primitiveRoot_or_deep_of_torsion_point") 1,
    .citation "[ABC3]" "veluQuotientFull_tate_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_deep") 1,
    .citation "[ABC3]" "exists_veluW_deep(第 1415、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_veluW_deep") 1,
    .citation "[ABC3]" "isUnit_c4_add_240_deep(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_add_240_deep") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1416）**——第 1404 から `hcop`（`l ∤ v(q_E)`）が落ちた。" ++
       "☆核が `μ_l` 型でも深くても `c₄(veluCurve) = c₄(E_q) + 240v` は単元である" ++
       "（前者は `l⁴·c₄(E_{q^l})`、後者は `v ∈ 𝔪`）。" ++
       "★★★これは `hcopDoesNotDescend2026_09_02`（底変換で `l ∤ jExp` が壊れる）を" ++
       "**根本から解消する**——悪い素点の半安定性は局所的には無条件になった。") 17 ]

end ABC3.Found.GaloisRep
