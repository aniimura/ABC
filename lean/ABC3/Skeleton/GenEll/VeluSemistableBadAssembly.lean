/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Found.GaloisRep.VeluSemistableBad
import ABC3.Meta.Claim

/-!
# 第 1388 ブロック —— **悪い素点での Vélu の商の半安定性（配管）**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——第 1327 に局所データを渡す

`semistableAt_velu_of_veluCurve_eq`（第 1327、証明済み）は

* `q`・`v`・`w`（Tate 母数と Vélu の和）
* `hunit : IsUnit ((tateCurveAt q hq).c₄ + 240·v)`
* `C₀` と `hEq : (C₀ ⊗ L_v) • (E′ ⊗ L_v) = veluCurve(E_q, v, w) ⊗ L_v`
* `hu : v(C₀.u) = 0`

を受けて `SemistableAt p E′` を出す。★これらを作る配管は
`vAdd_Delta_veluQuotient_tate`（第 1061、証明済み）の証明に**そのまま入っている**。

☆本ブロックはその配管を取り出して第 1327 に流す。
★`hunit` は `h4`（`c₄(veluCurve) = l⁴·c₄(E_{q^l})`）と
`tateCurveAt_c4_isUnit`（在庫）と `hlu`（`l` は単元）から出る。
-/

namespace ABC3.Skeleton.GenEll

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定**——★**無条件**（第 1388）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1061（`vAdd_Delta_veluQuotient_tate`）の配管を第 1327 に流しただけである。

★★★これで `semistableAt_veluQuotientFull` の**悪い素点側が閉じる**。 -/
theorem semistableAt_veluQuotient_bad {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [CharZero Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [CharZero R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = (HeightOneSpectrum.valuation L p) x)
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    [(E.baseChange Lv).IsElliptic] [(E.baseChange Lv).IsMinimal R]
    [(E'.baseChange Lv).IsElliptic]
    [WeierstrassCurve.IsIntegral (primeSubring p) E']
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hcop : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).v
      (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h)
        (tateParamR_mem (E.baseChange Lv) h)
        (tateParamR_ne_zero (E.baseChange Lv) h)).Q))
    (hlu : IsUnit ((l : R))) (h2 : (2 : R) ≠ 0) (h2K : (2 : Lv) ≠ 0)
    {Q : E.toAffine.Point} (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E
      (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))))
    (hΔL : E'.Δ ≠ 0) (hc4L : E'.c₄ ≠ 0) :
    SemistableAt p E' := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := ABC3.Found.GaloisRep.tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    ABC3.Found.GaloisRep.pow_mem_of_mem_ideal hq hl.pos
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop P hP hP0
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective R Lv)).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((ABC3.Found.GaloisRep.tateCurveAt_c4_isUnit
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu'u, hueq'⟩ := ABC3.Found.GaloisRep.evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
  haveI : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    ABC3.Found.GaloisRep.tateCurveAt_map_isElliptic _ hql hev' hc4T'
  obtain ⟨v, w, hv, hw, hell, h4, h6⟩ :=
    exists_vw_tate_mu (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hodd hlu hql h2 ζ hζ
  have hu := ABC3.Found.GenEll.isUnit_one_sub_pow_of_isUnit_natCast hl.pos hζ hlu
  have hquot := ABC3.Found.GaloisRep.veluQuotientFull_tate_mu
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    hl.pos ζ uζ hζu hζl hord hu v w h2K hv hw
  have hquot' : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ))))
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := hquot
  have hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (ABC3.Found.GenEll.veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := by
    rw [← hcurve, hPz]
    exact hquot'
  have hu0 := ABC3.Found.GaloisRep.vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  -- ★`hunit` は `h4` ＋ `tateCurveAt_c4_isUnit` ＋ `hlu` から
  have hunit : IsUnit ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ + 240 * v) := by
    rw [h4]
    exact (hlu.pow 4).mul
      (ABC3.Found.GaloisRep.tateCurveAt_c4_isUnit
        ((tateParamR (E.baseChange Lv) h) ^ l) hql)
  exact ABC3.Found.GaloisRep.semistableAt_velu_of_veluCurve_eq p hp E' hΔL hc4L
    (tateParamR (E.baseChange Lv) h) hq v w hunit C₀ hEq hu0

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def semistableAt_veluQuotient_bad.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_velu_of_veluCurve_eq(第 1327、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_velu_of_veluCurve_eq") 1,
    .citation "[ABC3]" "exists_variableChange_veluQuotient_tateModel(第 1058、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_variableChange_veluQuotient_tateModel") 1,
    .citation "[ABC3]" "veluQuotientFull_tate_mu(第 1057、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1388）**——第 1061 の配管を第 1327 に流しただけである。" ++
       "☆`hunit` は `h4`（`c₄(veluCurve) = l⁴·c₄(E_{q^l})`）から出る。" ++
       "★★★これで `semistableAt_veluQuotientFull` の**悪い素点側が閉じた**。") 17 ]

end ABC3.Skeleton.GenEll
