/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Found.GaloisRep.VeluSemistableBad
import ABC3.Found.GaloisRep.RamifiedValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1404 ブロック —— **悪い素点の半安定性（分岐版）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`hp` を `^e` に弱める

第 1327（`semistableAt_velu_of_veluCurve_eq`）と第 1388
（`semistableAt_veluQuotient_bad`）は**不分岐**な局所体
（`v_{Lv}(x) = v_p(x)`）を仮定していた。

★しかし実際に手に入る局所体は `L_p(ζ_l)` であり、一般に**分岐**する
（`exists_locCyc_package`、第 1377 は `v_{Lv}(x) = v_p(x)^e` の形で出す）。

☆本ブロックは仮定を `^e` に弱める。★中身は 1 行である——
`vAdd_algebraMap_eq_valAdd` を `vAdd_algebraMap_eq_mul_valAdd`（第 1370、在庫）に替えると
結論が `e·valAdd = 0` になり、`e ≥ 1` から `valAdd = 0` が出る。

★★★α の葉（第 1370-1384）で同じ形の一般化をしたのと**同型の作業**である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

/-- ☆Tate モデルの Vélu の商と変数変換で結ばれていれば `c₄ ≠ 0`（第 1408）。

★`c₄(veluCurve) = c₄(E_q) + 240v` が単元だからである。 -/
theorem c4_ne_zero_of_veluCurve_eq {Lv : Type} [Field Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (v w : R)
    (hunit : IsUnit ((tateCurveAt q hq).c₄ + 240 * v))
    (C₀ : WeierstrassCurve.VariableChange R) (X : WeierstrassCurve Lv)
    (hEq : (C₀.map (algebraMap R Lv)) • X
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv)) :
    X.c₄ ≠ 0 := by
  intro h0
  have h1 := congrArg WeierstrassCurve.c₄ hEq
  rw [WeierstrassCurve.variableChange_c₄, WeierstrassCurve.map_c₄, veluCurve_c₄, h0,
    mul_zero] at h1
  refine hunit.ne_zero ?_
  refine (IsFractionRing.injective R Lv) ?_
  rw [map_zero]
  exact h1.symm

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Vélu の商は悪い素点で半安定（分岐版）**——★**無条件**（第 1404、第 1408 で `hc4` を落とした）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1327 の `hp` を `v_{Lv}(x) = v_p(x)^e`（`e ≥ 1`）に弱めた形である。
★★`c₄ ≠ 0` は `hEq` と `hunit` から出るので仮定しない（第 1408）。 -/
theorem semistableAt_velu_of_veluCurve_eq_ram {L : Type} [Field L] [NumberField L]
    {Lv : Type} [Field Lv] [Algebra L Lv]
    {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E']
    (hΔ : E'.Δ ≠ 0)
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (v w : R)
    (hunit : IsUnit ((tateCurveAt q hq).c₄ + 240 * v))
    (C₀ : WeierstrassCurve.VariableChange R)
    (hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R Lv))
    (hu : vAdd (tateDvrVal R Lv) ((C₀.map (algebraMap R Lv)).u) = 0) :
    SemistableAt p E' := by
  have hinj : Function.Injective (algebraMap L Lv) := (algebraMap L Lv).injective
  have hc4' : (E'.baseChange Lv).c₄ ≠ 0 :=
    c4_ne_zero_of_veluCurve_eq q hq v w hunit C₀ (E'.baseChange Lv) hEq
  have hc4 : E'.c₄ ≠ 0 := by
    intro h0
    refine hc4' ?_
    rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_c₄, h0, map_zero]
  have hloc := vAdd_c4_of_veluCurve_eq (R := R) (K := Lv) q hq v w hunit C₀
    (E'.baseChange Lv) hc4' hEq hu
  have hne2 : algebraMap L Lv E'.c₄ ≠ 0 := (map_ne_zero_iff _ hinj).2 hc4
  have hEqU : (Units.mk0 ((E'.baseChange Lv).c₄) hc4')
      = Units.mk0 (algebraMap L Lv E'.c₄) hne2 := by
    refine Units.ext ?_
    show (E'.baseChange Lv).c₄ = algebraMap L Lv E'.c₄
    rw [WeierstrassCurve.baseChange, WeierstrassCurve.map_c₄]
  rw [hEqU, vAdd_algebraMap_eq_mul_valAdd (R := R) p hpe E'.c₄ hc4 hne2] at hloc
  have hzero : valAdd p (Units.mk0 E'.c₄ hc4) = 0 := by
    rcases mul_eq_zero.1 hloc with h | h
    · exact absurd h (by positivity)
    · exact h
  exact semistableAt_of_c4_valAdd_zero p E' hΔ hc4 hzero

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点では Vélu の商は半安定（分岐版）**——★**無条件**（第 1404）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1388 の `hp` を `v_{Lv}(x) = v_p(x)^e`（`e ≥ 1`）に弱めた形である。

★★★これで `L_p(ζ_l)`（一般に分岐する）の上の Tate モデルが使える。 -/
theorem semistableAt_veluQuotient_bad_ram {L : Type} [Field L] [NumberField L]
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
    (hΔL : E'.Δ ≠ 0) :
    SemistableAt p E' := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    pow_mem_of_mem_ideal hq hl.pos
  obtain ⟨C₀, P, hP, hP0, hcurve, hCE⟩ :=
    exists_variableChange_veluQuotient_tateModel E E' h hl hQ h2K hE'
  obtain ⟨ζ, uζ, hζ, hζu, hζl, hord, hPz⟩ :=
    ABC3.Found.GenEll.exists_primitiveRoot_of_torsion_point
      (tateParamR (E.baseChange Lv) h) hq hq0 hΔ hl hcop P hP hP0
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
  have hquot' : veluQuotientFull
      ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map (algebraMap R Lv))
      (((range l).erase 0).image (fun k : ℕ => pointCoords
        (k • tatePhi (mkTateSetup (K := Lv) (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
          (QuotientGroup.mk uζ))))
      = (veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := hquot
  have hEq : (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv)
      = (veluCurve
        (tateCurveAt (tateParamR (E.baseChange Lv) h) hq) v w).map (algebraMap R Lv) := by
    rw [← hcurve, hPz]
    exact hquot'
  have hu0 := vAdd_tateModel_u_eq_zero (E.baseChange Lv) h C₀ hCE hq0
  have hunit : IsUnit ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).c₄ + 240 * v) := by
    rw [h4]
    exact (hlu.pow 4).mul
      (tateCurveAt_c4_isUnit ((tateParamR (E.baseChange Lv) h) ^ l) hql)
  exact semistableAt_velu_of_veluCurve_eq_ram he p hpe E' hΔL
    (tateParamR (E.baseChange Lv) h) hq v w hunit C₀ hEq hu0

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_velu_of_veluCurve_eq_ram.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商は悪い素点で半安定——分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_ram.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(悪い素点では Vélu の商は半安定——分岐版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_ram.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_veluQuotient_bad(第 1388、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.semistableAt_veluQuotient_bad") 1,
    .citation "[ABC3]" "vAdd_algebraMap_eq_mul_valAdd(第 1370、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_algebraMap_eq_mul_valAdd") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1404）**——実際に手に入る局所体は `L_p(ζ_l)` であり" ++
       "一般に**分岐**するので、第 1327・第 1388 の不分岐の仮定を `^e` に弱めた。" ++
       "☆中身は 1 行——`vAdd_algebraMap_eq_valAdd` を第 1370 の `e` 倍版に替えると" ++
       "結論が `e·valAdd = 0` になり、`e ≥ 1` から `valAdd = 0` が出る。" ++
       "★α の葉（第 1370-1384）と同型の作業である。") 17 ]

end ABC3.Found.GaloisRep
