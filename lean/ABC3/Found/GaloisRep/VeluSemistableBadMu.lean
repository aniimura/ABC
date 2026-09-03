/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluSemistableBadDeep
import ABC3.Found.GaloisRep.SemistableScaled
import ABC3.Found.GaloisRep.VAddScaled
import ABC3.Found.GaloisRep.JVeluTateMuK
import ABC3.Meta.Claim

/-!
# 第 1427 ブロック —— **`μ_l` 型の核でも `p ∣ l` で半安定**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——最後の場合

第 1424 で残ったのは **`p ∣ l` かつ核が `μ_l` 型**の場合である。
☆そこでは `1 − ζ^i` が単元でないので Vélu の `v`・`w` は `R` に入らないが、
**`K` の水準の道具はすべて `hlu` なしで揃っている**（第 1129-1138）:

* `veluQuotientFull_tate_mu_K`（第 1135）——商 `= veluCurve (E_q ⊗ K) v w`
* `c4_c6_veluCurve_tate_field`——`c₄ = l⁴·c₄(E_{q^l})`、`c₆ = l⁶·c₆(E_{q^l})`

★したがって第 1426 で `valAdd p (c₄ E′) = 4·valAdd p (l)`、
`valAdd p (c₆ E′) = 6·valAdd p (l)` が読め、
第 1425（短 Weierstrass 形）が **`p ∤ 6` だけを条件に**半安定性を与える。

☆`p ∤ 6` は `p ∣ l` の場合 `l ≥ 5` を意味する。★`l = 3` は別扱いである。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta

open scoped Classical

/-! ## ☆道具 -/

section Tools

variable {L : Type} [Field L] {Lv : Type} [Field Lv] [Algebra L Lv]

/-- ☆`Units.mk0` の付け替え。 -/
theorem vAdd_mk0_congr {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Algebra R Lv] [IsFractionRing R Lv] {a b : Lv} (ha : a ≠ 0) (hb : b ≠ 0) (hab : a = b) :
    vAdd (tateDvrVal R Lv) (Units.mk0 a ha) = vAdd (tateDvrVal R Lv) (Units.mk0 b hb) := by
  subst hab; rfl

/-- ☆変数変換で結ばれていれば `c₄ ≠ 0` は移る。 -/
theorem c4_ne_zero_of_variableChange_eq (X V : WeierstrassCurve Lv)
    (C : WeierstrassCurve.VariableChange Lv) (hEq : C • X = V) (hV : V.c₄ ≠ 0) : X.c₄ ≠ 0 := by
  intro h0
  refine hV ?_
  rw [← hEq, WeierstrassCurve.variableChange_c₄, h0, mul_zero]

/-- ☆変数変換で結ばれていれば `c₆ ≠ 0` は移る。 -/
theorem c6_ne_zero_of_variableChange_eq (X V : WeierstrassCurve Lv)
    (C : WeierstrassCurve.VariableChange Lv) (hEq : C • X = V) (hV : V.c₆ ≠ 0) : X.c₆ ≠ 0 := by
  intro h0
  refine hV ?_
  rw [← hEq, WeierstrassCurve.variableChange_c₆, h0, mul_zero]

/-- ☆底変換で `c₄ ≠ 0` は降りる。 -/
theorem c4_ne_zero_of_baseChange (E' : WeierstrassCurve L)
    (h : (E'.baseChange Lv).c₄ ≠ 0) : E'.c₄ ≠ 0 := by
  intro h0
  refine h ?_
  show (E'.map (algebraMap L Lv)).c₄ = 0
  rw [WeierstrassCurve.map_c₄, h0, map_zero]

/-- ☆底変換で `c₆ ≠ 0` は降りる。 -/
theorem c6_ne_zero_of_baseChange (E' : WeierstrassCurve L)
    (h : (E'.baseChange Lv).c₆ ≠ 0) : E'.c₆ ≠ 0 := by
  intro h0
  refine h ?_
  show (E'.map (algebraMap L Lv)).c₆ = 0
  rw [WeierstrassCurve.map_c₆, h0, map_zero]

end Tools

section Descend

variable {L : Type} [Field L] [NumberField L] {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★★★★`c₄` の付値を `p` の側で読む（`n` 倍版）。 -/
theorem valAdd_c4_scaled_of_baseChange {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) (hc4ne : E'.c₄ ≠ 0) {n : L} (hn : n ≠ 0)
    (hc4ne' : (E'.baseChange Lv).c₄ ≠ 0) (hn' : algebraMap L Lv n ≠ 0)
    (hloc : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).c₄) hc4ne')
      = 4 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv n) hn')) :
    valAdd p (Units.mk0 E'.c₄ hc4ne) = 4 * valAdd p (Units.mk0 n hn) := by
  have hmap : (E'.baseChange Lv).c₄ = algebraMap L Lv E'.c₄ := WeierstrassCurve.map_c₄ _ _
  have hne2 : algebraMap L Lv E'.c₄ ≠ 0 := by rw [← hmap]; exact hc4ne'
  have hloc' : vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv E'.c₄) hne2)
      = 4 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv n) hn') := by
    rw [← hloc]
    exact (vAdd_mk0_congr (R := R) hc4ne' hne2 hmap).symm
  exact valAdd_scaled_of_vAdd he p hpe E'.c₄ n hc4ne hn 4 hne2 hn' (by exact_mod_cast hloc')

/-- ★★★★`c₆` の付値を `p` の側で読む（`n` 倍版）。 -/
theorem valAdd_c6_scaled_of_baseChange {e : ℕ} (he : 1 ≤ e) (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = ((HeightOneSpectrum.valuation L p) x) ^ e)
    (E' : WeierstrassCurve L) (hc6ne : E'.c₆ ≠ 0) {n : L} (hn : n ≠ 0)
    (hc6ne' : (E'.baseChange Lv).c₆ ≠ 0) (hn' : algebraMap L Lv n ≠ 0)
    (hloc : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).c₆) hc6ne')
      = 6 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv n) hn')) :
    valAdd p (Units.mk0 E'.c₆ hc6ne) = 6 * valAdd p (Units.mk0 n hn) := by
  have hmap : (E'.baseChange Lv).c₆ = algebraMap L Lv E'.c₆ := WeierstrassCurve.map_c₆ _ _
  have hne2 : algebraMap L Lv E'.c₆ ≠ 0 := by rw [← hmap]; exact hc6ne'
  have hloc' : vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv E'.c₆) hne2)
      = 6 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv n) hn') := by
    rw [← hloc]
    exact (vAdd_mk0_congr (R := R) hc6ne' hne2 hmap).symm
  exact valAdd_scaled_of_vAdd he p hpe E'.c₆ n hc6ne hn 6 hne2 hn' (by exact_mod_cast hloc')

end Descend

/-! ## ★★★★★★★★★★★★★★★★★★★★本体 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `μ_l` 型の核でも Vélu の商は半安定**——★**`p ∣ l` を許す**（第 1427）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1424 で残った最後の場合である。`p ∣ l` でも次のように通る:

1. `K` の水準の恒等式（第 1129-1138、`hlu` なし）で
   `c₄(veluCurve) = l⁴·c₄(E_{q^l})`・`c₆(veluCurve) = l⁶·c₆(E_{q^l})`
2. 第 1426 で `valAdd p (c₄ E′) = 4·valAdd p (l)`・`valAdd p (c₆ E′) = 6·valAdd p (l)`
3. 第 1425（短 Weierstrass 形）で半安定性

★条件は **`p ∤ 6`**（`h48`・`h864`）だけである。 -/
theorem semistableAt_veluQuotient_bad_mu {L : Type} [Field L] [NumberField L]
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
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) (h2K : (2 : Lv) ≠ 0)
    (h48 : valAdd p (Units.mk0 (48 : L) (by norm_num)) = 0)
    (h864 : valAdd p (Units.mk0 (864 : L) (by norm_num)) = 0)
    (C₀ : WeierstrassCurve.VariableChange R)
    (P : ((tateCurveAt (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv)).toAffine.Point)
    (hcurve : veluQuotientFull
        ((tateCurveAt (tateParamR (E.baseChange Lv) h)
          (tateParamR_mem (E.baseChange Lv) h)).map (algebraMap R Lv))
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • P)))
      = (C₀.map (algebraMap R Lv)) • (E'.baseChange Lv))
    (hu0 : vAdd (tateDvrVal R Lv) ((C₀.map (algebraMap R Lv)).u) = 0)
    (ζ : R) (uζ : Lvˣ) (hζ : IsPrimitiveRoot ζ l)
    (hζu : algebraMap R Lv ζ = (uζ : Lv)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hPz : P = tatePhi (mkTateSetup (tateParamR (E.baseChange Lv) h)
      (tateParamR_mem (E.baseChange Lv) h)
      (tateParamR_ne_zero (E.baseChange Lv) h))
      (tateModel_map_Delta_ne_zero (E.baseChange Lv) h) (QuotientGroup.mk uζ)) :
    SemistableAt p E' := by
  have hq := tateParamR_mem (E.baseChange Lv) h
  have hq0 := tateParamR_ne_zero (E.baseChange Lv) h
  have hΔ := tateModel_map_Delta_ne_zero (E.baseChange Lv) h
  have hinj : Function.Injective (algebraMap R Lv) := IsFractionRing.injective R Lv
  have hql : (tateParamR (E.baseChange Lv) h) ^ l ∈ IsLocalRing.maximalIdeal R :=
    pow_mem_of_mem_ideal hq hl.pos
  have hqlne : algebraMap R Lv ((tateParamR (E.baseChange Lv) h) ^ l) ≠ 0 :=
    (map_ne_zero_iff _ hinj).2 (pow_ne_zero l hq0)
  have hc4T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ ≠ 0 :=
    ((tateCurveAt_c4_isUnit _ hql).map (algebraMap R Lv)).ne_zero
  have hc6T' : algebraMap R Lv
      (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₆ ≠ 0 :=
    ((tateCurveAt_c6_isUnit _ hql).map (algebraMap R Lv)).ne_zero
  obtain ⟨u', hu'u, hueq'⟩ := evalAdic_tateJinvSeries_eq_mul_unit
    (I := IsLocalRing.maximalIdeal R) ((tateParamR (E.baseChange Lv) h) ^ l) hql
  have hev' : algebraMap R Lv
      (evalAdic tateJinvSeries ((tateParamR (E.baseChange Lv) h) ^ l) hql) ≠ 0 := by
    rw [hueq', map_mul]
    exact mul_ne_zero hqlne ((hu'u.map (algebraMap R Lv)).ne_zero)
  haveI hellQl : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).IsElliptic :=
    tateCurveAt_map_isElliptic _ hql hev' hc4T'
  have hΔl : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).Δ ≠ 0 :=
    ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).isUnit_Δ.ne_zero
  -- ★段 1: `K` の水準の `v`・`w` と恒等式
  obtain ⟨v, w, hv, hw, hell⟩ := exists_vw_tate_mu_K hl hζ
    (tateParamR (E.baseChange Lv) h) hq hql hΔl
  have hne : ∀ i ∈ (range l).erase 0, algebraMap R Lv (1 - ζ ^ i) ≠ 0 :=
    fun i hi => zeta_pow_sub_ne_zero_K hinj hζ hi
  have hquot := veluQuotientFull_tate_mu_K
    (mkTateSetup (tateParamR (E.baseChange Lv) h) hq hq0) hΔ
    (dvrTatePhiAddEquiv (tateParamR (E.baseChange Lv) h) hq hq0 hΔ) (fun _ => rfl)
    hl.pos ζ uζ hζu hζl hord hne v w h2K hv hw
  have hvS : (l : Lv) ^ 6 * v
      = algebraMap R Lv (veluSV l ζ (tateParamR (E.baseChange Lv) h) hq) := by
    rw [hv]; exact sum_natCast_pow_mul_veluV2_K hinj hζ _ _ hq
  have hwS : (l : Lv) ^ 8 * (2 * w)
      = algebraMap R Lv (veluSW l ζ (tateParamR (E.baseChange Lv) h) hq) := by
    rw [hw]; exact sum_natCast_pow_mul_veluW_K hinj hζ _ _ hq
  obtain ⟨h4, h6⟩ := c4_c6_veluCurve_tate_field hl hζ
    (tateParamR (E.baseChange Lv) h) hq hql v w hvS hwS
  rw [hPz] at hcurve
  have hEq := hcurve.symm.trans hquot
  -- ★段 2: `l ≠ 0`
  have hlLv : ((l : Lv)) ≠ 0 := by
    simpa using (Nat.cast_ne_zero (R := Lv)).2 hl.ne_zero
  have hlL : ((l : L)) ≠ 0 := Nat.cast_ne_zero.2 hl.ne_zero
  have hlmap : ((l : Lv)) = algebraMap L Lv ((l : L)) := by rw [map_natCast]
  have hlLv' : algebraMap L Lv ((l : L)) ≠ 0 := by rw [← hlmap]; exact hlLv
  -- ★段 3: `c₄`・`c₆` の付値（`Lv` の側）
  have hzc4 : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₄ = algebraMap R Lv
        (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₄ :=
    WeierstrassCurve.map_c₄ _ _
  have hzc4ne : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₄ ≠ 0 := by rw [hzc4]; exact hc4T'
  have hz4 : vAdd (tateDvrVal R Lv) (Units.mk0 (((tateCurveAt
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).c₄) hzc4ne) = 0 := by
    rw [vAdd_mk0_congr (R := R) hzc4ne hc4T' hzc4]
    exact tateDvrVal_eq_zero_of_isUnit _ (tateCurveAt_c4_isUnit _ hql) hc4T'
  have hzc6 : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₆ = algebraMap R Lv
        (tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).c₆ :=
    WeierstrassCurve.map_c₆ _ _
  have hzc6ne : ((tateCurveAt ((tateParamR (E.baseChange Lv) h) ^ l) hql).map
      (algebraMap R Lv)).c₆ ≠ 0 := by rw [hzc6]; exact hc6T'
  have hz6 : vAdd (tateDvrVal R Lv) (Units.mk0 (((tateCurveAt
      ((tateParamR (E.baseChange Lv) h) ^ l) hql).map (algebraMap R Lv)).c₆) hzc6ne) = 0 := by
    rw [vAdd_mk0_congr (R := R) hzc6ne hc6T' hzc6]
    exact tateDvrVal_eq_zero_of_isUnit _ (tateCurveAt_c6_isUnit _ hql) hc6T'
  have hVc4ne : (veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)) v w).c₄ ≠ 0 := by
    rw [h4]
    exact mul_ne_zero (pow_ne_zero 4 hlLv) hzc4ne
  have hVc6ne : (veluCurve ((tateCurveAt (tateParamR (E.baseChange Lv) h) hq).map
      (algebraMap R Lv)) v w).c₆ ≠ 0 := by
    rw [h6]
    exact mul_ne_zero (pow_ne_zero 6 hlLv) hzc6ne
  have hXc4 : (E'.baseChange Lv).c₄ ≠ 0 :=
    c4_ne_zero_of_variableChange_eq _ _ _ hEq hVc4ne
  have hXc6 : (E'.baseChange Lv).c₆ ≠ 0 :=
    c6_ne_zero_of_variableChange_eq _ _ _ hEq hVc6ne
  have hloc4 := vAdd_c4_scaled_of_eq (R := R) (E'.baseChange Lv) hXc4
    (C₀.map (algebraMap R Lv)) _ hEq hu0 hlLv hzc4ne hz4 h4
  have hloc6 := vAdd_c6_scaled_of_eq (R := R) (E'.baseChange Lv) hXc6
    (C₀.map (algebraMap R Lv)) _ hEq hu0 hlLv hzc6ne hz6 h6
  -- ★段 4: `p` の側へ降ろす
  have hc4ne : E'.c₄ ≠ 0 := c4_ne_zero_of_baseChange (Lv := Lv) E' hXc4
  have hc6ne : E'.c₆ ≠ 0 := c6_ne_zero_of_baseChange (Lv := Lv) E' hXc6
  have hloc4' : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).c₄) hXc4)
      = 4 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv ((l : L))) hlLv') := by
    rw [hloc4, vAdd_mk0_congr (R := R) hlLv hlLv' hlmap]
  have hloc6' : vAdd (tateDvrVal R Lv) (Units.mk0 ((E'.baseChange Lv).c₆) hXc6)
      = 6 * vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv ((l : L))) hlLv') := by
    rw [hloc6, vAdd_mk0_congr (R := R) hlLv hlLv' hlmap]
  have hv4 := valAdd_c4_scaled_of_baseChange (R := R) he p hpe E' hc4ne hlL hXc4 hlLv' hloc4'
  have hv6 := valAdd_c6_scaled_of_baseChange (R := R) he p hpe E' hc6ne hlL hXc6 hlLv' hloc6'
  -- ★段 5: 第 1425
  exact semistableAt_of_valAdd_c4_c6_scaled p E' hlL hc4ne hc6ne hv4 hv6 h48 h864

/-! ## ★出典の紐付け(`.src`) -/

def semistableAt_veluQuotient_bad_mu.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l 型の核でも Vélu の商は半安定。★p ∣ l を許す)",
    sectionId := "genell-lemma-3-5" }

def semistableAt_veluQuotient_bad_mu.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_tate_mu_K(第 1135、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_tate_mu_K") 1,
    .citation "[ABC3]" "c4_c6_veluCurve_tate_field(第 1136 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.c4_c6_veluCurve_tate_field") 1,
    .citation "[ABC3]" "semistableAt_of_valAdd_c4_c6_scaled(第 1425、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_of_valAdd_c4_c6_scaled") 1,
    .implicitStep
      ("★★★★★**2026-09-02（第 1427）**——第 1424 で残った最後の場合" ++
       "（`p ∣ l` かつ核が `μ_l` 型）を閉じた。" ++
       "☆`K` の水準の恒等式（`hlu` なし）で `c₄ = l⁴·c₄(E_{q^l})` が読め、" ++
       "第 1426 で `p` の側の付値に降ろし、第 1425（短 Weierstrass 形）で半安定性が出る。" ++
       "★条件は `p ∤ 6` だけ——`p ∣ l` の場合は `l ≥ 5` を意味する。" ++
       "☆`l = 3` は別扱いである。") 17 ]

end ABC3.Found.GaloisRep
