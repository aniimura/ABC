import ABC3.Found.GaloisRep.DegInf
import ABC3.Found.GaloisRep.TateCurveWitness

/-!
# Galois (G8) 第 320 ブロック —— **★★★★★★★局所高さは極小判別式の付値である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★到達点

> **`localHeightOf W h = vAdd (Δ_W)`**(`W` が極小のとき)

★★★これが第 319 の `deg∞` と第 310 の局所高さを**繋ぐ**唯一の環である。
`deg∞` は `v_p(Δ_min)` の和として定義したので、この一致がなければ
`degInf_ge_localHeight` は書けない。

## ★★★★★証明の筋——「Tate モデルは極小である」

第 310 の `tateParamR_spec` は、`W` の整モデルを変数変換 `C` で
**Tate 標準形 `tateCurveAt q`** に移せることを言う。そこで

1. `tateCurveAt q` の `c₄` は**単元**(第 304 の `tateCurveAt_c4_isUnit`)
2. ★★★★**`v(c₄) = 0` なら極小**(`isMinimal_of_c4_vAdd_eq_zero`、本ブロック)
   ——`c₄` の付値は変数変換で `-4v(u)`、`Δ` は `-12v(u)` だけ動く。
   `v(c₄) = 0` かつ整なら `v(u) ≤ 0`、一方 `Δ` の極小性は `v(u) ≥ 0` を強いる。
3. したがって `C.map` の `u` の付値は `0`(第 312 の `minimal_u_vAdd_eq`)
4. ★★`Δ(C•W) = u⁻¹²Δ(W)` より `v(Δ_q) = v(Δ_W)`
5. ★★★`Δ_q = q·(単元)`(第 101)より `v(Δ_q) = v(q) = 局所高さ`

★★★★★**完備化も一意性も要らない**——極小性の判定を `c₄` の付値だけで済ませたのが鍵。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `integral_c4_vAdd_nonneg` | ★★整モデルの `c₄` は付値非負 |
| `vAdd_c4_variableChange` | ★★★`c₄` の付値は `-4v(u)` だけ動く |
| `variableChange_c4_ne_zero` | ★`c₄ ≠ 0` は変数変換で保たれる |
| `isMinimal_of_c4_vAdd_eq_zero` | ★★★★★**`v(c₄) = 0` ⟹ 極小** |
| `tateDvrVal_eq_zero_of_isUnit` | ★★単元の付値は `0` |
| `isIntegral_baseChange` | ★係数環から来た曲線は整 |
| `localHeight_eq_vAdd_Delta` | ★★★★★★★**局所高さ = `v(Δ_min)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Interface.GaloisRep

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-! ## ★★★★★`v(c₄) = 0` ならば極小 -/

/-- ★★整モデルの `c₄` は付値が非負である。 -/
theorem integral_c4_vAdd_nonneg (W : WeierstrassCurve K) [IsIntegral R W] (hc4 : W.c₄ ≠ 0) :
    0 ≤ vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4) := by
  refine tateDvrVal_nonneg_of_mem _ ⟨(integralModel R W).c₄, ?_⟩
  exact WeierstrassCurve.integralModel_c₄_eq R W

/-- ★★★`c₄` の付値は変数変換で `-4·v(u)` だけ動く。 -/
theorem vAdd_c4_variableChange (W : WeierstrassCurve K) (hc4 : W.c₄ ≠ 0) (C : VariableChange K)
    (hc4' : (C • W).c₄ ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 ((C • W).c₄) hc4')
      = vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4) - 4 * vAdd (tateDvrVal R K) C.u := by
  have hu : (Units.mk0 ((C • W).c₄) hc4') = C.u⁻¹ ^ 4 * Units.mk0 W.c₄ hc4 := by
    refine Units.ext ?_
    show (C • W).c₄ = _
    rw [WeierstrassCurve.variableChange_c₄]
    push_cast
    simp
  rw [hu, vAdd_mul, vAdd_pow, vAdd_inv]
  omega

/-- ★`c₄ ≠ 0` は変数変換で保たれる。 -/
theorem variableChange_c4_ne_zero (W : WeierstrassCurve K) (hc4 : W.c₄ ≠ 0) (C : VariableChange K) :
    (C • W).c₄ ≠ 0 := by
  rw [WeierstrassCurve.variableChange_c₄]
  exact mul_ne_zero (pow_ne_zero _ (Units.ne_zero _)) hc4

set_option maxHeartbeats 1600000 in
/-- ★★★★★**`c₄` の付値が `0` ならば極小モデルである**。

★★★乗法還元の判定条件そのもの——`v(c₄) = 0` は「`c₄` が単元」であり、
そのとき変数変換の `u` は付値 `0` しか許されないので `v(Δ)` はこれ以上下がらない。 -/
theorem isMinimal_of_c4_vAdd_eq_zero (W : WeierstrassCurve K) [hint : IsIntegral R W]
    (hΔ : W.Δ ≠ 0) (hc4 : W.c₄ ≠ 0)
    (h : vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4) = 0) : IsMinimal R W := by
  rw [WeierstrassCurve.isMinimal_iff]
  constructor
  · show IsIntegral R ((1 : VariableChange K) • W)
    rw [one_smul]
    exact hint
  · intro C hCint hle
    have hCint' : IsIntegral R (C • W) := hCint
    have hCΔ : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    have hCc4 : (C • W).c₄ ≠ 0 := variableChange_c4_ne_zero W hc4 C
    have h0 : WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W)
        ≤ WeierstrassCurve.valuation_Δ_aux R (C • W) := hle
    rw [one_smul] at h0
    rw [WeierstrassCurve.valuation_Δ_aux, dif_pos hint, WeierstrassCurve.valuation_Δ_aux,
      dif_pos hCint'] at h0
    have hle2 : (IsDiscreteValuationRing.maximalIdeal R).valuation K W.Δ
        ≤ (IsDiscreteValuationRing.maximalIdeal R).valuation K ((C • W).Δ) := h0
    have hbridge : vAdd (tateDvrVal R K) (Units.mk0 ((C • W).Δ) hCΔ)
        ≤ vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) := by
      rw [← valuation_le_iff_vAdd_le]
      exact hle2
    have hnn : 0 ≤ vAdd (tateDvrVal R K) (Units.mk0 ((C • W).c₄) hCc4) :=
      integral_c4_vAdd_nonneg (R := R) (C • W) hCc4
    have hchg4 := vAdd_c4_variableChange (R := R) W hc4 C hCc4
    have hchgΔ := vAdd_Delta_variableChange (R := R) W hΔ C
    show WeierstrassCurve.valuation_Δ_aux R (C • W)
      ≤ WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W)
    rw [one_smul, WeierstrassCurve.valuation_Δ_aux, dif_pos hCint',
      WeierstrassCurve.valuation_Δ_aux, dif_pos hint]
    show (IsDiscreteValuationRing.maximalIdeal R).valuation K ((C • W).Δ)
      ≤ (IsDiscreteValuationRing.maximalIdeal R).valuation K W.Δ
    exact (valuation_le_iff_vAdd_le (R := R) (Units.mk0 ((C • W).Δ) hCΔ)
      (Units.mk0 W.Δ hΔ)).2 (by omega)

/-! ## ★★単元と係数環 -/

/-- ★★係数環の単元は付値 `0` を持つ。 -/
theorem tateDvrVal_eq_zero_of_isUnit (r : R) (hr : IsUnit r) (h0 : algebraMap R K r ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K r) h0) = 0 := by
  have h1 : 0 ≤ vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K r) h0) :=
    tateDvrVal_nonneg_of_mem _ ⟨r, rfl⟩
  have h2 : 0 ≤ vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K r) h0)⁻¹ := by
    refine tateDvrVal_nonneg_of_mem _ ⟨((hr.unit⁻¹ : Rˣ) : R), ?_⟩
    show algebraMap R K ((hr.unit⁻¹ : Rˣ) : R) = _
    rw [Units.val_inv_eq_inv_val]
    show algebraMap R K ((hr.unit⁻¹ : Rˣ) : R) = (algebraMap R K r)⁻¹
    have h3 : algebraMap R K ((hr.unit⁻¹ : Rˣ) : R) * algebraMap R K r = 1 := by
      rw [← map_mul]
      have h4 : ((hr.unit⁻¹ : Rˣ) : R) * ((hr.unit : Rˣ) : R) = 1 := by
        rw [← Units.val_mul, inv_mul_cancel]
        rfl
      rw [hr.unit_spec] at h4
      rw [h4, map_one]
    field_simp at h3 ⊢
    linear_combination h3
  rw [vAdd_inv] at h2
  omega

/-- ★係数環から係数拡大した曲線は整である。 -/
theorem isIntegral_baseChange (E : WeierstrassCurve R) :
    WeierstrassCurve.IsIntegral R (E.baseChange K) := ⟨⟨E, rfl⟩⟩

/-! ## ★★★★★★★局所高さ = 極小判別式の付値 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**局所高さは極小判別式の付値に等しい**。

★★★第 319 の `deg∞`(`v_p(Δ_min)` の和)と第 310 の局所高さ(`v(q_E)`)を繋ぐ。 -/
theorem localHeight_eq_vAdd_Delta [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (W : WeierstrassCurve K) [hell : W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) :
    (localHeightOf W h : ℤ)
      = vAdd (tateDvrVal R K) (Units.mk0 W.Δ hell.isUnit.ne_zero) := by
  have hΔ : W.Δ ≠ 0 := hell.isUnit.ne_zero
  obtain ⟨hq, C, hqne, hCE⟩ := tateParamR_spec W h
  have h1 : (localHeightOf W h : ℤ) = vAdd (tateDvrVal R K) (tateParamK W h) := by
    have hpos := vAdd_tateParamK_pos W h
    rw [localHeightOf]
    omega
  have hWmap : (integralModel R W).map (algebraMap R K) = W :=
    WeierstrassCurve.baseChange_integralModel_eq R W
  have hkey : (C.map (algebraMap R K)) • W
      = (tateCurveAt (tateParamR W h) hq).map (algebraMap R K) := by
    conv_lhs => rw [← hWmap]
    rw [WeierstrassCurve.map_variableChange, hCE]
  -- ★★★Tate モデルは極小である
  have hTΔ : ((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).Δ ≠ 0 := by
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
      (tateCurveAt_Delta_ne_zero hq hqne)
  have hTc4ne : ((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).c₄ ≠ 0 := by
    rw [WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
      (tateCurveAt_c4_isUnit _ hq).ne_zero
  have hTc4 : vAdd (tateDvrVal R K)
      (Units.mk0 (((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).c₄) hTc4ne) = 0 := by
    have hveq : (Units.mk0 (((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)).c₄) hTc4ne)
        = Units.mk0 (algebraMap R K ((tateCurveAt (tateParamR W h) hq).c₄))
          (by rw [← WeierstrassCurve.map_c₄]; exact hTc4ne) := by
      refine Units.ext ?_
      exact WeierstrassCurve.map_c₄ _ _
    rw [hveq]
    exact tateDvrVal_eq_zero_of_isUnit _ (tateCurveAt_c4_isUnit _ hq) _
  haveI hTint : WeierstrassCurve.IsIntegral R
      ((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)) := isIntegral_baseChange _
  have hTmin : IsMinimal R ((tateCurveAt (tateParamR W h) hq).map (algebraMap R K)) :=
    isMinimal_of_c4_vAdd_eq_zero _ hTΔ hTc4ne hTc4
  -- ★★極小モデルどうしの `u` は付値 `0`
  have hmin1 : IsMinimal R ((1 : VariableChange K) • W) := by
    rw [one_smul]
    infer_instance
  have hmin2 : IsMinimal R ((C.map (algebraMap R K)) • W) := by
    rw [hkey]
    exact hTmin
  have hu0 := minimal_u_vAdd_eq W hΔ (1 : VariableChange K) (C.map (algebraMap R K)) hmin1 hmin2
  have hone : vAdd (tateDvrVal R K) ((1 : VariableChange K).u) = 0 := by
    show vAdd (tateDvrVal R K) 1 = 0
    rw [vAdd, map_one]
    rfl
  have hΔchg := vAdd_Delta_variableChange (R := R) W hΔ (C.map (algebraMap R K))
  -- ★★★`Δ_q = q·(単元)`
  obtain ⟨v, hv, hΔq⟩ := tateCurveAt_Delta_eq_mul_unit (tateParamR W h) hq
  have hveq2 : (Units.mk0 (((C.map (algebraMap R K)) • W).Δ)
      (variableChange_Delta_ne_zero W hΔ (C.map (algebraMap R K))))
      = tateParamK W h * Units.mk0 (algebraMap R K v)
        ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hv.ne_zero) := by
    refine Units.ext ?_
    show ((C.map (algebraMap R K)) • W).Δ = _
    rw [hkey, WeierstrassCurve.map_Δ, hΔq, map_mul]
    rfl
  rw [h1, ← hu0, hone] at *
  rw [hveq2, vAdd_mul, tateDvrVal_eq_zero_of_isUnit _ hv] at hΔchg
  omega

/-! ## ★★★★★★★★`v(q)` の一意性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Tate モデルの母数の付値は一意**

    `C • (整モデル) = tateCurveAt q` なら `v(Δ_W) = v(q)`

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★第 320 の `localHeight_eq_vAdd_Delta` は**選んだ** `tateParamR W h` について
述べていたが、本補題は**任意の** Tate モデルについて述べる。
したがって `v(q)` は `choose` の取り方に依らない。

★★★★☆**これで `Lemma 3.5` の最後の入力が最も鋭い形になる**——
`E′` が母数 `q^l` の Tate モデルを持つことさえ言えれば
`v(Δ_min(E′)) = v(q^l) = l·v(q) = l·v(Δ_min(E))` が出る。 -/
theorem vAdd_Delta_eq_of_tateModel [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (W : WeierstrassCurve K) [hell : W.IsElliptic] [W.IsMinimal R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0)
    (C : VariableChange R) (hC : C • integralModel R W = tateCurveAt q hq)
    (hqK : algebraMap R K q ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 W.Δ hell.isUnit.ne_zero)
      = vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K q) hqK) := by
  have hΔ : W.Δ ≠ 0 := hell.isUnit.ne_zero
  have hWmap : (integralModel R W).map (algebraMap R K) = W :=
    WeierstrassCurve.baseChange_integralModel_eq R W
  have hkey : (C.map (algebraMap R K)) • W
      = (tateCurveAt q hq).map (algebraMap R K) := by
    conv_lhs => rw [← hWmap]
    rw [WeierstrassCurve.map_variableChange, hC]
  have hTΔ : ((tateCurveAt q hq).map (algebraMap R K)).Δ ≠ 0 := by
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
      (tateCurveAt_Delta_ne_zero hq hq0)
  have hTc4ne : ((tateCurveAt q hq).map (algebraMap R K)).c₄ ≠ 0 := by
    rw [WeierstrassCurve.map_c₄]
    exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
      (tateCurveAt_c4_isUnit _ hq).ne_zero
  have hTc4 : vAdd (tateDvrVal R K)
      (Units.mk0 (((tateCurveAt q hq).map (algebraMap R K)).c₄) hTc4ne) = 0 := by
    have hveq : (Units.mk0 (((tateCurveAt q hq).map (algebraMap R K)).c₄) hTc4ne)
        = Units.mk0 (algebraMap R K ((tateCurveAt q hq).c₄))
          (by rw [← WeierstrassCurve.map_c₄]; exact hTc4ne) := by
      refine Units.ext ?_
      exact WeierstrassCurve.map_c₄ _ _
    rw [hveq]
    exact tateDvrVal_eq_zero_of_isUnit _ (tateCurveAt_c4_isUnit _ hq) _
  haveI hTint : WeierstrassCurve.IsIntegral R
      ((tateCurveAt q hq).map (algebraMap R K)) := isIntegral_baseChange _
  have hTmin : IsMinimal R ((tateCurveAt q hq).map (algebraMap R K)) :=
    isMinimal_of_c4_vAdd_eq_zero _ hTΔ hTc4ne hTc4
  have hmin1 : IsMinimal R ((1 : VariableChange K) • W) := by
    rw [one_smul]
    infer_instance
  have hmin2 : IsMinimal R ((C.map (algebraMap R K)) • W) := by
    rw [hkey]
    exact hTmin
  have hu0 := minimal_u_vAdd_eq W hΔ (1 : VariableChange K) (C.map (algebraMap R K))
    hmin1 hmin2
  have hone : vAdd (tateDvrVal R K) ((1 : VariableChange K).u) = 0 := by
    show vAdd (tateDvrVal R K) 1 = 0
    rw [vAdd, map_one]
    rfl
  have hΔchg := vAdd_Delta_variableChange (R := R) W hΔ (C.map (algebraMap R K))
  obtain ⟨v, hv, hΔq⟩ := tateCurveAt_Delta_eq_mul_unit q hq
  have hveq2 : (Units.mk0 (((C.map (algebraMap R K)) • W).Δ)
      (variableChange_Delta_ne_zero W hΔ (C.map (algebraMap R K))))
      = Units.mk0 (algebraMap R K q) hqK * Units.mk0 (algebraMap R K v)
        ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hv.ne_zero) := by
    refine Units.ext ?_
    show ((C.map (algebraMap R K)) • W).Δ = _
    rw [hkey, WeierstrassCurve.map_Δ, hΔq, map_mul]
    rfl
  rw [← hu0, hone, mul_zero, sub_zero] at hΔchg
  rw [hveq2, vAdd_mul, tateDvrVal_eq_zero_of_isUnit _ hv, add_zero] at hΔchg
  omega

def vAdd_Delta_eq_of_tateModel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate モデルの母数の付値は一意。★無条件)",
    sectionId := "genell-def-3-3" }

/-! ## ★★★★★★★★★★★★`primeSubring` の言葉での極小性判定 -/

section PrimeSubring

variable {L : Type} [Field L] [NumberField L]

open IsDedekindDomain NumberField

/-- ☆ 2 つの加法付値は同時に `0` になる——付値が同値だから。 -/
theorem vAdd_eq_zero_iff_valAdd_eq_zero (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    vAdd (tateDvrVal (primeSubring p) L) x = 0 ↔ valAdd p x = 0 := by
  rw [vAdd_eq_zero_iff, valAdd_eq_zero_iff]
  exact Valuation.isEquiv_iff_val_eq_one.1 (valuation_isEquiv p)

/-- ★★★★★★★★★★★★★★★★
**`v_p(c₄) = 0` なら `primeSubring p` 上極小**——★**無条件**（第 1190）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`isMinimal_of_c4_vAdd_eq_zero`（第 320）は完備 DVR と `tateDvrVal` の言葉だったが、
`primeSubring p` は DVR で分数体が `L` なので、そのまま当てはまる。
★2 つの加法付値が同時に `0` になることだけが橋である。

★★★これが `Skeleton/GenEll/LCyclicReading.lean` の節点 2d-1 で
第 1189 が名指しした**最後の葉**である。 -/
theorem isMinimal_of_c4_valAdd_eq_zero (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    (hΔ : W.Δ ≠ 0) (hc4 : W.c₄ ≠ 0)
    (h : valAdd p (Units.mk0 W.c₄ hc4) = 0) :
    WeierstrassCurve.IsMinimal (primeSubring p) W :=
  isMinimal_of_c4_vAdd_eq_zero W hΔ hc4 ((vAdd_eq_zero_iff_valAdd_eq_zero p _).2 h)

def isMinimal_of_c4_valAdd_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(c₄) = 0 なら primeSubring p 上極小。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_of_c4_valAdd_eq_zero.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isMinimal_of_c4_vAdd_eq_zero(第 320、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_of_c4_vAdd_eq_zero") 1,
    .citation "[mathlib]" "Valuation.isEquiv_iff_val_eq_one(同値な付値は同時に 1)"
      (.inMathlib "Valuation.isEquiv_iff_val_eq_one") 1 ]

end PrimeSubring

/-! ## ★出典の紐付け(`.src`) -/

def localHeight_eq_vAdd_Delta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(局所高さ = 極小判別式の付値)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
