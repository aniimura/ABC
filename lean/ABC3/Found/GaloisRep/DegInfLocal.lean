import ABC3.Found.GaloisRep.LocalHeightDelta

/-!
# Galois (G8) 第 321 ブロック —— **★★★★★★★`d·deg∞ ≥ (局所高さ)·log 2`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★到達点

> **`d · deg∞(E) ≥ (局所高さ) · log 2`**(`degInfOf_ge_localHeight`)

★★★これが `FaltingsHeightData` の `degInf_ge_localHeight` 欄の中身である
(第 318 で入れた**不分岐の仮定つき**の形)。

## ★★★★★証明の筋——3 段

1. ★★★★**局所高さ ≤ `minDeltaExp p E`**(`localHeight_le_minDeltaExp`)
   ——`E` を `p` で極小化する変数変換 `C` を取り、それを `Lv` に移す。
   `C•E` は `p` 上整だから、不分岐の仮定で `R` 上も整(`isIntegral_baseChange_of_isIntegral`)。
   `E⁄Lv` は `R` 上極小なので `v_R(Δ(E⁄Lv)) ≤ v_R(Δ((C•E)⁄Lv)) = v_p(Δ(C•E)) = minDeltaExp p E`。
   ★左辺は第 320 で局所高さに等しい。
2. ★★`log N(p) ≥ log 2`(`log_two_le_log_absNorm`)——素イデアルのノルムは `2` 以上。
3. ★★`minDeltaExp` は各項非負だから、**単項 ≤ 有限和**(`single_le_finsum`)。

## ★★★★★★不分岐の仮定が効く場所は 2 つだけ

    hp : ∀ x : L, v_R(algebraMap L Lv x) = v_p(x)

* ★`p` 上整 ⟹ `R` 上整(`isIntegral_baseChange_of_isIntegral`)
* ★`v_R` の指数 = `v_p` の指数(`vAdd_algebraMap_eq_valAdd`)

★★★第 318 で見つけた「分岐で偽になる」という欠陥は、この 2 箇所で塞がれている。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_preimage_of_valuation_le_one` | ★付値 `≤ 1` の元は係数環から来る |
| `exists_preimage_of_mem_primeSubring` | ★★`p` の整数環の元は `R` から来る |
| `isIntegral_baseChange_of_isIntegral` | ★★★★**`p` 上整 ⟹ `R` 上整** |
| `withZero_log_coe_ofAdd` | ★`log ∘ ofAdd = id` |
| `vAdd_algebraMap_eq_valAdd` | ★★★★**付値指数の一致** |
| `minimal_vAdd_Delta_le` | ★★★★極小モデルの `v(Δ)` は最小 |
| `localHeight_le_minDeltaExp` | ★★★★★★**局所高さ ≤ `v_p(Δ_min)`** |
| `two_le_absNorm`・`log_two_le_log_absNorm` | ★★`N(p) ≥ 2` |
| `finrank_mul_degInfOf`・`hasFiniteSupport_degInf` | ★★配管 |
| `degInfOf_ge_localHeight` | ★★★★★★★**`d·deg∞ ≥ (局所高さ)·log 2`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField
open scoped Classical

/-! ## ★★★★極小モデルの判別式の付値は最小である -/

section Minimal

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★極小モデルの判別式の付値指数は、どの整モデルのそれよりも小さい。

★`IsMinimal` は `MaximalFor`——「上に出るものはない」という形なので、
値が**線形順序**であることを使って「常に下」に読み替える。 -/
theorem minimal_vAdd_Delta_le (W : WeierstrassCurve K) [hmin : IsMinimal R W]
    (hΔ : W.Δ ≠ 0) (C : VariableChange K) (hC : IsIntegral R (C • W)) :
    vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ)
      ≤ vAdd (tateDvrVal R K) (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C)) := by
  haveI hint : IsIntegral R W := inferInstance
  obtain ⟨h1, h2⟩ := (WeierstrassCurve.isMinimal_iff R W).1 hmin
  rw [← valuation_le_iff_vAdd_le]
  have hstep : WeierstrassCurve.valuation_Δ_aux R (C • W)
      ≤ WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W) := by
    by_cases hcase : WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W)
        ≤ WeierstrassCurve.valuation_Δ_aux R (C • W)
    · exact h2 hC hcase
    · exact (not_le.1 hcase).le
  rw [one_smul] at hstep
  rw [WeierstrassCurve.valuation_Δ_aux, dif_pos hC, WeierstrassCurve.valuation_Δ_aux,
    dif_pos hint] at hstep
  exact hstep

end Minimal

/-! ## ★★★★不分岐の仮定——`p` 上整ならば `R` 上整 -/

variable {L : Type} [Field L] [NumberField L]
  {Lv : Type} [Field Lv] [Algebra L Lv]
  {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [Algebra R Lv] [IsFractionRing R Lv]

/-- ★付値が `1` 以下の元は係数環から来る(`0` も含めて)。 -/
theorem exists_preimage_of_valuation_le_one (x : Lv)
    (h : (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R)) x ≤ 1) :
    ∃ r : R, algebraMap R Lv r = x := by
  by_cases hx : x = 0
  · exact ⟨0, by rw [map_zero, hx]⟩
  · refine dvr_mem_of_nonneg (Units.mk0 x hx) ?_
    rw [tateDvrVal_nonneg_iff]
    exact h

/-- ★★不分岐の仮定のもとで、`p` の整数環の元は `R` から来る。 -/
theorem exists_preimage_of_mem_primeSubring (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (x : L) (hx : x ∈ primeSubring p) :
    ∃ r : R, algebraMap R Lv r = algebraMap L Lv x := by
  refine exists_preimage_of_valuation_le_one (R := R) _ ?_
  rw [hp x]
  exact (mem_primeSubring_iff p x).1 hx

set_option maxHeartbeats 1600000 in
/-- ★★★★**不分岐の仮定のもとで、`p` 上整な曲線は `R` 上整になる**。

★係数を 1 つずつ `R` に持ち上げて、`R` 上の Weierstrass 曲線を組み立てる。 -/
theorem isIntegral_baseChange_of_isIntegral (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (V : WeierstrassCurve L) [hV : WeierstrassCurve.IsIntegral (primeSubring p) V] :
    WeierstrassCurve.IsIntegral R (V.baseChange Lv) := by
  obtain ⟨r1, hr1⟩ := exists_preimage_of_mem_primeSubring (Lv := Lv) (R := R) p hp V.a₁
    (by rw [← WeierstrassCurve.integralModel_a₁_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₁).2)
  obtain ⟨r2, hr2⟩ := exists_preimage_of_mem_primeSubring (Lv := Lv) (R := R) p hp V.a₂
    (by rw [← WeierstrassCurve.integralModel_a₂_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₂).2)
  obtain ⟨r3, hr3⟩ := exists_preimage_of_mem_primeSubring (Lv := Lv) (R := R) p hp V.a₃
    (by rw [← WeierstrassCurve.integralModel_a₃_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₃).2)
  obtain ⟨r4, hr4⟩ := exists_preimage_of_mem_primeSubring (Lv := Lv) (R := R) p hp V.a₄
    (by rw [← WeierstrassCurve.integralModel_a₄_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₄).2)
  obtain ⟨r6, hr6⟩ := exists_preimage_of_mem_primeSubring (Lv := Lv) (R := R) p hp V.a₆
    (by rw [← WeierstrassCurve.integralModel_a₆_eq (primeSubring p) V]
        exact ((WeierstrassCurve.integralModel (primeSubring p) V).a₆).2)
  refine ⟨⟨⟨r1, r2, r3, r4, r6⟩, ?_⟩⟩
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_
  · show algebraMap L Lv V.a₁ = algebraMap R Lv r1
    exact hr1.symm
  · show algebraMap L Lv V.a₂ = algebraMap R Lv r2
    exact hr2.symm
  · show algebraMap L Lv V.a₃ = algebraMap R Lv r3
    exact hr3.symm
  · show algebraMap L Lv V.a₄ = algebraMap R Lv r4
    exact hr4.symm
  · show algebraMap L Lv V.a₆ = algebraMap R Lv r6
    exact hr6.symm

/-! ## ★★★★付値指数の一致 -/

/-- ★`WithZero` の `log` は `ofAdd` を戻す。 -/
theorem withZero_log_coe_ofAdd (k : ℤ) :
    (((Multiplicative.ofAdd k : Multiplicative ℤ) : WithZero (Multiplicative ℤ))).log = k := rfl

/-- ★★★★不分岐の仮定のもとで、`R` の付値指数と `p` の付値指数は一致する。 -/
theorem vAdd_algebraMap_eq_valAdd (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (x : L) (hx : x ≠ 0) (hx' : algebraMap L Lv x ≠ 0) :
    vAdd (tateDvrVal R Lv) (Units.mk0 (algebraMap L Lv x) hx') = valAdd p (Units.mk0 x hx) := by
  have h1 := valuation_eq_ofAdd_neg_vAdd (R := R) (Units.mk0 (algebraMap L Lv x) hx')
  have h1' : (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
      (algebraMap L Lv x) = _ := h1
  rw [valAdd_eq_neg_log]
  have h2 : ((Units.mk0 x hx : Lˣ) : L) = x := rfl
  rw [h2, ← hp x, h1', withZero_log_coe_ofAdd]
  omega

/-! ## ★★★★★★局所高さ ≤ `v_p(Δ_min)` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**局所高さは極小判別式の指数で上から抑えられる**。

★★`Lv` の側での極小性は `L` の側での極小性より強い(変数変換が多いから)ので、
不等式は「局所高さ ≤ `v_p(Δ_min)`」の向きに出る。 -/
theorem localHeight_le_minDeltaExp [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E : WeierstrassCurve L) [hell : (E.baseChange Lv).IsElliptic]
    [hmin : (E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x) :
    (localHeightOf (E.baseChange Lv) h : ℤ) ≤ minDeltaExp p E := by
  have hΔv : (E.baseChange Lv).Δ ≠ 0 := hell.isUnit.ne_zero
  have hmapΔ : (E.baseChange Lv).Δ = algebraMap L Lv E.Δ := WeierstrassCurve.map_Δ _ _
  have hΔL : E.Δ ≠ 0 := by
    intro h0
    exact hΔv (by rw [hmapΔ, h0, map_zero])
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  have hCΔ : (C • E).Δ ≠ 0 := variableChange_Delta_ne_zero E hΔL C
  have hmde := minDeltaExp_eq p E hΔL C hC
  have h320 := localHeight_eq_vAdd_Delta (R := R) (E.baseChange Lv) h
  have hsm : (C.map (algebraMap L Lv)) • (E.baseChange Lv) = (C • E).baseChange Lv :=
    WeierstrassCurve.map_variableChange E C (algebraMap L Lv)
  haveI hCint : WeierstrassCurve.IsIntegral (primeSubring p) (C • E) := inferInstance
  have hint : WeierstrassCurve.IsIntegral R ((C.map (algebraMap L Lv)) • (E.baseChange Lv)) := by
    rw [hsm]
    exact isIntegral_baseChange_of_isIntegral p hp (C • E)
  have hle := minimal_vAdd_Delta_le (R := R) (E.baseChange Lv) hΔv
    (C.map (algebraMap L Lv)) hint
  have hΔeq : (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).Δ)
      = algebraMap L Lv ((C • E).Δ) := by
    rw [hsm]
    exact WeierstrassCurve.map_Δ _ _
  have hne' : algebraMap L Lv ((C • E).Δ) ≠ 0 := by
    rw [← hΔeq]
    exact variableChange_Delta_ne_zero (E.baseChange Lv) hΔv _
  have hunit : (Units.mk0 (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).Δ)
      (variableChange_Delta_ne_zero (E.baseChange Lv) hΔv (C.map (algebraMap L Lv))))
      = Units.mk0 (algebraMap L Lv ((C • E).Δ)) hne' := Units.ext hΔeq
  rw [hunit, vAdd_algebraMap_eq_valAdd (R := R) p hp ((C • E).Δ) hCΔ hne'] at hle
  rw [h320, hmde]
  exact hle

/-! ## ★★★★★★★★局所高さ = `v_p(Δ_min)`（半安定な場合） -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定なら局所高さは `v_p(Δ_min)` に等しい**

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★★不等式（`localHeight_le_minDeltaExp`）は「`Lv` の側の極小性が強い」から出るが、
**乗法還元では `v(c₄) = 0` が極小性を決める**（`isMinimal_of_c4_vAdd_eq_zero`）ので、
`p`-極小なモデルはそのまま `R`-極小でもある。したがって等号になる。

★★★★☆**これが `Lemma 3.5` の最後の入力
（`v_p(Δ_min(E′)) = l·v_p(Δ_min(E))`）を局所高さの言葉へ移す環である**——
局所高さの側は `Lemma 3.2, (ii)`（`q_{E′} = q_E^l`）が与える。 -/
theorem localHeight_eq_minDeltaExp [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E : WeierstrassCurve L) [hell : (E.baseChange Lv).IsElliptic]
    [hmin : (E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x)
    (C : VariableChange L) (hC : IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0)
    (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0) :
    (localHeightOf (E.baseChange Lv) h : ℤ) = minDeltaExp p E := by
  have hΔv : (E.baseChange Lv).Δ ≠ 0 := hell.isUnit.ne_zero
  have hmapΔ : (E.baseChange Lv).Δ = algebraMap L Lv E.Δ := WeierstrassCurve.map_Δ _ _
  have hΔL : E.Δ ≠ 0 := by
    intro h0
    exact hΔv (by rw [hmapΔ, h0, map_zero])
  have hCΔ : (C • E).Δ ≠ 0 := variableChange_Delta_ne_zero E hΔL C
  have hmde := minDeltaExp_eq p E hΔL C hC
  have hsm : (C.map (algebraMap L Lv)) • (E.baseChange Lv) = (C • E).baseChange Lv :=
    WeierstrassCurve.map_variableChange E C (algebraMap L Lv)
  haveI hCint : WeierstrassCurve.IsIntegral (primeSubring p) (C • E) := inferInstance
  haveI hint : WeierstrassCurve.IsIntegral R
      ((C.map (algebraMap L Lv)) • (E.baseChange Lv)) := by
    rw [hsm]
    exact isIntegral_baseChange_of_isIntegral p hp (C • E)
  -- ★`c₄` の付値が `0` であること（`Lv` の側）
  have hc4eq : ((C.map (algebraMap L Lv)) • (E.baseChange Lv)).c₄
      = algebraMap L Lv ((C • E).c₄) := by
    rw [hsm]
    exact WeierstrassCurve.map_c₄ _ _
  have hne2 : algebraMap L Lv ((C • E).c₄) ≠ 0 :=
    (map_ne_zero_iff _ (algebraMap L Lv).injective).2 hc4ne
  have hc4ne' : ((C.map (algebraMap L Lv)) • (E.baseChange Lv)).c₄ ≠ 0 := by
    rw [hc4eq]; exact hne2
  have hc4v : vAdd (tateDvrVal R Lv)
      (Units.mk0 (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).c₄) hc4ne') = 0 := by
    have hu : (Units.mk0 (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).c₄) hc4ne')
        = Units.mk0 (algebraMap L Lv ((C • E).c₄)) hne2 := Units.ext hc4eq
    rw [hu, vAdd_algebraMap_eq_valAdd (R := R) p hp ((C • E).c₄) hc4ne hne2]
    exact hc4
  have hΔne' : ((C.map (algebraMap L Lv)) • (E.baseChange Lv)).Δ ≠ 0 :=
    variableChange_Delta_ne_zero (E.baseChange Lv) hΔv _
  -- ★★★`v(c₄) = 0` なので `R`-極小
  have hRmin : IsMinimal R ((C.map (algebraMap L Lv)) • (E.baseChange Lv)) :=
    isMinimal_of_c4_vAdd_eq_zero _ hΔne' hc4ne' hc4v
  have hmin1 : IsMinimal R ((1 : VariableChange Lv) • (E.baseChange Lv)) := by
    rw [one_smul]
    infer_instance
  have hu0 := minimal_u_vAdd_eq (E.baseChange Lv) hΔv (1 : VariableChange Lv)
    (C.map (algebraMap L Lv)) hmin1 hRmin
  have hone : vAdd (tateDvrVal R Lv) ((1 : VariableChange Lv).u) = 0 := by
    show vAdd (tateDvrVal R Lv) 1 = 0
    rw [vAdd, map_one]
    rfl
  have hΔchg := vAdd_Delta_variableChange (R := R) (E.baseChange Lv) hΔv
    (C.map (algebraMap L Lv))
  rw [← hu0, hone, mul_zero, sub_zero] at hΔchg
  -- ★`v(Δ)` を `p` の側へ
  have hΔeq : (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).Δ)
      = algebraMap L Lv ((C • E).Δ) := by
    rw [hsm]; exact WeierstrassCurve.map_Δ _ _
  have hne3 : algebraMap L Lv ((C • E).Δ) ≠ 0 := by rw [← hΔeq]; exact hΔne'
  have hΔu : (Units.mk0 (((C.map (algebraMap L Lv)) • (E.baseChange Lv)).Δ) hΔne')
      = Units.mk0 (algebraMap L Lv ((C • E).Δ)) hne3 := Units.ext hΔeq
  have h320 := localHeight_eq_vAdd_Delta (R := R) (E.baseChange Lv) h
  rw [h320, hmde, ← hΔchg, hΔu,
    vAdd_algebraMap_eq_valAdd (R := R) p hp ((C • E).Δ) hCΔ hne3]

def localHeight_eq_minDeltaExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(半安定なら局所高さ = v_p(Δ_min)。★無条件)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★素イデアルのノルムは `2` 以上 -/

/-- ★★`N(p) ≥ 2`——`p` は `⊥` でも `⊤` でもないから。 -/
theorem two_le_absNorm (p : HeightOneSpectrum (𝓞 L)) : 2 ≤ Ideal.absNorm p.asIdeal := by
  have h0 : Ideal.absNorm p.asIdeal ≠ 0 := by
    rw [Ne, Ideal.absNorm_eq_zero_iff]
    exact p.ne_bot
  have h1 : Ideal.absNorm p.asIdeal ≠ 1 := by
    rw [Ne, Ideal.absNorm_eq_one_iff]
    exact p.isPrime.ne_top'
  omega

/-- ★★`log 2 ≤ log N(p)`。 -/
theorem log_two_le_log_absNorm (p : HeightOneSpectrum (𝓞 L)) :
    Real.log 2 ≤ Real.log (Ideal.absNorm p.asIdeal) := by
  refine Real.log_le_log (by norm_num) ?_
  exact_mod_cast two_le_absNorm p

/-! ## ★★配管 -/

set_option maxHeartbeats 1600000 in
theorem finrank_mul_degInfOf (E : WeierstrassCurve L) :
    (Module.finrank ℚ L : ℝ) * degInfOf L E
      = ∑ᶠ q : HeightOneSpectrum (𝓞 L),
          (minDeltaExp q E : ℝ) * Real.log (Ideal.absNorm q.asIdeal) := by
  have hd : (Module.finrank ℚ L : ℝ) ≠ 0 := by
    have := Module.finrank_pos (R := ℚ) (M := L)
    positivity
  rw [degInfOf, mul_div_cancel₀ _ hd]

theorem hasFiniteSupport_degInf (E : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) :
    Function.HasFiniteSupport
      (fun q : HeightOneSpectrum (𝓞 L) =>
        (minDeltaExp q E : ℝ) * Real.log (Ideal.absNorm q.asIdeal)) := by
  refine Set.Finite.subset (minDeltaExp_finite E hΔ) ?_
  intro q hq
  simp only [Function.mem_support, ne_eq, Set.mem_setOf_eq] at hq ⊢
  intro h0
  exact hq (by rw [h0]; push_cast; ring)

/-! ## ★★★★★★★`d·deg∞ ≥ (局所高さ)·log 2` -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**`d · deg∞(E) ≥ (局所高さ) · log 2`**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★★★`FaltingsHeightData` の `degInf_ge_localHeight` 欄の中身である。 -/
theorem degInfOf_ge_localHeight [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (E : WeierstrassCurve L) [hell : (E.baseChange Lv).IsElliptic]
    [hmin : (E.baseChange Lv).IsMinimal R]
    (h : (E.baseChange Lv).HasSplitMultiplicativeReduction R)
    (p : HeightOneSpectrum (𝓞 L))
    (hp : ∀ x : L, (HeightOneSpectrum.valuation Lv (IsDiscreteValuationRing.maximalIdeal R))
        (algebraMap L Lv x) = (HeightOneSpectrum.valuation L p) x) :
    (Module.finrank ℚ L : ℝ) * degInfOf L E
      ≥ (localHeightOf (E.baseChange Lv) h : ℝ) * Real.log 2 := by
  have hΔv : (E.baseChange Lv).Δ ≠ 0 := hell.isUnit.ne_zero
  have hmapΔ : (E.baseChange Lv).Δ = algebraMap L Lv E.Δ := WeierstrassCurve.map_Δ _ _
  have hΔL : E.Δ ≠ 0 := by
    intro h0
    exact hΔv (by rw [hmapΔ, h0, map_zero])
  have hkey : (localHeightOf (E.baseChange Lv) h : ℤ) ≤ minDeltaExp p E :=
    localHeight_le_minDeltaExp E h p hp
  rw [ge_iff_le, finrank_mul_degInfOf E]
  calc (localHeightOf (E.baseChange Lv) h : ℝ) * Real.log 2
      ≤ (minDeltaExp p E : ℝ) * Real.log 2 := by
        refine mul_le_mul_of_nonneg_right ?_ (Real.log_nonneg (by norm_num))
        exact_mod_cast hkey
    _ ≤ (minDeltaExp p E : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
        refine mul_le_mul_of_nonneg_left (log_two_le_log_absNorm p) ?_
        exact_mod_cast minDeltaExp_nonneg p E
    _ ≤ ∑ᶠ q : HeightOneSpectrum (𝓞 L),
          (minDeltaExp q E : ℝ) * Real.log (Ideal.absNorm q.asIdeal) := by
        refine single_le_finsum p (hasFiniteSupport_degInf E hΔL) (fun q => ?_)
        exact mul_nonneg (by exact_mod_cast minDeltaExp_nonneg q E) (log_absNorm_nonneg q)

/-! ## ★出典の紐付け(`.src`) -/

def degInfOf_ge_localHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(d·deg∞ ≥ (局所高さ)·log 2)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
