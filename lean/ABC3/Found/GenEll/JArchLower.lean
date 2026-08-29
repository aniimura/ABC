/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.JArchBound
import Mathlib.NumberTheory.ModularForms.EisensteinSeries.QExpansion
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★逆向き——`‖j‖·‖Δ‖` は下に有界（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★これは何か —— `Prop 3.4` の**第 3 の `≲`**

`§9-1000〜1001` は **上から**の評価

    `‖j(τ)‖·(‖Δ‖_Pet(τ))^{1+ϵ} ≤ M`

を取った。★★本ファイルは**逆向き**——`ϵ = 0` での**下から**の評価

    **`max(‖j(τ)‖, 1) · ‖Δ‖_Pet(τ) ≥ m > 0`**

を取る。★★★これが `Prop 3.4` の第 3 の `≲`
（`12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht_∞`）のアルキメデス側である。

## ★なぜ `max(‖j‖, 1)` なのか —— `log⁺` だから

高さは `log⁺|j| = max(log|j|, 0) = log(max(|j|, 1))` を足す。
★`τ = ρ = e^{2πi/3}` では `j = 0` なので `‖j‖·pD = 0` になってしまうが、
`max(‖j‖,1)·pD = pD(ρ) > 0` で破綻しない。

## ★★機構 —— カスプで `E₄ → 1`

`‖j(τ)‖·pD(τ) = ‖E₄(τ)‖³·(Im τ)⁶`（`ϵ = 0` の分解）。

★カスプでは `E₄ → 1`（正規化 Eisenstein 級数の定数項が `1`——
mathlib の `EisensteinSeries.E_qExpansion_coeff_zero`）なので、
`Im τ ≥ 1` かつ `‖E₄‖ > 1/2` の領域では `‖j‖·pD > 1/8`。
★★コンパクト側（切り詰めた基本領域）では連続な正値関数の**最小値**を取る。
★★★あとは `SL(2,ℤ)` 不変性。

## ★★上からの評価との対称性

| 向き | 主張 | どこ |
|---|---|---|
| 上から | `‖j‖·pD^{1+ϵ} ≤ M`（`ϵ > 0` が要る） | `JPeterssonBound.lean`（§9-1000） |
| 下から | `max(‖j‖,1)·pD ≥ m > 0`（`ϵ = 0` でよい） | ★本ファイル |

☆`ϵ` の非対称性が原文の `12(1+ϵ)` の出どころである。
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane Filter Topology Asymptotics

/-! ## ★★★★★★★カスプでの `E₄` -/

/-- ★★★★★**`E₄` はカスプで `valueAtInfty` に収束する**。

★保型形式は「値 − カスプ値」が指数的に減衰する
（`ModularFormClass.exp_decay_sub_atImInfty'`）。 -/
theorem tendsto_E4_atImInfty :
    Tendsto (fun τ : ℍ => (ModularForm.E₄ τ)) atImInfty
      (𝓝 (valueAtInfty (ModularForm.E₄ : ℍ → ℂ))) := by
  obtain ⟨c, hc, hO⟩ := ModularFormClass.exp_decay_sub_atImInfty' (ModularForm.E₄)
  have h0 := tendsto_rpow_mul_exp_neg_mul_atTop_nhds_zero 0 c hc
  simp only [Real.rpow_zero, one_mul] at h0
  have hlim : Tendsto (fun τ : ℍ => Real.exp (-c * τ.im)) atImInfty (𝓝 0) :=
    h0.comp tendsto_comap
  have h := hO.trans_tendsto hlim
  simpa using h.add_const (valueAtInfty (ModularForm.E₄ : ℍ → ℂ))

/-- ★★★★★★**`E₄` のカスプ値は `1`**——正規化 Eisenstein 級数の定数項。 -/
theorem valueAtInfty_E4 : valueAtInfty (ModularForm.E₄ : ℍ → ℂ) = 1 := by
  have hper := SlashInvariantFormClass.periodic_comp_ofComplex
    (f := (ModularForm.E₄ : ModularForm _ 4)) (h := (1:ℝ)) one_mem_strictPeriods_SL
  have hana := UpperHalfPlane.analyticAt_cuspFunction_zero (h := (1:ℝ)) one_pos hper
    (ModularFormClass.holo ModularForm.E₄) (ModularFormClass.bdd_at_infty ModularForm.E₄)
  rw [← UpperHalfPlane.qExpansion_coeff_zero (h := 1) one_pos hana hper]
  exact EisensteinSeries.E_qExpansion_coeff_zero (by norm_num) (by decide)

/-- ★★★★★★★**`‖E₄‖ → 1`（カスプ）**。 -/
theorem tendsto_norm_E4 :
    Tendsto (fun τ : ℍ => ‖ModularForm.E₄ τ‖) atImInfty (𝓝 1) := by
  have h := tendsto_E4_atImInfty
  rw [valueAtInfty_E4] at h
  simpa using h.norm

/-! ## ★★★★`ϵ = 0` の分解 -/

/-- ★★★★**`‖j‖·pD = ‖E₄‖³·(Im τ)⁶`**（`ϵ = 0`）。 -/
theorem norm_jFun_mul_peterssonDelta (τ : ℍ) :
    ‖jFun τ‖ * peterssonDelta τ = ‖ModularForm.E₄ τ‖ ^ 3 * τ.im ^ 6 := by
  have hΔ : ‖ModularForm.discriminant τ‖ ≠ 0 :=
    norm_ne_zero_iff.2 (ModularForm.discriminant_ne_zero τ)
  have hj : ‖jFun τ‖ = ‖ModularForm.E₄ τ‖ ^ 3 / ‖ModularForm.discriminant τ‖ := by
    rw [jFun, norm_div, norm_pow]
  rw [hj, peterssonDelta]
  field_simp

/-! ## ★★★★★不変性と連続性 -/

theorem maxJ_pD_smul (g : Matrix.SpecialLinearGroup (Fin 2) ℤ) (τ : ℍ) :
    max ‖jFun (g • τ)‖ 1 * peterssonDelta (g • τ)
      = max ‖jFun τ‖ 1 * peterssonDelta τ := by
  rw [jFun_smul, peterssonDelta_smul]

theorem continuous_maxJ_pD : Continuous (fun τ : ℍ => max ‖jFun τ‖ 1 * peterssonDelta τ) :=
  (continuous_jFun.norm.max continuous_const).mul continuous_peterssonDelta

/-! ## ★★★★★★★★★★★★★★★★★★★★下からの評価 -/

/-- ★★★★★★★★★★★★★★★★★★★★**`max(‖j(τ)‖,1)·‖Δ‖_Pet(τ) ≥ m > 0`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★カスプ側（`Im τ ≥ 1` かつ `‖E₄‖ > 1/2`）では `‖j‖·pD = ‖E₄‖³(Im τ)⁶ > 1/8`。
★★コンパクト側（切り詰めた基本領域）では連続な正値関数の最小値。
★★★あとは `SL(2,ℤ)` 不変性。 -/
theorem exists_lower_bound_maxJ_pD :
    ∃ m : ℝ, 0 < m ∧ ∀ τ : ℍ, m ≤ max ‖jFun τ‖ 1 * peterssonDelta τ := by
  obtain ⟨y₀, hy₀⟩ := (UpperHalfPlane.atImInfty_mem _).1
    (tendsto_norm_E4.eventually_const_lt (by norm_num : (1:ℝ)/2 < 1))
  set y₁ := max y₀ 1 with hy₁
  have hcpt := ModularGroup.isCompact_truncatedFundamentalDomain y₁
  have hne : (ModularGroup.truncatedFundamentalDomain y₁).Nonempty := by
    refine ⟨UpperHalfPlane.I, ModularGroup.I_mem_fd, ?_⟩
    simp [hy₁]
  obtain ⟨x₀, hx₀mem, hx₀min⟩ :=
    hcpt.exists_isMinOn hne continuous_maxJ_pD.continuousOn
  have hx₀pos : 0 < max ‖jFun x₀‖ 1 * peterssonDelta x₀ :=
    mul_pos (lt_of_lt_of_le zero_lt_one (le_max_right _ _)) (peterssonDelta_pos x₀)
  refine ⟨min (max ‖jFun x₀‖ 1 * peterssonDelta x₀) (1/8),
    lt_min hx₀pos (by norm_num), fun τ => ?_⟩
  obtain ⟨g, hg⟩ := ModularGroup.exists_smul_mem_fd τ
  have hinv := maxJ_pD_smul g τ
  rcases le_or_gt ((g • τ).im) y₁ with hle | hgt
  · have hm : max ‖jFun x₀‖ 1 * peterssonDelta x₀
        ≤ max ‖jFun (g • τ)‖ 1 * peterssonDelta (g • τ) := hx₀min ⟨hg, hle⟩
    rw [hinv] at hm
    exact le_trans (min_le_left _ _) hm
  · have hE : (1:ℝ)/2 < ‖ModularForm.E₄ (g • τ)‖ :=
      hy₀ (g • τ) (le_of_lt (lt_of_le_of_lt (le_max_left y₀ 1) hgt))
    have hy : (1:ℝ) ≤ (g • τ).im := le_of_lt (lt_of_le_of_lt (le_max_right y₀ 1) hgt)
    have hkey : (1:ℝ)/8 ≤ max ‖jFun (g • τ)‖ 1 * peterssonDelta (g • τ) := by
      have h1 : ‖jFun (g • τ)‖ * peterssonDelta (g • τ)
          ≤ max ‖jFun (g • τ)‖ 1 * peterssonDelta (g • τ) :=
        mul_le_mul_of_nonneg_right (le_max_left _ _) (peterssonDelta_pos _).le
      refine le_trans ?_ h1
      rw [norm_jFun_mul_peterssonDelta]
      have h2 : (1:ℝ)/8 ≤ ‖ModularForm.E₄ (g • τ)‖ ^ 3 := by
        have h := pow_le_pow_left₀ (by norm_num : (0:ℝ) ≤ 1/2) hE.le 3
        norm_num at h ⊢
        linarith
      have h3 : (1:ℝ) ≤ (g • τ).im ^ 6 := one_le_pow₀ hy
      nlinarith
    rw [hinv] at hkey
    exact le_trans (min_le_right _ _) hkey

/-! ## ★★★★★★★★★★曲線の言葉へ -/

/-- ★`log⁺` を `max` の中へ入れる。 -/
theorem max_log_eq_log_max {x : ℝ} (hx : 0 ≤ x) : max (Real.log x) 0 = Real.log (max x 1) := by
  rcases le_or_gt x 1 with h | h
  · rw [max_eq_right h, Real.log_one, max_eq_right (Real.log_nonpos hx h)]
  · rw [max_eq_left h.le, max_eq_left (Real.log_nonneg h.le)]

/-- ★★★★★★★★★★**`log⁺‖j(W)‖ + log((2π)¹²·curveArchInv W) ≥ −C`**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`§9-1001` の `exists_logPlus_j_bound`（上から、`(1+ϵ)` つき）の**逆向き**である
——★★こちらは `ϵ` が要らない。 -/
theorem exists_logPlus_j_lower :
    ∃ C : ℝ, ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      -C ≤ max (Real.log ‖W.j‖) 0
        + Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) := by
  obtain ⟨m, hmpos, hm⟩ := exists_lower_bound_maxJ_pD
  set A : ℝ := (2 * Real.pi) ^ 12 * (4096 * Real.pi ^ 12) with hAdef
  have hApos : (0:ℝ) < A := by
    have := Real.pi_pos
    rw [hAdef]; positivity
  refine ⟨-(Real.log m + Real.log A), fun W hW => ?_⟩
  obtain ⟨τ, hjτ, hci⟩ := exists_curveArchInv_eq_petersson W hW
  have hpd : (0:ℝ) < peterssonDelta τ := peterssonDelta_pos τ
  have hkey : (2 * Real.pi) ^ 12 * curveArchInv W = A * peterssonDelta τ := by
    rw [hci, hAdef]; ring
  have hjn : ‖W.j‖ = ‖jFun τ‖ := by rw [hjτ]
  have hmaxpos : (0:ℝ) < max ‖jFun τ‖ 1 := lt_of_lt_of_le zero_lt_one (le_max_right _ _)
  rw [neg_neg, hkey, hjn, max_log_eq_log_max (norm_nonneg _),
    Real.log_mul (ne_of_gt hApos) (ne_of_gt hpd)]
  have hlm : Real.log m ≤ Real.log (max ‖jFun τ‖ 1 * peterssonDelta τ) :=
    Real.log_le_log hmpos (hm τ)
  rw [Real.log_mul (ne_of_gt hmaxpos) (ne_of_gt hpd)] at hlm
  linarith

/-- ★★★★★**埋め込みごとの形**——`archNorm` の言葉で。 -/
theorem logPlus_archNorm_lower {L : Type*} [Field L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ) (C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      -C ≤ max (Real.log ‖W.j‖) 0 + Real.log ((2 * Real.pi) ^ 12 * curveArchInv W)) :
    -C ≤ max (Real.log ‖σ E.j‖) 0
      + Real.log ((2 * Real.pi) ^ 12 * archNorm E σ) := by
  simpa [archNorm, WeierstrassCurve.map_j] using hC (E.map σ)

/-! ## ★出典の紐付け(`.src`) -/

def valueAtInfty_E4.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(E₄ のカスプ値は 1)",
    sectionId := "genell-prop-3-4" }

def exists_lower_bound_maxJ_pD.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ ——max(‖j‖,1)·‖Δ‖_Pet は下に有界)",
    sectionId := "genell-prop-3-4" }

def exists_logPlus_j_lower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ のアルキメデス側——下からの形)",
    sectionId := "genell-prop-3-4" }

def exists_logPlus_j_lower.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "EisensteinSeries.E_qExpansion_coeff_zero(正規化 Eisenstein 級数の定数項は 1)"
      (.inMathlib "EisensteinSeries.E_qExpansion_coeff_zero") 3,
    .citation "[mathlib]"
      "ModularFormClass.exp_decay_sub_atImInfty'(保型形式はカスプ値へ指数的に収束)"
      (.inMathlib "ModularFormClass.exp_decay_sub_atImInfty'") 3,
    .otherPaper "[Silv2]"
      ("Proposition 2.1——★原文が Prop 3.4 の根拠として引く。" ++
       "★★**その第 3 の ≲ のアルキメデス側が本ファイルである**") 5,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 557): §9-1000〜1001 の**逆向き**を取った。" ++
       "★上からは ‖j‖·pD^{1+ϵ} ≤ M(ϵ > 0 が要る)、" ++
       "下からは max(‖j‖,1)·pD ≥ m > 0(ϵ は要らない)。" ++
       "☆この非対称性が原文の 12(1+ϵ) の ϵ の出どころである。" ++
       "★★機構はカスプで E₄ → 1(正規化 Eisenstein 級数の定数項が 1)") 8 ]

end ABC3.Found.GenEll
