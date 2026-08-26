import ABC3.Found.GenEll.LatticeDiscNonzero
import Mathlib.Analysis.Complex.CauchyIntegral

/-!
# GenEll 第 347 ブロック —— **★★★★★★モジュラー `j` 関数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜ `j` 関数を建てるのか——一意化 (iii)

第 346 で (i)(判別式の非消失)が閉じ、`j(ℤ + τℤ) = E₄(τ)³/Δ(τ)` まで出た。
★スケルトン `Skeleton/GenEll/Uniformization.lean` の `exists_periodPair`
(任意の `E/ℂ` が格子曲線と同型)は、mathlib の
**`WeierstrassCurve.exists_variableChange_of_j_eq`**(分離閉体上で `j` が等しければ同型)
と合わせると **`j : ℍ → ℂ` の全射性**だけに帰着する。

★★★**mathlib にモジュラー `j` 関数は無い**(2026-08-26 実測)。本ブロックで建てる。

## ★★★★★★全射性の筋——**開かつ閉**

古典的には Sturm 境界や valence 公式(留数計算)を使うが、
★★★**開写像定理 + 基本領域の切り詰めのコンパクト性**で足りる:

| 段 | 内容 |
|---|---|
| 1 | `j` は上半平面で正則で**定数でない** ⟹ 像は**開**(開写像定理) |
| 2 | `j` は `SL(2,ℤ)` 不変 ⟹ 任意の点は基本領域 `𝒟` に移せる |
| 3 | `‖j‖ → ∞`(カスプ)⟹ 有界な値は `𝒟 ∩ {im ≤ M}` から来る |
| 4 | `𝒟 ∩ {im ≤ M}` は**コンパクト**(mathlib `isCompact_truncatedFundamentalDomain`)⟹ 像は**閉** |
| 5 | `ℂ` は連結 ⟹ 像は `ℂ` 全体 |

★★★★★**valence 公式(留数計算)を一切使わない**。段 4 の在庫が効いた。

## ★★★★★2026-08-26 の実測(mathlib の在庫)

| 段 | mathlib |
|---|---|
| 基本領域 `𝒟`・`∃ g, g • z ∈ 𝒟` | ✅ `ModularGroup.fd`・`exists_smul_mem_fd` |
| **切り詰めた基本領域のコンパクト性** | ✅ `isCompact_truncatedFundamentalDomain` |
| 開写像定理(連結開集合版) | ✅ `AnalyticOnNhd.is_constant_or_isOpen` |
| `f(γ•z) = (cz+d)^k f(z)` | ✅ `SlashInvariantForm.slash_action_eqn''` |
| `Δ` はカスプ形式(`→ 0`) | ✅ `CuspFormClass.zero_at_infty` |
| モジュラー形式の**カスプでの極限** | ★無い(`cuspFunction_apply_zero` から 5 行で出る) |
| **`j` 関数そのもの** | ★無い |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `jFun` | ★★★★**`j = E₄³/Δ`** |
| `jFun_smul` | ★★★★`SL(2,ℤ)` 不変性 |
| `jFun_eq_latticeCurve_j` | ★★★格子曲線の `j` 不変量と一致 |
| `jC`・`analyticOnNhd_jC` | ★★★★★上半平面での正則性 |
| `tendsto_atImInfty_coeff_zero` | ★★★★★カスプでの極限は定数項 |
| `tendsto_E₄_atImInfty`・`tendsto_discriminant_atImInfty` | ★★★★`E₄ → 1`・`Δ → 0` |
| `tendsto_norm_jFun_atImInfty` | ★★★★★**`‖j‖ → ∞`** |
| `jC_not_constant` | ★★★★**`j` は定数でない** |
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane EisensteinSeries ABC3.Found.GaloisRep
open Filter Topology

/-! ## ★★★★モジュラー `j` 関数 -/

/-- ★★★★★★**モジュラー `j` 関数** `j = E₄³/Δ`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
noncomputable def jFun (τ : ℍ) : ℂ := ModularForm.E₄ τ ^ 3 / ModularForm.discriminant τ

/-- ★★★★`𝒮ℒ` 不変性——重さ 4 の 3 乗と重さ 12 で `(cτ+d)¹²` が消える。 -/
theorem jFun_smul_GL {γ : GL (Fin 2) ℝ} (hγ : γ ∈ 𝒮ℒ) (τ : ℍ) : jFun (γ • τ) = jFun τ := by
  have hd := UpperHalfPlane.denom_ne_zero γ τ
  have hΔ := ModularForm.discriminant_ne_zero τ
  have h4 := SlashInvariantForm.slash_action_eqn'' (ModularForm.E₄) hγ τ
  have h12 := SlashInvariantForm.slash_action_eqn'' CuspForm.discriminant hγ τ
  rw [CuspForm.coe_discriminant] at h12
  rw [jFun, jFun, h4, h12]
  field_simp

/-- ★★★★**`j` は `SL(2,ℤ)` 不変**。 -/
theorem jFun_smul (g : Matrix.SpecialLinearGroup (Fin 2) ℤ) (τ : ℍ) : jFun (g • τ) = jFun τ := by
  rw [show (g • τ) = (Matrix.SpecialLinearGroup.mapGL ℝ g) • τ from
    (MulAction.compHom_smul_def (Matrix.SpecialLinearGroup.mapGL ℝ) g τ)]
  exact jFun_smul_GL ⟨g, rfl⟩ τ

/-- ★★★**`j` は格子 `ℤ + τℤ` の曲線の `j` 不変量**(第 346)。 -/
theorem jFun_eq_latticeCurve_j (τ : ℍ) : jFun τ = (latticeCurve (tauPair τ)).j :=
  (latticeCurve_j_tauPair τ).symm

/-! ## ★★★★★上半平面での正則性 -/

/-- ★★★★`j` を `ℂ` 上に持ち上げたもの(上半平面の外では junk 値)。 -/
noncomputable def jC : ℂ → ℂ := jFun ∘ UpperHalfPlane.ofComplex

theorem jC_apply {z : ℂ} (hz : 0 < z.im) : jC z = jFun ⟨z, hz⟩ := by
  simp [jC, UpperHalfPlane.ofComplex_apply_of_im_pos hz]

theorem jC_comp_coe : jC ∘ UpperHalfPlane.coe = jFun := by
  funext σ
  exact jC_apply σ.2

/-- ★★★★★**`j` は上半平面で正則**——`E₄` と `Δ` の商で `Δ ≠ 0`。 -/
theorem analyticOnNhd_jC : AnalyticOnNhd ℂ jC UpperHalfPlane.upperHalfPlaneSet := by
  refine DifferentiableOn.analyticOnNhd ?_ UpperHalfPlane.isOpen_upperHalfPlaneSet
  intro z hz
  have hz' : 0 < z.im := hz
  have hE : DifferentiableAt ℂ (⇑ModularForm.E₄ ∘ UpperHalfPlane.ofComplex) z :=
    UpperHalfPlane.mdifferentiableAt_iff.mp (ModularFormClass.holo ModularForm.E₄ ⟨z, hz'⟩)
  have hD : DifferentiableAt ℂ (ModularForm.discriminant ∘ UpperHalfPlane.ofComplex) z := by
    have := UpperHalfPlane.mdifferentiableAt_iff.mp
      (ModularFormClass.holo CuspForm.discriminant ⟨z, hz'⟩)
    simpa [CuspForm.coe_discriminant] using this
  have hne : (ModularForm.discriminant ∘ UpperHalfPlane.ofComplex) z ≠ 0 :=
    ModularForm.discriminant_ne_zero _
  exact (((hE.pow 3).div hD hne)).differentiableWithinAt

/-- ★★★`j` は連続。 -/
theorem continuous_jFun : Continuous jFun := by
  rw [← jC_comp_coe, continuous_iff_continuousAt]
  intro τ
  exact ((analyticOnNhd_jC (τ : ℂ) τ.2).continuousAt).comp
    UpperHalfPlane.continuous_coe.continuousAt

/-! ## ★★★★★カスプでの極限 -/

/-- ★★★★★**レベル 1 のモジュラー形式はカスプで q 展開の定数項に収束する**。

★mathlib は `cuspFunction_apply_zero`(値の同定)までで、極限の形は持っていない。 -/
theorem tendsto_atImInfty_coeff_zero {k : ℤ} {F : Type*} [FunLike F ℍ ℂ]
    [ModularFormClass F 𝒮ℒ k] (f : F) :
    Tendsto (f : ℍ → ℂ) atImInfty (𝓝 ((PowerSeries.coeff 0) (qExpansion 1 (f : ℍ → ℂ)))) := by
  have hper := SlashInvariantFormClass.periodic_comp_ofComplex f
    (h := (1:ℝ)) one_mem_strictPeriods_SL
  have hana := ModularFormClass.analyticAt_cuspFunction_zero f (h := (1:ℝ))
    one_pos one_mem_strictPeriods_SL
  have hc : (PowerSeries.coeff 0) (qExpansion 1 (f : ℍ → ℂ)) = cuspFunction 1 (f : ℍ → ℂ) 0 := by
    rw [qExpansion_coeff]
    simp
  rw [hc]
  have heq : (cuspFunction 1 (f : ℍ → ℂ) ∘ fun τ : ℍ => Function.Periodic.qParam 1 (τ : ℂ))
      = (f : ℍ → ℂ) := by
    funext τ
    simpa using eq_cuspFunction τ one_ne_zero hper
  simpa [heq] using hana.continuousAt.tendsto.comp (qParam_tendsto_atImInfty (h := (1:ℝ)) one_pos)

/-- ★★★★**`E₄ → 1`**(カスプ)。 -/
theorem tendsto_E₄_atImInfty : Tendsto (ModularForm.E₄ : ℍ → ℂ) atImInfty (𝓝 1) := by
  have h := tendsto_atImInfty_coeff_zero (ModularForm.E₄)
  rwa [EisensteinSeries.E_qExpansion_coeff_zero (by norm_num) (by decide)] at h

/-- ★★★★**`Δ → 0`**(カスプ)——`Δ` はカスプ形式。 -/
theorem tendsto_discriminant_atImInfty : Tendsto ModularForm.discriminant atImInfty (𝓝 0) := by
  have h := CuspFormClass.zero_at_infty (CuspForm.discriminant)
  rw [CuspForm.coe_discriminant] at h
  exact h

/-- ★★★★★**`‖j‖ → ∞`**(カスプ)——`‖E₄‖³ → 1` を `‖Δ‖⁻¹ → ∞` に掛ける。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★カスプでは `E₄ → 1`・`Δ → 0` なので商は発散する。 -/
theorem tendsto_norm_jFun_atImInfty : Tendsto (fun τ : ℍ => ‖jFun τ‖) atImInfty atTop := by
  have hE : Tendsto (fun τ : ℍ => ‖ModularForm.E₄ τ‖ ^ 3) atImInfty (𝓝 1) := by
    have h1 : Tendsto (fun τ : ℍ => ‖ModularForm.E₄ τ‖) atImInfty (𝓝 1) := by
      simpa using (tendsto_E₄_atImInfty).norm
    simpa using h1.pow 3
  have hD0 : Tendsto (fun τ : ℍ => ‖ModularForm.discriminant τ‖) atImInfty (𝓝[>] 0) := by
    refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _
      (by simpa using tendsto_discriminant_atImInfty.norm) ?_
    filter_upwards with τ
    exact norm_pos_iff.2 (ModularForm.discriminant_ne_zero τ)
  have hmul : Tendsto (fun τ : ℍ => ‖ModularForm.E₄ τ‖ ^ 3 * ‖ModularForm.discriminant τ‖⁻¹)
      atImInfty atTop :=
    Filter.Tendsto.pos_mul_atTop one_pos hE hD0.inv_tendsto_nhdsGT_zero
  refine hmul.congr (fun τ => ?_)
  rw [jFun, norm_div, norm_pow, div_eq_mul_inv]

/-- ★★★★**`j` は定数でない**——カスプで発散するから。 -/
theorem jC_not_constant : ¬ (∃ w, ∀ z ∈ UpperHalfPlane.upperHalfPlaneSet, jC z = w) := by
  rintro ⟨w, hw⟩
  have hconst : ∀ τ : ℍ, ‖jFun τ‖ = ‖w‖ := by
    intro τ
    rw [← jC_apply τ.2, hw (τ : ℂ) τ.2]
  obtain ⟨τ, hτ⟩ := (tendsto_norm_jFun_atImInfty.eventually_ge_atTop (‖w‖ + 1)).exists
  rw [hconst τ] at hτ
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def jFun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def tendsto_norm_jFun_atImInfty.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def tendsto_atImInfty_coeff_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
