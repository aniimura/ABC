import ABC3.Found.GenEll.JSurjective

/-!
# GenEll 第 349 ブロック —— **★★★★★★周期束の共体積とアルキメデス不変量**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜ共体積が要るのか——界面の欠陥 #6

`Interface/GaloisRep/Reduction.lean` の `FaltingsHeightData` は `htFalt` に
**変数変換不変性**しか課していないので、第 329 の witness は `htFalt := deg∞/12` で
埋まってしまい、`prop_3_4` が恒等的に成り立つ形になった
(`Check/GaloisRep/HtFaltNotPinned.lean`)。

★★塞ぐには**アルキメデス素点での計量**が要る。真の Faltings 高さは

    12·d·ht^Falt(E) = Σ_p v_p(Δ_min)·log N(p) − Σ_{v|∞} log( (2π)^{12}·‖Δ‖_arch(E_v) )

で、`‖Δ‖_arch` に**周期束の共体積**が入る。★★★本ブロックはその共体積を建てる。

## ★★★★★★★核心——**`‖Δ_lat‖·covol⁶` は相似で変わらない**

| 量 | `Λ ↦ cΛ` での変わり方 |
|---|---|
| `Δ_lat = g₂³ − 27g₃²` | `|c|⁻¹²` 倍(第 333 `latticeDisc_scalePair`) |
| `covol` | `|c|²` 倍(本ブロック `covol_scalePair`) |
| **`‖Δ_lat‖·covol⁶`** | ★**不変**(`archInv_scalePair`) |

★★★★★これが「曲線に付く」量である——周期束は相似の分だけ選べるが、`archInv` は選び方に依らない。
★★これこそ Faltings 高さのアルキメデス局所因子である。

## ★★★★★モジュラーな姿——`Δ` の Petersson ノルム

正規化された束 `ℤ + τℤ` では `covol = Im τ` なので

    archInv(ℤ + τℤ) = 4096·π¹²·‖Δ(τ)‖·(Im τ)⁶

で、`‖Δ(τ)‖(Im τ)⁶` は **`Δ` の Petersson ノルム**である。
★★★これは `SL(2,ℤ)` 不変(`peterssonDelta_smul`)——重さ 12 の `|cτ+d|¹²` と
`Im(γτ) = Im τ/|cτ+d|²` の 6 乗がちょうど打ち消す。

★★★★★★**任意の周期束は Petersson ノルムに帰着する**(`exists_archInv_eq_petersson`)。

## ★★2026-08-26 の実測(mathlib の在庫)

★mathlib は `ZLattice.covolume`(Haar 測度による基本領域の体積)を持つ
(`Algebra/Module/ZLattice/Covolume.lean`)。
★★本ブロックは **2 次元の行列式そのもの**(`|Im(ω̄₁ω₂)|`)で定義した
——`PeriodPair` に `IsZLattice` を付ける配管を通さずに済み、
必要な**スケール則と正値性**は 3 行で出るからである。
★★★これは**逸脱**にあたるので記録する:`ZLattice.covolume` との一致
(`covolume_eq_det`)は将来必要になったら足す。下流が使うのは
`covol_scalePair`・`covol_pos`・`covol_tauPair` の 3 本だけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `covol` | ★★★★共体積 `|Im(ω̄₁ω₂)|` |
| `covol_pos` | ★★★正値(一次独立から) |
| `covol_scalePair` | ★★★★★**`covol(cΛ) = |c|²covol(Λ)`** |
| `archInv` | ★★★★★**`‖Δ_lat‖·covol⁶`** |
| `archInv_scalePair` | ★★★★★★**相似不変** |
| `covol_tauPair`・`archInv_tauPair` | ★★★★`covol = Im τ`、`archInv = 4096π¹²‖Δ‖(Im τ)⁶` |
| `peterssonDelta`・`peterssonDelta_smul` | ★★★★★**`SL(2,ℤ)` 不変** |
| `exists_archInv_eq_petersson` | ★★★★★★**任意の束は Petersson ノルムに帰着** |
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane ABC3.Found.GaloisRep

/-! ## ★★★★共体積 -/

/-- ★★★★**周期束の共体積** `|Im(ω̄₁ω₂)| = |ω₁.re·ω₂.im − ω₁.im·ω₂.re|`。

★★`ℂ ≅ ℝ²` での基本平行四辺形の面積である。 -/
noncomputable def covol (L : PeriodPair) : ℝ := |L.ω₁.re * L.ω₂.im - L.ω₁.im * L.ω₂.re|

/-- ★★★共体積は正——`ω₁, ω₂` は `ℝ` 上一次独立だから。 -/
theorem covol_pos (L : PeriodPair) : 0 < covol L := by
  rw [covol, abs_pos]
  intro h0
  refine ratio_im_ne_zero L ?_
  have h2 := omega₂_ne_zero L
  rw [Complex.div_im]
  have hns : Complex.normSq L.ω₂ ≠ 0 := by
    simpa [Complex.normSq_eq_zero] using h2
  field_simp
  linarith [h0]

/-- ★★★★★**`covol(cΛ) = |c|²·covol(Λ)`**——`c` 倍は `ℝ` 上の行列式 `|c|²` の写像。 -/
theorem covol_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    covol (scalePair L c hc) = Complex.normSq c * covol L := by
  have h : (scalePair L c hc).ω₁.re * (scalePair L c hc).ω₂.im
      - (scalePair L c hc).ω₁.im * (scalePair L c hc).ω₂.re
      = Complex.normSq c * (L.ω₁.re * L.ω₂.im - L.ω₁.im * L.ω₂.re) := by
    show (c * L.ω₁).re * (c * L.ω₂).im - (c * L.ω₁).im * (c * L.ω₂).re = _
    simp only [Complex.mul_re, Complex.mul_im, Complex.normSq_apply]
    ring
  rw [covol, covol, h, abs_mul, abs_of_nonneg (Complex.normSq_nonneg c)]

/-- ★周期の入れ替えで共体積は変わらない(行列式の符号だけが変わる)。 -/
theorem covol_swap (L : PeriodPair) :
    covol (⟨L.ω₂, L.ω₁, indep_swap L⟩ : PeriodPair) = covol L := by
  rw [covol, covol]
  show |L.ω₂.re * L.ω₁.im - L.ω₂.im * L.ω₁.re| = _
  rw [show L.ω₂.re * L.ω₁.im - L.ω₂.im * L.ω₁.re
      = -(L.ω₁.re * L.ω₂.im - L.ω₁.im * L.ω₂.re) by ring, abs_neg]

/-! ## ★★★★★★アルキメデス不変量 -/

/-- ★★★★★**アルキメデス不変量** `‖Δ_lat‖·covol⁶`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これが Faltings 高さのアルキメデス局所因子である。 -/
noncomputable def archInv (L : PeriodPair) : ℝ := ‖latticeDisc L‖ * covol L ^ 6

/-- ★★★★★★**`archInv` は相似で変わらない**——
`‖Δ‖` は `|c|⁻¹²` 倍、`covol⁶` は `|c|¹²` 倍でちょうど打ち消す。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★これが「周期束の選び方に依らない」ことの内容である。 -/
theorem archInv_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    archInv (scalePair L c hc) = archInv L := by
  have hn : Complex.normSq c = ‖c‖ ^ 2 := (Complex.normSq_eq_norm_sq c)
  have hcn : ‖c‖ ≠ 0 := norm_ne_zero_iff.2 hc
  rw [archInv, archInv, latticeDisc_scalePair L c hc, covol_scalePair L c hc, hn]
  rw [norm_mul, norm_inv, norm_pow]
  field_simp

theorem archInv_swap (L : PeriodPair) :
    archInv (⟨L.ω₂, L.ω₁, indep_swap L⟩ : PeriodPair) = archInv L := by
  rw [archInv, archInv, covol_swap, latticeDisc_congr (swap_lattice L (indep_swap L))]

/-- ★★★`archInv` は正。 -/
theorem archInv_pos (L : PeriodPair) : 0 < archInv L :=
  mul_pos (norm_pos_iff.2 (latticeDisc_ne_zero L)) (pow_pos (covol_pos L) 6)

/-! ## ★★★★正規化された束 -/

/-- ★★★**`covol(ℤ + τℤ) = Im τ`**。 -/
theorem covol_tauPair (τ : UpperHalfPlane) : covol (tauPair τ) = τ.im := by
  have h : (tauPair τ).ω₁.re * (tauPair τ).ω₂.im - (tauPair τ).ω₁.im * (tauPair τ).ω₂.re
      = -(τ : ℂ).im := by
    show (τ : ℂ).re * (1 : ℂ).im - (τ : ℂ).im * (1 : ℂ).re = _
    simp
  rw [covol, h, abs_neg, abs_of_pos τ.2]
  rfl

/-- ★★★★★**`archInv(ℤ + τℤ) = 4096π¹²·‖Δ(τ)‖·(Im τ)⁶`**。 -/
theorem archInv_tauPair (τ : UpperHalfPlane) :
    archInv (tauPair τ) = 4096 * Real.pi ^ 12 * ‖ModularForm.discriminant τ‖ * τ.im ^ 6 := by
  rw [archInv, covol_tauPair, latticeDisc_tauPair, norm_mul, norm_mul]
  simp only [Complex.norm_ofNat, norm_pow, Complex.norm_real, Real.norm_eq_abs,
    abs_of_nonneg Real.pi_nonneg]

/-! ## ★★★★★`Δ` の Petersson ノルム -/

/-- ★★★★**`Δ` の Petersson ノルム** `‖Δ(τ)‖·(Im τ)⁶`。 -/
noncomputable def peterssonDelta (τ : UpperHalfPlane) : ℝ :=
  ‖ModularForm.discriminant τ‖ * τ.im ^ 6

theorem peterssonDelta_pos (τ : UpperHalfPlane) : 0 < peterssonDelta τ :=
  mul_pos (norm_pos_iff.2 (ModularForm.discriminant_ne_zero τ)) (pow_pos τ.2 6)

/-- ★★★★★**Petersson ノルムは `𝒮ℒ` 不変**——重さ 12 の `|cτ+d|¹²` と
`Im(γτ) = Im τ/|cτ+d|²` の 6 乗がちょうど打ち消す。 -/
theorem peterssonDelta_smul_GL {γ : GL (Fin 2) ℝ} (hγ : γ ∈ 𝒮ℒ) (τ : UpperHalfPlane) :
    peterssonDelta (γ • τ) = peterssonDelta τ := by
  have hd : ‖UpperHalfPlane.denom γ (τ : ℂ)‖ ≠ 0 :=
    norm_ne_zero_iff.2 (UpperHalfPlane.denom_ne_zero γ τ)
  have hdet : ((Matrix.GeneralLinearGroup.det γ : ℝˣ) : ℝ) = 1 := by
    rw [Subgroup.HasDetOne.det_eq (Γ := 𝒮ℒ) hγ]; rfl
  have h12 := SlashInvariantForm.slash_action_eqn'' CuspForm.discriminant hγ τ
  rw [CuspForm.coe_discriminant] at h12
  have him : (γ • τ).im = τ.im / Complex.normSq (UpperHalfPlane.denom γ (τ : ℂ)) := by
    rw [UpperHalfPlane.im_smul_eq_div_normSq, hdet]
    simp
  have hns : Complex.normSq (UpperHalfPlane.denom γ (τ : ℂ))
      = ‖UpperHalfPlane.denom γ (τ : ℂ)‖ ^ 2 := Complex.normSq_eq_norm_sq _
  rw [peterssonDelta, peterssonDelta, h12, him, hns, norm_mul, norm_zpow]
  rw [show ((12 : ℤ)) = ((12 : ℕ) : ℤ) from rfl, zpow_natCast]
  field_simp

/-- ★★★★★**Petersson ノルムは `SL(2,ℤ)` 不変**。 -/
theorem peterssonDelta_smul (g : Matrix.SpecialLinearGroup (Fin 2) ℤ) (τ : UpperHalfPlane) :
    peterssonDelta (g • τ) = peterssonDelta τ := by
  rw [show (g • τ) = (Matrix.SpecialLinearGroup.mapGL ℝ g) • τ from
    (MulAction.compHom_smul_def (Matrix.SpecialLinearGroup.mapGL ℝ) g τ)]
  exact peterssonDelta_smul_GL ⟨g, rfl⟩ τ

/-! ## ★★★★★★束から Petersson ノルムへ -/

theorem archInv_eq_petersson_of_im_pos (M : PeriodPair) (h : 0 < (M.ω₁ / M.ω₂).im) :
    archInv M = 4096 * Real.pi ^ 12 * peterssonDelta ⟨M.ω₁ / M.ω₂, h⟩ := by
  have hE := eq_scalePair_tauPair M h
  calc archInv M
      = archInv (scalePair (tauPair ⟨M.ω₁ / M.ω₂, h⟩) M.ω₂ (omega₂_ne_zero M)) := by rw [hE]
    _ = archInv (tauPair ⟨M.ω₁ / M.ω₂, h⟩) := archInv_scalePair _ _ _
    _ = 4096 * Real.pi ^ 12 * peterssonDelta ⟨M.ω₁ / M.ω₂, h⟩ := by
        rw [archInv_tauPair, peterssonDelta]; ring

/-- ★★★★★★**任意の周期束のアルキメデス不変量は `Δ` の Petersson ノルムである**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★向きが負のときは `ω₁ ↔ ω₂` を入れ替える(束も `archInv` も変わらない)。 -/
theorem exists_archInv_eq_petersson (L : PeriodPair) :
    ∃ τ : UpperHalfPlane, archInv L = 4096 * Real.pi ^ 12 * peterssonDelta τ := by
  rcases lt_or_gt_of_ne (ratio_im_ne_zero L) with hneg | hpos
  · have hind := indep_swap L
    have him : 0 < ((⟨L.ω₂, L.ω₁, hind⟩ : PeriodPair).ω₁
        / (⟨L.ω₂, L.ω₁, hind⟩ : PeriodPair).ω₂).im := by
      show 0 < (L.ω₂ / L.ω₁).im
      rw [show L.ω₂ / L.ω₁ = (L.ω₁ / L.ω₂)⁻¹ from (inv_div _ _).symm]
      have hz : L.ω₁ / L.ω₂ ≠ 0 := fun h0 => by simp [h0] at hneg
      rw [Complex.inv_im]
      exact div_pos (by linarith) (Complex.normSq_pos.2 hz)
    exact ⟨_, by rw [← archInv_swap L, archInv_eq_petersson_of_im_pos _ him]⟩
  · exact ⟨_, archInv_eq_petersson_of_im_pos L hpos⟩

/-! ## ★出典の紐付け(`.src`) -/

def covol.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def archInv_scalePair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_archInv_eq_petersson.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
