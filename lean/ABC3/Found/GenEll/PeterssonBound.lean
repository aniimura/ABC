import ABC3.Found.GenEll.CurveArchInv

/-!
# GenEll 第 355 ブロック —— **★★★★★★アルキメデス因子は有界**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★なぜ有界性が要るのか——`prop_3_4` を保つため

第 354 で `archNorm` が建った。★これを `ht^Falt` の界面条件に入れると、
witness は `ht^Falt := deg∞/12 − A`(`A` はアルキメデス項)という形になる。
★★このとき `prop_3_4`(`deg∞/(12(1+ε)) ≤ ht^Falt + C`)を保つには、
**`A` が普遍定数で上から抑えられる**ことが要る。

★★★`A` は `log((2π)¹²·archNorm)` の平均なので、**`archNorm` の上界**があればよい。

## ★★★★★★段取り——**カスプでの減衰 + 基本領域の切り詰めのコンパクト性**

`archNorm = 4096π¹²·‖Δ(τ)‖(Im τ)⁶`(第 354)なので、
`peterssonDelta τ := ‖Δ(τ)‖(Im τ)⁶` の有界性に帰着する。

| 段 | 内容 |
|---|---|
| 1 | `Δ` はカスプ形式なので `Δ =O[atImInfty] exp(−2π·Im τ)`(mathlib) |
| 2 | `exp(−2πy)·y⁶ → 0`(`Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero` を `x = 2πy` で) |
| 3 | ⟹ `peterssonDelta → 0`(カスプ)、したがって `Im τ` が大きい所では `< 1` |
| 4 | 残りは `𝒟 ∩ {Im ≤ y₁}` で、これは**コンパクト**なので最大値を取る |
| 5 | 任意の `τ` は `SL(2,ℤ)` で `𝒟` に移せ、`peterssonDelta` は不変(第 349) |

★★★★★第 348(`j` の全射性)で使ったのと**同じ 2 つの在庫**——
指数減衰と `isCompact_truncatedFundamentalDomain`——がここでも効いた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tendsto_exp_mul_pow` | ★`exp(−2πy)y⁶ → 0` |
| `continuous_discriminant`・`continuous_peterssonDelta` | ★連続性 |
| `tendsto_peterssonDelta_atImInfty` | ★★★★カスプで 0 |
| `exists_bound_peterssonDelta` | ★★★★★★**Petersson ノルムは有界** |
| `exists_bound_curveArchInv`・`archNorm_le` | ★★★★★★**アルキメデス因子は有界** |
-/

namespace ABC3.Found.GenEll

open Complex Real ModularForm MatrixGroups UpperHalfPlane Filter Topology ABC3.Found.GaloisRep

/-! ## ★指数減衰が多項式に勝つ -/

theorem tendsto_exp_mul_pow : Tendsto (fun y : ℝ => Real.exp (-2 * Real.pi * y / 1) * y ^ 6)
    atTop (𝓝 0) := by
  have hpi : (0:ℝ) < 2 * Real.pi := by positivity
  have hbase := Real.tendsto_pow_mul_exp_neg_atTop_nhds_zero 6
  have hcomp : Tendsto (fun y : ℝ => (2 * Real.pi * y) ^ 6 * Real.exp (-(2 * Real.pi * y)))
      atTop (𝓝 0) := hbase.comp (Filter.Tendsto.const_mul_atTop hpi tendsto_id)
  have h := hcomp.div_const ((2 * Real.pi) ^ 6)
  rw [zero_div] at h
  refine h.congr (fun y => ?_)
  have h6 : ((2 * Real.pi) ^ 6 : ℝ) ≠ 0 := by positivity
  field_simp

/-! ## ★連続性 -/

theorem continuous_discriminant : Continuous ModularForm.discriminant := by
  have h := ModularFormClass.holo (CuspForm.discriminant)
  rw [CuspForm.coe_discriminant] at h
  exact h.continuous

theorem continuous_peterssonDelta : Continuous peterssonDelta :=
  continuous_discriminant.norm.mul (UpperHalfPlane.continuous_im.pow 6)

/-! ## ★★★★カスプでの減衰 -/

/-- ★★★★**Petersson ノルムはカスプで 0 に収束する**——`Δ` の指数減衰が `(Im τ)⁶` に勝つ。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Δ` がカスプ形式であることが効く。 -/
theorem tendsto_peterssonDelta_atImInfty : Tendsto peterssonDelta atImInfty (𝓝 0) := by
  have hO := CuspFormClass.exp_decay_atImInfty (CuspForm.discriminant) (h := (1:ℝ))
    one_pos one_mem_strictPeriods_SL
  rw [CuspForm.coe_discriminant] at hO
  have hO2 : peterssonDelta
      =O[atImInfty] (fun τ : ℍ => Real.exp (-2 * Real.pi * τ.im / 1) * τ.im ^ 6) :=
    hO.norm_left.mul (Asymptotics.isBigO_refl (fun τ : ℍ => τ.im ^ 6) atImInfty)
  refine hO2.trans_tendsto ?_
  exact tendsto_exp_mul_pow.comp tendsto_comap

/-! ## ★★★★★★有界性 -/

/-- ★★★★★★**Petersson ノルムは有界**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★カスプの外は切り詰めた基本領域(コンパクト)に収まり、`SL(2,ℤ)` 不変性で全体に及ぶ。 -/
theorem exists_bound_peterssonDelta : ∃ M : ℝ, ∀ τ : ℍ, peterssonDelta τ ≤ M := by
  obtain ⟨y₀, hy₀⟩ := (UpperHalfPlane.atImInfty_mem _).1
    (tendsto_peterssonDelta_atImInfty.eventually_lt_const (zero_lt_one))
  set y₁ := max y₀ 1 with hy₁
  have hcpt := ModularGroup.isCompact_truncatedFundamentalDomain y₁
  have hne : (ModularGroup.truncatedFundamentalDomain y₁).Nonempty := by
    refine ⟨UpperHalfPlane.I, ModularGroup.I_mem_fd, ?_⟩
    simp [hy₁]
  obtain ⟨x₀, hx₀mem, hx₀max⟩ :=
    hcpt.exists_isMaxOn hne (continuous_peterssonDelta.continuousOn)
  refine ⟨max (peterssonDelta x₀) 1, fun τ => ?_⟩
  obtain ⟨g, hg⟩ := ModularGroup.exists_smul_mem_fd τ
  have hinv : peterssonDelta (g • τ) = peterssonDelta τ := peterssonDelta_smul g τ
  rcases le_or_gt ((g • τ).im) y₁ with hle | hgt
  · have hm : peterssonDelta (g • τ) ≤ peterssonDelta x₀ := hx₀max ⟨hg, hle⟩
    rw [hinv] at hm
    exact le_trans hm (le_max_left _ _)
  · have h1 : peterssonDelta (g • τ) < 1 :=
      hy₀ (g • τ) (le_of_lt (lt_of_le_of_lt (le_max_left y₀ 1) hgt))
    rw [hinv] at h1
    exact le_trans h1.le (le_max_right _ _)

/-- ★★★★★★**`curveArchInv` は有界**——`prop_3_4` を保つのに要る。 -/
theorem exists_bound_curveArchInv :
    ∃ M : ℝ, 0 < M ∧ ∀ W : WeierstrassCurve ℂ, W.IsElliptic → curveArchInv W ≤ M := by
  obtain ⟨M₀, hM₀⟩ := exists_bound_peterssonDelta
  refine ⟨max (4096 * Real.pi ^ 12 * M₀) 1, lt_of_lt_of_le zero_lt_one (le_max_right _ _), ?_⟩
  intro W hW
  obtain ⟨τ, _, hτ⟩ := exists_curveArchInv_eq_petersson W hW
  rw [hτ]
  refine le_trans ?_ (le_max_left _ _)
  have hpos : (0:ℝ) ≤ 4096 * Real.pi ^ 12 := by positivity
  exact mul_le_mul_of_nonneg_left (hM₀ τ) hpos

theorem archNorm_le {L : Type*} [Field L] (M : ℝ)
    (hM : ∀ W : WeierstrassCurve ℂ, W.IsElliptic → curveArchInv W ≤ M)
    (E : WeierstrassCurve L) [E.IsElliptic] (σ : L →+* ℂ) : archNorm E σ ≤ M :=
  hM (E.map σ) inferInstance

/-! ## ★出典の紐付け(`.src`) -/

def exists_bound_peterssonDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def tendsto_peterssonDelta_atImInfty.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_bound_curveArchInv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
