/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.Ring.Finset

/-!
# Elementary —— `[GenEll] Lemma 3.6` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real

/-! ## Lemma 3.6 —— An Elementary Estimate -/

/-- **[GenEll] Lemma 3.6**(An Elementary Estimate)。

原文 (GenEll p.17):
> Lemma 3.6. (An Elementary Estimate) Let

`ϵ ∈ ℝ_{>0}` に対し**定数 `C₀ ∈ ℝ_{>0}` が存在して**、`y ≥ 1` かつ `x ≥ C₀·y^{1+ϵ}` なる
全ての `x, y ∈ ℝ` について **`x ≥ y·log(x)`**。

## ★定数を明示的に構成した

原文 (GenEll p.18):
> This follows immediately from the well-known elementary fact that x

——続きは `x^{1/(1+ϵ)}·log(x)/x = log(x)·x^{−ϵ/(1+ϵ)} → 0`(`x → ∞`)。
★原文は**極限**で片づけているが、**それでは定数が出ない**。

★我々は `δ ≝ ϵ/(1+ϵ)` と置き、**`C₀ ≝ (1/δ)^{1+ϵ}` を明示的に取った**。
道筋は次の 3 つの不等式である:

1. `log(x) ≤ x^δ/δ` —— `log(x^δ) ≤ x^δ − 1`(`Real.log_le_sub_one_of_pos`)から。
2. `y ≤ (x/C₀)^{1/(1+ϵ)}` —— 仮定を `1/(1+ϵ)` 乗して。
3. `(x/C₀)^{1/(1+ϵ)}·x^δ/δ = x` —— `C₀^{1/(1+ϵ)} = 1/δ` と
   `1/(1+ϵ) + δ = 1` から**きっかり `x` になる**。

★★**3 の等号がきっかり成り立つように `C₀` を選んである**。
`δ < 1` ゆえ `C₀ > 1`、したがって `x > 1` で `log x > 0` も従う——
この符号が無いと 2 の両辺に `log x` を掛けられない。 -/
theorem lemma_3_6 (eps : ℝ) (heps : 0 < eps) :
    ∃ C₀ : ℝ, 0 < C₀ ∧ ∀ x y : ℝ, 1 ≤ y → C₀ * y ^ (1 + eps) ≤ x → y * Real.log x ≤ x := by
  have h1e : (0 : ℝ) < 1 + eps := by linarith
  set δ : ℝ := eps / (1 + eps) with hδdef
  have hδ : 0 < δ := div_pos heps h1e
  have hδ1 : δ < 1 := by
    rw [hδdef, div_lt_one h1e]; linarith
  have hinvδ : 1 < 1 / δ := by
    rw [lt_div_iff₀ hδ]; linarith
  -- ★`C₀ ≝ (1/δ)^{1+ϵ}`
  set C₀ : ℝ := (1 / δ) ^ (1 + eps) with hC₀def
  have hC₀1 : 1 < C₀ := by
    have := Real.rpow_lt_rpow (by norm_num) hinvδ h1e
    rwa [Real.one_rpow] at this
  have hC₀0 : 0 < C₀ := lt_trans one_pos hC₀1
  have hmul : (1 + eps) * (1 / (1 + eps)) = 1 := by field_simp
  -- ★`C₀^{1/(1+ϵ)} = 1/δ`
  have hC₀root : C₀ ^ (1 / (1 + eps)) = 1 / δ := by
    rw [hC₀def, ← Real.rpow_mul (by positivity), hmul, Real.rpow_one]
  refine ⟨C₀, hC₀0, ?_⟩
  intro x y hy hx
  have hy0 : (0 : ℝ) < y := lt_of_lt_of_le one_pos hy
  -- `1 ≤ y^{1+ϵ}`
  have hyp1 : (1 : ℝ) ≤ y ^ (1 + eps) := by
    have := Real.rpow_le_rpow (by norm_num) hy h1e.le
    rwa [Real.one_rpow] at this
  -- `1 < C₀ ≤ C₀·y^{1+ϵ} ≤ x`
  have hx1 : 1 < x := by
    have : C₀ ≤ C₀ * y ^ (1 + eps) := le_mul_of_one_le_right hC₀0.le hyp1
    linarith
  have hx0 : (0 : ℝ) < x := lt_trans one_pos hx1
  have hlogpos : 0 ≤ Real.log x := Real.log_nonneg hx1.le
  -- ★① `log(x) ≤ x^δ/δ`
  have hlog : Real.log x ≤ x ^ δ / δ := by
    have hxδ : (0 : ℝ) < x ^ δ := Real.rpow_pos_of_pos hx0 δ
    have h := Real.log_le_sub_one_of_pos hxδ
    rw [Real.log_rpow hx0] at h
    rw [le_div_iff₀ hδ]
    nlinarith
  -- ★② `y ≤ (x/C₀)^{1/(1+ϵ)}`
  have hyle : y ≤ (x / C₀) ^ (1 / (1 + eps)) := by
    have hdiv : y ^ (1 + eps) ≤ x / C₀ := by
      rw [le_div_iff₀ hC₀0]; linarith [hx]
    have := Real.rpow_le_rpow (by positivity) hdiv (by positivity : (0:ℝ) ≤ 1 / (1 + eps))
    rwa [← Real.rpow_mul hy0.le, hmul, Real.rpow_one] at this
  -- ★③ `(x/C₀)^{1/(1+ϵ)} · x^δ/δ = x`
  have hkey : (x / C₀) ^ (1 / (1 + eps)) * (x ^ δ / δ) = x := by
    rw [Real.div_rpow hx0.le hC₀0.le, hC₀root]
    have hsum : 1 / (1 + eps) + δ = 1 := by
      rw [hδdef]; field_simp
    have : x ^ (1 / (1 + eps)) * x ^ δ = x := by
      rw [← Real.rpow_add hx0, hsum, Real.rpow_one]
    field_simp
    linarith [this]
  calc y * Real.log x
      ≤ (x / C₀) ^ (1 / (1 + eps)) * (x ^ δ / δ) := by
        exact mul_le_mul hyle hlog hlogpos (by positivity)
    _ = x := hkey

/-! ## ★出典の紐付け(`.src`) -/

def lemma_3_6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17, item := "Lemma 3.6",
    sectionId := "genell-lemma-3-6" }

end ABC3.Found.GenEll
