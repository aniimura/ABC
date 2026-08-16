import ABC3.Meta.Claim
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Analysis.SpecialFunctions.Pow.Real

/-!
# [NCBelyi] Lemma 2.1 —— Belyi 写像の分離性(計算核)(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.2。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.2):
> Lemma 2.1. (Separating Properties of Belyi Maps) Let C ∈ R be such

## ★★本補題は [NCBelyi] の計算の中心である

`Theorem 2.5`(= [GenEll] `Theorem 2.1` が引く非臨界 Belyi 写像の存在)は
`Lemma 2.2` の `|S|` についての帰納法に帰着し、その各段が本補題である。

## ★本ファイルが取るもの

`f(x) ≝ x^m·(x−1)^n`(`m, n ≥ 1`)について:

- ★**導関数の恒等式** `f′(x) = x^{m−1}·(x−1)^{n−1}·{(m+n)x − m}`(原文の (b) の根拠)
- ★★**性質 (e)** `f(β)/f(α) ≥ (β/α)² ≥ β/α`(`1 < α ≤ β`)——(d) の両場合を貫く
- `x > 1` での**単調性**
- `[0,1]` 上で `|f(x)| ≤ 1`、`β ≥ 2` で `f(β) ≥ β`(原文の (c) の根拠)

## ★★★原文は単調性を導関数から出すが、要らなかった

原文は「`f′(x) > 0` for real `x > 1`、ゆえに `f` は単調増加」と進む。
★**性質 (e) から直接出る**——`β ≥ α > 1` なら `f(β)/f(α) ≥ (β/α)² ≥ 1`。
微分も中間値定理も要らない。
★「正面から要ると思ったものが要らなかった」の型がまた出た。
(導関数の恒等式そのものは (b)——臨界点の位置——に要るので別に取る。)

## ★★pdftotext が落とすもの(2026-08-17 目視で発見)

- 条件 (iii) の `α ≠ 0, r, 1, ∞` の **`≠` が `=` に化ける**
- **`f′(x)` のプライム `′` が消える**(`f (x) = 0` に見える)。条件 (d) の `≠` も同様

★逐語照合は pdftotext の出力に対して行うので、
`1_Structured` では **verbatim を pdftotext の綴りに合わせ、reading に真の文を書いた**。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★導関数の恒等式 -/

/-- ★**`f′(x) = x^{m−1}·(x−1)^{n−1}·{(m+n)x − m}`**(`m = a+1`, `n = b+1`)。

原文 (NCBelyi p.2):
> mediately from the fact that:

★`m, n ≥ 1` を `a+1, b+1` と書いて自然数の引き算を避けた。 -/
theorem derivative_belyi (a b : ℕ) :
    derivative ((X : ℝ[X]) ^ (a + 1) * (X - 1) ^ (b + 1))
      = X ^ a * (X - 1) ^ b
          * (C ((a + b + 2 : ℕ) : ℝ) * X - C ((a + 1 : ℕ) : ℝ)) := by
  rw [derivative_mul, derivative_X_pow, derivative_pow, derivative_sub,
    derivative_X, derivative_one]
  push_cast
  simp only [C_add, C_1, map_ofNat]
  ring

/-! ## ★★性質 (e) -/

section Ratio

variable {α β : ℝ}

/-- ★`1 < α ≤ β` なら `(β−1)/(α−1) ≥ β/α`。

★`α(β−1) ≥ β(α−1)` ⟺ `β ≥ α` である。 -/
theorem sub_one_div_ge (hα : 1 < α) (hαβ : α ≤ β) :
    β / α ≤ (β - 1) / (α - 1) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hα1 : (0 : ℝ) < α - 1 := sub_pos.2 hα
  rw [div_le_div_iff₀ hα0 hα1]
  nlinarith

/-- ★`1 < α ≤ β` なら `β/α ≥ 1`。 -/
theorem one_le_div (hα : 1 < α) (hαβ : α ≤ β) : (1 : ℝ) ≤ β / α := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  rw [le_div_iff₀ hα0]
  linarith

/-- ★★**性質 (e)** —— `f(β)/f(α) ≥ (β/α)²`(`1 < α ≤ β`)。

原文 (NCBelyi p.3):
> (e) If α ∈ S\{∞} satisﬁes α > 1, then (β − 1)/(α − 1) ≥ β/α ≥ 1;

★`f(β)/f(α) = (β/α)^m·{(β−1)/(α−1)}^n ≥ (β/α)^{m+n} ≥ (β/α)²`。 -/
theorem belyi_ratio_ge (a b : ℕ) (hα : 1 < α) (hαβ : α ≤ β) :
    (β / α) ^ 2
      ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1)) / (α ^ (a + 1) * (α - 1) ^ (b + 1)) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hα1 : (0 : ℝ) < α - 1 := sub_pos.2 hα
  have hβ0 : (0 : ℝ) < β := lt_of_lt_of_le hα0 hαβ
  have hβ1 : (0 : ℝ) < β - 1 := lt_of_lt_of_le hα1 (by linarith)
  have hq1 : (1 : ℝ) ≤ β / α := one_le_div hα hαβ
  have hq0 : (0 : ℝ) < β / α := lt_of_lt_of_le zero_lt_one hq1
  -- 商を分解する
  have hsplit : (β ^ (a + 1) * (β - 1) ^ (b + 1)) / (α ^ (a + 1) * (α - 1) ^ (b + 1))
      = (β / α) ^ (a + 1) * ((β - 1) / (α - 1)) ^ (b + 1) := by
    rw [div_pow, div_pow, div_mul_div_comm]
  rw [hsplit]
  calc (β / α) ^ 2 ≤ (β / α) ^ (a + 1 + (b + 1)) := by
        refine pow_le_pow_right₀ hq1 ?_
        omega
    _ = (β / α) ^ (a + 1) * (β / α) ^ (b + 1) := pow_add _ _ _
    _ ≤ (β / α) ^ (a + 1) * ((β - 1) / (α - 1)) ^ (b + 1) := by
        refine mul_le_mul_of_nonneg_left ?_ (by positivity)
        exact pow_le_pow_left₀ hq0.le (sub_one_div_ge hα hαβ) _

/-- ★**`x > 1` での単調性** —— 性質 (e) の系。

★★**原文は導関数から出すが、要らなかった**——`(β/α)² ≥ 1` で足りる。 -/
theorem belyi_mono (a b : ℕ) (hα : 1 < α) (hαβ : α ≤ β) :
    α ^ (a + 1) * (α - 1) ^ (b + 1) ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hα1 : (0 : ℝ) < α - 1 := sub_pos.2 hα
  have hfα : (0 : ℝ) < α ^ (a + 1) * (α - 1) ^ (b + 1) := by positivity
  have hq1 : (1 : ℝ) ≤ β / α := one_le_div hα hαβ
  have h := belyi_ratio_ge a b hα hαβ
  have h1 : (1 : ℝ) ≤ (β / α) ^ 2 := one_le_pow₀ hq1
  rw [le_div_iff₀ hfα] at h
  nlinarith

end Ratio

/-! ## ★性質 (c) の根拠 -/

/-- ★`x ∈ [0,1]` では `|f(x)| ≤ 1`。 -/
theorem belyi_abs_le_one (a b : ℕ) {x : ℝ} (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    |x ^ (a + 1) * (x - 1) ^ (b + 1)| ≤ 1 := by
  rw [abs_mul, abs_pow, abs_pow]
  have hx : |x| ≤ 1 := abs_le.2 ⟨by linarith, h1⟩
  have hx1 : |x - 1| ≤ 1 := abs_le.2 ⟨by linarith, by linarith⟩
  calc |x| ^ (a + 1) * |x - 1| ^ (b + 1) ≤ 1 ^ (a + 1) * 1 ^ (b + 1) := by
        gcongr
    _ = 1 := by norm_num

/-- ★`β ≥ 2` なら `f(β) ≥ β`(とくに `f(β) > 1`)。 -/
theorem belyi_ge_self (a b : ℕ) {β : ℝ} (hβ : 2 ≤ β) :
    β ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := by
  have hβ0 : (0 : ℝ) < β := lt_of_lt_of_le zero_lt_two hβ
  have h1 : β ≤ β ^ (a + 1) := le_self_pow₀ (by linarith) (by omega)
  have h2 : (1 : ℝ) ≤ (β - 1) ^ (b + 1) := one_le_pow₀ (by linarith)
  nlinarith

/-! ## ★出典の紐付け(`.src`) -/

def belyi_ratio_ge.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 2,
    item := "Lemma 2.1(計算核——性質 (e) と導関数の恒等式のみ)",
    sectionId := "ncbelyi-lemma-2-1" }

def derivative_belyi.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 2,
    item := "Lemma 2.1(計算核——性質 (e) と導関数の恒等式のみ)",
    sectionId := "ncbelyi-lemma-2-1" }

end ABC3.Found.NCBelyi
