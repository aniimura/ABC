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

/-! ## ★★`f₀ ≤ 1/4` —— 原文の 4 場合分けは要らなかった -/

/-- ★`r = m/(m+n)` での値の絶対値は `(m/(m+n))^m·(n/(m+n))^n`。

原文 (NCBelyi p.3):
> f0 = |f(r)| = {m/(m + n)}m · {n/(m + n)}n
-/
theorem abs_belyi_at_r (a b : ℕ) :
    |((a + 1 : ℝ) / (a + 1 + (b + 1))) ^ (a + 1)
        * ((a + 1 : ℝ) / (a + 1 + (b + 1)) - 1) ^ (b + 1)|
      = ((a + 1 : ℝ) / (a + 1 + (b + 1))) ^ (a + 1)
          * ((b + 1 : ℝ) / (a + 1 + (b + 1))) ^ (b + 1) := by
  have hd : (0 : ℝ) < (a : ℝ) + 1 + ((b : ℝ) + 1) := by positivity
  have hsub : ((a + 1 : ℝ) / (a + 1 + (b + 1)) - 1)
      = -((b + 1 : ℝ) / (a + 1 + (b + 1))) := by
    field_simp
    ring
  rw [hsub, abs_mul, abs_pow, abs_pow, abs_neg, abs_of_nonneg (by positivity),
    abs_of_nonneg (by positivity)]

/-- ★★**`f₀ ≤ 1/4`**。

原文 (NCBelyi p.3):
> Note,   moreover,          that     this   expression       for  f0   implies   that   0 < f0      ≤   1  .  [Indeed,  this

★★**原文は 4 つの場合分けで確かめている**
(`m, n ≥ 2` / `m = n = 1` / 片方が 1 で他方が `≥ 3` / 片方が 1 で他方が 2)。
★**要らなかった。** `p ≝ m/(m+n)`、`q ≝ n/(m+n)` と置けば `p + q = 1` で

> `p^m·q^n = (p·q)·(p^{m−1}·q^{n−1}) ≤ (1/4)·1`

——**AM–GM(`p·q ≤ ((p+q)/2)² = 1/4`)1 回**で済む。
★「正面から要ると思ったものが要らなかった」の型がまた出た。 -/
theorem belyi_f0_le_quarter (a b : ℕ) :
    ((a + 1 : ℝ) / (a + 1 + (b + 1))) ^ (a + 1)
        * ((b + 1 : ℝ) / (a + 1 + (b + 1))) ^ (b + 1)
      ≤ 1 / 4 := by
  set d : ℝ := (a : ℝ) + 1 + ((b : ℝ) + 1) with hd
  have hd0 : (0 : ℝ) < d := by rw [hd]; positivity
  set p : ℝ := ((a : ℝ) + 1) / d with hp
  set q : ℝ := ((b : ℝ) + 1) / d with hq
  have hp0 : (0 : ℝ) < p := by rw [hp]; positivity
  have hq0 : (0 : ℝ) < q := by rw [hq]; positivity
  have hsum : p + q = 1 := by
    rw [hp, hq, div_add_div _ _ hd0.ne' hd0.ne']
    field_simp
    ring
  have hp1 : p ≤ 1 := by linarith
  have hq1 : q ≤ 1 := by linarith
  have hpq : p * q ≤ 1 / 4 := by nlinarith [sq_nonneg (p - q)]
  have hpa : p ^ a ≤ 1 := pow_le_one₀ hp0.le hp1
  have hqb : q ^ b ≤ 1 := pow_le_one₀ hq0.le hq1
  calc p ^ (a + 1) * q ^ (b + 1) = (p * q) * (p ^ a * q ^ b) := by ring
    _ ≤ (1 / 4) * (1 * 1) := by
        refine mul_le_mul hpq ?_ (by positivity) (by norm_num)
        exact mul_le_mul hpa hqb (by positivity) zero_le_one
    _ = 1 / 4 := by norm_num

/-! ## ★性質 (c) —— `f(β) ∉ f(S)` の根拠 -/

/-- ★**`β ≥ 2` かつ `x ∈ [0,1]` なら `f(β) ≠ f(x)`**。

★`f(β) ≥ β ≥ 2 > 1` かつ `|f(x)| ≤ 1` である。 -/
theorem belyi_ne_of_mem_unitInterval (a b : ℕ) {β x : ℝ} (hβ : 2 ≤ β)
    (h0 : 0 ≤ x) (h1 : x ≤ 1) :
    β ^ (a + 1) * (β - 1) ^ (b + 1) ≠ x ^ (a + 1) * (x - 1) ^ (b + 1) := by
  have hfβ : β ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := belyi_ge_self a b hβ
  have hfx : |x ^ (a + 1) * (x - 1) ^ (b + 1)| ≤ 1 := belyi_abs_le_one a b h0 h1
  have hfx' : x ^ (a + 1) * (x - 1) ^ (b + 1) ≤ 1 := le_trans (le_abs_self _) hfx
  intro hcon
  rw [hcon] at hfβ
  linarith

/-- ★**`1 < α < β` なら `f(α) < f(β)`**(狭義単調)。 -/
theorem belyi_strictMono (a b : ℕ) {α β : ℝ} (hα : 1 < α) (hαβ : α < β) :
    α ^ (a + 1) * (α - 1) ^ (b + 1) < β ^ (a + 1) * (β - 1) ^ (b + 1) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hα1 : (0 : ℝ) < α - 1 := sub_pos.2 hα
  have hfα : (0 : ℝ) < α ^ (a + 1) * (α - 1) ^ (b + 1) := by positivity
  have hq : (1 : ℝ) < β / α := (one_lt_div hα0).2 hαβ
  have h := belyi_ratio_ge a b hα hαβ.le
  have h1 : (1 : ℝ) < (β / α) ^ 2 := one_lt_pow₀ hq (by norm_num)
  rw [le_div_iff₀ hfα] at h
  nlinarith

/-! ## ★★★性質 (d) —— 原文の 4 つの場合

原文 (NCBelyi p.3):
> [by property (e)] or f(α) ≤ f0, in which case

★**集合 `S` をモデル化せずに、純粋な実数の不等式として取る。**
原文の場合分け(`n` 偶/奇 × `α > 1` / `α ∈ {0,1}` / `α = r`)は、
下の 4 本のいずれかに帰着する。
★★こうすると `Lemma 2.2` の帰納法を書くときに**そのまま使える形**になる。
-/

section PropD

variable {C α β f0 : ℝ}

/-- ★**(d) の場合 1** —— `f₀ = 0`、`α > 1`(原文の `n` 偶かつ `α > 1`)。

`f(β)/f(α) ≥ (β/α)² ≥ β/α ≥ C`。 -/
theorem belyi_d_pos (a b : ℕ) (hC : 2 ≤ C) (hα : 1 < α) (hCα : C * α ≤ β) :
    C ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1)) / (α ^ (a + 1) * (α - 1) ^ (b + 1)) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hq : C ≤ β / α := (le_div_iff₀ hα0).2 hCα
  have hq1 : (1 : ℝ) ≤ β / α := le_trans (by linarith) hq
  have hαβ : α ≤ β := by
    have := (le_div_iff₀ hα0).1 hq1
    linarith
  have h := belyi_ratio_ge a b hα hαβ
  have h2 : β / α ≤ (β / α) ^ 2 := by nlinarith
  linarith

/-- ★**(d) の場合 2** —— `f₀ = 0`、`α = r`(`0 < f(r) ≤ 1`)。

`f(β)/f(r) ≥ f(β) ≥ β ≥ C`。 -/
theorem belyi_d_r (a b : ℕ) (hC : 2 ≤ C) (hβ : C ≤ β)
    (hfr0 : 0 < f0) (hfr1 : f0 ≤ 1) :
    C ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1)) / f0 := by
  have hβ2 : (2 : ℝ) ≤ β := le_trans hC hβ
  have hfβ : β ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := belyi_ge_self a b hβ2
  rw [le_div_iff₀ hfr0]
  nlinarith

/-- ★★**(d) の場合 3** —— `f₀ > 0` かつ `f(α) ≥ f₀`(原文の `n` 奇の第 1 の場合)。

`(f(β)+f₀)/(f(α)+f₀) ≥ f(β)/(2f(α)) ≥ (1/2)(β/α)² ≥ β/α ≥ C`。
★最後から 2 つ目の不等号に **`β/α ≥ 2`** が要る——原文が `C ≥ 2` を仮定する理由である。 -/
theorem belyi_d_ge (a b : ℕ) (hC : 2 ≤ C) (hα : 1 < α) (hCα : C * α ≤ β)
    (hf0 : 0 < f0) (hge : f0 ≤ α ^ (a + 1) * (α - 1) ^ (b + 1)) :
    C ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1) + f0)
          / (α ^ (a + 1) * (α - 1) ^ (b + 1) + f0) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hα1 : (0 : ℝ) < α - 1 := sub_pos.2 hα
  have hfα : (0 : ℝ) < α ^ (a + 1) * (α - 1) ^ (b + 1) := by positivity
  have hq : C ≤ β / α := (le_div_iff₀ hα0).2 hCα
  have hq2 : (2 : ℝ) ≤ β / α := le_trans hC hq
  have hαβ : α ≤ β := by
    have : (1 : ℝ) ≤ β / α := le_trans (by linarith) hq2
    have := (le_div_iff₀ hα0).1 this
    linarith
  have h := belyi_ratio_ge a b hα hαβ
  rw [le_div_iff₀ hfα] at h
  have hden : (0 : ℝ) < α ^ (a + 1) * (α - 1) ^ (b + 1) + f0 := by linarith
  rw [le_div_iff₀ hden]
  -- `f(β) ≥ (β/α)²·f(α)` と `(β/α)² ≥ 2·(β/α)` から
  have hsq : 2 * (β / α) ≤ (β / α) ^ 2 := by nlinarith
  nlinarith

/-- ★★**(d) の場合 4** —— `f₀ > 0` かつ `f(α) ≤ f₀ ≤ 1/4`(原文の `n` 奇の第 2 の場合)。

`(f(β)+f₀)/(f(α)+f₀) ≥ f(β)/(2f₀) ≥ 2f(β) ≥ 2β ≥ C`。
★**`f₀ ≤ 1/4`(`belyi_f0_le_quarter`)がここで効く。** -/
theorem belyi_d_le (a b : ℕ) (hC : 2 ≤ C) (hα : 1 < α) (hCα : C * α ≤ β)
    (hf0 : 0 < f0) (hf04 : f0 ≤ 1 / 4)
    (hle : α ^ (a + 1) * (α - 1) ^ (b + 1) ≤ f0)
    (hfα0 : 0 ≤ α ^ (a + 1) * (α - 1) ^ (b + 1)) :
    C ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1) + f0)
          / (α ^ (a + 1) * (α - 1) ^ (b + 1) + f0) := by
  have hα0 : (0 : ℝ) < α := lt_trans zero_lt_one hα
  have hβC : C ≤ β := by nlinarith
  have hβ2 : (2 : ℝ) ≤ β := le_trans hC hβC
  have hfβ : β ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := belyi_ge_self a b hβ2
  have hden : (0 : ℝ) < α ^ (a + 1) * (α - 1) ^ (b + 1) + f0 := by linarith
  rw [le_div_iff₀ hden]
  -- 分母は `2·f₀ ≤ 1/2` 以下、分子は `f(β) ≥ β ≥ C` 以上
  nlinarith

/-- ★★**(d) の場合 4(一般形)** —— `α > 1` を仮定しない版。

`f(α) ≤ f₀ ≤ 1/4` と `C ≤ β` だけから
`(f(β)+f₀)/(f(α)+f₀) ≥ C` が出る。

★★**`α ∈ {0,1}` の場合(原文の最後の場合)もこれで覆える**——
そこでは `f(α) = 0 ≤ f₀` だからである。
★分母は `f(α) + f₀ ≤ 2f₀ ≤ 1/2`、分子は `f(β) + f₀ ≥ β ≥ C` である。 -/
theorem belyi_d_le' (a b : ℕ) {fα : ℝ} (hC : 2 ≤ C) (hβ : C ≤ β)
    (hf0 : 0 < f0) (hf04 : f0 ≤ 1 / 4) (hfα : fα ≤ f0) (hpos : 0 < fα + f0) :
    C ≤ (β ^ (a + 1) * (β - 1) ^ (b + 1) + f0) / (fα + f0) := by
  have hβ2 : (2 : ℝ) ≤ β := le_trans hC hβ
  have hfβ : β ≤ β ^ (a + 1) * (β - 1) ^ (b + 1) := belyi_ge_self a b hβ2
  rw [le_div_iff₀ hpos]
  nlinarith

end PropD

/-! ## ★★性質 (c) を集合の水準で -/

/-- ★★**性質 (c)** —— `f(β) ∉ f(S)`。

原文 (NCBelyi p.3):
> that property (c) is satisﬁed.

★`S` の元は原文の (i)(ii)(iii) より `0`、`r`、`1`、または `> 1` である。
- `0, r, 1 ∈ [0,1]` では `|f| ≤ 1 < f(β)`
- `α > 1` では `α < β` より `f(α) < f(β)`(狭義単調)

★★**原文は `S ⊆ ℙ¹(ℚ)` だが、有理性は (c) の証明に使われない**——
`S : Finset ℝ` で述べておく(仮説が弱いので原文の場合を含む)。 -/
theorem belyi_c (a b : ℕ) {C r β : ℝ} (hC : 2 ≤ C) (S : Finset ℝ)
    (hr0 : 0 ≤ r) (hr1 : r ≤ 1)
    (hS : ∀ α ∈ S, α = 0 ∨ α = r ∨ α = 1 ∨ 1 < α)
    (hone : (1 : ℝ) ∈ S)
    (hβ : ∀ α ∈ S, α ≠ 0 → C * α ≤ β) :
    ∀ α ∈ S, β ^ (a + 1) * (β - 1) ^ (b + 1) ≠ α ^ (a + 1) * (α - 1) ^ (b + 1) := by
  have hβ2 : (2 : ℝ) ≤ β := by
    have h := hβ 1 hone one_ne_zero
    linarith
  intro α hα
  rcases hS α hα with h | h | h | hgt
  · rw [h]; exact belyi_ne_of_mem_unitInterval a b hβ2 le_rfl zero_le_one
  · rw [h]; exact belyi_ne_of_mem_unitInterval a b hβ2 hr0 hr1
  · rw [h]; exact belyi_ne_of_mem_unitInterval a b hβ2 zero_le_one le_rfl
  · have hαβ : α < β := by
      have h := hβ α hα (by linarith)
      nlinarith
    exact ne_of_gt (belyi_strictMono a b hgt hαβ)

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
