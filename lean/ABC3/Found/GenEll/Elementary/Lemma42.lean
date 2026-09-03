/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import ABC3.Found.GenEll.Elementary.Lemma36

/-!
# Elementary —— `[GenEll] Lemma 4.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open Real

/-! ## Lemma 4.2 —— Some Elementary Estimates -/

/-- 原文が使う初等的事実の**強い版**。

原文 (GenEll p.21):
> follows from the easily verified fact that log(H +1) ≤ (3H/2)·log(2) for any positive integer H.

原文は「任意の正整数 `H` について `log(H+1) ≤ (3H/2)·log(2)`」を使う。
★**我々は `H + 1 ≤ 2^H`(`Nat.lt_two_pow_self`)から `log(H+1) ≤ H·log 2` を出す。**
係数 `3/2` が要らない分だけ強い。 -/
theorem log_succ_le_mul_log_two (H : ℕ) (hH : 0 < H) :
    Real.log ((H : ℝ) + 1) ≤ (H : ℝ) * Real.log 2 := by
  have hlt : H < 2 ^ H := Nat.lt_two_pow_self
  have hle : (H : ℝ) + 1 ≤ (2 : ℝ) ^ H := by
    have : (H : ℝ) + 1 ≤ ((2 ^ H : ℕ) : ℝ) := by exact_mod_cast hlt
    simpa using this
  have h0 : (0 : ℝ) < (H : ℝ) + 1 := by positivity
  calc Real.log ((H : ℝ) + 1)
      ≤ Real.log ((2 : ℝ) ^ H) := Real.log_le_log h0 hle
    _ = (H : ℝ) * Real.log 2 := by rw [Real.log_pow]

/-- **[GenEll] Lemma 4.2**(Some Elementary Estimates)。

原文 (GenEll p.21):
> Lemma 4.2. (Some Elementary Estimates) Let n be a positive integer;

`n` 正整数、`p_1,…,p_n` 素数、`h_1,…,h_n` 正整数、`h ≝ Σ h_j·log(p_j)` とすると
**`Σ log(p_j) ≤ h`** かつ **`Σ log(h_j) ≤ Σ log(h_j+1) ≤ 3h/2`**。

★3 つの不等式の道筋:
1. `h_j ≥ 1` と `log(p_j) ≥ 0` から各項ごとに `log(p_j) ≤ h_j·log(p_j)`。
2. `log` の単調性(`h_j ≤ h_j + 1`)。
3. ★`log(h_j+1) ≤ h_j·log 2 ≤ h_j·log(p_j)`(`p_j ≥ 2`)。総和して `≤ h ≤ 3h/2`。

★**3 で得ているのは原文より強い `≤ h` である**が、
主張は原文どおり `≤ 3h/2` のまま書いてある(G5——主張を強めも弱めもしない)。 -/
theorem lemma_4_2 (n : ℕ) (_hn : 0 < n) (p : Fin n → ℕ) (hp : ∀ j, (p j).Prime)
    (hh : Fin n → ℕ) (hhpos : ∀ j, 0 < hh j) :
    (∑ j, Real.log (p j)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j))
  ∧ (∑ j, Real.log (hh j)) ≤ (∑ j, Real.log ((hh j : ℝ) + 1))
  ∧ (∑ j, Real.log ((hh j : ℝ) + 1))
      ≤ 3 / 2 * (∑ j, (hh j : ℝ) * Real.log (p j)) := by
  -- 各素数について `log 2 ≤ log(p j)`(とくに `0 ≤ log(p j)`)
  have hlog2 : ∀ j, Real.log 2 ≤ Real.log (p j) := fun j =>
    Real.log_le_log (by norm_num) (by exact_mod_cast (hp j).two_le)
  have hlogpos : ∀ j, 0 ≤ Real.log (p j) := fun j =>
    le_trans (Real.log_nonneg (by norm_num)) (hlog2 j)
  have hh1 : ∀ j, (1 : ℝ) ≤ (hh j : ℝ) := fun j => by exact_mod_cast hhpos j
  -- ★第 1 不等式
  have first : (∑ j, Real.log (p j)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_le_sum fun j _ => le_mul_of_one_le_left (hlogpos j) (hh1 j)
  -- ★第 3 不等式の主部: `Σ log(h_j+1) ≤ Σ h_j·log(p_j)`
  have third' : (∑ j, Real.log ((hh j : ℝ) + 1)) ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_le_sum fun j _ =>
      le_trans (log_succ_le_mul_log_two (hh j) (hhpos j))
        (mul_le_mul_of_nonneg_left (hlog2 j) (by positivity))
  have hsum0 : 0 ≤ (∑ j, (hh j : ℝ) * Real.log (p j)) :=
    Finset.sum_nonneg fun j _ => mul_nonneg (by positivity) (hlogpos j)
  refine ⟨first, ?_, ?_⟩
  · -- ★第 2 不等式
    exact Finset.sum_le_sum fun j _ =>
      Real.log_le_log (by exact_mod_cast hhpos j) (by linarith)
  · -- ★第 3 不等式(原文どおり `3h/2` で書く)
    linarith

def lemma_4_2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Lemma 4.2",
    sectionId := "genell-lemma-4-2" }

def log_succ_le_mul_log_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 21, item := "Lemma 4.2",
    sectionId := "genell-lemma-4-2" }

end ABC3.Found.GenEll
