/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.House

/-!
# Liouville 不等式 —— チェーン `sixexp` の葉 `liouville-ineq`

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーンの葉。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を作るか

★★**0 でない代数的整数は小さくなれない**:

  `1 ≤ ‖σ α‖ · house(α)^{d-1}`   (`d = [K:ℚ]`)

超越性証明では**この対偶**を使う —— 解析側で「値が小さい」ことを示し、
算術側の下界と突き合わせて **値が 0 である**と結論する
(`eq_zero_of_norm_embedding_lt`)。

## ★★測定(2026-08-18)—— mathlib はここまで持っていた

| 要るもの | mathlib | 判定 |
|---|---|---|
| `house`(代数的数の共役の最大値) | `NumberField.house` | ★**ある** |
| `‖σ α‖ ≤ house α` | `NumberField.norm_embedding_le_house` | ★**ある** |
| 0 でない代数的整数は house ≥ 1 | `NumberField.one_le_house_of_isIntegral` | ★**ある** |
| ★**`‖N(α)‖ ≤ ‖σ α‖ · house(α)^{d-1}`** | `NumberField.norm_norm_le_norm_mul_house_pow` | ★★**ある** |
| ノルムが 0 でない整数 | `Algebra.norm_ne_zero_iff` + `Algebra.coe_norm_int` | ★**ある** |
| **`1 ≤ ‖N(α)‖` から下界を組む** | —— | ★**無い**(ここが本ファイル) |

★★**すなわち残っていたのは最後の 1 段だけだった。**
2026-08-17 に「mathlib に超越論の在庫が無い」と書いたのは、
存在しないディレクトリ(`Analysis/Transcendental/`)を探した結果である
(`Gap/FrdI/Section6.lean` の訂正欄)。
-/

namespace ABC3.Found.SixExp

open NumberField

variable {K : Type*} [Field K] [NumberField K]

/-- 0 でない代数的整数のノルムは絶対値 1 以上(`N(α) ∈ ℤ \ {0}` だから)。 -/
theorem one_le_norm_norm {α : 𝓞 K} (hα : α ≠ 0) : 1 ≤ ‖Algebra.norm ℚ (α : K)‖ := by
  have hz : Algebra.norm ℤ α ≠ 0 := Algebra.norm_ne_zero_iff.mpr hα
  rw [← Algebra.coe_norm_int α]
  have h1 : (1 : ℤ) ≤ |Algebra.norm ℤ α| := Int.one_le_abs (by omega)
  have h2 : ‖((Algebra.norm ℤ α : ℤ) : ℚ)‖ = |((Algebra.norm ℤ α : ℤ) : ℝ)| := by
    simp [Int.norm_eq_abs]
  rw [h2, ← Int.cast_abs]
  exact_mod_cast h1

/-- ★★**Liouville 不等式** —— `1 ≤ ‖σ α‖ · house(α)^{d-1}`。

★`∏_τ ‖τ α‖ = ‖N(α)‖ ≥ 1` を `σ` の項と残りに分け、残りを `house` で抑える。
分け方は mathlib の `norm_norm_le_norm_mul_house_pow` が済ませている。 -/
theorem one_le_norm_embedding_mul_house_pow {α : 𝓞 K} (hα : α ≠ 0) (σ : K →+* ℂ) :
    1 ≤ ‖σ (α : K)‖ * house (α : K) ^ (Module.finrank ℚ K - 1) :=
  (one_le_norm_norm hα).trans (norm_norm_le_norm_mul_house_pow (α : K) σ)

theorem one_le_house_pow {α : 𝓞 K} (hα : α ≠ 0) (n : ℕ) : 1 ≤ house (α : K) ^ n :=
  one_le_pow₀ (one_le_house_of_isIntegral α.2 (RingOfIntegers.coe_ne_zero_iff.mpr hα))

/-- ★**割った形** —— `house(α)^{-(d-1)} ≤ ‖σ α‖`。 -/
theorem inv_house_pow_le_norm_embedding {α : 𝓞 K} (hα : α ≠ 0) (σ : K →+* ℂ) :
    (house (α : K) ^ (Module.finrank ℚ K - 1))⁻¹ ≤ ‖σ (α : K)‖ := by
  have hH : (0:ℝ) < house (α : K) ^ (Module.finrank ℚ K - 1) :=
    lt_of_lt_of_le zero_lt_one (one_le_house_pow hα _)
  rw [inv_le_iff_one_le_mul₀ hH]
  exact one_le_norm_embedding_mul_house_pow hα σ

/-- ★★★**外挿で実際に使う形** —— `house` の**上界 `H`** が分かっているとき、
値が `H^{-(d-1)}` より小さければ **その代数的整数は 0 である**。

★これが「解析側の小ささ」と「算術側の下界」を突き合わせる段そのものである。
超越性証明の帰納は毎回この形で回る。 -/
theorem eq_zero_of_norm_embedding_lt {α : 𝓞 K} {H : ℝ} (σ : K →+* ℂ) (hH1 : 1 ≤ H)
    (hHle : house (α : K) ≤ H)
    (h : ‖σ (α : K)‖ < (H ^ (Module.finrank ℚ K - 1))⁻¹) : α = 0 := by
  by_contra hα
  have hHpos : (0:ℝ) < H ^ (Module.finrank ℚ K - 1) := by positivity
  have h1 : 1 ≤ ‖σ (α : K)‖ * house (α : K) ^ (Module.finrank ℚ K - 1) :=
    one_le_norm_embedding_mul_house_pow hα σ
  have h2 : ‖σ (α : K)‖ * house (α : K) ^ (Module.finrank ℚ K - 1)
      ≤ ‖σ (α : K)‖ * H ^ (Module.finrank ℚ K - 1) := by
    gcongr
    exact norm_nonneg _
  have h3 : ‖σ (α : K)‖ * H ^ (Module.finrank ℚ K - 1)
      < (H ^ (Module.finrank ℚ K - 1))⁻¹ * H ^ (Module.finrank ℚ K - 1) :=
    mul_lt_mul_of_pos_right h hHpos
  rw [inv_mul_cancel₀ hHpos.ne'] at h3
  exact absurd (lt_of_le_of_lt (h1.trans h2) h3) (lt_irrefl 1)

/-- ★**負の対照** —— `α = 0` では下界は成り立たない(不等式が非空虚であることの確認)。 -/
theorem not_one_le_norm_embedding_mul_house_pow_zero (σ : K →+* ℂ) :
    ¬ (1 ≤ ‖σ ((0 : 𝓞 K) : K)‖ * house ((0 : 𝓞 K) : K) ^ (Module.finrank ℚ K - 1)) := by
  simp

end ABC3.Found.SixExp
