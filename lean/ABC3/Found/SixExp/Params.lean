/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Setup

/-!
# パラメータの選択 —— `L² > 2N³` と初段の不等式を同時に満たす

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を満たさねばならないか

* **Siegel** —— 未定係数の個数 `|s| = L²` が方程式の個数 `|box| = N³` の 2 倍以上:
  `2N³ < L²`。(2 倍にすると Siegel の指数 `p/(q−p)` が `≤ 1` になり、
  `SiegelConst.lean` の多項式形が使える。)
* **数え上げの初段** —— `Setup.lean` の `gap_of_base` に渡す不等式。
  house と `δ` の寄与をまとめると `L·(N+1)·γ + δ ≤ N³·log2 / 2` の形になる。

★この 2 つは**逆を向く** —— Siegel は `L` を大きく、数え上げは `L` を小さくしたい。
両立するのは `L ≈ N^{3/2}` のときで、そこが six exponentials が成り立つ理由である
(`N³ log2` が `N^{5/2}` に勝つ)。

## ★取り方

`L := ⌊√(2N³)⌋ + 1` と取る。すると
* `2N³ < L²` は `Nat.lt_succ_sqrt'` そのもの
* `(L−1)² ≤ 2N³` は `Nat.sqrt_le'` そのもの

で、後者から `L ≤ T`(`T := (N³log2/2 − δ)/((N+1)γ)`)が出る。
★**平方根を実数で扱わずに済む**のが要点 —— `(T−1)² > 2N³` を見ればよい
(`le_of_sq_pred_le`)。あとは `N` を十分大きく取る。
-/

namespace ABC3.Found.SixExp

/-- ★★`L ≤ T` の判定 —— `(L−1)² ≤ m` かつ `(T−1)² > m` かつ `1 ≤ T` なら `L ≤ T`。

★これで `Nat.sqrt` の評価を**平方根を取らずに**使える。 -/
theorem le_of_sq_pred_le {L : ℕ} {T m : ℝ}
    (hL : (((L : ℝ) - 1)) ^ 2 ≤ m) (hT : ((T - 1)) ^ 2 > m) (hT1 : 1 ≤ T)
    (hL1 : 1 ≤ (L : ℝ)) : (L : ℝ) ≤ T := by
  by_contra hlt
  push_neg at hlt
  have h1 : T - 1 < (L : ℝ) - 1 := by linarith
  have h2 : (0:ℝ) ≤ T - 1 := by linarith
  have h3 : (T - 1) ^ 2 < ((L : ℝ) - 1) ^ 2 := by nlinarith
  linarith

/-- ★★★★**パラメータの選択** —— `L² > 2N³`(Siegel)と初段の不等式を
**同時に満たす** `N`・`L` が存在する。

★`L := ⌊√(2N³)⌋ + 1` と取り、`N` を十分大きくする。
本質は `N³ log2` が `N^{5/2}` に勝つことである。 -/
theorem exists_params (γ δ : ℝ) (hγ : 0 < γ) (_hδ : 0 ≤ δ) :
    ∃ N L : ℕ, 0 < N ∧ 2 * N ^ 3 < L ^ 2 ∧
      (L : ℝ) * ((N : ℝ) + 1) * γ + δ ≤ (N : ℝ) ^ 3 * Real.log 2 / 2 := by
  have hlog2 : (0:ℝ) < Real.log 2 := Real.log_pos (by norm_num)
  set a : ℝ := Real.log 2 / (8 * γ) with ha
  have ha0 : 0 < a := by positivity
  obtain ⟨N, hN⟩ := exists_nat_gt
    (max (max 1 (4 * δ / Real.log 2)) (max ((2 + 2 * a) / a ^ 2) (1 / a)))
  have hN1 : (1:ℝ) < (N:ℝ) := lt_of_le_of_lt (le_trans (le_max_left _ _) (le_max_left _ _)) hN
  have hNpos : 0 < N := by exact_mod_cast lt_trans zero_lt_one hN1
  have hNR : (0:ℝ) < (N:ℝ) := by exact_mod_cast hNpos
  have hNd : 4 * δ / Real.log 2 < (N:ℝ) :=
    lt_of_le_of_lt (le_trans (le_max_right _ _) (le_max_left _ _)) hN
  have hNa : (2 + 2 * a) / a ^ 2 < (N:ℝ) :=
    lt_of_le_of_lt (le_trans (le_max_left _ _) (le_max_right _ _)) hN
  have hNinv : 1 / a < (N:ℝ) :=
    lt_of_le_of_lt (le_trans (le_max_right _ _) (le_max_right _ _)) hN
  have hN3 : (N:ℝ) ≤ (N:ℝ) ^ 3 := by
    have h := pow_le_pow_right₀ (le_of_lt hN1) (by norm_num : 1 ≤ 3)
    simpa using h
  have hN2 : (N:ℝ) ^ 2 ≤ (N:ℝ) ^ 3 := pow_le_pow_right₀ (le_of_lt hN1) (by norm_num)
  refine ⟨N, Nat.sqrt (2 * N ^ 3) + 1, hNpos, ?_, ?_⟩
  · have h := Nat.lt_succ_sqrt' (2 * N ^ 3)
    simpa using h
  · set L := Nat.sqrt (2 * N ^ 3) + 1 with hL
    have hLR1 : (1:ℝ) ≤ (L:ℝ) := by
      have h1 : 1 ≤ L := by omega
      exact_mod_cast h1
    have hLsq : (((L:ℝ)) - 1) ^ 2 ≤ 2 * (N:ℝ) ^ 3 := by
      have h1 : ((L:ℝ)) - 1 = ((Nat.sqrt (2 * N ^ 3) : ℕ) : ℝ) := by
        rw [hL]; push_cast; ring
      rw [h1]
      have h2 : (Nat.sqrt (2 * N ^ 3)) ^ 2 ≤ 2 * N ^ 3 := Nat.sqrt_le' _
      have h3 : (((Nat.sqrt (2 * N ^ 3) : ℕ) : ℝ)) ^ 2 ≤ (((2 * N ^ 3 : ℕ)) : ℝ) := by
        exact_mod_cast h2
      push_cast at h3
      exact h3
    set T : ℝ := ((N:ℝ) ^ 3 * Real.log 2 / 2 - δ) / (((N:ℝ) + 1) * γ) with hT
    have hden : (0:ℝ) < ((N:ℝ) + 1) * γ := by positivity
    have hδN : δ ≤ (N:ℝ) ^ 3 * Real.log 2 / 4 := by
      rw [div_lt_iff₀ hlog2] at hNd
      nlinarith
    have haN2 : (1:ℝ) ≤ a * (N:ℝ) ^ 2 := by
      rw [div_lt_iff₀ ha0] at hNinv
      nlinarith
    have hTge : a * (N:ℝ) ^ 2 ≤ T := by
      rw [hT, le_div_iff₀ hden]
      have hexp : a * (N:ℝ) ^ 2 * (((N:ℝ) + 1) * γ)
          = (N:ℝ) ^ 2 * ((N:ℝ) + 1) * Real.log 2 / 8 := by
        rw [ha]
        field_simp
      rw [hexp]
      nlinarith
    have hT1 : (1:ℝ) ≤ T := le_trans haN2 hTge
    have hTsq : (T - 1) ^ 2 > 2 * (N:ℝ) ^ 3 := by
      have hbase : (a * (N:ℝ) ^ 2 - 1) ^ 2 > 2 * (N:ℝ) ^ 3 := by
        rw [div_lt_iff₀ (by positivity : (0:ℝ) < a ^ 2)] at hNa
        nlinarith
      have h1 : (0:ℝ) ≤ a * (N:ℝ) ^ 2 - 1 := by linarith
      nlinarith
    have hLT : (L:ℝ) ≤ T := le_of_sq_pred_le hLsq hTsq hT1 hLR1
    calc (L:ℝ) * ((N:ℝ) + 1) * γ + δ = (L:ℝ) * (((N:ℝ) + 1) * γ) + δ := by ring
      _ ≤ T * (((N:ℝ) + 1) * γ) + δ := by
          have h := mul_le_mul_of_nonneg_right hLT (le_of_lt hden)
          linarith
      _ = ((N:ℝ) ^ 3 * Real.log 2 / 2 - δ) + δ := by
          rw [hT, div_mul_cancel₀ _ (ne_of_gt hden)]
      _ = (N:ℝ) ^ 3 * Real.log 2 / 2 := by ring

end ABC3.Found.SixExp
