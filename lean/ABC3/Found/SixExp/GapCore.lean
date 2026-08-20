/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.SiegelPoly

/-!
# 数え上げの本体 —— 対数で書いた初段の不等式から全段へ

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★形

外挿の `hgap` は、定数を具体に置くと

  `K₀ · Z₀^{Lb·(n+1)} < 2^{n³}`   (すべての `n ≥ N`)

の形になる(`Lb = 2L'` は `p+q` の上界、`Z₀` は house・解析側・分母の増大率をまとめたもの、
`K₀` は `n` に依らない定数)。

★`gap_core` はこれを**対数で書いた初段の不等式 1 本**

  `log K₀ + Lb·(N+1)·log Z₀ ≤ N³·log2 / 2`

から出す。`Params.lean` の `exists_params` がちょうどこの形を与える。

## ★2 段の作り

1. 初段 —— `lt_two_pow_of_log` で対数から積の形に戻す。
2. 全段 —— `Setup.lean` の `gap_of_base`。
   `Z := Z₀^{Lb} ≤ 2^{3N²}` は、上の不等式を `N+1` で割れば出る
   (`N³/(2(N+1)) ≤ 3N²` はつねに真)。
-/

namespace ABC3.Found.SixExp

/-- ★対数で書いた不等式を積の形に戻す。 -/
theorem lt_two_pow_of_log {A Z : ℝ} (hA : 0 < A) (hZ : 0 < Z) {m T : ℕ}
    (h : Real.log A + m * Real.log Z < T * Real.log 2) :
    A * Z ^ m < 2 ^ T := by
  have hAZ : (0:ℝ) < A * Z ^ m := by positivity
  have h2 : (0:ℝ) < (2:ℝ) ^ T := by positivity
  rw [← Real.log_lt_log_iff hAZ h2]
  rw [Real.log_mul (ne_of_gt hA) (by positivity), Real.log_pow, Real.log_pow]
  exact h

/-- ★対数で書いた不等式(`≤`)を冪の形に戻す。 -/
theorem le_two_pow_of_log {Z : ℝ} (hZ : 0 < Z) {T : ℕ}
    (h : Real.log Z ≤ T * Real.log 2) : Z ≤ 2 ^ T := by
  have h2 : (0:ℝ) < (2:ℝ) ^ T := by positivity
  rw [← Real.log_le_log_iff hZ h2, Real.log_pow]
  exact h

/-- ★★★★**数え上げの本体** —— 対数で書いた初段の不等式から、
**すべての段の不等式**が出る。 -/
theorem gap_core {K₀ Z₀ : ℝ} (hK₀ : 1 ≤ K₀) (hZ₀ : 1 ≤ Z₀) {Lb N : ℕ} (hN : 1 ≤ N)
    (hlog : Real.log K₀ + (Lb:ℝ) * ((N:ℝ) + 1) * Real.log Z₀
      ≤ (N:ℝ) ^ 3 * Real.log 2 / 2) :
    ∀ n, N ≤ n → K₀ * Z₀ ^ (Lb * (n + 1)) < 2 ^ (n ^ 3) := by
  have hlog2 : (0:ℝ) < Real.log 2 := Real.log_pos (by norm_num)
  have hK0 : (0:ℝ) < K₀ := lt_of_lt_of_le zero_lt_one hK₀
  have hZ0 : (0:ℝ) < Z₀ := lt_of_lt_of_le zero_lt_one hZ₀
  have hlZ : (0:ℝ) ≤ Real.log Z₀ := Real.log_nonneg hZ₀
  have hlK : (0:ℝ) ≤ Real.log K₀ := Real.log_nonneg hK₀
  have hNR : (1:ℝ) ≤ (N:ℝ) := by exact_mod_cast hN
  have hN3 : (0:ℝ) < (N:ℝ) ^ 3 := by positivity
  set Z : ℝ := Z₀ ^ Lb with hZdef
  have hZ1 : (1:ℝ) ≤ Z := one_le_pow₀ hZ₀
  have hKZ : (0:ℝ) < K₀ * Z := by positivity
  have hrwgen : ∀ n : ℕ, (K₀ * Z) * Z ^ n = K₀ * Z₀ ^ (Lb * (n + 1)) := by
    intro n
    rw [hZdef, ← pow_mul, mul_assoc, ← pow_add]
    congr 2
    ring
  have hbase : (K₀ * Z) * Z ^ N < 2 ^ (N ^ 3) := by
    rw [hrwgen N]
    refine lt_two_pow_of_log hK0 hZ0 ?_
    have hcast : ((Lb * (N + 1) : ℕ) : ℝ) = (Lb:ℝ) * ((N:ℝ) + 1) := by push_cast; ring
    have hcast2 : (((N ^ 3 : ℕ)) : ℝ) = (N:ℝ) ^ 3 := by push_cast; ring
    rw [hcast, hcast2]
    calc Real.log K₀ + (Lb:ℝ) * ((N:ℝ) + 1) * Real.log Z₀
        ≤ (N:ℝ) ^ 3 * Real.log 2 / 2 := hlog
      _ < (N:ℝ) ^ 3 * Real.log 2 := by nlinarith
  have hZle : Z ≤ 2 ^ (3 * N ^ 2) := by
    refine le_two_pow_of_log (by positivity) ?_
    rw [hZdef, Real.log_pow]
    have h1 : (Lb:ℝ) * ((N:ℝ) + 1) * Real.log Z₀ ≤ (N:ℝ) ^ 3 * Real.log 2 / 2 := by linarith
    have h2 : (Lb:ℝ) * Real.log Z₀ * ((N:ℝ) + 1) ≤ (N:ℝ) ^ 3 * Real.log 2 / 2 := by
      calc (Lb:ℝ) * Real.log Z₀ * ((N:ℝ) + 1) = (Lb:ℝ) * ((N:ℝ) + 1) * Real.log Z₀ := by ring
        _ ≤ _ := h1
    have hcast : ((3 * N ^ 2 : ℕ) : ℝ) = 3 * (N:ℝ) ^ 2 := by push_cast; ring
    rw [hcast]
    nlinarith
  intro n hn
  have h := gap_of_base hKZ hZ1 hbase hZle n hn
  rwa [hrwgen n] at h

end ABC3.Found.SixExp
