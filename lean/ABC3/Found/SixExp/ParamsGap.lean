/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.GapCore

/-!
# 数え上げの完成形 —— パラメータを選ぶと全段の不等式が出る

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★組み上がる形

外挿の `hgap` を具体的な定数で展開すると

  `Δ k · M k · (H k)^e = K₀ · Z₀^{2L'·(n+1)}`   (`n = N + k`)

になる。ここで

* `Z₀ = Bβ³ · exp(30·X·Y) · A₀^{3e}` —— **段ごとに 1 乗ずつ増える分**
  (分母の倍率・解析側の指数・house の増大をまとめたもの)
* `K₀ = (L'²·Cbound)^{1+e}`、`Cbound = cS²·L'²·A₀^{6L'N}` —— **段に依らない分**
  (Siegel が与える係数の house の上界)

★`exists_params_gap` は「`N`・`L'` を選べば **すべての段で**
`K₀ · Z₀^{2L'(n+1)} < 2^{n³}`」を与える。

## ★2 段の作り

1. `exists_params_for_gap` —— `Params.lean` の `exists_params` に
   `γ := (1+e)(4 + 6 log A₀) + 2 log Z₀`、`δ := 2(1+e) log cS` を当てる。
   ★`log L' ≤ L' ≤ L'(N+1)` と `L'N ≤ L'(N+1)` で、`log K₀` の全項が
   `L'(N+1)` の定数倍に収まるのが要点。
2. `exists_params_gap` —— それを `GapCore.lean` の `gap_core` に渡す。
-/

namespace ABC3.Found.SixExp

/-- ★★★★**`gap_core` に渡す形のパラメータ選択**。

`K₀ = (L'²·Cbound)^{1+e}`、`Cbound = cS²·L'²·A₀^{6L'N}` のとき、
`log K₀ + 2L'(N+1)·log Z₀ ≤ N³ log2 / 2` を満たす `N`・`L'` がある。 -/
theorem exists_params_for_gap (A₀ Z₀ cS : ℝ) (hA₀ : 1 ≤ A₀) (hZ₀ : 1 ≤ Z₀) (hcS : 1 ≤ cS)
    (e : ℕ) :
    ∃ N L' : ℕ, 0 < N ∧ 2 * N ^ 3 < L' ^ 2 ∧
      ((1:ℝ) + e) * Real.log (((L':ℝ) ^ 2) * (cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N)))
        + (2 * (L':ℝ)) * ((N:ℝ) + 1) * Real.log Z₀
      ≤ (N:ℝ) ^ 3 * Real.log 2 / 2 := by
  have hA0 : (0:ℝ) < A₀ := lt_of_lt_of_le zero_lt_one hA₀
  have hc0 : (0:ℝ) < cS := lt_of_lt_of_le zero_lt_one hcS
  have hlA : (0:ℝ) ≤ Real.log A₀ := Real.log_nonneg hA₀
  have hlZ : (0:ℝ) ≤ Real.log Z₀ := Real.log_nonneg hZ₀
  have hlc : (0:ℝ) ≤ Real.log cS := Real.log_nonneg hcS
  set γ : ℝ := ((1:ℝ) + e) * (4 + 6 * Real.log A₀) + 2 * Real.log Z₀ with hγdef
  set δ : ℝ := 2 * ((1:ℝ) + e) * Real.log cS with hδdef
  have he0 : (0:ℝ) ≤ (e:ℝ) := by positivity
  have hγ : 0 < γ := by rw [hγdef]; nlinarith
  have hδ : 0 ≤ δ := by rw [hδdef]; positivity
  obtain ⟨N, L', hN, hcard, hle⟩ := exists_params γ δ hγ hδ
  refine ⟨N, L', hN, hcard, ?_⟩
  have hNR : (1:ℝ) ≤ (N:ℝ) := by exact_mod_cast hN
  have hNR0 : (0:ℝ) ≤ (N:ℝ) := by linarith
  have hL'1 : 1 ≤ L' := by
    by_contra h
    push_neg at h
    interval_cases L'
    omega
  have hL0 : (0:ℝ) < (L':ℝ) := by exact_mod_cast hL'1
  have hL1R : (1:ℝ) ≤ (L':ℝ) := by exact_mod_cast hL'1
  have hexpand : Real.log (((L':ℝ) ^ 2) * (cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N)))
      = 4 * Real.log (L':ℝ) + 2 * Real.log cS + (6 * (L':ℝ) * (N:ℝ)) * Real.log A₀ := by
    rw [Real.log_mul (by positivity) (by positivity),
        Real.log_mul (by positivity) (by positivity),
        Real.log_mul (by positivity) (by positivity),
        Real.log_pow, Real.log_pow, Real.log_pow]
    push_cast
    ring
  rw [hexpand]
  have hlogL : Real.log (L':ℝ) ≤ (L':ℝ) := by
    have h := Real.log_le_sub_one_of_pos hL0
    linarith
  have hLN : Real.log (L':ℝ) ≤ (L':ℝ) * ((N:ℝ) + 1) := by nlinarith
  have hkey : (L':ℝ) * ((N:ℝ) + 1) * γ + δ
      - (((1:ℝ) + e) * (4 * Real.log (L':ℝ) + 2 * Real.log cS
          + (6 * (L':ℝ) * (N:ℝ)) * Real.log A₀) + (2 * (L':ℝ)) * ((N:ℝ) + 1) * Real.log Z₀)
      = ((1:ℝ) + e) * 4 * ((L':ℝ) * ((N:ℝ) + 1) - Real.log (L':ℝ))
        + ((1:ℝ) + e) * 6 * Real.log A₀ * (L':ℝ) := by
    rw [hγdef, hδdef]
    ring
  have hnn : (0:ℝ) ≤ ((1:ℝ) + e) * 4 * ((L':ℝ) * ((N:ℝ) + 1) - Real.log (L':ℝ))
      + ((1:ℝ) + e) * 6 * Real.log A₀ * (L':ℝ) := by
    have t1 : (0:ℝ) ≤ (L':ℝ) * ((N:ℝ) + 1) - Real.log (L':ℝ) := by linarith
    have t2 : (0:ℝ) ≤ ((1:ℝ) + e) * 4 := by linarith
    have t3 : (0:ℝ) ≤ ((1:ℝ) + e) * 6 * Real.log A₀ * (L':ℝ) := by positivity
    nlinarith
  linarith [hle, hkey, hnn]

/-- ★★★★★**数え上げの完成形** —— `N`・`L'` を選べば、すべての段で
`K₀ · Z₀^{2L'(n+1)} < 2^{n³}` が成り立つ。 -/
theorem exists_params_gap (A₀ Z₀ cS : ℝ) (hA₀ : 1 ≤ A₀) (hZ₀ : 1 ≤ Z₀) (hcS : 1 ≤ cS)
    (e : ℕ) :
    ∃ N L' : ℕ, 0 < N ∧ 2 * N ^ 3 < L' ^ 2 ∧
      ∀ n, N ≤ n →
        (((L':ℝ) ^ 2) * (cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N))) ^ (1 + e)
          * Z₀ ^ ((2 * L') * (n + 1)) < 2 ^ (n ^ 3) := by
  obtain ⟨N, L', hN, hcard, hlog⟩ := exists_params_for_gap A₀ Z₀ cS hA₀ hZ₀ hcS e
  refine ⟨N, L', hN, hcard, ?_⟩
  have hA0 : (0:ℝ) < A₀ := lt_of_lt_of_le zero_lt_one hA₀
  have hc0 : (0:ℝ) < cS := lt_of_lt_of_le zero_lt_one hcS
  have hL'1 : 1 ≤ L' := by
    by_contra h
    push_neg at h
    interval_cases L'
    omega
  have hL1R : (1:ℝ) ≤ (L':ℝ) := by exact_mod_cast hL'1
  set Cb : ℝ := cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N) with hCb
  have h1 : (1:ℝ) ≤ cS ^ 2 := one_le_pow₀ hcS
  have h2 : (1:ℝ) ≤ ((L':ℝ) ^ 2) := one_le_pow₀ hL1R
  have h3 : (1:ℝ) ≤ A₀ ^ (6 * L' * N) := one_le_pow₀ hA₀
  have hCb1 : (1:ℝ) ≤ Cb := by
    rw [hCb]
    calc (1:ℝ) = 1 * 1 * 1 := by ring
      _ ≤ cS ^ 2 * ((L':ℝ) ^ 2) * A₀ ^ (6 * L' * N) := by gcongr
  have hbase1 : (1:ℝ) ≤ ((L':ℝ) ^ 2) * Cb := by
    calc (1:ℝ) = 1 * 1 := by ring
      _ ≤ ((L':ℝ) ^ 2) * Cb := by gcongr
  have hK₀1 : (1:ℝ) ≤ (((L':ℝ) ^ 2) * Cb) ^ (1 + e) := one_le_pow₀ hbase1
  have hlogK : Real.log ((((L':ℝ) ^ 2) * Cb) ^ (1 + e))
      = ((1:ℝ) + e) * Real.log (((L':ℝ) ^ 2) * Cb) := by
    rw [Real.log_pow]
    push_cast
    ring
  refine gap_core hK₀1 hZ₀ hN ?_
  rw [hlogK]
  have hcast : ((2 * L' : ℕ) : ℝ) = 2 * (L':ℝ) := by push_cast; ring
  rw [hcast]
  exact hlog

end ABC3.Found.SixExp
