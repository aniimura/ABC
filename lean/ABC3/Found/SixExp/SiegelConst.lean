/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.DenomClear

/-!
# Siegel の定数を名前で取り出す

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何が問題だったか

数え上げ(`Counting.lean`)を回すには、Siegel の補題が与える係数の house の上界

  `C ≈ c₁·(c₁·q·A)^{p/(q−p)}`

を**明示的に**押さえて `log C = O(L·N)` を言う必要がある。
★しかし mathlib の `c₁ K` は **private** で、外から名前で呼べない。

## ★取り出し方

`c₁ K` は **`p`・`q`・`A`・行列に依らない**。したがって

  `⟨_, fun p q … => NumberField.house.exists_ne_zero_int_vec_house_le …⟩`

と書けば、★**`_` が単一化で `c₁ K` に定まる**(`exists_siegel_const`)。
あとは `Exists.choose` で名前を付ければよい(`siegelC`)。

★符号は分からない(`c₁ K` の非負性は private な補題)ので、
`|c₁ K|` に丸め(`Real.abs_rpow_le_abs_rpow` で負の底の `rpow` を避ける)、
さらに `1` 以上に丸める(`siegelCpos`)。

## ★使える形

★`q ≥ 2p` なら指数 `p/(q−p) ≤ 1` なので、`x ≥ 1` に対し `x^e ≤ x` で

  **`house(ξ l) ≤ c²·q·A`**   (`siegelC_bound_poly`)

という多項式の形になる。★これが数え上げに渡せる形である ——
`log C ≤ 2 log c + log q + log A` で、`A` が `A₀^{L·N}` なら `O(L·N)`。
-/

namespace ABC3.Found.SixExp

open NumberField Complex

variable {K : Type*} [Field K] [NumberField K]

/-- ★★**Siegel の定数を名前で取り出す** ——
mathlib の評価に現れる定数 `c₁ K` は private だが、
**`p`・`q`・`A`・行列に依らない**ので、存在量化の証人として一度だけ取り出せる。

★`⟨_, …⟩` の `_` が単一化で `c₁ K` に定まるのが要点。 -/
theorem exists_siegel_const [DecidableEq (K →+* ℂ)] :
    ∃ c : ℝ, ∀ (p q : ℕ), 0 < p → p < q → ∀ (a : Matrix (Fin p) (Fin q) (𝓞 K)), a ≠ 0 →
      ∀ (A : ℝ), (∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) →
      ∃ ξ : Fin q → 𝓞 K, ξ ≠ 0 ∧ a.mulVec ξ = 0 ∧
        ∀ l, house ((ξ l : 𝓞 K) : K) ≤ c * ((c * q * A) ^ ((p : ℝ) / (q - p))) :=
  ⟨_, fun p q hp hpq a ha _A habs =>
    NumberField.house.exists_ne_zero_int_vec_house_le K a ha hp hpq (Fintype.card_fin q) habs
      (Fintype.card_fin p)⟩

/-- ★Siegel の定数(名前つき)。 -/
noncomputable def siegelC (K : Type*) [Field K] [NumberField K] [DecidableEq (K →+* ℂ)] : ℝ :=
  (exists_siegel_const (K := K)).choose

theorem siegelC_spec [DecidableEq (K →+* ℂ)] (p q : ℕ) (hp : 0 < p) (hpq : p < q)
    (a : Matrix (Fin p) (Fin q) (𝓞 K)) (ha : a ≠ 0) (A : ℝ)
    (habs : ∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) :
    ∃ ξ : Fin q → 𝓞 K, ξ ≠ 0 ∧ a.mulVec ξ = 0 ∧
      ∀ l, house ((ξ l : 𝓞 K) : K)
        ≤ siegelC K * ((siegelC K * q * A) ^ ((p : ℝ) / (q - p))) :=
  (exists_siegel_const (K := K)).choose_spec p q hp hpq a ha A habs

/-- ★★**符号を気にしない上界** —— `|c|` で置き換えた形。

★`|x^y| ≤ |x|^y`(`Real.abs_rpow_le_abs_rpow`)で、負の底の `rpow` の面倒を避ける。 -/
theorem siegelC_spec_abs [DecidableEq (K →+* ℂ)] (p q : ℕ) (hp : 0 < p) (hpq : p < q)
    (a : Matrix (Fin p) (Fin q) (𝓞 K)) (ha : a ≠ 0) {A : ℝ} (hA0 : 0 ≤ A)
    (habs : ∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) :
    ∃ ξ : Fin q → 𝓞 K, ξ ≠ 0 ∧ a.mulVec ξ = 0 ∧
      ∀ l, house ((ξ l : 𝓞 K) : K)
        ≤ |siegelC K| * ((|siegelC K| * q * A) ^ ((p : ℝ) / (q - p))) := by
  obtain ⟨ξ, hne, hmul, hb⟩ := siegelC_spec p q hp hpq a ha A habs
  refine ⟨ξ, hne, hmul, fun l => le_trans (hb l) ?_⟩
  set c := siegelC K with hc
  set e := ((p : ℝ) / (q - p)) with he
  have habs1 : |c * ((c * q * A) ^ e)| ≤ |c| * ((|c| * q * A) ^ e) := by
    rw [abs_mul]
    refine mul_le_mul_of_nonneg_left ?_ (abs_nonneg c)
    refine le_trans (Real.abs_rpow_le_abs_rpow _ _) ?_
    have habs2 : |c * (q : ℝ) * A| = |c| * (q : ℝ) * A := by
      rw [abs_mul, abs_mul, abs_of_nonneg (by positivity : (0:ℝ) ≤ (q:ℝ)), abs_of_nonneg hA0]
    rw [habs2]
  exact le_trans (le_abs_self _) habs1

/-- ★`1` 以上に丸めた Siegel の定数。 -/
noncomputable def siegelCpos (K : Type*) [Field K] [NumberField K] [DecidableEq (K →+* ℂ)] : ℝ :=
  max 1 |siegelC K|

theorem one_le_siegelCpos [DecidableEq (K →+* ℂ)] : 1 ≤ siegelCpos K := le_max_left _ _

theorem siegelC_spec_pos [DecidableEq (K →+* ℂ)] (p q : ℕ) (hp : 0 < p) (hpq : p < q)
    (a : Matrix (Fin p) (Fin q) (𝓞 K)) (ha : a ≠ 0) {A : ℝ} (hA0 : 0 ≤ A)
    (habs : ∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) :
    ∃ ξ : Fin q → 𝓞 K, ξ ≠ 0 ∧ a.mulVec ξ = 0 ∧
      ∀ l, house ((ξ l : 𝓞 K) : K)
        ≤ siegelCpos K * ((siegelCpos K * q * A) ^ ((p : ℝ) / (q - p))) := by
  obtain ⟨ξ, hne, hmul, hb⟩ := siegelC_spec_abs p q hp hpq a ha hA0 habs
  refine ⟨ξ, hne, hmul, fun l => le_trans (hb l) ?_⟩
  have hle : |siegelC K| ≤ siegelCpos K := le_max_right _ _
  have he0 : (0:ℝ) ≤ (p : ℝ) / (q - p) := by
    have hlt : (0:ℝ) < (q : ℝ) - (p : ℝ) := by
      have : (p:ℝ) < (q:ℝ) := by exact_mod_cast hpq
      linarith
    positivity
  have hb0 : (0:ℝ) ≤ |siegelC K| * (q : ℝ) * A := by positivity
  refine mul_le_mul hle (Real.rpow_le_rpow hb0 ?_ he0) (Real.rpow_nonneg hb0 _)
    (le_trans (abs_nonneg _) hle)
  have hq0 : (0:ℝ) ≤ (q:ℝ) := by positivity
  gcongr

/-- ★★★**多項式の形の Siegel 評価** —— `q ≥ 2p` なら指数は `≤ 1` なので、
`house(ξ l) ≤ c²·q·A` という**扱いやすい形**になる。

★これが数え上げ(`Counting.lean`)に渡せる形である ——
`log C ≤ 2 log c + log q + log A` で、`A` が `A₀^{L·N}` なら `O(L·N)`。 -/
theorem siegelC_bound_poly [DecidableEq (K →+* ℂ)] (p q : ℕ) (hp : 0 < p) (hpq : 2 * p ≤ q)
    (a : Matrix (Fin p) (Fin q) (𝓞 K)) (ha : a ≠ 0) {A : ℝ} (hA1 : 1 ≤ A)
    (habs : ∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) :
    ∃ ξ : Fin q → 𝓞 K, ξ ≠ 0 ∧ a.mulVec ξ = 0 ∧
      ∀ l, house ((ξ l : 𝓞 K) : K) ≤ siegelCpos K ^ 2 * q * A := by
  have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one hA1
  have hpq' : p < q := by omega
  obtain ⟨ξ, hne, hmul, hb⟩ := siegelC_spec_pos p q hp hpq' a ha hA0 habs
  refine ⟨ξ, hne, hmul, fun l => le_trans (hb l) ?_⟩
  set c := siegelCpos K with hc
  have hc1 : 1 ≤ c := one_le_siegelCpos
  have hqR : (1:ℝ) ≤ (q : ℝ) := by
    have h1 : 1 ≤ q := by omega
    exact_mod_cast h1
  have hx1 : (1:ℝ) ≤ c * (q : ℝ) * A := by
    have h1 : (1:ℝ) * 1 * 1 ≤ c * (q:ℝ) * A := by
      refine mul_le_mul (mul_le_mul hc1 hqR zero_le_one (le_trans zero_le_one hc1)) hA1
        zero_le_one ?_
      positivity
    simpa using h1
  have he1 : (p : ℝ) / ((q : ℝ) - (p : ℝ)) ≤ 1 := by
    have hpR : (p:ℝ) ≤ (q:ℝ) - (p:ℝ) := by
      have h2 : (2 * p : ℝ) ≤ (q : ℝ) := by exact_mod_cast hpq
      linarith
    have hpos : (0:ℝ) < (q:ℝ) - (p:ℝ) := by
      have h3 : (p:ℝ) < (q:ℝ) := by exact_mod_cast hpq'
      linarith
    rw [div_le_one hpos]
    exact hpR
  calc c * ((c * (q:ℝ) * A) ^ ((p : ℝ) / ((q:ℝ) - (p:ℝ))))
      ≤ c * ((c * (q:ℝ) * A) ^ (1:ℝ)) := by
        refine mul_le_mul_of_nonneg_left ?_ (le_trans zero_le_one hc1)
        exact Real.rpow_le_rpow_of_exponent_le hx1 he1
    _ = c * (c * (q:ℝ) * A) := by rw [Real.rpow_one]
    _ = c ^ 2 * (q:ℝ) * A := by ring

end ABC3.Found.SixExp
