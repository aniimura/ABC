/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Extrapolation

/-!
# 格子点での値の house 評価 —— 外挿の `halg` を供給する

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を作るか

`Extrapolation.lean` の `extrapolation_induction` は 2 種類の入力を要求する:

* **`halg`** —— 各格子点での値が、house が `H k` で抑えられた代数的整数の像であること
* **`hgap`** —— 解析側の上界が Liouville の下界を下回ること(パラメータの数え上げ)

★本ファイルは **`halg` を完全に供給する**。

## ★段取り

`E j i ∈ 𝓞_K` を `exp(x_j y_i)`(代数的整数と仮定)とし、係数 `c ∈ 𝓞_K` とすると、

  `latValO E m s c = ∑_{(p,q) ∈ s} (∏_i (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}) · c_{pq} ∈ 𝓞_K`

は**代数的整数**であり、`σ` で送ると補助関数の格子点での値そのものになる
(`map_latValK`)。house は

  `house ≤ |s| · A^{L·N} · C`   (`house_latValK_le`)

で抑えられる(`house(E j i) ≤ A`、`p+q ≤ L`、`∑ mᵢ ≤ N`、`house(c) ≤ C`)。

★mathlib の house の劣加法性・劣乗法性(`house_sum_le_sum_house` / `house_mul_le` /
`house_pow_le`)がそのまま効く。添字つき積の版だけ足した(`house_prod_le'`)。
-/

namespace ABC3.Found.SixExp

open NumberField Complex

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★1. house の補助 -/

theorem house_one : house (1 : K) = 1 := by
  have h := house_intCast (K := K) 1
  simpa using h

/-- ★添字つき積についての house の劣乗法性。 -/
theorem house_prod_le' {ι : Type*} (t : Finset ι) (f : ι → K) :
    house (∏ i ∈ t, f i) ≤ ∏ i ∈ t, house (f i) := by
  classical
  induction t using Finset.induction_on with
  | empty => simp [house_one]
  | insert a t ha ih =>
    rw [Finset.prod_insert ha, Finset.prod_insert ha]
    refine le_trans (house_mul_le _ _) ?_
    gcongr
    exact house_nonneg _

/-! ## ★2. 数体の側の格子点の値 -/

/-- ★数体の側の `auxMatrix`。 -/
noncomputable def auxMatrixK (E : Fin 2 → Fin 3 → K) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) : K :=
  ∏ i : Fin 3, (E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i)

/-- ★★`auxMatrix` の house 評価 —— 指数は `(p+q)·∑mᵢ`。 -/
theorem house_auxMatrixK_le (E : Fin 2 → Fin 3 → K) {A : ℝ} (hA1 : 1 ≤ A)
    (hA : ∀ j i, house (E j i) ≤ A) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    house (auxMatrixK E m pq) ≤ A ^ ((pq.1 + pq.2) * ∑ i : Fin 3, m i) := by
  have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one hA1
  have hstep : ∀ i : Fin 3, house ((E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i))
      ≤ A ^ ((pq.1 + pq.2) * m i) := by
    intro i
    have h1 : house (E 0 i ^ pq.1 * E 1 i ^ pq.2) ≤ A ^ pq.1 * A ^ pq.2 := by
      refine le_trans (house_mul_le _ _) ?_
      have h01 : house (E 0 i ^ pq.1) ≤ A ^ pq.1 :=
        le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) (hA 0 i) _)
      have h02 : house (E 1 i ^ pq.2) ≤ A ^ pq.2 :=
        le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) (hA 1 i) _)
      exact mul_le_mul h01 h02 (house_nonneg _) (by positivity)
    calc house ((E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i))
        ≤ (house (E 0 i ^ pq.1 * E 1 i ^ pq.2)) ^ (m i) := house_pow_le _ _
      _ ≤ (A ^ pq.1 * A ^ pq.2) ^ (m i) := pow_le_pow_left₀ (house_nonneg _) h1 _
      _ = A ^ ((pq.1 + pq.2) * m i) := by rw [← pow_add, ← pow_mul]
  calc house (auxMatrixK E m pq)
      ≤ ∏ i : Fin 3, house ((E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i)) := house_prod_le' _ _
    _ ≤ ∏ i : Fin 3, A ^ ((pq.1 + pq.2) * m i) :=
        Finset.prod_le_prod (fun i _ => house_nonneg _) (fun i _ => hstep i)
    _ = A ^ (∑ i : Fin 3, (pq.1 + pq.2) * m i) := Finset.prod_pow_eq_pow_sum _ _ _
    _ = A ^ ((pq.1 + pq.2) * ∑ i : Fin 3, m i) := by rw [← Finset.mul_sum]

/-- ★数体の側の格子点での値。 -/
noncomputable def latValK (E : Fin 2 → Fin 3 → K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → K) : K := ∑ pq ∈ s, auxMatrixK E m pq * c pq

/-- ★★★**格子点での値の house 評価** —— 外挿の `H k` を与える段。

`house(E j i) ≤ A`、`p+q ≤ L`、`∑ mᵢ ≤ N`、`house(c) ≤ C` なら
`house ≤ |s| · A^{L·N} · C`。 -/
theorem house_latValK_le (E : Fin 2 → Fin 3 → K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → K) {A C : ℝ} (hA1 : 1 ≤ A) (hA : ∀ j i, house (E j i) ≤ A)
    (hC : ∀ pq ∈ s, house (c pq) ≤ C)
    {L N : ℕ} (hsL : ∀ pq ∈ s, pq.1 + pq.2 ≤ L) (hm : ∑ i : Fin 3, m i ≤ N) :
    house (latValK E m s c) ≤ s.card * (A ^ (L * N) * C) := by
  have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one hA1
  have hterm : ∀ pq ∈ s, house (auxMatrixK E m pq * c pq) ≤ A ^ (L * N) * C := by
    intro pq hpq
    have h1 : house (auxMatrixK E m pq) ≤ A ^ (L * N) := by
      refine le_trans (house_auxMatrixK_le E hA1 hA m pq) ?_
      exact pow_le_pow_right₀ hA1 (Nat.mul_le_mul (hsL pq hpq) hm)
    have h2 : house (c pq) ≤ C := hC pq hpq
    refine le_trans (house_mul_le _ _) ?_
    exact mul_le_mul h1 h2 (house_nonneg _) (by positivity)
  calc house (latValK E m s c)
      ≤ ∑ pq ∈ s, house (auxMatrixK E m pq * c pq) := house_sum_le_sum_house _ _
    _ ≤ ∑ _pq ∈ s, A ^ (L * N) * C := Finset.sum_le_sum hterm
    _ = s.card * (A ^ (L * N) * C) := by rw [Finset.sum_const, nsmul_eq_mul]

/-- ★★数体の側の格子点の値を `σ` で送ると、解析側の補助関数の値になる。 -/
theorem map_latValK (σ : K →+* ℂ) (E : Fin 2 → Fin 3 → K) (x : Fin 2 → ℂ) (y : Fin 3 → ℂ)
    (hE : ∀ j i, σ (E j i) = sixExpVals x y j i) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → K) :
    σ (latValK E m s c) = auxFun x s (fun pq => σ (c pq)) (latticePt y m) := by
  rw [auxFun_latticePt, latValK, map_sum]
  refine Finset.sum_congr rfl (fun pq _ => ?_)
  rw [map_mul]
  congr 1
  rw [auxMatrixK, auxMatrix, map_prod]
  refine Finset.prod_congr rfl (fun i _ => ?_)
  rw [map_pow, map_mul, map_pow, map_pow, hE, hE]

/-! ## ★3. 整数環の側 -/

/-- ★整数環の側の `auxMatrix`。 -/
noncomputable def auxMatrixO (E : Fin 2 → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) : 𝓞 K :=
  ∏ i : Fin 3, (E 0 i ^ pq.1 * E 1 i ^ pq.2) ^ (m i)

/-- ★整数環の側の格子点の値 —— **代数的整数である**。 -/
noncomputable def latValO (E : Fin 2 → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → 𝓞 K) : 𝓞 K := ∑ pq ∈ s, auxMatrixO E m pq * c pq

omit [NumberField K] in
theorem coe_latValO (E : Fin 2 → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → 𝓞 K) :
    ((latValO E m s c : 𝓞 K) : K)
      = latValK (fun j i => ((E j i : 𝓞 K) : K)) m s (fun pq => ((c pq : 𝓞 K) : K)) := by
  simp only [latValO, latValK, auxMatrixO, auxMatrixK]
  push_cast
  rfl

/-- ★★★★**外挿の `halg` を供給する** —— 格子点での値は
**house が明示的に抑えられた代数的整数**である。

★これで `extrapolation_induction` の 2 つの入力のうち片方が完全に埋まった。
残るのは `hgap`(パラメータの数え上げ)である。 -/
theorem exists_algInt_latticePt (σ : K →+* ℂ) (E : Fin 2 → Fin 3 → 𝓞 K)
    (x : Fin 2 → ℂ) (y : Fin 3 → ℂ)
    (hE : ∀ j i, σ ((E j i : 𝓞 K) : K) = sixExpVals x y j i)
    (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → 𝓞 K)
    {A C : ℝ} (hA1 : 1 ≤ A) (hA : ∀ j i, house ((E j i : 𝓞 K) : K) ≤ A)
    (hC : ∀ pq ∈ s, house ((c pq : 𝓞 K) : K) ≤ C)
    {L N : ℕ} (hsL : ∀ pq ∈ s, pq.1 + pq.2 ≤ L)
    (m : Fin 3 → ℕ) (hm : ∑ i : Fin 3, m i ≤ N) :
    ∃ α : 𝓞 K, σ ((α : 𝓞 K) : K)
        = auxFun x s (fun pq => σ ((c pq : 𝓞 K) : K)) (latticePt y m)
      ∧ house ((α : 𝓞 K) : K) ≤ s.card * (A ^ (L * N) * C) := by
  refine ⟨latValO E m s c, ?_, ?_⟩
  · rw [coe_latValO]
    exact map_latValK σ _ x y hE m s _
  · rw [coe_latValO]
    exact house_latValK_le _ m s _ hA1 hA hC hsL hm

end ABC3.Found.SixExp
