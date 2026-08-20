/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Assembly

/-!
# 分母を払う —— 代数的数から代数的整数へ

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何が問題か

原文の仮定は `exp(x_j y_i)` が**代数的数**であることだけで、代数的**整数**とは限らない。
一方 Liouville 不等式(`Liouville.lean`)は `𝓞_K` の元にしか使えない。

★共通分母 `b ∈ 𝓞_K`(`b·E j i ∈ 𝓞_K`)を取っても、素朴に払うと

  `∏ᵢ ((b·E₀ᵢ)^p (b·E₁ᵢ)^q)^{mᵢ} = b^{(p+q)∑mᵢ}·∏ᵢ (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}`

で**指数が `(p,q)` に依存する**ので、係数ベクトルを固定したままでは扱えない。

## ★正しい払い方

★`b^{L−p−q}` を掛けて **`b^{L·∑mᵢ}` に揃える**:

  `denomBase b Ê L pq i := Ê₀ᵢ^p · Ê₁ᵢ^q · b^{L−p−q} ∈ 𝓞_K`  (`p+q ≤ L`)

これは `LatticeGeneral.lean` の一般形 `auxMatrixG` の底ベクトルとしてそのまま使え、

* `map_latValG_denomBase` —— ★`σ(latValG) = β^{L·∑mᵢ} · F(格子点)`(`β = σ b`)。
  ★これが `Assembly.lean` の**倍率つき外挿**の `δ` を与える。
* `house_denomBase_le` —— `house(Ê j i) ≤ A` かつ `house(b) ≤ A`(`A ≥ 1`)なら
  **`house(denomBase) ≤ A^L`**。指数がちょうど `L` に揃うのが要点。
-/

namespace ABC3.Found.SixExp

open Complex NumberField

variable {K : Type*} [Field K] [NumberField K]

theorem auxMatrix_eq_prod_latVec (E : Fin 2 → Fin 3 → ℂ) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    auxMatrix E m pq = ∏ i : Fin 3, (latVec E pq i) ^ (m i) := rfl

/-- ★★**分母を払った底ベクトル** `Ê₀ᵢ^p · Ê₁ᵢ^q · b^{L−p−q}`。 -/
noncomputable def denomBase (b : 𝓞 K) (Eh : Fin 2 → Fin 3 → 𝓞 K) (L : ℕ)
    (pq : ℕ × ℕ) (i : Fin 3) : 𝓞 K :=
  Eh 0 i ^ pq.1 * Eh 1 i ^ pq.2 * b ^ (L - pq.1 - pq.2)

omit [NumberField K] in
theorem map_denomBase (σ : K →+* ℂ) (b : 𝓞 K) (Eh : Fin 2 → Fin 3 → 𝓞 K)
    (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (L : ℕ)
    (hE : ∀ j i, σ ((Eh j i : 𝓞 K) : K) = σ ((b : 𝓞 K) : K) * sixExpVals x y j i)
    {pq : ℕ × ℕ} (hpq : pq.1 + pq.2 ≤ L) (i : Fin 3) :
    σ ((denomBase b Eh L pq i : 𝓞 K) : K)
      = (σ ((b : 𝓞 K) : K)) ^ L * latVec (sixExpVals x y) pq i := by
  set β := σ ((b : 𝓞 K) : K) with hb
  have hcoe : ((denomBase b Eh L pq i : 𝓞 K) : K)
      = ((Eh 0 i : 𝓞 K) : K) ^ pq.1 * ((Eh 1 i : 𝓞 K) : K) ^ pq.2
        * ((b : 𝓞 K) : K) ^ (L - pq.1 - pq.2) := by
    simp only [denomBase]
    push_cast
    rfl
  rw [hcoe, map_mul, map_mul, map_pow, map_pow, map_pow, hE, hE, latVec]
  have hL : pq.1 + pq.2 + (L - pq.1 - pq.2) = L := by omega
  calc (β * sixExpVals x y 0 i) ^ pq.1 * (β * sixExpVals x y 1 i) ^ pq.2 * β ^ (L - pq.1 - pq.2)
      = β ^ (pq.1 + pq.2 + (L - pq.1 - pq.2))
        * (sixExpVals x y 0 i ^ pq.1 * sixExpVals x y 1 i ^ pq.2) := by
        rw [mul_pow, mul_pow, pow_add, pow_add]
        ring
    _ = β ^ L * (sixExpVals x y 0 i ^ pq.1 * sixExpVals x y 1 i ^ pq.2) := by rw [hL]

/-- ★★★★**分母を払った格子点の値の `σ` 像** ——
`σ(latValG) = β^{L·∑mᵢ} · F(格子点)`。★これが倍率つき外挿の `δ` を与える。 -/
theorem map_latValG_denomBase (σ : K →+* ℂ) (b : 𝓞 K) (Eh : Fin 2 → Fin 3 → 𝓞 K)
    (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (L : ℕ)
    (hE : ∀ j i, σ ((Eh j i : 𝓞 K) : K) = σ ((b : 𝓞 K) : K) * sixExpVals x y j i)
    (s : Finset (ℕ × ℕ)) (hsL : ∀ pq ∈ s, pq.1 + pq.2 ≤ L) (c : ℕ × ℕ → 𝓞 K)
    (m : Fin 3 → ℕ) :
    σ ((latValG (denomBase b Eh L) m s c : 𝓞 K) : K)
      = (σ ((b : 𝓞 K) : K)) ^ (L * ∑ i : Fin 3, m i)
        * auxFun x s (fun pq => σ ((c pq : 𝓞 K) : K)) (latticePt y m) := by
  set β := σ ((b : 𝓞 K) : K) with hb
  rw [map_latValG, auxFun_latticePt, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun pq hpq => ?_)
  have hprod : (∏ i : Fin 3, (σ ((denomBase b Eh L pq i : 𝓞 K) : K)) ^ (m i))
      = β ^ (L * ∑ i : Fin 3, m i) * auxMatrix (sixExpVals x y) m pq := by
    have hstep : ∀ i : Fin 3, (σ ((denomBase b Eh L pq i : 𝓞 K) : K)) ^ (m i)
        = (β ^ L) ^ (m i) * (latVec (sixExpVals x y) pq i) ^ (m i) := by
      intro i
      rw [map_denomBase σ b Eh x y L hE (hsL pq hpq) i, mul_pow]
    rw [Finset.prod_congr rfl (fun i _ => hstep i), Finset.prod_mul_distrib,
      auxMatrix_eq_prod_latVec]
    congr 1
    rw [Finset.prod_pow_eq_pow_sum, ← pow_mul]
  rw [hprod]
  ring

/-- ★★分母を払った底ベクトルの house は **`A^L`** で抑えられる。

★指数がちょうど `L` に揃うのが要点(`p + q + (L−p−q) = L`)。 -/
theorem house_denomBase_le (b : 𝓞 K) (Eh : Fin 2 → Fin 3 → 𝓞 K) (L : ℕ) {A : ℝ}
    (hA1 : 1 ≤ A) (hAb : house ((b : 𝓞 K) : K) ≤ A)
    (hAE : ∀ j i, house ((Eh j i : 𝓞 K) : K) ≤ A)
    {pq : ℕ × ℕ} (hpq : pq.1 + pq.2 ≤ L) (i : Fin 3) :
    house ((denomBase b Eh L pq i : 𝓞 K) : K) ≤ A ^ L := by
  have hA0 : (0:ℝ) ≤ A := le_trans zero_le_one hA1
  have hcoe : ((denomBase b Eh L pq i : 𝓞 K) : K)
      = ((Eh 0 i : 𝓞 K) : K) ^ pq.1 * ((Eh 1 i : 𝓞 K) : K) ^ pq.2
        * ((b : 𝓞 K) : K) ^ (L - pq.1 - pq.2) := by
    simp only [denomBase]
    push_cast
    rfl
  rw [hcoe]
  have h1 : house (((Eh 0 i : 𝓞 K) : K) ^ pq.1) ≤ A ^ pq.1 :=
    le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) (hAE 0 i) _)
  have h2 : house (((Eh 1 i : 𝓞 K) : K) ^ pq.2) ≤ A ^ pq.2 :=
    le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) (hAE 1 i) _)
  have h3 : house (((b : 𝓞 K) : K) ^ (L - pq.1 - pq.2)) ≤ A ^ (L - pq.1 - pq.2) :=
    le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) hAb _)
  have h12 : house (((Eh 0 i : 𝓞 K) : K) ^ pq.1 * ((Eh 1 i : 𝓞 K) : K) ^ pq.2)
      ≤ A ^ pq.1 * A ^ pq.2 :=
    le_trans (house_mul_le _ _) (mul_le_mul h1 h2 (house_nonneg _) (by positivity))
  have hL : pq.1 + pq.2 + (L - pq.1 - pq.2) = L := by omega
  calc house (((Eh 0 i : 𝓞 K) : K) ^ pq.1 * ((Eh 1 i : 𝓞 K) : K) ^ pq.2
        * ((b : 𝓞 K) : K) ^ (L - pq.1 - pq.2))
      ≤ house (((Eh 0 i : 𝓞 K) : K) ^ pq.1 * ((Eh 1 i : 𝓞 K) : K) ^ pq.2)
        * house (((b : 𝓞 K) : K) ^ (L - pq.1 - pq.2)) := house_mul_le _ _
    _ ≤ (A ^ pq.1 * A ^ pq.2) * A ^ (L - pq.1 - pq.2) :=
        mul_le_mul h12 h3 (house_nonneg _) (by positivity)
    _ = A ^ (pq.1 + pq.2 + (L - pq.1 - pq.2)) := by rw [pow_add, pow_add]
    _ = A ^ L := by rw [hL]

end ABC3.Found.SixExp
