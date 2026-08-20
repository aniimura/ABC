/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Counting

/-!
# Siegel の補題で初段の零点を作る

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を作るか

超越性証明は「**未定係数の個数を方程式の個数より多く取る**」ことから始まる。
係数 `c_{pq}`(`(p,q) ∈ s`)を未知数、格子点 `m ∈ box` での消滅を方程式とすると、
`|s| > |box|` なら **0 でない代数的整数の解**が取れる —— これが Siegel の補題である。

★mathlib は数体版(house で評価する形)を持っている:
`NumberField.house.exists_ne_zero_int_vec_house_le`。
★ただし評価に現れる定数 `c₁ K` は **private** なので、
`exists_siegel_solution` で存在量化して受け取る。

## ★得られるもの

`exists_siegel_coeffs` は、`extrapolation_induction` の**初段の仮定 `h0`**
(`T 0` で `F` が消える)と、`LatticeHouse.lean` の `hC`(係数の house の上界)を
同時に供給する。

★行列が 0 でないことは **`m = 0` の行が全部 `1`** であることから出る
(`auxMatrixO_zero`)。
-/

namespace ABC3.Found.SixExp

open NumberField Complex

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★1. 整数環の側の `auxMatrix` の道具 -/

omit [NumberField K] in
theorem coe_auxMatrixO (E : Fin 2 → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    ((auxMatrixO E m pq : 𝓞 K) : K)
      = auxMatrixK (fun j i => ((E j i : 𝓞 K) : K)) m pq := by
  simp only [auxMatrixO, auxMatrixK]
  push_cast
  rfl

theorem house_auxMatrixO_le (E : Fin 2 → Fin 3 → 𝓞 K) {A : ℝ} (hA1 : 1 ≤ A)
    (hA : ∀ j i, house ((E j i : 𝓞 K) : K) ≤ A) {L N : ℕ}
    (m : Fin 3 → ℕ) (hm : ∑ i : Fin 3, m i ≤ N) (pq : ℕ × ℕ) (hpq : pq.1 + pq.2 ≤ L) :
    house ((auxMatrixO E m pq : 𝓞 K) : K) ≤ A ^ (L * N) := by
  rw [coe_auxMatrixO]
  refine le_trans (house_auxMatrixK_le _ hA1 hA m pq) ?_
  exact pow_le_pow_right₀ hA1 (Nat.mul_le_mul hpq hm)

omit [NumberField K] in
/-- ★`m = 0` の行はすべて `1` —— 行列が 0 でないことの証人。 -/
theorem auxMatrixO_zero (E : Fin 2 → Fin 3 → 𝓞 K) (pq : ℕ × ℕ) :
    auxMatrixO E (fun _ => 0) pq = 1 := by
  simp [auxMatrixO]

/-! ## ★2. Siegel の補題 -/

/-- ★★**Siegel の補題(数体版)の使いやすい形** ——
mathlib の評価に現れる定数は private なので、存在量化して受け取る。 -/
theorem exists_siegel_solution {α β : Type*} [Fintype α] [Fintype β] [DecidableEq (K →+* ℂ)]
    (a : Matrix α β (𝓞 K)) (ha : a ≠ 0) {p q : ℕ} (h0p : 0 < p) (hpq : p < q)
    (cardβ : Fintype.card β = q) (cardα : Fintype.card α = p)
    {A : ℝ} (habs : ∀ k l, house ((algebraMap (𝓞 K) K) (a k l)) ≤ A) :
    ∃ (C : ℝ) (ξ : β → 𝓞 K), ξ ≠ 0 ∧ a.mulVec ξ = 0
      ∧ ∀ l, house ((ξ l : 𝓞 K) : K) ≤ C := by
  obtain ⟨ξ, hne, hmul, hb⟩ :=
    NumberField.house.exists_ne_zero_int_vec_house_le K a ha h0p hpq cardβ habs cardα
  exact ⟨_, ξ, hne, hmul, hb⟩

/-- ★★★★**Siegel の補題で初段の零点を作る** ——
係数の個数 `|s|` が方程式の個数 `|box|` より多ければ、
**0 でない代数的整数の係数 `c`** が取れて、`box` の各格子点で補助関数が消え、
house は一様に抑えられる。

★これが `extrapolation_induction` の初段の仮定 `h0` と、
`LatticeHouse.lean` の `hC` を同時に供給する。 -/
theorem exists_siegel_coeffs [DecidableEq (K →+* ℂ)] (E : Fin 2 → Fin 3 → 𝓞 K)
    (s : Finset (ℕ × ℕ)) (box : Finset (Fin 3 → ℕ))
    (hzero_mem : (fun _ => 0) ∈ box) (hcard : box.card < s.card)
    {A : ℝ} (hA1 : 1 ≤ A) (hA : ∀ j i, house ((E j i : 𝓞 K) : K) ≤ A)
    {L N : ℕ} (hsL : ∀ pq ∈ s, pq.1 + pq.2 ≤ L)
    (hboxN : ∀ m ∈ box, ∑ i : Fin 3, m i ≤ N) :
    ∃ (C : ℝ) (c : ℕ × ℕ → 𝓞 K), (∃ pq ∈ s, c pq ≠ 0)
      ∧ (∀ m ∈ box, latValO E m s c = 0)
      ∧ (∀ pq, house ((c pq : 𝓞 K) : K) ≤ max C 0) := by
  classical
  set a : Matrix ↥box ↥s (𝓞 K) := fun m pq => auxMatrixO E (m : Fin 3 → ℕ) (pq : ℕ × ℕ)
    with ha_def
  have hbox0 : 0 < box.card := Finset.card_pos.mpr ⟨_, hzero_mem⟩
  have hs0 : s.Nonempty := Finset.card_pos.mp (lt_trans hbox0 hcard)
  obtain ⟨pq0, hpq0⟩ := hs0
  have hane : a ≠ 0 := by
    intro h
    have h1 : a ⟨_, hzero_mem⟩ ⟨pq0, hpq0⟩ = 0 := by rw [h]; rfl
    rw [ha_def] at h1
    simp only at h1
    rw [auxMatrixO_zero] at h1
    exact one_ne_zero h1
  have habs : ∀ (m : ↥box) (pq : ↥s),
      house ((algebraMap (𝓞 K) K) (a m pq)) ≤ A ^ (L * N) := by
    intro m pq
    exact house_auxMatrixO_le E hA1 hA _ (hboxN _ m.2) _ (hsL _ pq.2)
  obtain ⟨C, ξ, hξne, hmul, hbound⟩ :=
    exists_siegel_solution a hane hbox0 hcard (Fintype.card_coe s) (Fintype.card_coe box) habs
  have hzero : house ((0 : 𝓞 K) : K) = 0 := by
    have h := house_intCast (K := K) 0
    simpa using h
  refine ⟨C, fun pq => if h : pq ∈ s then ξ ⟨pq, h⟩ else 0, ?_, ?_, ?_⟩
  · by_contra hc
    push_neg at hc
    refine hξne (funext fun l => ?_)
    have hl := hc (l : ℕ × ℕ) l.2
    simpa using hl
  · intro m hm
    have hm0 : (a.mulVec ξ) ⟨m, hm⟩ = 0 := by rw [hmul]; rfl
    have hm1 : ∑ pq : ↥s, a ⟨m, hm⟩ pq * ξ pq = 0 := by
      rw [← hm0]
      simp [Matrix.mulVec, dotProduct]
    rw [latValO, ← Finset.sum_attach s (fun pq => auxMatrixO E m pq *
      (if h : pq ∈ s then ξ ⟨pq, h⟩ else 0)), ← hm1]
    refine Finset.sum_congr rfl (fun pq _ => ?_)
    simp [ha_def, pq.2]
  · intro pq
    by_cases h : pq ∈ s
    · simp only [h, dif_pos]
      exact le_trans (hbound _) (le_max_left _ _)
    · simp only [h, dif_neg, not_false_iff]
      rw [hzero]
      exact le_max_right _ _

end ABC3.Found.SixExp
