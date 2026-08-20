/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Params

/-!
# Siegel の補題(一般形・明示定数)

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★なぜもう一度書くか

`LatticeGeneral.lean` の `exists_siegel_coeffsG` は係数の house の上界を
**匿名の実数**として返す。★これでは数え上げ(`Params.lean`)に渡せない ——
`log C = O(L·N)` を言うのに上界の**形**が要るからである。

★`SiegelConst.lean` で Siegel の定数に名前を付け、`q ≥ 2p` のとき
`house ≤ c²·q·A` という多項式形を得た。本ファイルはそれを
`Finset` で添字づけた行列(`box` × `s`)に**移す**。

★移し方は `Fintype.equivFinOfCardEq` で `↥box ≃ Fin |box|`、`↥s ≃ Fin |s|` を取り、
行列を引き戻して `siegelC_bound_poly` を当て、解を押し戻すだけである。

## ★得られるもの

  `house(c pq) ≤ c²·|s|·B^N`

★これで数え上げに渡せる: `log` を取ると `2 log c + log|s| + N log B` で、
`B = A₀^{L}` なら `O(L·N)`、`log|s| = log L²` は `O(log L)`。
-/

namespace ABC3.Found.SixExp

open NumberField Complex

variable {K : Type*} [Field K] [NumberField K]

/-- ★★★★**Siegel の補題(一般形・明示定数)** —— `|s| ≥ 2|box|` なら
係数の house は **`c²·|s|·B^N`** で抑えられる。

★`exists_siegel_coeffsG` は定数が匿名だったので数え上げに渡せなかった。
こちらは `SiegelConst.lean` の多項式形を使う。 -/
theorem exists_siegel_coeffsG_poly [DecidableEq (K →+* ℂ)] (G : ℕ × ℕ → Fin 3 → 𝓞 K)
    (s : Finset (ℕ × ℕ)) (box : Finset (Fin 3 → ℕ))
    (hzero_mem : (fun _ => 0) ∈ box) (hcard : 2 * box.card ≤ s.card)
    {B : ℝ} (hB1 : 1 ≤ B) (hB : ∀ pq i, house ((G pq i : 𝓞 K) : K) ≤ B)
    {N : ℕ} (hboxN : ∀ m ∈ box, ∑ i : Fin 3, m i ≤ N) :
    ∃ c : ℕ × ℕ → 𝓞 K, (∃ pq ∈ s, c pq ≠ 0)
      ∧ (∀ m ∈ box, latValG G m s c = 0)
      ∧ (∀ pq, house ((c pq : 𝓞 K) : K) ≤ siegelCpos K ^ 2 * s.card * B ^ N) := by
  classical
  have hbox0 : 0 < box.card := Finset.card_pos.mpr ⟨_, hzero_mem⟩
  have hs0 : 0 < s.card := by omega
  obtain ⟨pq0, hpq0⟩ := Finset.card_pos.mp hs0
  set eα : ↥box ≃ Fin box.card := Fintype.equivFinOfCardEq (Fintype.card_coe box) with heα
  set eβ : ↥s ≃ Fin s.card := Fintype.equivFinOfCardEq (Fintype.card_coe s) with heβ
  set a : Matrix (Fin box.card) (Fin s.card) (𝓞 K) :=
    fun k l => auxMatrixG G ((eα.symm k : ↥box) : Fin 3 → ℕ) ((eβ.symm l : ↥s) : ℕ × ℕ)
    with ha
  have hane : a ≠ 0 := by
    intro h
    have h1 : a (eα ⟨_, hzero_mem⟩) (eβ ⟨pq0, hpq0⟩) = 0 := by rw [h]; rfl
    rw [ha] at h1
    simp only [Equiv.symm_apply_apply] at h1
    rw [auxMatrixG_zero] at h1
    exact one_ne_zero h1
  have habs : ∀ (k : Fin box.card) (l : Fin s.card),
      house ((algebraMap (𝓞 K) K) (a k l)) ≤ B ^ N := fun k l =>
    le_trans (house_auxMatrixG_le G hB1 hB _ _)
      (pow_le_pow_right₀ hB1 (hboxN _ (eα.symm k).2))
  have hBN : (1:ℝ) ≤ B ^ N := one_le_pow₀ hB1
  obtain ⟨ξ, hξne, hmul, hbound⟩ :=
    siegelC_bound_poly box.card s.card hbox0 hcard a hane hBN habs
  have hzero : house ((0 : 𝓞 K) : K) = 0 := by
    have h := house_intCast (K := K) 0
    simpa using h
  set c : ℕ × ℕ → 𝓞 K := fun pq => if h : pq ∈ s then ξ (eβ ⟨pq, h⟩) else 0 with hc
  refine ⟨c, ?_, ?_, ?_⟩
  · by_contra hcc
    push_neg at hcc
    refine hξne (funext fun l => ?_)
    have hl := hcc ((eβ.symm l : ↥s) : ℕ × ℕ) (eβ.symm l).2
    rw [hc] at hl
    simp only [(eβ.symm l).2, dif_pos] at hl
    rw [Subtype.coe_eta, Equiv.apply_symm_apply] at hl
    exact hl
  · intro m hm
    have hm0 : (a.mulVec ξ) (eα ⟨m, hm⟩) = 0 := by rw [hmul]; rfl
    have hm1 : ∑ l : Fin s.card, a (eα ⟨m, hm⟩) l * ξ l = 0 := by
      rw [← hm0]
      simp [Matrix.mulVec, dotProduct]
    have hstep : latValG G m s c = ∑ l : Fin s.card, a (eα ⟨m, hm⟩) l * ξ l := by
      rw [latValG, ← Finset.sum_attach s (fun pq => auxMatrixG G m pq * c pq)]
      exact Fintype.sum_equiv eβ
        (fun t : ↥s => auxMatrixG G m ((t : ↥s) : ℕ × ℕ) * c ((t : ↥s) : ℕ × ℕ))
        (fun l : Fin s.card => a (eα ⟨m, hm⟩) l * ξ l)
        (fun t => by simp only [ha, hc, Equiv.symm_apply_apply, t.2, dif_pos, Subtype.coe_eta])
    rw [hstep]
    exact hm1
  · intro pq
    rw [hc]
    by_cases h : pq ∈ s
    · simp only [h, dif_pos]
      exact hbound _
    · simp only [h, dif_neg, not_false_iff]
      rw [hzero]
      positivity

end ABC3.Found.SixExp
