/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.Siegel

/-!
# 一般形の格子行列と、半径が動く外挿

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★なぜ一般形が要るか

`LatticeHouse.lean` の `auxMatrixO` は `∏ᵢ (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}` という形で、
`E j i` が**代数的整数**であることを要求する。しかし原文の仮定は
`exp(x_j y_i)` が**代数的数**であることだけで、整数とは限らない。

★分母 `D` を払うと `Ê j i := D·E j i ∈ 𝓞_K` になるが、
`∏ᵢ (Ê₀ᵢ^p Ê₁ᵢ^q)^{mᵢ} = D^{(p+q)∑mᵢ}·∏ᵢ (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}` であり、
**指数が `(p,q)` に依存する**ので、係数ベクトルを固定したままでは扱えない。

★★正しい払い方は `D^{L·∑mᵢ}` を掛けることで、そのとき行列の成分は

  `∏ᵢ (Ê₀ᵢ^p · Ê₁ᵢ^q · D^{L−p−q})^{mᵢ}`

になる。★これは「各 `(p,q)` に**底ベクトル** `G pq : Fin 3 → 𝓞_K` を与えて
`∏ᵢ (G pq i)^{mᵢ}` を作る」という一般形に収まる。本ファイルはその一般形を扱う。

★house の評価も一般形のほうが素直になる —— `house(G pq i) ≤ B` なら
`house(auxMatrixG) ≤ B^{∑mᵢ}`。

## ★半径が動く外挿

実際の証明では格子点が段ごとに遠くへ行くので、`r`・`R`・`M` はすべて段に依存する。
`extrapolation_induction'` はその形。
-/

namespace ABC3.Found.SixExp

open Metric Complex NumberField

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★1. 半径が動く外挿 -/

/-- ★★★★★**外挿の帰納(半径が段ごとに動く版)**。

★`Schwarz` の多零点版は「零点も新しい点も `|z| ≤ r k` の中、
上界は `|z| = R k` で `M k`、`5·r k ≤ R k`」という形で使う。 -/
theorem extrapolation_induction' (σ : K →+* ℂ) {F : ℂ → ℂ} (r R M H : ℕ → ℝ)
    (hr : ∀ k, 0 ≤ r k) (h5 : ∀ k, 5 * r k ≤ R k) (hrR : ∀ k, r k < R k)
    (hF : ∀ k, DifferentiableOn ℂ F (closedBall 0 (R k)))
    (hMk : ∀ k, ∀ ζ : ℂ, ‖ζ‖ = R k → ‖F ζ‖ ≤ M k)
    (T : ℕ → Finset ℂ)
    (hTr : ∀ k, ∀ w ∈ T (k + 1), ‖w‖ ≤ r k)
    (hmono : ∀ k, T k ⊆ T (k + 1))
    (h0 : ∀ w ∈ T 0, F w = 0)
    (hH1 : ∀ k, 1 ≤ H k)
    (halg : ∀ k, ∀ w ∈ T (k + 1),
      ∃ α : 𝓞 K, σ ((α : 𝓞 K) : K) = F w ∧ house ((α : 𝓞 K) : K) ≤ H k)
    (hgap : ∀ k, M k * (1/2 : ℝ) ^ (T k).card < (H k ^ (Module.finrank ℚ K - 1))⁻¹) :
    ∀ k, ∀ w ∈ T k, F w = 0 := by
  intro k
  induction k with
  | zero => exact h0
  | succ n ih =>
    intro w hw
    obtain ⟨α, hα, hHle⟩ := halg n w hw
    have hsr : ∀ v ∈ T n, ‖v‖ ≤ r n := fun v hv => hTr n v (hmono n hv)
    exact extrapolation_step σ (hr n) (h5 n) (hrR n) (hF n) hsr ih (hMk n)
      (hTr n w hw) hα (hH1 n) hHle (hgap n)

/-! ## ★2. 一般形の格子行列 -/

/-- ★★**一般化した格子行列** —— 各 `(p,q)` に「底ベクトル」`G pq : Fin 3 → 𝓞 K` を与える。

★`auxMatrixO` は `G pq i = E₀ᵢ^p · E₁ᵢ^q` の場合。
★分母を払うには `G pq i = Ê₀ᵢ^p · Ê₁ᵢ^q · D^{L−p−q}` と取る。 -/
noncomputable def auxMatrixG (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) : 𝓞 K :=
  ∏ i : Fin 3, (G pq i) ^ (m i)

/-- ★一般化した格子点の値 —— **代数的整数である**。 -/
noncomputable def latValG (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → 𝓞 K) : 𝓞 K := ∑ pq ∈ s, auxMatrixG G m pq * c pq

omit [NumberField K] in
/-- ★`m = 0` の行はすべて `1` —— 行列が 0 でないことの証人。 -/
theorem auxMatrixG_zero (G : ℕ × ℕ → Fin 3 → 𝓞 K) (pq : ℕ × ℕ) :
    auxMatrixG G (fun _ => 0) pq = 1 := by
  simp [auxMatrixG]

omit [NumberField K] in
theorem coe_auxMatrixG (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    ((auxMatrixG G m pq : 𝓞 K) : K) = ∏ i : Fin 3, ((G pq i : 𝓞 K) : K) ^ (m i) := by
  simp only [auxMatrixG]
  push_cast
  rfl

omit [NumberField K] in
theorem coe_latValG (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → 𝓞 K) :
    ((latValG G m s c : 𝓞 K) : K)
      = ∑ pq ∈ s, ((auxMatrixG G m pq : 𝓞 K) : K) * ((c pq : 𝓞 K) : K) := by
  simp only [latValG]
  push_cast
  rfl

/-- ★★`auxMatrixG` の house 評価 —— 底が `B` で抑えられれば `B^{∑mᵢ}`。 -/
theorem house_auxMatrixG_le (G : ℕ × ℕ → Fin 3 → 𝓞 K) {B : ℝ} (_hB1 : 1 ≤ B)
    (hB : ∀ pq i, house ((G pq i : 𝓞 K) : K) ≤ B) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    house ((auxMatrixG G m pq : 𝓞 K) : K) ≤ B ^ (∑ i : Fin 3, m i) := by
  rw [coe_auxMatrixG]
  calc house (∏ i : Fin 3, ((G pq i : 𝓞 K) : K) ^ (m i))
      ≤ ∏ i : Fin 3, house (((G pq i : 𝓞 K) : K) ^ (m i)) := house_prod_le' _ _
    _ ≤ ∏ i : Fin 3, B ^ (m i) := by
        refine Finset.prod_le_prod (fun i _ => house_nonneg _) (fun i _ => ?_)
        exact le_trans (house_pow_le _ _) (pow_le_pow_left₀ (house_nonneg _) (hB pq i) _)
    _ = B ^ (∑ i : Fin 3, m i) := Finset.prod_pow_eq_pow_sum _ _ _

/-- ★★★格子点の値の house 評価(一般形)—— 外挿の `H k` を与える。 -/
theorem house_latValG_le (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → 𝓞 K) {B C : ℝ} (hB1 : 1 ≤ B)
    (hB : ∀ pq i, house ((G pq i : 𝓞 K) : K) ≤ B)
    (hC : ∀ pq ∈ s, house ((c pq : 𝓞 K) : K) ≤ C)
    {N : ℕ} (hm : ∑ i : Fin 3, m i ≤ N) :
    house ((latValG G m s c : 𝓞 K) : K) ≤ s.card * (B ^ N * C) := by
  have hB0 : (0:ℝ) ≤ B := le_trans zero_le_one hB1
  rw [coe_latValG]
  have hterm : ∀ pq ∈ s,
      house (((auxMatrixG G m pq : 𝓞 K) : K) * ((c pq : 𝓞 K) : K)) ≤ B ^ N * C := by
    intro pq hpq
    have h1 : house ((auxMatrixG G m pq : 𝓞 K) : K) ≤ B ^ N :=
      le_trans (house_auxMatrixG_le G hB1 hB m pq) (pow_le_pow_right₀ hB1 hm)
    refine le_trans (house_mul_le _ _) ?_
    exact mul_le_mul h1 (hC pq hpq) (house_nonneg _) (by positivity)
  calc house (∑ pq ∈ s, ((auxMatrixG G m pq : 𝓞 K) : K) * ((c pq : 𝓞 K) : K))
      ≤ ∑ pq ∈ s, house (((auxMatrixG G m pq : 𝓞 K) : K) * ((c pq : 𝓞 K) : K)) :=
        house_sum_le_sum_house _ _
    _ ≤ ∑ _pq ∈ s, B ^ N * C := Finset.sum_le_sum hterm
    _ = s.card * (B ^ N * C) := by rw [Finset.sum_const, nsmul_eq_mul]

omit [NumberField K] in
/-- ★★一般形の格子点の値を `σ` で送った姿。 -/
theorem map_latValG (σ : K →+* ℂ) (G : ℕ × ℕ → Fin 3 → 𝓞 K) (m : Fin 3 → ℕ)
    (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → 𝓞 K) :
    σ ((latValG G m s c : 𝓞 K) : K)
      = ∑ pq ∈ s, (∏ i : Fin 3, (σ ((G pq i : 𝓞 K) : K)) ^ (m i)) * σ ((c pq : 𝓞 K) : K) := by
  rw [coe_latValG, map_sum]
  refine Finset.sum_congr rfl (fun pq _ => ?_)
  rw [map_mul, coe_auxMatrixG, map_prod]
  congr 1
  exact Finset.prod_congr rfl (fun i _ => by rw [map_pow])

/-! ## ★3. Siegel の補題(一般形) -/

/-- ★★★★**Siegel の補題(一般形)** —— 係数の個数 `|s|` が方程式の個数 `|box|` より
多ければ、**0 でない代数的整数の係数**が取れて `box` の各格子点で値が消える。 -/
theorem exists_siegel_coeffsG [DecidableEq (K →+* ℂ)] (G : ℕ × ℕ → Fin 3 → 𝓞 K)
    (s : Finset (ℕ × ℕ)) (box : Finset (Fin 3 → ℕ))
    (hzero_mem : (fun _ => 0) ∈ box) (hcard : box.card < s.card)
    {B : ℝ} (hB1 : 1 ≤ B) (hB : ∀ pq i, house ((G pq i : 𝓞 K) : K) ≤ B)
    {N : ℕ} (hboxN : ∀ m ∈ box, ∑ i : Fin 3, m i ≤ N) :
    ∃ (C : ℝ) (c : ℕ × ℕ → 𝓞 K), (∃ pq ∈ s, c pq ≠ 0)
      ∧ (∀ m ∈ box, latValG G m s c = 0)
      ∧ (∀ pq, house ((c pq : 𝓞 K) : K) ≤ max C 0) := by
  classical
  set a : Matrix ↥box ↥s (𝓞 K) := fun m pq => auxMatrixG G (m : Fin 3 → ℕ) (pq : ℕ × ℕ)
    with ha_def
  have hbox0 : 0 < box.card := Finset.card_pos.mpr ⟨_, hzero_mem⟩
  have hs0 : s.Nonempty := Finset.card_pos.mp (lt_trans hbox0 hcard)
  obtain ⟨pq0, hpq0⟩ := hs0
  have hane : a ≠ 0 := by
    intro h
    have h1 : a ⟨_, hzero_mem⟩ ⟨pq0, hpq0⟩ = 0 := by rw [h]; rfl
    rw [ha_def] at h1
    simp only at h1
    rw [auxMatrixG_zero] at h1
    exact one_ne_zero h1
  have habs : ∀ (m : ↥box) (pq : ↥s),
      house ((algebraMap (𝓞 K) K) (a m pq)) ≤ B ^ N := fun m pq =>
    le_trans (house_auxMatrixG_le G hB1 hB _ _) (pow_le_pow_right₀ hB1 (hboxN _ m.2))
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
    rw [latValG, ← Finset.sum_attach s (fun pq => auxMatrixG G m pq *
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
