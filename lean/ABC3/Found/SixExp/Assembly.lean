/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.LatticeBox

/-!
# 組み立てに要る 2 群

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★1. 倍率つき外挿

分母を払うと、格子点での値は `F(w)` そのものではなく **`D^{L·∑mᵢ}·F(w)`** の形の
代数的整数になる。したがって外挿の 1 段も

  `σ α = δ · F w`、`‖δ‖ ≤ Δ`、`δ ≠ 0`

の形で使えなければならない。★`Δ` は解析側の上界に掛かるだけなので、
`hgap` が `Δ·M·(1/2)^{|T|} < H^{-(d-1)}` になる。
`Counting.lean` の数え上げ(3 乗は 1 乗に勝つ)は `Δ` が `N` の**指数の 1 乗**で
増えても通る。

## ★2. 指標の相異性

`coeff_eq_zero_of_auxFun_latticePt`(背理法の閉じ目)は
`Set.InjOn (latVec (sixExpVals x y)) s` を要求する。★これは仮定から出る:

* `latVec (sixExpVals x y) pq i = exp((p·x₀ + q·x₁)·yᵢ)`(`latVec_sixExpVals`)
* `x` が ℚ 上一次独立 ⟹ `(p,q) ↦ p·x₀ + q·x₁` は単射(`auxExp_injective`)
* ★`latVec pq = latVec pq'` なら `λ := 差 ≠ 0` について `exp(λ·yᵢ) = 1`、
  すなわち `λ·yᵢ = nᵢ·2πi`。ここから `n₁·y₀ = n₀·y₁` という **ℚ 上の関係**が出て、
  `y` の一次独立性に反する(`latVec_injective`)。

★★**これで six exponentials の仮定(`x`・`y` の ℚ 上一次独立性)が
どこでどう効くかが全部見えた** ——
`y` は「格子点が相異なる」(`latticePt_injective`、零点の個数 `N³`)と
「指標が相異なる」の 2 か所、`x` は「指数が相異なる」で効く。
-/

namespace ABC3.Found.SixExp

open Metric Complex NumberField

/-! ## ★1. 倍率つき外挿 -/

variable {K : Type*} [Field K] [NumberField K]

/-- ★★★★**外挿の 1 段(倍率つき)** —— 格子点での値が `δ·F(w)` の形の
代数的整数になる場合(分母を払ったとき)。 -/
theorem extrapolation_step_scaled (σ : K →+* ℂ) {R r M H Δ : ℝ}
    (hr : 0 ≤ r) (h5 : 5 * r ≤ R) (hrR : r < R)
    {F : ℂ → ℂ} (hF : DifferentiableOn ℂ F (closedBall 0 R))
    {s : Finset ℂ} (hs : ∀ w ∈ s, ‖w‖ ≤ r) (hzero : ∀ w ∈ s, F w = 0)
    (hM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖F ζ‖ ≤ M)
    {w : ℂ} (hw : ‖w‖ ≤ r)
    {α : 𝓞 K} {δ : ℂ} (hδ0 : δ ≠ 0) (hδ : ‖δ‖ ≤ Δ) (hΔ0 : 0 ≤ Δ)
    (hα : σ ((α : 𝓞 K) : K) = δ * F w)
    (hH1 : 1 ≤ H) (hHle : house ((α : 𝓞 K) : K) ≤ H)
    (hsmall : Δ * (M * (1/2 : ℝ) ^ s.card) < (H ^ (Module.finrank ℚ K - 1))⁻¹) :
    F w = 0 := by
  have hb := norm_le_half_pow hr h5 hrR s F M hF hs hzero hM w hw
  have hlt : ‖σ ((α : 𝓞 K) : K)‖ < (H ^ (Module.finrank ℚ K - 1))⁻¹ := by
    rw [hα, norm_mul]
    calc ‖δ‖ * ‖F w‖ ≤ Δ * (M * (1/2 : ℝ) ^ s.card) :=
          mul_le_mul hδ hb (norm_nonneg _) hΔ0
      _ < _ := hsmall
  have hz : α = 0 := eq_zero_of_norm_embedding_lt σ hH1 hHle hlt
  have hzz : δ * F w = 0 := by
    rw [← hα, hz]
    simp
  exact (mul_eq_zero.mp hzz).resolve_left hδ0

/-- ★★★★★**外挿の帰納(倍率つき・半径が動く版)** ——
実際の six exponentials の証明で使う形。 -/
theorem extrapolation_induction_scaled (σ : K →+* ℂ) {F : ℂ → ℂ}
    (r R M H Δ : ℕ → ℝ)
    (hr : ∀ k, 0 ≤ r k) (h5 : ∀ k, 5 * r k ≤ R k) (hrR : ∀ k, r k < R k)
    (hF : ∀ k, DifferentiableOn ℂ F (closedBall 0 (R k)))
    (hMk : ∀ k, ∀ ζ : ℂ, ‖ζ‖ = R k → ‖F ζ‖ ≤ M k)
    (T : ℕ → Finset ℂ)
    (hTr : ∀ k, ∀ w ∈ T (k + 1), ‖w‖ ≤ r k)
    (hmono : ∀ k, T k ⊆ T (k + 1))
    (h0 : ∀ w ∈ T 0, F w = 0)
    (hH1 : ∀ k, 1 ≤ H k) (hΔ0 : ∀ k, 0 ≤ Δ k)
    (halg : ∀ k, ∀ w ∈ T (k + 1), ∃ (α : 𝓞 K) (δ : ℂ), δ ≠ 0 ∧ ‖δ‖ ≤ Δ k ∧
      σ ((α : 𝓞 K) : K) = δ * F w ∧ house ((α : 𝓞 K) : K) ≤ H k)
    (hgap : ∀ k, Δ k * (M k * (1/2 : ℝ) ^ (T k).card)
      < (H k ^ (Module.finrank ℚ K - 1))⁻¹) :
    ∀ k, ∀ w ∈ T k, F w = 0 := by
  intro k
  induction k with
  | zero => exact h0
  | succ n ih =>
    intro w hw
    obtain ⟨α, δ, hδ0, hδ, hα, hHle⟩ := halg n w hw
    have hsr : ∀ v ∈ T n, ‖v‖ ≤ r n := fun v hv => hTr n v (hmono n hv)
    exact extrapolation_step_scaled σ (hr n) (h5 n) (hrR n) (hF n) hsr ih (hMk n)
      (hTr n w hw) hδ0 hδ (hΔ0 n) hα (hH1 n) hHle (hgap n)

/-! ## ★2. 指標の相異性 -/

theorem latVec_sixExpVals (x : Fin 2 → ℂ) (y : Fin 3 → ℂ) (pq : ℕ × ℕ) (i : Fin 3) :
    latVec (sixExpVals x y) pq i = Complex.exp (auxExp x pq * y i) := by
  rw [latVec, sixExpVals, sixExpVals, auxExp, ← Complex.exp_nat_mul, ← Complex.exp_nat_mul,
    ← Complex.exp_add]
  ring_nf

/-- ★`x` が ℚ 上一次独立なら、相異なる `(p,q)` は相異なる指数を与える。 -/
theorem auxExp_injective {x : Fin 2 → ℂ} (hx : LinearIndependent ℚ x) :
    Function.Injective (auxExp x) := by
  rintro ⟨p, q⟩ ⟨p', q'⟩ h
  set a : Fin 2 → ℚ := ![(p : ℚ) - (p' : ℚ), (q : ℚ) - (q' : ℚ)] with ha
  have ha0 : a 0 = (p : ℚ) - (p' : ℚ) := rfl
  have ha1 : a 1 = (q : ℚ) - (q' : ℚ) := rfl
  have hsum : ∑ j : Fin 2, a j • x j = 0 := by
    rw [Fin.sum_univ_two, ha0, ha1]
    rw [auxExp, auxExp] at h
    simp only at h
    have h2 : ((p : ℂ) * x 0 + (q : ℂ) * x 1) - ((p' : ℂ) * x 0 + (q' : ℂ) * x 1) = 0 := by
      rw [h, sub_self]
    rw [Rat.smul_def, Rat.smul_def]
    push_cast
    linear_combination h2
  have hz := Fintype.linearIndependent_iff.mp hx a hsum
  have h0 := hz 0
  have h1 := hz 1
  rw [ha0] at h0
  rw [ha1] at h1
  rw [sub_eq_zero] at h0 h1
  have hp : p = p' := by exact_mod_cast h0
  have hq : q = q' := by exact_mod_cast h1
  rw [hp, hq]

/-- ★★★**`x`・`y` がともに ℚ 上一次独立なら、`latVec` は単射**。

★これが `coeff_eq_zero_of_auxFun_latticePt`(背理法の閉じ目)の仮定である。 -/
theorem latVec_injective {x : Fin 2 → ℂ} {y : Fin 3 → ℂ}
    (hx : LinearIndependent ℚ x) (hy : LinearIndependent ℚ y) :
    Function.Injective (latVec (sixExpVals x y)) := by
  intro pq pq' h
  by_contra hne
  set l := auxExp x pq - auxExp x pq' with hl
  have hl0 : l ≠ 0 := sub_ne_zero.mpr (fun heq => hne (auxExp_injective hx heq))
  have hexp : ∀ i, Complex.exp (l * y i) = 1 := by
    intro i
    have h1 : Complex.exp (auxExp x pq * y i) = Complex.exp (auxExp x pq' * y i) := by
      have hi := congrFun h i
      rwa [latVec_sixExpVals, latVec_sixExpVals] at hi
    have h2 : l * y i = auxExp x pq * y i - auxExp x pq' * y i := by rw [hl]; ring
    rw [h2, Complex.exp_sub, h1, div_self (Complex.exp_ne_zero _)]
  choose n hn using fun i => Complex.exp_eq_one_iff.mp (hexp i)
  have hy0 : y 0 ≠ 0 := hy.ne_zero 0
  have hn0 : (n 0 : ℂ) ≠ 0 := by
    intro hz
    have h0 := hn 0
    rw [hz, zero_mul] at h0
    exact hy0 ((mul_eq_zero.mp h0).resolve_left hl0)
  have hrel : -((n 1 : ℂ)) * y 0 + (n 0 : ℂ) * y 1 = 0 := by
    have e0 := hn 0
    have e1 := hn 1
    have hz : l * (-((n 1 : ℂ)) * y 0 + (n 0 : ℂ) * y 1) = 0 := by
      have h1 : l * (-((n 1 : ℂ)) * y 0) = -((n 1 : ℂ)) * (l * y 0) := by ring
      have h2 : l * ((n 0 : ℂ) * y 1) = ((n 0 : ℂ)) * (l * y 1) := by ring
      rw [mul_add, h1, h2, e0, e1]
      ring
    exact (mul_eq_zero.mp hz).resolve_left hl0
  set a : Fin 3 → ℚ := ![-((n 1 : ℚ)), ((n 0 : ℚ)), 0] with ha
  have ha0 : a 0 = -((n 1 : ℚ)) := rfl
  have ha1 : a 1 = ((n 0 : ℚ)) := rfl
  have ha2 : a 2 = 0 := rfl
  have hsum : ∑ i : Fin 3, a i • y i = 0 := by
    rw [Fin.sum_univ_three, ha0, ha1, ha2]
    rw [Rat.smul_def, Rat.smul_def, Rat.smul_def]
    push_cast
    linear_combination hrel
  have hz := Fintype.linearIndependent_iff.mp hy a hsum
  have h1 := hz 1
  rw [ha1] at h1
  exact hn0 (by exact_mod_cast h1)

end ABC3.Found.SixExp
