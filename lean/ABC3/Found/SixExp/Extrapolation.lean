/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.SchwarzZeros
import ABC3.Found.SixExp.Liouville
import ABC3.Found.SixExp.AuxFunction
import Mathlib.LinearAlgebra.LinearIndependent.Defs
import Mathlib.FieldTheory.Finite.Basic

/-!
# 外挿と締め —— チェーン `sixexp` の節点 `extrapolation`

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★この段で何をするか

超越性証明の骨は次の 3 つで、前 2 つは既に在庫にある:

| 段 | 場所 |
|---|---|
| 解析側の小ささ(多零点版 Schwarz) | `SchwarzZeros.lean` |
| 算術側の下界(Liouville 不等式) | `Liouville.lean` |
| **その 2 つを突き合わせる**(外挿) | ★本ファイル |

★**外挿の 1 段**(`extrapolation_step`)——
`F` が `s` の各点で消え `|z|=R` で `|F| ≤ M`、`5r ≤ R` なら
`|z| ≤ r` で `|F| ≤ M·(1/2)^{|s|}`。新しい点 `w` での値が house `≤ H` の
代数的整数の像なら、`M·(1/2)^{|s|} < H^{-(d-1)}` の瞬間に **`F w = 0`**。

★**外挿の帰納**(`extrapolation_induction`)——
零点集合 `T k` が段ごとに増え、各段でこの不等式が成り立てば、
**すべての段で `F` は消える**。

## ★★締め —— 「すべての格子点で消える ⟹ 係数が 0」

外挿が終わると `F` は**すべての格子点で消える**。そこから係数が 0 になることは、
★**乗法的指標の一次独立性**(Dedekind–Artin。mathlib の
`linearIndependent_monoidHom`)で出る:

* `coeff_eq_zero_of_auxFun_eq_zero` —— `z ↦ exp(l z)` を `ℂ`(加法)上の指標と見る版
* `coeff_eq_zero_of_latticeVals_eq_zero` —— `m ↦ ∏ᵢ (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}` を
  **`ℕ³` 上の指標**と見る版。★格子点だけで消える情報から係数を潰せるのはこちら。

★★**これが背理法の閉じ目である** —— Siegel の補題は「0 でない係数」を与えるので、
「すべての格子点で消える」と矛盾する。

## ★残っているもの(記録)

★本ファイルは**外挿の骨組みと締め**であって、six exponentials theorem そのものでは
まだない。残るのは **パラメータの数え上げ**(`M`・`H k`・`|T k|` を具体的に選んで
`hgap` を満たすこと)と、格子点での値が実際に house `≤ H k` の代数的整数であること
(`AuxFunction.lean` の `isIntegral_auxMatrix` と Siegel の補題の係数評価をつなぐ)。
-/

namespace ABC3.Found.SixExp

open Metric Complex NumberField

/-! ## ★1. 外挿 -/

variable {K : Type*} [Field K] [NumberField K]

/-- ★★★★**外挿の 1 段** —— 解析側の小ささと算術側の下界を突き合わせる。

`F` が半径 `R` の閉円板で正則、`s ⊆ {|z| ≤ r}` の各点で消え、`|z| = R` で `|F| ≤ M`、
`5r ≤ R` とする。新しい点 `w`(`|w| ≤ r`)での値 `F w` が代数的整数 `α` の像で
`house α ≤ H` なら、**解析側の上界 `M·(1/2)^{|s|}` が Liouville の下界 `H^{-(d-1)}` を
下回った瞬間に `F w = 0`** である。 -/
theorem extrapolation_step (σ : K →+* ℂ) {R r M H : ℝ}
    (hr : 0 ≤ r) (h5 : 5 * r ≤ R) (hrR : r < R)
    {F : ℂ → ℂ} (hF : DifferentiableOn ℂ F (closedBall 0 R))
    {s : Finset ℂ} (hs : ∀ w ∈ s, ‖w‖ ≤ r) (hzero : ∀ w ∈ s, F w = 0)
    (hM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖F ζ‖ ≤ M)
    {w : ℂ} (hw : ‖w‖ ≤ r)
    {α : 𝓞 K} (hα : σ (α : K) = F w)
    (hH1 : 1 ≤ H) (hHle : house (α : K) ≤ H)
    (hsmall : M * (1/2 : ℝ) ^ s.card < (H ^ (Module.finrank ℚ K - 1))⁻¹) :
    F w = 0 := by
  have hb := norm_le_half_pow hr h5 hrR s F M hF hs hzero hM w hw
  have hlt : ‖σ (α : K)‖ < (H ^ (Module.finrank ℚ K - 1))⁻¹ := by
    rw [hα]
    exact lt_of_le_of_lt hb hsmall
  have hz : α = 0 := eq_zero_of_norm_embedding_lt σ hH1 hHle hlt
  rw [← hα, hz]
  simp

/-- ★★★★★**外挿の帰納** —— `sixexp` チェーンの節点 `extrapolation`。

各段の格子点の集合 `T k` について、
* `T 0` では `F` が消える(Siegel の補題で係数を選んだ結果)
* `T (k+1)` の各点での値は、house が `H k` で抑えられる代数的整数の像
* 解析側の上界 `M·(1/2)^{|T k|}` が Liouville の下界 `H k^{-(d-1)}` を下回る

なら、**すべての段で `F` は消える**。 -/
theorem extrapolation_induction (σ : K →+* ℂ) {R r M : ℝ}
    (hr : 0 ≤ r) (h5 : 5 * r ≤ R) (hrR : r < R)
    {F : ℂ → ℂ} (hF : DifferentiableOn ℂ F (closedBall 0 R))
    (hM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖F ζ‖ ≤ M)
    (T : ℕ → Finset ℂ) (hT : ∀ k, ∀ w ∈ T k, ‖w‖ ≤ r)
    (h0 : ∀ w ∈ T 0, F w = 0)
    (H : ℕ → ℝ) (hH1 : ∀ k, 1 ≤ H k)
    (halg : ∀ k, ∀ w ∈ T (k + 1), ∃ α : 𝓞 K, σ (α : K) = F w ∧ house (α : K) ≤ H k)
    (hgap : ∀ k, M * (1/2 : ℝ) ^ (T k).card < (H k ^ (Module.finrank ℚ K - 1))⁻¹) :
    ∀ k, ∀ w ∈ T k, F w = 0 := by
  intro k
  induction k with
  | zero => exact h0
  | succ n ih =>
    intro w hw
    obtain ⟨α, hα, hHle⟩ := halg n w hw
    exact extrapolation_step σ hr h5 hrR hF (hT n) ih hM (hT (n + 1) w hw) hα (hH1 n) hHle
      (hgap n)

/-! ## ★2. 締め(その 1) —— `ℂ` 上の指標 -/

/-- ★指数関数 `z ↦ exp(l·z)` を、加法群 `ℂ` 上の乗法的指標として見る。 -/
noncomputable def expChar (l : ℂ) : Multiplicative ℂ →* ℂ where
  toFun z := Complex.exp (l * Multiplicative.toAdd z)
  map_one' := by simp
  map_mul' a b := by
    show Complex.exp (l * (Multiplicative.toAdd a + Multiplicative.toAdd b))
      = Complex.exp (l * Multiplicative.toAdd a) * Complex.exp (l * Multiplicative.toAdd b)
    rw [mul_add, Complex.exp_add]

@[simp] theorem expChar_apply (l z : ℂ) :
    expChar l (Multiplicative.ofAdd z) = Complex.exp (l * z) := rfl

/-- ★★**相異なる `l` は相異なる指標を与える**。

★`exp(l z) = exp(m z)` を `z = (l−m)⁻¹` で見ると `exp 1 = 1` になり、
`exp_eq_one_iff` の実部で矛盾する。 -/
theorem expChar_injective : Function.Injective expChar := by
  intro l m h
  by_contra hne
  have hd : l - m ≠ 0 := sub_ne_zero.mpr hne
  have h1 : Complex.exp (l * (l - m)⁻¹) = Complex.exp (m * (l - m)⁻¹) := by
    have h0 := congrArg (fun f : Multiplicative ℂ →* ℂ => f (Multiplicative.ofAdd ((l - m)⁻¹))) h
    simpa using h0
  have h2 : l * (l - m)⁻¹ = m * (l - m)⁻¹ + 1 := by
    field_simp
    ring
  rw [h2, Complex.exp_add] at h1
  have h4 : Complex.exp 1 = 1 := by
    field_simp at h1
    exact h1
  obtain ⟨n, hn⟩ := Complex.exp_eq_one_iff.mp h4
  have h5 : (1 : ℂ).re = ((n : ℂ) * (2 * Real.pi * Complex.I)).re := by rw [← hn]
  simp at h5

/-- ★★★**指数和が恒等的に 0 なら係数はすべて 0**(指数が相異なるとき)。 -/
theorem coeff_eq_zero_of_auxFun_eq_zero (x : Fin 2 → ℂ) (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℂ)
    (hinj : Set.InjOn (auxExp x) s)
    (h : ∀ z : ℂ, auxFun x s c z = 0) : ∀ pq ∈ s, c pq = 0 := by
  classical
  set w : ↥s → (Multiplicative ℂ →* ℂ) := fun i => expChar (auxExp x (i : ℕ × ℕ)) with hw
  have hwinj : Function.Injective w := by
    intro a b hab
    exact Subtype.ext (hinj a.2 b.2 (expChar_injective hab))
  have hli : LinearIndependent ℂ ((fun f : Multiplicative ℂ →* ℂ => ⇑f) ∘ w) :=
    (linearIndependent_monoidHom (Multiplicative ℂ) ℂ).comp w hwinj
  have hsum : ∑ i : ↥s, c (i : ℕ × ℕ) • ((fun f : Multiplicative ℂ →* ℂ => ⇑f) ∘ w) i = 0 := by
    funext z
    have hz := h (Multiplicative.toAdd z)
    simp only [Finset.sum_apply, Pi.smul_apply, Function.comp_apply, smul_eq_mul, Pi.zero_apply]
    rw [← hz, auxFun, ← Finset.sum_attach s
      (fun pq => c pq * Complex.exp (auxExp x pq * Multiplicative.toAdd z))]
    rfl
  intro pq hpq
  exact Fintype.linearIndependent_iff.mp hli (fun i => c (i : ℕ × ℕ)) hsum ⟨pq, hpq⟩

/-! ## ★3. 締め(その 2) —— `ℕ³` 上の指標

★**格子点だけで消える情報から係数を潰せるのはこちら**である。
`m ↦ ∏ᵢ (E₀ᵢ^p E₁ᵢ^q)^{mᵢ}` は `ℕ³` 上の乗法的指標だから。 -/

/-- ★格子点の指標が決まるベクトル `i ↦ E₀ᵢ^p · E₁ᵢ^q`。 -/
noncomputable def latVec (E : Fin 2 → Fin 3 → ℂ) (pq : ℕ × ℕ) (i : Fin 3) : ℂ :=
  E 0 i ^ pq.1 * E 1 i ^ pq.2

/-- ★★格子点の値が定める `ℕ³` 上の乗法的指標。 -/
noncomputable def latChar (E : Fin 2 → Fin 3 → ℂ) (pq : ℕ × ℕ) :
    Multiplicative (Fin 3 → ℕ) →* ℂ where
  toFun m := ∏ i : Fin 3, latVec E pq i ^ (Multiplicative.toAdd m i)
  map_one' := by simp
  map_mul' a b := by
    show (∏ i : Fin 3, latVec E pq i ^ ((Multiplicative.toAdd a + Multiplicative.toAdd b) i))
      = (∏ i : Fin 3, latVec E pq i ^ (Multiplicative.toAdd a i))
        * (∏ i : Fin 3, latVec E pq i ^ (Multiplicative.toAdd b i))
    rw [← Finset.prod_mul_distrib]
    exact Finset.prod_congr rfl (fun i _ => by rw [Pi.add_apply, pow_add])

@[simp] theorem latChar_apply (E : Fin 2 → Fin 3 → ℂ) (pq : ℕ × ℕ) (m : Fin 3 → ℕ) :
    latChar E pq (Multiplicative.ofAdd m) = ∏ i : Fin 3, latVec E pq i ^ (m i) := rfl

theorem auxMatrix_eq_latChar (E : Fin 2 → Fin 3 → ℂ) (m : Fin 3 → ℕ) (pq : ℕ × ℕ) :
    auxMatrix E m pq = latChar E pq (Multiplicative.ofAdd m) := rfl

/-- ★指標が等しければベクトルも等しい —— 標準基底で見る。 -/
theorem latVec_of_latChar_eq {E : Fin 2 → Fin 3 → ℂ} {pq pq' : ℕ × ℕ}
    (h : latChar E pq = latChar E pq') : latVec E pq = latVec E pq' := by
  classical
  funext i
  have h1 := congrArg (fun f : Multiplicative (Fin 3 → ℕ) →* ℂ =>
    f (Multiplicative.ofAdd (Pi.single i 1))) h
  simp only [latChar_apply] at h1
  have h2 : ∀ (v : Fin 3 → ℂ), (∏ j : Fin 3, v j ^ ((Pi.single i 1 : Fin 3 → ℕ) j)) = v i := by
    intro v
    rw [Finset.prod_eq_single i]
    · simp
    · intro j _ hj
      simp [hj]
    · intro hi
      exact absurd (Finset.mem_univ i) hi
  rw [h2, h2] at h1
  exact h1

/-- ★★★★**格子点でつねに消えるなら係数はすべて 0** ——
外挿の帰納が終わったあとの「締め」の段。 -/
theorem coeff_eq_zero_of_latticeVals_eq_zero (E : Fin 2 → Fin 3 → ℂ) (s : Finset (ℕ × ℕ))
    (c : ℕ × ℕ → ℂ) (hinj : Set.InjOn (latVec E) s)
    (h : ∀ m : Fin 3 → ℕ, ∑ pq ∈ s, auxMatrix E m pq * c pq = 0) :
    ∀ pq ∈ s, c pq = 0 := by
  classical
  set w : ↥s → (Multiplicative (Fin 3 → ℕ) →* ℂ) := fun i => latChar E (i : ℕ × ℕ) with hw
  have hwinj : Function.Injective w := by
    intro a b hab
    exact Subtype.ext (hinj a.2 b.2 (latVec_of_latChar_eq hab))
  have hli : LinearIndependent ℂ ((fun f : Multiplicative (Fin 3 → ℕ) →* ℂ => ⇑f) ∘ w) :=
    (linearIndependent_monoidHom (Multiplicative (Fin 3 → ℕ)) ℂ).comp w hwinj
  have hsum : ∑ i : ↥s, c (i : ℕ × ℕ) • ((fun f : Multiplicative (Fin 3 → ℕ) →* ℂ => ⇑f) ∘ w) i
      = 0 := by
    funext m
    have hm := h (Multiplicative.toAdd m)
    simp only [Finset.sum_apply, Pi.smul_apply, Function.comp_apply, smul_eq_mul, Pi.zero_apply]
    rw [← hm, ← Finset.sum_attach s
      (fun pq => auxMatrix E (Multiplicative.toAdd m) pq * c pq)]
    exact Finset.sum_congr rfl (fun i _ => by
      rw [auxMatrix_eq_latChar]
      exact mul_comm _ _)
  intro pq hpq
  exact Fintype.linearIndependent_iff.mp hli (fun i => c (i : ℕ × ℕ)) hsum ⟨pq, hpq⟩

/-- ★★★★★**外挿が終わったあとの締め** —— 補助関数が**すべての格子点で消える**なら、
係数はすべて 0 である。

★これは Siegel の補題が「0 でない係数」を与えることと矛盾する ——
すなわち six exponentials theorem の背理法の閉じ目である。 -/
theorem coeff_eq_zero_of_auxFun_latticePt (x : Fin 2 → ℂ) (y : Fin 3 → ℂ)
    (s : Finset (ℕ × ℕ)) (c : ℕ × ℕ → ℂ)
    (hinj : Set.InjOn (latVec (sixExpVals x y)) s)
    (h : ∀ m : Fin 3 → ℕ, auxFun x s c (latticePt y m) = 0) :
    ∀ pq ∈ s, c pq = 0 := by
  refine coeff_eq_zero_of_latticeVals_eq_zero (sixExpVals x y) s c hinj (fun m => ?_)
  rw [← auxFun_latticePt]
  exact h m

end ABC3.Found.SixExp
