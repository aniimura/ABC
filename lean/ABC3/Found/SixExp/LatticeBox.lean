/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.LatticeGeneral

/-!
# 格子の箱と解析側の上界

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を作るか

外挿の帰納 `extrapolation_induction'` に渡す **零点集合の鎖 `T k`** と、
**解析側の上界 `M k`** を用意する。

* `latticePt_injective` —— ★`y` が **ℚ 上一次独立**なら格子点は相異なる。
  これが `|T k|` を数える根拠であり、six exponentials の仮定が効く一点目である。
* `natBox N` = `[0,N)³`、`latBox y N` = その像。★`|latBox y N| = N³`。
* `norm_le_of_mem_latBox` —— 箱の点は `|z| ≤ 3·N·Y` の中(`Y = max‖yᵢ‖`)。
* `norm_auxFun_le_of_box` —— `|z| = R` で `|F| ≤ (∑‖c‖)·exp(2·L·X·R)`
  (`X = max‖xⱼ‖`、`p,q ≤ L`)。★これが `M k` である。

★★**数え上げの構図がここで見える** ——
零点の個数は `N³`、`M k` の指数は `2·L·X·R ~ L·N`、house は `B^{3N}`。
`Counting.lean` の「3 乗は 1 乗に勝つ」がこの 3 つの比較である。
-/

namespace ABC3.Found.SixExp

open Complex Finset

/-! ## ★1. 格子点は相異なる -/

/-- ★★`y` が ℚ 上一次独立なら、格子点は相異なる。

★これが `|T k| = N³` を数える根拠であり、`y` が 3 個ある効き目でもある。 -/
theorem latticePt_injective {y : Fin 3 → ℂ} (hy : LinearIndependent ℚ y) :
    Function.Injective (latticePt y) := by
  intro m m' h
  set a : Fin 3 → ℚ := fun i => (m i : ℚ) - (m' i : ℚ) with ha
  have hsum : ∑ i : Fin 3, a i • y i = 0 := by
    have h1 : ∑ i : Fin 3, ((a i : ℚ) : ℂ) * y i = 0 := by
      have h2 : (∑ i : Fin 3, (m i : ℂ) * y i) - (∑ i : Fin 3, (m' i : ℂ) * y i) = 0 := by
        rw [latticePt, latticePt] at h
        rw [h, sub_self]
      rw [← Finset.sum_sub_distrib] at h2
      rw [← h2]
      refine Finset.sum_congr rfl (fun i _ => ?_)
      rw [ha]
      push_cast
      ring
    rw [← h1]
    exact Finset.sum_congr rfl (fun i _ => by rw [Rat.smul_def])
  have hz := Fintype.linearIndependent_iff.mp hy a hsum
  funext i
  have hi := hz i
  rw [ha] at hi
  simp only [sub_eq_zero] at hi
  exact_mod_cast hi

/-! ## ★2. 格子の箱 -/

/-- ★格子の箱 `[0,N)³` の添字。 -/
def natBox (N : ℕ) : Finset (Fin 3 → ℕ) := Fintype.piFinset (fun _ => Finset.range N)

theorem mem_natBox {N : ℕ} {m : Fin 3 → ℕ} : m ∈ natBox N ↔ ∀ i, m i < N := by
  simp [natBox, Fintype.mem_piFinset]

theorem card_natBox (N : ℕ) : (natBox N).card = N ^ 3 := by
  rw [natBox, Fintype.card_piFinset]
  simp [pow_succ]

theorem zero_mem_natBox {N : ℕ} (hN : 0 < N) : (fun _ => 0) ∈ natBox N :=
  mem_natBox.mpr (fun _ => hN)

theorem natBox_mono {N M : ℕ} (h : N ≤ M) : natBox N ⊆ natBox M := by
  intro m hm
  exact mem_natBox.mpr (fun i => lt_of_lt_of_le (mem_natBox.mp hm i) h)

theorem sum_le_of_mem_natBox {N : ℕ} {m : Fin 3 → ℕ} (h : m ∈ natBox N) :
    ∑ i : Fin 3, m i ≤ 3 * N := by
  have hlt : ∀ i, m i ≤ N := fun i => le_of_lt (mem_natBox.mp h i)
  calc ∑ i : Fin 3, m i ≤ ∑ _i : Fin 3, N := Finset.sum_le_sum (fun i _ => hlt i)
    _ = 3 * N := by simp [mul_comm]

/-- ★複素平面の中の格子の箱。 -/
noncomputable def latBox (y : Fin 3 → ℂ) (N : ℕ) : Finset ℂ :=
  (natBox N).image (latticePt y)

/-- ★★**`|latBox y N| = N³`** —— 零点の個数。 -/
theorem card_latBox {y : Fin 3 → ℂ} (hy : LinearIndependent ℚ y) (N : ℕ) :
    (latBox y N).card = N ^ 3 := by
  rw [latBox, Finset.card_image_of_injective _ (latticePt_injective hy), card_natBox]

theorem latBox_mono (y : Fin 3 → ℂ) {N M : ℕ} (h : N ≤ M) : latBox y N ⊆ latBox y M :=
  Finset.image_subset_image (natBox_mono h)

theorem norm_latticePt_le {y : Fin 3 → ℂ} {Y : ℝ} (hY : ∀ i, ‖y i‖ ≤ Y)
    {N : ℕ} {m : Fin 3 → ℕ} (hm : m ∈ natBox N) :
    ‖latticePt y m‖ ≤ 3 * N * Y := by
  have h1 : ‖latticePt y m‖ ≤ ∑ i : Fin 3, ‖(m i : ℂ) * y i‖ := by
    rw [latticePt]
    exact norm_sum_le _ _
  refine le_trans h1 ?_
  have h2 : ∀ i : Fin 3, ‖(m i : ℂ) * y i‖ ≤ (N : ℝ) * Y := by
    intro i
    rw [norm_mul]
    have hmi : ‖((m i : ℂ))‖ ≤ (N : ℝ) := by
      rw [Complex.norm_natCast]
      exact_mod_cast le_of_lt (mem_natBox.mp hm i)
    exact mul_le_mul hmi (hY i) (norm_nonneg _) (by positivity)
  calc ∑ i : Fin 3, ‖(m i : ℂ) * y i‖ ≤ ∑ _i : Fin 3, (N : ℝ) * Y :=
        Finset.sum_le_sum (fun i _ => h2 i)
    _ = 3 * ((N : ℝ) * Y) := by simp [mul_comm]
    _ = 3 * N * Y := by ring

/-- ★箱の点は `|z| ≤ 3·N·Y` の中。 -/
theorem norm_le_of_mem_latBox {y : Fin 3 → ℂ} {Y : ℝ} (hY : ∀ i, ‖y i‖ ≤ Y)
    {N : ℕ} {w : ℂ} (hw : w ∈ latBox y N) : ‖w‖ ≤ 3 * N * Y := by
  rw [latBox, Finset.mem_image] at hw
  obtain ⟨m, hm, rfl⟩ := hw
  exact norm_latticePt_le hY hm

/-! ## ★3. 解析側の上界 -/

theorem norm_auxExp_le {x : Fin 2 → ℂ} {X : ℝ} (hX : ∀ j, ‖x j‖ ≤ X)
    {L : ℕ} {pq : ℕ × ℕ} (hp : pq.1 ≤ L) (hq : pq.2 ≤ L) :
    ‖auxExp x pq‖ ≤ 2 * L * X := by
  have h1 : ‖auxExp x pq‖ ≤ ‖((pq.1 : ℂ)) * x 0‖ + ‖((pq.2 : ℂ)) * x 1‖ := by
    rw [auxExp]
    exact norm_add_le _ _
  refine le_trans h1 ?_
  have h2 : ‖((pq.1 : ℂ)) * x 0‖ ≤ (L : ℝ) * X := by
    rw [norm_mul, Complex.norm_natCast]
    exact mul_le_mul (by exact_mod_cast hp) (hX 0) (norm_nonneg _) (by positivity)
  have h3 : ‖((pq.2 : ℂ)) * x 1‖ ≤ (L : ℝ) * X := by
    rw [norm_mul, Complex.norm_natCast]
    exact mul_le_mul (by exact_mod_cast hq) (hX 1) (norm_nonneg _) (by positivity)
  calc ‖((pq.1 : ℂ)) * x 0‖ + ‖((pq.2 : ℂ)) * x 1‖ ≤ (L : ℝ) * X + (L : ℝ) * X :=
        add_le_add h2 h3
    _ = 2 * L * X := by ring

/-- ★★**解析側の上界 `M`** —— `|z| = R` での `|F|` の評価。 -/
theorem norm_auxFun_le_of_box {x : Fin 2 → ℂ} {X : ℝ} (hX : ∀ j, ‖x j‖ ≤ X)
    {L : ℕ} (s : Finset (ℕ × ℕ)) (hsL : ∀ pq ∈ s, pq.1 ≤ L ∧ pq.2 ≤ L)
    (c : ℕ × ℕ → ℂ) {R : ℝ} {z : ℂ} (hz : ‖z‖ = R) :
    ‖auxFun x s c z‖ ≤ (∑ pq ∈ s, ‖c pq‖) * Real.exp (2 * L * X * R) := by
  have hB : ∀ pq ∈ s, ‖auxExp x pq‖ ≤ 2 * L * X := fun pq hpq =>
    norm_auxExp_le hX (hsL pq hpq).1 (hsL pq hpq).2
  have h := norm_auxFun_le x s c hB z
  rw [hz] at h
  exact h

end ABC3.Found.SixExp
