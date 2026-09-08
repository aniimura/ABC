import ABC3.Found.PGC.TotallyRamifiedLayer

/-!
# 骨組み(作業中) —— 跳びの下界 `1 ≤ u`
-/

namespace ABC3.Found.PGC

namespace WildBreak

open Finset

/-! ## §1 抽象核(超距離だけ) -/

section Ultra

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- 各項が `r` 以下なら有限和も `r` 以下(超距離)。 -/
theorem norm_sum_le_of_forall_le {ι : Type*} (s : Finset ι) (x : ι → M) {r : ℝ}
    (hr : 0 ≤ r) (h : ∀ i ∈ s, ‖x i‖ ≤ r) : ‖∑ i ∈ s, x i‖ ≤ r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hr
  | insert a s ha ih =>
      rw [Finset.sum_insert ha]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans ?_
      exact max_le (h a (Finset.mem_insert_self a s))
        (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

/-- 各項が `r` 未満なら有限和も `r` 未満(超距離、`0 < r`)。 -/
theorem norm_sum_lt_of_forall_lt {ι : Type*} (s : Finset ι) (x : ι → M) {r : ℝ}
    (hr : 0 < r) (h : ∀ i ∈ s, ‖x i‖ < r) : ‖∑ i ∈ s, x i‖ < r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hr
  | insert a s ha ih =>
      rw [Finset.sum_insert ha]
      refine lt_of_le_of_lt (IsUltrametricDist.norm_add_le_max _ _) ?_
      exact max_lt (h a (Finset.mem_insert_self a s))
        (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

/-- ★★**抽象核** —— ノルムが**相異なる**族では、各項のノルムは**和のノルム以下**。

★超距離で「打ち消しが起きない」ことの中身。分岐・付値の語彙が 1 語も出ない。 -/
theorem norm_le_norm_sum_of_pairwise_ne {ι : Type*} [Fintype ι] [DecidableEq ι] (x : ι → M)
    (hpair : ∀ i j, i ≠ j → ‖x i‖ ≠ ‖x j‖) (i₀ : ι) : ‖x i₀‖ ≤ ‖∑ i, x i‖ := by
  obtain ⟨i₁, -, hmax⟩ :=
    Finset.exists_max_image (Finset.univ : Finset ι) (fun i => ‖x i‖) ⟨i₀, Finset.mem_univ i₀⟩
  rcases (norm_nonneg (x i₁)).lt_or_eq with hpos | hzero
  · have hlt : ‖∑ i ∈ Finset.univ.erase i₁, x i‖ < ‖x i₁‖ := by
      refine norm_sum_lt_of_forall_lt _ _ hpos (fun i hi => ?_)
      have hne : i ≠ i₁ := Finset.ne_of_mem_erase hi
      exact lt_of_le_of_ne (hmax i (Finset.mem_univ i)) (hpair i i₁ hne)
    have hsum : ‖∑ i, x i‖ = ‖x i₁‖ := by
      rw [← Finset.add_sum_erase _ x (Finset.mem_univ i₁),
        IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_gt hlt)]
      exact max_eq_left (le_of_lt hlt)
    rw [hsum]
    exact hmax i₀ (Finset.mem_univ i₀)
  · have : ‖x i₀‖ = 0 :=
      le_antisymm (by rw [hzero]; exact hmax i₀ (Finset.mem_univ i₀)) (norm_nonneg _)
    rw [this]
    exact norm_nonneg _

end Ultra

end WildBreak

end ABC3.Found.PGC
