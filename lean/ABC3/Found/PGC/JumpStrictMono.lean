import ABC3.Found.PGC.WildBreakUpperBound

/-!
# 骨組み(作業中) —— 跳びの狭義単調 `u m < u (m+1)`
-/

namespace ABC3.Found.PGC

namespace JumpMono

open Finset WildBreak

/-! ## §1 抽象核(超距離) —— 積の 1 次の項を取り出す -/

section Ultra

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★**抽象核** —— `‖η i‖ ≤ r ≤ 1` なら
`‖∏(1 + η i) − 1 − Σ η i‖ ≤ r²`(★1 次の項を取り出した残りは 2 次)。 -/
theorem norm_prod_one_add_sub_one_sub_sum_le {ι : Type*} (s : Finset ι) (η : ι → M) {r : ℝ}
    (hr0 : 0 ≤ r) (hr1 : r ≤ 1) (h : ∀ i ∈ s, ‖η i‖ ≤ r) :
    ‖(∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i‖ ≤ r * r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using mul_nonneg hr0 hr0
  | insert a s ha ih =>
      have hha : ‖η a‖ ≤ r := h a (Finset.mem_insert_self a s)
      have hsub : ∀ i ∈ s, ‖η i‖ ≤ r := fun i hi => h i (Finset.mem_insert_of_mem hi)
      have hP : ‖(∏ i ∈ s, (1 + η i)) - 1‖ ≤ r :=
        norm_prod_one_add_sub_one_le s η hr0 hr1 hsub
      have hIH : ‖(∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i‖ ≤ r * r := ih hsub
      rw [Finset.prod_insert ha, Finset.sum_insert ha,
        show (1 + η a) * (∏ i ∈ s, (1 + η i)) - 1 - (η a + ∑ i ∈ s, η i)
          = ((∏ i ∈ s, (1 + η i)) - 1 - ∑ i ∈ s, η i)
            + η a * ((∏ i ∈ s, (1 + η i)) - 1) by ring]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le hIH ?_)
      rw [norm_mul]
      exact mul_le_mul hha hP (norm_nonneg _) hr0

end Ultra

/-! ## §2 「`h` は元を `‖π‖^t` だけしか動かさない」 -/

section Expansion

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- 素元冪の基底展開。★各項のノルムは全体のノルム以下(§1 の核)。

★`WildBreak.exists_sub_algebraMap_norm_le` の前半と同じ構成だが、
そちらは `c 0` だけを取り出す形なので、ここでは**係数の族ごと**返す形に切り直した。 -/
theorem exists_coeff_norm_le [FiniteDimensional K M] {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (z : M) :
    ∃ c : Fin n → K, (∑ l, c l • π ^ (l : ℕ)) = z ∧ ∀ l, ‖c l • π ^ (l : ℕ)‖ ≤ ‖z‖ := by
  classical
  have hπne : ‖π‖ ≠ 1 := ne_of_lt hπ1
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine TotallyRamified.linearIndependent_of_ne_mod hπ0 hπne hvalK _
      (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  have hzmem : z ∈ Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) := by
    rw [hli.span_eq_top_of_card_eq_finrank' (by simp [hn])]; exact Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp hzmem
  refine ⟨c, hc, ?_⟩
  set x : Fin n → M := fun l => c l • π ^ (l : ℕ) with hxdef
  have hxval : ∀ l : Fin n, x l = c l • π ^ (l : ℕ) := fun _ => rfl
  have hsum : ∑ l, x l = z := hc
  have hexp : ∀ l : Fin n, c l ≠ 0 →
      ∃ m : ℤ, ‖x l‖ = ‖π‖ ^ ((n : ℤ) * m + ((l : ℕ) : ℤ)) := by
    intro l hl
    obtain ⟨m, hm⟩ := hvalK (c l) hl
    exact ⟨m, by rw [hxval, Algebra.smul_def, norm_mul, hm, norm_pow,
      ← zpow_natCast ‖π‖ (l : ℕ), ← zpow_add₀ (ne_of_gt hπ0)]⟩
  have hczero : ∀ l : Fin n, x l ≠ 0 → c l ≠ 0 := by
    intro l hx hcl
    exact hx (by rw [hxval, hcl, zero_smul])
  have hpair : ∀ i j : Fin n, i ≠ j → x i ≠ 0 → x j ≠ 0 → ‖x i‖ ≠ ‖x j‖ := by
    intro i j hij hi hj
    obtain ⟨mi, hmi⟩ := hexp i (hczero i hi)
    obtain ⟨mj, hmj⟩ := hexp j (hczero j hj)
    rw [hmi, hmj]
    intro hcon
    have h2 : (n : ℤ) * mi + ((i : ℕ) : ℤ) = (n : ℤ) * mj + ((j : ℕ) : ℤ) :=
      (zpow_right_inj₀ hπ0 hπne).mp hcon
    refine TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt
      (fun h => hij (Fin.ext h)) ⟨mj - mi, ?_⟩
    linarith [mul_sub (n : ℤ) mj mi]
  intro l
  have := norm_le_norm_sum_of_pairwise_ne' x hpair l
  rwa [hsum] at this

/-- ★★★★**`h` は元を `‖π‖^t` だけしか動かさない** —— `‖h z − z‖ ≤ ‖z‖·‖π‖^t`。

★これが Serre の `v(σz − z) ≥ v(z) + i_σ` のノルム版であり、★狭義単調(§3)の要である。
★証明は基底展開(`exists_coeff_norm_le`)と超距離だけ。

★在庫の測定: `‖w^m − π^m‖ ≤ ‖π‖^{m−1}·‖w − π‖` は
★**`GainedTowerStep.norm_pow_sub_pow_le`(159 行)に在った**(自作しかけて衝突で気づいた)。 -/
theorem norm_sub_apply_le_mul [FiniteDimensional K M] {π : M} {n t : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) (z : M) :
    ‖h z - z‖ ≤ ‖z‖ * ‖π‖ ^ t := by
  classical
  obtain ⟨c, hc, hnorm⟩ := exists_coeff_norm_le hπ0 hπ1 hn hvalK z
  have hzz : h z - z = ∑ l : Fin n, (c l • ((h π) ^ (l : ℕ) - π ^ (l : ℕ))) := by
    rw [← hc, map_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def, Algebra.smul_def, map_mul, AlgEquiv.commutes,
      map_pow, mul_sub]
  rw [hzz]
  refine norm_sum_le_of_forall_le _ _ (by positivity) (fun l _ => ?_)
  rcases Nat.eq_zero_or_pos (l : ℕ) with hl0 | hl1
  · rw [hl0]
    simp [mul_nonneg (norm_nonneg z) (by positivity : (0:ℝ) ≤ ‖π‖ ^ t)]
  · have hstep : ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖ ≤ ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) := by
      rw [← hbr]
      exact ABC3.Found.PGC.GainedTowerStep.norm_pow_sub_pow_le hπ0 (le_of_eq (hiso π)) _
    have hcomb : ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) = ‖π‖ ^ (l : ℕ) * ‖π‖ ^ t := by
      rw [← pow_add, ← pow_add]
      congr 1
      omega
    rw [Algebra.smul_def, norm_mul]
    calc ‖algebraMap K M (c l)‖ * ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖
        ≤ ‖algebraMap K M (c l)‖ * (‖π‖ ^ (l : ℕ) * ‖π‖ ^ t) := by
          rw [← hcomb]; exact mul_le_mul_of_nonneg_left hstep (norm_nonneg _)
      _ = (‖algebraMap K M (c l)‖ * ‖π‖ ^ (l : ℕ)) * ‖π‖ ^ t := by ring
      _ ≤ ‖z‖ * ‖π‖ ^ t := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          have := hnorm l
          rwa [Algebra.smul_def, norm_mul, norm_pow] at this

end Expansion

end JumpMono

end ABC3.Found.PGC
