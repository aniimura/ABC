import ABC3.Found.PGC.HasseArfCongruence

/-!
# 骨組み(作業中) —— 分岐群のノルム言語での特徴づけ
-/

namespace ABC3.Found.PGC

namespace RamNormBridge

open Finset WildBreak JumpMono

section Core

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- ★**指数の議論** —— `0 < l < n` で `‖c·π^l‖ ≤ 1` なら `‖c·π^l‖ ≤ ‖π‖`。

`‖c·π^l‖ = ‖π‖^{n·m + l}` で、`≤ 1` から指数 `≥ 0`、`0 < l < n` から `n ∤ l` ゆえ指数 `≠ 0`。
★これが「単元の展開では `l = 0` の項だけが大きい」ことの中身である。 -/
theorem norm_le_norm_pi_of_le_one {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {c : K} {l : ℕ} (hl0 : 0 < l) (hln : l < n)
    (hle : ‖algebraMap K M c * π ^ l‖ ≤ 1) :
    ‖algebraMap K M c * π ^ l‖ ≤ ‖π‖ := by
  rcases eq_or_ne c 0 with hc0 | hc0
  · rw [hc0, map_zero, zero_mul, norm_zero]; exact le_of_lt hπ0
  · obtain ⟨m, hm⟩ := hvalK c hc0
    have hval : ‖algebraMap K M c * π ^ l‖ = ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := by
      rw [norm_mul, hm, norm_pow, ← zpow_natCast ‖π‖ l, ← zpow_add₀ (ne_of_gt hπ0)]
    have hE0 : 0 ≤ (n : ℤ) * m + (l : ℤ) := by
      by_contra hcon
      rw [not_le] at hcon
      have hgt : (1:ℝ) < ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := one_lt_zpow_of_neg₀ hπ0 hπ1 hcon
      rw [← hval] at hgt
      linarith
    have hEne : (n : ℤ) * m + (l : ℤ) ≠ 0 := by
      intro hE
      refine TotallyRamifiedLayer.not_dvd_sub_of_lt hln (by omega : 0 < n)
        (by omega : l ≠ 0) ⟨-m, ?_⟩
      push_cast
      linarith [mul_neg (n : ℤ) m]
    rw [hval]
    calc ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) ≤ ‖π‖ ^ (1 : ℤ) :=
          zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)
      _ = ‖π‖ := zpow_one _

/-- ★★★★**(a) 単元での 1 つ分の改良** —— `‖z‖ ≤ 1` なら `‖h z − z‖ ≤ ‖π‖^{t+1}`。

★前波 `JumpMono.norm_sub_apply_le_mul` は `‖z‖·‖π‖^t` までしか出さず、
`‖z‖ ≤ 1` では `‖π‖^t`(1 つ足りない)であった。★`l = 0` の項が消えることと、
`l ≥ 1` の項が `‖π‖` 以下であること(上の指数の議論)で 1 つ取り戻す。

★これで `h ∈ G_t`(真の値)がノルムの言葉で言える。 -/
theorem norm_sub_apply_le_of_norm_le_one [FiniteDimensional K M] {π : M} {n t : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) {z : M} (hz : ‖z‖ ≤ 1) :
    ‖h z - z‖ ≤ ‖π‖ ^ (t + 1) := by
  classical
  obtain ⟨c, hc, hnorm⟩ := exists_coeff_norm_le hπ0 hπ1 hn hvalK z
  have hzz : h z - z = ∑ l : Fin n, (c l • ((h π) ^ (l : ℕ) - π ^ (l : ℕ))) := by
    rw [← hc, map_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def, Algebra.smul_def, map_mul, AlgEquiv.commutes, map_pow, mul_sub]
  rw [hzz]
  refine norm_sum_le_of_forall_le _ _ (by positivity) (fun l _ => ?_)
  rcases Nat.eq_zero_or_pos (l : ℕ) with hl0 | hl1
  · rw [hl0]
    simp [pow_nonneg (norm_nonneg π) (t + 1)]
  · have hle1 : ‖algebraMap K M (c l) * π ^ (l : ℕ)‖ ≤ 1 := by
      have h0 := hnorm l
      rw [Algebra.smul_def] at h0
      exact le_trans h0 hz
    have hcl : ‖algebraMap K M (c l) * π ^ (l : ℕ)‖ ≤ ‖π‖ :=
      norm_le_norm_pi_of_le_one hπ0 hπ1 hvalK hl1 l.isLt hle1
    have hstep : ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖ ≤ ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) := by
      rw [← hbr]
      exact GainedTowerStep.norm_pow_sub_pow_le hπ0 (le_of_eq (hiso π)) _
    have hcomb : ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) = ‖π‖ ^ (l : ℕ) * ‖π‖ ^ t := by
      rw [← pow_add, ← pow_add]
      congr 1
      omega
    rw [Algebra.smul_def, norm_mul]
    calc ‖algebraMap K M (c l)‖ * ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖
        ≤ ‖algebraMap K M (c l)‖ * (‖π‖ ^ (l : ℕ) * ‖π‖ ^ t) := by
          rw [← hcomb]; exact mul_le_mul_of_nonneg_left hstep (norm_nonneg _)
      _ = (‖algebraMap K M (c l)‖ * ‖π‖ ^ (l : ℕ)) * ‖π‖ ^ t := by ring
      _ ≤ ‖π‖ * ‖π‖ ^ t := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          rwa [← norm_pow, ← norm_mul]
      _ = ‖π‖ ^ (t + 1) := by rw [pow_succ]; ring

/-- ★★★★★**分岐群のノルム言語での特徴づけ**(両向き)。

    (∀ z, ‖z‖ ≤ 1 → ‖h z − z‖ ≤ ‖π‖^{i+1})  ↔  i ≤ t        (`‖hπ − π‖ = ‖π‖^{t+1}`)

★左辺は `LowerRamificationGroup.lean:270` の
`σ ∈ lowerRamificationGroup B G i ↔ ∀ x : B, σ•x − x ∈ 𝔪^{i+1}` の**ノルム版**である。
⇒ ★**`h ∈ G_i ⟺ i ≤ t`**。前波の「1 つずれる」は解消した。 -/
theorem mem_ramification_iff [FiniteDimensional K M] {π : M} {n t i : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) :
    (∀ z : M, ‖z‖ ≤ 1 → ‖h z - z‖ ≤ ‖π‖ ^ (i + 1)) ↔ i ≤ t := by
  constructor
  · intro hall
    have hpi := hall π (le_of_lt hπ1)
    rw [hbr] at hpi
    by_contra hcon
    rw [not_le] at hcon
    have hlt : t + 1 < i + 1 := by omega
    have hstrict : ‖π‖ ^ (i + 1) < ‖π‖ ^ (t + 1) :=
      pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 hlt
    linarith
  · intro hit z hz
    refine le_trans (norm_sub_apply_le_of_norm_le_one h hiso hπ0 hπ1 hn hvalK hbr hz) ?_
    exact pow_le_pow_of_le_one (le_of_lt hπ0) (le_of_lt hπ1) (by omega)

end Core

end RamNormBridge

end ABC3.Found.PGC
