import ABC3.Found.PGC.WildBreakLowerBound
import ABC3.Found.PGC.RamificationJumpBound

/-!
# 骨組み(作業中) —— 跳びの上界 `(p−1)t ≤ p·e`
-/

namespace ABC3.Found.PGC

namespace WildBreakUpper

open Finset WildBreak

/-! ## §1 抽象核 —— 共役はすべて等距離 -/

section Equidistant

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★**抽象核** —— `g` が等長なら `‖g^j π − π‖ ≤ ‖g π − π‖`(超距離＋群の性質だけ)。 -/
theorem norm_pow_apply_sub_le (g : M ≃ₐ[K] M) (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) :
    ∀ j : ℕ, ‖(g ^ j) π - π‖ ≤ ‖g π - π‖ := by
  intro j
  induction j with
  | zero => simp
  | succ j ih =>
      have hrw : (g ^ (j + 1)) π - π = g ((g ^ j) π - π) + (g π - π) := by
        rw [map_sub, ← AlgEquiv.mul_apply, ← pow_succ']
        ring
      rw [hrw]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ le_rfl)
      rw [hiso]
      exact ih

/-- ★★**抽象核** —— `orderOf g = p`(素数)で `p ∤ j` なら `‖g^j π − π‖ = ‖g π − π‖`。

★`g^j` も生成元なので上の不等式を**両向き**に使える。★分岐・付値の語彙は 1 語も出ない。 -/
theorem norm_pow_apply_sub_eq {p : ℕ} (hp : p.Prime) (g : M ≃ₐ[K] M) (hg : orderOf g = p)
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) {j : ℕ} (hj : ¬ (p ∣ j)) :
    ‖(g ^ j) π - π‖ = ‖g π - π‖ := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine le_antisymm (norm_pow_apply_sub_le g hiso π j) ?_
  -- `j` の法 `p` 逆元 `k` を取る
  have hjne : ((j : ZMod p)) ≠ 0 := by
    rw [Ne, ZMod.natCast_eq_zero_iff]
    exact hj
  set k : ℕ := ((j : ZMod p)⁻¹).val with hkdef
  have hmod : j * k ≡ 1 [MOD p] := by
    rw [← ZMod.natCast_eq_natCast_iff]
    push_cast
    rw [hkdef, ZMod.natCast_val, ZMod.cast_id]
    field_simp
  have hgjk : (g ^ j) ^ k = g := by
    rw [← pow_mul, ← pow_one g]
    exact (pow_eq_pow_iff_modEq).mpr (by simpa [hg] using hmod)
  have hisoj : ∀ z : M, ‖(g ^ j) z‖ = ‖z‖ := norm_pow_apply g hiso j
  have := norm_pow_apply_sub_le (g ^ j) hisoj π k
  rwa [hgjk] at this

end Equidistant

end WildBreakUpper

end ABC3.Found.PGC
