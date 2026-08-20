import Mathlib.Data.Fintype.EquivFin
import Mathlib.SetTheory.Cardinal.Finite
import Mathlib.Tactic.Ring

/-!
# 負の対照 —— **`≃` では群構造を一切課さない**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

## ★★★★何を確かめるか

(G1) `TorsionStructureData` の `structure_eq` は 2026-08-17 以前

    Nonempty (torsionPoints W n ≃ (ZMod n × ZMod n))

——すなわち**型の全単射**——を要求していた。

★★**有限型では、濃度が等しいだけで全単射は存在する。**
したがって旧 `structure_eq` は `#E[n] = n²`(数え上げ)だけで埋まり、
★★★**群構造は一切課されていなかった**。

★★★★原文 `Theorem 3.8` が主張するのは**群同型**であり、
`GL₂` が表現の行き先になるのはそのためである。★2026-08-17 に `≃+` へ強めた(§9-400)。
-/

namespace ABC3.Check.GaloisRep

/-- ★★★★★**濃度が等しければ全単射は在る**——群構造は要らない。

★この定理が通ること自体が、旧 `structure_eq`(`≃`)が
**数え上げだけで埋まる**ことの証明である。 -/
theorem equiv_of_card_eq {A B : Type} [Finite A] [Finite B] (h : Nat.card A = Nat.card B) :
    Nonempty (A ≃ B) := by
  haveI := Fintype.ofFinite A
  haveI := Fintype.ofFinite B
  refine ⟨Fintype.equivOfCardEq ?_⟩
  simpa [Nat.card_eq_fintype_card] using h

/-- ★★**したがって `(ZMod n × ZMod n) ≃ A` は `#A = n²` から出る**。

★★★これが「弱い posit」の姿である——`≃+` にして初めて群の情報が入る。 -/
theorem equiv_zmod_sq_of_card {A : Type} [Finite A] (n : ℕ) [NeZero n]
    (h : Nat.card A = n ^ 2) : Nonempty ((ZMod n × ZMod n) ≃ A) := by
  refine equiv_of_card_eq ?_
  rw [Nat.card_prod, Nat.card_zmod, h]
  ring

end ABC3.Check.GaloisRep
