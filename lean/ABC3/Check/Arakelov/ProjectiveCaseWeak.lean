import Mathlib.Topology.Compactness.Compact
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# 負の対照 —— **「連続かつ単射」では埋め込みにならない**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

## ★★★★何を確かめるか

(C2) `ProjectiveModelData` の `projectiveCase` は 2026-08-20 以前

    連続かつ単射で像が閉有界 ⟹ 定義域はコンパクト

を要求していた。★★**これは偽である。**

★★★本ファイルは反例を**機械で**作る:

| 部品 | 取り方 |
|---|---|
| `X` | `Option ℕ` に**離散位相** |
| `f` | `none ↦ 0`、`some n ↦ 1/(n+1)` |

★`f` は連続(離散だから)で単射、像は `{0} ∪ {1/(n+1)}` でコンパクト。
★★しかし `X` は無限離散なのでコンパクトでない。

★★★★したがって `projectiveCase` は **`IsInducing`**(埋め込み)を要求せねばならない
——`Found/GenEll/ArcModel.lean` が実際に持っているのもそれである(§9-406)。
-/

namespace ABC3.Check.Arakelov

open Filter Topology

/-- ★★★★★**連続単射で像がコンパクトでも、定義域はコンパクトとは限らない**。

★この定理が通ること自体が、旧 `projectiveCase`(連続+単射)が**偽**であったことの証明である。 -/
theorem exists_cont_inj_compact_range_not_compact :
    ∃ (X : Type) (t : TopologicalSpace X) (f : X → ℝ),
      @Continuous X ℝ t _ f ∧ Function.Injective f ∧ IsCompact (Set.range f) ∧
      ¬ @CompactSpace X t := by
  classical
  refine ⟨Option ℕ, ⊥, fun x => Option.elim x 0 (fun n => 1 / (n + 1 : ℝ)), ?_, ?_, ?_, ?_⟩
  · letI : TopologicalSpace (Option ℕ) := ⊥
    haveI : DiscreteTopology (Option ℕ) := ⟨rfl⟩
    exact continuous_of_discreteTopology
  · rintro (_ | n) (_ | m) h <;> simp only [Option.elim] at h
    · rfl
    · exact absurd h.symm (ne_of_gt (by positivity))
    · exact absurd h (ne_of_gt (by positivity))
    · have hpos1 : (0 : ℝ) < (n : ℝ) + 1 := by positivity
      have hpos2 : (0 : ℝ) < (m : ℝ) + 1 := by positivity
      have hnm : (n : ℝ) = (m : ℝ) := by
        field_simp at h
        linarith
      exact congrArg some (Nat.cast_injective hnm)
  · have hr : Set.range (fun x : Option ℕ => Option.elim x (0 : ℝ) (fun n => 1 / (n + 1 : ℝ)))
        = insert 0 (Set.range (fun n : ℕ => 1 / (n + 1 : ℝ))) := by
      ext y
      constructor
      · rintro ⟨(_ | n), rfl⟩
        · exact Set.mem_insert _ _
        · exact Set.mem_insert_of_mem _ ⟨n, rfl⟩
      · rintro (rfl | ⟨n, rfl⟩)
        · exact ⟨none, rfl⟩
        · exact ⟨some n, rfl⟩
    rw [hr]
    refine Filter.Tendsto.isCompact_insert_range_of_cofinite ?_
    rw [Nat.cofinite_eq_atTop]
    exact tendsto_one_div_add_atTop_nhds_zero_nat
  · letI : TopologicalSpace (Option ℕ) := ⊥
    haveI : DiscreteTopology (Option ℕ) := ⟨rfl⟩
    intro hc
    have hfin : Finite (Option ℕ) := finite_of_compact_of_discrete
    have hfin2 : Finite ℕ :=
      Finite.of_injective (fun n => (some n : Option ℕ)) (by intro a b h; simpa using h)
    exact (Infinite.not_finite hfin2)

end ABC3.Check.Arakelov
