/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Topology.Algebra.Category.ProfiniteGrp.Limits
import Mathlib.GroupTheory.PGroup
import Mathlib.Data.ZMod.Basic

/-!
# pro-`l` 群 —— チェーン `prol` の葉 `prol-def`

★`ResearchPaper/frdi-decomposition.json` の `prol` チェーンの葉。
最終目標は [FrdI] `Definition 2.8, (ii)`。

## ★測定(2026-08-18、探索範囲つき)

★**mathlib に pro-`l` 群は無い**。
探索: `lean/.lake/packages/mathlib/Mathlib/` 全体を
`IsProPGroup|ProPGroup|pro-p group` で grep して **0 件**。
★`ProfiniteGrp`(圏・極限表示)と `IsPGroup` は**ある**ので、
定義は「すべての有限商が `l` 群」と書けば済む。

## ★中身

| | 宣言 |
|---|---|
| 定義 | `IsProL` |
| 自明群 | `isProL_of_subsingleton` |
| 離散なら通常の `l` 群と一致 | `isProL_iff_isPGroup_of_discrete` |
| ★**負の対照** | `not_isProL_three_Z2`(`ℤ/2` は pro-3 でない) |
| 正の対照 | `isProL_two_Z2` |

★負の対照を置くのは、定義が**空虚でない**ことを機械で示すためである
(`tools/check.mjs` 冒頭 B 群の趣旨)。
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

/-- ★★**pro-`l` 群** —— すべての有限商(開正規部分群による商)が `l` 群であるような
副有限群。

★mathlib に無いので我々が定義する(上の測定)。 -/
def IsProL (l : ℕ) (M : ProfiniteGrp.{u}) : Prop :=
  ∀ U : OpenNormalSubgroup M, IsPGroup l (M ⧸ U.toSubgroup)

theorem isProL_iff (l : ℕ) (M : ProfiniteGrp.{u}) :
    IsProL l M ↔ ∀ (U : OpenNormalSubgroup M) (x : M ⧸ U.toSubgroup), ∃ k, x ^ l ^ k = 1 :=
  Iff.rfl

/-- 自明群は任意の `l` について pro-`l`。 -/
theorem isProL_of_subsingleton (l : ℕ) (M : ProfiniteGrp.{u}) [Subsingleton M] : IsProL l M := by
  intro U x
  obtain ⟨y, rfl⟩ := QuotientGroup.mk_surjective x
  rw [Subsingleton.elim y 1]
  exact ⟨0, by simp⟩

/-- 離散な副有限群では `⊥` が開正規部分群。 -/
def botONS (M : ProfiniteGrp.{u}) [DiscreteTopology M] : OpenNormalSubgroup M :=
  { toSubgroup := ⊥, isOpen' := isOpen_discrete _, isNormal' := inferInstance }

set_option maxHeartbeats 1000000 in
/-- ★**離散な場合は通常の `l` 群と一致する** —— 定義が正しい一般化であることの確認。 -/
theorem isProL_iff_isPGroup_of_discrete (l : ℕ) (M : ProfiniteGrp.{u}) [DiscreteTopology M] :
    IsProL l M ↔ IsPGroup l M := by
  constructor
  · intro h x
    obtain ⟨k, hk⟩ := h (botONS M) (QuotientGroup.mk x)
    refine ⟨k, ?_⟩
    rw [← QuotientGroup.mk_pow] at hk
    exact Subgroup.mem_bot.mp
      ((QuotientGroup.eq_one_iff (N := (botONS M).toSubgroup) (x ^ l ^ k)).mp hk)
  · intro h U
    exact h.to_quotient U.toSubgroup

/-! ## ★負の対照 —— 定義は空虚ではない -/

section NegControl

local instance instTopZ2 : TopologicalSpace (Multiplicative (ZMod 2)) := ⊥
local instance instDiscZ2 : DiscreteTopology (Multiplicative (ZMod 2)) := ⟨rfl⟩

/-- 対照に使う副有限群 `ℤ/2`(離散)。 -/
noncomputable def Z2 : ProfiniteGrp.{0} := ProfiniteGrp.of (Multiplicative (ZMod 2))

instance : DiscreteTopology Z2 := instDiscZ2

theorem not_isPGroup_three_Z2 : ¬ IsPGroup 3 (Multiplicative (ZMod 2)) := by
  intro h
  obtain ⟨k, hk⟩ := h (Multiplicative.ofAdd 1)
  have h1 : ((3:ℕ) ^ k : ZMod 2) = 0 := by
    have := congrArg Multiplicative.toAdd hk
    simpa [nsmul_eq_mul] using this
  have h3 : ((3:ℕ) ^ k : ZMod 2) = 1 := by
    push_cast
    rw [show ((3 : ZMod 2)) = 1 from by decide, one_pow]
  rw [h3] at h1
  exact one_ne_zero h1

/-- ★★**負の対照** —— `ℤ/2` は pro-3 ではない。 -/
theorem not_isProL_three_Z2 : ¬ IsProL 3 Z2 := fun h =>
  not_isPGroup_three_Z2 ((isProL_iff_isPGroup_of_discrete 3 Z2).mp h)

/-- ★正の対照 —— `ℤ/2` は pro-2 である。 -/
theorem isProL_two_Z2 : IsProL 2 Z2 := by
  refine (isProL_iff_isPGroup_of_discrete 2 Z2).mpr ?_
  intro x
  refine ⟨1, ?_⟩
  show x ^ (2:ℕ) = 1
  have h : ∀ y : Multiplicative (ZMod 2), y ^ (2:ℕ) = 1 := by decide
  exact h x

end NegControl

end ABC3.Found.ProL
