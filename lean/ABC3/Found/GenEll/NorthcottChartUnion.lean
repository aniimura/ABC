/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGlobalChart
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★チャートごとに有限なら全体でも有限（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★これは何か —— 段 E2 を Northcott に効かせる形

`§9-888`（段 E3h）は**チャート 1 枚分**の Northcott であった。
★`X = ⋃_i X_{s_i}`（段 E2）で全体へ広げるとき要るのは、
集合論の 1 行——**有限個の有限集合の和は有限**——だけである。

    `{p | ht p ≤ C} = ⋃_i {p | チャート = i ∧ ht p ≤ C}`

★★添字が `Fin (N+1)`（有限）なので `Set.finite_iUnion` が効く。

## ★★★これで `Proposition 1.4, (iv)` の骨格が揃う

| 段 | 内容 | 場所 |
|---|---|---|
| チャート 1 枚 | `X_{s_i}` の点について Northcott | `§9-888` |
| **全体** | **有限個の和で `X` へ** | ★★**本ファイル** |

★★★★残るのは「各点をどのチャートに割り当てるか」の配管で、
それは `§9-C2b`（`exists_chart_range`——体値の点はどれかのチャートに入る）が与える。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★★★★★有限個の有限集合の和 -/

/-- ★★★★★★★★★★**チャートごとに有限なら全体でも有限**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`X = ⋃_i X_{s_i}`（段 E2）を Northcott に効かせる形である。
★★機構は `Set.finite_iUnion`（添字が有限）だけ。 -/
theorem finite_of_finite_charts {P : Type} {ι : Type} [Finite ι]
    (ht : P → ℝ) (C : ℝ) (chart : P → ι)
    (hfin : ∀ i : ι, {p : P | chart p = i ∧ ht p ≤ C}.Finite) :
    {p : P | ht p ≤ C}.Finite := by
  have hcover : {p : P | ht p ≤ C} = ⋃ i : ι, {p : P | chart p = i ∧ ht p ≤ C} := by
    ext p
    simp only [Set.mem_setOf_eq, Set.mem_iUnion]
    constructor
    · intro hp; exact ⟨chart p, rfl, hp⟩
    · rintro ⟨i, -, hp⟩; exact hp
  rw [hcover]
  exact Set.finite_iUnion hfin

/-- ★★**部分集合の版** —— チャートごとの有限性を部分集合で受ける形。 -/
theorem finite_of_finite_charts_subset {P : Type} {ι : Type} [Finite ι]
    (ht : P → ℝ) (C : ℝ) (chart : P → ι) (S : ι → Set P)
    (hmem : ∀ p, p ∈ S (chart p))
    (hfin : ∀ i : ι, {p : P | p ∈ S i ∧ ht p ≤ C}.Finite) :
    {p : P | ht p ≤ C}.Finite := by
  refine finite_of_finite_charts ht C chart (fun i => (hfin i).subset ?_)
  rintro p ⟨hchart, hp⟩
  exact ⟨hchart ▸ hmem p, hp⟩

/-! ## ★出典の紐付け(`.src`) -/

def finite_of_finite_charts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートごとに有限なら全体でも有限)",
    sectionId := "genell-prop-1-4" }

def finite_of_finite_charts_subset.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートごとの有限性を部分集合で受ける形)",
    sectionId := "genell-prop-1-4" }

def finite_of_finite_charts.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_globalChart(チャート 1 枚分の Northcott、段 E3h、§9-888)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_globalChart") 2,
    .citation "[mathlib]" "Set.finite_iUnion(有限個の有限集合の和は有限)"
      (.inMathlib "Set.finite_iUnion") 1,
    .implicitStep
      ("★X = ⋃_i X_{s_i}(段 E2)で全体へ広げるとき要るのは、" ++
       "集合論の 1 行——**有限個の有限集合の和は有限**——だけである。" ++
       "添字が Fin (N+1)(有限)なので Set.finite_iUnion が効く") 1,
    .implicitStep
      ("★★残るのは「各点をどのチャートに割り当てるか」の配管で、" ++
       "それは §9-C2b(exists_chart_range——体値の点はどれかのチャートに入る)が与える") 3 ]

end ABC3.Found.GenEll
