/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPerfFactorial
import ABC3.Found.Divisor.ArithPhiPrime

/-!
# `Φ(L)` は perf-factorial —— [FrdI] `Example 6.3` の「one verifies immediately」

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★何をするか

`ArithPerfFactorial.lean` は「`Γ ⊆ (S →₀ ℝ)` が
**座標ごとに閉じ**、各素点で局所群が**離散か全体**」という枠で
`IsPerfFactorialWith` を組んだ。ここに `Γ = arithDivGroup L` を当てはめる。

★確かめるのは `IsCoordwiseR` だけである(他の 2 条は `ArithPhiPrime.lean` で済んでいる)——
`arithDivGroup L` の条件は「各非アルキメデス成分が `log(N v)` の整数倍」という
**座標ごとの条件**なので、成分を取り出しても保たれる。

★★これで原文の

> one verifies immediately that Φ(L) ̸= 0 is perf-factorial

が閉じる。
-/

namespace ABC3.Found.Divisor

open NumberField Finsupp ABC3.Found.FrdI

universe u

variable {L : Type u} [Field L] [NumberField L]

/-- ★★**`arithDivGroup L` は座標ごとに閉じている**。

★条件が「各非アルキメデス成分が `log(N v)` の整数倍」という座標ごとの形だからである。 -/
theorem isCoordwiseR_arithDivGroup : IsCoordwiseR (arithDivGroup L) := by
  intro y hy s
  cases s with
  | inl w =>
      obtain ⟨n, hn⟩ := hy w
      exact single_inl_mem_arithDivGroup w hn
  | inr w => exact single_inr_mem_arithDivGroup w _

/-- ★★★★★★**[FrdI] Example 6.3 —— `Φ(L)` は perf-factorial**。

★原文の「one verifies immediately that `Φ(L) ≠ 0` is perf-factorial」。

★★幾何(`Example 6.1`)との違いは `Definition 2.4, (i)` の `realScale` である ——
幾何では `M_p ≃+ ℕ` なので空虚だったが、算術ではアルキメデス素点で
`M_p ≃+ ℝ≥0` なので**実際に確かめる必要がある**
(`ArithPerfFactorial.lean` の `image_iotaR_univ_of_full`)。 -/
theorem isPerfFactorial_arithEff : IsPerfFactorial (effR (arithDivGroup L)) :=
  isPerfFactorial_effR isCoordwiseR_arithDivGroup isLocallyMonoprimeR_arithDivGroup

/-- ★★★★★**族 `ι` を明示した形**。 -/
theorem isPerfFactorialWith_arithEff :
    IsPerfFactorialWith (effR (arithDivGroup L))
      (iotaR (isGenSubgroupR_of_isLocallyMonoprimeR
        (isLocallyMonoprimeR_arithDivGroup (L := L)))) :=
  isPerfFactorialWith_effR isCoordwiseR_arithDivGroup isLocallyMonoprimeR_arithDivGroup

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Example 6.3` の「`Φ(L)` は perf-factorial」。 -/
def isPerfFactorial_arithEff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ(L) は perf-factorial",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
