/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartToProj
import ABC3.Meta.Claim

/-!
# ★★★★★★★ample なら**有限個**の切断で覆える —— 段 E2 の前半（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★これは台帳の段 E2 の前半である

`IsAmple M`（`AmpleDef.lean`、[Stacks] 28.27.1）は

> `X` が準コンパクトで、各点 `x` に `n ≥ 1` と `s ∈ Γ(X, M^{⊗n})` があって
> `x ∈ X_s` かつ `X_s` がアフィン

という**各点ごとの**条件である。★射 `X ⟶ ℙᴺ_R` を作るには
**有限個**の切断で `X` を覆う必要がある。

★★本ファイルは準コンパクト性からその有限化を取る。

## ★★★機構

各点の `X_s` は開なので `{X_s}` は `X` の開被覆である。
★`isCompact_univ.elim_finite_subcover` で有限部分被覆が取れる。

## ★残っている段（明示）

★★**次数を揃える段（段 E2 の後半）は本ファイルに無い**。
有限被覆の切断はそれぞれ別の `M^{⊗n_j}` に属するので、
`L = lcm(n_j)` として `s_j^{⊗(L/n_j)} ∈ Γ(X, M^{⊗L})` に置き換える必要がある。

★★★そこで要るのは **`X_{s^{⊗k}} = X_s`** である
——`TrivTensor.lean` の `X_{s⊗t} = X_s ∩ X_t` の冪版であり、
`presheafTensorPow` の結合の整合（コヒーレンス）が絡む。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-- ★`IsAmple` の証拠がくれるチャート——テンソル冪の次数と切断の対。 -/
structure AmpleChart (M : X.PresheafOfModules) where
  /-- テンソル冪の次数。 -/
  n : ℕ
  /-- 次数は正である。 -/
  hn : 0 < n
  /-- `M^{⊗n}` の大域切断。 -/
  s : ((presheafTensorPow M n).obj (op ⊤) : Type)

/-- ★チャートの非消失軌跡 `X_s`。 -/
noncomputable def AmpleChart.open' {M : X.PresheafOfModules} (c : AmpleChart M) : X.Opens :=
  nonVanishing (presheafTensorPow M c.n) c.s

/-- ★★★★★★★**ample なら有限個の切断で `X` を覆える**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は「各点の `X_s` は開だから開被覆、`X` は準コンパクト」だけである。
★★これが射 `X ⟶ ℙᴺ_R` の**添字集合**を有限にする段である。 -/
theorem exists_finite_cover_of_isAmple (M : X.PresheafOfModules) (h : IsAmple M) :
    ∃ T : Set (AmpleChart M), T.Finite ∧ (⨆ c ∈ T, (c.open' : Set X)) = Set.univ := by
  obtain ⟨hcpt, hpt⟩ := h
  have hcover : (Set.univ : Set X) ⊆ ⋃ c : AmpleChart M, ((c.open' : Set X)) := by
    intro x _
    obtain ⟨n, hn, s, hx, -⟩ := hpt x
    exact Set.mem_iUnion.2 ⟨⟨n, hn, s⟩, hx⟩
  obtain ⟨T, hT⟩ := isCompact_univ.elim_finite_subcover
    (fun c : AmpleChart M => ((c.open' : Set X))) (fun c => c.open'.2) hcover
  refine ⟨(T : Set (AmpleChart M)), T.finite_toSet, ?_⟩
  refine Set.eq_univ_of_univ_subset (le_trans hT ?_)
  simp

/-! ## ★出典の紐付け(`.src`) -/

def AmpleChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——ample の証拠がくれるチャート)",
    sectionId := "genell-prop-1-4" }

def exists_finite_cover_of_isAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ample なら有限個の切断で X を覆える。次数揃えは含まない)",
    sectionId := "genell-prop-1-4" }

def exists_finite_cover_of_isAmple.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "IsAmple(ample の定義、[Stacks] 28.27.1)"
      (.inProject "ABC3" "ABC3.Found.GenEll.IsAmple") 6,
    .citation "[mathlib]" "IsCompact.elim_finite_subcover"
      (.inMathlib "IsCompact.elim_finite_subcover") 6,
    .implicitStep
      ("★**次数を揃える段(段 E2 の後半)は本ファイルに無い**。" ++
       "有限被覆の切断はそれぞれ別の M^{⊗n_j} に属するので、L = lcm(n_j) として " ++
       "s_j^{⊗(L/n_j)} ∈ Γ(X, M^{⊗L}) に置き換える必要がある。" ++
       "★★そこで要るのは X_{s^{⊗k}} = X_s である" ++
       "——TrivTensor.lean の X_{s⊗t} = X_s ∩ X_t の冪版であり、" ++
       "presheafTensorPow の結合の整合が絡む") 6 ]

end ABC3.Found.GenEll
