import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [BC] Theorem 2.2.3 (Faltings) の posit(`Interface`)

原典: O. Brinon, B. Conrad, *CMI Summer School Notes on p-adic Hodge Theory
(Preliminary Version)*, 2009(`papers.json` 短縮タグ `BC`)。物理 p.13。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

★[LocProP] Proposition 2.4 が引く「Comparison Theorem」(原文は [Falt2]
Theorem 6.2 を指すが、[Falt2] は paywall で未取得。[BC] の Theorem 2.2.3
は同じ内容([Falt2] の元の結果そのもの)を述べているので、等価な代替
典拠として `Interface/BC` を立てる——`Skeleton/LocProP/Section2.lean` の
`hodgeTateExactSeq.needs` 冒頭のコメント参照)。

スキーム論・エタールコホモロジー・Hodge コホモロジーそのものは posit
しない(Falt1 と同じ流儀)。左辺 `C_K ⊗ H^n_et(X_K̄,Qp)` と右辺
`⊕_q(C_K(-q) ⊗ H^{n-q}(X,Ω^q_{X/K}))` を抽象的な加法群の対とみなし、
定理が主張する同型だけを posit する。
-/

namespace ABC3.Interface.BC

universe u

structure ComparisonIsoSetup where
  /-- 左辺: `C_K ⊗_{Qp} H^n_et(X_K̄, Qp)`。 -/
  Het : Type u
  HetGrp : AddCommGroup Het
  /-- 右辺: `⊕_q (C_K(-q) ⊗_K H^{n-q}(X, Ω^q_{X/K}))`。 -/
  HHodgeSum : Type u
  HHodgeSumGrp : AddCommGroup HHodgeSum
  /-- **[BC] Theorem 2.2.3 (Faltings)**、式 (2.2.1) の同型。 -/
  thm223 : Het ≃+ HHodgeSum

@[reducible] def ComparisonIsoSetup.example : ComparisonIsoSetup.{0} where
  Het := ℤ
  HetGrp := inferInstance
  HHodgeSum := ℤ
  HHodgeSumGrp := inferInstance
  thm223 := AddEquiv.refl ℤ

@[reducible] def ComparisonIsoSetup.nonvacuous : Nonempty (ComparisonIsoSetup.{0}) :=
  ⟨ComparisonIsoSetup.example⟩

end ABC3.Interface.BC
