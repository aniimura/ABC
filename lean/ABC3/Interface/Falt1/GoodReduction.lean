import ABC3.Meta.Claim
import ABC3.Interface.Falt1.AlmostEtale

/-!
# [Falt1] Chapter I §3, Theorem 3.1・3.2 の posit(`Interface`)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)。物理 p.10・12(印字 p.263・265)。
**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

★注: OCR 層は数式記号を激しく壊すため地の文で写す
(Interface/Falt1/Ramification.lean 冒頭の注記と同じ運用)。

Theorem 3.1・3.2 はどちらも「(specific construction) は almost étale
covering である」という結論を持つ——Definition 2.1 と同じ
`AlmostEtaleSetup.isAlmostEtale` 述語を再利用し、それぞれ別インスタンス
(具体的な `R,S,R_∞,S_∞` 等が異なる)として posit する。
-/

namespace ABC3.Interface.Falt1

structure GoodReductionSetup where
  /-- **Theorem 3.1** 用のインスタンス: `S[1/p]` が `R[1/p]` 上 étale なとき、
  `S⊗_R R_∞` の正規化 `S_∞` は `R_∞` の almost étale covering。 -/
  setup31 : AlmostEtaleSetup
  thm31 : setup31.isAlmostEtale
  /-- **Theorem 3.2** 用のインスタンス: `S[1/(pu_1...u_d)]` が
  `R[1/(pu_1...u_d)]` 上 étale なとき、`S_n` たちの合併 `S_∞` は
  `R_∞` の almost étale covering。 -/
  setup32 : AlmostEtaleSetup
  thm32 : setup32.isAlmostEtale

@[reducible] def GoodReductionSetup.example : GoodReductionSetup where
  setup31 := AlmostEtaleSetup.example
  thm31 := trivial
  setup32 := AlmostEtaleSetup.example
  thm32 := trivial

@[reducible] def GoodReductionSetup.nonvacuous : Nonempty GoodReductionSetup :=
  ⟨GoodReductionSetup.example⟩

end ABC3.Interface.Falt1
