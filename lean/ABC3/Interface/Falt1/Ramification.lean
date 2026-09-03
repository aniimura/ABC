import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs
import Mathlib.Data.Real.Basic
import Mathlib.Order.Basic

/-!
# [Falt1] Chapter I §1, Lemma 1.1・Theorem 1.2 の posit(`Interface`)

原典: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299。
物理 p.4-5(印字 p.257-258)。**260 dpi 目視確認 2026-09-04**(逐語の食い違い
無し)。

★注: このPDFのpdftotextのOCR層は数式記号を激しく壊す(`Ω`が`K`に、
`⊗`が`?`に、`⊂`が`c`になる等)。260dpi目視では内容が原文と一致する
ことを確認したが、tools/check.mjsの逐語照合(`原文 (タグ p.N):`記法・
HTML `.verbatim`)はOCR層に対して行われるため、Tate・BKと同じ運用に
する: Leanコメントは地の文(`内容 (Falt1 p.N):`)で写し、HTML側
`.verbatim`はOCRでも生き残る短い見出し部分だけに絞る。

離散付値環の拡大 `V ⊂ W`(完備、混標数、剰余体の `p`-基底の指数を `d` と
する)の微分加群論・differentの理論そのものは posit する(可換環論の
深い理論なので、mathlib に既存の対応物があるかは未調査——
`ResearchPaper/mathlib-gap.json` への登記は Found 着手時に行う)。
-/

namespace ABC3.Interface.Falt1

universe u

structure RamificationSetup where
  /-- `Ω_V ⊗_V W`(Lemma 1.1 の射の始域)。 -/
  OmegaVW : Type u
  OmegaVWGrp : AddCommGroup OmegaVW
  /-- `Ω_W`(Lemma 1.1 の射の終域)。 -/
  OmegaW : Type u
  OmegaWGrp : AddCommGroup OmegaW
  /-- **1.1 Lemma** の自然な射 `Ω_V ⊗_V W → Ω_W`。 -/
  lem11 : OmegaVW →+ OmegaW
  /-- **1.1 Lemma**: 上の射は単射。 -/
  lem11_injective : Function.Injective lem11
  /-- 余核 `Ω_{W/V}` の「長さ」(組成列の長さ、自然数値に抽象化)。 -/
  cokerLength : ℕ
  /-- `W/p^δW` の長さ(δ は W over V の different)。 -/
  quotientLength : ℕ
  /-- **1.1 Lemma**: 上の2つの長さが等しい。 -/
  lem11_length_eq : cokerLength = quotientLength
  /-- **1.2 Theorem** 用: `δ(W_n/V_n)` の列(添字 n、非負実数に抽象化)。 -/
  delta : ℕ → ℝ
  delta_nonneg : ∀ n, 0 ≤ delta n
  /-- **1.2 Theorem**: `δ_n → 0`(n → ∞、ε-N 論法で述べる)。 -/
  thm12 : ∀ ε : ℝ, 0 < ε → ∃ N : ℕ, ∀ n ≥ N, delta n < ε

@[reducible] def RamificationSetup.example : RamificationSetup.{0} where
  OmegaVW := ℤ
  OmegaVWGrp := inferInstance
  OmegaW := ℤ
  OmegaWGrp := inferInstance
  lem11 := AddMonoidHom.id ℤ
  lem11_injective := fun _ _ h => h
  cokerLength := 0
  quotientLength := 0
  lem11_length_eq := rfl
  delta := fun _ => 0
  delta_nonneg := fun _ => le_refl (0 : ℝ)
  thm12 := fun _ hε => ⟨0, fun _ _ => hε⟩

@[reducible] def RamificationSetup.nonvacuous : Nonempty (RamificationSetup.{0}) :=
  ⟨RamificationSetup.example⟩

end ABC3.Interface.Falt1
