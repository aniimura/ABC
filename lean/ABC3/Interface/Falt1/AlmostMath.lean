import ABC3.Meta.Claim
import Mathlib.Algebra.Group.Equiv.Basic
import Mathlib.Algebra.Group.Int.Defs

/-!
# [Falt1] "almost mathematics" の核心概念の共有部品

原典物理p.6(印字p.259)、§2(a)冒頭: "By m we denote the ideal generated
by all powers of p^ε, ε > 0."(m は全ての `p^ε`(ε>0)が生成するイデアル)。
論文全体で頻出する「m で零化される」("annihilated by m")を、この定義に
沿って「ある `e : ℕ` に対し `p^e` 倍写像が零写像になる」と具体化し、
Chapter I の複数の項目(§2 Theorem 2.4、§4 等)で共有する。

almost mathematics の理論そのもの(A[1/p]-加群の圏、almost 同型の圏論的
扱い等)は posit しない——あくまで個々の定理の主張を書き下すための
最小限の道具として使う。
-/

namespace ABC3.Interface.Falt1

universe u

/-- `X` が「m で零化される」ことの具体化。`pPow e` を「p^e 倍写像」の
posit とし、ある `e` でこれが零写像になることを要求する。 -/
structure MTorsionWitness (X : Type u) [AddCommGroup X] where
  pPow : ℕ → X →+ X
  ann : ∃ e : ℕ, pPow e = 0

/-- 自明な(零写像の場合の)証拠: `pPow` を恒等的に零写像とすれば `e = 0` で成立。 -/
@[reducible] def MTorsionWitness.trivial : MTorsionWitness ℤ where
  pPow := fun _ => 0
  ann := ⟨0, rfl⟩

@[reducible] def MTorsionWitness.nonvacuous : Nonempty (MTorsionWitness ℤ) :=
  ⟨MTorsionWitness.trivial⟩

end ABC3.Interface.Falt1
