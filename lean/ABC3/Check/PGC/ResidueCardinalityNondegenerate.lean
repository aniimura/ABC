import ABC3.Interface.PGC.LocalFieldData
import ABC3.Check.PGC.Section1Discriminating

/-!
# `ResidueCardinality` の非退化検査

`Interface/PGC/LocalFieldData.lean` の `ResidueCardinality` の docstring は
「`isPrimePow` が無いと `card := fun _ => 0` でも通ってしまい、内容が消える」と
主張している。これは**散文のまま**では検査されていなかったので、機械にかける。

**これは原典の主張ではない**(我々のモデルについての事実)ので `.src` を持たない。

## なぜ discharge の直後にこれが要るか

`isPrimePow : ∀ K, ∃ f, 0 < f ∧ card K = p ^ f` は **`PAdicLocalField p` 上の全称**である。
もし `PAdicLocalField p` が空なら、この条件は空虚に真になり、
`Found` が供給した witness も「何も言っていない」ことになる——
`ABC3/Meta/Calibration.lean` が実演した失敗と同型。

`Check/PGC/Section1Discriminating.lean` の `base : PAdicLocalField p`(= `ℚ_[p]` 自身)が
その空虚化を塞ぐ。以下の2つは `base` を実際に消費しているので、
`PAdicLocalField p` が空なら**証明できない**。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- 退化した `card`(恒等的に 0)は `ResidueCardinality` を満たさない。
すなわち `isPrimePow` は飾りではない。 -/
theorem no_zero_residueCardinality (r : ResidueCardinality p) : r.card base ≠ 0 := by
  obtain ⟨f, hf, hcard⟩ := r.isPrimePow base
  rw [hcard]
  exact pow_ne_zero f (Fact.out : p.Prime).pos.ne'

/-- より強く、`card` は 1 にもなれない(`0 < f` が効いている箇所)。
`0 < f` を落とすと `f = 0` で `p ^ 0 = 1` が通ってしまう——原文の
「q = p^f 元からなる体」は自明体を含意しないので、ここは主張の一部。 -/
theorem no_one_residueCardinality (r : ResidueCardinality p) : r.card base ≠ 1 := by
  obtain ⟨f, hf, hcard⟩ := r.isPrimePow base
  rw [hcard]
  exact Nat.ne_of_gt (Nat.one_lt_pow hf.ne' (Fact.out : p.Prime).one_lt)

end ABC3.Check.PGC
