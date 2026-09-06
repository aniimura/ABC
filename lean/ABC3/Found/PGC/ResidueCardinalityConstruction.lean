import ABC3.Found.PGC.AbsoluteRamification
import ABC3.Found.PGC.ResidueCardinality
import ABC3.Interface.PGC.LocalFieldData

/-!
# `Interface.PGC.ResidueCardinality` を**構成する**

`Interface/PGC/LocalFieldData.lean` の `ResidueCardinality p` は

* `card : PAdicLocalField p → ℕ`(剰余体の元の個数 `q`)
* `isPrimePow : ∀ K, ∃ f, 0 < f ∧ card K = p ^ f`

という自由なデータだった。`Skeleton/PGC/Section1Cor13.lean` の設計メモは

> 退化は消えていない——移動した …… Track B が本物の `ResidueCardinality` を
> 構成した時点で、ここに依存する全ての statement が一斉に非空虚性の検査を受ける。

と書いている。本ファイルはその**本物**を与える:

* `card K := Nat.card 𝓀[K.carrier]`(剰余体の元の個数——`Found/ResidueFieldFinite.lean`
  で有限性は確立済み)
* `isPrimePow` は `Found/PGC/AbsoluteRamification.lean::residueCard_eq_pow`
  (`q = p^{f(K/ℚ_p)}`)と `e·f = [K:ℚ_p] > 0`(したがって `f > 0`)から。
* `card_congr`(2026-09-05 追加、第 1012 の修理)は
  `Found/PGC/ResidueCardinality.lean::residueCard_congr`——
  ℚ_p-代数同型はスペクトルノルムを保つので整数環の環同型を誘導し、剰余体の濃度が一致する。

これで `Proposition 1.2`(`residueCard_and_degree_recoverable`)と
`Corollary 1.3`(`inertia_recoverable`)の仮説 `RD : ResidueCardinality p` に
**具体的な値**が入る——両者はもはや「自由なデータについての条件付き主張」では
なく、実在の量についての主張として読める。

## ★2026-09-06(D13 の実行): 仮説は消えた

上の「仮説 `RD` に具体的な値が入る」は、2026-09-06 に**仮説そのものを外す**形で決着した。
`Skeleton/PGC/Section1.lean` の Proposition 1.2 は `∀ RD : ResidueCardinality p` をやめ
`Found/PGC/ResidueCardinality.lean::realResidueCardinality` に固定してある
(理由: `∀ RD` 版は原典が偽と述べている命題と同値だった。
`Check/PGC/Prop12ForallRD.lean::forall_RD_recoverable_iff_algEquiv`)。
Corollary 1.3 の主語は 2026-09-05 から本ファイルの `residueCardinality` のままである。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- 絶対慣性次数は正——`e·f = [K:ℚ_p] > 0` だから。 -/
theorem absoluteInertiaDegree_pos (K : PAdicLocalField p) : 0 < absoluteInertiaDegree K := by
  have h := absoluteRamificationIndex_mul_absoluteInertiaDegree K
  have hpos : 0 < Module.finrank ℚ_[p] K.carrier := Module.finrank_pos
  rcases Nat.eq_zero_or_pos (absoluteInertiaDegree K) with h0 | h0
  · rw [h0, Nat.mul_zero] at h
    omega
  · exact h0

/-- **★★`ResidueCardinality p` の本物の構成**——剰余体の元の個数。 -/
noncomputable def residueCardinality (p : ℕ) [Fact p.Prime] : ResidueCardinality p where
  card K := Nat.card 𝓀[K.carrier]
  isPrimePow K := ⟨absoluteInertiaDegree K, absoluteInertiaDegree_pos K, residueCard_eq_pow K⟩
  card_congr e := residueCard_congr e

@[simp] theorem residueCardinality_card (K : PAdicLocalField p) :
    (residueCardinality p).card K = Nat.card 𝓀[K.carrier] := rfl

end ABC3.Found.PGC
