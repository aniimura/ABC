import ABC3.Interface.PGC.LocalFieldData
import ABC3.Found.PGC.LocalFieldNorm
import ABC3.Found.PGC.QpResidueField

/-!
# Track B — `Interface.PGC.ResidueCardinality` の実装

`Interface/PGC/LocalFieldData.lean` の `ResidueCardinality` に対する、
名前の付いた実装。**`sorry` 無し**。

材料は全て `Found/` にある:

| 成分 | 実体 | 置き場所 |
|---|---|---|
| `card` | `residueCard K = Nat.card 𝓀[K]` | `Found/PGC/LocalFieldNorm.lean` |
| `isPrimePow` | `residueCard_isPrimePow` | 同上(一般部分は `Found/ResidueFieldFinite.lean`) |

## ★G2 witness をここに置く理由

`ResidueCardinality.nonvacuous` は `structure` と同じファイル(`Interface/`)ではなく
**ここ**にある。`Interface/` から `Found/` を import すると、`Interface` を import する
`Skeleton/` が `Found` を推移的に引くことになり、「実装が無くても statement を書ける」
という条件付き形式化の要点が壊れるため(PLAN §3)。

`check.mjs` の G2 は宣言名を木全体から探すので、この配置で通る。
`Interface → Found` の import は `check.mjs` が禁止している(fixture D23/D24)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable (p : ℕ) [Fact p.Prime]

/-- **`ResidueCardinality` の本物**。`Interface` の仮説を落とすための witness。 -/
noncomputable def realResidueCardinality : ResidueCardinality p where
  card := residueCard
  isPrimePow := residueCard_isPrimePow

@[simp] theorem realResidueCardinality_card (K : PAdicLocalField p) :
    (realResidueCardinality p).card K = residueCard K := rfl

/-- ★**値が出る最初の点**: `ℚ_[p]` において `card = p`。

`Found/PGC/QpResidueField.lean` の計算。これで `realResidueCardinality` は
「型が付くだけ」の対象ではなくなった。 -/
@[simp] theorem realResidueCardinality_card_selfField :
    (realResidueCardinality p).card (selfField p) = p :=
  residueCard_selfField p

end ABC3.Found.PGC

namespace ABC3.Interface.PGC

variable (p : ℕ) [Fact p.Prime]

/-- **G2 非空虚 witness**。`sorry` 無し。

`structure` は `Interface/PGC/LocalFieldData.lean` にあるが、witness は
実装を要するのでこちら(`Found/`)で `Interface.PGC` 名前空間に足す。
向きの理由は本ファイル冒頭の docstring。

これにより、`ResidueCardinality` を仮説に持つ `Skeleton/` の主張
(`Skeleton/PGC/Section1.lean`・`Section1Cor13.lean`)は
**空虚でないことが確定した**——PLAN §3 が言う「Track B が実例を1つ供給した瞬間、
それに依存する Track A の statement 群が一斉に非空虚性の検査を受ける」の一例目。

`∀ K` の側が空虚でないこと(= `PAdicLocalField p` が空でないこと)の検査は
`Check/PGC/ResidueCardinalityNondegenerate.lean`。 -/
theorem ResidueCardinality.nonvacuous : Nonempty (ResidueCardinality p) :=
  ⟨ABC3.Found.PGC.realResidueCardinality p⟩

end ABC3.Interface.PGC
