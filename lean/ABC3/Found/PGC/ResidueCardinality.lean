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
| `card_congr` | `residueCard_congr` | 本ファイル(2026-09-05 追加) |

## ★2026-09-05: `card_congr`(同型不変性)を追加した

`Check/PGC/Prop12Degenerate.lean`(第 1012)が、`card`・`isPrimePow` だけの旧形では
[pGC] Proposition 1.2 の形式化が偽であることを示した。落ちていた条件は**同型不変性**。
`Interface` 側に場を足したので、実物側でそれを証明する必要が生じた——それが本ファイルの
`norm_algEquiv_carrier` → `integerRingEquiv` → `residueCard_congr` の3段である。

筋(`Found/PGC/UnramifiedExtension.lean::norm_algEquiv` の中間体版と同じ形):

1. ℚ_p-代数同型 `e : K ≃ₐ[ℚ_[p]] K'` は**スペクトルノルムを保つ**。
   スペクトルノルムは `spectralValue (minpoly ℚ_[p] x)` であり、`minpoly` は
   ℚ_p-代数同型で不変(`minpoly.algEquiv_eq`)だから。
   ★中間体版が使う `spectralNorm_eq_of_equiv` は `Gal(L/K)`(同じ `L`)専用なので
   ここでは使えない。`minpoly` に降りると2つの型をまたげる。
2. ノルムを保つので、整数環(= 閉単位球、`Valued.integer.mem_iff`)を
   整数環に写す環同型 `integerRingEquiv` を誘導する。
3. 環同型は剰余体の濃度を保つ(`Found/ResidueFieldFinite.lean::card_residueField_eq_of_ringEquiv`)。

極大イデアルの対応を明示的に作る必要は無かった——`IsLocalRing.ResidueField.mapEquiv`
(`card_residueField_eq_of_ringEquiv` の中身)が局所環の同型から直に剰余体の同型を作る。

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
open scoped NormedField Valued

attribute [local instance] Algebra.IsAlgebraic.of_finite

variable (p : ℕ) [Fact p.Prime]

/-! ## 同型不変性(`Interface` の `card_congr` を満たすため) -/

/-- **ℚ_p-代数同型はスペクトルノルムを保つ**。

スペクトルノルムは `spectralValue (minpoly ℚ_[p] x)` であり、最小多項式は
ℚ_p-代数同型で不変(`minpoly.algEquiv_eq`——これは**始域と終域の型が違ってよい**)。

`Found/PGC/UnramifiedExtension.lean::norm_algEquiv` は同じ主張の中間体版だが、
そちらが使う `spectralNorm_eq_of_equiv` は `Gal(L/K)`(同じ `L` の自己同型)専用で、
2つの `PAdicLocalField p` をまたげない。 -/
theorem norm_algEquiv_carrier {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
    (e : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) (x : K.carrier) : ‖e x‖ = ‖x‖ := by
  rw [NormedAlgebra.norm_eq_spectralNorm ℚ_[p] (e x), NormedAlgebra.norm_eq_spectralNorm ℚ_[p] x]
  show spectralValue (minpoly ℚ_[p] (e x)) = spectralValue (minpoly ℚ_[p] x)
  rw [minpoly.algEquiv_eq e x]

/-- ノルムを保つので、`e` は整数環(= 閉単位球)の環同型を誘導する。 -/
noncomputable def integerRingEquiv {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
    (e : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) : 𝒪[K.carrier] ≃+* 𝒪[K'.carrier] where
  toFun z := ⟨e z, Valued.integer.mem_iff.mpr
    (by rw [norm_algEquiv_carrier e]; exact Valued.integer.mem_iff.mp z.2)⟩
  invFun z := ⟨e.symm z, Valued.integer.mem_iff.mpr
    (by rw [norm_algEquiv_carrier e.symm]; exact Valued.integer.mem_iff.mp z.2)⟩
  left_inv z := by apply Subtype.ext; simp
  right_inv z := by apply Subtype.ext; simp
  map_mul' a b := by apply Subtype.ext; simp
  map_add' a b := by apply Subtype.ext; simp

/-- **★★★剰余体の元の個数は ℚ_p-代数同型で不変**。

`Interface` の `ResidueCardinality.card_congr` の実体。
整数環の環同型(`integerRingEquiv`)から `IsLocalRing.ResidueField.mapEquiv` 経由で
剰余体の同型が出る(`Found/ResidueFieldFinite.lean::card_residueField_eq_of_ringEquiv`)。 -/
theorem residueCard_congr {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
    (e : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) : residueCard K = residueCard K' :=
  card_residueField_eq_of_ringEquiv (integerRingEquiv e)

/-- **`ResidueCardinality` の本物**。`Interface` の仮説を落とすための witness。 -/
noncomputable def realResidueCardinality : ResidueCardinality p where
  card := residueCard
  isPrimePow := residueCard_isPrimePow
  card_congr := residueCard_congr

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
