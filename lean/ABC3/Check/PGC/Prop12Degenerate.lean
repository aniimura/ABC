import ABC3.Skeleton.PGC.Section1
import ABC3.Found.PGC.GaloisTransferContinuous
import ABC3.Found.PGC.QpResidueField
import Mathlib.Algebra.Algebra.TransferInstance
import Mathlib.Algebra.Field.TransferInstance
import Mathlib.Algebra.Module.TransferInstance

/-!
# [pGC] Proposition 1.2 の現在の形は偽だった——`ResidueCardinality` は同型不変でない

原文 (pGC p.3):

> The number q of elements in the residue field of O[scr]_K, and well as the absolute
> degree [K : Q[bb]_p] of K, can be recovered entirely group-theoretically from Γ_K.

`Skeleton/PGC/Section1.lean` の現在の形は

```
∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal
```

であり、`RD : ResidueCardinality p`(`Interface/PGC/LocalFieldData.lean`)は
`card : PAdicLocalField p → ℕ` と `isPrimePow`(`card K = p^f`, `f > 0`)だけを持つ
自由なデータである。**`card` が同型不変であることは要求されていない。**

`card` の定義域は p進局所体の同型類ではなく `PAdicLocalField p` の**項**なので、
同じ体を別の項として2つ作れば、`card` はその2つに別の値を割り当ててよい。
本ファイルはそれを実行する。

## 反例(`prop_1_2_statement_false`、`sorry` 無し)

台を `ℚ_[p]` のままにして、体構造だけを全単射 `x ↦ -x` に沿って移送する:

* `TwistedQp p := ℚ_[p]`(型としては `rfl` で等しい)に、
  `Equiv.field (Equiv.neg ℚ_[p])` で捻った `Field` を入れる。
  演算は `x *' y = -(x*y)`、単位元は `1' = -1`(加法は変わらない)。
* `Equiv.algebra` / `Equiv.algEquiv` で `ℚ_[p]`-代数構造と
  `twistedAlgEquiv : TwistedQp p ≃ₐ[ℚ_[p]] ℚ_[p]` を得る。これは**本物の p進局所体**である
  (`ℚ_[p]` と ℚ_p-代数として同型、次数 1)。
* 既存の `Found/PGC/GaloisTransferContinuous.lean::galContinuousMulEquiv` に
  その代数同型を渡して `α : Γ_{TwistedQp} ≅ Γ_{ℚ_p}`(位相群の同型)を得る。
* 病的な `RD` は「台が `ℚ_[p]` に等しいとき、その `1` が標準の `1` か」で分岐する
  (`OneIsStandard`)。捻った側は `1' = -1 ≠ 1` なので**判定が異なる**。
  `card := if OneIsStandard K then p else p^2` は `isPrimePow` を両分岐で満たす(f = 1 と f = 2)。

α を当てると `p^2 = p` になり、`p` が素数なので矛盾。次数成分は両方 1 なので効かない。

## ★これは原典 Proposition 1.2 の反証ではない

反例の2体は ℚ_p-代数として**同型**であり、q も次数も本来は等しい。
落ちているのは原典の数学ではなく、**我々の形式化**である:

> `ResidueCardinality.card` は `PAdicLocalField p` の項の関数であって、
> 同型類の関数ではない。同型不変性を課していないので、
> 「Γ_K から回復できる」は最初から成り立ちようがない。

すなわち「落とした条件は主張を偽にするか自明にする」の7例目。
落とした条件は**同型不変性**である。

## 修理の方向(本体セッションへ)

`ResidueCardinality` に不変性の場を1つ足せばこの反例は消える:

```
card_congr : ∀ {K K' : PAdicLocalField p}, (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) → card K = card K'
```

実物 `Found/PGC/ResidueCardinality.lean::realResidueCardinality` はこれを満たすはずである
(剰余体の元の個数は ℚ_p-代数同型で保たれる)ので、非空虚性は失われない。
★ただし `card_congr` を足すと Proposition 1.2 は「α : Γ_K ≅ Γ_K′ から体の同型を作れ」
という原典本来の内容に戻る——**易しくはならない**。

## ★2026-08-14 の記録の訂正

`Check/PGC/RefutationAttempts.lean` は Proposition 1.2 について
「反証できなかった」と記録し、その理由を2つ挙げていた。両方とも今日の反例で更新される:

1. 「`ULift` による同型コピーは `Field (ULift α)` が mathlib に無いので作れない」
   ——`Mathlib/Algebra/Field/ULift.lean` に `ULift.field` が在るので、この理由は今は誤り。
   ただし同ファイルが続けて挙げる「`ℚ_[p] ≠ ULift ℚ_[p]`(型の非等号)は Lean で証明できない」
   の方は**今も正しい**。本ファイルが `ULift` を使わないのはそのためである
   ——台の型は `ℚ_[p]` のまま**構造だけ**を捻れば、2つの項の非等号は
   `1' = -1 ≠ 1` から**証明できる**。
2. 「`refutation_reduces_to_alpha` により、反証には次数の違う2体の間の α が要る」
   ——これは**次数成分**についての正しい観察だが、`RD` 成分には効かない。
   本反例の2体は次数が等しく(ともに 1)、落ちるのは `card` 成分である。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Found.PGC ABC3.Interface.PGC
open scoped Classical

variable (p : ℕ) [Fact p.Prime]

/-! ## 1. 体構造だけを捻った p進局所体 -/

/-- 捻りに使う全単射 `x ↦ -x`。**環準同型ではない**(積を保たない)ので、
これに沿って移送した体構造は元の体構造と異なる。 -/
noncomputable abbrev negEquiv : ℚ_[p] ≃ ℚ_[p] := Equiv.neg ℚ_[p]

/-- `negEquiv` に沿って移送した `Field ℚ_[p]`。`1' = -1`、`x *' y = -(x*y)`。 -/
@[implicit_reducible] noncomputable def negField : Field ℚ_[p] := Equiv.field (negEquiv p)

/-- 同じ移送で得られる `ℚ_[p]`-代数構造。

★`Semiring ℚ_[p]` のインスタンスが標準と捻りの2つ在る状況なので、
`Equiv.algebra` の引数はすべて明示する(`_` に任せると単一化が捻った方を拾い、
「synthesized type class instance is not definitionally equal」で落ちる)。 -/
@[implicit_reducible] noncomputable def negAlgebra :
    @Algebra ℚ_[p] ℚ_[p] (by infer_instance) (negField p).toSemiring :=
  @Equiv.algebra ℚ_[p] ℚ_[p] ℚ_[p] (by infer_instance) (negEquiv p) (by infer_instance)
    (by infer_instance)

/-- 移送が与える ℚ_p-代数同型(台の型を付け替える前の生の形)。 -/
noncomputable def negAlgEquivRaw :
    @AlgEquiv ℚ_[p] ℚ_[p] ℚ_[p] (by infer_instance) (negField p).toSemiring (by infer_instance)
      (negAlgebra p) (by infer_instance) :=
  @Equiv.algEquiv ℚ_[p] ℚ_[p] ℚ_[p] (by infer_instance) (negEquiv p) (by infer_instance)
    (by infer_instance)

/-- 捻った体の台。**型としては `ℚ_[p]` そのもの**(`rfl` で等しい)。

★型同義語にするのが要点である。`ℚ_[p]` のままだとインスタンス探索が常に標準構造を拾い、
`AlgEquiv.symm` や `Module.Finite.equiv` が
「synthesized type class instance is not definitionally equal」で落ちる。
一方 `structure` で包む(`ULift` など)と `TwistedQp p = ℚ_[p]` が証明できなくなり、
判定関数 `OneIsStandard` が使えなくなる。**`def` による型同義語だけが両立する。** -/
def TwistedQp (p : ℕ) [Fact p.Prime] : Type := ℚ_[p]

noncomputable instance instFieldTwistedQp : Field (TwistedQp p) := negField p

noncomputable instance instAlgebraTwistedQp : Algebra ℚ_[p] (TwistedQp p) := negAlgebra p

/-- **捻った体は `ℚ_[p]` と ℚ_p-代数として同型**——本物の p進局所体である。 -/
noncomputable def twistedAlgEquiv : TwistedQp p ≃ₐ[ℚ_[p]] ℚ_[p] := negAlgEquivRaw p

instance instFiniteTwistedQp : Module.Finite ℚ_[p] (TwistedQp p) :=
  Module.Finite.equiv (twistedAlgEquiv p).symm.toLinearEquiv

/-- 台は `ℚ_[p]` そのもの、体構造だけを捻った p進局所体。 -/
noncomputable def twistedField : PAdicLocalField p where
  carrier := TwistedQp p
  isField := instFieldTwistedQp p
  isAlgebra := instAlgebraTwistedQp p
  isFinite := instFiniteTwistedQp p

/-- **2つの体の絶対 Galois 群の間の連続同型**。

`Found/PGC/GaloisTransferContinuous.lean::galContinuousMulEquiv`(第 969)に
上の代数同型を渡すだけ。★`RefutationAttempts.lean` が「構成できない」と書いていた α が、
ここでは**体の同型から**得られる——2体が同型だからである。 -/
noncomputable def twistedGalEquiv :
    ContinuousMulEquiv (twistedField p).absGal (selfField p).absGal :=
  galContinuousMulEquiv (twistedAlgEquiv p)

/-! ## 2. 2つの項を区別する述語 -/

/-- 「台が `ℚ_[p]` に等しいなら、その `1` は標準の `1` である」。

★同型不変ではない述語である(同型な2体で値が変わる)。それが要点で、
`ResidueCardinality.card` にはそういう関数を禁じる条件が無い。 -/
def OneIsStandard (K : PAdicLocalField p) : Prop :=
  ∀ h : K.carrier = ℚ_[p], cast h (1 : K.carrier) = (1 : ℚ_[p])

theorem oneIsStandard_selfField : OneIsStandard p (selfField p) := by
  intro h; exact cast_eq h 1

/-- 捻った体の `1` は `-1` なので、判定は逆に出る。 -/
theorem not_oneIsStandard_twisted : ¬ OneIsStandard p (twistedField p) := by
  intro hc
  have h1 : @Eq ℚ_[p] (-(1 : ℚ_[p])) 1 := hc rfl
  have h2 : (2 : ℚ_[p]) = 0 := by linear_combination -h1
  exact two_ne_zero h2

/-- **2つの項は等しくない**——`PAdicLocalField p` の項の等号は台の型だけでなく
体構造の等号も要求するから。 -/
theorem twistedField_ne_selfField : twistedField p ≠ selfField p := by
  intro h
  refine not_oneIsStandard_twisted p ?_
  rw [h]
  exact oneIsStandard_selfField p

/-! ## 3. 病的な `ResidueCardinality` -/

/-- 病的な剰余体の元の個数——**同型な2体に別の値を返す**。
`isPrimePow`(原文が課している唯一の条件)は両分岐で満たしている。 -/
noncomputable def badRD : ResidueCardinality p where
  card K := if OneIsStandard p K then p else p ^ 2
  isPrimePow K := by
    by_cases hK : OneIsStandard p K
    · exact ⟨1, one_pos, by rw [if_pos hK, pow_one]⟩
    · exact ⟨2, two_pos, by rw [if_neg hK]⟩

theorem badRD_selfField : (badRD p).card (selfField p) = p :=
  if_pos (oneIsStandard_selfField p)

theorem badRD_twisted : (badRD p).card (twistedField p) = p ^ 2 :=
  if_neg (not_oneIsStandard_twisted p)

/-! ## 4. 反証 -/

/-- **★★★★★★★[pGC] Proposition 1.2 の現在の形(自由な `RD`)は偽**。

`Skeleton/PGC/Section1.lean::residueCard_and_degree_recoverable` は
`(RD : ResidueCardinality p)` を仮説に取る条件付き形式化だが、
`RD` に同型不変性を課していないため、**どの `RD` でも成り立つ**という形にはならない。 -/
theorem prop_1_2_statement_false (p : ℕ) [Fact p.Prime] :
    ¬ (∀ RD : ResidueCardinality p, (residueCardAndDegreeObject RD).RecoverableFromAbsGal) := by
  intro h
  have key := h (badRD p) (twistedGalEquiv p)
  have h1 : (badRD p).card (twistedField p) = (badRD p).card (selfField p) :=
    congrArg Prod.fst key
  rw [badRD_twisted, badRD_selfField] at h1
  have hp : 1 < p := (Fact.out : p.Prime).one_lt
  have h2 : p * p = p * 1 := by rw [mul_one, ← pow_two]; exact h1
  have h3 : p = 1 := Nat.eq_of_mul_eq_mul_left (by omega) h2
  omega

#print axioms prop_1_2_statement_false

end ABC3.Check.PGC
