import ABC3.Skeleton.PGC.Setup
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic
import Mathlib.LinearAlgebra.FiniteDimensional.Defs

/-!
# 単純拡大 `K.carrier⟮x⟯` それ自体を新たな `PAdicLocalField p` として再構成

`sorry` 無し。

## 動機

`Found/PGC/LubinTateDistinguishedSeparable.lean` で `Λ_n`(Lubin-Tate の
`n`-段の全捩れ点)の元 `x` について、`x` を添加した単純拡大 `K.carrier⟮x⟯`
が `K.carrier` 上有限次であること・完備であることまで確立した。次の
目標は `PowerSeries.aeval` で `[a]_f` を `x` へ実際に評価することだが、
そのためには評価先の環に `IsLinearTopology S S` が要る——これは**体**
`K.closure`・`K.carrier⟮x⟯` それ自身では成り立たない(体のイデアルは
`{0}`と全体のみ)。古典的な Lubin-Tate 理論でも `[a]_f` は「体」でなく
「付値環」の間の写像として評価される。

★鍵となる観察: `K.carrier⟮x⟯` は `K.carrier` 上有限次(`x` の整性
から)であり、推移律(`FiniteDimensional.trans`)で `ℚ_[p]` 上も有限次。
つまり **`K.carrier⟮x⟯` 自身が `PAdicLocalField p` の条件を満たす**——
`Field`・`Algebra ℚ_[p] _`は`AlgebraicClosure`/`IntermediateField`の
一般論から自動的に手に入り、`FiniteDimensional ℚ_[p] _` だけ
`FiniteDimensional.trans` で明示的に組み立てればよい。

これにより、`K` に対して `Found/PGC/LocalFieldNorm.lean`・
`Found/PGC/ValuationRingDVR.lean`・`ValuationRingComplete.lean` が
与えていた機構(`𝒪[K.carrier]`・`Valued K.carrier`・その adic 完備性・
`IsLinearTopology` の土台となる `Ideal.isLinearTopology`)が、**この
新しい `PAdicLocalField p` インスタンスへそのまま流用できる**。

## 未解決の注意点(2026-09-04 時点、次に確かめるべきこと)

`IntermediateField.adjoin K.carrier {x}` は `K.closure` の部分体なので、
mathlib の一般論により(`K.closure` に載せた `NormedField` から)
**自動的に**独自の `NormedField` 構造を持つ。一方、ここで作る
`adjoinPAdicLocalField K x` を実際に `PAdicLocalField p` として使うと、
`Found/PGC/LocalFieldNorm.lean` の `normedField` scoped instance が
`spectralNorm.normedField ℚ_[p] _` で**別の**ノルム構造を与える。両者は
(スペクトルノルムの延長の一意性から)数学的には一致するはずだが、
**definitionally 一致するとは限らない**——両方が同時にスコープに
入ると instance diamond になりうる。使う側で `open scoped
ABC3.Found.PGC` する際は、この型の上でどちらのノルムを使っているかを
明示的に管理する必要がある(または両者の一致を橋渡しする補題を先に
用意する)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

/-- `x : K.closure` を `K.carrier` に添加した単純拡大 `K.carrier⟮x⟯`
自身を、新たな `PAdicLocalField p` として組み立てる——`Field`・
`Algebra ℚ_[p] _` は一般論から自動的に手に入り、
`FiniteDimensional ℚ_[p] _` だけ `FiniteDimensional.trans` で明示的に
組み立てる。 -/
noncomputable def adjoinPAdicLocalField {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    PAdicLocalField p where
  carrier := IntermediateField.adjoin K.carrier ({x} : Set K.closure)
  isFinite :=
    FiniteDimensional.trans ℚ_[p] K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))

end ABC3.Found.PGC
