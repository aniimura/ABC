import ABC3.Found.PGC.ResidueCardinalityConstruction

/-!
# 有限拡大を `PAdicLocalField p` として構成する / 不分岐性の判定条件

## `adjoinField`

`PAdicLocalField p` は `carrier : Type` に `Field`・`Algebra ℚ_[p]`・
`FiniteDimensional ℚ_[p]` を載せただけの構造(`Skeleton/PGC/Setup.lean`)。
したがって単項拡大 `K(x) ⊆ K̄` はそのまま `PAdicLocalField p` になる
——`Algebra ℚ_[p] ↥K⟮x⟯` は自動、`FiniteDimensional ℚ_[p] ↥K⟮x⟯` は塔
(`Module.Finite.trans`)。

これは `Interface/PGC/LocalFieldData.lean` の `SubgroupCorrespondence`
(開部分群 ↔ 中間体、しかも `PAdicLocalField p` として)へ向かう第一歩。

★**未解決の配管**(次に戻るときの出発点): `(adjoinField K x).carrier` に
`Found/PGC/LocalFieldNorm.lean` の `normedField` インスタンス(ℚ_p 上の
スペクトルノルム)が載るのに対し、`adjoinIntegers K x` が使っているのは
`K.closure` の部分体として `↥K⟮x⟯` が継ぐノルム。両者は
**命題としては等しい**(完備体上のノルム延長の一意性)が **definitional
には別物**なので、`Nat.card 𝓀[(adjoinField K x).carrier] = residueDegree K x`
は `rfl` で通らない。橋渡しには `spectralNorm_unique_field_norm_ext`
(完備体からの体ノルムの一意性)を使う必要がある。

## 不分岐性の判定条件

原文 (pGC p.3) の

> L is unramified over K if and only q_L = q^[Γ_K:H]

を、本プロジェクトの量で書いたもの:`q_L`(= `residueDegree K x`、
`𝒪_{K(x)}` の剰余体の元の個数)が `q_K^{[K(x):K]}` に等しいことと、
`e = 1` は同値。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- **単項拡大を `PAdicLocalField p` として見る**。 -/
noncomputable def adjoinField (K : PAdicLocalField p) (x : K.closure) : PAdicLocalField p where
  carrier := IntermediateField.adjoin K.carrier ({x} : Set K.closure)
  isFinite := by
    haveI : IsScalarTower ℚ_[p] K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by infer_instance
    exact Module.Finite.trans (R := ℚ_[p]) K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))

@[simp] theorem adjoinField_carrier (K : PAdicLocalField p) (x : K.closure) :
    (adjoinField K x).carrier = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := rfl

/-- **★不分岐性の判定条件**——`q_L = q_K^{[L:K]}` ⟺ 不分岐。
原文 p.3 の "L is unramified over K if and only q_L = q^[Γ_K:H]" を、
単項拡大について本プロジェクトの量で書いたもの。 -/
theorem isUnramifiedAdjoin_iff_residueDegree (K : PAdicLocalField p) (x : K.closure) :
    IsUnramifiedAdjoin K x ↔ residueDegree K x
      = Nat.card 𝓀[K.carrier]
        ^ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  constructor
  · intro hu
    rw [residueDegree_eq_residueCard_pow K x, inertiaDegree_eq_finrank_of_isUnramified K x hu]
  · exact isUnramifiedAdjoin_of_residueDegree K x

/-- 判定条件を `Interface` の `residueCardinality` の言葉で書いたもの
(基礎体側だけ——拡大体側は上記の配管が未解決)。 -/
theorem isUnramifiedAdjoin_iff_residueDegree_interface (K : PAdicLocalField p) (x : K.closure) :
    IsUnramifiedAdjoin K x ↔ residueDegree K x
      = (residueCardinality p).card K
        ^ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  isUnramifiedAdjoin_iff_residueDegree K x

end ABC3.Found.PGC
