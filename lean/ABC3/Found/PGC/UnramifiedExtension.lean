import ABC3.Found.PGC.AdjoinIntegers

/-!
# 不分岐拡大への出発点(スケルトン、`sorry` 無し・現時点では定義のみ)

節目(5)(古典的Lubin-Tate理論の射影極限、完全分岐アーベル拡大の部分は
完成済み——`memory/pgc-lubin-tate-existence-progress.md`)を`/goal`の
`Proposition 1.2`(および実は`Proposition 2.1`も、`memory/pgc-prop12-
reciprocity-gap.md`参照)へ接続するには、古典的局所類体論の相互律
`Γ_K^ab≅(K^×)^∧`の**もう半分**——`K`の最大不分岐拡大`K^ur`と
`Gal(K^ur/K)≅Ẑ`(Frobenius元が生成する)——が要る。この部分は
mathlib・本プロジェクトともに完全に未着手(2026-09-05実測、
`pgc-prop12-reciprocity-gap.md`に記録済み)。

本ファイルは、その最初のスケルトン(CLAUDE.mdの進め方: スケルトンで
依存グラフを作成後、葉から形式化する)——**定義のみ**を置く。証明を
要する部分(存在性・一意性・Galois群の構造)はまだ`sorry`すら書いて
いない(定義には証明が要らないため、G8/G9のnonvacuous検査の対象外)。

## 数学的な見取り図(証明はまだ無い、方針のみ)

- `x∈K.closure`(`K(x)/K`有限)が**不分岐**であるとは、剰余次数
  `f=[残余体(adjoinIntegers K x):残余体(𝒪_K)]`が拡大次数`[K(x):K]`
  に一致すること(古典的な`e·f=[K(x):K]`のうち`e=1`と同値、
  `Skeleton/PGC/Section1Cor13.lean::IsUnramifiedAt`と同じ特徴づけを
  `Interface`経由ではなく実在の`K`・`adjoinIntegers`で与える)。
- **存在性**(未着手): 各`n≥1`に対し、次数`n`の不分岐拡大が存在する
  ——候補となるmathlibの道具: `HenselianLocalRing`(`𝒪_K`は
  `IsAdicComplete`から完備、ゆえにHenselian——のはずだが本プロジェクト
  でのinstance確認はまだ)・`FiniteField.isSplittingField_sub`
  (`X^{q^n}-X`の分解体が`F_{q^n}`)・分離的多項式の分解体は不分岐
  という一般論(`Algebra.isUnramifiedAt_iff_map_eq`の分離性条件)。
- **一意性**(未着手): 次数`n`の不分岐拡大は(`K.closure`内で)ただ
  1つ——残余体の拡大`F_{q^n}/F_q`の一意性(`GaloisField`・
  `FiniteField.algEquivOfCardEq`)から。
- **`K^ur`とそのGalois群**(未着手): `K^ur:=⋃n,(次数nの不分岐拡大)`、
  `Gal(K^ur/K)≅Gal(F̄_q/F_q)`(剰余体還元による)、
  `Gal(F̄_q/F_q)≅Ẑ`(Frobenius`x↦x^q`が位相的生成元)。

## 現時点で確立したもの

`residueDegree`(剰余次数、定義のみ)と、退化した基点(`x=0`、
`K(0)=K`)での有限次拡大性の基本的な健全性チェック(`FiniteDimensional`
instanceが取れること)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `x:K.closure`(`K(x)/K`有限)の**剰余次数**——`adjoinIntegers K x`
の剰余体の(有限とは限らない)濃度。`𝒪_K`自身の剰余体の濃度の
`[K(x):K]`乗に一致するとき(かつそのときに限り)`x`は不分岐——この
判定条件自体はまだ`Prop`として切り出していない(次の一歩)。 -/
noncomputable def residueDegree {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] : ℕ :=
  haveI := isLocalRing_adjoinIntegers K x
  Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x))

/-- 退化した基点`x=0`では`K(x)=K`(有限次拡大性の健全性チェック、
`finrank=1`のケース)——不分岐性の議論を始める前の最小の確認。 -/
theorem finiteDimensional_adjoin_zero {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure)) := by
  rw [show IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure) = ⊥ by simp]
  infer_instance

end ABC3.Found.PGC
