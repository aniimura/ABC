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

`residueDegree`(剰余次数、定義のみ)、退化した基点(`x=0`、
`K(0)=K`)での有限次拡大性の基本的な健全性チェック、および
**`adjoinIntegers K 0≃+*𝒪[K.carrier]`という環同型**(`adjoinIntegers
ZeroEquiv`)——これで`residueDegree K 0`が`𝒪_K`自身の剰余体の濃度に
一致すること(`residueDegree_zero`)が示せた。退化した基点は
(`[K(0):K]=1`と合わせ)不分岐性の条件を自明に満たす、という最初の
健全性チェックが完成した。

環同型の構成: `IntermediateField.adjoin K.carrier {0}=⊥`
(`IntermediateField.mem_bot`)と、`K.closure`上のノルムが`K.carrier`
の元では元のノルムに一致すること(mathlibの`norm_algebraMap'`、
探索だけで見つかった)を組み合わせ、`k↦algebraMap K.carrier K.closure
k`という単純な対応が`𝒪[K.carrier]→adjoinIntegers K 0`の環同型を
与えることを直接示した(`RingEquiv.ofBijective`)。剰余体の対応は
`IsLocalRing.ResidueField.mapEquiv`(mathlib)で得る。
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

/-- **`𝒪[K.carrier]→adjoinIntegers K 0`という環準同型**——`k↦algebraMap
K.carrier K.closure k`。値が`adjoinIntegers K 0`(`K.carrier⟮0⟯=⊥`の
ノルム≤1部分環)に入ることは`IntermediateField.mem_bot`+`norm_
algebraMap'`(mathlib、`K.closure`上のノルムが`K.carrier`の元では
元のノルムに一致する)から。 -/
noncomputable def adjoinIntegersZeroRingHom {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    𝒪[K.carrier] →+* adjoinIntegers K 0 where
  toFun k := ⟨⟨algebraMap K.carrier K.closure (k : K.carrier), by
      rw [show IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure) = ⊥ by simp,
        IntermediateField.mem_bot]
      exact ⟨(k : K.carrier), rfl⟩⟩, by
    show ‖algebraMap K.carrier K.closure (k : K.carrier)‖ ≤ 1
    rw [norm_algebraMap']
    exact k.2⟩
  map_one' := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure (1 : K.carrier) = 1
    exact map_one (algebraMap K.carrier K.closure)
  map_mul' a b := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((a : K.carrier) * (b : K.carrier)) =
      algebraMap K.carrier K.closure (a : K.carrier) * algebraMap K.carrier K.closure (b : K.carrier)
    exact map_mul (algebraMap K.carrier K.closure) (a : K.carrier) (b : K.carrier)
  map_zero' := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure (0 : K.carrier) = 0
    exact map_zero (algebraMap K.carrier K.closure)
  map_add' a b := by
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure ((a : K.carrier) + (b : K.carrier)) =
      algebraMap K.carrier K.closure (a : K.carrier) + algebraMap K.carrier K.closure (b : K.carrier)
    exact map_add (algebraMap K.carrier K.closure) (a : K.carrier) (b : K.carrier)

/-- **`adjoinIntegersZeroRingHom`は全単射**——単射性は`algebraMap
K.carrier K.closure`自身の単射性(体の準同型)から、全射性は
`adjoinIntegers K 0`の元が(`IntermediateField.mem_bot`により)
`Set.range (algebraMap K.carrier K.closure)`に入ることと、そのノルム
条件が`norm_algebraMap'`で引き戻せることから。 -/
theorem adjoinIntegersZeroRingHom_bijective {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Function.Bijective (adjoinIntegersZeroRingHom K) := by
  constructor
  · intro a b hab
    apply Subtype.ext
    have h1 : algebraMap K.carrier K.closure (a : K.carrier) = algebraMap K.carrier K.closure (b : K.carrier) :=
      congrArg (fun z : adjoinIntegers K 0 =>
        ((z : IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure)) : K.closure)) hab
    exact (algebraMap K.carrier K.closure).injective h1
  · intro y
    have hmem : (((y : adjoinIntegers K 0) : IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure)) :
        K.closure) ∈ (⊥ : IntermediateField K.carrier K.closure) :=
      (show IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure) = ⊥ by simp) ▸
        (y : IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure)).2
    obtain ⟨k', hk'⟩ := IntermediateField.mem_bot.mp hmem
    have hnorm : ‖k'‖ ≤ 1 := by
      rw [← norm_algebraMap' K.closure k', hk']
      exact y.2
    refine ⟨⟨k', hnorm⟩, ?_⟩
    apply Subtype.ext; apply Subtype.ext
    show algebraMap K.carrier K.closure k' = _
    exact hk'

/-- **`𝒪[K.carrier]≃+*adjoinIntegers K 0`という環同型**——退化した
基点`x=0`では、`K(x)`の整数環が`𝒪_K`自身と(自然に)一致する。 -/
noncomputable def adjoinIntegersZeroEquiv {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    𝒪[K.carrier] ≃+* adjoinIntegers K 0 :=
  RingEquiv.ofBijective (adjoinIntegersZeroRingHom K) (adjoinIntegersZeroRingHom_bijective K)

/-- **退化した基点`x=0`での剰余次数**: `residueDegree K 0`は`𝒪_K`
自身の剰余体の濃度に一致する——`adjoinIntegersZeroEquiv`と
`IsLocalRing.ResidueField.mapEquiv`(mathlib)を組み合わせるだけ。
`[K(0):K]=1`(`finiteDimensional_adjoin_zero`)と合わせ、退化した
基点が不分岐性の条件(剰余次数`=`剰余体濃度の`[K(x):K]`乗)を自明に
満たすことの、最初の健全性チェック。 -/
theorem residueDegree_zero {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure))] :
    residueDegree K 0 = Nat.card (IsLocalRing.ResidueField (𝒪[K.carrier])) := by
  haveI := isLocalRing_adjoinIntegers K 0
  unfold residueDegree
  rw [(Nat.card_congr (IsLocalRing.ResidueField.mapEquiv (adjoinIntegersZeroEquiv K)).toEquiv).symm]

/-! ## 拡大体の整数環の分岐理論への地ならし

`adjoinIntegers K x`は、実は**拡大体`K.carrier⟮x⟯`に同じ`𝒪[...]`記法を
適用したものと`rfl`で一致する**(実測済み)——これにより`Valued`系の
mathlib機構がそのまま使える。以下はその活用で、mathlibの分岐理論
(`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing`、基本等式
`e·f=[L:K]`)を`𝒪[K.carrier]→adjoinIntegers K x`に適用するために要る
instanceを順に埋めていく作業(Hensel's lemmaの自作——mathlibの
`ℤ_[p]`版は約450行——を回避する筋)。 -/

/-- **拡大体`K.carrier⟮x⟯`の付値は非自明**——`p`の像のノルムが`1`未満
(`norm_natCast_p_lt_one`+`norm_algebraMap'`)かつ`0`でない
(`Valuation.ne_zero_iff`+標数`0`)ことから。基礎体では
`(rankOne K).toIsNontrivial`が自動的に見つかるが、拡大体では
instance探索が届かないため明示的に構成する。 -/
theorem isNontrivial_valued_adjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    (Valued.v (R := IntermediateField.adjoin K.carrier ({x} : Set K.closure))).IsNontrivial := by
  haveI : CharZero K.carrier := charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : CharZero (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    charZero_of_injective_algebraMap
      (algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))).injective
  have hpcast : ((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) =
      algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ((p : ℕ) : K.carrier) := by
    push_cast; rfl
  have hnormlt : ‖((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ < 1 := by
    rw [hpcast,
      show ‖algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ((p : ℕ) : K.carrier)‖ =
        ‖((p : ℕ) : K.carrier)‖ from
        norm_algebraMap' (↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ((p : ℕ) : K.carrier)]
    exact norm_natCast_p_lt_one K
  constructor
  refine ⟨((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure)), ?_, ?_⟩
  · rw [Valuation.ne_zero_iff]
    exact_mod_cast (Nat.cast_ne_zero (R := IntermediateField.adjoin K.carrier ({x} : Set K.closure))).mpr
      (Nat.Prime.ne_zero Fact.out)
  · have hv : Valued.v ((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) =
        ‖((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖₊ := NNReal.eq rfl
    rw [hv]
    intro hcon
    have hone : ‖((p : ℕ) : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1 := by
      have := congrArg NNReal.toReal hcon
      simpa using this
    rw [hone] at hnormlt
    exact absurd hnormlt (lt_irrefl 1)

/-- **`adjoinIntegers K x`は離散付値環**——基礎体での`valuationRing_isDVR`
と同じ`Valued.integer.isDiscreteValuationRing_of_compactSpace`に、
`compactSpace_adjoinIntegers`(既出)と上の`isNontrivial_valued_adjoin`
を与えるだけ。これで`IsDedekindDomain (adjoinIntegers K x)`も従い、
分岐理論の基本等式を使う準備の半分が整う(残るは
`Module.Finite (𝒪[K.carrier]) (adjoinIntegers K x)`)。

★配管(記録): `compactSpace_adjoinIntegers K x`は`CompactSpace
(adjoinIntegers K x)`という形だが、補題側が要求するのは
`CompactSpace 𝒪[K.carrier⟮x⟯]`——両者は`rfl`で一致するのに**instance
探索は`def`の壁を越えられない**(`tools/lean-idioms.md` #23/#31 の
類型)。`haveI : CompactSpace 𝒪[...] := compactSpace_adjoinIntegers K x`
と**補題側の形で書いて`:=`で渡す**(defeqは`exact`が受け入れる)ことで
解決する。 -/
theorem isDiscreteValuationRing_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    IsDiscreteValuationRing (adjoinIntegers K x) := by
  haveI : CompactSpace (𝒪[IntermediateField.adjoin K.carrier ({x} : Set K.closure)]) :=
    compactSpace_adjoinIntegers K x
  haveI := isNontrivial_valued_adjoin K x
  exact Valued.integer.isDiscreteValuationRing_of_compactSpace
    (K := IntermediateField.adjoin K.carrier ({x} : Set K.closure))

end ABC3.Found.PGC
