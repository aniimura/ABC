import ABC3.Found.PGC.AdjoinIntegers
import Mathlib.RingTheory.DedekindDomain.IntegralClosure
import Mathlib.NumberTheory.RamificationInertia.Basic
import Mathlib.RingTheory.Henselian
import Mathlib.Algebra.Polynomial.Eval.Irreducible
import Mathlib.RingTheory.Polynomial.GaussLemma

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


/-! ## `Module.Finite` への道——スペクトルノルムで「ノルム≤1 ⟺ `𝒪_K`上整」を得る

`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing`(`e·f=[L:K]`)を
使うために最後に残っていた前提が `Module.Finite 𝒪[K.carrier]
(adjoinIntegers K x)`。これは `IsIntegralClosure (adjoinIntegers K x)
𝒪[K.carrier] K.carrier⟮x⟯`(=「`adjoinIntegers K x` はちょうど
`𝒪[K.carrier]` の `K.carrier⟮x⟯` における整閉包」)から
`IsIntegralClosure.finite` で出る。その核心が下の同値。

★発見(2026-09-05): mathlib の **spectral norm** の道具立てで両方向とも
数行で書ける。Hensel の補題も共役の対称式も自分で組む必要が無かった。

* 難しい向き(ノルム≤1 ⟹ 整): `NormedAlgebra.norm_eq_spectralNorm`
  (`‖y‖ = spectralNorm K L y`、`K`が完備なら成立)と、`spectralNorm`の
  定義が `spectralValue (minpoly K y)` であること、そして
  `spectralValue_le_one_iff`(monic なら `spectralValue P ≤ 1 ↔
  ∀n, ‖P.coeff n‖ ≤ 1`)を繋ぐ。すなわち **`‖y‖≤1` なら最小多項式の
  係数がすべて `𝒪_K` に入る**——古典的な「共役のノルムが等しい」議論が
  mathlib 側に既に畳み込まれている。あとは `Polynomial.toSubring` で
  係数を `𝒪_K` に落とすだけ。
* 易しい向き(整 ⟹ ノルム≤1): `norm_root_le_spectralValue` に
  `f := spectralAlgNorm K L` を与える(`spectralAlgNorm_isPowMul`・
  `isNonarchimedean_spectralNorm` が揃っている)。

★配管(記録): `rw [← NormedAlgebra.norm_eq_spectralNorm ...]` は
**instance 経路の違いで発火しない**(`Valued`由来の`NormedField`と
補題側のそれが syntactic に一致しない)。`exact`/`le_of_eq_of_le` に
すると defeq で通る(`tools/lean-idioms.md` #37 の類型)。 -/

/-- **ノルム≤1 ⟹ `𝒪[K.carrier]` 上整**(難しい向き)。最小多項式の
係数がすべてノルム≤1 であることを `spectralValue_le_one_iff` から得て、
`Polynomial.toSubring` で `𝒪[K.carrier]` 係数の monic 多項式に落とす。 -/
theorem isIntegral_of_norm_le_one {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (hy : ‖y‖ ≤ 1) :
    IsIntegral 𝒪[K.carrier] y := by
  have hmonic : (minpoly K.carrier y).Monic := minpoly.monic (IsIntegral.of_finite K.carrier y)
  have hcoeff : ∀ n : ℕ, ‖(minpoly K.carrier y).coeff n‖ ≤ 1 := by
    rw [← spectralValue_le_one_iff hmonic]
    have h1 : spectralNorm K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) y
        = spectralValue (minpoly K.carrier y) := rfl
    rw [← h1]
    exact le_of_eq_of_le (NormedAlgebra.norm_eq_spectralNorm K.carrier y).symm hy
  have hsub : ↑(minpoly K.carrier y).coeffs ⊆ (𝒪[K.carrier] : Set K.carrier) := by
    intro c hc
    obtain ⟨n, _, rfl⟩ := Polynomial.mem_coeffs_iff.mp hc
    show _ ∈ 𝒪[K.carrier]
    rw [Valuation.mem_integer_iff]
    have hv : Valued.v ((minpoly K.carrier y).coeff n)
        = (‖(minpoly K.carrier y).coeff n‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    exact_mod_cast hcoeff n
  refine ⟨(minpoly K.carrier y).toSubring 𝒪[K.carrier] hsub, ?_, ?_⟩
  · exact (Polynomial.monic_toSubring _ _ _).mpr hmonic
  · have h2 : (algebraMap 𝒪[K.carrier]
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        = (algebraMap K.carrier
            (IntermediateField.adjoin K.carrier ({x} : Set K.closure))).comp
          (Subring.subtype 𝒪[K.carrier]) := rfl
    rw [h2, ← Polynomial.eval₂_map, Polynomial.map_toSubring]
    exact minpoly.aeval K.carrier y

/-- **`𝒪[K.carrier]` 上整 ⟹ ノルム≤1**(易しい向き)。整方程式を
`K.carrier` 係数に写し、`norm_root_le_spectralValue` を
`spectralAlgNorm` に適用する。 -/
theorem norm_le_one_of_isIntegral {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hy : IsIntegral 𝒪[K.carrier] y) : ‖y‖ ≤ 1 := by
  obtain ⟨q, hqm, hq⟩ := hy
  have hmap : (q.map (Subring.subtype 𝒪[K.carrier])).Monic := hqm.map _
  have hroot : Polynomial.aeval y (q.map (Subring.subtype 𝒪[K.carrier])) = 0 := by
    rw [Polynomial.aeval_def, Polynomial.eval₂_map]
    exact hq
  have hcoeff : ∀ n : ℕ, ‖(q.map (Subring.subtype 𝒪[K.carrier])).coeff n‖ ≤ 1 := by
    intro n
    rw [Polynomial.coeff_map]
    have h0 : Valued.v ((q.coeff n : K.carrier)) ≤ 1 := (q.coeff n).2
    have hv : Valued.v ((q.coeff n : K.carrier))
        = (‖(q.coeff n : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [hv] at h0
    exact_mod_cast h0
  have hsv : spectralValue (q.map (Subring.subtype 𝒪[K.carrier])) ≤ 1 :=
    (spectralValue_le_one_iff hmap).mpr hcoeff
  have hle : (spectralAlgNorm K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) y
      ≤ spectralValue (q.map (Subring.subtype 𝒪[K.carrier])) :=
    norm_root_le_spectralValue spectralAlgNorm_isPowMul isNonarchimedean_spectralNorm hmap hroot
  rw [spectralAlgNorm_def] at hle
  calc ‖y‖ = spectralNorm K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) y :=
      NormedAlgebra.norm_eq_spectralNorm K.carrier y
    _ ≤ _ := hle
    _ ≤ 1 := hsv

/-- **`adjoinIntegers K x` の元＝`𝒪[K.carrier]` 上整な元**——上の二つを
まとめた同値。`IsIntegralClosure` を作るための材料。 -/
theorem isIntegral_iff_norm_le_one {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (y : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    IsIntegral 𝒪[K.carrier] y ↔ ‖y‖ ≤ 1 :=
  ⟨norm_le_one_of_isIntegral K x y, isIntegral_of_norm_le_one K x y⟩


/-! ## 分岐理論の基本等式 `e·f = [L:K]` ——不分岐性を定義できる地点

上の同値から `IsIntegralClosure (adjoinIntegers K x) 𝒪[K.carrier]
K.carrier⟮x⟯` が出て、`IsIntegralClosure.finite` が
`Module.Finite 𝒪[K.carrier] (adjoinIntegers K x)` を与える——これが
`Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing` に残っていた
最後の前提だった。これで **`e·f = [K(x):K]`** が本プロジェクトの
設定でそのまま使える。

★配管(記録): `Ideal.ramificationIdx`/`inertiaDeg` は
`IsLocalRing (adjoinIntegers K x)` を**主張の型の段階で**要求する
(`haveI` を証明の中に置いても遅い)。`residueDegree` と同じく
`haveI := isLocalRing_adjoinIntegers K x` を **`def` の本体**に置いた
薄いラッパー(`ramificationIndex`・`inertiaDegree`)を挟むことで、
利用側にインスタンス束縛を波及させずに済む。 -/

/-- 基礎体 `K.carrier` の付値も自明でない(`isNontrivial_valued_adjoin`
の基礎体版)。 -/
theorem isNontrivial_valued_carrier {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    (Valued.v (R := K.carrier)).IsNontrivial := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  constructor
  refine ⟨((p : ℕ) : K.carrier), ?_, ?_⟩
  · rw [Valuation.ne_zero_iff]
    exact_mod_cast (Nat.cast_ne_zero (R := K.carrier)).mpr (Nat.Prime.ne_zero Fact.out)
  · have hv : Valued.v ((p : ℕ) : K.carrier) = (‖((p : ℕ) : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    intro hcon
    have hone : ‖((p : ℕ) : K.carrier)‖ = 1 := by
      have := congrArg NNReal.toReal hcon
      simpa using this
    have hlt := norm_natCast_p_lt_one K
    rw [hone] at hlt
    exact absurd hlt (lt_irrefl 1)

/-- 基礎体の整数環 `𝒪[K.carrier]` は離散付値環。`IsNoetherianRing`・
`IsDedekindDomain` はここから従う。 -/
theorem isDiscreteValuationRing_carrierIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsDiscreteValuationRing 𝒪[K.carrier] :=
  haveI := isNontrivial_valued_carrier K
  Valued.integer.isDiscreteValuationRing_of_compactSpace (K := K.carrier)

/-- **`adjoinIntegers K x` は `𝒪[K.carrier]` の `K.carrier⟮x⟯` における
整閉包**——`isIntegral_iff_norm_le_one` をそのまま `IsIntegralClosure`
の形に組み替えたもの。 -/
theorem isIntegralClosure_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    IsIntegralClosure (adjoinIntegers K x) 𝒪[K.carrier]
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  constructor
  · exact Subtype.val_injective
  · intro z
    rw [isIntegral_iff_norm_le_one K x z]
    exact ⟨fun hz => ⟨⟨z, hz⟩, rfl⟩, fun ⟨w, hw⟩ => hw ▸ w.2⟩

/-- **`adjoinIntegers K x` は `𝒪[K.carrier]` 上有限**——長らく残って
いた最後の前提。`IsIntegralClosure.finite`(整閉包が有限であること、
`A` が整閉Noetherianで `L/Frac(A)` が有限分離のとき)に上の
`IsIntegralClosure` を与えるだけ。分離性は標数0から自動。 -/
theorem module_finite_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Module.Finite 𝒪[K.carrier] (adjoinIntegers K x) := by
  haveI : IsScalarTower 𝒪[K.carrier] K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  haveI : IsScalarTower 𝒪[K.carrier] (adjoinIntegers K x)
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  haveI := isIntegralClosure_adjoinIntegers K x
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := isDiscreteValuationRing_carrierIntegers K
  exact IsIntegralClosure.finite 𝒪[K.carrier] K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (adjoinIntegers K x)

/-- **分岐指数 `e`**——`𝒪[K.carrier]` の極大イデアルが
`adjoinIntegers K x` の極大イデアルの何乗まで入るか。 -/
noncomputable def ramificationIndex {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] : ℕ :=
  haveI := isLocalRing_adjoinIntegers K x
  (IsLocalRing.maximalIdeal 𝒪[K.carrier]).ramificationIdx
    (IsLocalRing.maximalIdeal (adjoinIntegers K x))

/-- **慣性次数 `f`**——剰余体の拡大次数。`residueDegree`(剰余体の
**元の個数** `q^f`)とは別物なので注意。 -/
noncomputable def inertiaDegree {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] : ℕ :=
  haveI := isLocalRing_adjoinIntegers K x
  (IsLocalRing.maximalIdeal 𝒪[K.carrier]).inertiaDeg
    (IsLocalRing.maximalIdeal (adjoinIntegers K x))

/-- **分岐理論の基本等式 `e·f = [K(x):K]`**。局所体(完備離散付値体・
剰余体有限)における古典的な等式が、本プロジェクトの
`PAdicLocalField` 設定でそのまま使えるようになった。不分岐性
(`e = 1`)・完全分岐(`f = 1`)を定義し、次数から一方を他方で決める
ための土台。 -/
theorem ramificationIndex_mul_inertiaDegree {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    ramificationIndex K x * inertiaDegree K x
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  haveI := isLocalRing_adjoinIntegers K x
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI : IsScalarTower 𝒪[K.carrier] K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  haveI : IsScalarTower 𝒪[K.carrier] (adjoinIntegers K x)
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  haveI : IsFractionRing (adjoinIntegers K x)
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    ValuationRing.instIsFractionRingInteger
      (K := IntermediateField.adjoin K.carrier ({x} : Set K.closure)) Valued.v
  haveI := module_finite_adjoinIntegers K x
  exact Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing (adjoinIntegers K x) K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (IsDiscreteValuationRing.not_a_field 𝒪[K.carrier])

/-- **不分岐**——分岐指数が `1`。 -/
def IsUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Prop :=
  ramificationIndex K x = 1

/-- 不分岐なら慣性次数が拡大次数そのもの——`e·f=[L:K]` の直接の帰結。 -/
theorem inertiaDegree_eq_finrank_of_isUnramified {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (h : IsUnramifiedAdjoin K x) :
    inertiaDegree K x
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  have := ramificationIndex_mul_inertiaDegree K x
  rw [show ramificationIndex K x = 1 from h, one_mul] at this
  exact this


/-! ## 剰余体の側から不分岐性を判定する

`e·f = [L:K]` が手に入ったので、**`f` を下から押さえれば `e = 1` が
出る**。これは不分岐拡大の存在を示す際に Hensel の補題を迂回できる
可能性のある経路——「剰余体が十分大きい」ことさえ言えれば
`f = [L:K]` となり、`e·f = [L:K]` から `e = 1` が従う。

そのために本節では
* `adjoinIntegers K x` の剰余体が**有限**であること
  (`Found/ResidueFieldFinite.lean` の一般結果を拡大体に適用)、
* `inertiaDegree K x` が**剰余体の拡大次数**そのものであること、
* したがって `residueDegree K x`(剰余体の**元の個数**)が `q^f` で
  あること
を示し、最後に判定条件
`residueDegree K x = q^[L:K] ⟹ IsUnramifiedAdjoin K x` を得る。 -/

/-- `IsLocalRing (adjoinIntegers K x)` をインスタンスとして登録する。
`Ideal.inertiaDeg` や `IsLocalRing.ResidueField` は**主張の型の段階**で
これを要求するので、`theorem` のままだと利用側が毎回
インスタンス束縛を書く羽目になる。 -/
instance instIsLocalRingAdjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) : IsLocalRing (adjoinIntegers K x) :=
  isLocalRing_adjoinIntegers K x

/-- **`𝒪[K.carrier] → adjoinIntegers K x` は局所準同型**——ノルムが
`algebraMap` で保たれる(`norm_algebraMap'`)ことから、非単元は非単元へ
写る。`Ideal.LiesOver`(したがって `Ideal.inertiaDeg_algebraMap`)と
剰余体の間の `Algebra` 構造がこれで出る。 -/
instance isLocalHom_adjoinIntegersAlgebraMap {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) :
    IsLocalHom (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := by
  constructor
  intro a ha
  by_contra hcon
  rw [Valuation.Integer.not_isUnit_iff_valuation_lt_one] at hcon
  refine absurd ha ((Valuation.Integer.not_isUnit_iff_valuation_lt_one
    (v := (Valued.v : Valuation (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) NNReal))
    (x := (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a))).mpr ?_)
  show Valued.v ((algebraMap K.carrier
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (a : K.carrier)) < 1
  have e1 : Valued.v ((algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (a : K.carrier))
      = (‖(algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (a : K.carrier)‖₊ : NNReal) :=
    NNReal.eq rfl
  have e2 : Valued.v ((a : K.carrier)) = (‖(a : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
  have e3 : (‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (a : K.carrier)‖₊ : NNReal)
      = (‖(a : K.carrier)‖₊ : NNReal) := NNReal.eq (norm_algebraMap' _ _)
  rw [e1, e3, ← e2]
  exact hcon

/-- 拡大体 `K.carrier⟮x⟯` も自明でないノルム体。基礎体の
`NormedField.exists_one_lt_norm` の元を `algebraMap` で送るだけ。

★配管(記録): `@[implicit_reducible]` を**必ず付ける**。付けないと
`letI` で入れた瞬間に `IsUltrametricDist ↥K.carrier⟮x⟯` の
インスタンス探索が(ノルム構造の経路が変わって)失敗する。 -/
@[implicit_reducible] noncomputable def nontriviallyNormedField_adjoin {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure) :
    NontriviallyNormedField (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  { (inferInstance : NormedField (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) with
    non_trivial := by
      obtain ⟨y, hy⟩ := NormedField.exists_one_lt_norm K.carrier
      exact ⟨algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) y,
        lt_of_lt_of_eq hy (norm_algebraMap' _ y).symm⟩ }

/-- **拡大体の剰余体も有限**——`Found/ResidueFieldFinite.lean` の
一般結果 `finite_residueField`(完備・proper な非アルキメデスノルム体)
を `K.carrier⟮x⟯` に適用する。`ProperSpace` は有限次元性から。 -/
instance finite_residueField_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Finite (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
  letI := nontriviallyNormedField_adjoin K x
  haveI : ProperSpace (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    FiniteDimensional.proper K.carrier _
  exact finite_residueField (K := IntermediateField.adjoin K.carrier ({x} : Set K.closure))

/-- **慣性次数は剰余体の拡大次数**。`Ideal.inertiaDeg_algebraMap` に
`Ideal.LiesOver`(局所準同型から自動)を与えるだけ。 -/
theorem inertiaDegree_eq_finrank_residueField {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    inertiaDegree K x = Module.finrank (IsLocalRing.ResidueField 𝒪[K.carrier])
      (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
  unfold inertiaDegree
  exact Ideal.inertiaDeg_algebraMap _ _

/-- **`q_L = q^f`**——拡大体の剰余体の元の個数は、基礎体のそれの
慣性次数乗。有限体上の有限次元ベクトル空間の元の個数
(`Module.card_eq_pow_finrank`)そのもの。 -/
theorem residueDegree_eq_residueCard_pow {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    residueDegree K x
      = Nat.card (IsLocalRing.ResidueField 𝒪[K.carrier]) ^ inertiaDegree K x := by
  haveI : Fintype (IsLocalRing.ResidueField 𝒪[K.carrier]) := Fintype.ofFinite _
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
  unfold residueDegree
  rw [inertiaDegree_eq_finrank_residueField K x, Nat.card_eq_fintype_card,
    Nat.card_eq_fintype_card, Module.card_eq_pow_finrank
      (K := IsLocalRing.ResidueField 𝒪[K.carrier])
      (V := IsLocalRing.ResidueField (adjoinIntegers K x))]

/-- **不分岐性の判定(慣性次数版)**——`f = [L:K]` なら `e = 1`。
`e·f = [L:K]` と `[L:K] > 0` から。 -/
theorem isUnramifiedAdjoin_of_inertiaDegree_eq_finrank {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (h : inertiaDegree K x
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    IsUnramifiedAdjoin K x := by
  have hmul := ramificationIndex_mul_inertiaDegree K x
  rw [h] at hmul
  have hpos : 0 < Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
  unfold IsUnramifiedAdjoin
  nlinarith [hmul, hpos]

/-- 剰余体は有限体なので元の個数は `2` 以上。指数の一意性
(`Nat.pow_right_injective`)に要る。 -/
theorem one_lt_residueCard {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    1 < Nat.card (IsLocalRing.ResidueField 𝒪[K.carrier]) := by
  haveI : Finite (IsLocalRing.ResidueField 𝒪[K.carrier]) := residueField_finite K
  exact Finite.one_lt_card

/-- **不分岐性の判定(剰余体の元の個数版)**——`q_L = q^{[L:K]}` なら
不分岐。**Hensel の補題を使わずに不分岐拡大を作る**ときの目標形:
剰余体が `q^n` 個の元を持つことさえ言えれば、`e=1` が自動的に出る。 -/
theorem isUnramifiedAdjoin_of_residueDegree {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (h : residueDegree K x = Nat.card (IsLocalRing.ResidueField 𝒪[K.carrier])
      ^ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    IsUnramifiedAdjoin K x := by
  apply isUnramifiedAdjoin_of_inertiaDegree_eq_finrank
  have h2 := residueDegree_eq_residueCard_pow K x
  rw [h] at h2
  exact Nat.pow_right_injective (one_lt_residueCard K) h2.symm


/-! ## ★訂正: Hensel の補題は mathlib に**ある**——`𝒪[K.carrier]` は Henselian

`memory/pgc-unramified-extension-progress.md` に「完備局所環が
Henselian であることに相当するインスタンスは mathlib に無い」と
記録していたが、**これは誤りだった**。探し方が悪く、
`HenselianLocalRing` を結論とする宣言(`Field.henselian` しか無い)
だけを見ていた。正しい入口は **`HenselianRing R I`** の方で、

```
IsAdicComplete.henselianRing (R) (I) [IsAdicComplete I R] : HenselianRing R I
```

が存在する(`Mathlib/RingTheory/Henselian.lean`)。そして本プロジェクトは
すでに `ABC3.Found.PGC.isAdicComplete_valuationRing`
(`Found/PGC/ValuationRingComplete.lean`)を持っている——
`𝒪[K.carrier]` は `maximalIdeal`-進完備。

したがって `HenselianLocalRing 𝒪[K.carrier]` は**数行**で出る。
`HenselianRing` の仮定が「`f'(a₀)` の剰余体での像が単元」なのに対し
`HenselianLocalRing` は「`f'(a₀)` 自身が単元」なので、
`IsUnit.map` を一つ挟むだけ。

★測定の教訓(`Found/ResidueFieldFinite.lean` の docstring と同じ轍):
「無い」という測定は**探索範囲を書かないと再現できない**。今回は
`HenselianLocalRing` だけを引いて `HenselianRing` を引かなかった。 -/

/-- **`𝒪[K.carrier]` は Henselian 局所環**——`isAdicComplete_
valuationRing`(既存)に `IsAdicComplete.henselianRing` を当てるだけ。
剰余体での単根が `𝒪[K.carrier]` へ持ち上がる。 -/
instance henselianLocalRing_carrierIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    HenselianLocalRing 𝒪[K.carrier] := by
  haveI := isAdicComplete_valuationRing K
  haveI := IsAdicComplete.henselianRing 𝒪[K.carrier] (IsLocalRing.maximalIdeal 𝒪[K.carrier])
  constructor
  intro f hf a₀ h0 hu
  exact HenselianRing.is_henselian f hf a₀ h0 (hu.map (Ideal.Quotient.mk _))


/-! ## ★不分岐拡大の構成——剰余体の既約多項式を持ち上げる

`isUnramifiedAdjoin_of_inertiaDegree_eq_finrank`(`f=[L:K] ⟹ e=1`)が
あるので、**Hensel の補題を使わずに**不分岐拡大が作れる:

1. 剰余体 `𝓀[K.carrier]` 上のモニック既約多項式 `g` を取る。
2. `Polynomial.lifts_and_degree_eq_and_monic` で `𝒪[K.carrier]` 上の
   モニック `f`(同次数、`f ↦ g`)へ持ち上げる。
3. `Polynomial.Monic.irreducible_of_irreducible_map` で `f` も既約。
4. Gauss(`Monic.irreducible_iff_irreducible_map_fraction_map`、
   `IsIntegrallyClosed 𝒪[K.carrier]` は付値環だから自動)で
   `K.carrier` 上でも既約。
5. `K.closure` は代数閉体なので根 `x` が取れ、`minpoly K.carrier x`
   がその既約多項式そのものだから `[K(x):K] = deg g`。
6. `f` がモニックだから `x` は `𝒪[K.carrier]` 上整、よって `‖x‖ ≤ 1`
   (`norm_le_one_of_isIntegral`)——`x` は `adjoinIntegers K x` の元。
7. 剰余をとると `x̄` は `g` の根。`g` は既約モニックだから
   `minpoly 𝓀 x̄ = g`、よって `deg g ≤ f = inertiaDegree K x`。
8. `e·f = [K(x):K] = deg g` と `e ≥ 1` から `f ≤ deg g`。
   両方合わせて `f = deg g = [K(x):K]`、すなわち **`e = 1`**。

残るのは「有限体上の任意次数のモニック既約多項式の存在」だけ
(これは純粋に有限体の話で、この節とは独立に足せる)。 -/

/-- `K.closure` は `K.carrier` の代数閉包なので、単項生成の中間体は
つねに有限次元。インスタンスにしておくと以降の
`[FiniteDimensional ...]` 束縛が自動で埋まる。 -/
instance finiteDimensional_adjoin_closure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) :
    FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  IntermediateField.adjoin.finiteDimensional
    (IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic x))

/-- **★不分岐拡大の構成**——剰余体上のモニック既約多項式 `g` から、
次数 `deg g` の**不分岐**単項拡大 `K(x)/K` を作る。Hensel の補題は
使わない(上の見取り図を参照)。 -/
theorem exists_isUnramifiedAdjoin_of_irreducible {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (g : Polynomial 𝓀[K.carrier]) (hgm : g.Monic) (hgi : Irreducible g) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        = g.natDegree ∧ IsUnramifiedAdjoin K x := by
  have hlifts : g ∈ Polynomial.lifts (IsLocalRing.residue 𝒪[K.carrier]) :=
    Polynomial.lifts_iff_coeff_lifts g |>.mpr (fun n => ⟨_, Quotient.out_eq' _⟩)
  obtain ⟨f, hfmap, hfdeg, hfm⟩ := Polynomial.lifts_and_degree_eq_and_monic hlifts hgm
  have hgnd : g.natDegree = f.natDegree := by rw [← hfmap, hfm.natDegree_map]
  have hfi : Irreducible f := hfm.irreducible_of_irreducible_map _ f (by rw [hfmap]; exact hgi)
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  have hfKi : Irreducible (f.map (algebraMap 𝒪[K.carrier] K.carrier)) :=
    (hfm.irreducible_iff_irreducible_map_fraction_map (K := K.carrier)).mp hfi
  have hFm : (f.map (algebraMap 𝒪[K.carrier] K.carrier)).Monic := hfm.map _
  have hFCm : ((f.map (algebraMap 𝒪[K.carrier] K.carrier)).map
      (algebraMap K.carrier K.closure)).Monic := hFm.map _
  have hndpos : 0 < f.natDegree := hgnd ▸ Irreducible.natDegree_pos hgi
  obtain ⟨x, hx⟩ := IsAlgClosed.exists_root
    ((f.map (algebraMap 𝒪[K.carrier] K.carrier)).map (algebraMap K.carrier K.closure))
    (by
      rw [Polynomial.degree_eq_natDegree hFCm.ne_zero, hFm.natDegree_map, hfm.natDegree_map]
      exact_mod_cast Nat.pos_iff_ne_zero.mp hndpos)
  have haev : Polynomial.aeval x (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0 := by
    rw [Polynomial.aeval_def, ← Polynomial.eval_map]; exact hx
  have hint : IsIntegral K.carrier x := IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic x)
  have hmp : f.map (algebraMap 𝒪[K.carrier] K.carrier) = minpoly K.carrier x :=
    minpoly.eq_of_irreducible_of_monic hfKi haev hFm
  have hrank : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = g.natDegree := by
    rw [IntermediateField.adjoin.finrank hint, ← hmp, hfm.natDegree_map, hgnd]
  refine ⟨x, hrank, ?_⟩
  have hval : (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val
      ⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ = x := rfl
  have hkey := Polynomial.aeval_algHom_apply
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val
    (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (f.map (algebraMap 𝒪[K.carrier] K.carrier))
  rw [hval, haev] at hkey
  have h1 : Polynomial.aeval
      (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0 :=
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val.injective
      (by rw [map_zero]; exact hkey.symm)
  have halg : (algebraMap 𝒪[K.carrier]
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = (algebraMap K.carrier
          (IntermediateField.adjoin K.carrier ({x} : Set K.closure))).comp
        (Subring.subtype 𝒪[K.carrier]) := rfl
  have h3 : Polynomial.eval₂ (algebraMap 𝒪[K.carrier]
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) f = 0 := by
    rw [halg, ← Polynomial.eval₂_map]; exact h1
  have hintL : IsIntegral 𝒪[K.carrier]
      (⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨f, hfm, h3⟩
  have hnorm : ‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1 :=
    norm_le_one_of_isIntegral K x _ hintL
  have hsub : ((Subring.subtype (adjoinIntegers K x)).comp
      (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)))
      = algebraMap 𝒪[K.carrier]
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := rfl
  have h2 := Polynomial.hom_eval₂ f (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
    (Subring.subtype (adjoinIntegers K x))
    (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩ : adjoinIntegers K x)
  rw [hsub] at h2
  have hrootO : Polynomial.eval₂ (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩ : adjoinIntegers K x)
      f = 0 := Subtype.ext (h2.trans h3)
  have hbar : Polynomial.aeval (IsLocalRing.residue (adjoinIntegers K x)
      ⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩)
      (f.map (IsLocalRing.residue 𝒪[K.carrier])) = 0 := by
    rw [Polynomial.aeval_def, Polynomial.eval₂_map]
    have hcomp : (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))).comp
        (IsLocalRing.residue 𝒪[K.carrier])
        = (IsLocalRing.residue (adjoinIntegers K x)).comp
          (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := RingHom.ext (congrFun rfl)
    rw [hcomp, ← Polynomial.hom_eval₂ f (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (IsLocalRing.residue (adjoinIntegers K x)) _, hrootO, map_zero]
  rw [hfmap] at hbar
  haveI := module_finite_adjoinIntegers K x
  have hge : g.natDegree ≤ inertiaDegree K x := by
    have hmpbar : g = minpoly 𝓀[K.carrier] (IsLocalRing.residue (adjoinIntegers K x)
        ⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩) :=
      minpoly.eq_of_irreducible_of_monic hgi hbar hgm
    rw [inertiaDegree_eq_finrank_residueField K x, hmpbar]
    exact minpoly.natDegree_le _
  apply isUnramifiedAdjoin_of_inertiaDegree_eq_finrank
  rw [hrank]
  have hmul := ramificationIndex_mul_inertiaDegree K x
  rw [hrank] at hmul
  have hnpos : 0 < g.natDegree := Irreducible.natDegree_pos hgi
  have he : ramificationIndex K x ≠ 0 := by
    intro h; rw [h, zero_mul] at hmul; omega
  have hle : inertiaDegree K x ≤ ramificationIndex K x * inertiaDegree K x :=
    Nat.le_mul_of_pos_left _ (Nat.pos_of_ne_zero he)
  omega

end ABC3.Found.PGC
