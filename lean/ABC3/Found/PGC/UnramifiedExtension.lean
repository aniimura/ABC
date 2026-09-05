import ABC3.Found.PGC.AdjoinIntegers
import Mathlib.RingTheory.DedekindDomain.IntegralClosure
import Mathlib.NumberTheory.RamificationInertia.Basic
import Mathlib.RingTheory.Henselian
import Mathlib.Algebra.Polynomial.Eval.Irreducible
import Mathlib.RingTheory.Polynomial.GaussLemma
import ABC3.Found.FiniteFieldIrreducible
import ABC3.Found.HenselianSplits
import ABC3.Found.Teichmuller
import ABC3.Found.PGC.ValuationRingComplete

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


/-- `K.closure` は `K.carrier` の代数閉包なので、単項生成の中間体は
つねに有限次元。インスタンスにしておくと以降の
`[FiniteDimensional ...]` 束縛が自動で埋まる。 -/
instance finiteDimensional_adjoin_closure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) :
    FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
  IntermediateField.adjoin.finiteDimensional
    (IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic x))

/-! ## 拡大体側も Henselian——一意性へ向けた準備

不分岐拡大の**一意性**(同じ次数なら一致)を示すには、
「`𝒪_{K(y)}` の剰余体 `𝔽_{q^n}` にある根を `𝒪_{K(y)}` へ持ち上げる」
段で `HenselianLocalRing (adjoinIntegers K y)` が要る。基礎体側
(`henselianLocalRing_carrierIntegers`)と同じ筋——`IsAdic`(進位相=
距離位相)+ コンパクト性からの完備性 → `IsAdicComplete` →
`IsAdicComplete.henselianRing`——が拡大体でもそのまま通る。

★配管(記録): `Found/PGC/ValuationRingComplete.lean` の基礎体版は
`rw [show maximalIdeal ... = Valued.maximalIdeal ... from rfl, ...]` で
書けたが、拡大体では `↥(adjoinIntegers K x)` と `↥𝒪[K.carrier⟮x⟯]` の
`CommRing` インスタンス経路が違う(`CommRing.toCommSemiring` と
`SubsemiringClass.toCommSemiring`)ため `rw` が型検査で落ちる。
**`▸` の項レベルのキャスト**にすると defeq で通る
(`tools/lean-idioms.md` #37 の類型)。また `Valued.integer.norm_
irreducible_pos` 等は `[NontriviallyNormedField K]` を要求するので、
`letI := nontriviallyNormedField_adjoin K x` を先に置かないと
インスタンス探索が `whnf` でタイムアウトする。 -/

/-- 拡大体の整数環でも「`maximalIdeal`-進位相 = 距離位相」。 -/
theorem isAdic_maximalIdeal_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) :
    IsAdic (IsLocalRing.maximalIdeal (adjoinIntegers K x)) := by
  letI := nontriviallyNormedField_adjoin K x
  haveI : IsDiscreteValuationRing 𝒪[IntermediateField.adjoin K.carrier ({x} : Set K.closure)] :=
    isDiscreteValuationRing_adjoinIntegers K x
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible
    𝒪[IntermediateField.adjoin K.carrier ({x} : Set K.closure)]
  have hϖnorm0 : ‖ϖ‖ ≠ 0 := (Valued.integer.norm_irreducible_pos hϖ).ne'
  have hϖlt1 : ‖ϖ‖ < 1 := Valued.integer.norm_irreducible_lt_one hϖ
  have hball := Irreducible.maximalIdeal_pow_eq_closedBall_pow hϖ
  rw [isAdic_iff]
  refine ⟨fun n => ?_, fun s hs => ?_⟩
  · have h1 : IsOpen (Metric.closedBall
        (0 : 𝒪[IntermediateField.adjoin K.carrier ({x} : Set K.closure)]) (‖ϖ‖ ^ n)) :=
      IsUltrametricDist.isOpen_closedBall 0 (pow_ne_zero n hϖnorm0)
    exact (hball n) ▸ h1
  · rw [Metric.mem_nhds_iff] at hs
    obtain ⟨ε, hεpos, hεsub⟩ := hs
    obtain ⟨n, hn⟩ := exists_pow_lt_of_lt_one hεpos hϖlt1
    refine ⟨n, ?_⟩
    have h2 : Metric.closedBall
        (0 : 𝒪[IntermediateField.adjoin K.carrier ({x} : Set K.closure)]) (‖ϖ‖ ^ n) ⊆ s :=
      (Metric.closedBall_subset_ball hn).trans hεsub
    exact (hball n) ▸ h2

/-- 拡大体の整数環も `maximalIdeal`-進完備。 -/
theorem isAdicComplete_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) :
    IsAdicComplete (IsLocalRing.maximalIdeal (adjoinIntegers K x)) (adjoinIntegers K x) := by
  haveI : CompactSpace (adjoinIntegers K x) := compactSpace_adjoinIntegers K x
  exact (isAdic_maximalIdeal_adjoinIntegers K x).isAdicComplete_iff.mpr
    ⟨inferInstance, inferInstance⟩

/-- **拡大体の整数環も Henselian 局所環**。剰余体での単根が
`adjoinIntegers K x` へ持ち上がる——不分岐拡大の一意性の要。 -/
instance henselianLocalRing_adjoinIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) : HenselianLocalRing (adjoinIntegers K x) := by
  haveI := isAdicComplete_adjoinIntegers K x
  haveI := IsAdicComplete.henselianRing (adjoinIntegers K x)
    (IsLocalRing.maximalIdeal (adjoinIntegers K x))
  constructor
  intro f hf a₀ h0 hu
  exact HenselianRing.is_henselian f hf a₀ h0 (hu.map (Ideal.Quotient.mk _))


/-! ## Galois 群から剰余体の Galois 群への還元射

不分岐拡大の理論の核心は
**`Gal(K(x)/K) ≅ Gal(𝓀_{K(x)}/𝓀)`**(そして右辺は Frobenius が生成する
巡回群)——これが `Gal(K^ur/K) ≅ Ẑ` の出どころ。本節ではその射
`residueGalHom` を構成する。

構成は 3 段:
1. `σ : K(x) ≃ₐ[K] K(x)` は**ノルムを保つ**(`norm_algEquiv`)——
   スペクトルノルムが Galois 共役で不変であること
   (`spectralNorm_eq_of_equiv`)と `NormedAlgebra.norm_eq_spectralNorm`
   から。したがって `σ` は整数環 `adjoinIntegers K x` を保つ。
2. よって `σ` は `adjoinIntegers K x` の環同型 `algEquivIntegers` を誘導し、
3. さらに剰余体の `𝓀`-代数同型 `residueAlgEquiv` を誘導する
   (`σ` が `K.carrier` を固定するので `𝓀` を固定する)。

これを群準同型にまとめたものが `residueGalHom`。 -/

/-- `K(x)` の `K`-自己同型はノルムを保つ。 -/
theorem norm_algEquiv {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖σ z‖ = ‖z‖ := by
  rw [NormedAlgebra.norm_eq_spectralNorm K.carrier (σ z),
    NormedAlgebra.norm_eq_spectralNorm K.carrier z]
  exact (spectralNorm_eq_of_equiv σ z).symm

/-- ノルムを保つので、`σ` は整数環 `adjoinIntegers K x` の環同型を誘導する。 -/
noncomputable def algEquivIntegers {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    adjoinIntegers K x ≃+* adjoinIntegers K x where
  toFun z := ⟨σ (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
    by show ‖σ (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1
       rw [norm_algEquiv K x σ]; exact z.2⟩
  invFun z := ⟨σ.symm (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
    by show ‖σ.symm (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1
       rw [norm_algEquiv K x σ.symm]; exact z.2⟩
  left_inv z := by apply Subtype.ext; simp
  right_inv z := by apply Subtype.ext; simp
  map_mul' a b := by apply Subtype.ext; simp
  map_add' a b := by apply Subtype.ext; simp

/-- さらに剰余体の `𝓀[K.carrier]`-代数同型を誘導する。 -/
noncomputable def residueAlgEquiv {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :
    IsLocalRing.ResidueField (adjoinIntegers K x)
      ≃ₐ[𝓀[K.carrier]] IsLocalRing.ResidueField (adjoinIntegers K x) :=
  { IsLocalRing.ResidueField.mapEquiv (algEquivIntegers K x σ) with
    commutes' := by
      intro r
      obtain ⟨a, rfl⟩ := IsLocalRing.residue_surjective r
      have h1 : algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))
          (IsLocalRing.residue 𝒪[K.carrier] a)
          = IsLocalRing.residue (adjoinIntegers K x)
            (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) a) := rfl
      show (IsLocalRing.ResidueField.mapEquiv (algEquivIntegers K x σ)) _ = _
      rw [h1, IsLocalRing.ResidueField.mapEquiv_apply, IsLocalRing.ResidueField.map_residue]
      congr 1
      apply Subtype.ext
      show σ (algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (a : K.carrier)) = _
      rw [AlgEquiv.commutes]
      rfl }

/-- `residueAlgEquiv` の剰余元での値。 -/
theorem residueAlgEquiv_apply {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (b : adjoinIntegers K x) :
    residueAlgEquiv K x σ (IsLocalRing.residue (adjoinIntegers K x) b)
      = IsLocalRing.residue (adjoinIntegers K x) (algEquivIntegers K x σ b) := by
  show (IsLocalRing.ResidueField.mapEquiv (algEquivIntegers K x σ)) _ = _
  rw [IsLocalRing.ResidueField.mapEquiv_apply, IsLocalRing.ResidueField.map_residue]
  rfl

/-- **還元射 `Gal(K(x)/K) →* Gal(𝓀_{K(x)}/𝓀)`**。 -/
noncomputable def residueGalHom {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure) :
    ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    →* (IsLocalRing.ResidueField (adjoinIntegers K x)
      ≃ₐ[𝓀[K.carrier]] IsLocalRing.ResidueField (adjoinIntegers K x)) where
  toFun := residueAlgEquiv K x
  map_one' := by
    apply AlgEquiv.ext
    intro z
    obtain ⟨b, rfl⟩ := IsLocalRing.residue_surjective z
    rw [residueAlgEquiv_apply, AlgEquiv.one_apply]
    rfl
  map_mul' σ τ := by
    apply AlgEquiv.ext
    intro z
    obtain ⟨b, rfl⟩ := IsLocalRing.residue_surjective z
    rw [residueAlgEquiv_apply, AlgEquiv.mul_apply, residueAlgEquiv_apply, residueAlgEquiv_apply]
    rfl


/-! ### 還元射の単射性と全単射性

単射性の要は **Hensel の一意性**(`IsLocalRing.eq_of_eval_eq_zero_of_
not_isUnit_sub`): `σ` が剰余体で恒等なら `σ(x)` と `x` は同じ剰余を持つ
`f` の 2 根、`f'(x)` が単元(= `f̄` が分離的で `x̄` が単根)なので
`σ(x) = x`、`x` は生成元だから `σ = 1`。

全単射性は位数の比較——`|Gal(K(x)/K)| = [K(x):K]`(Galois)と
`|Gal(𝓀_{K(x)}/𝓀)| = f`(有限体は Galois)、不分岐なので両者は等しい。 -/

/-- 生成元 `x` を固定する `K`-自己同型は恒等。 -/
theorem algEquiv_eq_one_of_fixes_gen {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (h : σ ⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩
      = ⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩) : σ = 1 := by
  have hext : (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) →ₐ[K.carrier]
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = AlgHom.id K.carrier _ := by
    refine IntermediateField.algHom_ext_of_eq_adjoin K.carrier rfl ?_
    intro y hy
    simp only [Set.mem_singleton_iff] at hy
    subst hy
    exact h
  apply AlgEquiv.ext
  intro z
  have := AlgHom.congr_fun hext z
  simpa using this

/-- `f̄` が既約(有限体上なので分離的)なら、その根で `f'` は単元。 -/
theorem isUnit_eval_derivative {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (f : Polynomial 𝒪[K.carrier])
    (hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])))
    (b : adjoinIntegers K x)
    (hbar : Polynomial.aeval (IsLocalRing.residue (adjoinIntegers K x) b)
      (f.map (IsLocalRing.residue 𝒪[K.carrier])) = 0) :
    IsUnit (Polynomial.eval b (Polynomial.derivative
      (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))))) := by
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  have hcomp : (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))).comp
      (IsLocalRing.residue 𝒪[K.carrier])
      = (IsLocalRing.residue (adjoinIntegers K x)).comp
        (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := RingHom.ext (congrFun rfl)
  have hFAres : (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))
      = (f.map (IsLocalRing.residue 𝒪[K.carrier])).map
        (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))) := by
    rw [Polynomial.map_map, Polynomial.map_map, ← hcomp]
  refine (IsLocalRing.residue_ne_zero_iff_isUnit _).mp ?_
  have hev := Polynomial.hom_eval₂ (Polynomial.derivative
      (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)))) (RingHom.id _)
    (IsLocalRing.residue (adjoinIntegers K x)) b
  simp only [Polynomial.eval₂_id, RingHom.comp_id] at hev
  rw [hev, ← Polynomial.eval_map, ← Polynomial.derivative_map, hFAres]
  have hsep : ((f.map (IsLocalRing.residue 𝒪[K.carrier])).map
      (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)))).Separable :=
    (PerfectField.separable_of_irreducible hgi).map
  have hroot : Polynomial.eval (IsLocalRing.residue (adjoinIntegers K x) b)
      ((f.map (IsLocalRing.residue 𝒪[K.carrier])).map
        (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)))) = 0 := by
    rw [Polynomial.eval_map, ← Polynomial.aeval_def]
    exact hbar
  obtain ⟨u, v, huv⟩ := hsep
  intro hcon
  have hev2 := congrArg (Polynomial.eval (IsLocalRing.residue (adjoinIntegers K x) b)) huv
  simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_one, hcon, hroot,
    mul_zero, zero_add] at hev2
  exact one_ne_zero hev2.symm

/-- **還元射は単射**——Hensel の一意性から。 -/
theorem residueGalHom_injective {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (f : Polynomial 𝒪[K.carrier])
    (hnorm : ‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1)
    (hroot : Polynomial.eval₂ (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩ : adjoinIntegers K x)
      f = 0)
    (hunit : IsUnit (Polynomial.eval
      (⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩ : adjoinIntegers K x)
      (Polynomial.derivative
        (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)))))) :
    Function.Injective (residueGalHom K x) := by
  rw [injective_iff_map_eq_one]
  intro σ hσ
  set xO : adjoinIntegers K x :=
    ⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩ with hxO
  have hres : IsLocalRing.residue (adjoinIntegers K x) (algEquivIntegers K x σ xO)
      = IsLocalRing.residue (adjoinIntegers K x) xO := by
    have h1 := congrArg (fun (e : IsLocalRing.ResidueField (adjoinIntegers K x)
      ≃ₐ[𝓀[K.carrier]] IsLocalRing.ResidueField (adjoinIntegers K x)) =>
        e (IsLocalRing.residue (adjoinIntegers K x) xO)) hσ
    have h0 : residueGalHom K x σ = residueAlgEquiv K x σ := rfl
    simp only [h0] at h1
    rw [residueAlgEquiv_apply] at h1
    simpa using h1
  have hgψ : ((algEquivIntegers K x σ : adjoinIntegers K x →+* adjoinIntegers K x).comp
      (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)))
      = algebraMap 𝒪[K.carrier] (adjoinIntegers K x) := by
    refine RingHom.ext fun a => ?_
    apply Subtype.ext
    show σ (algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (a : K.carrier)) = _
    rw [AlgEquiv.commutes]
    rfl
  have hb : Polynomial.eval xO (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))) = 0 := by
    rw [Polynomial.eval_map]; exact hroot
  have ha : Polynomial.eval (algEquivIntegers K x σ xO)
      (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))) = 0 := by
    have h2 := Polynomial.hom_eval₂ f (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (algEquivIntegers K x σ : adjoinIntegers K x →+* adjoinIntegers K x) xO
    rw [hgψ, hroot, map_zero] at h2
    rw [Polynomial.eval_map]
    exact h2.symm
  have hnu : ¬ IsUnit (xO - algEquivIntegers K x σ xO) := by
    rw [← mem_nonunits_iff, ← IsLocalRing.mem_maximalIdeal,
      ← IsLocalRing.residue_eq_zero_iff, map_sub, hres, sub_self]
  have hkey : xO = algEquivIntegers K x σ xO :=
    IsLocalRing.eq_of_eval_eq_zero_of_not_isUnit_sub hb ha hnu hunit
  refine algEquiv_eq_one_of_fixes_gen K x σ ?_
  have h3 := congrArg (fun (z : adjoinIntegers K x) =>
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) hkey
  exact h3.symm

/-- **還元射は全単射**——単射性と位数の一致(不分岐なので `[K(x):K] = f`)から。 -/
theorem residueGalHom_bijective {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (hnor : Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (hu : IsUnramifiedAdjoin K x)
    (hinj : Function.Injective (residueGalHom K x)) :
    Function.Bijective (residueGalHom K x) := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := hnor
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : IsGalois 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
    infer_instance
  rw [Fintype.bijective_iff_injective_and_card]
  refine ⟨hinj, ?_⟩
  rw [← Nat.card_eq_fintype_card, ← Nat.card_eq_fintype_card,
    IsGalois.card_aut_eq_finrank, IsGalois.card_aut_eq_finrank,
    ← inertiaDegree_eq_finrank_residueField K x,
    inertiaDegree_eq_finrank_of_isUnramified K x hu]


/-! ## Hensel による分裂の持ち上げ——不分岐拡大は normal

`Found/HenselianSplits.lean`(一般の結果)を本設定に流し込む。剰余体
`𝓀_{K(x)}` は `𝓀` の**有限体の有限次拡大**なので normal——したがって
`x̄` の最小多項式 `ḡ` は `𝓀_{K(x)}` で分裂する。分離性は有限体が完全体
であることから。これを Hensel で `𝒪_{K(x)}` へ持ち上げると、定義多項式
`f` が `K(x)` で分裂する。あとは `K(x)` が `f` の分裂体であることを見て
`Normal.of_isSplittingField` で結論する。

標数 0 なので `Normal` から直ちに `IsGalois` が従う——不分岐拡大の
Galois 性・一意性・`Gal(K^ur/K) ≅ Ẑ` へ向かう足場。 -/

/-- 剰余体の既約多項式 `ḡ` の持ち上げ `f` は、`K(x)` の中で分裂する。 -/
theorem splits_adjoin_of_lift {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (f : Polynomial 𝒪[K.carrier]) (hfm : f.Monic)
    (hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])))
    (xbar : IsLocalRing.ResidueField (adjoinIntegers K x))
    (hmpbar : f.map (IsLocalRing.residue 𝒪[K.carrier]) = minpoly 𝓀[K.carrier] xbar) :
    ((f.map (algebraMap 𝒪[K.carrier] K.carrier)).map (algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))).Splits := by
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  have hcomp : (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))).comp
      (IsLocalRing.residue 𝒪[K.carrier])
      = (IsLocalRing.residue (adjoinIntegers K x)).comp
        (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := RingHom.ext (congrFun rfl)
  have hFAres : (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))
      = (f.map (IsLocalRing.residue 𝒪[K.carrier])).map
        (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))) := by
    rw [Polynomial.map_map, Polynomial.map_map, ← hcomp]
  have hsep : ((f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))).Separable := by
    rw [hFAres]; exact (PerfectField.separable_of_irreducible hgi).map
  have hspl : ((f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))).Splits := by
    rw [hFAres, hmpbar]; exact Normal.splits (by infer_instance) xbar
  have hmain := ABC3.Found.splits_map_of_residue_splits
    (Subring.subtype (adjoinIntegers K x)) Subtype.val_injective
    (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))) (hfm.map _) hsep hspl
  rw [Polynomial.map_map] at hmain
  rw [Polynomial.map_map]
  exact hmain

/-- `K(x)` の中で `F` が分裂し `x` が `F` の根なら、`K(x)/K` は normal
——`K(x)` が `F` の分裂体そのものになる(`K.closure` の中で `F` の根は
すべて `K(x)` に収まる)。 -/
theorem normal_of_splits_in_adjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (F : Polynomial K.carrier) (hFne : F ≠ 0)
    (hFsplitsL : (F.map (algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))).Splits)
    (hxroot : Polynomial.aeval x F = 0) :
    Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  haveI := IsAlgClosure.normal K.carrier K.closure
  have hsplitsC : (F.map (algebraMap K.carrier K.closure)).Splits := IsAlgClosed.splits _
  have hsplitfield := IntermediateField.adjoin_rootSet_isSplittingField hsplitsC
  have htower : (algebraMap K.carrier K.closure)
      = (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
        (algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    IsScalarTower.algebraMap_eq _ _ _
  have heq : IntermediateField.adjoin K.carrier (F.rootSet K.closure)
      = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    apply le_antisymm
    · rw [IntermediateField.adjoin_le_iff]
      intro z hz
      have hz' : z ∈ (F.map (algebraMap K.carrier K.closure)).roots := by
        rw [Polynomial.mem_roots (by
          simpa [Polynomial.map_eq_zero_iff (algebraMap K.carrier K.closure).injective] using hFne)]
        obtain ⟨h1, h2⟩ := Polynomial.mem_rootSet'.mp hz
        rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
        exact h2
      rw [htower, ← Polynomial.map_map, Polynomial.Splits.roots_map hFsplitsL] at hz'
      obtain ⟨r, _, hr⟩ := Multiset.mem_map.mp hz'
      exact hr ▸ r.2
    · rw [IntermediateField.adjoin_le_iff]
      intro z hz
      simp only [Set.mem_singleton_iff] at hz
      subst hz
      refine IntermediateField.subset_adjoin K.carrier (F.rootSet K.closure) ?_
      rw [Polynomial.mem_rootSet']
      exact ⟨by simpa [Polynomial.map_eq_zero_iff (algebraMap K.carrier K.closure).injective]
        using hFne, hxroot⟩
  rw [heq] at hsplitfield
  exact Normal.of_isSplittingField F

/-- `F` が `K(x)` で分裂するなら、`K.closure` にある `F` の根はすべて
`K(x)` に入る(`normal_of_splits_in_adjoin` の「⊆」の部分を切り出したもの)。 -/
theorem mem_adjoin_of_root_of_splits {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) (F : Polynomial K.carrier) (hFne : F ≠ 0)
    (hsplits : (F.map (algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))).Splits)
    {z : K.closure} (hz : Polynomial.aeval z F = 0) :
    z ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  have htower : (algebraMap K.carrier K.closure)
      = (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
        (algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    IsScalarTower.algebraMap_eq _ _ _
  have hz' : z ∈ (F.map (algebraMap K.carrier K.closure)).roots := by
    rw [Polynomial.mem_roots (by
      simpa [Polynomial.map_eq_zero_iff (algebraMap K.carrier K.closure).injective] using hFne)]
    rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
    exact hz
  rw [htower, ← Polynomial.map_map, Polynomial.Splits.roots_map hsplits] at hz'
  obtain ⟨r, _, hr⟩ := Multiset.mem_map.mp hz'
  exact hr ▸ r.2

/-! ## ★不分岐拡大の判定子——「剰余体で既約になるモニック多項式の根」

以下が本節の中心定理。`𝒪[K.carrier]` 上のモニック多項式 `f` で、
その**剰余が既約**なものの根 `x` を取ると、`K(x)/K` は自動的に

* 次数 `deg f`、
* **不分岐**、
* **normal**(標数 0 なので Galois)、
* `Gal(K(x)/K) ≃ Gal(𝓀_{K(x)}/𝓀)`

を満たす。存在(各次数)も一意性も、この判定子に流し込むだけで出る。 -/

/-- **★不分岐拡大の判定子**。 -/
theorem isUnramifiedAdjoin_of_lift {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) (f : Polynomial 𝒪[K.carrier]) (hfm : f.Monic)
    (hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])))
    (haev : Polynomial.aeval x (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0) :
    Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        = f.natDegree
      ∧ IsUnramifiedAdjoin K x
      ∧ Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ∧ Function.Bijective (residueGalHom K x) := by
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  have hfKi : Irreducible (f.map (algebraMap 𝒪[K.carrier] K.carrier)) :=
    (hfm.irreducible_iff_irreducible_map_fraction_map (K := K.carrier)).mp
      (hfm.irreducible_of_irreducible_map _ f hgi)
  have hFm : (f.map (algebraMap 𝒪[K.carrier] K.carrier)).Monic := hfm.map _
  have hint : IsIntegral K.carrier x := IsAlgebraic.isIntegral (Algebra.IsAlgebraic.isAlgebraic x)
  have hmp : f.map (algebraMap 𝒪[K.carrier] K.carrier) = minpoly K.carrier x :=
    minpoly.eq_of_irreducible_of_monic hfKi haev hFm
  have hrank : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = f.natDegree := by
    rw [IntermediateField.adjoin.finrank hint, ← hmp, hfm.natDegree_map]
  have hnpos : 0 < f.natDegree := by
    have := Irreducible.natDegree_pos hgi
    rwa [hfm.natDegree_map] at this
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
  have hmpbar : f.map (IsLocalRing.residue 𝒪[K.carrier])
      = minpoly 𝓀[K.carrier] (IsLocalRing.residue (adjoinIntegers K x)
        ⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, hnorm⟩) :=
    minpoly.eq_of_irreducible_of_monic hgi hbar (hfm.map _)
  have hunram : IsUnramifiedAdjoin K x := by
    haveI := module_finite_adjoinIntegers K x
    have hge : f.natDegree ≤ inertiaDegree K x := by
      rw [inertiaDegree_eq_finrank_residueField K x, ← hfm.natDegree_map
        (IsLocalRing.residue 𝒪[K.carrier]), hmpbar]
      exact minpoly.natDegree_le _
    apply isUnramifiedAdjoin_of_inertiaDegree_eq_finrank
    rw [hrank]
    have hmul := ramificationIndex_mul_inertiaDegree K x
    rw [hrank] at hmul
    have he : ramificationIndex K x ≠ 0 := by
      intro h; rw [h, zero_mul] at hmul; omega
    have hle : inertiaDegree K x ≤ ramificationIndex K x * inertiaDegree K x :=
      Nat.le_mul_of_pos_left _ (Nat.pos_of_ne_zero he)
    omega
  have hnormal : Normal K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    normal_of_splits_in_adjoin K x (f.map (algebraMap 𝒪[K.carrier] K.carrier)) hFm.ne_zero
      (splits_adjoin_of_lift K x f hfm hgi _ hmpbar) haev
  have hunit := isUnit_eval_derivative K x f hgi _ hbar
  exact ⟨hrank, hunram, hnormal,
    residueGalHom_bijective K x hnormal hunram
      (residueGalHom_injective K x f hnorm hrootO hunit)⟩


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

/-- **★不分岐拡大の構成**——剰余体上のモニック既約多項式 `g` から、
次数 `deg g` の**不分岐**単項拡大 `K(x)/K` を作る。`g` を `𝒪[K.carrier]`
へ持ち上げて `K.closure` の中の根を取り、判定子
`isUnramifiedAdjoin_of_lift` に流し込むだけ。 -/
theorem exists_isUnramifiedAdjoin_of_irreducible {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (g : Polynomial 𝓀[K.carrier]) (hgm : g.Monic) (hgi : Irreducible g) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        = g.natDegree ∧ IsUnramifiedAdjoin K x
        ∧ Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ∧ Function.Bijective (residueGalHom K x) := by
  have hlifts : g ∈ Polynomial.lifts (IsLocalRing.residue 𝒪[K.carrier]) :=
    Polynomial.lifts_iff_coeff_lifts g |>.mpr (fun n => ⟨_, Quotient.out_eq' _⟩)
  obtain ⟨f, hfmap, hfdeg, hfm⟩ := Polynomial.lifts_and_degree_eq_and_monic hlifts hgm
  have hgnd : g.natDegree = f.natDegree := by rw [← hfmap, hfm.natDegree_map]
  have hgi' : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])) := by
    rw [hfmap]; exact hgi
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
  obtain ⟨hrank, hunram, hnormal, hbij⟩ := isUnramifiedAdjoin_of_lift K x f hfm hgi' haev
  exact ⟨x, by rw [hrank, hgnd], hunram, hnormal, hbij⟩

/-- **★各次数の不分岐拡大の存在**——`n ≥ 1` に対し、`K` の次数 `n` の
**不分岐**単項拡大 `K(x)/K` が `K.closure` の中に存在し、しかも
`Normal`(標数 0 なので Galois)である。
`Found/FiniteFieldIrreducible.lean`(有限体上の任意次数のモニック既約
多項式)を上の構成に流し込むだけ。

これで「不分岐拡大の理論」の第一の柱——**存在**——が立った。
残るのは一意性(同じ次数の不分岐拡大は一致)と、`K^ur` および
`Gal(K^ur/K) ≅ Ẑ`。 -/
theorem exists_isUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ∧ Function.Bijective (residueGalHom K x) := by
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  obtain ⟨g, hgm, hgi, hgd⟩ := ABC3.Found.exists_monic_irreducible_natDegree_eq 𝓀[K.carrier] n hn
  obtain ⟨x, hrank, hu, hnor, hbij⟩ := exists_isUnramifiedAdjoin_of_irreducible K g hgm hgi
  exact ⟨x, hgd ▸ hrank, hu, hnor, hbij⟩

/-- 標数 0 なので、上で得た不分岐拡大は **Galois**。 -/
theorem exists_isGalois_isUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  obtain ⟨x, hrank, hu, hnor, _⟩ := exists_isUnramifiedAdjoin K n hn
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := hnor
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  exact ⟨x, hrank, hu, ⟨⟩⟩


/-- **★★不分岐拡大の Galois 群は剰余体の Galois 群と同型**——各次数 `n`
について、`Gal(K(x)/K) ≃* Gal(𝓀_{K(x)}/𝓀)`。右辺は有限体の Galois 群
なので Frobenius が生成する巡回群 `ℤ/n`——これが `Gal(K^ur/K) ≅ Ẑ` の
出どころになる。 -/
theorem exists_mulEquiv_residueGal {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ Nonempty (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
          ≃* (IsLocalRing.ResidueField (adjoinIntegers K x)
            ≃ₐ[𝓀[K.carrier]] IsLocalRing.ResidueField (adjoinIntegers K x))) := by
  obtain ⟨x, hrank, hu, _, hbij⟩ := exists_isUnramifiedAdjoin K n hn
  exact ⟨x, hrank, hu, ⟨MulEquiv.ofBijective (residueGalHom K x) hbij⟩⟩


/-! ## 不分岐拡大の Galois 群は `ℤ/n`

`Gal(K(x)/K) ≃* Gal(𝓀_{K(x)}/𝓀)` と、**有限体の Galois 群は巡回群**
(mathlib のインスタンス `IsCyclic (E ≃ₐ[F] E)`——Frobenius が生成する)
を合わせると、次数 `n` の不分岐拡大の Galois 群は位数 `n` の巡回群、
すなわち `ℤ/n` と同型。

`Gal(K^ur/K) ≅ Ẑ` は、これらを `n` について射影極限に組み上げたもの
——残る仕事は一意性(次数 `n` の不分岐拡大が `K.closure` の中で一意)と
極限の構成。 -/

/-- **次数 `n` の不分岐拡大の Galois 群は位数 `n` の巡回群**。 -/
theorem exists_isCyclic_gal {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ IsCyclic ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
        ∧ Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = n := by
  obtain ⟨x, hrank, hu, hnor, hbij⟩ := exists_isUnramifiedAdjoin K n hn
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := hnor
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  refine ⟨x, hrank, hu, ?_, ?_⟩
  · exact (MulEquiv.isCyclic (MulEquiv.ofBijective (residueGalHom K x) hbij).symm).mp
      inferInstance
  · rw [IsGalois.card_aut_eq_finrank, hrank]

/-- **★★★次数 `n` の不分岐拡大の Galois 群は `ℤ/n`**。 -/
theorem exists_gal_mulEquiv_zmod {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ Nonempty (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
            ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
          ≃* Multiplicative (ZMod n)) := by
  obtain ⟨x, hrank, hu, hcyc, hcard⟩ := exists_isCyclic_gal K n hn
  haveI := hcyc
  refine ⟨x, hrank, hu, ⟨?_⟩⟩
  exact (hcard ▸ (zmodCyclicMulEquiv hcyc)).symm


/-- 剰余体にある `f̄` の根は、Hensel で `adjoinIntegers K x` の中の `f` の根に
持ち上がる(`Found/HenselianSplits.lean` の一般補題を本設定に流し込んだもの)。 -/
theorem exists_root_lift_of_residue_root {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) (f : Polynomial 𝒪[K.carrier]) (hfm : f.Monic)
    (hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])))
    (b : IsLocalRing.ResidueField (adjoinIntegers K x))
    (hb : Polynomial.aeval b (f.map (IsLocalRing.residue 𝒪[K.carrier])) = 0) :
    ∃ z : adjoinIntegers K x,
      Polynomial.aeval (((z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        : K.closure)) (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0 := by
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  have hcomp : (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))).comp
      (IsLocalRing.residue 𝒪[K.carrier])
      = (IsLocalRing.residue (adjoinIntegers K x)).comp
        (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) := RingHom.ext (congrFun rfl)
  have hFAres : (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))
      = (f.map (IsLocalRing.residue 𝒪[K.carrier])).map (algebraMap 𝓀[K.carrier]
        (IsLocalRing.ResidueField (adjoinIntegers K x))) := by
    rw [Polynomial.map_map, ← hcomp, ← Polynomial.map_map]
  have hr0 : Polynomial.eval b ((f.map (IsLocalRing.residue 𝒪[K.carrier])).map
      (algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)))) = 0 := by
    rw [Polynomial.eval_map, ← Polynomial.aeval_def]; exact hb
  have hroot0 : Polynomial.eval b ((f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
      (IsLocalRing.residue (adjoinIntegers K x))) = 0 := by rw [hFAres]; exact hr0
  have hder0 : Polynomial.eval b (Polynomial.derivative
      ((f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))).map
        (IsLocalRing.residue (adjoinIntegers K x)))) ≠ 0 := by
    rw [hFAres]
    obtain ⟨u, v, huv⟩ := (PerfectField.separable_of_irreducible hgi).map
      (f := algebraMap 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)))
    intro hcon
    have hev2 := congrArg (Polynomial.eval b) huv
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_one, hcon, hr0,
      mul_zero, zero_add] at hev2
    exact one_ne_zero hev2.symm
  obtain ⟨zA, hzA, hzAres⟩ := ABC3.Found.exists_root_of_residue_root
    (f.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))) (hfm.map _) b hroot0 hder0
  have hzA' : Polynomial.eval₂ (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) zA f = 0 := by
    rw [← Polynomial.eval_map]; exact hzA
  refine ⟨zA, ?_⟩
  have hι : (((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val.toRingHom).comp
      (Subring.subtype (adjoinIntegers K x))).comp
      (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      = algebraMap 𝒪[K.carrier] K.closure := rfl
  have h2 := Polynomial.hom_eval₂ f (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
    (((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val.toRingHom).comp
      (Subring.subtype (adjoinIntegers K x))) zA
  rw [hι, hzA', map_zero] at h2
  rw [Polynomial.aeval_def, Polynomial.eval₂_map]
  exact h2.symm

/-! ## ★不分岐拡大の一意性

「不分岐 ⟹ 剰余体の原始元を Hensel で持ち上げると、それが `K(x)` 全体の
生成元になる」(`exists_integral_generator`)——不分岐性からその
**定義多項式**を取り出す操作。これがあると、任意の不分岐拡大に判定子
`isUnramifiedAdjoin_of_lift` を適用でき、次が出る:

* `normal_of_isUnramifiedAdjoin` : 不分岐 ⟹ normal(標数 0 なので Galois)
* `adjoin_eq_of_isUnramified` : **同じ次数の不分岐拡大は `K.closure` の
  中で一致する**(一意性)

一意性の筋: `K(y)` の定義多項式 `f`(その剰余 `f̄` は既約 `n` 次)を取る。
`𝓀_{K(x)}` も `𝓀` の `n` 次拡大だから `f̄` の根 `b` を含む
(`exists_root_of_finrank_eq`)。`b` を Hensel で `𝒪_{K(x)}` へ持ち上げて
`f` の根 `z ∈ K(x)` を得る。`z` は `f_K` の根で `f_K` は既約 `n` 次だから
`[K(z):K] = n`、`K(z) ⊆ K(x)` と次数一致で `K(z) = K(x)`。一方 `K(y)` は
normal で `f_K` は `K(y)` で分裂するから `z ∈ K(y)`、同様に `K(z) = K(y)`。 -/

/-- **不分岐拡大の「整な生成元」**——剰余体の原始元を Hensel で持ち上げると
`K(x)` の生成元になり、その最小多項式は剰余体で既約な多項式の持ち上げ。 -/
theorem exists_integral_generator {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x : K.closure)
    (hu : IsUnramifiedAdjoin K x) :
    ∃ (θ : K.closure) (f : Polynomial 𝒪[K.carrier]),
      f.Monic ∧ Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier]))
        ∧ Polynomial.aeval θ (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0
        ∧ IntermediateField.adjoin K.carrier ({θ} : Set K.closure)
          = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  obtain ⟨tb, htb⟩ := Field.exists_primitive_element 𝓀[K.carrier]
    (IsLocalRing.ResidueField (adjoinIntegers K x))
  have htbint : IsIntegral 𝓀[K.carrier] tb := IsIntegral.of_finite _ _
  have hgbm : (minpoly 𝓀[K.carrier] tb).Monic := minpoly.monic htbint
  have hgbi : Irreducible (minpoly 𝓀[K.carrier] tb) := minpoly.irreducible htbint
  have hgbdeg : (minpoly 𝓀[K.carrier] tb).natDegree
      = Module.finrank 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
    rw [← IntermediateField.adjoin.finrank htbint, htb, IntermediateField.finrank_top']
  have hlifts : (minpoly 𝓀[K.carrier] tb) ∈ Polynomial.lifts (IsLocalRing.residue 𝒪[K.carrier]) :=
    Polynomial.lifts_iff_coeff_lifts _ |>.mpr (fun n => ⟨_, Quotient.out_eq' _⟩)
  obtain ⟨f, hfmap, hfdeg, hfm⟩ := Polynomial.lifts_and_degree_eq_and_monic hlifts hgbm
  have hgi : Irreducible (f.map (IsLocalRing.residue 𝒪[K.carrier])) := by rw [hfmap]; exact hgbi
  have hbar : Polynomial.aeval tb (f.map (IsLocalRing.residue 𝒪[K.carrier])) = 0 := by
    rw [hfmap]; exact minpoly.aeval _ _
  obtain ⟨θA, haev⟩ := exists_root_lift_of_residue_root K x f hfm hgi tb hbar
  refine ⟨((θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure), f,
    hfm, hgi, haev, ?_⟩
  obtain ⟨hrankθ, -, -, -⟩ := isUnramifiedAdjoin_of_lift K
    (((θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) f hfm hgi haev
  have hnd : (minpoly 𝓀[K.carrier] tb).natDegree = f.natDegree := by
    rw [← hfmap, hfm.natDegree_map]
  refine IntermediateField.eq_of_le_of_finrank_eq ?_ ?_
  · rw [IntermediateField.adjoin_simple_le_iff]
    exact (θA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2
  · rw [hrankθ, ← hnd, hgbdeg, ← inertiaDegree_eq_finrank_residueField K x,
      inertiaDegree_eq_finrank_of_isUnramified K x hu]

/-- **不分岐 ⟹ normal**(標数 0 なので Galois)。 -/
theorem normal_of_isUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x : K.closure) (hu : IsUnramifiedAdjoin K x) :
    Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  obtain ⟨θ, f, hfm, hgi, haev, hadj⟩ := exists_integral_generator K x hu
  obtain ⟨-, -, hnor, -⟩ := isUnramifiedAdjoin_of_lift K θ f hfm hgi haev
  exact hadj ▸ hnor

/-- **★★不分岐拡大の一意性**——同じ次数の不分岐拡大は `K.closure` の中で
一致する。 -/
theorem adjoin_eq_of_isUnramified {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (x y : K.closure) (hux : IsUnramifiedAdjoin K x) (huy : IsUnramifiedAdjoin K y)
    (hdeg : Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      = IntermediateField.adjoin K.carrier ({y} : Set K.closure) := by
  obtain ⟨θ, f, hfm, hgi, haev, hadj⟩ := exists_integral_generator K y huy
  obtain ⟨hrankθ, huθ, hnorθ, -⟩ := isUnramifiedAdjoin_of_lift K θ f hfm hgi haev
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  have hfKi : Irreducible (f.map (algebraMap 𝒪[K.carrier] K.carrier)) :=
    (hfm.irreducible_iff_irreducible_map_fraction_map (K := K.carrier)).mp
      (hfm.irreducible_of_irreducible_map _ f hgi)
  have hFm : (f.map (algebraMap 𝒪[K.carrier] K.carrier)).Monic := hfm.map _
  have hdegres : Module.finrank 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K x))
      = (f.map (IsLocalRing.residue 𝒪[K.carrier])).natDegree := by
    rw [← inertiaDegree_eq_finrank_residueField K x,
      inertiaDegree_eq_finrank_of_isUnramified K x hux, hdeg, ← hadj, hrankθ, hfm.natDegree_map]
  have hrankx : Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = f.natDegree := by
    rw [← inertiaDegree_eq_finrank_of_isUnramified K x hux,
      inertiaDegree_eq_finrank_residueField K x, hdegres, hfm.natDegree_map]
  obtain ⟨b, hb⟩ := ABC3.Found.exists_root_of_finrank_eq 𝓀[K.carrier]
    (f.map (IsLocalRing.residue 𝒪[K.carrier])) hgi
    (IsLocalRing.ResidueField (adjoinIntegers K x)) hdegres
  obtain ⟨zA, haevz⟩ := exists_root_lift_of_residue_root K x f hfm hgi b hb
  obtain ⟨hrankz, -, -, -⟩ := isUnramifiedAdjoin_of_lift K
    ((zA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) f hfm hgi haevz
  have hzx : IntermediateField.adjoin K.carrier
      ({((zA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)}
        : Set K.closure)
      = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    refine IntermediateField.eq_of_le_of_finrank_eq ?_ ?_
    · rw [IntermediateField.adjoin_simple_le_iff]
      exact (zA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2
    · rw [hrankz, hrankx]
  have hvalθ : (IntermediateField.adjoin K.carrier ({θ} : Set K.closure)).val
      ⟨θ, IntermediateField.mem_adjoin_simple_self K.carrier θ⟩ = θ := rfl
  have hkeyθ := Polynomial.aeval_algHom_apply
    (IntermediateField.adjoin K.carrier ({θ} : Set K.closure)).val
    (⟨θ, IntermediateField.mem_adjoin_simple_self K.carrier θ⟩ :
      IntermediateField.adjoin K.carrier ({θ} : Set K.closure))
    (f.map (algebraMap 𝒪[K.carrier] K.carrier))
  rw [hvalθ, haev] at hkeyθ
  have h1θ : Polynomial.aeval
      (⟨θ, IntermediateField.mem_adjoin_simple_self K.carrier θ⟩ :
        IntermediateField.adjoin K.carrier ({θ} : Set K.closure))
      (f.map (algebraMap 𝒪[K.carrier] K.carrier)) = 0 :=
    (IntermediateField.adjoin K.carrier ({θ} : Set K.closure)).val.injective
      (by rw [map_zero]; exact hkeyθ.symm)
  have hmpθL : f.map (algebraMap 𝒪[K.carrier] K.carrier)
      = minpoly K.carrier (⟨θ, IntermediateField.mem_adjoin_simple_self K.carrier θ⟩ :
        IntermediateField.adjoin K.carrier ({θ} : Set K.closure)) :=
    minpoly.eq_of_irreducible_of_monic hfKi h1θ hFm
  have hsplitsθ : ((f.map (algebraMap 𝒪[K.carrier] K.carrier)).map
      (algebraMap K.carrier
        (IntermediateField.adjoin K.carrier ({θ} : Set K.closure)))).Splits := by
    rw [hmpθL]; exact Normal.splits hnorθ _
  have hz_mem := mem_adjoin_of_root_of_splits K θ (f.map (algebraMap 𝒪[K.carrier] K.carrier))
    hFm.ne_zero hsplitsθ haevz
  have hzθ : IntermediateField.adjoin K.carrier
      ({((zA : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)}
        : Set K.closure)
      = IntermediateField.adjoin K.carrier ({θ} : Set K.closure) := by
    refine IntermediateField.eq_of_le_of_finrank_eq ?_ ?_
    · rw [IntermediateField.adjoin_simple_le_iff]; exact hz_mem
    · rw [hrankz, hrankθ]
  exact hzx.symm.trans (hzθ.trans hadj)

/-- **★★★各次数の不分岐拡大は存在して一意**——`Gal(K^ur/K) ≅ Ẑ` へ向けた
足場がこれで揃う。 -/
theorem exists_unique_adjoin_isUnramified {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ x : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) = n
        ∧ IsUnramifiedAdjoin K x
        ∧ ∀ y : K.closure, IsUnramifiedAdjoin K y →
            Module.finrank K.carrier
              (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) = n →
            IntermediateField.adjoin K.carrier ({y} : Set K.closure)
              = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  obtain ⟨x, hrank, hu, -, -⟩ := exists_isUnramifiedAdjoin K n hn
  refine ⟨x, hrank, hu, fun y huy hranky => ?_⟩
  exact adjoin_eq_of_isUnramified K y x huy hu (by rw [hranky, hrank])


/-! ## ★最大不分岐拡大 `K^ur`

存在・一意性が揃ったので、不分岐拡大の全体は**有向系**になる:
次数 `m`・`n` の不分岐拡大は、次数 `m*n` の不分岐拡大に両方とも含まれる
(`adjoin_le_of_dvd`——次数が割り切れば包含する)。したがって
`K^ur := ⨆ {K(x) | x は不分岐}` は「有向和」であり、`z ∈ K^ur` は
「ある不分岐 `x` について `z ∈ K(x)`」と同値。

`Gal(K^ur/K) ≅ Ẑ` は、`Gal(K_n/K) ≅ ℤ/n`(`exists_gal_mulEquiv_zmod`)を
`n` について射影極限に組み上げたもの——次の段。 -/

/-- **次数が割り切れば包含する**——`[K(x):K] ∣ [K(y):K]` なら `K(x) ⊆ K(y)`。
`f̄` の根が `𝓀_{K(y)}` にあること(`exists_root_of_natDegree_dvd`)を
Hensel で持ち上げ、一意性で `K(x)` と同定する。 -/
theorem adjoin_le_of_dvd {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x y : K.closure)
    (hux : IsUnramifiedAdjoin K x) (huy : IsUnramifiedAdjoin K y)
    (hdvd : Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ∣ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure) := by
  obtain ⟨θ, f, hfm, hgi, haev, hadj⟩ := exists_integral_generator K x hux
  obtain ⟨hrankθ, -, -, -⟩ := isUnramifiedAdjoin_of_lift K θ f hfm hgi haev
  haveI := module_finite_adjoinIntegers K y
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  have hdvd' : (f.map (IsLocalRing.residue 𝒪[K.carrier])).natDegree
      ∣ Module.finrank 𝓀[K.carrier] (IsLocalRing.ResidueField (adjoinIntegers K y)) := by
    rw [hfm.natDegree_map, ← inertiaDegree_eq_finrank_residueField K y,
      inertiaDegree_eq_finrank_of_isUnramified K y huy, ← hrankθ, hadj]
    exact hdvd
  obtain ⟨b, hb⟩ := ABC3.Found.exists_root_of_natDegree_dvd 𝓀[K.carrier]
    (f.map (IsLocalRing.residue 𝒪[K.carrier])) hgi
    (IsLocalRing.ResidueField (adjoinIntegers K y)) hdvd'
  obtain ⟨zA, haevz⟩ := exists_root_lift_of_residue_root K y f hfm hgi b hb
  obtain ⟨hrankz, huz, -, -⟩ := isUnramifiedAdjoin_of_lift K
    ((zA : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) : K.closure) f hfm hgi haevz
  have hzx : IntermediateField.adjoin K.carrier
      ({((zA : IntermediateField.adjoin K.carrier ({y} : Set K.closure)) : K.closure)}
        : Set K.closure)
      = IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    refine adjoin_eq_of_isUnramified K _ x huz hux ?_
    rw [hrankz, ← hrankθ, hadj]
  rw [← hzx, IntermediateField.adjoin_simple_le_iff]
  exact (zA : IntermediateField.adjoin K.carrier ({y} : Set K.closure)).2

/-- **不分岐拡大は有向系をなす**——次数 `m*n` の不分岐拡大が両方を含む。 -/
theorem exists_isUnramified_ge {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x y : K.closure)
    (hux : IsUnramifiedAdjoin K x) (huy : IsUnramifiedAdjoin K y) :
    ∃ z : K.closure, IsUnramifiedAdjoin K z
      ∧ IntermediateField.adjoin K.carrier ({x} : Set K.closure)
          ≤ IntermediateField.adjoin K.carrier ({z} : Set K.closure)
      ∧ IntermediateField.adjoin K.carrier ({y} : Set K.closure)
          ≤ IntermediateField.adjoin K.carrier ({z} : Set K.closure) := by
  have hmpos : 0 < Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := Module.finrank_pos
  have hnpos : 0 < Module.finrank K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) := Module.finrank_pos
  obtain ⟨z, hrankz, huz, -, -⟩ := exists_isUnramifiedAdjoin K
    (Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      * Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)))
    (by positivity)
  exact ⟨z, huz,
    adjoin_le_of_dvd K x z hux huz (by rw [hrankz]; exact Dvd.intro _ rfl),
    adjoin_le_of_dvd K y z huy huz (by rw [hrankz]; exact Dvd.intro_left _ rfl)⟩

/-- **`K` の最大不分岐拡大 `K^ur`**——不分岐な単項拡大すべての上限。 -/
noncomputable def unramifiedClosure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IntermediateField K.carrier K.closure :=
  ⨆ x : {x : K.closure // IsUnramifiedAdjoin K x},
    IntermediateField.adjoin K.carrier ({(x : K.closure)} : Set K.closure)

/-- 不分岐な元は存在する(次数 1、すなわち `K` 自身)。 -/
theorem nonempty_isUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Nonempty {x : K.closure // IsUnramifiedAdjoin K x} := by
  obtain ⟨x, -, hu, -, -⟩ := exists_isUnramifiedAdjoin K 1 one_ne_zero
  exact ⟨⟨x, hu⟩⟩

/-- 不分岐な単項拡大の族は有向。 -/
theorem directed_isUnramifiedAdjoin {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Directed (· ≤ ·) (fun x : {x : K.closure // IsUnramifiedAdjoin K x} =>
      IntermediateField.adjoin K.carrier ({(x : K.closure)} : Set K.closure)) := by
  rintro ⟨a, ha⟩ ⟨b, hb⟩
  obtain ⟨z, huz, h1, h2⟩ := exists_isUnramified_ge K a b ha hb
  exact ⟨⟨z, huz⟩, h1, h2⟩

/-- **`K^ur` の元の特徴づけ**——`z ∈ K^ur` ⟺ ある不分岐 `x` について
`z ∈ K(x)`。有向性(`directed_isUnramifiedAdjoin`)から上限が和集合になる。 -/
theorem mem_unramifiedClosure_iff {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (z : K.closure) :
    z ∈ unramifiedClosure K ↔ ∃ x : K.closure, IsUnramifiedAdjoin K x
      ∧ z ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  haveI := nonempty_isUnramifiedAdjoin K
  constructor
  · intro hz
    have hcoe := IntermediateField.coe_iSup_of_directed (directed_isUnramifiedAdjoin K)
    have hz' : z ∈ (⋃ i : {x : K.closure // IsUnramifiedAdjoin K x},
        ((IntermediateField.adjoin K.carrier ({(i : K.closure)} : Set K.closure) :
          IntermediateField K.carrier K.closure) : Set K.closure)) := hcoe ▸ hz
    obtain ⟨i, hzi⟩ := Set.mem_iUnion.mp hz'
    exact ⟨(i : K.closure), i.2, hzi⟩
  · rintro ⟨x, hx, hz⟩
    exact le_iSup (fun x : {x : K.closure // IsUnramifiedAdjoin K x} =>
      IntermediateField.adjoin K.carrier ({(x : K.closure)} : Set K.closure)) ⟨x, hx⟩ hz

/-- 不分岐な単項拡大は `K^ur` に含まれる。 -/
theorem adjoin_le_unramifiedClosure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {x : K.closure} (hx : IsUnramifiedAdjoin K x) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤ unramifiedClosure K :=
  le_iSup (fun x : {x : K.closure // IsUnramifiedAdjoin K x} =>
    IntermediateField.adjoin K.carrier ({(x : K.closure)} : Set K.closure)) ⟨x, hx⟩


/-! ## `K^ur/K` は Galois、そして `Γ_K ↠ ℤ/n`

`K^ur` は normal な単項拡大の上限なので normal(`IntermediateField.
normal_iSup`)、標数 0 なので **Galois**。

さらに、次数 `n` の不分岐拡大 `K_n` は normal だから制限射
`Γ_K = Gal(K̄/K) →* Gal(K_n/K)` は**全射**
(`AlgEquiv.restrictNormalHom_surjective`)。`Gal(K_n/K) ≃ ℤ/n` と合わせて

**`Γ_K` は任意の `n ≥ 1` に対して `ℤ/n` へ全射する**

——これが `Γ_K ↠ Gal(K^ur/K) ≅ Ẑ` の具体的な中身。`Ẑ` 自体は
mathlib に見当たらない(2026-09-05 実測)ので、射影極限としての
`Ẑ` を作る代わりに、各 `n` での全射という形で述べる。 -/

/-- `K^ur/K` は normal——normal な単項拡大の上限だから。 -/
theorem normal_unramifiedClosure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Normal K.carrier (unramifiedClosure K) := by
  haveI : ∀ i : {x : K.closure // IsUnramifiedAdjoin K x},
      Normal K.carrier (IntermediateField.adjoin K.carrier ({(i : K.closure)} : Set K.closure)) :=
    fun i => normal_of_isUnramifiedAdjoin K (i : K.closure) i.2
  exact IntermediateField.normal_iSup (F := K.carrier) (K := K.closure) _

/-- **`K^ur/K` は Galois**(標数 0 なので分離性は自動)。 -/
theorem isGalois_unramifiedClosure {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    IsGalois K.carrier (unramifiedClosure K) := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI := normal_unramifiedClosure K
  haveI : Algebra.IsSeparable K.carrier (unramifiedClosure K) :=
    IntermediateField.isSeparable_tower_bot K.carrier (unramifiedClosure K)
  exact ⟨⟩

/-- **★★★`Γ_K` は任意の `n ≥ 1` に対して `ℤ/n` へ全射する**——次数 `n` の
不分岐拡大 `K_n` への制限射が全射(`K_n` は normal)で、
`Gal(K_n/K) ≃ ℤ/n` だから。古典的局所類体論の不分岐部分
`Gal(K^ur/K) ≅ Ẑ` の、`Ẑ` を経由しない具体的な言い換え。 -/
theorem exists_surjective_absGal_to_zmod {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (n : ℕ)
    (hn : n ≠ 0) :
    ∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →* Multiplicative (ZMod n),
      Function.Surjective φ := by
  obtain ⟨x, hrank, hu, hnor, hbij⟩ := exists_isUnramifiedAdjoin K n hn
  haveI := hnor
  haveI := IsAlgClosure.normal K.carrier K.closure
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  have hcyc : IsCyclic ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    (MulEquiv.isCyclic (MulEquiv.ofBijective (residueGalHom K x) hbij).symm).mp inferInstance
  have hcard : Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = n := by
    rw [IsGalois.card_aut_eq_finrank, hrank]
  haveI := hcyc
  have e := (hcard ▸ (zmodCyclicMulEquiv hcyc)).symm
  refine ⟨(e : _ ≃* Multiplicative (ZMod n)).toMonoidHom.comp
    (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))), ?_⟩
  exact (e.surjective).comp
    (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (E := K.closure))


/-! ## Teichmüller 持ち上げ——`𝒪_K^×` は `𝓀^×`(位数 `q-1` の巡回群)を含む

`Found/Teichmuller.lean`(一般の結果)を `𝒪[K.carrier]` に適用する。
`𝒪[K.carrier]` は Henselian(`henselianLocalRing_carrierIntegers`)で
剰余体は有限(`residueField_finite`)なので、剰余写像
`𝒪_K^× → 𝓀^×` は群としての切断を持つ。

古典的局所類体論の `𝒪_K^× ≅ μ_{q-1} × (1+𝔪_K)` の**第一因子**——
[pGC] Proposition 1.2 が `Γ_K^ab` の群構造から `q` を読み取るときの入口。 -/

/-- **`𝒪_K^×` の Teichmüller 持ち上げ**——剰余写像 `𝒪_K^× → 𝓀^×` の
群としての切断が存在する(したがって `𝓀^×` は `𝒪_K^×` の部分群として
実現される)。 -/
theorem exists_teichmuller_section {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    ∃ ω : (𝓀[K.carrier])ˣ →* (𝒪[K.carrier])ˣ,
      (∀ b, Units.map (IsLocalRing.residue 𝒪[K.carrier] :
        𝒪[K.carrier] →* 𝓀[K.carrier]) (ω b) = b) ∧ Function.Injective ω := by
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
  exact ⟨ABC3.Found.teichmuller 𝒪[K.carrier],
    ABC3.Found.residue_teichmuller 𝒪[K.carrier],
    ABC3.Found.teichmuller_injective 𝒪[K.carrier]⟩


/-! ## `𝒪_K^× ≅ 𝓀^× × (1+𝔪_K)`——単数群の直積分解

Teichmüller 切断(`exists_teichmuller_section`)と
`ABC3.Found.prodKerOfRightInverse`(切断を持つ準同型は直積分解を与える)を
合わせるだけ。第二因子(剰余が `1` の単数=主単数 `1+𝔪_K`)は pro-p 群で、
その `ℤ_p`-階数が `[K:ℚ_p]` に対応する——[pGC] Proposition 1.2 が
`Γ_K^ab` の群構造から `q` と `[K:ℚ_p]` を読み取る二つの入口のうち、
第一因子から `q`(=`|𝓀|`、`𝓀^×` の位数 +1)が読める。 -/

/-- **単数群の直積分解** `𝒪_K^× ≅ 𝓀^× × ker(還元)`。 -/
theorem nonempty_units_mulEquiv_prod {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) :
    Nonempty ((𝒪[K.carrier])ˣ ≃* (𝓀[K.carrier])ˣ ×
      (Units.map (IsLocalRing.residue 𝒪[K.carrier] :
        𝒪[K.carrier] →* 𝓀[K.carrier])).ker) := by
  obtain ⟨ω, hsec, -⟩ := exists_teichmuller_section K
  exact ⟨ABC3.Found.prodKerOfRightInverse
    (Units.map (IsLocalRing.residue 𝒪[K.carrier] : 𝒪[K.carrier] →* 𝓀[K.carrier])) ω hsec⟩

/-- 第二因子は「剰余が `1` の単数」——すなわち主単数 `1+𝔪_K`。 -/
theorem mem_ker_units_map_residue_iff {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (u : (𝒪[K.carrier])ˣ) :
    u ∈ (Units.map (IsLocalRing.residue 𝒪[K.carrier] :
        𝒪[K.carrier] →* 𝓀[K.carrier])).ker
      ↔ IsLocalRing.residue 𝒪[K.carrier] (u : 𝒪[K.carrier]) = 1 := by
  rw [MonoidHom.mem_ker, Units.ext_iff]
  rfl


/-! ## `Γ_K ↠ Gal(K^ur/K) ↠ ℤ/n`

`K^ur/K` は Galois(`isGalois_unramifiedClosure`)なので、絶対 Galois 群
からの制限射は全射(`AlgEquiv.restrictNormalHom_surjective`)。さらに
次数 `n` の不分岐拡大 `K_n ⊆ K^ur` への制限も全射で、
`Gal(K_n/K) ≃ ℤ/n`。したがって

```
Γ_K ↠ Gal(K^ur/K) ↠ ℤ/n   (任意の n ≥ 1)
```

`Ẑ = lim ℤ/n` を作れば右側は `Gal(K^ur/K) ≅ Ẑ` になる(`Ẑ` は mathlib
に不在、2026-09-05 実測)。 -/

/-- **`Γ_K ↠ Gal(K^ur/K)`**——`K^ur/K` が normal だから。 -/
theorem exists_surjective_absGal_to_unramifiedClosureGal {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) :
    ∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →*
      ((unramifiedClosure K) ≃ₐ[K.carrier] (unramifiedClosure K)), Function.Surjective φ := by
  haveI := normal_unramifiedClosure K
  haveI := IsAlgClosure.normal K.carrier K.closure
  exact ⟨AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K),
    AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (unramifiedClosure K)) (E := K.closure)⟩

/-- **`Gal(K^ur/K) ↠ ℤ/n`**——次数 `n` の不分岐拡大 `K_n ⊆ K^ur` への
制限が全射で、`Gal(K_n/K) ≃ ℤ/n` だから。 -/
theorem exists_surjective_unramifiedClosureGal_to_zmod {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) (n : ℕ) (hn : n ≠ 0) :
    ∃ ψ : ((unramifiedClosure K) ≃ₐ[K.carrier] (unramifiedClosure K))
      →* Multiplicative (ZMod n), Function.Surjective ψ := by
  obtain ⟨x, hrank, hu, hnor, hbij⟩ := exists_isUnramifiedAdjoin K n hn
  haveI := hnor
  haveI := normal_unramifiedClosure K
  haveI := module_finite_adjoinIntegers K x
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  have hcyc : IsCyclic ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    (MulEquiv.isCyclic (MulEquiv.ofBijective (residueGalHom K x) hbij).symm).mp inferInstance
  have hcard : Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) = n := by
    rw [IsGalois.card_aut_eq_finrank, hrank]
  haveI := hcyc
  have e := (hcard ▸ (zmodCyclicMulEquiv hcyc)).symm
  letI : Algebra (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (unramifiedClosure K) :=
    (IntermediateField.inclusion (adjoin_le_unramifiedClosure K hu)).toRingHom.toAlgebra
  haveI : IsScalarTower K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (unramifiedClosure K) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  refine ⟨(e : _ ≃* Multiplicative (ZMod n)).toMonoidHom.comp
    (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := (unramifiedClosure K))
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))), ?_⟩
  exact (e.surjective).comp
    (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
      (K₁ := (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      (E := (unramifiedClosure K)))


/-- 単項拡大の包含は次数の整除を導く(塔の公式)。不分岐性は不要。 -/
theorem finrank_dvd_of_adjoin_le {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {x y : K.closure}
    (hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ∣ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) := by
  letI : Algebra (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    (IntermediateField.inclusion hle).toRingHom.toAlgebra
  haveI : IsScalarTower K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    IsScalarTower.of_algebraMap_eq (fun _ => rfl)
  exact ⟨Module.finrank (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)),
    (Module.finrank_mul_finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure))).symm⟩

/-- **★不分岐拡大の塔は次数の整除で決まる**——`K_m ⊆ K_n ⟺ m ∣ n`。
`⟸` は `adjoin_le_of_dvd`(Hensel + 一意性)、`⟹` は塔の公式。
`K^ur` の内部構造が `(ℕ, ∣)` の順序で完全に記述されることを意味する
——`Gal(K^ur/K) ≅ Ẑ = lim ℤ/n` の「順序系」の側。 -/
theorem adjoin_le_iff_dvd {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (x y : K.closure)
    (hux : IsUnramifiedAdjoin K x) (huy : IsUnramifiedAdjoin K y) :
    IntermediateField.adjoin K.carrier ({x} : Set K.closure)
        ≤ IntermediateField.adjoin K.carrier ({y} : Set K.closure)
      ↔ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ∣ Module.finrank K.carrier
          (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
  ⟨fun hle => finrank_dvd_of_adjoin_le K hle, adjoin_le_of_dvd K x y hux huy⟩

end ABC3.Found.PGC
