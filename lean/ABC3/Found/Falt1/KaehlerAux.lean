import Mathlib.RingTheory.Kaehler.Basic
import Mathlib.RingTheory.Kaehler.Polynomial
import Mathlib.LinearAlgebra.Span.Basic
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.RingTheory.DedekindDomain.Different
import Mathlib.RingTheory.IsAdjoinRoot

/-!
# [Falt1] Lemma 1.1 に向けた補助補題(`Found`、sorry 無し)

★★これは `/goal Falt1 Chapter I Found` の一部として、Lemma 1.1
(`Ω_V ⊗_V W → Ω_W` の単射性・余核の長さ)を実際に **証明** する試みの
第一歩。**Lemma 1.1 そのものはまだ証明できていない**——ここにあるのは
その途中で使う再利用可能な補題群である(正直な進捗として記録する)。

## 完成した部分(2026-09-04、すべて sorry 無し)

**① 余核 `Ω_{W/V}` の計算**(`R := V`・`B := W` と代入すれば Lemma 1.1
の余核の主張そのもの): `B = Polynomial R ⧸ (f)` のとき
`Ω_{B/R} ≅ B ⧸ (f'(root))` ——`omega_quotient_eq_derivative_span`
(補助: `kerQuotEquivOmega`・`tensorPolynomialOmegaEquiv(_D)`・
`range_kerCotangentToTensor_span`・`map_ker_mapBaseChange_eq_span`)。
Faltings の原文「`Ω_W` は `Ω_V⊗W⊕WdT` を `f'(w)dT` で割った商」の
核心計算そのものを検証した。

**② differentIdeal への接続**: `w` が `W` を `V`-代数として生成する
(かつ拡大体レベルでも生成する)なら `differentIdeal V W = span{f'(w)}`
——`differentIdeal_eq_span_derivative`(`conductor_mul_differentIdeal`+
`conductor_eq_top_of_adjoin_eq_top`)。①と合わせれば
`Ω_{W/V} ≅ W ⧸ differentIdeal V W` が見える。

**③ 一般の `W` への橋渡し(1本目)**: `V` が整閉整域で `w` が `V` 上整
かつ `W` を生成するなら `AdjoinRoot (minpoly V w) ≃ₐ[V] W`
——`adjoinRootMinpolyEquiv`(`IsAdjoinRootMonic.mkOfAdjoinEqTop`、
体限定の `AlgEquiv.adjoinSingletonEquivAdjoinRootMinpoly` ではなく
一般整閉整域版を使った)。`AdjoinRoot f = Polynomial V ⧸ (f)` なので、
①の具体形と Falt1 の一般の `W` を繋ぐ最初の橋。

**④ 代数同型に沿った `Ω` の転送(完成)**: 環同型 `e : A ≃ₐ[R] B` は
`Ω[A⁄R] ≃+ Ω[B⁄R]` を誘導する——`omegaCongr`(mathlib には代数の**塔**
用の `KaehlerDifferential.map` はあるが、代数**同型**用の transport の
既製品が無かったため、`e`・`e.symm` から作った局所インスタンスで
両方向の `map` を組み立てて互いに逆写像であることを示した、
`omegaCongr_leftInv`/`_rightInv`・`isScalarTower_of_algEquiv`)。

**①②③④の貼り合わせ(完成)**: `falt1CokernelIso` として1つの定理に
組み立てた——`V` が Dedekind 整域、`w:W` が `V` 上整・`W` を生成し
(かつ拡大体レベルでも生成する)なら
`Ω_{W/V} ≅ W ⧸ differentIdeal V W`(`≃+`)。Lemma 1.1 の「`Ω_{W/V}` が
何であるか」の特定が完全に証明された。途中で見つけた鍵:
`AdjoinRoot f` は `Polynomial V ⧸ (f)` と定義的に等しいのに
`Algebra (Polynomial V) (AdjoinRoot f)` が自動では見つからない
(`AdjoinRoot` が `def` で simp-reducible ではないため)——
`AdjoinRoot.instAlgebraPolynomial`(`inferInstanceAs` で明示登記)で
解決した。

## 残っている作業(正直な記録、Lemma 1.1 完成にはまだ遠い)

1. **「長さ」への接続**: `Module.length`(`RingTheory/Length.lean`)は
   `LinearEquiv.length_eq (e : M ≃ₗ[R] N) : length R M = length R N`
   という形で不変性を持つ——**`LinearEquiv`(線形同型)を要求する**、
   `AddEquiv`(加法群としての同型)では直接使えない。`falt1CokernelIso`
   はシグネチャでの instance 解決順の問題を避けるため意図的に `≃+` に
   弱めていた(`omegaCongr` 等の内部の `KaehlerDifferential.map` 自体は
   本来 A-線形だが、型としては落としてある)——`≃ₗ[W]` まで持ち上げる
   作業がまだ残っている(2026-09-04 に `Module.length` の存在と
   `LinearEquiv.length_eq` の要件を確認)。
2. **単射性(Lemma 1.1 の本体の主張、第一完全列)**: `Ω_V ⊗_V W → Ω_W`
   (絶対微分、Z=絶対基底とする塔 Z→V→W への「第一完全列」)の単射性。
   Faltings の議論は `Ω_{V[T]/Z} ≅ Ω_V⊗V[T]⊕V[T]dT` という直和分解と
   f'(w) が非零因子であることを使うが、この直和分解の既製品も mathlib
   に無い(`polynomialEquiv`・`mvPolynomialBasis` は自明底のみ)——
   Lemma 1.1 の残り部分の中でも最も骨が折れそうな箇所。
3. Falt1 の `V` は完備離散付値環(すべて `IsDedekindDomain` 等の一般論
   でカバーされる保証はまだ無い、個別確認が必要)。

★見積り: 上記4点だけでもさらに多くの補題を要し、Lemma 1.1 の完成には
まだ日数単位の作業が必要と見られる。§2-4 は Faltings の
"almost mathematics" 自体が mathlib に無い(`Almost`/`Malcev`/
`Tannakian` を mathlib 全体検索して確認済み、2026-09-04)ため、
Lemma 1.1 よりさらに大きな、ライブラリ新設に近い規模の作業になる。
-/

namespace ABC3.Found.Falt1

/-- 非零因子 `f` が生成する単項イデアル `(f)` の cotangent 加群
`(f).Cotangent = (f)/(f)²` は `A' ⧸ (f)` に同型である。

証明: `A' ≃ₗ[A'] (f)`(`x ↦ x*f`、`f` が非零因子なので単射、値域は
定義から `(f)`)を作り、この同型の下で `(f)•⊤`(= `(f)²`)が
`(f) ⊂ A'` に対応することを見る。 -/
noncomputable def nzdCotangentEquivQuot {A' : Type*} [CommRing A'] (f : A')
    (hf : f ∈ nonZeroDivisorsRight A') :
    (Ideal.span ({f} : Set A')).Cotangent ≃ₗ[A'] A' ⧸ (Ideal.span ({f} : Set A')) := by
  set I := Ideal.span ({f} : Set A') with hI
  have hinj : Function.Injective (LinearMap.toSpanSingleton A' A' f) :=
    LinearMap.ker_eq_bot.mp (LinearMap.ker_toSpanSingleton_eq_bot_iff.mpr hf)
  set e0 := LinearEquiv.ofInjective (LinearMap.toSpanSingleton A' A' f) hinj
  have hrange : (LinearMap.toSpanSingleton A' A' f).range = I := by
    rw [LinearMap.range_toSpanSingleton]
  set e1 : A' ≃ₗ[A'] I := hrange ▸ e0
  refine Submodule.Quotient.equiv (I • (⊤ : Submodule A' I)) I e1.symm ?_
  rw [Submodule.map_smul'' I ⊤ (e1.symm : ↥I →ₗ[A'] A')]
  have h2 : Submodule.map (e1.symm : ↥I →ₗ[A'] A') ⊤ = ⊤ := by
    rw [Submodule.map_top]; exact LinearEquiv.range e1.symm
  rw [h2]
  simp

/-- `AdjoinRoot f`(`= Polynomial R ⧸ (f)`)への商写像の核は`(f)`そのもの。 -/
theorem ker_adjoinRoot_mk {R : Type*} [CommRing R] (f : Polynomial R) :
    RingHom.ker (AdjoinRoot.mk f) = Ideal.span ({f} : Set (Polynomial R)) := by
  show RingHom.ker (Ideal.Quotient.mk (Ideal.span {f})) = Ideal.span {f}
  exact Ideal.mk_ker

/-- 環の塔 `R → A → B`(`A → B` が全射)に対し、Kähler 微分の「第二完全列」
`B ⊗[A] Ω[A⁄R] → Ω[B⁄R] → 0`(`mapBaseChange` の core exactness)を
商同型の形で述べたもの: `(B⊗_A Ω_{A/R}) ⧸ ker(mapBaseChange) ≃ₗ[B] Ω_{B/R}`。 -/
noncomputable def kerQuotEquivOmega {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    (hsurj : Function.Surjective (algebraMap A B)) :
    letI M := TensorProduct A B Ω[A⁄R]
    (M ⧸ (LinearMap.ker (KaehlerDifferential.mapBaseChange R A B) : Submodule B M)) ≃ₗ[B] Ω[B⁄R] :=
  LinearMap.quotKerEquivOfSurjective _ (KaehlerDifferential.mapBaseChange_surjective R A B hsurj)

/-- 上の `ker(mapBaseChange)` の台集合は `kerCotangentToTensor` の値域に一致する
(Jacobi-Zariski 型の完全列 `exact_kerCotangentToTensor_mapBaseChange` の
集合レベルの言い換え)。`kerQuotEquivOmega` と組み合わせて、`Ω_{B/R}` を
`I/I²`(`I = ker(A→B)`)を経由して具体的に計算する足場になる。 -/
theorem ker_mapBaseChange_eq_range_kerCotangentToTensor {R A B : Type*}
    [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra A B] [Algebra R B] [IsScalarTower R A B]
    (hsurj : Function.Surjective (algebraMap A B)) :
    (LinearMap.ker (KaehlerDifferential.mapBaseChange R A B) : Set (TensorProduct A B Ω[A⁄R])) =
      Set.range (KaehlerDifferential.kerCotangentToTensor R A B) := by
  have hex := KaehlerDifferential.exact_kerCotangentToTensor_mapBaseChange R A B hsurj
  ext y
  exact hex y

/-- `A = Polynomial R` のとき、`B ⊗[A] Ω[A⁄R] ≃ₗ[B] B`(`Ω[A⁄R] ≃ₗ[A] A` を
`KaehlerDifferential.polynomialEquiv` で示した上で base change + rid)。
`Ω_{B/R} ≅ B/(f'(root))` の計算(進捗ノート参照)で使う。 -/
noncomputable def tensorPolynomialOmegaEquiv {R B : Type*} [CommRing R] [CommRing B]
    [Algebra (Polynomial R) B] [Algebra R B] [IsScalarTower R (Polynomial R) B] :
    TensorProduct (Polynomial R) B Ω[Polynomial R⁄R] ≃ₗ[B] B :=
  (LinearEquiv.baseChange (Polynomial R) B Ω[Polynomial R⁄R] (Polynomial R)
    (KaehlerDifferential.polynomialEquiv R)).trans
  (TensorProduct.AlgebraTensorModule.rid (Polynomial R) B B)

/-- **主要な中間結果**: `B = Polynomial R ⧸ (f)` のとき、
`kerCotangentToTensor` の値域は `{1 ⊗ₜ D(f)}` が張る `Polynomial R`-部分
加群に一致する。

証明の構造: `I = ker(A→B) = (f)` は(`toSpanSingleton A A f` の値域として)
`A` から生成される巡回加群であり、`I.toCotangent` は全射なので、
`kerCotangentToTensor` の値域は「生成元 `⟨f,_⟩` の像」`A`-生成する。
その像は `kerCotangentToTensor_toCotangent` より `1⊗ₜD(f)`。 -/
theorem range_kerCotangentToTensor_span {R B : Type*} [CommRing R] [CommRing B]
    [Algebra (Polynomial R) B] [Algebra R B] [IsScalarTower R (Polynomial R) B] (f : Polynomial R)
    (hB : RingHom.ker (algebraMap (Polynomial R) B) = Ideal.span ({f} : Set (Polynomial R))) :
    LinearMap.range (KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B) =
      Submodule.span (Polynomial R) {(1 : B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f)} := by
  have hmem : ∀ x : Polynomial R, x * f ∈ RingHom.ker (algebraMap (Polynomial R) B) := by
    intro x; rw [hB]; exact Ideal.mem_span_singleton'.mpr ⟨x, rfl⟩
  set φ := (LinearMap.toSpanSingleton (Polynomial R) (Polynomial R) f).codRestrict _ hmem with hφdef
  have hφsurj : Function.Surjective φ := by
    intro y
    have hy : (y : Polynomial R) ∈ Ideal.span ({f} : Set (Polynomial R)) := hB ▸ y.2
    obtain ⟨x, hx⟩ := Ideal.mem_span_singleton'.mp hy
    exact ⟨x, Subtype.ext hx⟩
  have hcomp : Function.Surjective
      ((RingHom.ker (algebraMap (Polynomial R) B)).toCotangent ∘ₗ φ) :=
    (RingHom.ker (algebraMap (Polynomial R) B)).toCotangent_surjective.comp hφsurj
  have hrange1 : LinearMap.range (KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B) =
      LinearMap.range (KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B ∘ₗ
        ((RingHom.ker (algebraMap (Polynomial R) B)).toCotangent ∘ₗ φ)) := by
    rw [LinearMap.range_comp, LinearMap.range_eq_top.mpr hcomp, Submodule.map_top]
  rw [hrange1]
  have hfB : algebraMap (Polynomial R) B f = 0 := by
    have : f ∈ RingHom.ker (algebraMap (Polynomial R) B) := hB ▸ Ideal.mem_span_singleton_self f
    exact this
  have heq : (KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B ∘ₗ
      ((RingHom.ker (algebraMap (Polynomial R) B)).toCotangent ∘ₗ φ)) =
      LinearMap.toSpanSingleton (Polynomial R) _
        ((1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f)) := by
    apply LinearMap.ext
    intro x
    show KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B
      ((RingHom.ker (algebraMap (Polynomial R) B)).toCotangent (φ x)) = _
    rw [KaehlerDifferential.kerCotangentToTensor_toCotangent]
    show (1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) (x*f))
        = x • ((1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f))
    rw [Derivation.leibniz, TensorProduct.tmul_add]
    have hzero : (1:B) ⊗ₜ[Polynomial R] (f • KaehlerDifferential.D R (Polynomial R) x) = 0 := by
      rw [TensorProduct.tmul_smul]
      show (f • (1:B)) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) x) = 0
      rw [Algebra.smul_def, hfB, zero_mul, TensorProduct.zero_tmul]
    rw [hzero, add_zero, TensorProduct.tmul_smul]
  rw [heq, LinearMap.range_toSpanSingleton]

/-- `tensorPolynomialOmegaEquiv` の下で `1 ⊗ₜ D(f)` は `f` の(A→Bを経由した)
微分 `algebraMap A B (derivative f)` に写る(= Faltings の `f'(w)`)。 -/
theorem tensorPolynomialOmegaEquiv_D {R B : Type*} [CommRing R] [CommRing B]
    [Algebra (Polynomial R) B] [Algebra R B] [IsScalarTower R (Polynomial R) B] (f : Polynomial R) :
    tensorPolynomialOmegaEquiv ((1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f)) =
      algebraMap (Polynomial R) B (Polynomial.derivative f) := by
  show (TensorProduct.AlgebraTensorModule.rid (Polynomial R) B B)
    ((LinearEquiv.baseChange (Polynomial R) B Ω[Polynomial R⁄R] (Polynomial R)
      (KaehlerDifferential.polynomialEquiv R)) ((1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f))) = _
  rw [LinearEquiv.baseChange_tmul, KaehlerDifferential.polynomialEquiv_D,
    TensorProduct.AlgebraTensorModule.rid_tmul, Algebra.smul_def, mul_one]

/-- `ker(mapBaseChange)` を `tensorPolynomialOmegaEquiv` で `B` に転送すると、
`f'(root)`(= δ)が生成するイデアルにちょうど一致する。`range_kerCotangentToTensor_span`
(値域が `{1⊗ₜD(f)}` の A-張る空間)と `tensorPolynomialOmegaEquiv_D`
(その像が δ)を、A→B が全射という事実で B-張る空間に架け替える。 -/
theorem map_ker_mapBaseChange_eq_span {R B : Type*} [CommRing R] [CommRing B]
    [Algebra (Polynomial R) B] [Algebra R B] [IsScalarTower R (Polynomial R) B] (f : Polynomial R)
    (hB : RingHom.ker (algebraMap (Polynomial R) B) = Ideal.span ({f} : Set (Polynomial R)))
    (hsurj : Function.Surjective (algebraMap (Polynomial R) B)) :
    Submodule.map (tensorPolynomialOmegaEquiv (R:=R) (B:=B)).toLinearMap
      (LinearMap.ker (KaehlerDifferential.mapBaseChange R (Polynomial R) B))
      = Ideal.span ({algebraMap (Polynomial R) B (Polynomial.derivative f)} : Set B) := by
  set δ := algebraMap (Polynomial R) B (Polynomial.derivative f) with hδ
  set e := (tensorPolynomialOmegaEquiv (R:=R) (B:=B))
  have hsmul : ∀ (a : Polynomial R) (x : TensorProduct (Polynomial R) B Ω[Polynomial R⁄R]),
      e (a • x) = a • e x := fun a x => LinearMap.map_smul_of_tower (e.toLinearMap) a x
  apply SetLike.ext'
  show (e : TensorProduct (Polynomial R) B Ω[Polynomial R⁄R] → B) '' _ = _
  rw [ker_mapBaseChange_eq_range_kerCotangentToTensor hsurj]
  have hrange : Set.range (KaehlerDifferential.kerCotangentToTensor R (Polynomial R) B)
      = (Submodule.span (Polynomial R)
          {(1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f)} : Set _) := by
    rw [← range_kerCotangentToTensor_span f hB]; rfl
  rw [hrange]
  ext b
  simp only [Set.mem_image, SetLike.mem_coe, Submodule.mem_span_singleton]
  constructor
  · rintro ⟨y, ⟨a, rfl⟩, rfl⟩
    refine ⟨algebraMap (Polynomial R) B a, ?_⟩
    rw [hsmul, tensorPolynomialOmegaEquiv_D, hδ]
    simp [Algebra.smul_def]
  · rintro ⟨c, rfl⟩
    obtain ⟨a, rfl⟩ := hsurj c
    refine ⟨a • ((1:B) ⊗ₜ[Polynomial R] (KaehlerDifferential.D R (Polynomial R) f)), ⟨a, rfl⟩, ?_⟩
    rw [hsmul, tensorPolynomialOmegaEquiv_D, hδ]
    simp [Algebra.smul_def]

/-- **[Falt1] Lemma 1.1 の "第二完全列" 側の核心計算(完成)**: `B = R[X] ⧸ (f)`
のとき、`Ω_{B/R} ≅ B ⧸ (f'(root))`。これは Faltings の原文の議論
「`Ω_W` は `Ω_V⊗W⊕WdT` を `f'(w)dT` で割った商」の mathlib 版。

★これ自体は Lemma 1.1 の**主張そのもの**ではない——Falt1 の `V,W` は
一般の完備離散付値環の拡大で、ここでの `R,B` はその中でも
`B = R[X]/(f)` という「モノジェニックな塔の1段」の場合。Lemma 1.1 は
`V→W` の塔(降下引数で結局モノジェニックな場合に帰着するが、その帰着
自体は Faltings の議論の一部で、まだ形式化していない)に対する主張。
それでも本定理は Lemma 1.1 の証明で実際に使われる計算そのものであり、
「原文の証明の心臓部を実際に検証した」という意味で重要な一歩。 -/
noncomputable def omega_quotient_eq_derivative_span {R B : Type*} [CommRing R] [CommRing B]
    [Algebra (Polynomial R) B] [Algebra R B] [IsScalarTower R (Polynomial R) B] (f : Polynomial R)
    (hB : RingHom.ker (algebraMap (Polynomial R) B) = Ideal.span ({f} : Set (Polynomial R)))
    (hsurj : Function.Surjective (algebraMap (Polynomial R) B)) :
    Ω[B⁄R] ≃ₗ[B]
      B ⧸ Ideal.span ({algebraMap (Polynomial R) B (Polynomial.derivative f)} : Set B) :=
  (kerQuotEquivOmega hsurj).symm.trans
    (Submodule.Quotient.equiv _ _ (tensorPolynomialOmegaEquiv (R:=R) (B:=B))
      (map_ker_mapBaseChange_eq_span f hB hsurj))

/-- **differentIdeal への接続(完成)**: `w` が `W` を `V`-代数として生成し、
かつ(拡大体レベルでも)`L` を `K`-代数として生成するなら、
`differentIdeal V W = span{f'(w)}`(`f = minpoly V w`)。
`conductor_mul_differentIdeal`(mathlib)と、生成元による conductor の
自明化(`conductor_eq_top_of_adjoin_eq_top`)を組み合わせる。

これで Lemma 1.1 の「長さ」の主張(`Ω_{W/V}` の長さ `= length(W/p^δW)`、
`p^δ` は differentIdeal の生成元)への接続が完成した:
`omega_quotient_eq_derivative_span` と本定理を貼り合わせれば
`Ω_{W/V} ≅ W ⧸ differentIdeal V W` が得られる(まだ貼り合わせていない
——`omega_quotient_eq_derivative_span` は `Polynomial R ⧸ (f)` という
具体形を使うが、本定理は `AdjoinRoot`/一般の `W` で述べており、両者を
繋ぐには `W ≃ Polynomial V ⧸ (minpoly V w)` という同型がさらに要る)。 -/
theorem differentIdeal_eq_span_derivative {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤)
    (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤) :
    differentIdeal V W = Ideal.span {(Polynomial.aeval w) (Polynomial.derivative (minpoly V w))} := by
  have hc : conductor V w = ⊤ := conductor_eq_top_of_adjoin_eq_top hadjoin
  have := conductor_mul_differentIdeal V K L w hw
  rw [hc, Ideal.top_mul] at this
  exact this

/-- **一般の W への橋渡し(1本目)**: `V` が整閉整域(例: Dedekind整域)で
`w : W` が `V` 上整かつ `W` を生成するなら、`AdjoinRoot (minpoly V w) ≃ₐ[V] W`。
`IsAdjoinRootMonic.mkOfAdjoinEqTop`(体に限らず一般の整閉整域 `V` で成立、
`AlgEquiv.adjoinSingletonEquivAdjoinRootMinpoly` は `V` が体の場合のみ
なので使えない)経由。

`AdjoinRoot f` は定義から `Polynomial V ⧸ (f)` なので、これで
`omega_quotient_eq_derivative_span`(`Polynomial R ⧸ (f)` の具体形)と
Falt1 の一般の `W` を繋ぐ最初の橋になる。②の橋渡しは `omegaCongr`
(下記)で完成した。 -/
noncomputable def adjoinRootMinpolyEquiv {V W : Type*} [CommRing V] [CommRing W] [Algebra V W]
    [IsDomain V] [IsDomain W] [Module.IsTorsionFree V W] [IsIntegrallyClosed V]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤) :
    AdjoinRoot (minpoly V w) ≃ₐ[V] W :=
  IsAdjoinRoot.adjoinRootAlgEquiv (IsAdjoinRootMonic.mkOfAdjoinEqTop hint hadjoin).toIsAdjoinRoot

/-- `e : A ≃ₐ[R] B` から `Algebra A B`(`e` 自身を経由)を作ると
`IsScalarTower R A B` が成り立つ(`e` が `R`-線形であることから)。
`omegaCongr` の下準備。 -/
theorem isScalarTower_of_algEquiv {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra R B] (e : A ≃ₐ[R] B) :
    letI : Algebra A B := e.toRingHom.toAlgebra
    IsScalarTower R A B := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  constructor
  intro r a b
  show algebraMap A B (r • a) * b = r • (algebraMap A B a * b)
  have h : algebraMap A B (r • a) = r • algebraMap A B a := by
    show e (r • a) = r • e a
    exact map_smul e r a
  rw [h, smul_mul_assoc]

/-- `omegaCongr` の左逆性(`map R R B A ∘ map R R A B = id`)。生成元の集合
`Set.range (D R A)` 上で `map_D` により具体的に計算し(`e.symm(e a)=a`)、
`Submodule.span_induction` で加法・A-スカラー倍に持ち上げる。 -/
theorem omegaCongr_leftInv {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra R B] (e : A ≃ₐ[R] B) :
    letI : Algebra A B := e.toRingHom.toAlgebra
    letI : IsScalarTower R A B := isScalarTower_of_algEquiv e
    letI : Algebra B A := e.symm.toRingHom.toAlgebra
    letI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
    ∀ x : Ω[A⁄R], (KaehlerDifferential.map R R B A) ((KaehlerDifferential.map R R A B) x) = x := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower R A B := isScalarTower_of_algEquiv e
  letI : Algebra B A := e.symm.toRingHom.toAlgebra
  haveI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
  have hspan := KaehlerDifferential.span_range_derivation (R := R) (S := A)
  intro x
  have hx : x ∈ Submodule.span A (Set.range (KaehlerDifferential.D R A)) := hspan ▸ Submodule.mem_top
  induction hx using Submodule.span_induction with
  | mem y hy =>
      obtain ⟨a, rfl⟩ := hy
      rw [KaehlerDifferential.map_D, KaehlerDifferential.map_D]
      congr 1
      show e.symm (e a) = a
      exact e.symm_apply_apply a
  | zero => simp
  | add y z hy hz ihy ihz => rw [map_add, map_add, ihy, ihz]
  | smul a y hy ih =>
      show (KaehlerDifferential.map R R B A) ((KaehlerDifferential.map R R A B) (a • y)) = a • y
      rw [map_smul]
      show (KaehlerDifferential.map R R B A)
        (algebraMap A B a • (KaehlerDifferential.map R R A B) y) = a • y
      rw [map_smul]
      show algebraMap B A (algebraMap A B a) •
        (KaehlerDifferential.map R R B A) ((KaehlerDifferential.map R R A B) y) = a • y
      rw [ih]
      congr 1
      show e.symm (e a) = a
      exact e.symm_apply_apply a

/-- `omegaCongr` の右逆性(対称な議論)。 -/
theorem omegaCongr_rightInv {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra R B] (e : A ≃ₐ[R] B) :
    letI : Algebra A B := e.toRingHom.toAlgebra
    letI : IsScalarTower R A B := isScalarTower_of_algEquiv e
    letI : Algebra B A := e.symm.toRingHom.toAlgebra
    letI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
    ∀ x : Ω[B⁄R], (KaehlerDifferential.map R R A B) ((KaehlerDifferential.map R R B A) x) = x := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower R A B := isScalarTower_of_algEquiv e
  letI : Algebra B A := e.symm.toRingHom.toAlgebra
  haveI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
  have hspan := KaehlerDifferential.span_range_derivation (R := R) (S := B)
  intro x
  have hx : x ∈ Submodule.span B (Set.range (KaehlerDifferential.D R B)) := hspan ▸ Submodule.mem_top
  induction hx using Submodule.span_induction with
  | mem y hy =>
      obtain ⟨b, rfl⟩ := hy
      rw [KaehlerDifferential.map_D, KaehlerDifferential.map_D]
      congr 1
      show e (e.symm b) = b
      exact e.apply_symm_apply b
  | zero => simp
  | add y z hy hz ihy ihz => rw [map_add, map_add, ihy, ihz]
  | smul b y hy ih =>
      show (KaehlerDifferential.map R R A B) ((KaehlerDifferential.map R R B A) (b • y)) = b • y
      rw [map_smul]
      show (KaehlerDifferential.map R R A B)
        (algebraMap B A b • (KaehlerDifferential.map R R B A) y) = b • y
      rw [map_smul]
      show algebraMap A B (algebraMap B A b) •
        (KaehlerDifferential.map R R A B) ((KaehlerDifferential.map R R B A) y) = b • y
      rw [ih]
      congr 1
      show e (e.symm b) = b
      exact e.apply_symm_apply b

/-- **②の橋渡し(完成)**: 環同型 `e : A ≃ₐ[R] B` は `Ω[A⁄R] ≃+ Ω[B⁄R]`
(加法群としての同型)を誘導する。mathlib には代数の**塔**
(`KaehlerDifferential.map`)用の道具はあるが、代数**同型**に沿った
transport の既製品は無かった——`map R R A B` と `map R R B A` を
(`e`・`e.symm` から作った局所インスタンス `Algebra A B`・`Algebra B A`
経由で)組み立て、`omegaCongr_leftInv`/`_rightInv` で互いに逆写像である
ことを示して構成した。★加法群としての同型に留めた(`≃+`、`A`-線形性は
求めない)——長さ・濃度の比較には十分。 -/
noncomputable def omegaCongr {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra R B] (e : A ≃ₐ[R] B) :
    Ω[A⁄R] ≃+ Ω[B⁄R] := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower R A B := isScalarTower_of_algEquiv e
  letI : Algebra B A := e.symm.toRingHom.toAlgebra
  haveI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
  exact AddEquiv.ofBijective (KaehlerDifferential.map R R A B).toAddMonoidHom
    ⟨Function.LeftInverse.injective (omegaCongr_leftInv e),
     fun y => ⟨KaehlerDifferential.map R R B A y, omegaCongr_rightInv e y⟩⟩

/-! ## ①②③④ の貼り合わせ(完成) -/

/-- `AdjoinRoot f`(`= Polynomial R ⧸ (f)`)には自動では見つからない
`Algebra (Polynomial R) (AdjoinRoot f)` を明示的に登記する(`AdjoinRoot`
が `def` で simp-reducible ではないため、instance 探索が定義を展開して
既存の `Ideal.Quotient.algebra` を見つけられない——2026-09-04 に発見)。 -/
noncomputable instance AdjoinRoot.instAlgebraPolynomial {R : Type*} [CommRing R] (f : Polynomial R) :
    Algebra (Polynomial R) (AdjoinRoot f) :=
  inferInstanceAs (Algebra (Polynomial R) (Polynomial R ⧸ Ideal.span ({f} : Set (Polynomial R))))

instance AdjoinRoot.instIsScalarTowerPolynomial {R : Type*} [CommRing R] (f : Polynomial R) :
    IsScalarTower R (Polynomial R) (AdjoinRoot f) :=
  inferInstanceAs (IsScalarTower R (Polynomial R) (Polynomial R ⧸ Ideal.span ({f} : Set (Polynomial R))))

/-- ①(`omega_quotient_eq_derivative_span`)を `AdjoinRoot f` に特殊化し、
`≃+`(加法群としての同型)に弱めたもの。 -/
noncomputable def omegaAdjoinRootQuot {R : Type*} [CommRing R] (f : Polynomial R) :
    Ω[(AdjoinRoot f)⁄R] ≃+
      (AdjoinRoot f) ⧸ Ideal.span
        ({algebraMap (Polynomial R) (AdjoinRoot f) (Polynomial.derivative f)} : Set (AdjoinRoot f)) :=
  (omega_quotient_eq_derivative_span f (ker_adjoinRoot_mk f) AdjoinRoot.mk_surjective).toAddEquiv

/-- `IsAdjoinRootMonic.mkOfAdjoinEqTop hint hadjoin` の根は `w` そのもの。 -/
theorem adjoinRootMinpolyEquiv_root {V W : Type*} [CommRing V] [CommRing W] [Algebra V W]
    [IsDomain V] [IsDomain W] [Module.IsTorsionFree V W] [IsIntegrallyClosed V]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤) :
    (IsAdjoinRootMonic.mkOfAdjoinEqTop hint hadjoin).toIsAdjoinRoot.root = w :=
  IsAdjoinRoot.mkOfAdjoinEqTop_root (hα := hint) (hα₂ := hadjoin)

/-- `adjoinRootMinpolyEquiv` は `Polynomial V` からの `algebraMap` の像を
`w` での evaluation(`aeval`)に写す——`e` の構成(`IsAdjoinRoot.map` =
`aeval` at `root`、`IsAdjoinRoot.aeval_root_eq_map`)と
`adjoinRootMinpolyEquiv_root`(根が `w`)から従う。 -/
theorem adjoinRootMinpolyEquiv_algebraMap {V W : Type*} [CommRing V] [CommRing W] [Algebra V W]
    [IsDomain V] [IsDomain W] [Module.IsTorsionFree V W] [IsIntegrallyClosed V]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤) (p : Polynomial V) :
    (adjoinRootMinpolyEquiv w hint hadjoin) (algebraMap (Polynomial V) (AdjoinRoot (minpoly V w)) p)
      = Polynomial.aeval w p := by
  unfold adjoinRootMinpolyEquiv
  rw [show algebraMap (Polynomial V) (AdjoinRoot (minpoly V w)) p = AdjoinRoot.mk (minpoly V w) p from rfl,
    IsAdjoinRoot.adjoinRootAlgEquiv_apply_mk]
  conv_rhs => rw [← adjoinRootMinpolyEquiv_root w hint hadjoin]
  exact (congrFun (congrArg _
    (IsAdjoinRoot.aeval_root_eq_map (h := (IsAdjoinRootMonic.mkOfAdjoinEqTop hint hadjoin).toIsAdjoinRoot))) p).symm

/-- **[Falt1] Lemma 1.1 の余核の主張(完成)**: `V` が Dedekind 整域、
`w : W` が `V` 上整・`W` を生成し(かつ拡大体レベルでも生成する)なら、
`Ω_{W/V} ≅ W ⧸ differentIdeal V W`(`≃+`、加法群としての同型)。

これで Lemma 1.1「`Ω_{W/V}` の長さは `length(W/p^δW)` に等しい」の
**余核側**(`Ω_{W/V}` が何であるかの特定)が完全に証明された——`p^δ` は
`differentIdeal V W` の生成元(δ は Faltings の記法での different の
指数)。①②③④(`omega_quotient_eq_derivative_span`・
`differentIdeal_eq_span_derivative`・`adjoinRootMinpolyEquiv`・
`omegaCongr`)を `adjoinRootMinpolyEquiv_algebraMap`(橋渡しの自然性)・
`Ideal.map_span`・`Ideal.quotientEquiv` で貼り合わせた。

★まだ残っている: (a) 「長さ」概念そのもの(`length(W/p^δW)`)を
mathlib のどの道具で表すか(本定理は「同型」を与えるので、あとは
「同型な加法群は同じ長さを持つ」という一般論を適用するだけのはず)、
(b) Lemma 1.1 の**単射性**の主張(`Ω_V⊗_VW→Ω_W` の単射性、絶対微分の
「第一完全列」——本定理とは独立な、別の議論)。 -/
noncomputable def falt1CokernelIso {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Ω[W⁄V] ≃+ (W ⧸ differentIdeal V W) := by
  set e := adjoinRootMinpolyEquiv w hint hadjoin
  set step1 := (omegaCongr e).symm.trans (omegaAdjoinRootQuot (minpoly V w))
  set J := Ideal.map (e.toRingEquiv : AdjoinRoot (minpoly V w) →+* W)
    (Ideal.span ({algebraMap (Polynomial V) (AdjoinRoot (minpoly V w))
      (Polynomial.derivative (minpoly V w))} : Set (AdjoinRoot (minpoly V w))))
  set step2 := Ideal.quotientEquiv _ J e.toRingEquiv rfl
  have hJ : J = differentIdeal V W := by
    show Ideal.map _ (Ideal.span _) = _
    rw [Ideal.map_span, Set.image_singleton]
    have heval : (e.toRingEquiv : AdjoinRoot (minpoly V w) →+* W)
        (algebraMap (Polynomial V) (AdjoinRoot (minpoly V w)) (Polynomial.derivative (minpoly V w)))
        = Polynomial.aeval w (Polynomial.derivative (minpoly V w)) :=
      adjoinRootMinpolyEquiv_algebraMap w hint hadjoin (Polynomial.derivative (minpoly V w))
    rw [heval, ← differentIdeal_eq_span_derivative w hw hadjoin]
  exact step1.trans (hJ ▸ step2.toAddEquiv)

end ABC3.Found.Falt1
