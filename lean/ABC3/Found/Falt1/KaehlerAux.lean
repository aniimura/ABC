import Mathlib.RingTheory.Kaehler.Basic
import Mathlib.RingTheory.Kaehler.Polynomial
import Mathlib.RingTheory.Kaehler.JacobiZariski
import Mathlib.LinearAlgebra.Span.Basic
import Mathlib.RingTheory.AdjoinRoot
import Mathlib.RingTheory.DedekindDomain.Different
import Mathlib.RingTheory.IsAdjoinRoot
import Mathlib.RingTheory.Smooth.Basic
import Mathlib.RingTheory.Kaehler.TensorProduct
import Mathlib.LinearAlgebra.TensorProduct.Prod
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# [Falt1] Chapter I §1 の補助補題群(`Found`、sorry 無し)

★★★これは `/goal Falt1 Chapter I Found` の一部。**Lemma 1.1(§1
1/2)は `Found/Falt1/Lemma11.lean`(`lemma_1_1_falt1`・
`falt1RamificationSetup`)で数学的内容+Interface統合ともに完成した**
——本ファイルはその土台(①〜⑤、下記)と、**Theorem 1.2(§1 2/2、
まだ未完成)** に向けた新たな補助補題群(`pushoutKaehlerSplit` 等)を
両方含む(正直な進捗として記録する)。

## Lemma 1.1 の土台(2026-09-04、すべて sorry 無し。完成)

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

**⑤ `≃ₗ[W]` への格上げと「長さ」への接続(完成)**: `falt1CokernelIso`
は`≃+`(加法群としての同型)だったが、`falt1CokernelIsoLinear` で
`Ω_{W/V} ≃ₗ[W] (W ⧸ differentIdeal V W)`(本物の `W`-線形同型)まで
証明した。3段の合成(`(omegaCongr e).symm`・`omegaAdjoinRootQuot`・
`Ideal.quotientEquiv`)それぞれのスムル整合性を辿って
`AddEquiv.toLinearEquiv` に渡した(鍵は `Ideal.quotientEquiv_mk` の
naturality——`Ideal.quotientEquivAlg` の重い `AlgEquiv→+*` 強制を
一切使わずに済んだ)。これで `LinearEquiv.length_eq` が直接使え、
`falt1CokernelLengthEq : Module.length W Ω[W⁄V] =
Module.length W (W ⧸ differentIdeal V W)` を得た——Lemma 1.1 の
「長さ」の主張そのものが完成した。

## Lemma 1.1(§1 1/2)は完成した(2026-09-04)

上の①〜⑤(余核の特定・differentIdeal への接続・一般の `W` への橋渡し・
`Ω` の代数同型に沿った転送・`≃ₗ[W]` への格上げ)に加え、本ファイル後半
で**単射性**(`falt1MapBaseChangeInjective`、`mapBaseChange_injective_
of_nzd`・`mapBaseChange_injective_transport` 等を経由)と
**`differentIdeal V W ≠ ⊥` の導出**(`falt1_differentIdeal_ne_bot`、
`FractionRing.liftAlgebra` 経由)も完成し、`Found/Falt1/Lemma11.lean`
で **`Interface.RamificationSetup` への正式な差し替え**
(`falt1RamificationSetup`、`Module.length` の `ℕ∞→ℕ` 有限性込み)まで
仕上げた。`lemma_1_1 (falt1RamificationSetup w hint hadjoin hw)` が
型検査を通ることを確認済み——Lemma 1.1 は数学的内容・Interface統合
ともに完成している。詳細な逸脱の記録は `Lemma11.lean` を参照。

## Theorem 1.2(§1 2/2)へ向けた新たな補助補題(2026-09-04、未完成)

原文(Falt1 p.5、260dpi 目視で読み直し済み、`ResearchPaper/falt1-goal.md`
参照)は「非常に分岐した拡大の列 `V_n`」に対する differentIdeal の
収束を Lemma 1.1 の繰り返し適用で示す——その核心である**複数生成元版
Lemma 1.1**(`Ω_{W_{n+1}/V_n}` が `d+1` 個の直和、という主張)の
**2生成元の場合**を `pushoutKaehlerSplit` として完成した(下記の節)。
★残っている作業(`ResearchPaper/falt1-goal.md` に詳細): `d+1` 個への
一般化・`mapBaseChange` 単射性の具体的な条件・「非常に分岐した」`V_n`
の族そのものの形式化・`Module.length` の漸化不等式。§2-4 は Faltings
の "almost mathematics" 自体が mathlib に無い(`Almost`/`Malcev`/
`Tannakian` を mathlib 全体検索して確認済み)ため、ライブラリ新設に
近い規模の作業になる。
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

/-- `(omegaCongr e).symm` は具体的に `map R R B A` そのもの(一意性から)。
`map R R B A` は `B`-線形(=元の `omegaCongr` の定義で `B := W` のとき
`W`-線形)なので、`falt1CokernelIso` を `≃ₗ[W]` へ持ち上げる際の第一段の
スムル整合性はこれで手に入る(`map_smul` を直接使えばよい)。 -/
theorem omegaCongr_symm_eq {R A B : Type*} [CommRing R] [CommRing A] [CommRing B]
    [Algebra R A] [Algebra R B] (e : A ≃ₐ[R] B) :
    letI : Algebra A B := e.toRingHom.toAlgebra
    letI : IsScalarTower R A B := isScalarTower_of_algEquiv e
    letI : Algebra B A := e.symm.toRingHom.toAlgebra
    letI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
    ∀ x, (omegaCongr e).symm x = KaehlerDifferential.map R R B A x := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower R A B := isScalarTower_of_algEquiv e
  letI : Algebra B A := e.symm.toRingHom.toAlgebra
  haveI : IsScalarTower R B A := isScalarTower_of_algEquiv e.symm
  intro x
  apply (omegaCongr e).injective
  rw [AddEquiv.apply_symm_apply]
  show x = (KaehlerDifferential.map R R A B) ((KaehlerDifferential.map R R B A) x)
  exact (omegaCongr_rightInv e x).symm

/-!
## `falt1CokernelIso` の `≃ₗ[W]` への格上げ(完成、2026-09-04)

★★**性能上の壁を回避できた実測**: 当初は `step2 = Ideal.quotientEquiv
...` の段を `Ideal.quotientEquivAlg` に差し替える計画だったが、
`Ideal.map (e : AdjoinRoot (minpoly V w) →+* W) ...` の型検査が
`maxHeartbeats 1000000`(既定の5倍、44秒)でも `timeout at whnf` で
止まった(`AlgEquiv` の `→+*` への強制が `RingEquiv` 版よりはるかに
重いため)。代わりに **`Ideal.quotientEquiv` 自身の naturality**
(`Ideal.quotientEquiv_mk : (I.quotientEquiv J f hIJ) (mk I x) = mk J (f x)`
——mathlib に既製品として存在、`f : R ≃+* S` という軽い `RingEquiv`
だけを使う)を直接使うことで、`AlgEquiv` 強制を一切経由せずに
`step2` のスムル整合性を 0.2 秒で証明できた。
-/

/-- **`falt1CokernelIso` の `≃ₗ[W]` 版(完成)**: `Ω_{W/V} ≃ₗ[W]
W ⧸ differentIdeal V W`——`falt1CokernelIso` と同じ内容だが本物の
`W`-線形同型として構成し直した。3段の合成それぞれのスムル整合性:

1. `(omegaCongr e).symm` の段: `omegaCongr_symm_eq` により
   `map V V W (AdjoinRoot ...)` そのもの——`map_smul` で `W`-線形性。
2. `omegaAdjoinRootQuot (minpoly V w)` の段: `omega_quotient_eq_derivative_span`
   は本物の `LinearEquiv`(`AdjoinRoot(minpoly V w)`-線形)なので
   `map_smul` がそのまま使える。
3. `step2 = Ideal.quotientEquiv ...` の段: `Ideal.quotientEquiv_mk`
   (naturality)を `Ideal.Quotient.mk_surjective` で商の代表元に
   帰着させて適用(`AlgEquiv→+*` 強制を経由しない、上のノート参照)。

3段を貼り合わせ、`AddEquiv.toLinearEquiv` にスムル整合性の証明を渡した。 -/
noncomputable def falt1CokernelIsoLinear {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Ω[W⁄V] ≃ₗ[W] (W ⧸ differentIdeal V W) := by
  set e := adjoinRootMinpolyEquiv w hint hadjoin with he
  letI : Algebra (AdjoinRoot (minpoly V w)) W := e.toRingHom.toAlgebra
  haveI : IsScalarTower V (AdjoinRoot (minpoly V w)) W := isScalarTower_of_algEquiv e
  letI : Algebra W (AdjoinRoot (minpoly V w)) := e.symm.toRingHom.toAlgebra
  haveI : IsScalarTower V W (AdjoinRoot (minpoly V w)) := isScalarTower_of_algEquiv e.symm
  set J := Ideal.map (e.toRingEquiv : AdjoinRoot (minpoly V w) →+* W)
    (Ideal.span ({algebraMap (Polynomial V) (AdjoinRoot (minpoly V w))
      (Polynomial.derivative (minpoly V w))} : Set (AdjoinRoot (minpoly V w)))) with hJdef
  set step1 := (omegaCongr e).symm.trans (omegaAdjoinRootQuot (minpoly V w)) with hstep1
  set step2 := Ideal.quotientEquiv _ J e.toRingEquiv rfl with hstep2
  have hJ : J = differentIdeal V W := by
    show Ideal.map _ (Ideal.span _) = _
    rw [Ideal.map_span, Set.image_singleton]
    have heval : (e.toRingEquiv : AdjoinRoot (minpoly V w) →+* W)
        (algebraMap (Polynomial V) (AdjoinRoot (minpoly V w)) (Polynomial.derivative (minpoly V w)))
        = Polynomial.aeval w (Polynomial.derivative (minpoly V w)) :=
      adjoinRootMinpolyEquiv_algebraMap w hint hadjoin (Polynomial.derivative (minpoly V w))
    rw [heval, ← differentIdeal_eq_span_derivative w hw hadjoin]
  have hsmul : ∀ (c : W) (x : Ω[W⁄V]), (step1.trans step2.toAddEquiv) (c • x) =
      c • (step1.trans step2.toAddEquiv) x := by
    intro c x
    show step2 (step1 (c • x)) = c • step2 (step1 x)
    have hA : step1 (c • x) = (algebraMap W (AdjoinRoot (minpoly V w)) c) • step1 x := by
      show (omegaAdjoinRootQuot (minpoly V w)) ((omegaCongr e).symm (c • x)) =
        (algebraMap W (AdjoinRoot (minpoly V w)) c) • (omegaAdjoinRootQuot (minpoly V w)) ((omegaCongr e).symm x)
      have h1cx := omegaCongr_symm_eq e (c • x)
      have h1x := omegaCongr_symm_eq e x
      rw [h1cx, h1x]
      have h2 := map_smul (KaehlerDifferential.map V V W (AdjoinRoot (minpoly V w))) c x
      rw [h2]
      show (omega_quotient_eq_derivative_span (minpoly V w) (ker_adjoinRoot_mk (minpoly V w))
              AdjoinRoot.mk_surjective).toAddEquiv
            ((algebraMap W (AdjoinRoot (minpoly V w)) c) •
              (KaehlerDifferential.map V V W (AdjoinRoot (minpoly V w))) x) =
          (algebraMap W (AdjoinRoot (minpoly V w)) c) •
            (omega_quotient_eq_derivative_span (minpoly V w) (ker_adjoinRoot_mk (minpoly V w))
              AdjoinRoot.mk_surjective).toAddEquiv ((KaehlerDifferential.map V V W (AdjoinRoot (minpoly V w))) x)
      exact map_smul (omega_quotient_eq_derivative_span (minpoly V w) (ker_adjoinRoot_mk (minpoly V w))
        AdjoinRoot.mk_surjective) (algebraMap W (AdjoinRoot (minpoly V w)) c)
        ((KaehlerDifferential.map V V W (AdjoinRoot (minpoly V w))) x)
    rw [hA]
    show step2 ((algebraMap W (AdjoinRoot (minpoly V w)) c) • step1 x) = c • step2 (step1 x)
    have hC : ∀ (a : AdjoinRoot (minpoly V w)) (z : AdjoinRoot (minpoly V w) ⧸
        Ideal.span ({algebraMap (Polynomial V) (AdjoinRoot (minpoly V w))
          (Polynomial.derivative (minpoly V w))} : Set (AdjoinRoot (minpoly V w)))),
        step2 (a • z) = (e.toRingEquiv a) • step2 z := by
      intro a z
      obtain ⟨z', rfl⟩ := Ideal.Quotient.mk_surjective z
      show step2 (Ideal.Quotient.mk _ (a * z')) = e.toRingEquiv a • step2 (Ideal.Quotient.mk _ z')
      rw [hstep2, Ideal.quotientEquiv_mk, Ideal.quotientEquiv_mk, map_mul]
      rfl
    rw [hC]
    congr 1
    show e.toRingEquiv (algebraMap W (AdjoinRoot (minpoly V w)) c) = c
    show e (e.symm c) = c
    exact e.apply_symm_apply c
  exact hJ ▸ AddEquiv.toLinearEquiv (step1.trans step2.toAddEquiv) hsmul

/-- **[Falt1] Lemma 1.1「長さ」の主張(完成)**: `Module.length` の
`LinearEquiv` に対する不変性(`LinearEquiv.length_eq`)を
`falt1CokernelIsoLinear` に適用しただけ——`Ω_{W/V}` の長さが
`W ⧸ differentIdeal V W` の長さ(Faltings の記法で `length(W/p^δW)`)
に等しいという、Lemma 1.1 の主張そのもの。 -/
theorem falt1CokernelLengthEq {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Module.length W Ω[W⁄V] = Module.length W (W ⧸ differentIdeal V W) :=
  LinearEquiv.length_eq (falt1CokernelIsoLinear w hint hadjoin hw)

/-!
## Lemma 1.1「単射性」への一歩(第一完全列の直和分解、完成、2026-09-04)

Faltings の議論は塔 `Z → V → V[T] → W(=V[T]/(f))` を使い、`Ω_{V[T]/Z} ≅
(V[T]⊗_VΩ_{V/Z}) ⊕ V[T]dT` という直和分解と `f'(w)` が非零因子であること
から `Ω_V⊗_VW → Ω_W` の単射性を導く。★★この直和分解そのもの
(`polynomialKaehlerSplit`)を今回 mathlib の道具だけで組み立てた——
mathlib には**既製品が無かった**(`polynomialEquiv`・`mvPolynomialBasis`は
`R[X]/R` の自明底のみで、2段の塔 `Z→V→V[T]` は扱わない)。

鍵になった発見:
1. `Function.Exact.splitSurjectiveEquiv`(`Algebra/Exact/Basic.lean`)——
   完全列 `f,g` で `f` 単射・`g` にセクション `l` があれば
   `N ≃ₗ[R] M × P` を直接与える既製品(探すまで知らなかった)。
2. **`mapBaseChange` の単射性(一般の塔 `Z→V→V[T]` で常に成立)**:
   `Algebra.H1Cotangent.exact_δ_mapBaseChange`(Jacobi-Zariski 完全列の
   境界写像 `δ`)+ `Subsingleton (H1Cotangent V (Polynomial V))` から従う。
   後者は mathlib に無かったので `subsingleton_H1Cotangent_self` として
   一般に証明した(**`FormallySmooth R S` ⟹ `Subsingleton (H1Cotangent R S)`**
   ——`S` 自身を「関係式 0 個の自明な presentation」`Extension.ofSurjective
   (AlgHom.id R S)` とみなし、その `Cotangent` が `(⊥:Ideal S).Cotangent`
   に一致して自明であることと `equivH1CotangentOfFormallySmooth` を貼り
   合わせた)。`Polynomial V` は `FormallySmooth V (Polynomial V)`
   (mathlib に既製品`Algebra.Extension.Algebra.FormallySmooth.polynomial`)
   なので直ちに使える。
3. **`g`(`Ω_{V[T]/Z}→Ω_{V[T]/V}`)のセクション**: `polynomialEquiv V :
   Ω_{V[T]/V}≃ₗ[V[T]]V[T]` と `D_Z(T)` から `y ↦ (polynomialEquiv V y)•D_Z(T)`
   という明示的な `V[T]`-線形写像を作り、`map_D`(自然性)で
   セクションであることを直接計算で確認した。

★これは Lemma 1.1 の単射性の**最難関の前提**(直和分解そのもの)だが、
Lemma 1.1 本体にはまだ届いていない——残る作業(未着手):
(a) `f'(w)` が `D_Z(f)` の「`dT` 成分」であることの明示計算
    (`polynomialKaehlerSplit` の下で `D_Z(f)` を書き下す)、
(b) `W` への base change(`V[T]⊗_{V[T]}W = W` 経由で直和分解を `W` に
    移す)、(c) 「自由な直和成分に非零因子で割ると他の成分の単射性は
    保たれる」という初等的だが未証明の議論。これらは (a)(b)(c) の
    3ステップとして構造は見えているが、まだ1行も書いていない。
-/

theorem subsingleton_H1Cotangent_self {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.FormallySmooth R S] : Subsingleton (Algebra.H1Cotangent R S) := by
  have h1 : Subsingleton (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id).H1Cotangent := by
    have hcot : Subsingleton (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id).Cotangent := by
      show Subsingleton (RingHom.ker (algebraMap
        (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id).Ring S)).Cotangent
      have hker : RingHom.ker (algebraMap
          (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id).Ring S) = ⊥ := by
        show RingHom.ker (algebraMap S S) = ⊥
        simp [RingHom.ker_eq_bot_iff_eq_zero]
      rw [hker]
      have hsurj := (⊥ : Ideal S).toCotangent_surjective
      constructor
      intro a b
      obtain ⟨a', rfl⟩ := hsurj a
      obtain ⟨b', rfl⟩ := hsurj b
      congr 1
      exact Subsingleton.elim a' b'
    exact Subsingleton.intro (fun a b => Subtype.ext (Subsingleton.elim a.1 b.1))
  haveI : Algebra.FormallySmooth R
      (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id).Ring := by
    show Algebra.FormallySmooth R S; infer_instance
  exact (Algebra.Extension.equivH1CotangentOfFormallySmooth
    (Algebra.Extension.ofSurjective (AlgHom.id R S) Function.surjective_id)).symm.injective.subsingleton

/-- 塔 `Z → V → V[T]` に対し、`mapBaseChange`(`V[T]⊗_VΩ_{V/Z}→Ω_{V[T]/Z}`)は
常に単射(`V→V[T]` が多項式拡大であることから、`H1Cotangent V (Polynomial V)`
が自明になるため)。 -/
theorem mapBaseChange_injective_polynomial {Z V : Type*} [CommRing Z] [CommRing V] [Algebra Z V] :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V (Polynomial V)) := by
  haveI := subsingleton_H1Cotangent_self (R := V) (S := Polynomial V)
  rw [← LinearMap.ker_eq_bot]
  have hex := Algebra.H1Cotangent.exact_δ_mapBaseChange Z V (Polynomial V)
  rw [LinearMap.exact_iff] at hex
  rw [hex, eq_bot_iff]
  rintro y ⟨x, rfl⟩
  simp [Subsingleton.elim x 0]

/-- **Falt1 Lemma 1.1 の「第一完全列」直和分解(完成)**: `Ω_{V[T]/Z} ≅
(V[T]⊗_V Ω_{V/Z}) × V[T]`(`dT` 方向を素朴な直積の第二成分として表す)。
`mapBaseChange_injective_polynomial`(単射性)と、`polynomialEquiv` から
作った明示的なセクション(`Function.Exact.splitSurjectiveEquiv` に渡す)
を貼り合わせた。 -/
noncomputable def polynomialKaehlerSplit (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    Ω[(Polynomial V)⁄Z] ≃ₗ[Polynomial V]
      (TensorProduct V (Polynomial V) Ω[V⁄Z]) × Polynomial V := by
  have hl : (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)) ∘ₗ
      (LinearMap.smulRight (KaehlerDifferential.polynomialEquiv V).toLinearMap
        (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = LinearMap.id := by
    apply LinearMap.ext
    intro y
    show (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V))
      ((KaehlerDifferential.polynomialEquiv V y) • (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = y
    rw [map_smul]
    have hmapD := KaehlerDifferential.map_D Z V (Polynomial V) (Polynomial V) Polynomial.X
    rw [show algebraMap (Polynomial V) (Polynomial V) Polynomial.X = Polynomial.X from rfl] at hmapD
    rw [hmapD, ← KaehlerDifferential.polynomialEquiv_symm]
    exact (KaehlerDifferential.polynomialEquiv V).symm_apply_apply y
  set e := (Function.Exact.splitSurjectiveEquiv (KaehlerDifferential.exact_mapBaseChange_map Z V (Polynomial V))
    (mapBaseChange_injective_polynomial (Z := Z) (V := V))
    ⟨_, hl⟩).1
  exact e.trans (LinearEquiv.refl (Polynomial V) _ |>.prodCongr (KaehlerDifferential.polynomialEquiv V))

/-- `kerCotangentToTensor` の生成元計算(`range_kerCotangentToTensor_span` の
一般化): 塔 `Z→V→(Polynomial V)→B` で(絶対基底 `Z` と多項式環の係数環 `V`
が別々でよい)、`B = Polynomial V ⧸ (f)` なら
`range(kerCotangentToTensor Z (Polynomial V) B) = (Polynomial V)-span{1⊗D_Z(f)}`。
`Ω_V⊗_VW→Ω_W` の単射性(Lemma 1.1 の第一完全列)を仕上げる際、
`W = AdjoinRoot f` の場合にこの一般化が必要になる
(`range_kerCotangentToTensor_span` 自身は係数環と絶対基底が同じ場合しか
扱えない)。 -/
theorem range_kerCotangentToTensor_span_tower {Z V B : Type*} [CommRing Z] [CommRing V] [CommRing B]
    [Algebra Z V] [Algebra (Polynomial V) B] [Algebra Z B] [Algebra V B]
    [IsScalarTower Z (Polynomial V) B] [IsScalarTower V (Polynomial V) B] [IsScalarTower Z V B]
    (f : Polynomial V)
    (hB : RingHom.ker (algebraMap (Polynomial V) B) = Ideal.span ({f} : Set (Polynomial V))) :
    LinearMap.range (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B) =
      Submodule.span (Polynomial V) {(1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f)} := by
  have hmem : ∀ x : Polynomial V, x * f ∈ RingHom.ker (algebraMap (Polynomial V) B) := by
    intro x; rw [hB]; exact Ideal.mem_span_singleton'.mpr ⟨x, rfl⟩
  set φ := (LinearMap.toSpanSingleton (Polynomial V) (Polynomial V) f).codRestrict _ hmem with hφdef
  have hφsurj : Function.Surjective φ := by
    intro y
    have hy : (y : Polynomial V) ∈ Ideal.span ({f} : Set (Polynomial V)) := hB ▸ y.2
    obtain ⟨x, hx⟩ := Ideal.mem_span_singleton'.mp hy
    exact ⟨x, Subtype.ext hx⟩
  have hcomp : Function.Surjective
      ((RingHom.ker (algebraMap (Polynomial V) B)).toCotangent ∘ₗ φ) :=
    (RingHom.ker (algebraMap (Polynomial V) B)).toCotangent_surjective.comp hφsurj
  have hrange1 : LinearMap.range (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B) =
      LinearMap.range (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B ∘ₗ
        ((RingHom.ker (algebraMap (Polynomial V) B)).toCotangent ∘ₗ φ)) := by
    rw [LinearMap.range_comp, LinearMap.range_eq_top.mpr hcomp, Submodule.map_top]
  rw [hrange1]
  have hfB : algebraMap (Polynomial V) B f = 0 := by
    have : f ∈ RingHom.ker (algebraMap (Polynomial V) B) := hB ▸ Ideal.mem_span_singleton_self f
    exact this
  have heq : (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B ∘ₗ
      ((RingHom.ker (algebraMap (Polynomial V) B)).toCotangent ∘ₗ φ)) =
      LinearMap.toSpanSingleton (Polynomial V) _
        ((1:B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f)) := by
    apply LinearMap.ext
    intro x
    show KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B
      ((RingHom.ker (algebraMap (Polynomial V) B)).toCotangent (φ x)) = _
    rw [KaehlerDifferential.kerCotangentToTensor_toCotangent]
    show (1:B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) (x*f))
        = x • ((1:B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f))
    rw [Derivation.leibniz, TensorProduct.tmul_add]
    have hzero : (1:B) ⊗ₜ[Polynomial V] (f • KaehlerDifferential.D Z (Polynomial V) x) = 0 := by
      rw [TensorProduct.tmul_smul]
      show (f • (1:B)) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) x) = 0
      rw [Algebra.smul_def, hfB, zero_mul, TensorProduct.zero_tmul]
    rw [hzero, add_zero, TensorProduct.tmul_smul]
  rw [heq, LinearMap.range_toSpanSingleton]

/-- **`KaehlerDifferential.map` の3段の塔での合成則(完成)**: 塔
`Z→V→A→B` で `map Z Z A B ∘ map Z Z V A = map Z Z V B`。mathlib には
この合成則の既製品が見当たらなかった——`Ω[V⁄Z]` が `D Z V` の像で
`V`-生成されること(`span_range_derivation`)を使い、生成元上で
`map_D` を2回・スカラータワー`IsScalarTower V A B`の
`algebraMap_apply`で確認し、`Submodule.span_induction` で加法・
`V`-スカラー倍に持ち上げた(`LinearMap.map_smul_of_tower` が鍵——
`map Z Z A B` は `A`-線形だが `V` を経由したスカラー倍とも両立する)。

Lemma 1.1 の単射性(`Ω_V⊗_VW→Ω_W`)を`polynomialKaehlerSplit`
(`Ω_{V[T]/Z}` の直和分解)から `W=AdjoinRoot f` へ橋渡しする際、
`mapBaseChange Z V W` を `mapBaseChange Z (Polynomial V) W` 経由に
factor する(`mapBaseChange_tmul` + 本補題)ために必要になる。 -/
theorem map_comp_map_tower {Z V A B : Type*} [CommRing Z] [CommRing V] [CommRing A] [CommRing B]
    [Algebra Z V] [Algebra V A] [Algebra Z A] [Algebra A B] [Algebra V B] [Algebra Z B]
    [IsScalarTower Z V A] [IsScalarTower Z V B] [IsScalarTower Z A B] [IsScalarTower V A B] :
    ∀ y : Ω[V⁄Z], (KaehlerDifferential.map Z Z A B) ((KaehlerDifferential.map Z Z V A) y) =
      (KaehlerDifferential.map Z Z V B) y := by
  have hspan := KaehlerDifferential.span_range_derivation (R := Z) (S := V)
  intro y
  have hy : y ∈ Submodule.span V (Set.range (KaehlerDifferential.D Z V)) := hspan ▸ Submodule.mem_top
  induction hy using Submodule.span_induction with
  | mem z hz =>
      obtain ⟨v, rfl⟩ := hz
      rw [KaehlerDifferential.map_D, KaehlerDifferential.map_D, KaehlerDifferential.map_D]
      congr 1
      show algebraMap A B (algebraMap V A v) = algebraMap V B v
      rw [← IsScalarTower.algebraMap_apply]
  | zero => simp
  | add y z hy hz ihy ihz => simp only [map_add, ihy, ihz]
  | smul c y hy ih =>
      show (KaehlerDifferential.map Z Z A B) ((KaehlerDifferential.map Z Z V A) (c • y)) =
        (KaehlerDifferential.map Z Z V B) (c • y)
      rw [LinearMap.map_smul_of_tower (KaehlerDifferential.map Z Z V A) c y,
        LinearMap.map_smul_of_tower (KaehlerDifferential.map Z Z A B) c
          ((KaehlerDifferential.map Z Z V A) y),
        LinearMap.map_smul_of_tower (KaehlerDifferential.map Z Z V B) c y, ih]

/-!
## Lemma 1.1「単射性」の完成(2026-09-04)

前節の見取り図の3ステップ(ψ の構成・split 性の base change による
単射性・NZD からの初等的な議論)をすべて実装できた。鍵になったのは
`ψ` を **`LinearMap.baseChange`(スカラー拡大)+ `cancelBaseChange`
(結合律)の合成として明示的に構成する**ことで、「`polynomialKaehlerSplit`
の split 性を base change する」という抽象的な議論を、`LinearMap.
baseChange_comp`/`baseChange_id` という具体的な関手性の等式に
還元できたこと。`ψ` 自身の像が「第一成分のみ」であることを直接
確かめる代わりに、**`ψ` の像の上で `baseChange(map Z V (Polynomial V)
(Polynomial V))` が恒等的に 0 になる**という(より弱いが十分な)事実を
使い、そこから「NZD で割ったら 0 しか残らない」という議論に繋げた。
-/

/-- `mapBaseChange Z V (Polynomial V)` の split 性から作った、
`Ω[(Polynomial V)⁄Z]` の直和分解の(`polynomialKaehlerSplit` と同じ
`Function.Exact.splitSurjectiveEquiv` の呼び出しから作った)版。
以降の補題群では `Ω[(Polynomial V)⁄V]` を`polynomialEquiv` で
`Polynomial V` に単純化する前の形が必要になるため、`polynomialKaehlerSplit`
とは別に用意した。 -/
noncomputable def polyKaehlerSplitEquiv (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    Ω[(Polynomial V)⁄Z] ≃ₗ[Polynomial V]
      (TensorProduct V (Polynomial V) Ω[V⁄Z]) × Ω[(Polynomial V)⁄V] := by
  have hl : (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)) ∘ₗ
      (LinearMap.smulRight (KaehlerDifferential.polynomialEquiv V).toLinearMap
        (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = LinearMap.id := by
    apply LinearMap.ext
    intro y
    show (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V))
      ((KaehlerDifferential.polynomialEquiv V y) • (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = y
    rw [map_smul]
    have hmapD := KaehlerDifferential.map_D Z V (Polynomial V) (Polynomial V) Polynomial.X
    rw [show algebraMap (Polynomial V) (Polynomial V) Polynomial.X = Polynomial.X from rfl] at hmapD
    rw [hmapD, ← KaehlerDifferential.polynomialEquiv_symm]
    exact (KaehlerDifferential.polynomialEquiv V).symm_apply_apply y
  exact (Function.Exact.splitSurjectiveEquiv (KaehlerDifferential.exact_mapBaseChange_map Z V (Polynomial V))
    (mapBaseChange_injective_polynomial (Z := Z) (V := V))
    ⟨_, hl⟩).1

/-- `polyKaehlerSplitEquiv` の定義性質: `mapBaseChange Z V (Polynomial V)`
は `polyKaehlerSplitEquiv` の下で「第一成分への埋め込み `inl`」そのもの
——`Function.Exact.splitSurjectiveEquiv` の返す組の1つ目の性質。 -/
theorem polyKaehlerSplitEquiv_prop (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    KaehlerDifferential.mapBaseChange Z V (Polynomial V) =
      (polyKaehlerSplitEquiv Z V).symm.toLinearMap ∘ₗ LinearMap.inl (Polynomial V) _ _ := by
  have hl : (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)) ∘ₗ
      (LinearMap.smulRight (KaehlerDifferential.polynomialEquiv V).toLinearMap
        (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = LinearMap.id := by
    apply LinearMap.ext
    intro y
    show (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V))
      ((KaehlerDifferential.polynomialEquiv V y) • (KaehlerDifferential.D Z (Polynomial V) Polynomial.X)) = y
    rw [map_smul]
    have hmapD := KaehlerDifferential.map_D Z V (Polynomial V) (Polynomial V) Polynomial.X
    rw [show algebraMap (Polynomial V) (Polynomial V) Polynomial.X = Polynomial.X from rfl] at hmapD
    rw [hmapD, ← KaehlerDifferential.polynomialEquiv_symm]
    exact (KaehlerDifferential.polynomialEquiv V).symm_apply_apply y
  exact (Function.Exact.splitSurjectiveEquiv (KaehlerDifferential.exact_mapBaseChange_map Z V (Polynomial V))
    (mapBaseChange_injective_polynomial (Z := Z) (V := V))
    ⟨_, hl⟩).2.1

/-- `mapBaseChange Z V (Polynomial V)` の明示的な**左逆**(retraction)——
`polyKaehlerSplitEquiv` の第一成分への射影。 -/
noncomputable def polyMapBaseChangeRetraction (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    Ω[(Polynomial V)⁄Z] →ₗ[Polynomial V] TensorProduct V (Polynomial V) Ω[V⁄Z] :=
  LinearMap.fst (Polynomial V) _ _ ∘ₗ (polyKaehlerSplitEquiv Z V).toLinearMap

theorem polyMapBaseChangeRetraction_comp (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    (polyMapBaseChangeRetraction Z V) ∘ₗ (KaehlerDifferential.mapBaseChange Z V (Polynomial V)) = LinearMap.id := by
  unfold polyMapBaseChangeRetraction
  rw [polyKaehlerSplitEquiv_prop]
  ext x
  simp

/-- `KaehlerDifferential.map` の合成 `Ω[V⁄Z] → Ω[(Polynomial V)⁄Z] →
Ω[(Polynomial V)⁄V]` は常に 0(`Ω[V⁄Z]` から来た元は「定数」なので
`Ω[(Polynomial V)⁄V]`(`V` 上の微分)では消える——`Derivation.
map_algebraMap`)。`mapBaseChange Z V B` の単射性の証明で、`ψ` の像が
`kerCotangentToTensor` の値域(NZD 方向)と横断的に交わることを示す
ときの鍵。 -/
theorem map_map_tower_eq_zero (Z V : Type*) [CommRing Z] [CommRing V] [Algebra Z V] :
    ∀ v : Ω[V⁄Z], (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V))
      ((KaehlerDifferential.map Z Z V (Polynomial V)) v) = 0 := by
  intro v
  have hspan := KaehlerDifferential.span_range_derivation (R := Z) (S := V)
  have hv : v ∈ Submodule.span V (Set.range (KaehlerDifferential.D Z V)) := hspan ▸ Submodule.mem_top
  induction hv using Submodule.span_induction with
  | mem z hz =>
      obtain ⟨w, rfl⟩ := hz
      rw [KaehlerDifferential.map_D, KaehlerDifferential.map_D]
      exact Derivation.map_algebraMap (KaehlerDifferential.D V (Polynomial V)) w
  | zero => simp
  | add y z hy hz ihy ihz => rw [map_add, map_add, ihy, ihz, add_zero]
  | smul c y hy ih =>
      show (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V))
        ((KaehlerDifferential.map Z Z V (Polynomial V)) (c • y)) = 0
      rw [LinearMap.map_smul_of_tower (KaehlerDifferential.map Z Z V (Polynomial V)) c y,
        LinearMap.map_smul_of_tower (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)) c
          ((KaehlerDifferential.map Z Z V (Polynomial V)) y), ih, smul_zero]

variable {Z V B : Type*} [CommRing Z] [CommRing V] [CommRing B] [Algebra Z V]
  [Algebra (Polynomial V) B] [Algebra V B] [Algebra Z B]
  [IsScalarTower V (Polynomial V) B] [IsScalarTower Z V B] [IsScalarTower Z (Polynomial V) B]
  [SMulCommClass (Polynomial V) B B]

/-- **`ψ`(`mapBaseChange Z V (Polynomial V)` の `B` への base change、
結合律で書き直したもの)**: `B⊗_VΩ_V → B⊗_{Polynomial V}Ω_{Polynomial V/Z}`。
`LinearMap.baseChange`(スカラー拡大)と `cancelBaseChange`
(`B⊗_{PolyV}((PolyV)⊗_VM) ≅ B⊗_VM`、結合律)の合成として構成した——
これにより「split 性が base change で保たれる」という主張が
`LinearMap.baseChange_comp`/`baseChange_id` という具体的な関手性の
等式に帰着する(下の `psiMap_injective` 参照)。 -/
noncomputable def psiMap :
    TensorProduct V B Ω[V⁄Z] →ₗ[B]
      TensorProduct (Polynomial V) B Ω[(Polynomial V)⁄Z] :=
  (LinearMap.baseChange B (KaehlerDifferential.mapBaseChange Z V (Polynomial V))) ∘ₗ
    (TensorProduct.AlgebraTensorModule.cancelBaseChange V (Polynomial V) B (M := B) (N := Ω[V⁄Z])).symm.toLinearMap

omit [Algebra Z B] [IsScalarTower Z V B] [IsScalarTower Z (Polynomial V) B] in
theorem psiMap_tmul (b : B) (v : Ω[V⁄Z]) :
    psiMap (Z := Z) (V := V) (B := B) (b ⊗ₜ v) =
      b ⊗ₜ[Polynomial V] (KaehlerDifferential.map Z Z V (Polynomial V) v) := by
  show (LinearMap.baseChange B (KaehlerDifferential.mapBaseChange Z V (Polynomial V)))
    ((TensorProduct.AlgebraTensorModule.cancelBaseChange V (Polynomial V) B (M := B) (N := Ω[V⁄Z])).symm (b ⊗ₜ v)) = _
  rw [TensorProduct.AlgebraTensorModule.cancelBaseChange_symm_tmul]
  show b ⊗ₜ[Polynomial V] (KaehlerDifferential.mapBaseChange Z V (Polynomial V) ((1:Polynomial V) ⊗ₜ[V] v)) = _
  rw [KaehlerDifferential.mapBaseChange_tmul, one_smul]

omit [Algebra Z B] [IsScalarTower Z V B] [IsScalarTower Z (Polynomial V) B] in
/-- `ψ` は単射——`mapBaseChange Z V (Polynomial V)` の split 性
(`polyMapBaseChangeRetraction`)を `LinearMap.baseChange` で `B` へ
運んでも split 性が保たれること(`baseChange_comp`/`baseChange_id`)
と、`cancelBaseChange` が同型であることから。 -/
theorem psiMap_injective : Function.Injective (psiMap (Z := Z) (V := V) (B := B)) := by
  have hretr : (LinearMap.baseChange B (polyMapBaseChangeRetraction Z V)) ∘ₗ
      (LinearMap.baseChange B (KaehlerDifferential.mapBaseChange Z V (Polynomial V))) = LinearMap.id := by
    rw [← LinearMap.baseChange_comp, polyMapBaseChangeRetraction_comp, LinearMap.baseChange_id]
  have h1 : Function.Injective (LinearMap.baseChange B (KaehlerDifferential.mapBaseChange Z V (Polynomial V))) := by
    intro a b hab
    have := congrArg (LinearMap.baseChange B (polyMapBaseChangeRetraction Z V)) hab
    simp only [← LinearMap.comp_apply, hretr, LinearMap.id_apply] at this
    exact this
  unfold psiMap
  exact h1.comp
    (TensorProduct.AlgebraTensorModule.cancelBaseChange V (Polynomial V) B (M := B) (N := Ω[V⁄Z])).symm.injective

/-- `mapBaseChange Z V B` は `mapBaseChange Z (Polynomial V) B ∘ ψ` に
factor する(`mapBaseChange_tmul` + `map_comp_map_tower` を単純テンソル
上で確認し、線形性で拡張)。 -/
theorem mapBaseChange_eq_comp_psiMap :
    KaehlerDifferential.mapBaseChange Z V B =
      (KaehlerDifferential.mapBaseChange Z (Polynomial V) B) ∘ₗ (psiMap (Z := Z) (V := V) (B := B)) := by
  apply LinearMap.ext
  intro x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul b v =>
      rw [LinearMap.comp_apply, psiMap_tmul, KaehlerDifferential.mapBaseChange_tmul,
        KaehlerDifferential.mapBaseChange_tmul, map_comp_map_tower]
  | add x y hx hy => simp [map_add, hx, hy]

omit [Algebra Z B] [IsScalarTower Z V B] [IsScalarTower Z (Polynomial V) B] in
/-- `ψ` の像の上では `map Z V (Polynomial V)(Polynomial V)`(の `B` への
base change)は恒等的に 0(`map_map_tower_eq_zero` を単純テンソル上で
適用し、線形性で拡張)。 -/
theorem baseChange_map_comp_psiMap_eq_zero (x : TensorProduct V B Ω[V⁄Z]) :
    (LinearMap.baseChange B (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)))
      (psiMap (Z := Z) (V := V) (B := B) x) = 0 := by
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul b v =>
      rw [psiMap_tmul]
      show b ⊗ₜ[Polynomial V] (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)
        (KaehlerDifferential.map Z Z V (Polynomial V) v)) = 0
      rw [map_map_tower_eq_zero, TensorProduct.tmul_zero]
  | add x y hx hy => rw [map_add, map_add, hx, hy, add_zero]

/-- `π : B⊗_{Polynomial V}Ω_{Polynomial V/V} →ₗ[B] B`(`polynomialEquiv`
を `B` へ base change し、`rid` で `B⊗_{PolyV}(PolyV)≅B` と同一視した
もの)。`1⊗D_Z(f)` の「`dT` 成分」を `f'(w)` として取り出す道具。 -/
noncomputable def piMap : TensorProduct (Polynomial V) B Ω[(Polynomial V)⁄V] →ₗ[B] B :=
  (TensorProduct.AlgebraTensorModule.rid (Polynomial V) B B).toLinearMap ∘ₗ
    (LinearMap.baseChange B (KaehlerDifferential.polynomialEquiv V).toLinearMap)

omit [Algebra V B] [IsScalarTower V (Polynomial V) B] in
theorem piMap_tmul (b : B) (y : Ω[(Polynomial V)⁄V]) :
    piMap (V := V) (B := B) (b ⊗ₜ y) = b * (algebraMap (Polynomial V) B (KaehlerDifferential.polynomialEquiv V y)) := by
  show (TensorProduct.AlgebraTensorModule.rid (Polynomial V) B B)
    (LinearMap.baseChange B (KaehlerDifferential.polynomialEquiv V).toLinearMap (b ⊗ₜ y)) = _
  rw [LinearMap.baseChange_tmul]
  show (TensorProduct.AlgebraTensorModule.rid (Polynomial V) B B) (b ⊗ₜ (KaehlerDifferential.polynomialEquiv V y)) = _
  rw [TensorProduct.AlgebraTensorModule.rid_tmul, Algebra.smul_def, mul_comm]

/-- **[Falt1] Lemma 1.1「単射性」(第一完全列、完成)**: `B = Polynomial V ⧸ (f)`
(`hB`・`hsurj`)で `f'` の像(Falt1 の `f'(w)`)が `B` の非零因子なら、
`mapBaseChange Z V B : B⊗_VΩ_{V/Z} → Ω_{B/Z}` は単射。

証明の流れ: `z ∈ ker(mapBaseChange Z V B)` とすると
`mapBaseChange_eq_comp_psiMap` より `ψ(z) ∈ ker(mapBaseChange Z(PolyV)B)`。
`range_kerCotangentToTensor_span_tower` + mathlib の
`range_kerCotangentToTensor`(`kerCotangentToTensor` の値域が `ker
(mapBaseChange)` の scalar 制限に一致)を貼り合わせて
`ψ(z) = c•(1⊗D_Zf)`(`c:Polynomial V`)と書く。`baseChange_map_comp_
psiMap_eq_zero` で左辺に `baseChange(map)` を当てると 0、右辺は
`π`(`piMap`)を経由して計算すると `algebraMap(c) * f'(w)`——NZD なので
`algebraMap(c)=0`、したがって `ψ(z)=c•(1⊗D_Zf)=0`。最後に `ψ` の単射性
(`psiMap_injective`)から `z=0`。 -/
theorem mapBaseChange_injective_of_nzd (f : Polynomial V)
    (hB : RingHom.ker (algebraMap (Polynomial V) B) = Ideal.span ({f} : Set (Polynomial V)))
    (hsurj : Function.Surjective (algebraMap (Polynomial V) B))
    (hnzd : algebraMap (Polynomial V) B (Polynomial.derivative f) ∈ nonZeroDivisors B) :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V B) := by
  rw [← LinearMap.ker_eq_bot, eq_bot_iff]
  intro z hz
  rw [LinearMap.mem_ker, mapBaseChange_eq_comp_psiMap, LinearMap.comp_apply] at hz
  have hrange_eq : LinearMap.range (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B) =
      Submodule.span (Polynomial V) {(1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f)} :=
    range_kerCotangentToTensor_span_tower f hB
  have hrange_eq2 : LinearMap.range (KaehlerDifferential.kerCotangentToTensor Z (Polynomial V) B) =
      Submodule.restrictScalars (Polynomial V) (KaehlerDifferential.mapBaseChange Z (Polynomial V) B).ker :=
    KaehlerDifferential.range_kerCotangentToTensor Z (Polynomial V) B hsurj
  have hmem : psiMap (Z := Z) (V := V) (B := B) z ∈
      Submodule.span (Polynomial V) {(1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f)} := by
    rw [← hrange_eq, hrange_eq2, Submodule.restrictScalars_mem]
    exact hz
  obtain ⟨c, hc⟩ := Submodule.mem_span_singleton.mp hmem
  have hc' : (algebraMap (Polynomial V) B c) • ((1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f))
      = psiMap (Z := Z) (V := V) (B := B) z := by
    rw [Algebra.algebraMap_eq_smul_one, smul_assoc, one_smul]
    exact hc
  have h0 : (0:B) = piMap (V := V) (B := B)
      ((LinearMap.baseChange B (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)))
        (psiMap (Z := Z) (V := V) (B := B) z)) := by
    rw [baseChange_map_comp_psiMap_eq_zero]; simp
  rw [← hc'] at h0
  have hstep : (LinearMap.baseChange B (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)))
      ((algebraMap (Polynomial V) B c) • ((1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f))) =
      (algebraMap (Polynomial V) B c) •
        (LinearMap.baseChange B (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)))
        ((1 : B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D Z (Polynomial V) f)) := map_smul _ _ _
  rw [hstep] at h0
  rw [LinearMap.baseChange_tmul] at h0
  have hval : (KaehlerDifferential.map Z V (Polynomial V) (Polynomial V)) (KaehlerDifferential.D Z (Polynomial V) f)
      = KaehlerDifferential.D V (Polynomial V) f := by
    have := KaehlerDifferential.map_D Z V (Polynomial V) (Polynomial V) f
    rwa [show algebraMap (Polynomial V) (Polynomial V) f = f from rfl] at this
  rw [hval] at h0
  have hpi : piMap (V := V) (B := B)
      ((algebraMap (Polynomial V) B c) • ((1:B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D V (Polynomial V) f)))
      = (algebraMap (Polynomial V) B c) •
        (piMap (V := V) (B := B) ((1:B) ⊗ₜ[Polynomial V] (KaehlerDifferential.D V (Polynomial V) f))) :=
    map_smul _ _ _
  rw [hpi] at h0
  rw [piMap_tmul, KaehlerDifferential.polynomialEquiv_D, one_mul, smul_eq_mul] at h0
  have hc0 : algebraMap (Polynomial V) B c = 0 :=
    (mul_right_mem_nonZeroDivisors_eq_zero_iff hnzd).mp h0.symm
  have hpsi0 : psiMap (Z := Z) (V := V) (B := B) z = 0 := by
    rw [← hc', hc0, zero_smul]
  exact psiMap_injective (by rw [hpsi0, map_zero])

/-!
## Falt1 の実際の `W` への橋渡し(完成、2026-09-04)

`mapBaseChange_injective_of_nzd` を `B := AdjoinRoot (minpoly V w)`・
`f := minpoly V w` に適用し、`adjoinRootMinpolyEquiv` に沿って
`mapBaseChange Z V (AdjoinRoot (minpoly V w))` の単射性を `mapBaseChange
Z V W` の単射性へ輸送する——`falt1CokernelIso` と同じ橋渡しパターンだが、
**輸送すべき対象が `Ω` ではなく `mapBaseChange`(2つの空間の間の写像)
なので、`omegaCongr` とは別の transport 補題が必要**だった
(`mapBaseChange_injective_transport`、下記)。
-/

/-- `KaehlerDifferential.map` の `B`-slot での自然性(点ごと): `e:A≃ₐ[V]B`
のとき、`mapBaseChange Z V A` と `mapBaseChange Z V B` は `e` を挟んで
両立する。`map_comp_map_tower`(3段の塔 `Z→V→A→B`)を鍵に使う。 -/
theorem mapBaseChange_transport {Z V A B : Type*} [CommRing Z] [CommRing V] [CommRing A] [CommRing B]
    [Algebra Z V] [Algebra V A] [Algebra Z A] [Algebra V B] [Algebra Z B]
    [IsScalarTower Z V A] [IsScalarTower Z V B]
    (e : A ≃ₐ[V] B) :
    letI : Algebra A B := e.toRingHom.toAlgebra
    letI : IsScalarTower V A B := isScalarTower_of_algEquiv e
    letI : IsScalarTower Z A B := isScalarTower_of_algEquiv (e.restrictScalars Z)
    ∀ (b : B) (v : Ω[V⁄Z]),
      (KaehlerDifferential.map Z Z A B)
        (KaehlerDifferential.mapBaseChange Z V A (e.symm b ⊗ₜ v)) =
      KaehlerDifferential.mapBaseChange Z V B (b ⊗ₜ v) := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower V A B := isScalarTower_of_algEquiv e
  haveI : IsScalarTower Z A B := isScalarTower_of_algEquiv (e.restrictScalars Z)
  intro b v
  rw [KaehlerDifferential.mapBaseChange_tmul (R := Z) (A := V) (B := A),
    KaehlerDifferential.mapBaseChange_tmul (R := Z) (A := V) (B := B)]
  have hmapsmul : (KaehlerDifferential.map Z Z A B) (e.symm b • KaehlerDifferential.map Z Z V A v) =
      e.symm b • (KaehlerDifferential.map Z Z A B) (KaehlerDifferential.map Z Z V A v) :=
    map_smul (KaehlerDifferential.map Z Z A B) (e.symm b) (KaehlerDifferential.map Z Z V A v)
  rw [hmapsmul, map_comp_map_tower]
  show algebraMap A B (e.symm b) • (KaehlerDifferential.map Z Z V B) v = b • (KaehlerDifferential.map Z Z V B) v
  congr 1
  show e (e.symm b) = b
  exact e.apply_symm_apply b

/-- **`mapBaseChange` の単射性は代数同型に沿って輸送できる**:
`e:A≃ₐ[V]B` で `mapBaseChange Z V A` が単射なら `mapBaseChange Z V B`
も単射。`φ := (e.symm)⊗id`(`TensorProduct.AlgebraTensorModule.congr`、
全単射)と `mapBaseChange_transport`・`omegaCongr` の単射性を組み合わせた。 -/
theorem mapBaseChange_injective_transport {Z V A B : Type*} [CommRing Z] [CommRing V] [CommRing A] [CommRing B]
    [Algebra Z V] [Algebra V A] [Algebra Z A] [Algebra V B] [Algebra Z B]
    [IsScalarTower Z V A] [IsScalarTower Z V B]
    (e : A ≃ₐ[V] B) (hA : Function.Injective (KaehlerDifferential.mapBaseChange Z V A)) :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V B) := by
  letI : Algebra A B := e.toRingHom.toAlgebra
  haveI : IsScalarTower V A B := isScalarTower_of_algEquiv e
  haveI : IsScalarTower Z A B := isScalarTower_of_algEquiv (e.restrictScalars Z)
  set φ : TensorProduct V B Ω[V⁄Z] ≃ₗ[V] TensorProduct V A Ω[V⁄Z] :=
    TensorProduct.AlgebraTensorModule.congr e.symm.toLinearEquiv (LinearEquiv.refl V Ω[V⁄Z]) with hφdef
  have hkey : ∀ x : TensorProduct V B Ω[V⁄Z],
      (KaehlerDifferential.map Z Z A B) (KaehlerDifferential.mapBaseChange Z V A (φ x)) =
        KaehlerDifferential.mapBaseChange Z V B x := by
    intro x
    induction x using TensorProduct.induction_on with
    | zero => simp
    | tmul b v =>
        show (KaehlerDifferential.map Z Z A B)
          (KaehlerDifferential.mapBaseChange Z V A (e.symm b ⊗ₜ v)) =
          KaehlerDifferential.mapBaseChange Z V B (b ⊗ₜ v)
        exact mapBaseChange_transport e b v
    | add x y hx hy => simp only [map_add, hx, hy]
  have hmapAB_inj : Function.Injective (KaehlerDifferential.map Z Z A B) :=
    (omegaCongr (e.restrictScalars Z)).injective
  intro x y hxy
  have hx0 : (KaehlerDifferential.map Z Z A B) (KaehlerDifferential.mapBaseChange Z V A (φ x)) =
      (KaehlerDifferential.map Z Z A B) (KaehlerDifferential.mapBaseChange Z V A (φ y)) := by
    rw [hkey, hkey]; exact hxy
  have := hA (hmapAB_inj hx0)
  exact φ.injective this

/-- **[Falt1] Lemma 1.1「単射性」、Falt1 の実際の `W` に対して(完成)**:
`V` が Dedekind 整域、`w:W` が `V` 上整・`W` を生成し(かつ拡大体レベル
でも生成する)、`differentIdeal V W ≠ ⊥`(有限可分拡大なら常に成立——
mathlib の `differentIdeal_ne_bot` から従うはずだが、本定理では
仮定として直接渡す)なら、`mapBaseChange Z V W : W⊗_VΩ_{V/Z} → Ω_{W/Z}`
は単射。

証明: `AdjoinRoot(minpoly V w)` は `e`(`adjoinRootMinpolyEquiv`)を
通じて `W`(整域)と同型なので整域。`mapBaseChange_injective_of_nzd`
を `f:=minpoly V w` に適用するための非零因子性は、`e` の下で
`f'(w)`(= `differentIdeal` の生成元、`differentIdeal_eq_span_derivative`)
が非零であることから(`differentIdeal ≠ ⊥` ⟹ 生成元 ≠ 0 ⟹ 整域なので
非零因子)。最後に `mapBaseChange_injective_transport` で `W` へ輸送。 -/
theorem falt1MapBaseChangeInjective {Z V K L W : Type*} [CommRing Z] [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W] [Algebra Z V]
    [Algebra Z W] [IsScalarTower Z V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤)
    (hdne : differentIdeal V W ≠ ⊥) :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V W) := by
  set e := adjoinRootMinpolyEquiv w hint hadjoin with he
  haveI hdom : IsDomain (AdjoinRoot (minpoly V w)) := Function.Injective.isDomain e.toRingHom e.injective
  haveI hZP : IsScalarTower Z (Polynomial V) (AdjoinRoot (minpoly V w)) := by
    apply IsScalarTower.of_algebraMap_eq
    intro x
    rw [IsScalarTower.algebraMap_apply Z V (AdjoinRoot (minpoly V w)),
      IsScalarTower.algebraMap_apply Z V (Polynomial V),
      IsScalarTower.algebraMap_apply V (Polynomial V) (AdjoinRoot (minpoly V w))]
  apply mapBaseChange_injective_transport e
  apply mapBaseChange_injective_of_nzd (minpoly V w) (ker_adjoinRoot_mk (minpoly V w)) AdjoinRoot.mk_surjective
  rw [mem_nonZeroDivisors_iff_ne_zero]
  intro hcontra
  have heval : e (algebraMap (Polynomial V) (AdjoinRoot (minpoly V w)) (Polynomial.derivative (minpoly V w))) = 0 := by
    rw [hcontra, map_zero]
  rw [adjoinRootMinpolyEquiv_algebraMap w hint hadjoin (Polynomial.derivative (minpoly V w))] at heval
  apply hdne
  rw [differentIdeal_eq_span_derivative w hw hadjoin]
  exact Ideal.span_singleton_eq_bot.mpr heval

/-!
★★これで **Lemma 1.1(余核の特定・「長さ」の等式・単射性、Interface
統合込み)は `Found/Falt1/Lemma11.lean` で完成した**。以下は
`Theorem 1.2`(§1 のもう1項目、未完成)へ向けた探索で見つけた新規の
一般補題——「非常に分岐した拡大の列」の各段で使う**複数生成元版**の
Lemma 1.1 に相当する「pushout(テンソル積)の Kähler 微分の直和分解」
である。
-/

/-- `Algebra.IsPushout R C D B`(`B` は `C`・`D` の `R` 上の pushout=
テンソル積)のとき、mathlib の `tensorKaehlerEquiv`(`D⊗_RΩ_{D... }`
ではなく`TensorProduct D B Ω[D⁄R] ≃ₗ[B] Ω[B⁄C]`、**同型**)は、
「`map R C B B` を経由した `mapBaseChange R D B` の像」と**一致する**
——`map R C B B ∘ mapBaseChange R D B = tensorKaehlerEquiv R C D B`。
この自然性が、`tensorKaehlerEquiv.symm` を経由した `mapBaseChange
R D B` が `map R C B B` の**明示的なセクション**であることの鍵になる
(`pushoutKaehlerSplit` で使う)。 -/
theorem tensorKaehlerEquiv_eq_map_mapBaseChange {R C D B : Type*} [CommRing R] [CommRing C] [CommRing D]
    [CommRing B] [Algebra R C] [Algebra R D] [Algebra R B] [Algebra C B] [Algebra D B]
    [IsScalarTower R C B] [IsScalarTower R D B] [Algebra.IsPushout R C D B]
    (x : TensorProduct D B Ω[D⁄R]) :
    (KaehlerDifferential.map R C B B) (KaehlerDifferential.mapBaseChange R D B x)
      = (KaehlerDifferential.tensorKaehlerEquiv R C D B) x := by
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul b y =>
      have hspan := KaehlerDifferential.span_range_derivation (R := R) (S := D)
      have hy : y ∈ Submodule.span D (Set.range (KaehlerDifferential.D R D)) := hspan ▸ Submodule.mem_top
      induction hy using Submodule.span_induction with
      | mem z hz =>
          obtain ⟨d, rfl⟩ := hz
          rw [KaehlerDifferential.mapBaseChange_tmul, KaehlerDifferential.tensorKaehlerEquiv_tmul_D]
          show (KaehlerDifferential.map R C B B) (b • (KaehlerDifferential.map R R D B) (KaehlerDifferential.D R D d)) = _
          rw [map_smul, KaehlerDifferential.map_D]
          show b • (KaehlerDifferential.map R C B B) (KaehlerDifferential.D R B (algebraMap D B d)) = _
          rw [KaehlerDifferential.map_D]
          congr 2
      | zero => simp
      | add y z hy hz ihy ihz =>
          rw [TensorProduct.tmul_add, map_add, map_add, ihy, ihz, map_add]
      | smul d y hy ih =>
          have h1 : (KaehlerDifferential.mapBaseChange R D B) (b ⊗ₜ[D] (d • y))
              = d • (KaehlerDifferential.mapBaseChange R D B) (b ⊗ₜ[D] y) := by
            rw [TensorProduct.tmul_smul]
            exact LinearMap.map_smul_of_tower (KaehlerDifferential.mapBaseChange R D B) d (b ⊗ₜ[D] y)
          have h2 : (KaehlerDifferential.tensorKaehlerEquiv R C D B) (b ⊗ₜ[D] (d • y))
              = d • (KaehlerDifferential.tensorKaehlerEquiv R C D B) (b ⊗ₜ[D] y) := by
            rw [TensorProduct.tmul_smul]
            exact LinearMap.map_smul_of_tower (KaehlerDifferential.tensorKaehlerEquiv R C D B).toLinearMap d (b ⊗ₜ[D] y)
          rw [h1, h2, LinearMap.map_smul_of_tower (KaehlerDifferential.map R C B B) d
            ((KaehlerDifferential.mapBaseChange R D B) (b ⊗ₜ[D] y)), ih]
  | add x y hx hy => simp only [map_add, hx, hy]

/-- **Theorem 1.2 で使う「複数生成元版 Lemma 1.1」の核心(完成)**:
`B` が `C`・`D`(ともに `R`-代数)の `R` 上の pushout(テンソル積)の
とき、`mapBaseChange R C B` が単射なら(`falt1MapBaseChangeInjective`
と同型の条件——`C` 側が `f'` 非零因子な monogenic 拡大なら成立)、
`Ω_{B/R} ≃ₗ[B] (B⊗_CΩ_{C/R}) × (B⊗_DΩ_{D/R})`——`polynomialKaehlerSplit`
と全く同じ「`Function.Exact.splitSurjectiveEquiv` + 明示的セクション」
のパターンだが、セクションは `polynomialEquiv` ではなく
`tensorKaehlerEquiv`(`tensorKaehlerEquiv_eq_map_mapBaseChange` の
自然性)から得た。★★これを `d+1` 個の生成元(帰納的に `pushout` を
繰り返す)へ拡張すれば Theorem 1.2 の「`Ω_{W_{n+1}/V_n}` は `d+1` 個の
直和」という主張になるはずだが、その帰納・`V_n` の族の形式化・
`Module.length` の漸化不等式は未着手。 -/
noncomputable def pushoutKaehlerSplit {R C D B : Type*} [CommRing R] [CommRing C] [CommRing D] [CommRing B]
    [Algebra R C] [Algebra R D] [Algebra R B] [Algebra C B] [Algebra D B]
    [IsScalarTower R C B] [IsScalarTower R D B] [Algebra.IsPushout R C D B]
    (hinj : Function.Injective (KaehlerDifferential.mapBaseChange R C B)) :
    Ω[B⁄R] ≃ₗ[B] (TensorProduct C B Ω[C⁄R]) × (TensorProduct D B Ω[D⁄R]) := by
  have hretr : (KaehlerDifferential.map R C B B) ∘ₗ
      ((KaehlerDifferential.mapBaseChange R D B).comp
        (KaehlerDifferential.tensorKaehlerEquiv R C D B).symm.toLinearMap) = LinearMap.id := by
    apply LinearMap.ext
    intro y
    show (KaehlerDifferential.map R C B B)
      ((KaehlerDifferential.mapBaseChange R D B) ((KaehlerDifferential.tensorKaehlerEquiv R C D B).symm y)) = y
    rw [tensorKaehlerEquiv_eq_map_mapBaseChange]
    exact (KaehlerDifferential.tensorKaehlerEquiv R C D B).apply_symm_apply y
  set e := (Function.Exact.splitSurjectiveEquiv (KaehlerDifferential.exact_mapBaseChange_map R C B)
    hinj ⟨_, hretr⟩).1
  exact e.trans (LinearEquiv.refl B _ |>.prodCongr (KaehlerDifferential.tensorKaehlerEquiv R C D B).symm)

/-!
## `pushoutKaehlerSplit` の仮定 `hinj` を満たす具体例(2026-09-04、完成)

`pushoutKaehlerSplit` を実際に使うには `mapBaseChange R C B` の単射性
(`hinj`)が要る。「非常に分岐した拡大の列」の各段では
`B = C ⊗[R] AdjoinRoot g`(`D := AdjoinRoot g`、`g` は次の生成元の
最小多項式)という形になるので、この特別な場合に `hinj` を
`mapBaseChange_injective_of_nzd`(Lemma 1.1 の証明で使った補題)へ
帰着させる——**単変数の場合の「テンソル積は商と可換」
(`Algebra.TensorProduct.tensorQuotientEquiv`)経由で
`C ⊗[R] AdjoinRoot g ≃ₐ[C] AdjoinRoot (g.map (algebraMap R C))` を
構成し、`mapBaseChange_injective_transport` で輸送する**。
falt1-goal.md で「3つ目(AdjoinRoot基底変換)」と呼んでいた経路。
-/

/-- `C ⊗[R] Polynomial R ≃ₐ[C] Polynomial C`(単変数版の
`MvPolynomial.algebraTensorAlgEquiv`)。`MvPolynomial PUnit` を経由して
構成する(`Polynomial R ≃ₐ[R] MvPolynomial PUnit R` の `uniqueAlgEquiv`
を挟む)。`PUnit` 自身の宇宙を `Type`(宇宙 0)に固定しないと `R C` を
`Type*` にした瞬間に解決不能な宇宙メタ変数が残る(実測済み)。 -/
noncomputable def tensorPolynomialAlgEquiv {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] :
    TensorProduct R C (Polynomial R) ≃ₐ[C] Polynomial C :=
  (Algebra.TensorProduct.congr (AlgEquiv.refl : C ≃ₐ[C] C)
    (MvPolynomial.uniqueAlgEquiv R (PUnit : Type)).symm).trans
  ((MvPolynomial.algebraTensorAlgEquiv R C (σ := (PUnit : Type))).trans
    (MvPolynomial.uniqueAlgEquiv C (PUnit : Type)))

theorem tensorPolynomialAlgEquiv_one_tmul_X {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] :
    (tensorPolynomialAlgEquiv (R := R) (C := C)) ((1:C) ⊗ₜ[R] (Polynomial.X : Polynomial R)) = Polynomial.X := by
  unfold tensorPolynomialAlgEquiv
  simp only [AlgEquiv.trans_apply, Algebra.TensorProduct.congr_apply,
    Algebra.TensorProduct.map_tmul, AlgEquiv.coe_algHom,
    MvPolynomial.algebraTensorAlgEquiv_tmul]
  rw [show (AlgEquiv.refl : C ≃ₐ[C] C) 1 = 1 from rfl, one_smul]
  have hX : (Polynomial.X : Polynomial R)
      = Polynomial.monomial (Finsupp.single (default : PUnit) 1 default) 1 := by
    rw [Finsupp.single_eq_same]; rfl
  rw [hX, MvPolynomial.uniqueAlgEquiv_symm_monomial, MvPolynomial.map_monomial, map_one,
    MvPolynomial.uniqueAlgEquiv_monomial, Finsupp.single_eq_same]
  rfl

/-- `tensorPolynomialAlgEquiv` は基底変換と両立する:
`1 ⊗ₜ g ↦ g.map (algebraMap R C)`。`X` での一致(上)から
`Polynomial.algHom_ext` で `R`-代数準同型として一意性を使う。 -/
theorem tensorPolynomialAlgEquiv_one_tmul {R C : Type*} [CommRing R] [CommRing C] [Algebra R C]
    (g : Polynomial R) :
    (tensorPolynomialAlgEquiv (R := R) (C := C)) ((1:C) ⊗ₜ[R] g) = g.map (algebraMap R C) := by
  set φ : Polynomial R →ₐ[R] Polynomial C :=
    ((tensorPolynomialAlgEquiv (R := R) (C := C)).toAlgHom.restrictScalars R).comp
      (Algebra.TensorProduct.includeRight) with hφ
  set ψ : Polynomial R →ₐ[R] Polynomial C := Polynomial.mapAlgHom (Algebra.ofId R C) with hψ
  have hφX : φ Polynomial.X = Polynomial.X := tensorPolynomialAlgEquiv_one_tmul_X
  have hψX : ψ Polynomial.X = Polynomial.X := by
    rw [hψ]
    show Polynomial.map (Algebra.ofId R C : R →+* C) Polynomial.X = Polynomial.X
    simp
  have heq : φ = ψ := Polynomial.algHom_ext (hφX.trans hψX.symm)
  have hg : φ g = ψ g := by rw [heq]
  rw [hψ] at hg
  show φ g = _
  rw [hg]
  rfl

theorem includeRight_tmul_eq {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    (Algebra.TensorProduct.includeRight : Polynomial R →ₐ[R] TensorProduct R C (Polynomial R)) g
      = (1:C) ⊗ₜ[R] g := rfl

/-- **`C ⊗[R] AdjoinRoot g ≃ₐ[C] AdjoinRoot (g.map (algebraMap R C))`**:
テンソル積と商の可換性(`Algebra.TensorProduct.tensorQuotientEquiv`)を
`AdjoinRoot g = Polynomial R ⧸ span {g}` に適用し、`tensorPolynomialAlgEquiv`
で `Polynomial R` の係数を `C` へ運ぶと、イデアルは
`span {g.map (algebraMap R C)}` に写る(`tensorPolynomialAlgEquiv_one_tmul`
が鍵)。`Ideal.quotientEquivAlg` で結合する。 -/
noncomputable def adjoinRootTensorEquiv {R C : Type*} [CommRing R] [CommRing C] [Algebra R C]
    (g : Polynomial R) :
    TensorProduct R C (AdjoinRoot g) ≃ₐ[C] AdjoinRoot (g.map (algebraMap R C)) := by
  have e1 := Algebra.TensorProduct.tensorQuotientEquiv (R := R) C (Polynomial R) C
    (Ideal.span ({g} : Set (Polynomial R)))
  have hpt : (tensorPolynomialAlgEquiv (R := R) (C := C) : TensorProduct R C (Polynomial R) →+* Polynomial C)
      ((1:C) ⊗ₜ[R] g) = g.map (algebraMap R C) := tensorPolynomialAlgEquiv_one_tmul g
  have hIJ : Ideal.span ({g.map (algebraMap R C)} : Set (Polynomial C))
      = (Ideal.map Algebra.TensorProduct.includeRight (Ideal.span ({g} : Set (Polynomial R)))).map
        (tensorPolynomialAlgEquiv (R := R) (C := C) : TensorProduct R C (Polynomial R) →+* Polynomial C) := by
    rw [Ideal.map_span, Set.image_singleton, Ideal.map_span, Set.image_singleton,
      includeRight_tmul_eq, hpt]
  exact e1.trans (Ideal.quotientEquivAlg _ _ (tensorPolynomialAlgEquiv (R := R) (C := C)) hIJ)

/-- **`pushoutKaehlerSplit` の `hinj` を `B = C ⊗[R] AdjoinRoot g` の場合に
Lemma 1.1 型の非零因子条件へ帰着**: `adjoinRootTensorEquiv` で `B` を
`AdjoinRoot (g.map (algebraMap R C))` に同一視し、
`mapBaseChange_injective_of_nzd`(`falt1MapBaseChangeInjective` と同じ
補題)を適用してから `mapBaseChange_injective_transport` で `B` へ戻す。
これで Theorem 1.2 の各帰納段の `hinj` は、Lemma 1.1 とまったく同じ形
——「その段の最小多項式の微分が(基底変換後も)非零因子」——の
仮定 1 つに帰着する。 -/
theorem mapBaseChange_injective_adjoinRoot_tensor {R C : Type*} [CommRing R] [CommRing C] [Algebra R C]
    (g : Polynomial R)
    (hnzd : algebraMap (Polynomial C) (AdjoinRoot (g.map (algebraMap R C)))
        (Polynomial.derivative (g.map (algebraMap R C)))
        ∈ nonZeroDivisors (AdjoinRoot (g.map (algebraMap R C)))) :
    Function.Injective (KaehlerDifferential.mapBaseChange R C (TensorProduct R C (AdjoinRoot g))) := by
  haveI hRP : IsScalarTower R (Polynomial C) (AdjoinRoot (g.map (algebraMap R C))) := by
    apply IsScalarTower.of_algebraMap_eq
    intro x
    rw [IsScalarTower.algebraMap_apply R C (AdjoinRoot (g.map (algebraMap R C))),
      IsScalarTower.algebraMap_apply R C (Polynomial C),
      IsScalarTower.algebraMap_apply C (Polynomial C) (AdjoinRoot (g.map (algebraMap R C)))]
  apply mapBaseChange_injective_transport (adjoinRootTensorEquiv g).symm
  exact mapBaseChange_injective_of_nzd (g.map (algebraMap R C)) (ker_adjoinRoot_mk _)
    AdjoinRoot.mk_surjective hnzd

/-!
## `pushoutKaehlerSplit` の `d+1` 項化(2026-09-04、3項の場合で実証)

Theorem 1.2 の証明本体(物理p.5=印字p.258、falt1-goal.md に逐語で記録)は
「`Ω_{V_{n+1}/V_n}⊗W_{n+1}` は `d+1` 個の巡回加群の直和」という事実を
使う。これは `pushoutKaehlerSplit` を `d+1` 回反復する(2 項版を
`B_k := C_0⊗_R⋯⊗_RC_k` の帰納法で繰り返す)ことで得られるはずと
falt1-goal.md 3a に記した。ここでは**3項(`d=2`)の場合を実際に
最後まで実行し、帰納の1ステップが機械的に閉じることを確認した**。 -/

/-- **`pushoutKaehlerSplit` の3項版**: `B1` が `C0`・`C1` の pushout、
`B` が `B1`・`C2` の pushout のとき(2段の反復)、
`Ω[B/R] ≃ₗ[B] (B⊗_{C0}Ω[C0/R]) × (B⊗_{C1}Ω[C1/R]) × (B⊗_{C2}Ω[C2/R])`。
`hinj1`(`C0,C1` 段)・`hinj2`(`B1,C2` 段)からの帰結。

証明の骨格(`d+1` 項への一般化もこの2つの道具の反復のみで閉じる):
1. `pushoutKaehlerSplit` を2回適用(`B1` の分解・`B` の分解)。
2. `B1` の分解(`B1`-加群の直積)を `TensorProduct.AlgebraTensorModule.congr`
   で `B` まで底変換し、`TensorProduct.prodRight` でテンソル積と直積の
   可換性を使う。
3. `TensorProduct.AlgebraTensorModule.cancelBaseChange`(「底変換の底変換
   =直接の底変換」)で `B⊗_{B1}(B1⊗_{C0}Ω[C0/R])` を `B⊗_{C0}Ω[C0/R]`
   まで潰す(`C1` も同様)。

`Algebra C0 B`・`Algebra C1 B`・その `IsScalarTower` は**仮定として
明示的に渡す**(`pushoutKaehlerSplit` 自身と同じ流儀)——`B1`・`B` を
具体的な `TensorProduct` として書いてしまうと、この宣言の**型
そのもの**(`TensorProduct C1 B Ω[C1/R]` に要る `Module C1 B`)が
`Algebra.TensorProduct.rightAlgebra`(`abbrev`、tools/lean-idioms.md
#23)を要求するのに対し、宣言の型は `by` ブロックの外(`letI` が使える
前)で先に確定してしまうため、**具体的な `TensorProduct` に特殊化した
補題を別途は書けない**(型それ自体が elaborate できない)ことを実測した
——具体例への適用は、使用箇所ごとに `letI := Algebra.TensorProduct.
rightAlgebra` を先に置いてから `pushoutKaehlerSplit3` を直接呼ぶ
(実際に小さな例で `Module C1 B`・`Algebra.IsPushout` 等すべてが
`inferInstance` で解決することを確認済み)。 -/
noncomputable def pushoutKaehlerSplit3 {R C0 C1 C2 B1 B : Type*} [CommRing R] [CommRing C0] [CommRing C1]
    [CommRing C2] [CommRing B1] [CommRing B]
    [Algebra R C0] [Algebra R C1] [Algebra R C2]
    [Algebra C0 B1] [Algebra C1 B1] [Algebra R B1] [IsScalarTower R C0 B1] [IsScalarTower R C1 B1]
    [Algebra.IsPushout R C0 C1 B1]
    [Algebra B1 B] [Algebra C2 B] [Algebra R B] [IsScalarTower R B1 B] [IsScalarTower R C2 B]
    [Algebra.IsPushout R B1 C2 B]
    [Algebra C0 B] [Algebra C1 B] [IsScalarTower C0 B1 B] [IsScalarTower C1 B1 B]
    (hinj1 : Function.Injective (KaehlerDifferential.mapBaseChange R C0 B1))
    (hinj2 : Function.Injective (KaehlerDifferential.mapBaseChange R B1 B)) :
    Ω[B⁄R] ≃ₗ[B]
      (TensorProduct C0 B Ω[C0⁄R]) × (TensorProduct C1 B Ω[C1⁄R]) × (TensorProduct C2 B Ω[C2⁄R]) := by
  set e1 := pushoutKaehlerSplit (R := R) (C := B1) (D := C2) (B := B) hinj2 with he1
  set e2 := pushoutKaehlerSplit (R := R) (C := C0) (D := C1) (B := B1) hinj1 with he2
  set eC0 := TensorProduct.AlgebraTensorModule.cancelBaseChange C0 B1 B B (Ω[C0⁄R]) with heC0
  set eC1 := TensorProduct.AlgebraTensorModule.cancelBaseChange C1 B1 B B (Ω[C1⁄R]) with heC1
  set e2' : TensorProduct B1 B Ω[B1⁄R] ≃ₗ[B]
      TensorProduct B1 B ((TensorProduct C0 B1 Ω[C0⁄R]) × (TensorProduct C1 B1 Ω[C1⁄R])) :=
    TensorProduct.AlgebraTensorModule.congr (LinearEquiv.refl B B) e2 with he2'
  set eProd := TensorProduct.prodRight B1 B B (TensorProduct C0 B1 Ω[C0⁄R]) (TensorProduct C1 B1 Ω[C1⁄R])
    with heProd
  set eStep : TensorProduct B1 B Ω[B1⁄R] ≃ₗ[B]
      (TensorProduct C0 B Ω[C0⁄R]) × (TensorProduct C1 B Ω[C1⁄R]) :=
    (e2'.trans eProd).trans (eC0.prodCongr eC1) with heStep
  exact e1.trans (eStep.prodCongr (LinearEquiv.refl B _) |>.trans (LinearEquiv.prodAssoc B _ _ _))

/-!
## Theorem 1.2 の4番目のピース: 長さの漸化不等式から `δ_n→0`(2026-09-04、完成)

falt1-goal.md に逐語で記録した証明本文(物理p.5=印字p.258)の最後の一段
落を、`V_n`・`W_n` の具体的な構成に一切依存しない、純粋な実数列の
不等式として抽出し形式化した:

> We derive that δ_n-δ_{n+1} ≥ β-(d+1)(δ_n-δ_{n+1}). So if δ_n ≥ d+1,
> then δ_{n+1} ≤ δ_n-1/(d+2), and otherwise δ_{n+1} ≤
> (1-1/((d+1)(d+2)))δ_n. In any case δ_n → 0 for n→∞.

原文の不等式 `δ_n-δ_{n+1} ≥ β-(d+1)(δ_n-δ_{n+1})`(`β=min{1,δ_n/(d+1)}`)
を整理すると `δ_{n+1} ≤ δ_n - β/(d+2)` になる(これが `hrec` の形)。
証明は2段: (1) `δ_n≥d+1` が続く限り固定量 `1/(d+2)` ずつ減るので、
非負性(`hδ0`)と合わせるとどこかで `δ_N<d+1` になる(有限回で脱出)。
(2) 一旦 `δ_n<d+1` になれば、比 `1-1/((d+1)(d+2))∈(0,1)` の等比的減衰
が続く(その領域に留まることも同時に示す)。以上を `squeeze_zero` と
`tendsto_pow_atTop_nhds_zero_of_lt_one` で結ぶ。 -/

/-- **Theorem 1.2 の長さの漸化不等式 → `δ_n→0`**(物理p.5の証明最終段落、
`V_n`・`W_n` の具体的構成に依存しない抽象形)。`δ n` は `W_n` の
(半局所環の各成分ごとの)different の「長さ」に相当する実数列
(非負・`hδ0`)で、`hrec` は原文の不等式を整理した形
`δ_{n+1} ≤ δ_n - min{1,δ_n/(d+1)}/(d+2)`。 -/
theorem delta_tendsto_zero (d : ℕ) (δ : ℕ → ℝ) (hδ0 : ∀ n, 0 ≤ δ n)
    (hrec : ∀ n, δ (n+1) ≤ δ n - (min 1 (δ n / (d+1))) / (d+2)) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  have hpos : (0:ℝ) < (d:ℝ)+2 := by positivity
  have hpos1 : (0:ℝ) < (d:ℝ)+1 := by positivity
  have hminnn : ∀ n, 0 ≤ min 1 (δ n / (d+1)) := fun n => le_min zero_le_one (div_nonneg (hδ0 n) hpos1.le)
  have hstepA : ∀ n, (d:ℝ)+1 ≤ δ n → δ (n+1) ≤ δ n - 1/((d:ℝ)+2) := by
    intro n hge
    have hmin1 : min (1:ℝ) (δ n / ((d:ℝ)+1)) = 1 := by
      apply min_eq_left; rw [le_div_iff₀ hpos1]; linarith
    have h1 := hrec n
    rw [hmin1] at h1
    exact h1
  set c : ℝ := 1 - 1/(((d:ℝ)+1)*((d:ℝ)+2)) with hc
  have hcpos : 0 < c := by
    rw [hc]; rw [sub_pos, div_lt_one (by positivity)]; nlinarith
  have hclt1 : c < 1 := by
    rw [hc]; have : 0 < 1/(((d:ℝ)+1)*((d:ℝ)+2)) := by positivity
    linarith
  have hstepB : ∀ n, δ n < (d:ℝ)+1 → δ (n+1) ≤ δ n * c := by
    intro n hlt
    have hmin2 : min (1:ℝ) (δ n / ((d:ℝ)+1)) = δ n / ((d:ℝ)+1) := by
      apply min_eq_right; rw [div_le_one hpos1]; linarith
    have h1 := hrec n
    rw [hmin2] at h1
    have heq : δ n - (δ n / ((d:ℝ)+1)) / ((d:ℝ)+2) = δ n * c := by
      rw [hc]; field_simp
    linarith [heq]
  -- 第1段: δ_n≥d+1 が永遠には続かない(有限回で δ_N<d+1 に到達)
  obtain ⟨N, hN⟩ : ∃ N, δ N < (d:ℝ)+1 := by
    by_contra hcon
    push Not at hcon
    have hbound : ∀ n, δ n ≤ δ 0 - n / ((d:ℝ)+2) := by
      intro n
      induction n with
      | zero => simp
      | succ k ih =>
          have hs := hstepA k (hcon k)
          have hsplit : ((k:ℝ)+1)/((d:ℝ)+2) = (k:ℝ)/((d:ℝ)+2) + 1/((d:ℝ)+2) := by ring
          push_cast
          linarith [hsplit]
    obtain ⟨n, hn⟩ := exists_nat_gt (δ 0 * ((d:ℝ)+2))
    have hb := hbound n
    have h0 := hδ0 n
    have hfinal : (n:ℝ) / ((d:ℝ)+2) > δ 0 := by
      rw [gt_iff_lt, lt_div_iff₀ hpos]; linarith
    linarith
  -- 第2段: δ_N<d+1 以降は比 c∈(0,1) の等比的減衰(その領域に留まることも同時に示す)
  have hbelow : ∀ k, δ (k+N) < (d:ℝ)+1 ∧ δ (k+N) ≤ δ N * c^k := by
    intro k
    induction k with
    | zero => constructor <;> simp [hN]
    | succ j ih =>
        obtain ⟨ihlt, ihle⟩ := ih
        have hstep := hstepB (j+N) ihlt
        constructor
        · have : δ (j+N) * c < (d:ℝ)+1 := by nlinarith [hδ0 (j+N)]
          calc δ (j+1+N) = δ (j+N+1) := by ring_nf
            _ ≤ δ (j+N) * c := hstep
            _ < (d:ℝ)+1 := this
        · calc δ (j+1+N) = δ (j+N+1) := by ring_nf
            _ ≤ δ (j+N) * c := hstep
            _ ≤ (δ N * c^j) * c := by
                apply mul_le_mul_of_nonneg_right ihle hcpos.le
            _ = δ N * c^(j+1) := by ring
  have htail : Filter.Tendsto (fun k => δ (k+N)) Filter.atTop (nhds 0) := by
    apply squeeze_zero (f := fun k => δ (k+N)) (g := fun k => δ N * c^k)
    · intro k; exact hδ0 (k+N)
    · intro k; exact (hbelow k).2
    · simpa using Filter.Tendsto.const_mul (δ N) (tendsto_pow_atTop_nhds_zero_of_lt_one hcpos.le hclt1)
  rwa [Filter.tendsto_add_atTop_iff_nat N] at htail

end ABC3.Found.Falt1
