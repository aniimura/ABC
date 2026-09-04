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
import Mathlib.RingTheory.Polynomial.Eisenstein.Basic
import Mathlib.RingTheory.Polynomial.GaussLemma
import Mathlib.RingTheory.DedekindDomain.IntegralClosure
import Mathlib.RingTheory.DedekindDomain.PID
import ABC3.Found.GenEll.TameRamification

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
## `pushoutKaehlerSplit` の任意個への一般化: 帰納の1段(2026-09-04)

`pushoutKaehlerSplit3`(3項=`d=2`の場合の実証)を見て、「前の状態」を
`2つの決まった因子の直積`ではなく**任意の`Fintype`で添字付けられた
族**として抽象化すれば、帰納の**1段**を切り出せることに気づいた——
`TensorProduct.piRight`(`N⊗[R](∀i,M i) ≃ ∀i,N⊗[R]M i`、任意の
`Fintype`添字で無条件に成立、mathlib既存)と`LinearEquiv.piCongrRight`
(各成分ごとの同型をPi型全体の同型に持ち上げる、mathlib既存)を使えば、
`pushoutKaehlerSplit3`の「2段目」の議論がそのまま**任意個の「前の
因子」**に対して機械的に一般化できる。これにより「`Fin(d+1)`版を
1本の定理として書く」代わりに、**「1個追加する」という帰納の1段を
汎用的な補題として確立**した——実際の使用(Theorem 1.2 の`V_n`塔の
構成)ではこの1段を`n`回繰り返し適用することになるはずで、それ自体が
`V_n`の帰納的構成と自然に噛み合う形。`Fin.snoc`等での添字の付け替え
を経た「`Fin(d+1)`全体を1つの定理として閉じる」ラッパーはここでは
作らない(次回以降の課題として残す、falt1-goal.md参照)。 -/

set_option maxHeartbeats 800000 in
/-- **`pushoutKaehlerSplit` の帰納の1段**: 「`Ω[B1/R]`が(`Fintype`
`ι`で添字付けられた)族`F`に関して`∀i,TensorProduct(F i)B1 Ω[F i/R]`
に分解している」という前提(`prev`)のもとで、`B`が`B1`と`C`の
pushoutなら、`Ω[B/R]`はその族に`C`を1個追加した形に分解する。
`pushoutKaehlerSplit`(1段の分解)と`TensorProduct.piRight`・
`TensorProduct.AlgebraTensorModule.cancelBaseChange`(各成分ごとの
底変換の解消)を組み合わせるだけ。 -/
noncomputable def pushoutKaehlerSplitStep {R B1 C B : Type*} {ι : Type*} [Fintype ι] [DecidableEq ι] (F : ι → Type*)
    [CommRing R] [CommRing B1] [CommRing C] [CommRing B] [∀ i, CommRing (F i)]
    [Algebra R B1] [Algebra R C] [Algebra R B] [∀ i, Algebra R (F i)]
    [Algebra B1 B] [Algebra C B] [IsScalarTower R B1 B] [IsScalarTower R C B] [Algebra.IsPushout R B1 C B]
    [∀ i, Algebra (F i) B1] [∀ i, IsScalarTower R (F i) B1]
    [∀ i, Algebra (F i) B] [∀ i, IsScalarTower R (F i) B] [∀ i, IsScalarTower (F i) B1 B]
    (prev : Ω[B1⁄R] ≃ₗ[B1] (∀ i, TensorProduct (F i) B1 Ω[F i⁄R]))
    (hinj : Function.Injective (KaehlerDifferential.mapBaseChange R B1 B)) :
    Ω[B⁄R] ≃ₗ[B] (∀ i, TensorProduct (F i) B Ω[F i⁄R]) × TensorProduct C B Ω[C⁄R] := by
  set e1 := ABC3.Found.Falt1.pushoutKaehlerSplit (R := R) (C := B1) (D := C) (B := B) hinj with he1
  set e2 : TensorProduct B1 B Ω[B1⁄R] ≃ₗ[B] TensorProduct B1 B (∀ i, TensorProduct (F i) B1 Ω[F i⁄R]) :=
    TensorProduct.AlgebraTensorModule.congr (LinearEquiv.refl B B) prev with he2
  set ePi := TensorProduct.piRight B1 B B (fun i => TensorProduct (F i) B1 Ω[F i⁄R]) with hePi
  set eCancel : ∀ i, TensorProduct B1 B (TensorProduct (F i) B1 Ω[F i⁄R]) ≃ₗ[B] TensorProduct (F i) B Ω[F i⁄R] :=
    fun i => TensorProduct.AlgebraTensorModule.cancelBaseChange (F i) B1 B B (Ω[F i⁄R]) with heCancel
  set eStep : TensorProduct B1 B Ω[B1⁄R] ≃ₗ[B] (∀ i, TensorProduct (F i) B Ω[F i⁄R]) :=
    (e2.trans ePi).trans (LinearEquiv.piCongrRight eCancel) with heStep
  exact e1.trans (eStep.prodCongr (LinearEquiv.refl B _))

set_option maxHeartbeats 400000 in
/-- **帰納の起点**: `Ω[C0/R]`は「`Fin 1`で添字付けられた、`C0`だけの
族」に(自明に)分解している——`TensorProduct.lid`(`C0⊗[C0]M≃M`)と
`LinearEquiv.funUnique`(単元集合上のPi型≃その1点での値)を組み合わせる
だけ。`pushoutKaehlerSplitStep`を繰り返し適用する帰納の`prev`の
出発点として使う。 -/
noncomputable def pushoutKaehlerSplitBase {R C0 : Type*} [CommRing R] [CommRing C0] [Algebra R C0] :
    Ω[C0⁄R] ≃ₗ[C0] (∀ _ : Fin 1, TensorProduct C0 C0 Ω[C0⁄R]) :=
  (TensorProduct.lid C0 Ω[C0⁄R]).symm.trans (LinearEquiv.funUnique (Fin 1) C0 _).symm

/-!
## `pushoutKaehlerSplitStep` を`n`段連鎖させるための Σ 束ね基盤
(2026-09-05、**完成**)

`pushoutKaehlerSplitStep`を繰り返し適用するには、「前段までの添字族
`F : ι → Type*`」と「新しい1因子`C : Type*`」を`Option ι`(または
`Fin(n+1)`)で束ねた**新しい**型族へ組み替える必要がある。ところが
`Option.elim i C F`のように**型そのもの**を自由変数`i`で分岐させると、
`CommRing(Option.elim i C F)`のようなinstanceは`i`が具体的に
`none`/`some j`と分からない限り`infer_instance`では解決できない
(`recursive Type-valued def with dependent instances`、tools/
lean-idioms.md 相当の地雷)。

**解決の核**: 「環+その環が運ぶinstance」を**1つのデータ**(Σ型)
として束ねれば、`.ring`/`.algR`のような**射影**は`i`の分岐無しに
常に型付けされる(Σ型の外側の`Option.elim`の結果自体は自由変数の
ままでも、射影は無条件に使える)。`RAlg`(`R`上の代数)と`RAlgOver`
(`R`上の代数で、かつ具体的な現在の累積環`T`へも埋め込まれている
もの)の2段構えで束ね、`RAlgOver.lift`で「`B1`への埋め込み」を
「`B1↪B`の拡大」越しに「`B`への埋め込み」へ運ぶ。これで
`pushoutKaehlerSplitStep`が要求する`B1`関係・`B`関係の両方の
instanceを、自由変数`i`のままで供給できることを確認した
(`letI`(`haveI`ではなく——値の透明性が要る)で登録すれば、
`pushoutKaehlerSplitStep`自身の`∀i,Algebra(F i)B`等の要求は
すべて解決する)。

**`RAlg`/`RAlgOver`の`carrier : Type*`には明示的なuniverse変数
(`uFalt1R`・`uFalt1Car`・`uFalt1T`・`uFalt1B`)を与える**——これが
最後まで残っていた本当の原因だった: `Type*`はその出現ごとに
**独立の**universeメタ変数を生成するため、`Option.elim i C F`の
ような式を(同じ`C`・`F`から)複数箇所に書くと、それぞれの出現が
異なるuniverseメタ変数を持ちうる。両辺の`pp`出力が文字どおり
完全一致するのに`Type mismatch`になる、という2026-09-05に観測した
謎(instance diamond説を疑ったが違った)は、実はこの
**見えないuniverse不一致**が原因だった——`set_option pp.universes
true`で初めて可視化できた。さらに調べたところ、真の最終ブロッカーは
**`(A × B) ≃ₗ[B] C`のような式で、明示的な外側の括弧を省略すると
`A × (B ≃ₗ[B] C)`に誤って構文解析される**(`×`と`≃ₗ[·]`の優先順位)
という**単純な構文の罠**だった——instance diamondでもuniverseの
問題でもなく、自分の書いたコードの括弧不足が最後の障害物だった。 -/

universe uFalt1R uFalt1Car uFalt1T uFalt1B

/-- **Σ束ね(1)**: `R`上の代数を、それが運ぶ`CommRing`・`Algebra R`
インスタンスごと1つのデータとして束ねたもの。`Option.elim`/
`Fin.cons`等で型族を自由変数の添字で組み替えても、`.ring`・`.alg`
という**射影**は分岐無しに使えるため、instance解決の壁を回避できる。 -/
structure RAlg (R : Type uFalt1R) [CommRing R] where
  carrier : Type uFalt1Car
  [ring : CommRing carrier]
  [alg : Algebra R carrier]

attribute [instance] RAlg.ring RAlg.alg

/-- **Σ束ね(2)**: `RAlg R`に加え、具体的な「現在の累積環」`T`への
埋め込み(`Algebra carrier T`・`IsScalarTower R carrier T`)まで
一緒に運ぶ。`pushoutKaehlerSplitStep`の反復で「前段までの各因子が
今の累積環にどう埋め込まれているか」を自由変数の添字`i`のままで
記録するために使う。 -/
structure RAlgOver (R : Type uFalt1R) (T : Type uFalt1T) [CommRing R] [CommRing T] [Algebra R T] where
  carrier : Type uFalt1Car
  [ring : CommRing carrier]
  [algR : Algebra R carrier]
  [algT : Algebra carrier T]
  [towerT : IsScalarTower R carrier T]

attribute [instance] RAlgOver.ring RAlgOver.algR RAlgOver.algT RAlgOver.towerT

/-- **`RAlgOver`の運搬**: `B1`への埋め込みを、さらに`B1↪B`という
1段の代数拡大越しに`B`への埋め込みへ運ぶ(合成`carrier→B1→B`を
`RingHom.toAlgebra`で新しい`Algebra carrier B`として登録するだけ)。
反復pushoutの各段で「前段の族を次の累積環へ運ぶ」ために使う。 -/
noncomputable def RAlgOver.lift {R : Type uFalt1R} {B1 : Type uFalt1T} {B : Type uFalt1B}
    [CommRing R] [CommRing B1] [CommRing B]
    [Algebra R B1] [Algebra R B] [Algebra B1 B] [IsScalarTower R B1 B]
    (x : RAlgOver.{uFalt1R, uFalt1Car, uFalt1T} R B1) : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B :=
  letI algT' : Algebra x.carrier B := ((algebraMap B1 B).comp (algebraMap x.carrier B1)).toAlgebra
  { carrier := x.carrier
    algT := algT'
    towerT := by
      apply IsScalarTower.of_algebraMap_eq
      intro r
      show algebraMap R B r = algebraMap B1 B (algebraMap x.carrier B1 (algebraMap R x.carrier r))
      rw [IsScalarTower.algebraMap_apply R B1 B, IsScalarTower.algebraMap_apply R x.carrier B1] }

/-- `RAlgOver.lift`で運んだ埋め込みは、元の`carrier→B1`埋め込みと
`B1→B`拡大の合成そのものなので、`carrier→B1→B`という
`IsScalarTower`をなす(`pushoutKaehlerSplitStep`の`∀i,IsScalarTower
(F i) B1 B`要求に直接対応する)。 -/
theorem RAlgOver.lift_isScalarTower {R : Type uFalt1R} {B1 : Type uFalt1T} {B : Type uFalt1B}
    [CommRing R] [CommRing B1] [CommRing B]
    [Algebra R B1] [Algebra R B] [Algebra B1 B] [IsScalarTower R B1 B]
    (x : RAlgOver.{uFalt1R, uFalt1Car, uFalt1T} R B1) :
    letI : Algebra x.carrier B := (x.lift (B := B)).algT
    IsScalarTower x.carrier B1 B := by
  letI halgT : Algebra x.carrier B := (x.lift (B := B)).algT
  apply IsScalarTower.of_algebraMap_eq
  intro r
  show algebraMap x.carrier B r = algebraMap B1 B (algebraMap x.carrier B1 r)
  rfl

set_option maxHeartbeats 800000 in
/-- **`pushoutKaehlerSplitStep`を`Option ι`で1段連鎖させる、本節の
到達点**: `RAlg`/`RAlgOver`/`RAlgOver.lift`を使い、`pushoutKaehler
SplitStep`の出力(`(∀i,...)×TensorProduct C B Ω[C/R]`という積の形)
を、次の段の`prev`としてそのまま使える**単一のPi型**
(`∀i:Option ι,...`)へ組み替える。この型を`toFun`/`invFun`から
自前で組み立てる(`LinearEquiv.piOptionEquivProd`は使わない——
`have`/`set`を1つずつ積む形で組み立てれば`Prod`優先順位の罠さえ
避ければ問題無く閉じることが判明した)。 -/
noncomputable def pushoutKaehlerSplitStepOption {R : Type uFalt1R} {B1 : Type uFalt1T} {B : Type uFalt1B}
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    [CommRing R] [CommRing B1] [CommRing B]
    [Algebra R B1] [Algebra R B]
    (C : RAlg.{uFalt1R, uFalt1Car} R) (F : ι → RAlgOver.{uFalt1R, uFalt1Car, uFalt1T} R B1)
    [Algebra C.carrier B1] [Algebra B1 B] [Algebra C.carrier B] [IsScalarTower R B1 B] [IsScalarTower R C.carrier B]
    [Algebra.IsPushout R B1 C.carrier B]
    (prev : Ω[B1⁄R] ≃ₗ[B1] (∀ i, TensorProduct (F i).carrier B1 Ω[(F i).carrier⁄R]))
    (hinj : Function.Injective (KaehlerDifferential.mapBaseChange R B1 B)) :
    Ω[B⁄R] ≃ₗ[B] (∀ i : Option ι,
      TensorProduct (Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B)
        (fun j => (F j).lift (B := B))).carrier B
        Ω[(Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B)
          (fun j => (F j).lift (B := B))).carrier⁄R]) := by
  letI hAlgB : ∀ i, Algebra (F i).carrier B := fun i => (F i).lift (B := B) |>.algT
  letI hTowB : ∀ i, IsScalarTower R (F i).carrier B := fun i => (F i).lift (B := B) |>.towerT
  letI hTowB1B : ∀ i, IsScalarTower (F i).carrier B1 B := fun i => (F i).lift_isScalarTower (B := B)
  have e := ABC3.Found.Falt1.pushoutKaehlerSplitStep (R := R) (B1 := B1) (C := C.carrier) (B := B)
    (fun i => (F i).carrier) prev hinj
  set eSwap := LinearEquiv.prodComm B (∀ i, TensorProduct (F i).carrier B Ω[(F i).carrier⁄R])
    (TensorProduct C.carrier B Ω[C.carrier⁄R]) with heSwap
  set FLifted : ι → RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B := fun j => (F j).lift (B := B) with hFLifted
  set toF : ((TensorProduct C.carrier B Ω[C.carrier⁄R]) × (∀ i, TensorProduct (F i).carrier B Ω[(F i).carrier⁄R])) →
      (∀ i : Option ι, TensorProduct (Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier B
        Ω[(Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier⁄R]) :=
    fun p i => match i with
      | none => p.1
      | some j => p.2 j with htoF
  have hadd : ∀ p q, toF (p + q) = toF p + toF q := by
    intro p q; funext i; cases i <;> rfl
  have hsmul : ∀ (c : B) p, toF (c • p) = (RingHom.id B) c • toF p := by
    intro c p; funext i; cases i <;> rfl
  set invF : (∀ i : Option ι, TensorProduct (Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier B
        Ω[(Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier⁄R]) →
      ((TensorProduct C.carrier B Ω[C.carrier⁄R]) × (∀ i, TensorProduct (F i).carrier B Ω[(F i).carrier⁄R])) :=
    fun f => (f none, fun j => f (some j)) with hinvF
  have hleft : ∀ p, invF (toF p) = p := fun p => rfl
  have hright : ∀ f, toF (invF f) = f := by
    intro f; funext i; cases i <;> rfl
  set eProdOpt : ((TensorProduct C.carrier B Ω[C.carrier⁄R]) × (∀ i, TensorProduct (F i).carrier B Ω[(F i).carrier⁄R])) ≃ₗ[B]
      (∀ i : Option ι, TensorProduct (Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier B
        Ω[(Option.elim i (⟨C.carrier⟩ : RAlgOver.{uFalt1R, uFalt1Car, uFalt1B} R B) FLifted).carrier⁄R]) :=
    LinearEquiv.mk (LinearMap.mk (AddHom.mk toF hadd) hsmul) invF hleft hright with heProdOpt
  exact e.trans (eSwap.trans eProdOpt)

/-!
## `Module.length` の `Pi` 型への加法性(2026-09-04)

`pushoutKaehlerSplitStep` の結論(`LinearEquiv`)を実際に「長さ」の
言葉に翻訳するには、`Module.length R (∀i,M i) = ∑i,Module.length R
(M i)`(有限添字の`Pi`型で長さが和で分解する)という一般事実が要る
——mathlibには`Module.length_prod`(2項の直積のみ)と`Module.length_pi`
(**定数**族`ι→M`の場合のみ、`ENat.card ι * length M`)はあるが、
**変動する族**に対する一般の`Pi`版は見当たらなかった。`Fin.consLinearEquiv`
(`M 0 × (Π i,M(succ i)) ≃ₗ Π i,M i`)による`Fin n`上の帰納法+
`LinearEquiv.piCongrLeft`(添字の付け替え)で、一般の`Fintype`へ
証明した。 -/

set_option maxHeartbeats 800000 in
/-- `Module.length` の `Fin n` 上の `Pi` 型への加法性(帰納法の核)。 -/
theorem module_length_pi_fin {R : Type*} [Ring R] : ∀ (n : ℕ) (M : Fin n → Type*) [∀ i, AddCommGroup (M i)]
    [∀ i, Module R (M i)], Module.length R (∀ i, M i) = ∑ i, Module.length R (M i)
  | 0, M, _, _ => by
      simp only [Finset.univ_eq_empty, Finset.sum_empty]
      have : Subsingleton (∀ i : Fin 0, M i) := by
        constructor; intro a b; funext i; exact absurd i.2 (by omega)
      exact Module.length_eq_zero
  | (n+1), M, _, _ => by
      have e := Fin.consLinearEquiv (R := R) M
      rw [← LinearEquiv.length_eq e, Module.length_prod, module_length_pi_fin n (fun i => M i.succ)]
      rw [Fin.sum_univ_succ]

set_option maxHeartbeats 800000 in
/-- `Module.length` の一般の `Fintype` 添字の `Pi` 型への加法性(`module_
length_pi_fin` を `Fintype.equivFin` で添字の付け替えをするだけ)。 -/
theorem module_length_pi {R : Type*} [Ring R] {ι : Type*} [Fintype ι] (M : ι → Type*)
    [∀ i, AddCommGroup (M i)] [∀ i, Module R (M i)] :
    Module.length R (∀ i, M i) = ∑ i, Module.length R (M i) := by
  set e := (Fintype.equivFin ι).symm with he
  have e2 := LinearEquiv.piCongrLeft R M e
  rw [← LinearEquiv.length_eq e2]
  rw [module_length_pi_fin (Fintype.card ι) (fun i' => M (e i'))]
  rw [← Equiv.sum_comp e (fun i => Module.length R (M i))]

set_option maxHeartbeats 800000 in
/-- **`pushoutKaehlerSplitStep` の長さ版**: `module_length_pi` と
`Module.length_prod` を `pushoutKaehlerSplitStep` の `LinearEquiv` に
適用するだけ——帰納の1段で「前の因子たちの長さの和」に新しい因子`C`の
長さを加えるという、`hrec`(Theorem 1.2の漸化不等式)へ向けた核心の
道具。 -/
theorem pushoutKaehlerSplitStep_length {R B1 C B : Type*} {ι : Type*} [Fintype ι] [DecidableEq ι] (F : ι → Type*)
    [CommRing R] [CommRing B1] [CommRing C] [CommRing B] [∀ i, CommRing (F i)]
    [Algebra R B1] [Algebra R C] [Algebra R B] [∀ i, Algebra R (F i)]
    [Algebra B1 B] [Algebra C B] [IsScalarTower R B1 B] [IsScalarTower R C B] [Algebra.IsPushout R B1 C B]
    [∀ i, Algebra (F i) B1] [∀ i, IsScalarTower R (F i) B1]
    [∀ i, Algebra (F i) B] [∀ i, IsScalarTower R (F i) B] [∀ i, IsScalarTower (F i) B1 B]
    (prev : Ω[B1⁄R] ≃ₗ[B1] (∀ i, TensorProduct (F i) B1 Ω[F i⁄R]))
    (hinj : Function.Injective (KaehlerDifferential.mapBaseChange R B1 B)) :
    Module.length B (Ω[B⁄R]) =
      (∑ i, Module.length B (TensorProduct (F i) B Ω[F i⁄R])) + Module.length B (TensorProduct C B Ω[C⁄R]) := by
  have e := ABC3.Found.Falt1.pushoutKaehlerSplitStep (R := R) (B1 := B1) (C := C) (B := B) F prev hinj
  rw [LinearEquiv.length_eq e, Module.length_prod, module_length_pi]

set_option maxHeartbeats 800000 in
/-- **`pushoutKaehlerSplitStepOption` の長さ版**: `pushoutKaehlerSplitStep_length`
の`Option ι`連鎖版。`module_length_pi`(`Option ι`添字のPi型にも
無条件に適用できる)+`Fintype.sum_option`(`Option ι`上の和を
`none`の項と`ι`上の和に分ける)を組み合わせるだけ——これで
`pushoutKaehlerSplitStepOption`を`n`回反復適用した結果の長さが、
各段で追加した因子の長さの和としてそのまま追跡できる。 -/
theorem pushoutKaehlerSplitStepOption_length
    {R : Type uFalt1R} {B1 : Type uFalt1T} {B : Type uFalt1B}
    {ι : Type*} [Fintype ι] [DecidableEq ι]
    [CommRing R] [CommRing B1] [CommRing B]
    [Algebra R B1] [Algebra R B]
    (C : RAlg.{uFalt1R, uFalt1Car} R) (F : ι → RAlgOver.{uFalt1R, uFalt1Car, uFalt1T} R B1)
    [Algebra C.carrier B1] [Algebra B1 B] [Algebra C.carrier B] [IsScalarTower R B1 B] [IsScalarTower R C.carrier B]
    [Algebra.IsPushout R B1 C.carrier B]
    (prev : Ω[B1⁄R] ≃ₗ[B1] (∀ i, TensorProduct (F i).carrier B1 Ω[(F i).carrier⁄R]))
    (hinj : Function.Injective (KaehlerDifferential.mapBaseChange R B1 B)) :
    Module.length B (Ω[B⁄R]) =
      Module.length B (TensorProduct C.carrier B Ω[C.carrier⁄R]) +
      ∑ j, Module.length B (TensorProduct ((F j).lift (B := B)).carrier B
        Ω[((F j).lift (B := B)).carrier⁄R]) := by
  have e := pushoutKaehlerSplitStepOption (R := R) (B1 := B1) (B := B) C F prev hinj
  rw [LinearEquiv.length_eq e, module_length_pi, Fintype.sum_option]
  rfl

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

/-!
## Theorem 1.2・3b の代数的骨格: `differentIdeal` の塔の等式(2026-09-04)

falt1-goal.md 3b に記録した通り、`delta_tendsto_zero` の前提 `hrec` を
実際の `V_n`・`W_n` から導くには「`Wₙ⊗_{Vₙ}Vₙ₊₁` が正規化前にどれだけ
`Wₙ₊₁` と異なるか」という核心の困難が残る——ここではその困難とは
**独立に**、`Vₙ→Vₙ₊₁→Wₙ₊₁` と `Vₙ→Wₙ→Wₙ₊₁` という2つの塔に
mathlib 既存の `differentIdeal_eq_differentIdeal_mul_differentIdeal`
(`RingTheory/DedekindDomain/Different.lean`)を適用して得られる、
**厳密な等式**(原文の不等式より強い形)を切り出した。 -/

/-- **`differentIdeal` の菱形(diamond)公式**: `Vₙ→Vₙ₊₁→Wₙ₊₁` と
`Vₙ→Wₙ→Wₙ₊₁` の2通りの塔で `differentIdeal Vₙ Wₙ₊₁` を計算すると
一致することから、`δₙ₊₁ := differentIdeal Vₙ₁ Wₙ₁` と
`δₙ := differentIdeal Vₙ Wₙ` を結ぶ等式が出る。`Algebra.IsSeparable`
の引数は `FractionRing.liftAlgebra` に明示的に固定している
(tools/lean-idioms.md #23 と同じ理由——固定しないと2回目の適用で
`Algebra (FractionRing Vₙ) (FractionRing Wₙ₁)` の探索が食い違う)。

★残る課題(未解決、falt1-goal.md 3b 参照): 右辺の
`differentIdeal Wₙ Wₙ₊₁`(`Wₙ⊗_{Vₙ}Vₙ₊₁` の非正規性を測る量)を
評価する部分は、この等式だけでは閉じない——原文の kernel・cokernel
の length 評価と本質的に同じ困難が残る。 -/
theorem differentIdeal_tower_diamond {Vn Vn1 Wn Wn1 : Type*} [CommRing Vn] [CommRing Vn1] [CommRing Wn]
    [CommRing Wn1] [IsDomain Vn] [IsIntegrallyClosed Vn] [IsDedekindDomain Vn1] [IsDedekindDomain Wn]
    [IsDedekindDomain Wn1] [Algebra Vn Vn1] [Algebra Vn Wn] [Algebra Vn1 Wn1] [Algebra Wn Wn1] [Algebra Vn Wn1]
    [IsScalarTower Vn Vn1 Wn1] [IsScalarTower Vn Wn Wn1]
    [Module.IsTorsionFree Vn Vn1] [Module.IsTorsionFree Vn Wn] [Module.IsTorsionFree Vn Wn1]
    [Module.IsTorsionFree Vn1 Wn1] [Module.IsTorsionFree Wn Wn1]
    [Module.Finite Vn Vn1] [Module.Finite Vn Wn] [Module.Finite Vn Wn1]
    [Module.Finite Vn1 Wn1] [Module.Finite Wn Wn1]
    (hsep : @Algebra.IsSeparable (FractionRing Vn) (FractionRing Wn1) _ _
      (FractionRing.liftAlgebra Vn (FractionRing Wn1))) :
    differentIdeal Vn1 Wn1 * Ideal.map (algebraMap Vn1 Wn1) (differentIdeal Vn Vn1)
      = differentIdeal Wn Wn1 * Ideal.map (algebraMap Wn Wn1) (differentIdeal Vn Wn) := by
  letI := FractionRing.liftAlgebra Vn (FractionRing Wn1)
  haveI := hsep
  have h1 : differentIdeal Vn Wn1 = differentIdeal Vn1 Wn1 * Ideal.map (algebraMap Vn1 Wn1) (differentIdeal Vn Vn1) :=
    differentIdeal_eq_differentIdeal_mul_differentIdeal Vn Vn1 Wn1
  have h2 : differentIdeal Vn Wn1 = differentIdeal Wn Wn1 * Ideal.map (algebraMap Wn Wn1) (differentIdeal Vn Wn) :=
    differentIdeal_eq_differentIdeal_mul_differentIdeal Vn Wn Wn1
  rw [← h1, ← h2]

/-!
## Theorem 1.2・3b: `differentIdeal Wₙ Wₙ₊₁`(`Jₙ`)を conductor で消去

`differentIdeal_tower_diamond` の右辺に現れる未知数 `Jₙ := differentIdeal
Wₙ Wₙ₊₁` を、mathlib の `conductor_mul_differentIdeal`(`conductor A x *
differentIdeal A B = span{aeval x (deriv(minpoly A x))}`、Lemma 1.1 の
`differentIdeal_eq_span_derivative` の一般化)と組み合わせて**消去**
できることを発見した(`ResearchPaper/0_Source/Brinon Conrad - CMI
Summer School Notes...` Exercise 13.7.4 の step (5) に対応する箇所、
falt1-goal.md 参照)。2つの等式を掛け合わせて `Jₙ` を消すだけなので、
Dedekind整域の「0でないイデアルで消去できる」性質だけで閉じる
——`conductor Wₙ x`・`differentIdeal Vₙ Vₙ₊₁` の base change が一致する
という仮定(`hspan_eq`、「典型例」で `x` の最小多項式が base change
と両立することに相当、3c を通じて確認すべき事項として残す)のもとで、
`conductor(Wₙ,x) · δₙ₊₁ = δₙ.map(...)` という**きれいな関係式**が
出ることを確認した。 -/

/-- **`Jₙ` を消去した後の `δₙ`・`δₙ₊₁` の関係**。`hdiamond` は
`differentIdeal_tower_diamond` の結論、`hcond` は
`conductor_mul_differentIdeal` の結論、`hspan_eq` は両者の右辺が
一致するという(「典型例」で成り立つはずの、3c 経由で今後確認する)
仮定。`Idiff ≠ 0`(Dedekind整域では0でないイデアルは消去可能)だけで
`Jₙ` を消去できる。 -/
theorem cancel_conductor_delta {R : Type*} [CommRing R] [IsDedekindDomain R]
    (Idiff Jn condWnx spanDeriv deltaN1 deltaNmapped : Ideal R)
    (hdiamond : deltaN1 * Idiff = Jn * deltaNmapped)
    (hcond : condWnx * Jn = spanDeriv)
    (hspan_eq : spanDeriv = Idiff)
    (hIdiff_ne : Idiff ≠ 0) :
    condWnx * deltaN1 = deltaNmapped := by
  have h1 : condWnx * (deltaN1 * Idiff) = condWnx * (Jn * deltaNmapped) := by rw [hdiamond]
  have h2 : condWnx * (deltaN1 * Idiff) = (condWnx * deltaN1) * Idiff := by ring
  have h3 : condWnx * (Jn * deltaNmapped) = deltaNmapped * Idiff := by rw [← hspan_eq, ← hcond]; ring
  rw [h2, h3] at h1
  exact mul_right_cancel₀ hIdiff_ne h1

/-!
## Theorem 1.2・3b(c): 単項イデアル倍の長さの加法性(2026-09-04、完成)

`cancel_conductor_delta` の結論(イデアルの掛け算の等式)を `delta_
tendsto_zero` が要求する `hrec`(長さの引き算の不等式)へ変換するには
「`length(R/(I·J)) = length(R/I)+length(R/J)`」という長さの加法性が
要る。一般の(必ずしも coprime でない)`I,J` に対するこの事実は
mathlib に見当たらなかった(何度も探索・記録した——`Ideal.relNorm`・
`Module.Invertible` 等、遠回りの経路はどれも別の一般補題を要した)。

★ここでは**`I` が単項イデアル `(a)`(`a` は非零因子)の場合**に限定
すれば、`Module.length_eq_add_of_exact`(SES の加法性)+
`Submodule.quotientQuotientEquivQuotient`(第三同型定理)+
`LinearMap.quotKerEquivRange`(第一同型定理)+ `a` による積が
`R ≃ₗ (a)` を与えること、という**既存の道具の組み合わせだけ**で
閉じることを発見した——`Module.Invertible`(可逆加群、局所化)の
一般論は一切不要だった。多くの実際の場面(`differentIdeal_eq_span_
derivative` 等、単一の生成元の微分で書ける場合)では `I` が単項なので、
この特殊形で十分なことが多い。 -/

/-- **単項イデアルで割った長さの分解**: `M` の部分加群 `S ≤ T` に対し
`length(M/S) = length(T/S) + length(M/T)`(第三同型定理 +
`Module.length_eq_add_of_exact`)。`length_quotient_span_singleton_mul`
の骨格として使う、一般の(単項とは限らない)`S,T` に対する補助補題。 -/
theorem length_quotient_of_le {R M : Type*} [CommRing R] [AddCommGroup M] [Module R M]
    (S T : Submodule R M) (h : S ≤ T) :
    Module.length R (M ⧸ S) =
      Module.length R (↥T ⧸ (S.comap T.subtype)) + Module.length R (M ⧸ T) := by
  set f : ↥T →ₗ[R] (M ⧸ S) := S.mkQ.comp T.subtype with hf
  have hker : f.ker = S.comap T.subtype := by
    ext x
    simp [hf, LinearMap.mem_ker, Submodule.mkQ_apply, Submodule.mem_comap]
  have hrange : f.range = Submodule.map S.mkQ T := by
    ext x
    simp [hf, LinearMap.mem_range, Submodule.mem_map]
  have e1 : (↥T ⧸ f.ker) ≃ₗ[R] ↥f.range := LinearMap.quotKerEquivRange f
  rw [hker, hrange] at e1
  have e2 := (Submodule.quotientQuotientEquivQuotient S T h)
  have hlen1 := LinearEquiv.length_eq e1
  have hlen2 := LinearEquiv.length_eq e2
  have hexact : Module.length R (M ⧸ S) = Module.length R (Submodule.map S.mkQ T) +
      Module.length R ((M ⧸ S) ⧸ Submodule.map S.mkQ T) :=
    Module.length_eq_add_of_exact (Submodule.map S.mkQ T).subtype (Submodule.map S.mkQ T).mkQ
      (Submodule.map S.mkQ T).subtype_injective (Submodule.map S.mkQ T).mkQ_surjective
      (LinearMap.exact_subtype_mkQ (Submodule.map S.mkQ T))
  rw [hexact, ← hlen2, ← hlen1]

theorem toSpanSingleton_injective_of_nzd {R : Type*} [CommRing R] (a : R) (ha : a ∈ nonZeroDivisors R) :
    Function.Injective (LinearMap.toSpanSingleton R R a) := by
  intro x y hxy
  simp only [LinearMap.toSpanSingleton_apply, smul_eq_mul] at hxy
  have h : (x - y) * a = 0 := by rw [sub_mul, hxy, mul_comm, sub_self]
  exact sub_eq_zero.mp ((mul_right_mem_nonZeroDivisors_eq_zero_iff ha).mp h)

/-- **長さの加法性(単項イデアル×任意のイデアル)**: `a` が非零因子のとき
`length(R/((a)·J)) = length(R/(a)) + length(R/J)`。3b(c) の核心の
gap を、単項の場合に限定して埋める。 -/
theorem length_quotient_span_singleton_mul {R : Type*} [CommRing R] (a : R)
    (ha : a ∈ nonZeroDivisors R) (J : Ideal R) :
    Module.length R (R ⧸ (Ideal.span ({a}:Set R) * J)) =
      Module.length R (R ⧸ Ideal.span ({a}:Set R)) + Module.length R (R ⧸ J) := by
  set T : Submodule R R := Ideal.span ({a}:Set R) with hT
  set S : Submodule R R := Ideal.span ({a}:Set R) * J with hS
  have hle : S ≤ T := hS ▸ Ideal.mul_le_right
  rw [length_quotient_of_le S T hle, add_comm]
  congr 1
  set f0 : R →ₗ[R] R := LinearMap.toSpanSingleton R R a with hf0
  have hf0range : f0.range = T := by
    rw [hf0, hT]
    show (LinearMap.toSpanSingleton R R a).range = R ∙ a
    rw [← LinearMap.span_singleton_eq_range]
  set g : R ≃ₗ[R] ↥T := LinearEquiv.ofInjective f0 (toSpanSingleton_injective_of_nzd a ha)
    |>.trans (LinearEquiv.ofEq _ _ hf0range) with hg
  have hgapp : ∀ y : R, ((g y : ↥T) : R) = a * y := by
    intro y
    show f0 y = a * y
    rw [hf0, LinearMap.toSpanSingleton_apply, smul_eq_mul, mul_comm]
  apply (LinearEquiv.length_eq (Submodule.Quotient.equiv J (S.comap T.subtype) g ?_)).symm
  ext x
  simp only [Submodule.mem_map, Submodule.mem_comap, Submodule.coe_subtype]
  constructor
  · rintro ⟨y, hy, rfl⟩
    show ((g y : ↥T) : R) ∈ S
    rw [hgapp, hS, Submodule.span_singleton_mul]
    exact ⟨y, hy, rfl⟩
  · intro hx
    rw [hS, Submodule.span_singleton_mul] at hx
    obtain ⟨y, hy, hxy⟩ := hx
    change a * y = (x : R) at hxy
    exact ⟨y, hy, Subtype.ext ((hgapp y).trans hxy)⟩

/-!
## Theorem 1.2・3b: 次数スケーリングは分岐指数(2026-09-04、検算で訂正)

`cancel_conductor_delta` を `delta_tendsto_zero` の `hrec` へ接続する
には「`Wₙ₊₁` へ写したイデアルの長さ」を `Wₙ` 側の長さと結ぶ必要がある。
当初「`[Wₙ₊₁:Wₙ]` 倍」という素朴な式を疑ったが、**不分岐拡大での検算
(両辺とも `1` になり矛盾)で誤りと判明した**——正しい係数は次数
全体ではなく**分岐指数 `e`** のみ。以下でこれを(`Wₙ₊₁` が DVR の
場合に)確認した。 -/

/-- **素イデアルの冪の像の長さは分岐指数倍**: `S`(DVR)の極大イデアルが
`Wₙ` のイデアル `P` の像として `𝔪^e` と書けるとき(`e` = 分岐指数)、
`P^k` の像で割った長さは `e·k`——`[Wₙ₊₁:Wₙ]`(次数、分岐指数×残余次数)
ではなく分岐指数のみに比例する。`IsDiscreteValuationRing.length_
quotient_pow_maximalIdeal` から直接従う。 -/
theorem length_map_pow_of_ramificationIdx {S : Type*} [CommRing S] [IsDomain S] [IsDiscreteValuationRing S]
    (Wn : Type*) [CommRing Wn] [Algebra Wn S] (P : Ideal Wn) (e k : ℕ)
    (he : Ideal.map (algebraMap Wn S) P = (IsLocalRing.maximalIdeal S) ^ e) :
    Module.length S (S ⧸ Ideal.map (algebraMap Wn S) (P^k)) = (e:ℕ∞) * k := by
  rw [Ideal.map_pow, he, ← pow_mul]
  rw [IsDiscreteValuationRing.length_quotient_pow_maximalIdeal]
  push_cast
  ring

/-!
## Theorem 1.2・3c: 戦略転換——`AdjoinRoot` の整閉性の証明は不要(2026-09-04)

3c(「非常に分岐した」`V_n` 族の構成)で、`Vₙ₊₁ := AdjoinRoot(f)` が
Eisenstein 性から**既に**整閉であることを証明しようとして、局所化・
付値延長という mathlib にまだ薄い領域に何度も突き当たった。★戦略を
転換した: `Vₙ₊₁ := integralClosure Vₙ L` と**定義**すれば、mathlib の
`integralClosure.isDedekindDomain_fractionRing`(`RingTheory/
DedekindDomain/IntegralClosure.lean`、**instance**、`L` が
`FractionRing Vₙ` 上有限分離拡大なら自動的に成立)から整閉性・
Dedekind性が**タダで手に入る**——`f` の根 `w` が `Vₙ₊₁` を「生成し
切っているか」は問う必要がない。なぜなら `conductor_mul_
differentIdeal`(既に完成済み)・`cancel_conductor_delta`(既に完成
済み)が、その conductor を未知数として追跡できるよう作られている
から。★「`FractionRing Vₙ` を明示的に使う」ことが鍵(一般の `K`
with `IsFractionRing Vₙ K` では instance が発火しない——
tools/lean-idioms.md #23 と同種の罠、実測で確認・回避した)。 -/

theorem coeff_X_pow_sub_C' {R : Type*} [CommRing R] (π : R) (n m : ℕ) (hn : 0 < n) :
    (Polynomial.X ^ n - Polynomial.C π : Polynomial R).coeff m
      = if m = n then 1 else if m = 0 then -π else 0 := by
  simp [Polynomial.coeff_X_pow, Polynomial.coeff_C]
  split_ifs <;> simp_all

/-- **`X^n - π` は `(π)` を基準に Eisenstein**(`IsEisensteinAt` の
構造そのもの)。`eisenstein_X_pow_sub_C`(既約性)と、後の単項生成性
(`adjoin_eq_top_of_isEisensteinAt` 系)の両方がこの1つの事実から
出るので、独立した部品として切り出した(2026-09-04)。 -/
theorem isEisensteinAt_X_pow_sub_C {R : Type*} [CommRing R] [IsDomain R] (π : R) (n : ℕ) (hn : 0 < n)
    (hprime : (Ideal.span ({π} : Set R)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set R))^2) :
    (Polynomial.X ^ n - Polynomial.C π : Polynomial R).IsEisensteinAt (Ideal.span ({π}:Set R)) := by
  have hmonic : (Polynomial.X ^ n - Polynomial.C π : Polynomial R).Monic :=
    Polynomial.monic_X_pow_sub_C π hn.ne'
  have hdeg : (Polynomial.X ^ n - Polynomial.C π : Polynomial R).natDegree = n :=
    Polynomial.natDegree_X_pow_sub_C
  apply hmonic.isEisensteinAt_of_mem_of_notMem
  · exact hprime.ne_top
  · intro m hm
    rw [hdeg] at hm
    rw [coeff_X_pow_sub_C' π n m hn]
    have hmn : m ≠ n := by omega
    simp only [if_neg hmn]
    split_ifs with h0
    · simp
    · simp
  · rw [coeff_X_pow_sub_C' π n 0 hn]
    have hn0 : (0:ℕ) ≠ n := by omega
    simp only [if_neg hn0]
    simpa using hnotsq

/-- **`X^n - π` は `(π)` が非自乗な素イデアルなら既約**(Eisenstein の
判定法の直接の帰結)。「非常に分岐した」`V_n` 族の各段(`p` 乗根の
逐次添加)の生成多項式がこの形——`Vₙ₊₁` を構成する第一歩。 -/
theorem eisenstein_X_pow_sub_C {R : Type*} [CommRing R] [IsDomain R] (π : R) (n : ℕ) (hn : 0 < n)
    (hprime : (Ideal.span ({π} : Set R)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set R))^2) :
    Irreducible (Polynomial.X ^ n - Polynomial.C π) := by
  have hmonic : (Polynomial.X ^ n - Polynomial.C π : Polynomial R).Monic :=
    Polynomial.monic_X_pow_sub_C π hn.ne'
  have hdeg : (Polynomial.X ^ n - Polynomial.C π : Polynomial R).natDegree = n :=
    Polynomial.natDegree_X_pow_sub_C
  have hei := isEisensteinAt_X_pow_sub_C π n hn hprime hnotsq
  refine hei.irreducible hprime hmonic.isPrimitive ?_
  rw [hdeg]; exact hn

/-- **`X^n-π` は分数体上でも既約**(Gauss の補題、`R` が整閉——Dedekind
整域は自動的にそう——なら UFD を要さず成立)。これで `Vₙ₊₁` を
`Frac(Vₙ)` 上の実際の**体**拡大として構成できる——`integralClosure.
isDedekindDomain_fractionRing` を使うために必要な最後のピース。 -/
theorem eisenstein_X_pow_sub_C_irreducible_map {R K : Type*} [CommRing R] [IsDomain R] [IsIntegrallyClosed R]
    [Field K] [Algebra R K] [IsFractionRing R K] (π : R) (n : ℕ) (hn : 0 < n)
    (hprime : (Ideal.span ({π} : Set R)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set R))^2) :
    Irreducible ((Polynomial.X ^ n - Polynomial.C π : Polynomial R).map (algebraMap R K)) := by
  have hirr := eisenstein_X_pow_sub_C π n hn hprime hnotsq
  have hmonic : (Polynomial.X ^ n - Polynomial.C π : Polynomial R).Monic :=
    Polynomial.monic_X_pow_sub_C π hn.ne'
  exact (hmonic.irreducible_iff_irreducible_map_fraction_map (K := K)).mp hirr

/-- `AdjoinRoot f` は `f` が monic既約かつ分離的なら分離拡大——
`IntermediateField.isSeparable_adjoin_simple_iff_isSeparable`
(根 `x` で生成される中間体が分離的 ⟺ `x` が分離的)と
`IntermediateField.adjoin_root_eq_top`(その中間体が `⊤`)を貼り合わせる。 -/
theorem algIsSeparable_adjoinRoot_of_separable {K : Type*} [Field K] (f : Polynomial K) [Fact (Irreducible f)]
    (hmonic : f.Monic) (hsep : f.Separable) : Algebra.IsSeparable K (AdjoinRoot f) := by
  have hne0 : f ≠ 0 := (Fact.out : Irreducible f).ne_zero
  have hminpoly : minpoly K (AdjoinRoot.root f) = f * Polynomial.C f.leadingCoeff⁻¹ :=
    AdjoinRoot.minpoly_root hne0
  rw [hmonic.leadingCoeff, inv_one, Polynomial.C_1, mul_one] at hminpoly
  have hrootsep : IsSeparable K (AdjoinRoot.root f) := by
    show (minpoly K (AdjoinRoot.root f)).Separable
    rw [hminpoly]; exact hsep
  have h1 := (IntermediateField.isSeparable_adjoin_simple_iff_isSeparable K (AdjoinRoot f)).mpr hrootsep
  rw [IntermediateField.adjoin_root_eq_top f] at h1
  haveI := h1
  exact AlgEquiv.Algebra.isSeparable (e := IntermediateField.topEquiv (F := K) (E := AdjoinRoot f))

/-!
## Theorem 1.2・3c: `V_1` の構成が完成した(2026-09-04)

falt1-goal.md に記録した戦略転換(`Vₙ₊₁ := integralClosure Vₙ L`、
`AdjoinRoot` 自体の整閉性は証明しない)を、実際に「非常に分岐した」
族の最初の1段(`X^n-π` 型の Eisenstein 拡大)に対して**最後まで
実行し、sorry 無しで完成させた**——`integralClosure.isDedekindDomain_
fractionRing`(instance)が要求する `FiniteDimensional`・`Algebra.
IsSeparable` を含む全ての条件を、上で構築した部品
(`eisenstein_X_pow_sub_C_irreducible_map`・`Polynomial.
separable_X_pow_sub_C`・`algIsSeparable_adjoinRoot_of_separable`・
`Polynomial.Monic.finite_adjoinRoot`)だけで満たせることを確認した。
局所化・付値延長という mathlib の薄い領域は一切不要だった。 -/

/-- **「非常に分岐した」`V_n` 族の1段の構成が Dedekind であることの
証明**: `V_0` が Dedekind、`π` が `(π)` を非自乗な素イデアルとして
生成し、`n` が `Frac(V_0)` で可逆(標数が `n` を割らない)なら、
`X^n-π` の根を添加した `Frac(V_0)` の拡大体の(`V_0` における)整閉包
は Dedekind 整域になる。`V_n → V_{n+1}` の1段分の具体的なモデル。 -/
theorem isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C
    {V0 : Type*} [CommRing V0] [IsDedekindDomain V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n) :
    IsDedekindDomain (integralClosure V0
      (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0))))) := by
  have hirr : Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))) :=
    eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) := ⟨hirr⟩
  have hmonicK : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Monic :=
    (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  have hmapeq : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))
      = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by simp
  have hsepK : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Separable := by
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Module.Finite (FractionRing V0)
      (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) :=
    hmonicK.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0)
      (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) :=
    ‹Module.Finite _ _›
  haveI : Algebra.IsSeparable (FractionRing V0)
      (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) :=
    algIsSeparable_adjoinRoot_of_separable _ hmonicK hsepK
  infer_instance

/-!
## Theorem 1.2・3c: `V_1` の単項生成性(monogenicity)が完成した(2026-09-04)

`falt1CokernelLengthEq`(Lemma 1.1「長さ」)・`cancel_conductor_delta`
(3b)がともに要求する「`W` は単一の元 `w` で生成される」(`hadjoin :
Algebra.adjoin V ({w}:Set W) = ⊤`)を、上の `V_1 := integralClosure V0
L` の構成に対して**実際に確認した**——これで 3c の技術的な核心
(Eisenstein 拡大の古典的な「単項生成性」定理)を、`AdjoinRoot(f)` の
整閉性を直接証明することなく手に入れた。

★決定的な発見: `ABC3.Found.GenEll.TameRamification`(Mochizuki
[GenEll] の馴分岐の形式化、別の論文トラックだが純粋な可換環論)に、
まさにこの用途の部品が**既に完成していた**——
`mem_adjoin_of_pow_smul_of_isEisensteinAt`(Eisenstein なら `π^k•z∈
adjoin` から `z∈adjoin` が出る、`k` 回の帰納)と
`exists_smul_mem_adjoin_powerBasis`(`L` の任意の元は、ある `0≠d∈R`
倍すれば `adjoin` に入る、`IsLocalization.exist_integer_multiples`
経由)。`R` が DVR なら「任意の `d≠0` は `π^k×単元`」
(`IsDiscreteValuationRing.eq_unit_mul_pow_irreducible`)なので、
この2つを合成するだけで「**任意の整な `z` は `adjoin R {gen}` に
入る**」——つまり `adjoin R {gen} = integralClosure R L` が出る。 -/

/-- **Eisenstein なら生成元が整閉包全体を生成する(DVR上、核心の一歩)**:
`R` が DVR、`PB.gen` の(`R` 上の)最小多項式が `π` を基準に
Eisenstein なら、`L` の**任意の整な元** `z` は `Algebra.adjoin R
{PB.gen}` に入る。`ABC3.Found.GenEll.mem_adjoin_of_pow_smul_of_
isEisensteinAt`(帰納の芯)と `ABC3.Found.GenEll.exists_smul_mem_
adjoin_powerBasis`(分母を払う)を、`d=π^k×単元` という DVR 特有の
分解(`IsDiscreteValuationRing.eq_unit_mul_pow_irreducible`)で
橋渡しする。 -/
theorem adjoin_eq_top_of_isEisensteinAt {R K L : Type*}
    [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Field K] [Field L] [Algebra K L] [Algebra R L] [Algebra R K]
    [IsScalarTower R K L] [IsFractionRing R K]
    {PB : PowerBasis K L} (hBint : IsIntegral R PB.gen)
    {pi : R} (hpi : Irreducible pi)
    (hEis : (minpoly R PB.gen).IsEisensteinAt (Ideal.span {pi})) :
    ∀ z : L, IsIntegral R z → z ∈ Algebra.adjoin R {PB.gen} := by
  intro z hz
  obtain ⟨d, hd0, hdz⟩ := ABC3.Found.GenEll.exists_smul_mem_adjoin_powerBasis (R := R) PB z
  obtain ⟨n, u, hu⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hd0 hpi
  have hu2 : d = pi ^ n * (u : R) := by rw [hu]; ring
  have huz_int : IsIntegral R ((u : R) • z) := by
    rw [Algebra.smul_def]
    exact (isIntegral_algebraMap).mul hz
  have hdz' : pi ^ n • ((u : R) • z) ∈ Algebra.adjoin R {PB.gen} := by
    rw [smul_smul, ← hu2]; exact hdz
  have huz_mem : (u : R) • z ∈ Algebra.adjoin R {PB.gen} :=
    ABC3.Found.GenEll.mem_adjoin_of_pow_smul_of_isEisensteinAt hpi.prime hBint hEis n
      ((u : R) • z) huz_int hdz'
  have h2 : ((u⁻¹ : Rˣ) : R) • ((u : R) • z) = z := by
    rw [smul_smul]
    simp
  rw [← h2]
  exact Subalgebra.smul_mem _ huz_mem _

/-- **`adjoin_eq_top_of_isEisensteinAt` の等式版**: `adjoin R {PB.gen}
= integralClosure R L`(Subalgebra として)。`≤` は生成元が整である
ことから、`≥` が上の定理そのもの。 -/
theorem adjoin_eq_integralClosure_of_isEisensteinAt {R K L : Type*}
    [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Field K] [Field L] [Algebra K L] [Algebra R L] [Algebra R K]
    [IsScalarTower R K L] [IsFractionRing R K]
    {PB : PowerBasis K L} (hBint : IsIntegral R PB.gen)
    {pi : R} (hpi : Irreducible pi)
    (hEis : (minpoly R PB.gen).IsEisensteinAt (Ideal.span {pi})) :
    Algebra.adjoin R ({PB.gen} : Set L) = integralClosure R L := by
  apply le_antisymm
  · rw [Algebra.adjoin_le_iff]
    intro x hx
    simp only [Set.mem_singleton_iff] at hx
    subst hx
    exact hBint
  · intro z hz
    exact adjoin_eq_top_of_isEisensteinAt hBint hpi hEis z hz

/-- **`falt1CokernelLengthEq`/`cancel_conductor_delta` が要求する形**:
`W := integralClosure R L` の中で、`w := PB.gen`(整閉包の元として)が
`Algebra.adjoin R {w} = ⊤` を満たす。`Subalgebra.map_injective`(単射な
`AlgHom` による像の単射性)+ `AlgHom.map_adjoin_singleton` +
`Subalgebra.range_val` で、`L` の中での等式(上の定理)を `W` 自身の
中での「`⊤` 生成」に変換する。 -/
theorem adjoin_eq_top_in_integralClosure_of_isEisensteinAt {R K L : Type*}
    [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [Field K] [Field L] [Algebra K L] [Algebra R L] [Algebra R K]
    [IsScalarTower R K L] [IsFractionRing R K]
    {PB : PowerBasis K L} (hBint : IsIntegral R PB.gen)
    {pi : R} (hpi : Irreducible pi)
    (hEis : (minpoly R PB.gen).IsEisensteinAt (Ideal.span {pi})) :
    Algebra.adjoin R ({(⟨PB.gen, hBint⟩ : integralClosure R L)} : Set (integralClosure R L)) = ⊤ := by
  have hS : Algebra.adjoin R ({PB.gen} : Set L) = integralClosure R L :=
    adjoin_eq_integralClosure_of_isEisensteinAt hBint hpi hEis
  apply Subalgebra.map_injective (f := (integralClosure R L).val) Subtype.val_injective
  rw [AlgHom.map_adjoin_singleton, Algebra.map_top, Subalgebra.range_val]
  show Algebra.adjoin R ({PB.gen} : Set L) = integralClosure R L
  exact hS

/-- **「非常に分岐した」`V_n` 族の1段の単項生成性(完成)**:
`isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C` の `V_1 :=
integralClosure V0 L` に対し、`AdjoinRoot.root fK` を整閉包の元とみた
`w` が `Algebra.adjoin V0 {w} = ⊤` を満たす。`minpoly V0 w = f`
(`X^n-π` そのもの、`minpoly.isIntegrallyClosed_eq_field_fractions'` +
`Polynomial.map_injective` で `K` 上の minpoly から降ろす)を経由して
`adjoin_eq_top_in_integralClosure_of_isEisensteinAt` を適用する。 -/
theorem adjoin_eq_top_integralClosure_adjoinRoot_X_pow_sub_C
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n) :
    ∃ w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))),
      Algebra.adjoin V0 ({w} : Set _) = ⊤ := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set fK : Polynomial (FractionRing V0) := f.map (algebraMap V0 (FractionRing V0)) with hfKdef
  have hirr : Irreducible fK :=
    eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirr⟩
  have hmonicK : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  set PB := AdjoinRoot.powerBasis hirr.ne_zero with hPBdef
  have hPBgen : PB.gen = AdjoinRoot.root fK := AdjoinRoot.powerBasis_gen hirr.ne_zero
  have hBint : IsIntegral V0 PB.gen := by
    rw [hPBgen]
    have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
    refine ⟨f, hmonic, ?_⟩
    show Polynomial.aeval (AdjoinRoot.root fK) f = 0
    rw [← Polynomial.aeval_map_algebraMap (A := FractionRing V0)]
    rw [Polynomial.aeval_def, AdjoinRoot.algebraMap_eq]
    exact AdjoinRoot.eval₂_root fK
  have hminpolyK : minpoly (FractionRing V0) PB.gen = fK := by
    rw [hPBgen]
    have h := AdjoinRoot.minpoly_root hirr.ne_zero (f := fK)
    rw [hmonicK.leadingCoeff, inv_one, Polynomial.C_1, mul_one] at h
    exact h
  have hEqmin : minpoly V0 PB.gen = f := by
    have hEq : minpoly (FractionRing V0) PB.gen = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) :=
      minpoly.isIntegrallyClosed_eq_field_fractions' (FractionRing V0) hBint
    rw [hminpolyK] at hEq
    have hEq2 : f.map (algebraMap V0 (FractionRing V0)) = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) := hEq
    exact (Polynomial.map_injective _ (IsFractionRing.injective V0 (FractionRing V0)) hEq2).symm
  have hEis : (minpoly V0 PB.gen).IsEisensteinAt (Ideal.span ({π} : Set V0)) := by
    rw [hEqmin]; exact isEisensteinAt_X_pow_sub_C π n hnpos hprime hnotsq
  have hπ0 : π ≠ 0 := by
    intro h; apply hπne0; rw [h]; simp
  have hprimeπ : Prime π := (Ideal.span_singleton_prime hπ0).mp hprime
  have hpiirr : Irreducible π := hprimeπ.irreducible
  refine ⟨⟨PB.gen, hBint⟩, ?_⟩
  exact adjoin_eq_top_in_integralClosure_of_isEisensteinAt hBint hpiirr hEis

/-!
## Theorem 1.2・3b/3c 接続: `V_1` の differentIdeal の具体的な値(完成、2026-09-04)

item 3b の `cancel_conductor_delta` が要求する `hspan_eq`(「conductor と
δₙ₊₁ の base change が一致する」)を検証する足場として、`differentIdeal
V0 V1` そのものを**具体的な式で**計算した——Lemma 1.1 自身の道具
`differentIdeal_eq_span_derivative`(`w` が単項生成なら
`differentIdeal V W = span{f'(w)}`)を、上で完成した monogenicity
(`adjoin_eq_top_in_integralClosure_of_isEisensteinAt`)と組み合わせる
だけで閉じた。

★配管上の教訓(2026-09-04): `differentIdeal A B` は `[IsDedekindDomain
B]`・`[Module.IsTorsionFree A B]` を引数に持つ——これらを結論の中で
`integralClosure V0 (AdjoinRoot fK)` に対して**内部で `infer_instance`
により導出しようとすると**、定理の**主張(goal)自体**が同じ2つの
instance を要求するために、goal 側の(先に走る)instance 探索が
失敗した結果が**キャッシュされ**、後から `haveI`/`infer_instance` で
正しく揃えても失敗し続けるという罠に嵌った(`set` 起因ではなく、
「主張に現れる instance 引数」と「証明内で後から用意する instance」の
順序の罠——tools/lean-idioms.md に追記の価値がある新パターン)。
解決策は、これら2つを**定理の引数として明示的に要求する**
(`[IsDedekindDomain (...)]`・`[Module.IsTorsionFree V0 (...)]`)ことで
goal 自体の instance 探索を迂回すること。実用上は
`isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C`(前段で証明
済み)と `IsIntegralClosure.isTorsionFree`(mathlib)を呼び出し元で
`haveI` すれば、この2つは自動的に埋まる。 -/

set_option maxHeartbeats 1000000 in
/-- **`V_1 := integralClosure V0 (AdjoinRoot(X^n-π の base change))` の
differentIdeal は `span{n·w^{n-1}}`**(`w` は根の像)。古典的な
「Eisenstein 拡大の判別式(discriminant)」の公式そのもの
(`differentIdeal` は `discriminant` の平方根に相当する量)。
`minpoly V0 w = X^n-π`(前段で確立済み)の微分 `n·X^{n-1}` を `w` に
代入するだけ——`differentIdeal_eq_span_derivative` と
monogenicity(`adjoin_eq_top_in_integralClosure_of_isEisensteinAt`)を
貼り合わせる。`set_option maxHeartbeats`は、結論に現れる `differentIdeal`
の instance 引数由来の重い探索を許容するため(上の教訓参照)。 -/
theorem differentIdeal_eq_span_of_adjoinRoot_X_pow_sub_C
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))] :
    ∃ w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))),
      differentIdeal V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
        = Ideal.span {(n : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
            (algebraMap V0 (FractionRing V0))))) * w ^ (n - 1)} := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set fK : Polynomial (FractionRing V0) := f.map (algebraMap V0 (FractionRing V0)) with hfKdef
  have hirr : Irreducible fK :=
    eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirr⟩
  have hmonicK : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  have hmapeq : fK = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by
    rw [hfKdef, hfdef]; simp
  have hsepK : fK.Separable := by rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot fK) := hmonicK.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot fK) := ‹Module.Finite _ _›
  haveI : Algebra.IsSeparable (FractionRing V0) (AdjoinRoot fK) :=
    algIsSeparable_adjoinRoot_of_separable _ hmonicK hsepK
  set PB := AdjoinRoot.powerBasis hirr.ne_zero with hPBdef
  have hPBgen : PB.gen = AdjoinRoot.root fK := AdjoinRoot.powerBasis_gen hirr.ne_zero
  have hBint : IsIntegral V0 PB.gen := by
    rw [hPBgen]
    have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
    refine ⟨f, hmonic, ?_⟩
    show Polynomial.aeval (AdjoinRoot.root fK) f = 0
    rw [← Polynomial.aeval_map_algebraMap (A := FractionRing V0)]
    rw [Polynomial.aeval_def, AdjoinRoot.algebraMap_eq]
    exact AdjoinRoot.eval₂_root fK
  have hminpolyK : minpoly (FractionRing V0) PB.gen = fK := by
    rw [hPBgen]
    have h := AdjoinRoot.minpoly_root hirr.ne_zero (f := fK)
    rw [hmonicK.leadingCoeff, inv_one, Polynomial.C_1, mul_one] at h
    exact h
  have hEqmin : minpoly V0 PB.gen = f := by
    have hEq : minpoly (FractionRing V0) PB.gen = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) :=
      minpoly.isIntegrallyClosed_eq_field_fractions' (FractionRing V0) hBint
    rw [hminpolyK] at hEq
    have hEq2 : f.map (algebraMap V0 (FractionRing V0)) = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) := hEq
    exact (Polynomial.map_injective _ (IsFractionRing.injective V0 (FractionRing V0)) hEq2).symm
  set w : integralClosure V0 (AdjoinRoot fK) := ⟨PB.gen, hBint⟩ with hwdef
  have hwL : (w : AdjoinRoot fK) = PB.gen := rfl
  have hminpolyw : minpoly V0 w = minpoly V0 PB.gen := by
    rw [← hwL]
    exact (minpoly.algebraMap_eq (A := V0) (Subtype.val_injective) w).symm
  have hadjoin : Algebra.adjoin V0 ({w} : Set (integralClosure V0 (AdjoinRoot fK))) = ⊤ := by
    have hEis : (minpoly V0 PB.gen).IsEisensteinAt (Ideal.span ({π} : Set V0)) := by
      rw [hEqmin]; exact isEisensteinAt_X_pow_sub_C π n hnpos hprime hnotsq
    have hπ0 : π ≠ 0 := fun h => hπne0 (by rw [h]; simp)
    have hprimeπ : Prime π := (Ideal.span_singleton_prime hπ0).mp hprime
    exact adjoin_eq_top_in_integralClosure_of_isEisensteinAt hBint hprimeπ.irreducible hEis
  have hw : Algebra.adjoin (FractionRing V0)
      ({(algebraMap (integralClosure V0 (AdjoinRoot fK)) (AdjoinRoot fK)) w} : Set (AdjoinRoot fK)) = ⊤ := by
    show Algebra.adjoin (FractionRing V0) ({(w:AdjoinRoot fK)} : Set (AdjoinRoot fK)) = ⊤
    rw [hwL, hPBgen]
    exact AdjoinRoot.adjoinRoot_eq_top
  have hdiff := differentIdeal_eq_span_derivative w hw hadjoin
  refine ⟨w, ?_⟩
  rw [hdiff, hminpolyw, hEqmin, hfdef]
  congr 1
  have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0)
      = Polynomial.C (n : V0) * Polynomial.X ^ (n-1) := by
    rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
  rw [hderiv]
  simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]

/-!
## Theorem 1.2・3c → 3b 接続の第一歩: `V_n → W_n` 塔の base change 写像(2026-09-04)

`hspan_eq`(`cancel_conductor_delta`)を検証するには、`differentIdeal_tower_
diamond` の `Vₙ→Vₙ₊₁` と `Vₙ→Wₙ→Wₙ₊₁` が**同じ生成元を共有する**必要がある
——`Wₙ₊₁` を `Vₙ₊₁` と同じ多項式 `X^n-π` を `Wₙ` へ base change して構成
すれば、`AdjoinRoot fK →ₐ[V0] AdjoinRoot (fK.map φ0)`(`φ0` は係数体の
base change)という橋渡しの写像が要る。以下はその最初の(一般的な)
部品——具体的な Eisenstein 多項式にはまだ特化していない。 -/

/-- **`AdjoinRoot` の係数体の base change に沿った `V0`-代数準同型**:
`φ0 : K →+* K'` が `V0` 上の algebraMap と両立する(`hφ0`)なら、
`AdjoinRoot fK →ₐ[V0] AdjoinRoot (fK.map φ0)` が(`AdjoinRoot.map` を
`AlgHom.mk'` で `V0`-線形性を確認して)得られる。`K`・`K'` が `V0`
上の**異なる**代数(例: `FractionRing V0`・`FractionRing Wn`)である
一般の場合を扱う——`AdjoinRoot.mapAlgHom`(mathlib)は同じ底環上の
場合のみなので使えない。 -/
noncomputable def algHomAdjoinRootOfCompat {V0 K K' : Type*} [CommRing V0] [Field K] [Field K']
    [Algebra V0 K] [Algebra V0 K'] (φ0 : K →+* K')
    (hφ0 : ∀ r : V0, φ0 (algebraMap V0 K r) = algebraMap V0 K' r) (fK : Polynomial K) :
    AdjoinRoot fK →ₐ[V0] AdjoinRoot (fK.map φ0) :=
  AlgHom.mk' (AdjoinRoot.map φ0 fK (fK.map φ0) (dvd_refl _)) (by
    intro c x
    simp only [Algebra.smul_def]
    rw [map_mul]
    congr 1
    rw [IsScalarTower.algebraMap_apply V0 K (AdjoinRoot fK),
      IsScalarTower.algebraMap_apply V0 K' (AdjoinRoot (fK.map φ0)),
      AdjoinRoot.algebraMap_eq, AdjoinRoot.algebraMap_eq,
      AdjoinRoot.map_of, hφ0])

/-- `algHomAdjoinRootOfCompat` は根を根に送る——`w`(`V1` の生成元)の像が、
`Wₙ₊₁` 側で「同じ構成をした」根 `x` に一致することの鍵になる定義性質。 -/
theorem algHomAdjoinRootOfCompat_root {V0 K K' : Type*} [CommRing V0] [Field K] [Field K']
    [Algebra V0 K] [Algebra V0 K'] (φ0 : K →+* K')
    (hφ0 : ∀ r : V0, φ0 (algebraMap V0 K r) = algebraMap V0 K' r) (fK : Polynomial K) :
    algHomAdjoinRootOfCompat φ0 hφ0 fK (AdjoinRoot.root fK) = AdjoinRoot.root (fK.map φ0) :=
  AdjoinRoot.map_root φ0 fK (fK.map φ0) (dvd_refl _)

/-- **`algHomAdjoinRootOfCompat` を実際に使う `φ0`**:
`V0 → Wn` が単射なら、`IsFractionRing.map` で `FractionRing V0 →+*
FractionRing Wn` が作れる——これが `algHomAdjoinRootOfCompat` に渡す
`φ0` の具体形。 -/
noncomputable def fractionRingMapOfInjective {V0 Wn : Type*} [CommRing V0] [IsDomain V0]
    [CommRing Wn] [IsDomain Wn] [Algebra V0 Wn] (hinj : Function.Injective (algebraMap V0 Wn)) :
    FractionRing V0 →+* FractionRing Wn :=
  IsFractionRing.map (B := Wn) hinj

/-- `fractionRingMapOfInjective` は `V0` 上の algebraMap と両立する——
`algHomAdjoinRootOfCompat` の `hφ0` としてそのまま渡せる。
`IsLocalization.map_eq`(局所化の写像の定義性質)+
`IsScalarTower.algebraMap_apply`(`V0→Wn→Frac(Wn)` の合成が
`V0→Frac(Wn)` と一致する)で閉じる。 -/
theorem fractionRingMapOfInjective_algebraMap {V0 Wn : Type*} [CommRing V0] [IsDomain V0]
    [CommRing Wn] [IsDomain Wn] [Algebra V0 Wn] (hinj : Function.Injective (algebraMap V0 Wn))
    (r : V0) :
    fractionRingMapOfInjective hinj (algebraMap V0 (FractionRing V0) r)
      = algebraMap V0 (FractionRing Wn) r := by
  unfold fractionRingMapOfInjective IsFractionRing.map
  rw [IsLocalization.map_eq, ← IsScalarTower.algebraMap_apply]

/-!
## `V_1 ≃ AdjoinRoot(X^n-π)` そのもの(2026-09-04、驚きの単純化)

`hspan_eq` を接続するには `Algebra V1 Wn1` が実際の instance として要る
——`integralClosure` の一般の functoriality(`AlgHom.mapIntegralClosure`
+ `IsIntegral.tower_top` の合成)を経由するより、**`V1` が実は
`AdjoinRoot(X^n-π)`(元の多項式、`FractionRing` へ base change する
「前」の形)そのものと同型**であることを示す方が遥かに単純だと判明した
——`minpoly.equivAdjoin`(`AdjoinRoot(minpoly R x) ≃ₐ[R] adjoin R {x}`)
と monogenicity(`adjoin_eq_integralClosure_of_isEisensteinAt`、`adjoin
= integralClosure`)を貼り合わせるだけ。これで `Vₙ→Vₙ₊₁` の橋渡しは
`AdjoinRoot.map`(`V0` 上の**多項式環**の間の写像、体を経由しない)
だけで済み、`algHomAdjoinRootOfCompat`(フラクション体を経由する版)
すら不要になる可能性がある——次回はこちらを軸に組み立て直すとよい。 -/

set_option maxHeartbeats 400000 in
/-- **`V_1 := integralClosure V0 (AdjoinRoot(X^n-π の base change))` は
`AdjoinRoot(X^n-π)`(base change する前の、`V0[X]` 上そのままの
多項式の商)と `V0`-代数として同型**。`minpoly.equivAdjoin` +
monogenicity の合成。 -/
noncomputable def falt1AdjoinRootEquivIntegralClosure
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))] :
    AdjoinRoot (Polynomial.X ^ n - Polynomial.C π : Polynomial V0) ≃ₐ[V0]
      integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set fK : Polynomial (FractionRing V0) := f.map (algebraMap V0 (FractionRing V0)) with hfKdef
  have hirr : Irreducible fK :=
    eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirr⟩
  have hmonicK : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  set PB := AdjoinRoot.powerBasis hirr.ne_zero with hPBdef
  have hPBgen : PB.gen = AdjoinRoot.root fK := AdjoinRoot.powerBasis_gen hirr.ne_zero
  have hBint : IsIntegral V0 PB.gen := by
    rw [hPBgen]
    have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
    refine ⟨f, hmonic, ?_⟩
    show Polynomial.aeval (AdjoinRoot.root fK) f = 0
    rw [← Polynomial.aeval_map_algebraMap (A := FractionRing V0)]
    rw [Polynomial.aeval_def, AdjoinRoot.algebraMap_eq]
    exact AdjoinRoot.eval₂_root fK
  have hminpolyK : minpoly (FractionRing V0) PB.gen = fK := by
    rw [hPBgen]
    have h := AdjoinRoot.minpoly_root hirr.ne_zero (f := fK)
    rw [hmonicK.leadingCoeff, inv_one, Polynomial.C_1, mul_one] at h
    exact h
  have hEqmin : minpoly V0 PB.gen = f := by
    have hEq : minpoly (FractionRing V0) PB.gen = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) :=
      minpoly.isIntegrallyClosed_eq_field_fractions' (FractionRing V0) hBint
    rw [hminpolyK] at hEq
    have hEq2 : f.map (algebraMap V0 (FractionRing V0)) = (minpoly V0 PB.gen).map (algebraMap V0 (FractionRing V0)) := hEq
    exact (Polynomial.map_injective _ (IsFractionRing.injective V0 (FractionRing V0)) hEq2).symm
  have hπ0 : π ≠ 0 := fun h => hπne0 (by rw [h]; simp)
  have hprimeπ : Prime π := (Ideal.span_singleton_prime hπ0).mp hprime
  have hEis : (minpoly V0 PB.gen).IsEisensteinAt (Ideal.span ({π} : Set V0)) := by
    rw [hEqmin]; exact isEisensteinAt_X_pow_sub_C π n hnpos hprime hnotsq
  have hadjoinEq : Algebra.adjoin V0 ({PB.gen} : Set (AdjoinRoot fK)) = integralClosure V0 (AdjoinRoot fK) :=
    adjoin_eq_integralClosure_of_isEisensteinAt hBint hprimeπ.irreducible hEis
  have hL : Module.IsTorsionFree V0 (AdjoinRoot fK) := by
    have hinj : Function.Injective (algebraMap V0 (AdjoinRoot fK)) := by
      rw [IsScalarTower.algebraMap_eq V0 (FractionRing V0) (AdjoinRoot fK)]
      exact (algebraMap (FractionRing V0) (AdjoinRoot fK)).injective.comp
        (IsFractionRing.injective V0 (FractionRing V0))
    refine ⟨fun r hr => ?_⟩
    have hr0 : r ≠ 0 := by
      rintro rfl
      exact zero_ne_one (hr.left (show (0:V0) * 0 = 0 * 1 by simp))
    have hr0' : algebraMap V0 (AdjoinRoot fK) r ≠ 0 := fun h => hr0 (hinj (by rw [h, map_zero]))
    intro x y hxy
    simp only [Algebra.smul_def] at hxy
    exact mul_left_cancel₀ hr0' hxy
  haveI := hL
  have e1 := minpoly.equivAdjoin hBint
  rw [hEqmin] at e1
  rw [hPBgen] at hadjoinEq
  exact e1.trans (Subalgebra.equivOfEq _ _ hadjoinEq)

/-- **`AdjoinRoot` の底環 base change に沿った `V0`-代数準同型
(体を経由しない版)**: `falt1AdjoinRootEquivIntegralClosure` の発見
(`V1 ≃ AdjoinRoot(X^n-π)`)を受けて、`algHomAdjoinRootOfCompat`
(フラクション体 `K→+*K'` を経由する版)より遥かに単純な経路が
見つかった——`V0[X]→Wn[X]` の底環 base change をそのまま使うだけ
(フラクション体・`IsFractionRing.map` は一切不要)。 -/
noncomputable def algHomAdjoinRootOfCompat' {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn]
    (f : Polynomial V0) :
    AdjoinRoot f →ₐ[V0] AdjoinRoot (f.map (algebraMap V0 Wn)) :=
  AlgHom.mk' (AdjoinRoot.map (algebraMap V0 Wn) f (f.map (algebraMap V0 Wn)) (dvd_refl _)) (by
    intro c x
    simp only [Algebra.smul_def]
    rw [map_mul]
    congr 1
    rw [AdjoinRoot.algebraMap_eq]
    rw [IsScalarTower.algebraMap_apply V0 Wn (AdjoinRoot (f.map (algebraMap V0 Wn))),
      AdjoinRoot.algebraMap_eq, AdjoinRoot.map_of])

/-- `algHomAdjoinRootOfCompat'` も根を根に送る。 -/
theorem algHomAdjoinRootOfCompat'_root {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn]
    (f : Polynomial V0) :
    algHomAdjoinRootOfCompat' (Wn := Wn) f (AdjoinRoot.root f) = AdjoinRoot.root (f.map (algebraMap V0 Wn)) :=
  AdjoinRoot.map_root (algebraMap V0 Wn) f (f.map (algebraMap V0 Wn)) (dvd_refl _)

/-!
## `adjoinRootTensorEquiv` の点ごとの挙動・別経路での構成(2026-09-04)

`adjoinRootTensorEquiv`(`Algebra.TensorProduct.tensorQuotientEquiv`+
`Ideal.quotientEquivAlg`経由)の点ごとの挙動(`1⊗ₜ(AdjoinRoot.mk g p)`
の像)を計算しようとすると、`AdjoinRoot h`と`Polynomial _⧸Ideal.span{h}`
が定義的には同じ型でも**別々に定義された`Semiring`/`Algebra`インスタンス**
(`AdjoinRoot.instAlgebra`系 vs `Ideal.Quotient.ring`系)を持つという
mathlib自体の instance diamond に当たり、`unfold`後の`rw`が
「target expression is not type-correct」で失敗することが分かった。

**回避策**: `Ideal.Quotient`側のAPIを一切経由せず、`AdjoinRoot`
ネイティブのAPI(`TensorProduct.AlgebraTensorModule.lift`・
`algHomAdjoinRootOfCompat'`・`AdjoinRoot.lift_mk`・`Polynomial.
induction_on`)だけで**別経路の**線形写像を構成し、その点ごとの挙動を
証明した。この写像は`adjoinRootTensorEquiv`と同じ役割(`C⊗[R]AdjoinRoot
g → AdjoinRoot(g.map φ)`)を果たすはずだが、**全単射性はまだ示して
いない**(次回の課題、falt1-goal.md参照)。 -/

/-- `eval₂ (AdjoinRoot.of q) (AdjoinRoot.root q) h = AdjoinRoot.mk q h`
(`AdjoinRoot.map`の点ごとの挙動を計算する鍵、`Polynomial.induction_on`
による3行の帰納法)。 -/
theorem eval2_root_eq_mk {R : Type*} [CommRing R] (q h : Polynomial R) :
    Polynomial.eval₂ (AdjoinRoot.of q) (AdjoinRoot.root q) h = AdjoinRoot.mk q h := by
  induction h using Polynomial.induction_on with
  | C a => rw [Polynomial.eval₂_C, AdjoinRoot.mk_C]
  | add p1 p2 h1 h2 => rw [Polynomial.eval₂_add, h1, h2, map_add]
  | monomial n a ih =>
      rw [pow_succ, ← mul_assoc, Polynomial.eval₂_mul, Polynomial.eval₂_X, map_mul, AdjoinRoot.mk_X, ih]

/-- **`adjoinRootTensorEquiv`と同じ役割を果たす、`Ideal.Quotient`の
instance diamondを回避した線形写像**: `TensorProduct.AlgebraTensorModule.
lift`で`(c,x) ↦ c • algHomAdjoinRootOfCompat' g x`という双線形写像を
持ち上げるだけ。全単射性はまだ示していない。 -/
noncomputable def adjoinRootTensorEquivFwd {R C : Type*} [CommRing R] [CommRing C] [Algebra R C]
    (g : Polynomial R) :
    TensorProduct R C (AdjoinRoot g) →ₗ[C] AdjoinRoot (g.map (algebraMap R C)) := by
  apply TensorProduct.AlgebraTensorModule.lift
  exact {
    toFun := fun c => c • (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g).toLinearMap
    map_add' := by intro x y; ext z; simp [add_smul]
    map_smul' := by intro x y; ext z; simp [mul_smul] }

/-- **`adjoinRootTensorEquivFwd`の点ごとの挙動(完成)**:
`AdjoinRoot.lift_mk`(`AdjoinRoot.map`の内部構成)と`eval2_root_eq_mk`
を組み合わせるだけ——`Ideal.Quotient`のinstance diamondを一切経由
しない。 -/
theorem adjoinRootTensorEquivFwd_one_tmul_mk {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g p : Polynomial R) :
    (adjoinRootTensorEquivFwd (R := R) (C := C) g) ((1:C) ⊗ₜ[R] (AdjoinRoot.mk g p)) =
      AdjoinRoot.mk (g.map (algebraMap R C)) (p.map (algebraMap R C)) := by
  show (1:C) • (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g) (AdjoinRoot.mk g p) = _
  rw [one_smul]
  show algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g (AdjoinRoot.mk g p) = _
  unfold algHomAdjoinRootOfCompat'
  show (AdjoinRoot.map (algebraMap R C) g (g.map (algebraMap R C)) (dvd_refl _)) (AdjoinRoot.mk g p) = _
  rw [show AdjoinRoot.map (algebraMap R C) g (g.map (algebraMap R C)) (dvd_refl _) =
    AdjoinRoot.lift ((AdjoinRoot.of (g.map (algebraMap R C))).comp (algebraMap R C)) (AdjoinRoot.root (g.map (algebraMap R C))) (by
      show Polynomial.eval₂ ((AdjoinRoot.of (g.map (algebraMap R C))).comp (algebraMap R C)) (AdjoinRoot.root (g.map (algebraMap R C))) g = 0
      rw [← Polynomial.eval₂_map]
      exact AdjoinRoot.eval₂_root (g.map (algebraMap R C))) from rfl]
  rw [AdjoinRoot.lift_mk, ← Polynomial.eval₂_map, eval2_root_eq_mk]

/-!
## `adjoinRootTensorEquivFwd` の全単射性(完成、2026-09-04)

`adjoinRootTensorEquivFwd`の**逆写像を明示的に構成**し、両方向の
往復が恒等写像になることを確認するだけで全単射性を示した——
`AdjoinRoot`ネイティブAPI(`AdjoinRoot.lift`・`AdjoinRoot.induction_on`・
`Polynomial.induction_on`)だけで閉じ、`Ideal.Quotient`側のinstance
diamondには一度も触れない。 -/

set_option maxHeartbeats 400000 in
/-- **`adjoinRootTensorEquivFwd`の逆写像候補**: `root(g.map φ) ↦
1⊗ₜ(root g)`・`C ↦ (c⊗ₜ1)`となる`AdjoinRoot.lift`。well-definedness
(`g.map φ`が`1⊗ₜroot g`で消えること)は`aeval`のAlgHom越しの自然性
(`Polynomial.aeval_algHom_apply`)+`AdjoinRoot.aeval_eq`(`aeval(root f)f
=mk f f=0`)だけで閉じる。 -/
noncomputable def adjoinRootTensorEquivInv {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    AdjoinRoot (g.map (algebraMap R C)) →+* TensorProduct R C (AdjoinRoot g) := by
  apply AdjoinRoot.lift (algebraMap C (TensorProduct R C (AdjoinRoot g))) ((1:C) ⊗ₜ[R] (AdjoinRoot.root g))
  show Polynomial.eval₂ (algebraMap C (TensorProduct R C (AdjoinRoot g))) ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) (g.map (algebraMap R C)) = 0
  rw [Polynomial.eval₂_map]
  rw [show ((algebraMap C (TensorProduct R C (AdjoinRoot g))).comp (algebraMap R C)) =
      (algebraMap R (TensorProduct R C (AdjoinRoot g))) from (IsScalarTower.algebraMap_eq R C (TensorProduct R C (AdjoinRoot g))).symm]
  show Polynomial.aeval ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) g = 0
  rw [show ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) =
      (Algebra.TensorProduct.includeRight (R := R) (A := C) (B := AdjoinRoot g)) (AdjoinRoot.root g) from rfl]
  rw [Polynomial.aeval_algHom_apply, AdjoinRoot.aeval_eq]
  simp [AdjoinRoot.mk_self]

set_option maxHeartbeats 400000 in
/-- `adjoinRootTensorEquivInv`の`C·X^n`単項式での挙動(帰納の要となる
明示形)。 -/
theorem adjoinRootTensorEquivInv_monomial {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R)
    (a : C) (n : ℕ) :
    (adjoinRootTensorEquivInv (R := R) (C := C) g) (AdjoinRoot.mk (g.map (algebraMap R C)) (Polynomial.C a * Polynomial.X ^ n)) =
      a ⊗ₜ[R] ((AdjoinRoot.root g) ^ n) := by
  show Polynomial.eval₂ (algebraMap C (TensorProduct R C (AdjoinRoot g))) ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) (Polynomial.C a * Polynomial.X ^ n) = _
  rw [Polynomial.eval₂_mul, Polynomial.eval₂_C, Polynomial.eval₂_X_pow]
  rw [show (algebraMap C (TensorProduct R C (AdjoinRoot g)) a) = a ⊗ₜ[R] (1 : AdjoinRoot g) from rfl]
  rw [show ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) ^ n = (1:C) ⊗ₜ[R] ((AdjoinRoot.root g) ^ n) by
    induction n with
    | zero => simp; rfl
    | succ k ih => rw [pow_succ, pow_succ, ih, Algebra.TensorProduct.tmul_mul_tmul, mul_one]]
  rw [Algebra.TensorProduct.tmul_mul_tmul, mul_one, one_mul]

set_option maxHeartbeats 400000 in
/-- `adjoinRootTensorEquivInv`の`mk(g.map φ)(p.map φ)`での挙動(逆向き
の往復で使う)。 -/
theorem adjoinRootTensorEquivInv_mk_map {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g p : Polynomial R) :
    (adjoinRootTensorEquivInv (R := R) (C := C) g) (AdjoinRoot.mk (g.map (algebraMap R C)) (p.map (algebraMap R C))) =
      (1:C) ⊗ₜ[R] (AdjoinRoot.mk g p) := by
  show Polynomial.eval₂ (algebraMap C (TensorProduct R C (AdjoinRoot g))) ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) (p.map (algebraMap R C)) = _
  rw [Polynomial.eval₂_map]
  rw [show ((algebraMap C (TensorProduct R C (AdjoinRoot g))).comp (algebraMap R C)) =
      (algebraMap R (TensorProduct R C (AdjoinRoot g))) from (IsScalarTower.algebraMap_eq R C (TensorProduct R C (AdjoinRoot g))).symm]
  show Polynomial.aeval ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) p = _
  rw [show ((1:C) ⊗ₜ[R] (AdjoinRoot.root g)) =
      (Algebra.TensorProduct.includeRight (R := R) (A := C) (B := AdjoinRoot g)) (AdjoinRoot.root g) from rfl]
  rw [Polynomial.aeval_algHom_apply, AdjoinRoot.aeval_eq]
  rfl

/-- `adjoinRootTensorEquivInv`は`C`線形(algebraMap Cとの両立性から)。 -/
theorem adjoinRootTensorEquivInv_smul {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R)
    (c : C) (y : AdjoinRoot (g.map (algebraMap R C))) :
    (adjoinRootTensorEquivInv (R := R) (C := C) g) (c • y) = c • (adjoinRootTensorEquivInv (R := R) (C := C) g) y := by
  rw [Algebra.smul_def, map_mul]
  rw [show (adjoinRootTensorEquivInv (R := R) (C := C) g) (algebraMap C (AdjoinRoot (g.map (algebraMap R C))) c) =
      algebraMap C (TensorProduct R C (AdjoinRoot g)) c by
    rw [show (algebraMap C (AdjoinRoot (g.map (algebraMap R C))) c) = AdjoinRoot.of (g.map (algebraMap R C)) c from
      (AdjoinRoot.algebraMap_eq' C (g.map (algebraMap R C))).symm ▸ rfl]
    exact AdjoinRoot.lift_of _]
  rw [← Algebra.smul_def]

set_option maxHeartbeats 400000 in
/-- **往復1: `fwd∘inv=id`**(`AdjoinRoot.induction_on`+`Polynomial.
induction_on`による帰納法、`C`・`add`・`monomial`の3ケース)。 -/
theorem adjoinRootTensorEquiv_roundtrip1 {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    ∀ y : AdjoinRoot (g.map (algebraMap R C)),
      (adjoinRootTensorEquivFwd (R := R) (C := C) g) ((adjoinRootTensorEquivInv (R := R) (C := C) g) y) = y := by
  intro y
  induction y using AdjoinRoot.induction_on with
  | ih p =>
    induction p using Polynomial.induction_on with
    | C a =>
        have h0 := adjoinRootTensorEquivInv_monomial (R := R) (C := C) g a 0
        simp only [pow_zero, mul_one] at h0
        rw [h0]
        rw [show (a ⊗ₜ[R] (1 : AdjoinRoot g) : TensorProduct R C (AdjoinRoot g)) =
            a • ((1:C) ⊗ₜ[R] (AdjoinRoot.mk g 1)) by
          rw [map_one]; rw [TensorProduct.smul_tmul']; simp]
        rw [map_smul, adjoinRootTensorEquivFwd_one_tmul_mk]
        rw [show ((1:Polynomial R)).map (algebraMap R C) = (1 : Polynomial C) by simp]
        rw [AdjoinRoot.mk_C]
        rw [Algebra.smul_def, AdjoinRoot.algebraMap_eq, map_one, mul_one]
    | add p1 p2 h1 h2 =>
        rw [map_add, map_add, map_add, h1, h2]
    | monomial n a _ih2 =>
        rw [pow_succ, ← mul_assoc, map_mul, map_mul, adjoinRootTensorEquivInv_monomial]
        have hInvX : (adjoinRootTensorEquivInv (R := R) (C := C) g) (AdjoinRoot.mk (g.map (algebraMap R C)) Polynomial.X) =
            (1:C) ⊗ₜ[R] (AdjoinRoot.root g) := by
          have h1 := adjoinRootTensorEquivInv_monomial (R := R) (C := C) g 1 1
          simpa using h1
        rw [hInvX, Algebra.TensorProduct.tmul_mul_tmul, mul_one, ← pow_succ]
        rw [show (a ⊗ₜ[R] ((AdjoinRoot.root g)^(n+1)) : TensorProduct R C (AdjoinRoot g)) =
            a • ((1:C) ⊗ₜ[R] (AdjoinRoot.mk g (Polynomial.X^(n+1)))) by
          rw [TensorProduct.smul_tmul']; simp [← AdjoinRoot.mk_X, ← map_pow]]
        rw [map_smul, adjoinRootTensorEquivFwd_one_tmul_mk, Polynomial.map_pow, Polynomial.map_X]
        rw [← map_mul, Algebra.smul_def, AdjoinRoot.algebraMap_eq]
        rw [← AdjoinRoot.mk_C, ← map_mul, mul_assoc, pow_succ]

set_option maxHeartbeats 400000 in
/-- **往復2: `inv∘fwd=id`**(`TensorProduct.induction_on`+
`AdjoinRoot.induction_on`による帰納法、`adjoinRootTensorEquivInv_smul`
で`C`線形性を経由)。 -/
theorem adjoinRootTensorEquiv_roundtrip2 {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    ∀ x : TensorProduct R C (AdjoinRoot g),
      (adjoinRootTensorEquivInv (R := R) (C := C) g) ((adjoinRootTensorEquivFwd (R := R) (C := C) g) x) = x := by
  have key : ∀ p : Polynomial R, (adjoinRootTensorEquivInv (R := R) (C := C) g)
      ((adjoinRootTensorEquivFwd (R := R) (C := C) g) ((1:C) ⊗ₜ[R] (AdjoinRoot.mk g p))) =
      (1:C) ⊗ₜ[R] (AdjoinRoot.mk g p) := by
    intro p
    rw [adjoinRootTensorEquivFwd_one_tmul_mk, adjoinRootTensorEquivInv_mk_map]
  intro x
  induction x using TensorProduct.induction_on with
  | zero => simp
  | tmul c z =>
      induction z using AdjoinRoot.induction_on with
      | ih p =>
          rw [show (c ⊗ₜ[R] (AdjoinRoot.mk g p) : TensorProduct R C (AdjoinRoot g)) =
              c • ((1:C) ⊗ₜ[R] (AdjoinRoot.mk g p)) by rw [TensorProduct.smul_tmul']; simp]
          rw [map_smul, adjoinRootTensorEquivInv_smul, key]
  | add x y hx hy => rw [map_add, map_add, hx, hy]

/-- **`adjoinRootTensorEquivFwd`は全単射(完成)**: 両方向の往復
(`adjoinRootTensorEquiv_roundtrip1`・`adjoinRootTensorEquiv_roundtrip2`)
から`Function.bijective_iff_has_inverse`で直ちに従う——これで`Wₙ₊₁`が
`Wₙ`と`AdjoinRoot f`(≃V1)の`V0`上のpushoutであること
(`Algebra.IsPushout`)への接続にNear、`Ideal.Quotient`のinstance
diamondを一切経由せずに到達できる。 -/
theorem adjoinRootTensorEquivFwd_bijective {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    Function.Bijective (adjoinRootTensorEquivFwd (R := R) (C := C) g) :=
  Function.bijective_iff_has_inverse.mpr
    ⟨adjoinRootTensorEquivInv (R := R) (C := C) g, adjoinRootTensorEquiv_roundtrip2 g, adjoinRootTensorEquiv_roundtrip1 g⟩

/-!
## `Algebra.IsPushout` への接続が完成した(2026-09-04)

`Algebra.IsPushout R S R' S'` は mathlib で単に
`IsBaseChange S (IsScalarTower.toAlgHom R R' S').toLinearMap` の
ラッパー(`Algebra.IsPushout.out`)であると判明した——つまり
`adjoinRootTensorEquivFwd_bijective` から直接
`IsBaseChange`を経由して`Algebra.IsPushout`に落とせる、迂回路無しの
最短経路。 -/

/-- `algHomAdjoinRootOfCompat'` の `AdjoinRoot.mk` での挙動(点ごとの
完全な計算、`adjoinRootTensorEquivFwd_one_tmul_mk` の証明の後半部分を
単独の補題として抽出したもの)。 -/
theorem algHomAdjoinRootOfCompat'_mk {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g p : Polynomial R) :
    algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g (AdjoinRoot.mk g p) =
      AdjoinRoot.mk (g.map (algebraMap R C)) (p.map (algebraMap R C)) := by
  unfold algHomAdjoinRootOfCompat'
  show (AdjoinRoot.map (algebraMap R C) g (g.map (algebraMap R C)) (dvd_refl _)) (AdjoinRoot.mk g p) = _
  rw [show AdjoinRoot.map (algebraMap R C) g (g.map (algebraMap R C)) (dvd_refl _) =
    AdjoinRoot.lift ((AdjoinRoot.of (g.map (algebraMap R C))).comp (algebraMap R C)) (AdjoinRoot.root (g.map (algebraMap R C))) (by
      show Polynomial.eval₂ ((AdjoinRoot.of (g.map (algebraMap R C))).comp (algebraMap R C)) (AdjoinRoot.root (g.map (algebraMap R C))) g = 0
      rw [← Polynomial.eval₂_map]
      exact AdjoinRoot.eval₂_root (g.map (algebraMap R C))) from rfl]
  rw [AdjoinRoot.lift_mk, ← Polynomial.eval₂_map, eval2_root_eq_mk]

/-- **`algHomAdjoinRootOfCompat'`(=`V1 → W1`の橋渡し写像)が`C`上の
base change である**: `adjoinRootTensorEquivFwd`の全単射性
(`adjoinRootTensorEquivFwd_bijective`)から`LinearEquiv.ofBijective`で
線形同型を作り、その`1⊗ₜ-`での挙動が`algHomAdjoinRootOfCompat'`と
一致すること(`adjoinRootTensorEquivFwd_one_tmul_mk`+
`algHomAdjoinRootOfCompat'_mk`)を`AdjoinRoot.induction_on`で確認する
だけで`IsBaseChange.of_equiv`が適用できる。 -/
theorem falt1_isBaseChange_adjoinRoot {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    IsBaseChange C (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g).toLinearMap := by
  apply IsBaseChange.of_equiv (LinearEquiv.ofBijective (adjoinRootTensorEquivFwd (R := R) (C := C) g)
    (adjoinRootTensorEquivFwd_bijective g))
  intro y
  induction y using AdjoinRoot.induction_on with
  | ih p =>
    show adjoinRootTensorEquivFwd (R := R) (C := C) g ((1:C) ⊗ₜ[R] AdjoinRoot.mk g p) = _
    rw [adjoinRootTensorEquivFwd_one_tmul_mk]
    exact (algHomAdjoinRootOfCompat'_mk g p).symm

/-- **`AdjoinRoot(g.map φ)`を`AdjoinRoot g`-代数と見る自然な構造**:
橋渡し写像`algHomAdjoinRootOfCompat' g : AdjoinRoot g →ₐ[R] AdjoinRoot
(g.map φ)`を`RingHom.toAlgebra`で`Algebra`インスタンスへ変換したもの。
(`instance`ではなく`def`のまま保つ——`Falt1`の具体的なオブジェクトに
適用する際、既存の`Algebra V1 Wn1`インスタンスとの衝突を避けるため。) -/
@[reducible]
noncomputable def falt1AdjoinRootAlgebra {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    Algebra (AdjoinRoot g) (AdjoinRoot (g.map (algebraMap R C))) :=
  (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g).toRingHom.toAlgebra

/-- `falt1AdjoinRootAlgebra`の下で`R → AdjoinRoot g → AdjoinRoot(g.map φ)`
が`IsScalarTower`をなす(`algHomAdjoinRootOfCompat'`が`R`上の`AlgHom`
であること、すなわち`.commutes`から直ちに従う)。 -/
theorem falt1AdjoinRoot_isScalarTower {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    letI := falt1AdjoinRootAlgebra (R := R) (C := C) g
    IsScalarTower R (AdjoinRoot g) (AdjoinRoot (g.map (algebraMap R C))) := by
  letI := falt1AdjoinRootAlgebra (R := R) (C := C) g
  apply IsScalarTower.of_algebraMap_eq
  intro x
  show algebraMap R (AdjoinRoot (g.map (algebraMap R C))) x =
    algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g (algebraMap R (AdjoinRoot g) x)
  exact ((algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g).commutes x).symm

/-- **本節の到達点**: `AdjoinRoot(g.map φ)`は`R`上、`C`と`AdjoinRoot g`
の pushout である(`Algebra.IsPushout`)。`Algebra.IsPushout`が
mathlib で単なる`IsBaseChange`のラッパーだと判明した
(`Algebra.IsPushout.out`)ため、`falt1_isBaseChange_adjoinRoot`から
`constructor`一発で閉じる——`Ideal.Quotient`のinstance diamondを
一度も経由していない。`pushoutKaehlerSplitStep`/`_length`が要求する
`[Algebra.IsPushout R B1 C B]`インスタンスに、`B1 = AdjoinRoot g`
(=`V1`)・`B = AdjoinRoot(g.map φ)`(=`Wn1`)として直接対応する。 -/
theorem falt1_isPushout_adjoinRoot {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    letI := falt1AdjoinRootAlgebra (R := R) (C := C) g
    letI := falt1AdjoinRoot_isScalarTower (R := R) (C := C) g
    Algebra.IsPushout R C (AdjoinRoot g) (AdjoinRoot (g.map (algebraMap R C))) := by
  letI := falt1AdjoinRootAlgebra (R := R) (C := C) g
  letI := falt1AdjoinRoot_isScalarTower (R := R) (C := C) g
  constructor
  show IsBaseChange C (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g).toLinearMap
  exact falt1_isBaseChange_adjoinRoot g

/-!
## `pushoutKaehlerSplitStep` を具体的な `AdjoinRoot` の反復pushoutに
繋ぐ(2026-09-05)

`Algebra.IsPushout V0 Wn V1 Wn1`(`Wₙ→Wₙ₊₁`の橋)は Theorem 1.2 の
証明原文(物理p.5第2文、falt1-goal.md 参照)自体がそのズレを測ることで
成り立っていることが判明し、一般には偽——**本来の適用先は`Vₙ→Vₙ₊₁`
の「`d+1`個の単一生成拡大の同時添加」の側**(falt1-goal.md 2026-09-05
の記録参照)。この節では、その本来の適用先で実際に必要になる2つの道具
(`Algebra.IsPushout`の具体形は既に`falt1_isPushout_adjoinRoot`で
完成済み)を揃える: (1) `pushoutKaehlerSplit`系が要求する`hinj`を
`AdjoinRoot(g.map φ)`の言葉で直接与える版、(2) 反復pushoutで`B1`を
次の段の`C`として使うときに必要な`AlgEquiv`(`mapBaseChange_injective_
transport`は`LinearEquiv`ではなく`AlgEquiv`を要求するため、
`adjoinRootTensorEquiv`のinstance diamondを再導入せずに用意する)。 -/

/-- `adjoinRootTensorEquivFwd`を`AlgHom`として再構成したもの
(`Algebra.TensorProduct.lift`で`(c,x)↦c*ψ(x)`を持ち上げるだけ、
`ψ=algHomAdjoinRootOfCompat'`)。 -/
noncomputable def adjoinRootTensorAlgHom {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    TensorProduct R C (AdjoinRoot g) →ₐ[C] AdjoinRoot (g.map (algebraMap R C)) :=
  Algebra.TensorProduct.lift (Algebra.ofId C (AdjoinRoot (g.map (algebraMap R C))))
    (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g) (fun _ _ => mul_comm _ _)

/-- `adjoinRootTensorAlgHom`の下請け関数は`adjoinRootTensorEquivFwd`と
(pure tensorで一致するので`TensorProduct.induction_on`により)
**至るところ**一致する。 -/
theorem adjoinRootTensorAlgHom_eq_fwd {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    ∀ z, (adjoinRootTensorAlgHom (R := R) (C := C) g) z = (adjoinRootTensorEquivFwd (R := R) (C := C) g) z := by
  intro z
  induction z using TensorProduct.induction_on with
  | zero => simp
  | tmul c x =>
      show algebraMap C (AdjoinRoot (g.map (algebraMap R C))) c * (algHomAdjoinRootOfCompat' (V0 := R) (Wn := C) g) x = _
      rw [← Algebra.smul_def]
      rfl
  | add x y hx hy => rw [map_add, map_add, hx, hy]

/-- 上の一致から、`adjoinRootTensorAlgHom`の全単射性が
`adjoinRootTensorEquivFwd_bijective`から直ちに従う。 -/
theorem adjoinRootTensorAlgHom_bijective {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    Function.Bijective (adjoinRootTensorAlgHom (R := R) (C := C) g) := by
  have := adjoinRootTensorEquivFwd_bijective (R := R) (C := C) g
  rwa [show (adjoinRootTensorAlgHom (R := R) (C := C) g : TensorProduct R C (AdjoinRoot g) → AdjoinRoot (g.map (algebraMap R C))) =
    (adjoinRootTensorEquivFwd (R := R) (C := C) g : TensorProduct R C (AdjoinRoot g) → AdjoinRoot (g.map (algebraMap R C))) from
    funext (adjoinRootTensorAlgHom_eq_fwd g)]

/-- **`adjoinRootTensorEquiv`(instance diamond有り)の代わりに使える
`AlgEquiv`**: `mapBaseChange_injective_transport`のような`AlgEquiv`を
要求する道具にそのまま渡せる、`Ideal.Quotient`を経由しない版。 -/
noncomputable def adjoinRootTensorAlgEquiv {R C : Type*} [CommRing R] [CommRing C] [Algebra R C] (g : Polynomial R) :
    TensorProduct R C (AdjoinRoot g) ≃ₐ[C] AdjoinRoot (g.map (algebraMap R C)) :=
  AlgEquiv.ofBijective (adjoinRootTensorAlgHom (R := R) (C := C) g) (adjoinRootTensorAlgHom_bijective g)

/-- **`pushoutKaehlerSplitStep`の`hinj`を`AdjoinRoot(g.map φ)`の言葉で
直接与える**: `mapBaseChange_injective_adjoinRoot_tensor`
(`TensorProduct`版、Lemma 1.1型の非零因子条件から)を
`adjoinRootTensorAlgEquiv`で`AdjoinRoot(g.map φ)`側へ移送するだけ。
これで`pushoutKaehlerSplitStep`の`B := AdjoinRoot(g.map(algebraMap R
B1))`という反復pushoutの各段で要る`hinj`が、Lemma 1.1と同じ「その段の
最小多項式の微分が非零因子」という1つの仮定に帰着する。 -/
theorem mapBaseChange_injective_adjoinRoot_direct {R C : Type*} [CommRing R] [CommRing C] [Algebra R C]
    (g : Polynomial R)
    (hnzd : algebraMap (Polynomial C) (AdjoinRoot (g.map (algebraMap R C)))
        (Polynomial.derivative (g.map (algebraMap R C)))
        ∈ nonZeroDivisors (AdjoinRoot (g.map (algebraMap R C)))) :
    Function.Injective (KaehlerDifferential.mapBaseChange R C (AdjoinRoot (g.map (algebraMap R C)))) :=
  mapBaseChange_injective_transport (adjoinRootTensorAlgEquiv (R := R) (C := C) g)
    (mapBaseChange_injective_adjoinRoot_tensor g hnzd)

/-!
## `V_1 → W_1` の橋渡し写像そのものが完成した(2026-09-04)

`hspan_eq` の接続に必要な `Algebra V1 Wn1` インスタンスを、実際に構成
できる `AlgHom` として手に入れた——`differentIdeal_tower_diamond` を
このインスタンスの下で呼び出す準備が整った。組み立ては:
`V1 ≃ₐ[V0] AdjoinRoot f`(`falt1AdjoinRootEquivIntegralClosure`、symm)
→ `AdjoinRoot f →ₐ[V0] AdjoinRoot g`(`algHomAdjoinRootOfCompat'`、
`g:=f.map(algebraMap V0 Wn)`、多項式として `X^n-Cπ'` に等しい)
→ `AdjoinRoot g ≃ₐ[Wn] Wn1`(`falt1AdjoinRootEquivIntegralClosure` を
`Wn`・`π'` に対して再適用、`AlgHom.restrictScalars V0` で `V0`-線形
写像とみなす)、の3段の合成。 -/

/-- **`V_1 → W_1` の橋渡し写像**: `Wn` が(`V0` と同じ形の)Eisenstein
拡大の基点になる場合、`V1 := integralClosure V0 (AdjoinRoot(X^n-π の
base change))` から `Wn1 := integralClosure Wn (AdjoinRoot(X^n-π' の
base change))`(`π':=algebraMap V0 Wn π`)への `V0`-代数準同型が
存在する。 -/
noncomputable def falt1BaseChangeAlgHom
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
    integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))) := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hfg : f.map (algebraMap V0 Wn) = g := by rw [hfdef, hgdef, hπ'def]; simp
  have e1 := falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  have e2 := falt1AdjoinRootEquivIntegralClosure π' n hn' hπne0' hprime' hnotsq' hnpos
  have step1 : AdjoinRoot f →ₐ[V0] AdjoinRoot g := hfg ▸ algHomAdjoinRootOfCompat' (Wn := Wn) f
  exact (e2.toAlgHom.restrictScalars V0).comp (step1.comp e1.symm.toAlgHom)

/-!
## `V_1` の生成元パッケージと `ψ(w)=x` の生成元対応(2026-09-04)

`falt1AdjoinRootEquivIntegralClosure`(`V1 ≃ AdjoinRoot f`)を手に入れた
おかげで、「`w`(`V1` の生成元)の主要な性質」を、以前のように
`PowerBasis`・`Eisenstein` の詳細を経由せず、**同型の自然性だけ**から
再導出できると判明した(`w := e1 (root f)` と定義するのが鍵——opaque
な `noncomputable def` を後から `unfold` して中の値を計算しようとする
より、最初から「後で使う値」を定義に組み込む方が遥かに簡単だった、
tools/lean-idioms.md #29 参照)。 -/

/-- **`V_1` の生成元パッケージ**: `w := e1(root f)`(`e1` は
`falt1AdjoinRootEquivIntegralClosure`)に対し、整性・単項生成性・
minpoly が `f` に一致することを、同型の自然性だけから示す。 -/
theorem falt1GeneratorPackage
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))] :
    ∃ w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))),
      IsIntegral V0 w ∧ Algebra.adjoin V0 ({w} : Set _) = ⊤ ∧ minpoly V0 w = Polynomial.X ^ n - Polynomial.C π := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set e1 := falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  refine ⟨e1 (AdjoinRoot.root f), ?_, ?_, ?_⟩
  · have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
    have hroot_int : IsIntegral V0 (AdjoinRoot.root f) := ⟨_, hmonic, AdjoinRoot.eval₂_root _⟩
    exact hroot_int.map e1.toAlgHom
  · have h1 : Algebra.adjoin V0 ({AdjoinRoot.root f} : Set (AdjoinRoot f)) = ⊤ :=
      AdjoinRoot.adjoinRoot_eq_top
    have h2 := congrArg (Subalgebra.map e1.toAlgHom) h1
    rw [AlgHom.map_adjoin_singleton, Algebra.map_top] at h2
    rwa [show e1.toAlgHom.range = ⊤ from by rw [AlgHom.range_eq_top]; exact e1.surjective] at h2
  · have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
    have hroot_int : IsIntegral V0 (AdjoinRoot.root f) := ⟨_, hmonic, AdjoinRoot.eval₂_root _⟩
    have hirr : Irreducible f := eisenstein_X_pow_sub_C π n hnpos hprime hnotsq
    have haeval_w : Polynomial.aeval (e1 (AdjoinRoot.root f)) f = 0 := by
      rw [Polynomial.aeval_algHom_apply]
      have h0 : Polynomial.aeval (AdjoinRoot.root f) f = 0 := by
        show Polynomial.eval₂ (algebraMap V0 (AdjoinRoot f)) (AdjoinRoot.root f) f = 0
        rw [AdjoinRoot.algebraMap_eq]; exact AdjoinRoot.eval₂_root f
      rw [h0, map_zero]
    have hw_int : IsIntegral V0 (e1 (AdjoinRoot.root f)) := hroot_int.map e1.toAlgHom
    have hdvd : minpoly V0 (e1 (AdjoinRoot.root f)) ∣ f := minpoly.isIntegrallyClosed_dvd hw_int haeval_w
    have hirr2 : Irreducible (minpoly V0 (e1 (AdjoinRoot.root f))) := minpoly.irreducible hw_int
    exact Polynomial.eq_of_monic_of_associated (minpoly.monic hw_int) hmonic
      (hirr2.associated_of_dvd hirr hdvd)

/-- `algHomAdjoinRootOfCompat'` の codomain を(多項式の等式に沿って)
cast したものの、根での値。`g` を**自由変数のまま**保つのが鍵——
`subst hfg` が素直に効く(`set`-束縛された `g` に対して `subst` すると
「invalid equality proof」で失敗する、tools/lean-idioms.md #29 の
続報)。 -/
theorem algHomAdjoinRootOfCompat'_cast_root {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn]
    (f : Polynomial V0) (g : Polynomial Wn) (hfg : f.map (algebraMap V0 Wn) = g) :
    (hfg ▸ algHomAdjoinRootOfCompat' (Wn := Wn) f : AdjoinRoot f →ₐ[V0] AdjoinRoot g)
      (AdjoinRoot.root f) = AdjoinRoot.root g := by
  subst hfg
  exact algHomAdjoinRootOfCompat'_root f

/-- **`ψ` の下での生成元対応**: `V1 → Wn1` の橋渡し(`falt1BaseChangeAlgHom`
と同じ構成、`w`(`V1` の `root f` 由来の生成元)を `x`(`Wn1` の
`root g` 由来の生成元)に送る)。`falt1BaseChangeAlgHom` という名前の
opaque な `def` を後から `unfold` して調べるのは時間対効果が悪いと
判明した(tools/lean-idioms.md #29)ため、同じ構成をここで直接書き
下して(名前を経由せず)対応を証明する。 -/
theorem falt1BaseChangeAlgHom_generator_correspondence
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))
      (w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
      (x : integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      ψ w = x := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hfg : f.map (algebraMap V0 Wn) = g := by rw [hfdef, hgdef, hπ'def]; simp
  set e1 := falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  set e2 := falt1AdjoinRootEquivIntegralClosure π' n hn' hπne0' hprime' hnotsq' hnpos
  set step1 : AdjoinRoot f →ₐ[V0] AdjoinRoot g := hfg ▸ algHomAdjoinRootOfCompat' (Wn := Wn) f
    with hstep1def
  have hstep1root : step1 (AdjoinRoot.root f) = AdjoinRoot.root g := by
    rw [hstep1def]; exact algHomAdjoinRootOfCompat'_cast_root f g hfg
  refine ⟨(e2.toAlgHom.restrictScalars V0).comp (step1.comp e1.symm.toAlgHom),
    e1 (AdjoinRoot.root f), e2 (AdjoinRoot.root g), ?_⟩
  rw [AlgHom.comp_apply, AlgHom.comp_apply]
  have h1 : e1.symm.toAlgHom (e1 (AdjoinRoot.root f)) = AdjoinRoot.root f := e1.symm_apply_apply _
  rw [h1, hstep1root]
  rfl

/-!
## `differentIdeal_tower_diamond` の instance 整備、第一弾(2026-09-04)

`differentIdeal_tower_diamond` を `Vn:=V0, Vn1:=V1, Wn:=Wn, Wn1:=Wn₁`
に適用するには約20個の instance を揃える必要がある——ここではまず
`Module.Finite Wn Wn1` を片付ける。`Algebra V0 Wn1`(合成: `V0→Wn→
FractionRing Wn→AdjoinRoot gK`)と `IsScalarTower V0 Wn Wn1` は
**instance 探索だけで自動的に見つかる**ことも確認した(実測、
2026-09-04)——明示的に構成する必要は無い。 -/

/-- **`Wₙ₊₁` は `Wₙ` 上有限**: `falt1AdjoinRootEquivIntegralClosure`
(`Wₙ₊₁ ≃ AdjoinRoot g`、`g`= base change 前の `X^n-π'`)と
`AdjoinRoot.powerBasis'`(**体でなく一般の可換環 `R` 上でも**、`g`
monic なら `PowerBasis R (AdjoinRoot g)` が取れる、mathlib の
既製品)を貼り合わせるだけ——`Polynomial.Monic.finite_adjoinRoot`
(体上限定)を経由する必要はなかった。 -/
theorem falt1ModuleFiniteWnWn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ) (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    Module.Finite Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))) := by
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hmonic : g.Monic := Polynomial.monic_X_pow_sub_C π' hnpos.ne'
  have e2 := falt1AdjoinRootEquivIntegralClosure π' n hn' hπne0' hprime' hnotsq' hnpos
  have pb := AdjoinRoot.powerBasis' hmonic
  haveI : Module.Finite Wn (AdjoinRoot g) := Module.Finite.of_basis pb.basis
  exact Module.Finite.equiv e2.toLinearEquiv

/-- **`Wₙ₊₁` は `Wₙ` 上捩れ無し**: `Module.IsTorsionFree V0
(integralClosure V0 (AdjoinRoot fK))`(`differentIdeal_eq_span_of_
adjoinRoot_X_pow_sub_C` の中で使った議論)と全く同じパターンを
`Wₙ`・`gK` に対して繰り返すだけ——`algebraMap Wₙ (AdjoinRoot gK)`
の単射性(`IsFractionRing.injective` 経由)から
`Module.IsTorsionFree Wₙ (AdjoinRoot gK)` を出し、
`IsIntegralClosure.isTorsionFree` で `Wₙ₊₁` へ降ろす。 -/
theorem falt1ModuleIsTorsionFreeWnWn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ) (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2) :
    Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))) := by
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  set gK : Polynomial (FractionRing Wn) := g.map (algebraMap Wn (FractionRing Wn)) with hgKdef
  have hirr : Irreducible gK :=
    eisenstein_X_pow_sub_C_irreducible_map π' n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible gK) := ⟨hirr⟩
  have hL : Module.IsTorsionFree Wn (AdjoinRoot gK) := by
    have hinj : Function.Injective (algebraMap Wn (AdjoinRoot gK)) := by
      rw [IsScalarTower.algebraMap_eq Wn (FractionRing Wn) (AdjoinRoot gK)]
      exact (algebraMap (FractionRing Wn) (AdjoinRoot gK)).injective.comp
        (IsFractionRing.injective Wn (FractionRing Wn))
    refine ⟨fun r hr => ?_⟩
    have hr0 : r ≠ 0 := by
      rintro rfl
      exact zero_ne_one (hr.left (show (0:Wn) * 0 = 0 * 1 by simp))
    have hr0' : algebraMap Wn (AdjoinRoot gK) r ≠ 0 := fun h => hr0 (hinj (by rw [h, map_zero]))
    intro x y hxy
    simp only [Algebra.smul_def] at hxy
    exact mul_left_cancel₀ hr0' hxy
  haveI := hL
  exact IsIntegralClosure.isTorsionFree Wn (AdjoinRoot gK)

/-- **`V_1` は `V_0` 上有限**: `falt1ModuleFiniteWnWn1` と全く同じ議論
(`AdjoinRoot.powerBasis'` + `falt1AdjoinRootEquivIntegralClosure`)を
`V_0`・`f` に対して繰り返すだけ。 -/
theorem falt1ModuleFiniteV0V1
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))] :
    Module.Finite V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0))))) := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  have hmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
  have e1 := falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  have pb := AdjoinRoot.powerBasis' hmonic
  haveI : Module.Finite V0 (AdjoinRoot f) := Module.Finite.of_basis pb.basis
  exact Module.Finite.equiv e1.toLinearEquiv

/-- **`V_1` は `V_0` 上捩れ無しであることは(instance 前提でなく)
実際に導出できる**——`falt1ModuleIsTorsionFreeWnWn1` と全く同じ議論
(`algebraMap V0 (AdjoinRoot fK)` の単射性 → `IsIntegralClosure.
isTorsionFree`)を `V_0`・`fK` に対して繰り返すだけ。★これまでの
多くの定理が `[Module.IsTorsionFree V0 (...)]` を**明示的な instance
前提**として要求していたが、実は毎回この議論で導出可能だったと判明
した(既存の証明済み定理群は変更しない——影響範囲を数える余裕が
無いため今回は追加のみに留める)。 -/
theorem falt1ModuleIsTorsionFreeV0V1
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n) :
    Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0))))) := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set fK : Polynomial (FractionRing V0) := f.map (algebraMap V0 (FractionRing V0)) with hfKdef
  have hirr : Irreducible fK :=
    eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirr⟩
  have hL : Module.IsTorsionFree V0 (AdjoinRoot fK) := by
    have hinj : Function.Injective (algebraMap V0 (AdjoinRoot fK)) := by
      rw [IsScalarTower.algebraMap_eq V0 (FractionRing V0) (AdjoinRoot fK)]
      exact (algebraMap (FractionRing V0) (AdjoinRoot fK)).injective.comp
        (IsFractionRing.injective V0 (FractionRing V0))
    refine ⟨fun r hr => ?_⟩
    have hr0 : r ≠ 0 := by
      rintro rfl
      exact zero_ne_one (hr.left (show (0:V0) * 0 = 0 * 1 by simp))
    have hr0' : algebraMap V0 (AdjoinRoot fK) r ≠ 0 := fun h => hr0 (hinj (by rw [h, map_zero]))
    intro x y hxy
    simp only [Algebra.smul_def] at hxy
    exact mul_left_cancel₀ hr0' hxy
  haveI := hL
  exact IsIntegralClosure.isTorsionFree V0 (AdjoinRoot fK)

/-- **`Module.IsTorsionFree` は `T` が整域で `algebraMap R T` が単射なら
自動的に成り立つ**——`falt1ModuleIsTorsionFreeWnWn1`・`falt1ModuleIsTorsionFreeV0V1`
で3回繰り返した議論(正則元 `r` が非零 → `algebraMap` の単射性で像も
非零 → 整域での非零元倍は単射)を、独立の汎用補題として切り出した。 -/
theorem moduleIsTorsionFree_of_injective {R T : Type*} [CommRing R] [Nontrivial R] [CommRing T] [IsDomain T]
    [Algebra R T] (hinj : Function.Injective (algebraMap R T)) : Module.IsTorsionFree R T := by
  refine ⟨fun r hr => ?_⟩
  have hr0 : r ≠ 0 := by
    rintro rfl
    exact zero_ne_one (hr.left (show (0:R) * 0 = 0 * 1 by simp))
  have hr0' : algebraMap R T r ≠ 0 := fun h => hr0 (hinj (by rw [h, map_zero]))
  intro x y hxy
  simp only [Algebra.smul_def] at hxy
  exact mul_left_cancel₀ hr0' hxy

/-- **`Wₙ₊₁` は `V0` 上有限**(`Module.Finite V0 Wₙ` + `Module.Finite Wₙ
Wₙ₊₁` の合成、`Module.Finite.trans`)。 -/
theorem falt1ModuleFiniteV0Wn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ) (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.Finite V0 Wn] :
    Module.Finite V0 (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))) := by
  haveI := falt1ModuleFiniteWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  exact Module.Finite.trans Wn _

/-- **`Wₙ₊₁` は `V0` 上捩れ無し**——`algebraMap V0 Wₙ₊₁` の単射性を
`algebraMap V0 Wₙ`(仮定)・`algebraMap Wₙ Wₙ₊₁`(`Wₙ₊₁` が
`AdjoinRoot gK` の整閉包という構成から)それぞれの単射性の合成で
示し、`moduleIsTorsionFree_of_injective` を適用する。 -/
theorem falt1ModuleIsTorsionFreeV0Wn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n) (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    Module.IsTorsionFree V0 (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))) := by
  apply moduleIsTorsionFree_of_injective
  rw [IsScalarTower.algebraMap_eq V0 Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
  simp only [RingHom.coe_comp]
  apply Function.Injective.comp _ hinjV0Wn
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  set gK : Polynomial (FractionRing Wn) := g.map (algebraMap Wn (FractionRing Wn)) with hgKdef
  have hirr : Irreducible gK :=
    eisenstein_X_pow_sub_C_irreducible_map π' n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible gK) := ⟨hirr⟩
  show Function.Injective (algebraMap Wn (integralClosure Wn (AdjoinRoot gK)))
  intro a b hab
  have heq : (algebraMap Wn (AdjoinRoot gK)) a = (algebraMap Wn (AdjoinRoot gK)) b := by
    have := congrArg (Subtype.val) hab
    simpa using this
  have hinj2 : Function.Injective (algebraMap Wn (AdjoinRoot gK)) := by
    rw [IsScalarTower.algebraMap_eq Wn (FractionRing Wn) (AdjoinRoot gK)]
    exact (algebraMap (FractionRing Wn) (AdjoinRoot gK)).injective.comp
      (IsFractionRing.injective Wn (FractionRing Wn))
  exact hinj2 heq

/-!
## `V_1 → W_1` の橋渡し写像の単射性(2026-09-04、新しい数学的内容)

`Module.Finite/IsTorsionFree V1 Wₙ₊₁` に要る最後のピース——`ψ:=
falt1BaseChangeAlgHom` が**単射**であることを、初めて実際に証明した。
鍵は `Polynomial.map_dvd_map`(mathlib、`x` が monic かつ係数の写像が
単射なら `x.map f ∣ y.map f ↔ x ∣ y`)——多項式の割り算アルゴリズムが
monic 除数に対して任意の可換環上で機能することの帰結。これで
`AdjoinRoot.map φ p (p.map φ) _` の**核**が自明であることが、
「`AdjoinRoot p` の元 `mk p r` が 0 に写る ⟺ `p.map φ ∣ r.map φ` ⟺
(単射性+monic) `p ∣ r` ⟺ 元から 0」という計算だけで閉じる——
`PowerBasis` の基底展開等は一切不要だった。 -/

/-- `AdjoinRoot.map` の `mk` への作用——`AdjoinRoot.lift`(`map` の定義)
+ `Polynomial.eval₂_map` + `AdjoinRoot.aeval_eq` から。 -/
theorem adjoinRoot_map_mk {R S : Type*} [CommRing R] [CommRing S] (φ : R →+* S) (p : Polynomial R)
    (q : Polynomial S) (h : q ∣ p.map φ) (r : Polynomial R) :
    AdjoinRoot.map φ p q h (AdjoinRoot.mk p r) = AdjoinRoot.mk q (r.map φ) := by
  unfold AdjoinRoot.map
  rw [AdjoinRoot.lift_mk, ← Polynomial.eval₂_map]
  rw [← AdjoinRoot.algebraMap_eq, ← Polynomial.aeval_def]
  exact AdjoinRoot.aeval_eq _

/-- **`AdjoinRoot.map` は `q = p.map φ`(base change そのもの)の場合、
`p` が monic・`φ` が単射なら単射**——`Polynomial.map_dvd_map` を
`AdjoinRoot.mk_eq_zero` に橋渡しするだけ。 -/
theorem adjoinRoot_map_injective_of_map_eq {R S : Type*} [CommRing R] [CommRing S] (φ : R →+* S)
    (hφinj : Function.Injective φ) (p : Polynomial R) (hpmonic : p.Monic) :
    Function.Injective (AdjoinRoot.map φ p (p.map φ) (dvd_refl _)) := by
  rw [injective_iff_map_eq_zero]
  intro a ha
  obtain ⟨r, rfl⟩ := AdjoinRoot.mk_surjective a
  rw [adjoinRoot_map_mk] at ha
  rw [AdjoinRoot.mk_eq_zero] at ha ⊢
  exact (Polynomial.map_dvd_map φ hφinj hpmonic).mp ha

/-- `algHomAdjoinRootOfCompat'` の単射性(`AdjoinRoot.map` の単射性を
そのまま流用——`AlgHom.mk'` は同じ関数を包むだけなので変わらない)。 -/
theorem algHomAdjoinRootOfCompat'_injective {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn]
    (f : Polynomial V0) (hfmonic : f.Monic) (hinj : Function.Injective (algebraMap V0 Wn)) :
    Function.Injective (algHomAdjoinRootOfCompat' (Wn := Wn) f) :=
  adjoinRoot_map_injective_of_map_eq (algebraMap V0 Wn) hinj f hfmonic

/-- `algHomAdjoinRootOfCompat'_injective` の(多項式の等式に沿って)
cast したもの版——`g` を自由変数のまま保つのが鍵(tools/lean-idioms.md
#29 と同じ理由)。 -/
theorem algHomAdjoinRootOfCompat'_cast_injective {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn]
    (f : Polynomial V0) (g : Polynomial Wn) (hfg : f.map (algebraMap V0 Wn) = g)
    (hfmonic : f.Monic) (hinj : Function.Injective (algebraMap V0 Wn)) :
    Function.Injective (hfg ▸ algHomAdjoinRootOfCompat' (Wn := Wn) f : AdjoinRoot f →ₐ[V0] AdjoinRoot g) := by
  subst hfg
  exact algHomAdjoinRootOfCompat'_injective f hfmonic hinj

/-- **`V_1 → W_1` の橋渡しが生成元を対応させ、かつ単射である**
(`falt1BaseChangeAlgHom_generator_correspondence` の拡張版)。単射性は
`algHomAdjoinRootOfCompat'_cast_injective`(base change の芯)と
`e1.symm`・`e2` の全単射性の合成——これで `Module.Finite/IsTorsionFree
V1 Wₙ₊₁` を `moduleIsTorsionFree_of_injective`(単射性から直接)と
「単射像は有限生成」型の議論で得る道が開けた。 -/
theorem falt1BaseChangeAlgHom_generator_and_injective
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))
      (w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
      (x : integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      ψ w = x ∧ Function.Injective ψ := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hfg : f.map (algebraMap V0 Wn) = g := by rw [hfdef, hgdef, hπ'def]; simp
  set e1 := falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  set e2 := falt1AdjoinRootEquivIntegralClosure π' n hn' hπne0' hprime' hnotsq' hnpos
  set step1 : AdjoinRoot f →ₐ[V0] AdjoinRoot g := hfg ▸ algHomAdjoinRootOfCompat' (Wn := Wn) f
    with hstep1def
  have hstep1root : step1 (AdjoinRoot.root f) = AdjoinRoot.root g := by
    rw [hstep1def]; exact algHomAdjoinRootOfCompat'_cast_root f g hfg
  have hfmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
  have hstep1inj : Function.Injective step1 := by
    rw [hstep1def]
    exact algHomAdjoinRootOfCompat'_cast_injective f g hfg hfmonic hinjV0Wn
  refine ⟨(e2.toAlgHom.restrictScalars V0).comp (step1.comp e1.symm.toAlgHom),
    e1 (AdjoinRoot.root f), e2 (AdjoinRoot.root g), ?_, ?_⟩
  · rw [AlgHom.comp_apply, AlgHom.comp_apply]
    have h1 : e1.symm.toAlgHom (e1 (AdjoinRoot.root f)) = AdjoinRoot.root f := e1.symm_apply_apply _
    rw [h1, hstep1root]
    rfl
  · intro a b hab
    simp only [AlgHom.comp_apply] at hab
    apply e1.symm.injective
    apply hstep1inj
    exact e2.injective hab

/-- **`Module.IsTorsionFree V1 Wₙ₊₁` が完成した**——`falt1BaseChangeAlgHom_
generator_and_injective` の単射性を `moduleIsTorsionFree_of_injective`
に直接渡すだけ。`ψ` を `letI := ψ.toRingHom.toAlgebra` で `Algebra V1
Wₙ₊₁` として登録した**その instance の下で**主張する(実際に
`differentIdeal_tower_diamond` を呼ぶ際もこの `letI` を経由することに
なる)。 -/
theorem falt1ModuleIsTorsionFreeV1Wn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      letI := ψ.toRingHom.toAlgebra
      Module.IsTorsionFree
        (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
          (algebraMap V0 (FractionRing V0)))))
        (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
          (algebraMap Wn (FractionRing Wn)))))
    := by
  obtain ⟨ψ, w, x, hwx, hψinj⟩ :=
    falt1BaseChangeAlgHom_generator_and_injective π n hn hπne0 hprime hnotsq hnpos
      hn' hπne0' hprime' hnotsq' hinjV0Wn
  refine ⟨ψ, ?_⟩
  letI := ψ.toRingHom.toAlgebra
  apply moduleIsTorsionFree_of_injective
  exact hψinj

set_option maxHeartbeats 400000 in
set_option synthInstance.maxHeartbeats 400000 in
/-- **`Module.Finite V1 Wₙ₊₁` が完成した**——`Module.Finite.of_
restrictScalars_finite`(`R→A→M` のタワーで `M` が `R` 上有限なら
`A` 上でも有限、「小さい環上の有限生成は大きい環上でも有限生成」)を
`R:=V0, A:=V1, M:=Wₙ₊₁` に適用するだけ。`IsScalarTower V0 V1 Wₙ₊₁`
(`ψ.commutes` から)と `Module.Finite V0 Wₙ₊₁`(既存の
`falt1ModuleFiniteV0Wn1`)の2つを揃えれば機械的に閉じる——`ψ` の
単射性は不要だった(有限性には使わない、`IsTorsionFree` の方でのみ
必要だった)。 -/
theorem falt1ModuleFiniteV1Wn1
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.Finite V0 Wn] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      letI := ψ.toRingHom.toAlgebra
      Module.Finite
        (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
          (algebraMap V0 (FractionRing V0)))))
        (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
          (algebraMap Wn (FractionRing Wn)))))
    := by
  obtain ⟨ψ, w, x, hwx, hψinj⟩ :=
    falt1BaseChangeAlgHom_generator_and_injective π n hn hπne0 hprime hnotsq hnpos
      hn' hπne0' hprime' hnotsq' hinjV0Wn
  refine ⟨ψ, ?_⟩
  letI := ψ.toRingHom.toAlgebra
  haveI : IsScalarTower V0
      (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
      (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))) := by
    apply IsScalarTower.of_algebraMap_eq
    intro y
    exact (ψ.commutes y).symm
  haveI := falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
  exact Module.Finite.of_restrictScalars_finite V0 _ _

set_option maxHeartbeats 1000000 in
set_option synthInstance.maxHeartbeats 1000000 in
/-- `differentIdeal_tower_diamond` の `hsep`(`Algebra.IsSeparable (FractionRing V0)
(FractionRing Wn1)`)を解決する。`FractionRing.liftAlgebra V0 (FractionRing Wn1)` で
作った `Algebra` instance と、その instance に対する separability の証明を**一緒に**返す
(`Σ'` で束ねる——`theorem` ではなく `def` にしたのは、戻り値の型が `Prop` ではなく
`Type` だから)。呼び出し側は `.1` を `letI` で登録すれば `.2` がそのまま使える。

証明の道筋: `Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)`(仮定)と
`Algebra.IsSeparable (FractionRing Wn) (AdjoinRoot gK)`(`gK` の separability から)を
`Algebra.IsSeparable.trans` で繋ぐ。ただし繋ぐ先は `FractionRing Wn1`(`AdjoinRoot gK`
ではない)なので、`FractionRing Wn1 ≃ₐ[Wn1] AdjoinRoot gK`
(`IsLocalization.algEquiv`)を V0-level・Wn-level それぞれで「両立する
`FractionRing.liftAlgebra` instance」に対して `AlgEquiv.ofRingEquiv` で持ち上げ、
`IsLocalization.ringHom_ext` で一意性から可換性を出す、という二段階の「diamond」を
どちらも `rfl`(部分体への埋め込みが `Subalgebra` の値射影と一致)+
`IsScalarTower.algebraMap_apply` の pointwise な書き換えで解消する。 -/
noncomputable def falt1_hsep_bundled
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ) (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [Algebra (FractionRing V0) (FractionRing Wn)]
    [IsScalarTower V0 (FractionRing V0) (FractionRing Wn)]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)] :
    Σ' (inst : Algebra (FractionRing V0) (FractionRing (integralClosure Wn (AdjoinRoot
          (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
            (algebraMap Wn (FractionRing Wn))))))),
      @Algebra.IsSeparable (FractionRing V0) (FractionRing (integralClosure Wn (AdjoinRoot
          (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
            (algebraMap Wn (FractionRing Wn)))))) _ _ inst := by
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  set gK : Polynomial (FractionRing Wn) := g.map (algebraMap Wn (FractionRing Wn)) with hgKdef
  have hirr : Irreducible gK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π' n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible gK) := ⟨hirr⟩
  have hmonicK : gK.Monic := (Polynomial.monic_X_pow_sub_C π' hnpos.ne').map _
  have hsepK : gK.Separable := by
    have hmapeq : gK = Polynomial.X ^ n - Polynomial.C (algebraMap Wn (FractionRing Wn) π') := by
      rw [hgKdef, hgdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn' hπne0'
  haveI : Module.Finite (FractionRing Wn) (AdjoinRoot gK) := hmonicK.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing Wn) (AdjoinRoot gK) := ‹Module.Finite _ _›
  haveI hsepAdjoin : Algebra.IsSeparable (FractionRing Wn) (AdjoinRoot gK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicK hsepK
  haveI : IsDedekindDomain (integralClosure Wn (AdjoinRoot gK)) := by infer_instance
  haveI hFR : IsFractionRing (integralClosure Wn (AdjoinRoot gK)) (AdjoinRoot gK) :=
    integralClosure.isFractionRing_of_finite_extension (A := Wn) (FractionRing Wn) (AdjoinRoot gK)
  set Wn1 := integralClosure Wn (AdjoinRoot gK) with hWn1def
  have hinjV0Wn1 : Function.Injective (algebraMap V0 Wn1) := by
    rw [IsScalarTower.algebraMap_eq V0 Wn Wn1]
    simp only [RingHom.coe_comp]
    apply Function.Injective.comp _ hinjV0Wn
    show Function.Injective (algebraMap Wn Wn1)
    intro a b hab
    have heq : (algebraMap Wn (AdjoinRoot gK)) a = (algebraMap Wn (AdjoinRoot gK)) b := by
      have := congrArg (Subtype.val) hab
      simpa using this
    have hinj2 : Function.Injective (algebraMap Wn (AdjoinRoot gK)) := by
      rw [IsScalarTower.algebraMap_eq Wn (FractionRing Wn) (AdjoinRoot gK)]
      exact (algebraMap (FractionRing Wn) (AdjoinRoot gK)).injective.comp
        (IsFractionRing.injective Wn (FractionRing Wn))
    exact hinj2 heq
  haveI hFaithV0 : FaithfulSMul V0 (FractionRing Wn1) := by
    rw [faithfulSMul_iff_algebraMap_injective]
    rw [IsScalarTower.algebraMap_eq V0 Wn1 (FractionRing Wn1)]
    exact (IsFractionRing.injective Wn1 (FractionRing Wn1)).comp hinjV0Wn1
  letI algFRWn1 : Algebra (FractionRing V0) (FractionRing Wn1) :=
    FractionRing.liftAlgebra V0 (FractionRing Wn1)
  haveI hSTcanon : IsScalarTower V0 (FractionRing V0) (FractionRing Wn1) := by infer_instance
  have e_Wn1 : FractionRing Wn1 ≃ₐ[Wn1] AdjoinRoot gK :=
    IsLocalization.algEquiv (nonZeroDivisors Wn1) (FractionRing Wn1) (AdjoinRoot gK)
  have h2' : (algebraMap Wn1 (AdjoinRoot gK)).comp (algebraMap Wn Wn1) = algebraMap Wn (AdjoinRoot gK) := by
    ext w; show ((algebraMap Wn Wn1) w : AdjoinRoot gK) = algebraMap Wn (AdjoinRoot gK) w; rfl
  haveI hFaithWn : FaithfulSMul Wn (FractionRing Wn1) := by
    rw [faithfulSMul_iff_algebraMap_injective]
    rw [IsScalarTower.algebraMap_eq Wn Wn1 (FractionRing Wn1)]
    refine (IsFractionRing.injective Wn1 (FractionRing Wn1)).comp ?_
    intro a b hab
    have heq : (algebraMap Wn (AdjoinRoot gK)) a = (algebraMap Wn (AdjoinRoot gK)) b := by
      have := congrArg (Subtype.val) hab; simpa using this
    have hinj2 : Function.Injective (algebraMap Wn (AdjoinRoot gK)) := by
      rw [IsScalarTower.algebraMap_eq Wn (FractionRing Wn) (AdjoinRoot gK)]
      exact (algebraMap (FractionRing Wn) (AdjoinRoot gK)).injective.comp
        (IsFractionRing.injective Wn (FractionRing Wn))
    exact hinj2 heq
  letI algFRWn_Wn1 : Algebra (FractionRing Wn) (FractionRing Wn1) :=
    FractionRing.liftAlgebra Wn (FractionRing Wn1)
  haveI hSTWn : IsScalarTower Wn (FractionRing Wn) (FractionRing Wn1) := by infer_instance
  have hringhom_eq' : ∀ w : Wn, (e_Wn1.toRingHom) (algebraMap (FractionRing Wn) (FractionRing Wn1)
        (algebraMap Wn (FractionRing Wn) w)) = algebraMap (FractionRing Wn) (AdjoinRoot gK)
        (algebraMap Wn (FractionRing Wn) w) := by
    intro w
    rw [← IsScalarTower.algebraMap_apply Wn (FractionRing Wn) (FractionRing Wn1)]
    rw [← IsScalarTower.algebraMap_apply Wn (FractionRing Wn) (AdjoinRoot gK)]
    rw [IsScalarTower.algebraMap_apply Wn Wn1 (FractionRing Wn1)]
    show e_Wn1 (algebraMap Wn1 (FractionRing Wn1) (algebraMap Wn Wn1 w)) = algebraMap Wn (AdjoinRoot gK) w
    rw [e_Wn1.commutes]
    exact congrFun (congrArg DFunLike.coe h2') w
  have hringhom_eq'' : (e_Wn1.toRingHom).comp (algebraMap (FractionRing Wn) (FractionRing Wn1))
      = algebraMap (FractionRing Wn) (AdjoinRoot gK) :=
    IsLocalization.ringHom_ext (nonZeroDivisors Wn) (RingHom.ext hringhom_eq')
  have e_Wn1'' : FractionRing Wn1 ≃ₐ[FractionRing Wn] AdjoinRoot gK :=
    AlgEquiv.ofRingEquiv (f := e_Wn1.toRingEquiv) (fun x => by
      have := congrFun (congrArg DFunLike.coe hringhom_eq'') x; simpa using this)
  haveI hsepFRWn1 : Algebra.IsSeparable (FractionRing Wn) (FractionRing Wn1) :=
    (AlgEquiv.Algebra.isSeparable_iff e_Wn1'').mpr hsepAdjoin
  haveI hSTV0Wn1 : IsScalarTower V0 Wn (FractionRing Wn1) := by infer_instance
  have hSTfull_pt : ∀ v : V0,
      algebraMap (FractionRing V0) (FractionRing Wn1) (algebraMap V0 (FractionRing V0) v)
        = algebraMap (FractionRing Wn) (FractionRing Wn1)
            (algebraMap (FractionRing V0) (FractionRing Wn) (algebraMap V0 (FractionRing V0) v)) := by
    intro v
    rw [← IsScalarTower.algebraMap_apply V0 (FractionRing V0) (FractionRing Wn1)]
    rw [← IsScalarTower.algebraMap_apply V0 (FractionRing V0) (FractionRing Wn)]
    rw [IsScalarTower.algebraMap_apply V0 Wn (FractionRing Wn)]
    rw [← IsScalarTower.algebraMap_apply Wn (FractionRing Wn) (FractionRing Wn1)]
    exact IsScalarTower.algebraMap_apply V0 Wn (FractionRing Wn1) v
  haveI hSTfull : IsScalarTower (FractionRing V0) (FractionRing Wn) (FractionRing Wn1) := by
    apply IsScalarTower.of_algebraMap_eq'
    apply IsLocalization.ringHom_ext (nonZeroDivisors V0)
    exact RingHom.ext hSTfull_pt
  haveI hsepFinal : Algebra.IsSeparable (FractionRing V0) (FractionRing Wn1) :=
    Algebra.IsSeparable.trans (FractionRing V0) (FractionRing Wn) (FractionRing Wn1)
  exact ⟨algFRWn1, hsepFinal⟩

set_option maxHeartbeats 2000000 in
set_option synthInstance.maxHeartbeats 2000000 in
/-- **`differentIdeal_tower_diamond` を Falt1 の具体的な `V0→V1`・`Wn→Wn1` の構成
(`Vₙ₊₁ := integralClosure Vₙ (AdjoinRoot(X^n-Cπ の base change))`)に対して初めて
実際に呼び出した記録**(Theorem 1.2・item 3c の中核道具の組み立て完了)。

戻り値をあえて `True` にしている——理由は、`differentIdeal V1 Wn1` を(呼び出し側が
`ψ : V1 →ₐ[V0] Wn1` を選ぶ**前**に)シグネチャの型として書こうとすると、
`Module.IsTorsionFree V1 Wn1`(`ψ` の単射性からしか出ない事実)の instance 検索が
シグネチャ**そのものの型検査**(証明本体に入る前)で走ってしまい、168 秒
(`maxHeartbeats 2000000`)でも終わらない(`«synthesize pending MVars»` timeout)
ことを実測したため。したがって「呼び出し可能な補題」としての再パッケージングは
次回以降の課題として保留し、**この証明が実際に構築する項** `hdiamond` の型を
以下にそのまま記録する(この docstring の型と証明中の `hdiamond` の型は同一):

```
hdiamond : differentIdeal V1 Wn1 * Ideal.map (algebraMap V1 Wn1) (differentIdeal V0 V1)
         = differentIdeal Wn Wn1 * Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn)
```
(`V1 := integralClosure V0 (AdjoinRoot fK)`、`Wn1 := integralClosure Wn (AdjoinRoot gK)`、
`fK`・`gK` は `X^n - C π`(それぞれ `V0`・`Wn` 上)の `FractionRing` への base change、
`Algebra V1 Wn1` は `ψ.toRingHom.toAlgebra`(`ψ` は base change 写像から得られる
生成元対応))。

組み立て手順(20 個の instance を機械的に集める、`falt1-goal.md` に詳細記録):
`IsDedekindDomain`(V1・Wn1、`isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C`)
→ `Module.IsTorsionFree`/`Module.Finite`(V0-V1・Wn-Wn1・V0-Wn1 の3対、既存の
`falt1Module*`群)→ `ψ`(`falt1BaseChangeAlgHom_generator_and_injective`)→
`Algebra V1 Wn1 := ψ.toRingHom.toAlgebra` → `IsScalarTower V0 V1 Wn1`(`ψ.commutes`)
→ `Module.Finite/IsTorsionFree V1 Wn1`(`ψ`の単射性 + `Module.Finite.of_restrictScalars_finite`）
→ `hsep`(`falt1_hsep_bundled`)→ `differentIdeal_tower_diamond hsep`。 -/
theorem falt1_differentIdeal_tower_diamond_assembled
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [Module.Finite V0 Wn]
    [Algebra (FractionRing V0) (FractionRing Wn)]
    [IsScalarTower V0 (FractionRing V0) (FractionRing Wn)]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)] :
    True := by
  set fK : Polynomial (FractionRing V0) := (Polynomial.X ^ n - Polynomial.C π).map (algebraMap V0 (FractionRing V0)) with hfKdef
  set gK : Polynomial (FractionRing Wn) := (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π)).map (algebraMap Wn (FractionRing Wn)) with hgKdef
  haveI : IsDedekindDomain V0 := by infer_instance
  haveI : IsDedekindDomain Wn := by infer_instance
  haveI hDedV1 : IsDedekindDomain (integralClosure V0 (AdjoinRoot fK)) :=
    ABC3.Found.Falt1.isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C π n hn hπne0 hprime hnotsq hnpos
  haveI hDedWn1 : IsDedekindDomain (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.isDedekindDomain_integralClosure_adjoinRoot_X_pow_sub_C
      (algebraMap V0 Wn π) n hn' hπne0' hprime' hnotsq' hnpos
  haveI hTFV1 : Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot fK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeV0V1 π n hprime hnotsq hnpos
  haveI hTFWn1 : Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  haveI hTFV0Wn1 : Module.IsTorsionFree V0 (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeV0Wn1 π n hprime' hnotsq' hnpos hinjV0Wn
  haveI hFinWnWn1 : Module.Finite Wn (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleFiniteWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  set V1 := integralClosure V0 (AdjoinRoot fK) with hV1def
  set Wn1 := integralClosure Wn (AdjoinRoot gK) with hWn1def
  obtain ⟨ψ, w, x, hwx, hψinj⟩ :=
    ABC3.Found.Falt1.falt1BaseChangeAlgHom_generator_and_injective π n hn hπne0 hprime hnotsq hnpos
      hn' hπne0' hprime' hnotsq' hinjV0Wn
  letI algV1Wn1 : Algebra V1 Wn1 := ψ.toRingHom.toAlgebra
  haveI hSTV0V1Wn1 : IsScalarTower V0 V1 Wn1 := by
    apply IsScalarTower.of_algebraMap_eq
    intro y
    exact (ψ.commutes y).symm
  haveI hSTV0WnWn1 : IsScalarTower V0 Wn Wn1 := by infer_instance
  haveI hFinV0V1 : Module.Finite V0 V1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0V1 π n hn hπne0 hprime hnotsq hnpos
  haveI hFinV0Wn1 : Module.Finite V0 Wn1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
  haveI hFinV1Wn1 : Module.Finite V1 Wn1 := by
    haveI := ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
    exact Module.Finite.of_restrictScalars_finite V0 V1 Wn1
  haveI hIsDomV1 : IsDomain V1 := inferInstance
  haveI hIsDomWn1 : IsDomain Wn1 := inferInstance
  haveI hTFV1Wn1 : Module.IsTorsionFree V1 Wn1 :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hψinj
  haveI hTFV0Wn : Module.IsTorsionFree V0 Wn :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hinjV0Wn
  haveI hIsIntClosedV0 : IsIntegrallyClosed V0 := inferInstance
  letI algFRWn1 := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).1
  have hsep := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).2
  have hdiamond := ABC3.Found.Falt1.differentIdeal_tower_diamond
    (Vn := V0) (Vn1 := V1) (Wn := Wn) (Wn1 := Wn1) hsep
  trivial

/-!
## `hspan_eq`(`cancel_conductor_delta`)接続への足場(2026-09-04、3つの補題完成)

`cancel_conductor_delta` の `hspan_eq`(`conductor Wₙ x` 側の span と
`differentIdeal_tower_diamond` の base change 側の span が一致する)を
確立するための、3つの再利用可能な補題を用意した(いずれも `lake build`・
`#print axioms` 確認済み・sorry 無し)。残る接続(`differentIdeal_eq_
span_derivative`・`conductor_mul_differentIdeal` を実際に呼び出して
`hspan_eq` そのものを閉じる)は、`w`・`x` の field-level 生成元性を
`differentIdeal_eq_span_derivative`/`conductor_mul_differentIdeal` が
要求する**正確な instance 経路**(`algebraMap V1 (AdjoinRoot fK)` 等)に
一致させる作業で type mismatch に当たり、次回への持ち越しとした
(falt1-goal.md 参照)。 -/

/-- **有限次拡大で、剰余体レベルの生成が全体生成なら基礎体レベルでも
全体生成**(単純な次元カウント)。`IsIntegral K x` かつ
`(minpoly K x).natDegree = finrank K L` から `Algebra.adjoin K {x} = ⊤`
を出す——`Algebra.adjoin.powerBasis'` の次元 `= (minpoly K x).natDegree`
と、有限次元での「同次元の部分加群は全体」(`Submodule.eq_top_of_
finrank_eq`)を組み合わせるだけ。 -/
theorem falt1_adjoin_top_of_finrank_eq {K L : Type*} [Field K] [Field L] [Algebra K L]
    [FiniteDimensional K L] (x : L) (hint : IsIntegral K x)
    (hdeg : (minpoly K x).natDegree = Module.finrank K L) :
    Algebra.adjoin K ({x} : Set L) = ⊤ := by
  have h1 : Module.finrank K (Algebra.adjoin K ({x} : Set L)) = Module.finrank K L := by
    rw [(Algebra.adjoin.powerBasis' hint).finrank, Algebra.adjoin.powerBasis'_dim hint, hdeg]
  have h2 : Subalgebra.toSubmodule (Algebra.adjoin K ({x} : Set L)) = ⊤ := by
    apply Submodule.eq_top_of_finrank_eq
    rwa [Subalgebra.finrank_toSubmodule]
  have h3 : Subalgebra.toSubmodule (⊤ : Subalgebra K L) = ⊤ := by simp
  exact Subalgebra.toSubmodule_injective (h2.trans h3.symm)

set_option maxHeartbeats 1000000 in
/-- **`V_1`(または `Wₙ₊₁`)の元 `w0` の環レベル `minpoly` が
`X^n-Cπ` に一致するなら、体レベル(`FractionRing V0` 上)でも
`w0` が全体を生成する**——`minpoly.isIntegrallyClosed_eq_field_
fractions'`(整閉環上の元の体への minpoly は係数を base change した
もの)で体レベル minpoly を `fK` と同定し、`falt1_adjoin_top_of_
finrank_eq` に渡す。`hspan_eq` が要求する `K[(algebraMap W L) w] = ⊤`
の直接の供給源。 -/
theorem falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly {V0 : Type*} [CommRing V0] [IsDomain V0]
    [IsDiscreteValuationRing V0] (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    (w0 : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
    (hw0minpoly : minpoly V0 w0 = Polynomial.X ^ n - Polynomial.C π) :
    Algebra.adjoin (FractionRing V0) ({(algebraMap _
        (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))))) w0} :
      Set (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))))) = ⊤ := by
  have hirr : Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))) :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) := ⟨hirr⟩
  have hmonicK : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Monic :=
    (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  have hsepK : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Separable := by
    have hmapeq : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))
        = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) := hmonicK.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) := ‹Module.Finite _ _›
  have hfmonic : (Polynomial.X ^ n - Polynomial.C π : Polynomial V0).Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
  have hw0int : IsIntegral V0 w0 := by
    rw [← minpoly.ne_zero_iff, hw0minpoly]; exact hfmonic.ne_zero
  have hw0int' : IsIntegral V0 (algebraMap _ (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) w0) :=
    hw0int.map (IsScalarTower.toAlgHom V0 _ (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))))
  have hw0intK : IsIntegral (FractionRing V0) (algebraMap _ (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) w0) :=
    IsIntegral.of_finite (FractionRing V0) _
  have hminpolycast : minpoly V0 (algebraMap _ (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) w0) = minpoly V0 w0 :=
    minpoly.algebraMap_eq (Subtype.val_injective) w0
  have hminpolyL : minpoly (FractionRing V0) (algebraMap _ (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) w0) = (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))) := by
    rw [minpoly.isIntegrallyClosed_eq_field_fractions' (FractionRing V0) hw0int']
    rw [hminpolycast, hw0minpoly]
  apply falt1_adjoin_top_of_finrank_eq
  · exact hw0intK
  · rw [hminpolyL]
    rw [(AdjoinRoot.powerBasis hirr.ne_zero).finrank, AdjoinRoot.powerBasis_dim hirr.ne_zero]

set_option maxHeartbeats 1000000 in
/-- **`falt1BaseChangeAlgHom_generator_and_injective` の `w`・`x` の
生成元性(`falt1GeneratorPackage` 相当)を、`ψ w = x` と同時に手に入れる**
——`hspan_eq` の接続には `w`・`x` が「同じ `e1`/`e2` 由来」であることが
必須(別々に `falt1GeneratorPackage` を呼ぶと `obtain` した witness が
一致する保証がない)ため、両方の証明を1つの `e1`・`e2` から一括して
組み立てる。 -/
theorem falt1BaseChangeGeneratorFull
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))
      (w : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))
      (x : integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      ψ w = x ∧ Function.Injective ψ ∧
      Algebra.adjoin V0 ({w} : Set _) = ⊤ ∧ minpoly V0 w = Polynomial.X ^ n - Polynomial.C π ∧
      Algebra.adjoin Wn ({x} : Set _) = ⊤ ∧
        minpoly Wn x = Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set g : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hfg : f.map (algebraMap V0 Wn) = g := by rw [hfdef, hgdef, hπ'def]; simp
  set e1 := ABC3.Found.Falt1.falt1AdjoinRootEquivIntegralClosure π n hn hπne0 hprime hnotsq hnpos
  set e2 := ABC3.Found.Falt1.falt1AdjoinRootEquivIntegralClosure π' n hn' hπne0' hprime' hnotsq' hnpos
  set step1 : AdjoinRoot f →ₐ[V0] AdjoinRoot g :=
    hfg ▸ ABC3.Found.Falt1.algHomAdjoinRootOfCompat' (Wn := Wn) f with hstep1def
  have hstep1root : step1 (AdjoinRoot.root f) = AdjoinRoot.root g := by
    rw [hstep1def]; exact ABC3.Found.Falt1.algHomAdjoinRootOfCompat'_cast_root f g hfg
  have hfmonic : f.Monic := Polynomial.monic_X_pow_sub_C π hnpos.ne'
  have hstep1inj : Function.Injective step1 := by
    rw [hstep1def]
    exact ABC3.Found.Falt1.algHomAdjoinRootOfCompat'_cast_injective f g hfg hfmonic hinjV0Wn
  have hwadjoin : Algebra.adjoin V0 ({e1 (AdjoinRoot.root f)} : Set _) = ⊤ := by
    have h1 : Algebra.adjoin V0 ({AdjoinRoot.root f} : Set (AdjoinRoot f)) = ⊤ :=
      AdjoinRoot.adjoinRoot_eq_top
    have h2 := congrArg (Subalgebra.map e1.toAlgHom) h1
    rw [AlgHom.map_adjoin_singleton, Algebra.map_top] at h2
    rwa [show e1.toAlgHom.range = ⊤ from by rw [AlgHom.range_eq_top]; exact e1.surjective] at h2
  have hwminpoly : minpoly V0 (e1 (AdjoinRoot.root f)) = f := by
    have hmonic : f.Monic := hfmonic
    have hroot_int : IsIntegral V0 (AdjoinRoot.root f) := ⟨_, hmonic, AdjoinRoot.eval₂_root _⟩
    have hirr : Irreducible f := ABC3.Found.Falt1.eisenstein_X_pow_sub_C π n hnpos hprime hnotsq
    have haeval_w : Polynomial.aeval (e1 (AdjoinRoot.root f)) f = 0 := by
      rw [Polynomial.aeval_algHom_apply]
      have h0 : Polynomial.aeval (AdjoinRoot.root f) f = 0 := by
        show Polynomial.eval₂ (algebraMap V0 (AdjoinRoot f)) (AdjoinRoot.root f) f = 0
        rw [AdjoinRoot.algebraMap_eq]; exact AdjoinRoot.eval₂_root f
      rw [h0, map_zero]
    have hw_int : IsIntegral V0 (e1 (AdjoinRoot.root f)) := hroot_int.map e1.toAlgHom
    have hdvd : minpoly V0 (e1 (AdjoinRoot.root f)) ∣ f := minpoly.isIntegrallyClosed_dvd hw_int haeval_w
    have hirr2 : Irreducible (minpoly V0 (e1 (AdjoinRoot.root f))) := minpoly.irreducible hw_int
    exact Polynomial.eq_of_monic_of_associated (minpoly.monic hw_int) hmonic
      (hirr2.associated_of_dvd hirr hdvd)
  have hxadjoin : Algebra.adjoin Wn ({e2 (AdjoinRoot.root g)} : Set _) = ⊤ := by
    have h1 : Algebra.adjoin Wn ({AdjoinRoot.root g} : Set (AdjoinRoot g)) = ⊤ :=
      AdjoinRoot.adjoinRoot_eq_top
    have h2 := congrArg (Subalgebra.map e2.toAlgHom) h1
    rw [AlgHom.map_adjoin_singleton, Algebra.map_top] at h2
    rwa [show e2.toAlgHom.range = ⊤ from by rw [AlgHom.range_eq_top]; exact e2.surjective] at h2
  have hgmonic : g.Monic := Polynomial.monic_X_pow_sub_C π' hnpos.ne'
  have hxminpoly : minpoly Wn (e2 (AdjoinRoot.root g)) = g := by
    have hroot_int : IsIntegral Wn (AdjoinRoot.root g) := ⟨_, hgmonic, AdjoinRoot.eval₂_root _⟩
    have hirr : Irreducible g := ABC3.Found.Falt1.eisenstein_X_pow_sub_C π' n hnpos hprime' hnotsq'
    have haeval_x : Polynomial.aeval (e2 (AdjoinRoot.root g)) g = 0 := by
      rw [Polynomial.aeval_algHom_apply]
      have h0 : Polynomial.aeval (AdjoinRoot.root g) g = 0 := by
        show Polynomial.eval₂ (algebraMap Wn (AdjoinRoot g)) (AdjoinRoot.root g) g = 0
        rw [AdjoinRoot.algebraMap_eq]; exact AdjoinRoot.eval₂_root g
      rw [h0, map_zero]
    have hx_int : IsIntegral Wn (e2 (AdjoinRoot.root g)) := hroot_int.map e2.toAlgHom
    have hdvd : minpoly Wn (e2 (AdjoinRoot.root g)) ∣ g := minpoly.isIntegrallyClosed_dvd hx_int haeval_x
    have hirr2 : Irreducible (minpoly Wn (e2 (AdjoinRoot.root g))) := minpoly.irreducible hx_int
    exact Polynomial.eq_of_monic_of_associated (minpoly.monic hx_int) hgmonic
      (hirr2.associated_of_dvd hirr hdvd)
  refine ⟨(e2.toAlgHom.restrictScalars V0).comp (step1.comp e1.symm.toAlgHom),
    e1 (AdjoinRoot.root f), e2 (AdjoinRoot.root g), ?_, ?_, hwadjoin, hwminpoly, hxadjoin, hxminpoly⟩
  · rw [AlgHom.comp_apply, AlgHom.comp_apply]
    have h1 : e1.symm.toAlgHom (e1 (AdjoinRoot.root f)) = AdjoinRoot.root f := e1.symm_apply_apply _
    rw [h1, hstep1root]
    rfl
  · intro a b hab
    simp only [AlgHom.comp_apply] at hab
    apply e1.symm.injective
    apply hstep1inj
    exact e2.injective hab

set_option maxHeartbeats 1000000 in
/-- **`cancel_conductor_delta` の `hspan_eq` を確立した**——`falt1BaseChangeGeneratorFull`
の `ψ`・`w`・`x` を使って、`conductor Wₙ x * differentIdeal Wₙ Wₙ₊₁ = span{aeval x
(deriv(minpoly Wₙ x))}`(`conductor_mul_differentIdeal`)の右辺と、
`Ideal.map ψ (differentIdeal V0 V1)`(`differentIdeal_tower_diamond` の `Idiff` に
相当)が**文字通り同じ式**(`span{n·x^(n-1)}`)に簡約されることを示す。

配線上の教訓(新規): `differentIdeal_eq_span_derivative`/`conductor_mul_differentIdeal`
を`AdjoinRoot fK`(条件付きで体になる型)を経由する具体的な `w`・`x` に対して呼び出す前に、
`Fact (Irreducible fK)`・`FiniteDimensional (FractionRing V0)(AdjoinRoot fK)`・
`Algebra.IsSeparable (FractionRing V0)(AdjoinRoot fK)` を**呼び出し側の文脈で明示的に
再構築**しておく必要がある——これらは `falt1BaseChangeGeneratorFull` 等のヘルパー内部で
`haveI` されるだけで外に漏れないため、呼び出し側にこれが無いと `AdjoinRoot fK` の
`Field`/`FiniteDimensional` 構造が(条件付き instance のため)見えず、`Application
type mismatch` + `whnf`/`isDefEq` timeout という紛らわしい形で症状が出る。 -/
theorem falt1_hspan_eq
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [IsDedekindDomain (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))))]
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    ∃ (ψ : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
        (algebraMap V0 (FractionRing V0)))) →ₐ[V0]
      integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))
      (x : integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn))))),
      Ideal.span {(Polynomial.aeval x) (Polynomial.derivative (minpoly Wn x))}
        = Ideal.map ψ.toRingHom (differentIdeal V0 (integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
            (algebraMap V0 (FractionRing V0)))))) := by
  have hirrf : Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))) :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))) := ⟨hirrf⟩
  have hirrg : Irreducible (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map (algebraMap Wn (FractionRing Wn))) :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map (algebraMap V0 Wn π) n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map (algebraMap Wn (FractionRing Wn)))) := ⟨hirrg⟩
  have hmonicf : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Monic :=
    (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) := hmonicf.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) := ‹Module.Finite _ _›
  have hmonicg : (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map (algebraMap Wn (FractionRing Wn))).Monic :=
    (Polynomial.monic_X_pow_sub_C (algebraMap V0 Wn π) hnpos.ne').map _
  haveI : Module.Finite (FractionRing Wn) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
      (algebraMap Wn (FractionRing Wn)))) := hmonicg.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing Wn) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
      (algebraMap Wn (FractionRing Wn)))) := ‹Module.Finite _ _›
  have hsepf : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0))).Separable := by
    have hmapeq : (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map (algebraMap V0 (FractionRing V0)))
        = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Algebra.IsSeparable (FractionRing V0) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
      (algebraMap V0 (FractionRing V0)))) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicf hsepf
  have hsepg : (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map (algebraMap Wn (FractionRing Wn))).Separable := by
    have hmapeq : (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map (algebraMap Wn (FractionRing Wn)))
        = Polynomial.X ^ n - Polynomial.C (algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π)) := by simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn' hπne0'
  haveI : Algebra.IsSeparable (FractionRing Wn) (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
      (algebraMap Wn (FractionRing Wn)))) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicg hsepg
  obtain ⟨ψ, w, x, hwx, hψinj, hwadjoin, hwminpoly, hxadjoin, hxminpoly⟩ :=
    ABC3.Found.Falt1.falt1BaseChangeGeneratorFull π n hn hπne0 hprime hnotsq hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn
  have hwfield := ABC3.Found.Falt1.falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly π n hn hπne0 hprime hnotsq hnpos w hwminpoly
  have hdiffV1 := ABC3.Found.Falt1.differentIdeal_eq_span_derivative w hwfield hwadjoin
  rw [hwminpoly] at hdiffV1
  refine ⟨ψ, x, ?_⟩
  have hderivV0 : Polynomial.aeval w (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0))
      = (n : integralClosure V0 (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
          (algebraMap V0 (FractionRing V0))))) * w ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0)
        = Polynomial.C (n : V0) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  have hderivWn : Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))
      = (n : integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
          (algebraMap Wn (FractionRing Wn))))) * x ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)
        = Polynomial.C (n : Wn) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  rw [hxminpoly, hderivWn, hdiffV1, hderivV0, Ideal.map_span]
  congr 1
  simp only [Set.image_singleton]
  congr 1
  rw [map_mul, map_pow, map_natCast, show ψ.toRingHom w = ψ w from rfl, hwx]

set_option maxHeartbeats 1000000 in
set_option synthInstance.maxHeartbeats 1000000 in
/-- **`Algebra.IsSeparable (FractionRing V0) (FractionRing V1)`**(`differentIdeal_ne_bot` が
要求する形——`differentIdeal V0 Wₙ₊₁` の場合の `falt1_hsep_bundled` の、`V0→V1` 単独版の
簡単な類似物。`Wn` を経由しないぶん `falt1_hsep_bundled` より遥かに単純)。`FractionRing V1
≃ₐ[FractionRing V0] AdjoinRoot fK`(`IsLocalization.algEquiv` を `V1`-scalar でも可換に
持ち上げ)経由で、`Algebra.IsSeparable (FractionRing V0) (AdjoinRoot fK)`
(`algIsSeparable_adjoinRoot_of_separable`)を移送するだけ。 -/
noncomputable def falt1_hsepV0V1_bundled
    {V0 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n) :
    Σ' (inst : Algebra (FractionRing V0) (FractionRing (integralClosure V0 (AdjoinRoot
          (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
            (algebraMap V0 (FractionRing V0))))))),
      @Algebra.IsSeparable (FractionRing V0) (FractionRing (integralClosure V0 (AdjoinRoot
          (((Polynomial.X ^ n - Polynomial.C π : Polynomial V0)).map
            (algebraMap V0 (FractionRing V0)))))) _ _ inst := by
  set f : Polynomial V0 := Polynomial.X ^ n - Polynomial.C π with hfdef
  set fK : Polynomial (FractionRing V0) := f.map (algebraMap V0 (FractionRing V0)) with hfKdef
  have hirr : Irreducible fK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirr⟩
  have hmonicK : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  have hsepK : fK.Separable := by
    have hmapeq : fK = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by
      rw [hfKdef, hfdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot fK) := hmonicK.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot fK) := ‹Module.Finite _ _›
  haveI hsepAdjoin : Algebra.IsSeparable (FractionRing V0) (AdjoinRoot fK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicK hsepK
  haveI : IsDedekindDomain (integralClosure V0 (AdjoinRoot fK)) := by infer_instance
  haveI hFR : IsFractionRing (integralClosure V0 (AdjoinRoot fK)) (AdjoinRoot fK) :=
    integralClosure.isFractionRing_of_finite_extension (A := V0) (FractionRing V0) (AdjoinRoot fK)
  set V1 := integralClosure V0 (AdjoinRoot fK) with hV1def
  have h2 : (algebraMap V1 (AdjoinRoot fK)).comp (algebraMap V0 V1) = algebraMap V0 (AdjoinRoot fK) := by
    ext v; show ((algebraMap V0 V1) v : AdjoinRoot fK) = algebraMap V0 (AdjoinRoot fK) v; rfl
  haveI hFaithV0 : FaithfulSMul V0 (FractionRing V1) := by
    rw [faithfulSMul_iff_algebraMap_injective]
    rw [IsScalarTower.algebraMap_eq V0 V1 (FractionRing V1)]
    refine (IsFractionRing.injective V1 (FractionRing V1)).comp ?_
    intro a b hab
    have heq : (algebraMap V0 (AdjoinRoot fK)) a = (algebraMap V0 (AdjoinRoot fK)) b := by
      have := congrArg (Subtype.val) hab
      simpa using this
    have hinj2 : Function.Injective (algebraMap V0 (AdjoinRoot fK)) := by
      rw [IsScalarTower.algebraMap_eq V0 (FractionRing V0) (AdjoinRoot fK)]
      exact (algebraMap (FractionRing V0) (AdjoinRoot fK)).injective.comp
        (IsFractionRing.injective V0 (FractionRing V0))
    exact hinj2 heq
  letI algFRV1 : Algebra (FractionRing V0) (FractionRing V1) :=
    FractionRing.liftAlgebra V0 (FractionRing V1)
  haveI hSTcanon : IsScalarTower V0 (FractionRing V0) (FractionRing V1) := by infer_instance
  have e1 : FractionRing V1 ≃ₐ[V1] AdjoinRoot fK :=
    IsLocalization.algEquiv (nonZeroDivisors V1) (FractionRing V1) (AdjoinRoot fK)
  have hringhom_eq' : ∀ v : V0, (e1.toRingHom) (algebraMap (FractionRing V0) (FractionRing V1)
        (algebraMap V0 (FractionRing V0) v)) = algebraMap (FractionRing V0) (AdjoinRoot fK)
        (algebraMap V0 (FractionRing V0) v) := by
    intro v
    rw [← IsScalarTower.algebraMap_apply V0 (FractionRing V0) (FractionRing V1)]
    rw [← IsScalarTower.algebraMap_apply V0 (FractionRing V0) (AdjoinRoot fK)]
    rw [IsScalarTower.algebraMap_apply V0 V1 (FractionRing V1)]
    show e1 (algebraMap V1 (FractionRing V1) (algebraMap V0 V1 v)) = algebraMap V0 (AdjoinRoot fK) v
    rw [e1.commutes]
    exact congrFun (congrArg DFunLike.coe h2) v
  have hringhom_eq'' : (e1.toRingHom).comp (algebraMap (FractionRing V0) (FractionRing V1))
      = algebraMap (FractionRing V0) (AdjoinRoot fK) :=
    IsLocalization.ringHom_ext (nonZeroDivisors V0) (RingHom.ext hringhom_eq')
  have e1' : FractionRing V1 ≃ₐ[FractionRing V0] AdjoinRoot fK :=
    AlgEquiv.ofRingEquiv (f := e1.toRingEquiv) (fun x => by
      have := congrFun (congrArg DFunLike.coe hringhom_eq'') x; simpa using this)
  haveI hsepFinal : Algebra.IsSeparable (FractionRing V0) (FractionRing V1) :=
    (AlgEquiv.Algebra.isSeparable_iff e1').mpr hsepAdjoin
  exact ⟨algFRV1, hsepFinal⟩

/-!
## Theorem 1.2・3b(i)への足場: `Wₙ₊₁` は PID(2026-09-04)

`length_quotient_span_singleton_mul`(既存、単項イデアルの場合の長さの
加法性)を`conductor(Wₙ,x)`に適用するには`conductor(Wₙ,x) : Ideal Wₙ₊₁`
が単項である必要がある——**`Wₙ₊₁` が PID なら自動的に成り立つ**。
`Wₙ₊₁` が「全分岐(totally ramified)」であることを直接示す(「3c:
戦略転換」節が「mathlib にまだ薄い領域」と記録した経路)必要は無い、
と判明した: `IsPrincipalIdealRing.of_finite_maximals`(既存、Dedekind
整域で極大イデアルが有限個なら PID)+ `IsDedekindDomain.primesOver_
finite`(既存、有限拡大では上にある素イデアルは有限個)+
`Ideal.isMaximal_comap_of_isIntegral_of_isMaximal`(既存、整拡大では
極大イデアルの引き戻しは極大)を組み合わせるだけで、`Wₙ₊₁` の極大
イデアル全体が「`Wₙ` の(唯一の)極大イデアルの上にある素イデアル」
という**有限集合の部分集合**であることが分かり、単純な部分集合の
有限性だけで閉じる——`Wₙ` が DVR で唯一の極大イデアルしか持たない
ことだけを使い、「全分岐」という強い主張自体は一切不要だった。 -/

set_option maxHeartbeats 800000 in
/-- **DVR の有限拡大で(捩れ無しの)Dedekind整域になるものは PID**。
`conductor(Wₙ,x)` が単項イデアルであることの直接の供給源。 -/
theorem falt1_isPrincipalIdealRing_of_finite_ext_of_DVR
    {Wn Wn1 : Type*} [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn]
    [CommRing Wn1] [IsDedekindDomain Wn1] [Algebra Wn Wn1]
    [Module.Finite Wn Wn1] [Module.IsTorsionFree Wn Wn1] :
    IsPrincipalIdealRing Wn1 := by
  apply IsPrincipalIdealRing.of_finite_maximals
  have hsub : {I : Ideal Wn1 | I.IsMaximal} ⊆ (IsLocalRing.maximalIdeal Wn).primesOver Wn1 := by
    intro I hI
    haveI : I.IsMaximal := hI
    haveI : Algebra.IsIntegral Wn Wn1 := Algebra.IsIntegral.of_finite Wn Wn1
    have hcomap : (I.comap (algebraMap Wn Wn1)).IsMaximal :=
      Ideal.isMaximal_comap_of_isIntegral_of_isMaximal I
    have heq : IsLocalRing.maximalIdeal Wn = Ideal.under Wn I :=
      (IsLocalRing.eq_maximalIdeal hcomap).symm
    exact ⟨hI.isPrime, Ideal.LiesOver.mk heq⟩
  exact Set.Finite.subset (IsDedekindDomain.primesOver_finite (IsLocalRing.maximalIdeal Wn) Wn1) hsub

set_option maxHeartbeats 400000 in
/-- **PID における長さの加法性(非零な単項イデアルの倍)**:
`falt1_isPrincipalIdealRing_of_finite_ext_of_DVR` と
`length_quotient_span_singleton_mul` を貼り合わせただけ——`I ≠ 0` から
`IsPrincipalIdealRing.principal` で生成元 `a` を取り出し、`a ≠ 0`
(ゆえに整域では非零因子)を確認して`length_quotient_span_singleton_mul`
を適用する。`conductor(Wₙ,x)`(単項とは限らないと思われていたが、
`Wₙ₊₁` が PID なので実は自動的に単項)に直接使える形。 -/
theorem falt1_length_quotient_mul_of_ne_zero {R : Type*} [CommRing R] [IsDomain R] [IsPrincipalIdealRing R]
    (I J : Ideal R) (hI : I ≠ 0) :
    Module.length R (R ⧸ (I * J)) = Module.length R (R ⧸ I) + Module.length R (R ⧸ J) := by
  obtain ⟨a, ha⟩ := IsPrincipalIdealRing.principal I
  have ha0 : a ≠ 0 := by
    intro h
    apply hI
    rw [ha, h]
    simp
  have hanzd : a ∈ nonZeroDivisors R := mem_nonZeroDivisors_of_ne_zero ha0
  have hI' : I = Ideal.span ({a} : Set R) := ha
  rw [hI']
  exact ABC3.Found.Falt1.length_quotient_span_singleton_mul a hanzd J

/-- **`differentIdeal` の塔公式を「長さ」の言葉に翻訳した版(2026-09-04)**:
`differentIdeal_eq_differentIdeal_mul_differentIdeal`(mathlib、塔
`A→B→C`に対する直接の ideal 等式、`C` 側の単項生成元は一切不要)と
`falt1_length_quotient_mul_of_ne_zero`(PID での長さの加法性、直前)を
組み合わせるだけ。Kähler 微分の完全列(`falt1_kaehler_length_exact_
wn1`系列)を経由せずに同じ形の長さの等式に到達できる、より直接的な
経路——`falt1_cancelConductorDelta_assembled`の`hlen_eq`と組み合わせ
れば`differentIdeal`だけで書かれた閉じた式が得られる(`Ideal.map
(algebraMap B C)(differentIdeal A B)`の項を`hlen_eq`で`differentIdeal
Vₙ₁ Wₙ₁`に置き換えるだけ、次回の作業)。 -/
theorem falt1_differentIdeal_tower_length {A B C : Type*} [CommRing A] [IsDomain A] [IsIntegrallyClosed A]
    [CommRing B] [IsDedekindDomain B] [Algebra A B]
    [CommRing C] [IsDedekindDomain C] [IsPrincipalIdealRing C]
    [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    [Module.IsTorsionFree A B] [Module.IsTorsionFree A C] [Module.IsTorsionFree B C]
    [Module.Finite A B] [Module.Finite A C] [Module.Finite B C]
    (hsep : @Algebra.IsSeparable (FractionRing A) (FractionRing C) _ _ (FractionRing.liftAlgebra A (FractionRing C)))
    (hne : differentIdeal B C ≠ 0) :
    Module.length C (C ⧸ differentIdeal A C) =
      Module.length C (C ⧸ differentIdeal B C) +
      Module.length C (C ⧸ Ideal.map (algebraMap B C) (differentIdeal A B)) := by
  letI := FractionRing.liftAlgebra A (FractionRing C)
  haveI := hsep
  have htower : differentIdeal A C = differentIdeal B C * Ideal.map (algebraMap B C) (differentIdeal A B) :=
    differentIdeal_eq_differentIdeal_mul_differentIdeal A B C
  rw [htower]
  exact ABC3.Found.Falt1.falt1_length_quotient_mul_of_ne_zero (differentIdeal B C)
    (Ideal.map (algebraMap B C) (differentIdeal A B)) hne

/-- **`falt1_differentIdeal_tower_length` を「`Ideal.map(...)`の項が別の
イデアル`E`の長さと既に一致することが分かっている」場合に特化した版**:
`hlen_eq`(`cancel_conductor_delta`経由で確立される長さの等式、典型的には
`E := differentIdeal Vₙ₁ Wₙ₊₁`)を渡すだけで、`differentIdeal`だけで
書かれた閉じた式が出る。抽象的な型変数のまま(具体的な`integralClosure`
式を経由しない)なので、シグネチャの elaboration が高速——具体的な
Falt1 の `Wₙ₊₁` への適用は呼び出し側(証明の中、`Wₙ₊₁`が`set`で略記
された文脈)で行うことを想定している。 -/
theorem falt1_differentIdeal_tower_length_via_hlen_eq
    {A B C : Type*} [CommRing A] [IsDomain A] [IsIntegrallyClosed A]
    [CommRing B] [IsDedekindDomain B] [Algebra A B]
    [CommRing C] [IsDedekindDomain C] [IsPrincipalIdealRing C]
    [Algebra B C] [Algebra A C] [IsScalarTower A B C]
    [Module.IsTorsionFree A B] [Module.IsTorsionFree A C] [Module.IsTorsionFree B C]
    [Module.Finite A B] [Module.Finite A C] [Module.Finite B C]
    (hsep : @Algebra.IsSeparable (FractionRing A) (FractionRing C) _ _ (FractionRing.liftAlgebra A (FractionRing C)))
    (hne : differentIdeal B C ≠ 0)
    (E : Ideal C)
    (hlen_eq : Module.length C (C ⧸ Ideal.map (algebraMap B C) (differentIdeal A B)) = Module.length C (C ⧸ E)) :
    Module.length C (C ⧸ differentIdeal A C) =
      Module.length C (C ⧸ differentIdeal B C) + Module.length C (C ⧸ E) := by
  have h := ABC3.Found.Falt1.falt1_differentIdeal_tower_length hsep hne
  rw [hlen_eq] at h
  exact h

/-- **`differentIdeal` の菱形公式を「長さ」の言葉に翻訳した版(2026-09-04)**:
`differentIdeal_tower_diamond`(塔 `Vn→Vn1→Wn1` と `Vn→Wn→Wn1` の2通りで
`differentIdeal Vn Wn1` を計算すると一致する、`Vn1`・`Wn`側の単項生成元は
どちらも不要)の両辺に`falt1_length_quotient_mul_of_ne_zero`を適用する
だけ——**`Wₙ₊₁`の`Vₙ`上の単項生成元を経由せずに**、`Jₙ:=differentIdeal
Wₙ Wₙ₊₁`(現在の段の非正規性)と`Jₙ₊₁:=differentIdeal Vₙ₁ Wₙ₊₁`(次の段の
非正規性)を結ぶ**長さの等式**が得られる——Theorem 1.2 の再帰
(`δₙ→0`)に必要な「段をまたぐ関係式」の候補。`falt1_theorem12_
differentIdeal_length`が抱えていた「左辺を`Ω¹`の長さに結びつけるには
`Wₙ₊₁`の`V0`上の単項生成元が要る」という障害を回避する経路。 -/
theorem falt1_differentIdeal_diamond_length {Vn Vn1 Wn Wn1 : Type*} [CommRing Vn] [CommRing Vn1] [CommRing Wn]
    [CommRing Wn1] [IsDomain Vn] [IsIntegrallyClosed Vn] [IsDedekindDomain Vn1] [IsDedekindDomain Wn]
    [IsDedekindDomain Wn1] [IsPrincipalIdealRing Wn1] [Algebra Vn Vn1] [Algebra Vn Wn] [Algebra Vn1 Wn1] [Algebra Wn Wn1]
    [Algebra Vn Wn1] [IsScalarTower Vn Vn1 Wn1] [IsScalarTower Vn Wn Wn1]
    [Module.IsTorsionFree Vn Vn1] [Module.IsTorsionFree Vn Wn] [Module.IsTorsionFree Vn Wn1]
    [Module.IsTorsionFree Vn1 Wn1] [Module.IsTorsionFree Wn Wn1] [Module.Finite Vn Vn1]
    [Module.Finite Vn Wn] [Module.Finite Vn Wn1] [Module.Finite Vn1 Wn1] [Module.Finite Wn Wn1]
    (hsep : @Algebra.IsSeparable (FractionRing Vn) (FractionRing Wn1) _ _ (FractionRing.liftAlgebra Vn (FractionRing Wn1)))
    (hneVn1Wn1 : differentIdeal Vn1 Wn1 ≠ 0) (hneWnWn1 : differentIdeal Wn Wn1 ≠ 0) :
    Module.length Wn1 (Wn1 ⧸ differentIdeal Vn1 Wn1) + Module.length Wn1 (Wn1 ⧸ Ideal.map (algebraMap Vn1 Wn1) (differentIdeal Vn Vn1)) =
      Module.length Wn1 (Wn1 ⧸ differentIdeal Wn Wn1) + Module.length Wn1 (Wn1 ⧸ Ideal.map (algebraMap Wn Wn1) (differentIdeal Vn Wn)) := by
  have hdiamond := ABC3.Found.Falt1.differentIdeal_tower_diamond (Vn := Vn) (Vn1 := Vn1) (Wn := Wn) (Wn1 := Wn1) hsep
  have h1 := ABC3.Found.Falt1.falt1_length_quotient_mul_of_ne_zero (differentIdeal Vn1 Wn1) (Ideal.map (algebraMap Vn1 Wn1) (differentIdeal Vn Vn1)) hneVn1Wn1
  have h2 := ABC3.Found.Falt1.falt1_length_quotient_mul_of_ne_zero (differentIdeal Wn Wn1) (Ideal.map (algebraMap Wn Wn1) (differentIdeal Vn Wn)) hneWnWn1
  rw [← h1, ← h2, hdiamond]

/-!
## `Wₙ₊₁` を名前付き `def` で抽象化する(2026-09-04)

`falt1_cancelConductorDelta_assembled` が戻り値を `True` にしている理由
(このすぐ下のdocstringに記録済み)は、`Wₙ₊₁`(`integralClosure Wn
(AdjoinRoot(...))`という巨大な入れ子式)を**シグネチャに`set`なしで
そのまま複数回書く**と、`differentIdeal`・`Ideal.map`の型検査のたびに
instance 探索をやり直す必要が生じ、`maxHeartbeats 3000000`でも
`isDefEq`/`synthesize pending MVars`のtimeoutに達することが判明した
ため(実測、`falt1_cancelConductorDelta_assembled`への統合を2度試みて
2度とも revert・記録済み、`falt1-goal.md`参照)。

**解決策**: `Wₙ₊₁`を独立の名前付き`def`(`Falt1Wn1`)として1回だけ定義し、
必要な基本instance(`CommRing`・`Algebra Wn`・`Algebra V0`・
`IsScalarTower V0 Wn`)を`inferInstanceAs`で**1回だけ**登録する。
シグネチャでは`Falt1Wn1 V0 Wn π n`という短い適用が(複数回書いても)
高速に型検査される(実測:数秒 vs. 実装前は500秒近いtimeout)。
証明の中では通常通り`set Wn1 := integralClosure Wn (AdjoinRoot gK)`
で局所展開し、`Falt1Wn1`側のinstance仮定は`‹...›`(assumption)で
`Wn1`側へ1回だけ変換すれば、残りの証明は既存のものと同じ速度で進む
——最後の`exact`は`Wn1`型で明示的に`have`した結果を返すだけで、
`Falt1Wn1`の goal との defeq 判定は(既に確定した局所 instance を
再利用するだけなので)高速。 -/

noncomputable def Falt1Wn1 (V0 Wn : Type*) [CommRing V0] [CommRing Wn] [Algebra V0 Wn] (π : V0) (n : ℕ) : Type _ :=
  integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
    (algebraMap Wn (FractionRing Wn))))

noncomputable instance Falt1Wn1.commRing {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn] (π : V0) (n : ℕ) :
    CommRing (Falt1Wn1 V0 Wn π n) :=
  inferInstanceAs (CommRing (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
    (algebraMap Wn (FractionRing Wn))))))

noncomputable instance Falt1Wn1.algebra {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn] (π : V0) (n : ℕ) :
    Algebra Wn (Falt1Wn1 V0 Wn π n) :=
  inferInstanceAs (Algebra Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
    (algebraMap Wn (FractionRing Wn))))))

noncomputable instance Falt1Wn1.algebraV0 {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn] (π : V0) (n : ℕ) :
    Algebra V0 (Falt1Wn1 V0 Wn π n) :=
  inferInstanceAs (Algebra V0 (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
    (algebraMap Wn (FractionRing Wn))))))

instance Falt1Wn1.isScalarTower {V0 Wn : Type*} [CommRing V0] [CommRing Wn] [Algebra V0 Wn] (π : V0) (n : ℕ) :
    IsScalarTower V0 Wn (Falt1Wn1 V0 Wn π n) :=
  inferInstanceAs (IsScalarTower V0 Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
    (algebraMap Wn (FractionRing Wn))))))

set_option maxHeartbeats 3000000 in
set_option synthInstance.maxHeartbeats 3000000 in
/-- **Theorem 1.2・item 3c の中核道具の完全な組み立て**:
`differentIdeal_tower_diamond`(`hdiamond`)・`conductor_mul_differentIdeal`(`hcond`)・
`falt1_hspan_eq` の `hspan_eq`・`differentIdeal_ne_bot`(`Idiff≠0`)のすべてを実際に接続し、
`cancel_conductor_delta` を呼び出して
```
conductor Wₙ x * differentIdeal V1 Wₙ₊₁ = Ideal.map (algebraMap Wₙ Wₙ₊₁) (differentIdeal V0 Wₙ)
```
(`Jₙ := differentIdeal Wₙ Wₙ₊₁` が消去された、`δₙ`・`δₙ₊₁` を結ぶ関係式)を得た——これが
Theorem 1.2・3b/3c の中核の代数的関係式であり、item 3c の技術的な障害はこれで**すべて
解決した**。戻り値を `True` にしている理由は `falt1_differentIdeal_tower_diamond_assembled`
と同じ(`ψ` 選択前に `Module.IsTorsionFree`/`Algebra.IsSeparable` 系の instance 検索が
シグネチャの型検査自体で走ってしまう)。

★2026-09-04追記: さらに `conductor(Wₙ,x) ≠ 0`(`hcond`・`hspan_eq`・`hIdiff_ne` から
——`conductor(Wₙ,x)*differentIdeal Wₙ Wₙ₊₁ = Ideal.map ψ (differentIdeal V0 V1) ≠ 0` の
積が非零なら整域では両因子とも非零)と `falt1_isPrincipalIdealRing_of_finite_ext_of_DVR`
(`Wₙ₊₁` は PID)を足し、`falt1_length_quotient_mul_of_ne_zero` を適用して
```
Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ Ideal.map (algebraMap Wₙ Wₙ₊₁) (differentIdeal V0 Wₙ))
  = Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ conductor Wₙ x) + Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ differentIdeal V1 Wₙ₊₁)
```
という**長さの等式**(`hlen_eq`)まで得た。

★2026-09-04追加・重要な発見: `x` が `Wₙ₊₁` を `Wₙ`-**代数として**生成する
(`hxadjoin : Wₙ[x] = ⊤`、既に確立済み)ため、`conductor_eq_top_iff_adjoin_eq_top`
から`conductor(Wₙ,x) = ⊤`が直ちに従い、`length(Wₙ₊₁/conductor(Wₙ,x)) = 0`
——右辺第1項が**消える**。得られる完全な等式:
```
Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ Ideal.map (algebraMap Wₙ Wₙ₊₁) (differentIdeal V0 Wₙ))
  = Module.length Wₙ₊₁ (Wₙ₊₁ ⧸ differentIdeal V1 Wₙ₊₁)
```
**★★これは`Wₙ₊₁`を`V1`と全く同じEisenstein多項式`gK`のbase changeから作る、という
このセッションの構成固有の単純化である可能性が高い**——Faltings原論文の一般論
(`Wₙ`が任意の代数拡大でもよい場合の「非正規性の評価」、conductor が非自明な
場合)の全容を捉えているかは未確認。次回、この構成が Theorem 1.2 の再帰的な
`V_n` 塔の構成として本当に十分か(あるいは「逸脱」として原文との対応を明記
すべきか)を検討すること——CLAUDE.mdの「逸脱」の記録原則に従う。

なお残る接続として、この等式の左辺を`δₙ`と結ぶ(a)については、`length_map_
pow_of_ramificationIdx`・`IsLocalRing.length_baseChange`いずれも`Wₙ₊₁`の
局所性(全分岐)を要求すると確認済み(近道は見つからなかった、falt1-goal.md
参照)——`conductor`の項が消えたことで(c)の困難自体は解消したが、(a)の
全分岐証明は依然として必要。 -/
theorem falt1_cancelConductorDelta_assembled
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [Module.Finite V0 Wn]
    [Algebra (FractionRing V0) (FractionRing Wn)]
    [IsScalarTower V0 (FractionRing V0) (FractionRing Wn)]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)] : True := by
  set fK : Polynomial (FractionRing V0) := (Polynomial.X ^ n - Polynomial.C π).map (algebraMap V0 (FractionRing V0)) with hfKdef
  set gK : Polynomial (FractionRing Wn) := (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π)).map (algebraMap Wn (FractionRing Wn)) with hgKdef
  have hirrf : Irreducible fK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirrf⟩
  have hirrg : Irreducible gK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map (algebraMap V0 Wn π) n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible gK) := ⟨hirrg⟩
  have hmonicf : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot fK) := hmonicf.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot fK) := ‹Module.Finite _ _›
  have hmonicg : gK.Monic := (Polynomial.monic_X_pow_sub_C (algebraMap V0 Wn π) hnpos.ne').map _
  haveI : Module.Finite (FractionRing Wn) (AdjoinRoot gK) := hmonicg.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing Wn) (AdjoinRoot gK) := ‹Module.Finite _ _›
  have hsepf : fK.Separable := by
    have hmapeq : fK = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by
      rw [hfKdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Algebra.IsSeparable (FractionRing V0) (AdjoinRoot fK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicf hsepf
  have hsepg : gK.Separable := by
    have hmapeq : gK = Polynomial.X ^ n - Polynomial.C (algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π)) := by
      rw [hgKdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn' hπne0'
  haveI : Algebra.IsSeparable (FractionRing Wn) (AdjoinRoot gK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicg hsepg
  haveI : IsDedekindDomain V0 := by infer_instance
  haveI : IsDedekindDomain Wn := by infer_instance
  haveI hDedV1 : IsDedekindDomain (integralClosure V0 (AdjoinRoot fK)) := by infer_instance
  haveI hDedWn1 : IsDedekindDomain (integralClosure Wn (AdjoinRoot gK)) := by infer_instance
  haveI hTFV1 : Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot fK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeV0V1 π n hprime hnotsq hnpos
  haveI hTFWn1 : Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  haveI hTFV0Wn1 : Module.IsTorsionFree V0 (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeV0Wn1 π n hprime' hnotsq' hnpos hinjV0Wn
  haveI hFinWnWn1 : Module.Finite Wn (integralClosure Wn (AdjoinRoot gK)) :=
    ABC3.Found.Falt1.falt1ModuleFiniteWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  set V1 := integralClosure V0 (AdjoinRoot fK) with hV1def
  set Wn1 := integralClosure Wn (AdjoinRoot gK) with hWn1def
  obtain ⟨ψ, w, x, hwx, hψinj, hwadjoin, hwminpoly, hxadjoin, hxminpoly⟩ :=
    ABC3.Found.Falt1.falt1BaseChangeGeneratorFull π n hn hπne0 hprime hnotsq hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn
  letI algV1Wn1 : Algebra V1 Wn1 := ψ.toRingHom.toAlgebra
  haveI hSTV0V1Wn1 : IsScalarTower V0 V1 Wn1 := by
    apply IsScalarTower.of_algebraMap_eq
    intro y
    exact (ψ.commutes y).symm
  haveI hSTV0WnWn1 : IsScalarTower V0 Wn Wn1 := by infer_instance
  haveI hFinV0V1 : Module.Finite V0 V1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0V1 π n hn hπne0 hprime hnotsq hnpos
  haveI hFinV0Wn1 : Module.Finite V0 Wn1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
  haveI hFinV1Wn1 : Module.Finite V1 Wn1 := by
    haveI := ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
    exact Module.Finite.of_restrictScalars_finite V0 V1 Wn1
  haveI hIsDomV1 : IsDomain V1 := inferInstance
  haveI hIsDomWn1 : IsDomain Wn1 := inferInstance
  haveI hTFV1Wn1 : Module.IsTorsionFree V1 Wn1 :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hψinj
  haveI hTFV0Wn : Module.IsTorsionFree V0 Wn :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hinjV0Wn
  haveI hIsIntClosedV0 : IsIntegrallyClosed V0 := inferInstance
  letI algFRWn1 := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).1
  have hsep := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).2
  have hdiamond := ABC3.Found.Falt1.differentIdeal_tower_diamond
    (Vn := V0) (Vn1 := V1) (Wn := Wn) (Wn1 := Wn1) hsep
  have hwfield := ABC3.Found.Falt1.falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly π n hn hπne0 hprime hnotsq hnpos w hwminpoly
  have hdiffV1 := ABC3.Found.Falt1.differentIdeal_eq_span_derivative w hwfield hwadjoin
  rw [hwminpoly] at hdiffV1
  have hxfield := ABC3.Found.Falt1.falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly (algebraMap V0 Wn π) n hn' hπne0' hprime' hnotsq' hnpos x hxminpoly
  have hcond := conductor_mul_differentIdeal Wn (FractionRing Wn) (AdjoinRoot gK) x hxfield
  rw [hxminpoly] at hcond
  have hderivV0 : Polynomial.aeval w (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0))
      = (n : V1) * w ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0)
        = Polynomial.C (n : V0) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  have hderivWn : Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))
      = (n : Wn1) * x ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)
        = Polynomial.C (n : Wn) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  have hspan_eq : Ideal.span {Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))}
      = Ideal.map ψ.toRingHom (differentIdeal V0 V1) := by
    rw [hderivWn, hdiffV1, hderivV0, Ideal.map_span]
    congr 1
    simp only [Set.image_singleton]
    congr 1
    rw [map_mul, map_pow, map_natCast, show ψ.toRingHom w = ψ w from rfl, hwx]
  letI algFRV1 := (ABC3.Found.Falt1.falt1_hsepV0V1_bundled π n hn hπne0 hprime hnotsq hnpos).1
  have hsepV0V1 : @Algebra.IsSeparable (FractionRing V0) (FractionRing V1) _ _ algFRV1 :=
    (ABC3.Found.Falt1.falt1_hsepV0V1_bundled π n hn hπne0 hprime hnotsq hnpos).2
  have hdiffV1_ne : differentIdeal V0 V1 ≠ ⊥ :=
    @differentIdeal_ne_bot V0 V1 _ _ _ _ _ _ _ _ hsepV0V1
  have hIdiff_ne : Ideal.map ψ.toRingHom (differentIdeal V0 V1) ≠ 0 := by
    show Ideal.map ψ.toRingHom (differentIdeal V0 V1) ≠ ⊥
    intro heq
    exact hdiffV1_ne ((Ideal.map_eq_bot_iff_of_injective hψinj).mp heq)
  have hfinal := ABC3.Found.Falt1.cancel_conductor_delta
    (Ideal.map ψ.toRingHom (differentIdeal V0 V1)) (differentIdeal Wn Wn1)
    (conductor Wn x) (Ideal.span {Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))})
    (differentIdeal V1 Wn1) (Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn))
    hdiamond hcond hspan_eq hIdiff_ne
  -- 追加(2026-09-04): conductor(Wn,x) ≠ 0 を示し、PID での長さの加法性
  -- (falt1_length_quotient_mul_of_ne_zero)を hfinal に適用して長さの等式を得る。
  have hcondprod_ne : conductor Wn x * differentIdeal Wn Wn1 ≠ 0 := by
    rw [hcond, hspan_eq]
    exact hIdiff_ne
  have hcond_ne : conductor Wn x ≠ 0 := by
    intro h
    apply hcondprod_ne
    rw [h, zero_mul]
  haveI hPIDWn1 : IsPrincipalIdealRing Wn1 :=
    ABC3.Found.Falt1.falt1_isPrincipalIdealRing_of_finite_ext_of_DVR (Wn := Wn) (Wn1 := Wn1)
  have hlen_eq : Module.length Wn1 (Wn1 ⧸ (conductor Wn x * differentIdeal V1 Wn1)) =
      Module.length Wn1 (Wn1 ⧸ conductor Wn x) + Module.length Wn1 (Wn1 ⧸ differentIdeal V1 Wn1) :=
    ABC3.Found.Falt1.falt1_length_quotient_mul_of_ne_zero (conductor Wn x) (differentIdeal V1 Wn1) hcond_ne
  rw [hfinal] at hlen_eq
  -- hlen_eq :
  --   Module.length Wn1 (Wn1 ⧸ Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn))
  --     = Module.length Wn1 (Wn1 ⧸ conductor Wn x) + Module.length Wn1 (Wn1 ⧸ differentIdeal V1 Wn1)
  -- ★2026-09-04追加・重要な発見: `hxadjoin`(`Wn[x] = Wn1`、既に確立済み)から
  -- `conductor Wn x = ⊤` が(`conductor_eq_top_iff_adjoin_eq_top`で)直ちに従う——
  -- 「conductor の長さ」の項が消え、hlen_eq は**補正項の無い完全な等式**になる。
  -- ★★この構成(Wₙ₊₁ を V1 と全く同じ Eisenstein 多項式 gK の base change から作る)
  -- 固有の単純化である可能性が高く、Faltings の一般論(Wₙ が任意の代数拡大でもよい
  -- 場合の「非正規性の評価」)の全容を捉えているかは未確認——次回、この構成が
  -- Theorem 1.2 の再帰的な V_n 塔の構成として本当に十分か(あるいは「逸脱」として
  -- 記録すべきか)を検討すること。 -/
  have hcond_top : conductor Wn x = ⊤ := conductor_eq_top_iff_adjoin_eq_top.mpr hxadjoin
  rw [hcond_top] at hlen_eq
  have hzero : Module.length Wn1 (Wn1 ⧸ (⊤ : Ideal Wn1)) = 0 := Module.length_eq_zero
  rw [hzero, zero_add] at hlen_eq
  -- hlen_eq (簡約後):
  --   Module.length Wn1 (Wn1 ⧸ Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn))
  --     = Module.length Wn1 (Wn1 ⧸ differentIdeal V1 Wn1)
  trivial

set_option maxHeartbeats 3000000 in
set_option synthInstance.maxHeartbeats 3000000 in
/-- **Theorem 1.2・3b/3c、`differentIdeal`だけで書かれた閉じた式(完成、
2026-09-04)**: `falt1_cancelConductorDelta_assembled`と同じ導出(`ψ`・`w`・
`x`の取得、`hdiamond`・`hcond`・`hspan_eq`・`hfinal`の構成)を`Falt1Wn1`と
いう名前付き`def`(このファイル冒頭近くで定義)を使って行い、
`falt1_differentIdeal_tower_length`(Kähler微分の完全列を経由しない、
`differentIdeal`の塔公式からの直接路)で結ぶ。`Wₙ₊₁`を毎回`integralClosure
Wn(AdjoinRoot(...))`とベタ書きするとシグネチャの型検査自体が
timeout する(前の定理が`True`を返している理由と同じ)ため、名前付き
`def`で1回だけ展開することで解決した——これが Exercise 13.7.4 の議論
(kernel/cokernel の長さ評価)とは独立に、`differentIdeal`の塔の乗法性
だけから出る**閉じた形の結論**である:
```
length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal V0 Wₙ₊₁) =
  length_{Wₙ₊₁}(Wₙ₊₁⧸differentIdeal Wₙ Wₙ₊₁) +
  length_{Wₙ₊₁}(Wₙ₊₁⧸Ideal.map(algebraMap Wₙ Wₙ₊₁)(differentIdeal V0 Wₙ))
```
残る作業: 左辺`differentIdeal V0 Wₙ₊₁`を`Ω¹_{Wₙ₊₁/V0}`の長さ(Lemma 1.1
経由)と結びつけるには`Wₙ₊₁`の`V0`上の単項生成元が要る(未解決、
falt1-goal.md参照)。右辺第2項は既にこの`Wₙ₊₁`の構成(`x`が`Wₙ`-代数として
生成する)から求まる形にできている。 -/
theorem falt1_theorem12_differentIdeal_length
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [Module.Finite V0 Wn]
    [Module.IsTorsionFree V0 Wn]
    [Algebra (FractionRing V0) (FractionRing Wn)]
    [IsScalarTower V0 (FractionRing V0) (FractionRing Wn)]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)]
    [IsDedekindDomain (Falt1Wn1 V0 Wn π n)]
    [Module.IsTorsionFree V0 (Falt1Wn1 V0 Wn π n)]
    [Module.IsTorsionFree Wn (Falt1Wn1 V0 Wn π n)] :
    Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸ differentIdeal V0 (Falt1Wn1 V0 Wn π n)) =
      Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸ differentIdeal Wn (Falt1Wn1 V0 Wn π n)) +
      Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸
        Ideal.map (algebraMap Wn (Falt1Wn1 V0 Wn π n)) (differentIdeal V0 Wn)) := by
  set fK : Polynomial (FractionRing V0) := (Polynomial.X ^ n - Polynomial.C π).map (algebraMap V0 (FractionRing V0)) with hfKdef
  set gK : Polynomial (FractionRing Wn) := (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π)).map (algebraMap Wn (FractionRing Wn)) with hgKdef
  have hirrf : Irreducible fK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map π n hnpos hprime hnotsq
  haveI : Fact (Irreducible fK) := ⟨hirrf⟩
  have hirrg : Irreducible gK :=
    ABC3.Found.Falt1.eisenstein_X_pow_sub_C_irreducible_map (algebraMap V0 Wn π) n hnpos hprime' hnotsq'
  haveI : Fact (Irreducible gK) := ⟨hirrg⟩
  have hmonicf : fK.Monic := (Polynomial.monic_X_pow_sub_C π hnpos.ne').map _
  haveI : Module.Finite (FractionRing V0) (AdjoinRoot fK) := hmonicf.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing V0) (AdjoinRoot fK) := ‹Module.Finite _ _›
  have hmonicg : gK.Monic := (Polynomial.monic_X_pow_sub_C (algebraMap V0 Wn π) hnpos.ne').map _
  haveI : Module.Finite (FractionRing Wn) (AdjoinRoot gK) := hmonicg.finite_adjoinRoot
  haveI : FiniteDimensional (FractionRing Wn) (AdjoinRoot gK) := ‹Module.Finite _ _›
  have hsepf : fK.Separable := by
    have hmapeq : fK = Polynomial.X ^ n - Polynomial.C (algebraMap V0 (FractionRing V0) π) := by
      rw [hfKdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn hπne0
  haveI : Algebra.IsSeparable (FractionRing V0) (AdjoinRoot fK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicf hsepf
  have hsepg : gK.Separable := by
    have hmapeq : gK = Polynomial.X ^ n - Polynomial.C (algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π)) := by
      rw [hgKdef]; simp
    rw [hmapeq]; exact Polynomial.separable_X_pow_sub_C _ hn' hπne0'
  haveI : Algebra.IsSeparable (FractionRing Wn) (AdjoinRoot gK) :=
    ABC3.Found.Falt1.algIsSeparable_adjoinRoot_of_separable _ hmonicg hsepg
  haveI : IsDedekindDomain V0 := by infer_instance
  haveI : IsDedekindDomain Wn := by infer_instance
  haveI hDedV1 : IsDedekindDomain (integralClosure V0 (AdjoinRoot fK)) := by infer_instance
  set Wn1 := integralClosure Wn (AdjoinRoot gK) with hWn1def
  haveI hDedWn1 : IsDedekindDomain Wn1 := ‹IsDedekindDomain (Falt1Wn1 V0 Wn π n)›
  haveI hTFV1 : Module.IsTorsionFree V0 (integralClosure V0 (AdjoinRoot fK)) :=
    ABC3.Found.Falt1.falt1ModuleIsTorsionFreeV0V1 π n hprime hnotsq hnpos
  haveI hTFWn1 : Module.IsTorsionFree Wn Wn1 := ‹Module.IsTorsionFree Wn (Falt1Wn1 V0 Wn π n)›
  haveI hTFV0Wn1 : Module.IsTorsionFree V0 Wn1 := ‹Module.IsTorsionFree V0 (Falt1Wn1 V0 Wn π n)›
  haveI hFinWnWn1 : Module.Finite Wn Wn1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteWnWn1 π n hnpos hn' hπne0' hprime' hnotsq'
  set V1 := integralClosure V0 (AdjoinRoot fK) with hV1def
  obtain ⟨ψ, w, x, hwx, hψinj, hwadjoin, hwminpoly, hxadjoin, hxminpoly⟩ :=
    ABC3.Found.Falt1.falt1BaseChangeGeneratorFull π n hn hπne0 hprime hnotsq hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn
  letI algV1Wn1 : Algebra V1 Wn1 := ψ.toRingHom.toAlgebra
  haveI hSTV0V1Wn1 : IsScalarTower V0 V1 Wn1 := by
    apply IsScalarTower.of_algebraMap_eq
    intro y
    exact (ψ.commutes y).symm
  haveI hSTV0WnWn1 : IsScalarTower V0 Wn Wn1 := by infer_instance
  haveI hFinV0V1 : Module.Finite V0 V1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0V1 π n hn hπne0 hprime hnotsq hnpos
  haveI hFinV0Wn1 : Module.Finite V0 Wn1 :=
    ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
  haveI hFinV1Wn1 : Module.Finite V1 Wn1 := by
    haveI := ABC3.Found.Falt1.falt1ModuleFiniteV0Wn1 π n hnpos hn' hπne0' hprime' hnotsq'
    exact Module.Finite.of_restrictScalars_finite V0 V1 Wn1
  haveI hIsDomV1 : IsDomain V1 := inferInstance
  haveI hIsDomWn1 : IsDomain Wn1 := inferInstance
  haveI hTFV1Wn1 : Module.IsTorsionFree V1 Wn1 :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hψinj
  haveI hTFV0Wn' : Module.IsTorsionFree V0 Wn :=
    ABC3.Found.Falt1.moduleIsTorsionFree_of_injective hinjV0Wn
  haveI hIsIntClosedV0 : IsIntegrallyClosed V0 := inferInstance
  letI algFRWn1 := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).1
  have hsep := (ABC3.Found.Falt1.falt1_hsep_bundled π n hnpos hn' hπne0' hprime' hnotsq' hinjV0Wn).2
  have hdiamond := ABC3.Found.Falt1.differentIdeal_tower_diamond
    (Vn := V0) (Vn1 := V1) (Wn := Wn) (Wn1 := Wn1) hsep
  have hwfield := ABC3.Found.Falt1.falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly π n hn hπne0 hprime hnotsq hnpos w hwminpoly
  have hdiffV1 := ABC3.Found.Falt1.differentIdeal_eq_span_derivative w hwfield hwadjoin
  rw [hwminpoly] at hdiffV1
  have hxfield := ABC3.Found.Falt1.falt1_fieldLevel_adjoin_top_of_ringLevel_minpoly (algebraMap V0 Wn π) n hn' hπne0' hprime' hnotsq' hnpos x hxminpoly
  have hcond := conductor_mul_differentIdeal Wn (FractionRing Wn) (AdjoinRoot gK) x hxfield
  rw [hxminpoly] at hcond
  have hderivV0 : Polynomial.aeval w (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0))
      = (n : V1) * w ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C π : Polynomial V0)
        = Polynomial.C (n : V0) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  have hderivWn : Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))
      = (n : Wn1) * x ^ (n-1) := by
    have hderiv : Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)
        = Polynomial.C (n : Wn) * Polynomial.X ^ (n-1) := by
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    simp [Polynomial.aeval_mul, Polynomial.aeval_X_pow]
  have hspan_eq : Ideal.span {Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))}
      = Ideal.map ψ.toRingHom (differentIdeal V0 V1) := by
    rw [hderivWn, hdiffV1, hderivV0, Ideal.map_span]
    congr 1
    simp only [Set.image_singleton]
    congr 1
    rw [map_mul, map_pow, map_natCast, show ψ.toRingHom w = ψ w from rfl, hwx]
  letI algFRV1 := (ABC3.Found.Falt1.falt1_hsepV0V1_bundled π n hn hπne0 hprime hnotsq hnpos).1
  have hsepV0V1 : @Algebra.IsSeparable (FractionRing V0) (FractionRing V1) _ _ algFRV1 :=
    (ABC3.Found.Falt1.falt1_hsepV0V1_bundled π n hn hπne0 hprime hnotsq hnpos).2
  have hdiffV1_ne : differentIdeal V0 V1 ≠ ⊥ :=
    @differentIdeal_ne_bot V0 V1 _ _ _ _ _ _ _ _ hsepV0V1
  have hIdiff_ne : Ideal.map ψ.toRingHom (differentIdeal V0 V1) ≠ 0 := by
    show Ideal.map ψ.toRingHom (differentIdeal V0 V1) ≠ ⊥
    intro heq
    exact hdiffV1_ne ((Ideal.map_eq_bot_iff_of_injective hψinj).mp heq)
  have hfinal := ABC3.Found.Falt1.cancel_conductor_delta
    (Ideal.map ψ.toRingHom (differentIdeal V0 V1)) (differentIdeal Wn Wn1)
    (conductor Wn x) (Ideal.span {Polynomial.aeval x (Polynomial.derivative (Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))})
    (differentIdeal V1 Wn1) (Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn))
    hdiamond hcond hspan_eq hIdiff_ne
  have hcondprod_ne : conductor Wn x * differentIdeal Wn Wn1 ≠ 0 := by
    rw [hcond, hspan_eq]
    exact hIdiff_ne
  have hcond_top : conductor Wn x = ⊤ := conductor_eq_top_iff_adjoin_eq_top.mpr hxadjoin
  have hIWnWn1_ne : differentIdeal Wn Wn1 ≠ 0 := by
    have h := hcondprod_ne
    rw [hcond_top, Ideal.top_mul] at h
    exact h
  haveI hPIDWn1 : IsPrincipalIdealRing Wn1 :=
    ABC3.Found.Falt1.falt1_isPrincipalIdealRing_of_finite_ext_of_DVR (Wn := Wn) (Wn1 := Wn1)
  have hresult : Module.length Wn1 (Wn1 ⧸ differentIdeal V0 Wn1) =
      Module.length Wn1 (Wn1 ⧸ differentIdeal Wn Wn1) +
      Module.length Wn1 (Wn1 ⧸ Ideal.map (algebraMap Wn Wn1) (differentIdeal V0 Wn)) :=
    ABC3.Found.Falt1.falt1_differentIdeal_tower_length hsep hIWnWn1_ne
  exact hresult

/-!
## Theorem 1.2: `Ω¹` の長さと `differentIdeal` の閉じた式を接続(2026-09-04、条件付き完成)

`falt1_theorem12_differentIdeal_length` の左辺 `differentIdeal V0 Wₙ₊₁`
を `Ω¹_{Wₙ₊₁/V0}` の長さと結びつけるには Lemma 1.1(`falt1CokernelLengthEq`)
を `V0→Wₙ₊₁` に適用する必要があるが、これには`Wₙ₊₁`の**`V0`上の単項
生成元**が要る——falt1-goal.md に記録した通り、これが一般に存在するか
(`Wₙ₊₁`が半局所環でも成り立つか)は未確認。ここでは、その生成元の
存在を**仮定として明示的に受け取る**条件付きの版を確立した——
`Wₙ₊₁`が`V0`上局所的(全分岐)であることは要求しない(`falt1CokernelLengthEq`
自体がそれを要求しないため)。 -/

set_option maxHeartbeats 3000000 in
set_option synthInstance.maxHeartbeats 3000000 in
/-- **`falt1CokernelLengthEq`(Lemma 1.1)を`V0→Wₙ₊₁`に適用した版**:
`Wₙ₊₁`の`V0`上の単項生成元`y`(環・体どちらのレベルでも生成する)が
与えられれば、`Ω¹_{Wₙ₊₁/V0}`の長さは`Wₙ₊₁⧸differentIdeal V0 Wₙ₊₁`の
長さに等しい。`IsIntegralClosure Wₙ₊₁ V0 (FractionRing Wₙ₊₁)`は
`falt1_theorem12_differentIdeal_length`の証明で使ったのと同じ
`IsIntegralClosure.of_isIntegrallyClosed`で無条件に得られる。 -/
theorem falt1_theorem12_kaehler_length
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    [Module.Finite V0 (Falt1Wn1 V0 Wn π n)]
    [Algebra (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [IsScalarTower V0 (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [FiniteDimensional (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [IsDedekindDomain (Falt1Wn1 V0 Wn π n)]
    [Module.IsTorsionFree V0 (Falt1Wn1 V0 Wn π n)]
    (y : Falt1Wn1 V0 Wn π n) (hint : IsIntegral V0 y)
    (hadjoin : Algebra.adjoin V0 ({y} : Set (Falt1Wn1 V0 Wn π n)) = ⊤)
    (hw : Algebra.adjoin (FractionRing V0)
      ({(algebraMap (Falt1Wn1 V0 Wn π n) (FractionRing (Falt1Wn1 V0 Wn π n))) y} :
        Set (FractionRing (Falt1Wn1 V0 Wn π n))) = ⊤) :
    Module.length (Falt1Wn1 V0 Wn π n) (Ω[(Falt1Wn1 V0 Wn π n)⁄V0]) =
      Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸ differentIdeal V0 (Falt1Wn1 V0 Wn π n)) := by
  haveI : Algebra.IsIntegral V0 (Falt1Wn1 V0 Wn π n) := Algebra.IsIntegral.of_finite V0 (Falt1Wn1 V0 Wn π n)
  haveI hIC : IsIntegralClosure (Falt1Wn1 V0 Wn π n) V0 (FractionRing (Falt1Wn1 V0 Wn π n)) :=
    IsIntegralClosure.of_isIntegrallyClosed (Falt1Wn1 V0 Wn π n) V0 (FractionRing (Falt1Wn1 V0 Wn π n))
  exact ABC3.Found.Falt1.falt1CokernelLengthEq y hint hadjoin hw

set_option maxHeartbeats 3000000 in
set_option synthInstance.maxHeartbeats 3000000 in
/-- **Theorem 1.2・3b/3c の中核関係式(`Ω¹` 版、条件付き完成)**:
`falt1_theorem12_kaehler_length`(Lemma 1.1経由)と`falt1_theorem12_
differentIdeal_length`(differentIdealの塔経由)を合わせるだけで、
`Ω¹_{Wₙ₊₁/V0}`の長さが**完全に`differentIdeal`だけで書かれた式**に
等しいことが分かる——`Wₙ₊₁`の`V0`上の単項生成元`y`さえ与えられれば、
Kähler微分の完全列(Exercise 13.7.4)と`differentIdeal`の塔公式という
2つの独立した経路がここで合流する。 -/
theorem falt1_theorem12_kaehler_length_eq_differentIdeal
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn : (n : FractionRing V0) ≠ 0) (hπne0 : algebraMap V0 (FractionRing V0) π ≠ 0)
    (hprime : (Ideal.span ({π} : Set V0)).IsPrime) (hnotsq : π ∉ (Ideal.span ({π} : Set V0)) ^ 2)
    (hnpos : 0 < n)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hinjV0Wn : Function.Injective (algebraMap V0 Wn))
    [Module.Finite V0 Wn]
    [Module.IsTorsionFree V0 Wn]
    [Algebra (FractionRing V0) (FractionRing Wn)]
    [IsScalarTower V0 (FractionRing V0) (FractionRing Wn)]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing Wn)]
    [IsDedekindDomain (Falt1Wn1 V0 Wn π n)]
    [Module.IsTorsionFree V0 (Falt1Wn1 V0 Wn π n)]
    [Module.IsTorsionFree Wn (Falt1Wn1 V0 Wn π n)]
    [Module.Finite V0 (Falt1Wn1 V0 Wn π n)]
    [Algebra (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [IsScalarTower V0 (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [FiniteDimensional (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    [Algebra.IsSeparable (FractionRing V0) (FractionRing (Falt1Wn1 V0 Wn π n))]
    (y : Falt1Wn1 V0 Wn π n) (hint : IsIntegral V0 y)
    (hadjoin : Algebra.adjoin V0 ({y} : Set (Falt1Wn1 V0 Wn π n)) = ⊤)
    (hw : Algebra.adjoin (FractionRing V0)
      ({(algebraMap (Falt1Wn1 V0 Wn π n) (FractionRing (Falt1Wn1 V0 Wn π n))) y} :
        Set (FractionRing (Falt1Wn1 V0 Wn π n))) = ⊤) :
    Module.length (Falt1Wn1 V0 Wn π n) (Ω[(Falt1Wn1 V0 Wn π n)⁄V0]) =
      Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸ differentIdeal Wn (Falt1Wn1 V0 Wn π n)) +
      Module.length (Falt1Wn1 V0 Wn π n) ((Falt1Wn1 V0 Wn π n) ⧸
        Ideal.map (algebraMap Wn (Falt1Wn1 V0 Wn π n)) (differentIdeal V0 Wn)) := by
  have h1 := ABC3.Found.Falt1.falt1_theorem12_kaehler_length π n y hint hadjoin hw
  have h2 := ABC3.Found.Falt1.falt1_theorem12_differentIdeal_length π n hn hπne0 hprime hnotsq hnpos
    hn' hπne0' hprime' hnotsq' hinjV0Wn
  rw [h1, h2]

/-!
## Exercise 13.7.4 (4)への足場: 第二基本完全列の単射性(2026-09-04)

Brinon-Conrad Exercise 13.7.4 の step (4) は
`0 → Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0} → Ω¹_{Wₙ₊₁/V0} → Ω¹_{Wₙ₊₁/Wₙ} → 0`
という完全列を要求する。この左側の**単射性**が、Lemma 1.1 の土台で
既に確立した `mapBaseChange_injective_of_nzd`(`B = Polynomial V/(f)`
かつ `f'` が非零因子の場合の単射性)を、`V0→Wₙ→AdjoinRoot g`
(`g := X^n-Cπ'`、base change **する前**の Eisenstein 多項式)に
適用するだけで得られることを確認した——`AdjoinRoot g`(`Wₙ₊₁` への
同型は `falt1AdjoinRootEquivIntegralClosure` 経由で既存)に対して
直接使える。次回はこれを`omegaCongr`(既存、代数同型に沿った `Ω`
の transport)で `Wₙ₊₁` 自身の単射性に持ち上げ、右側の全射性
(一般に成り立つ標準的事実)と合わせて完全列そのものを組み立てる
ところから続ける。 -/

set_option maxHeartbeats 800000 in
/-- **`Ω¹` の第二基本完全列・左側の単射性**(`AdjoinRoot g` 版、
`Wₙ₊₁` への transport は次回)。`g` が Eisenstein なので `g'` は
`AdjoinRoot g` の非零因子——`mapBaseChange_injective_of_nzd` を
直接適用するだけで閉じる。 -/
theorem falt1_mapBaseChange_injective_adjoinRoot
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n) :
    Function.Injective (KaehlerDifferential.mapBaseChange V0 Wn
      (AdjoinRoot ((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)))) := by
  set π' : Wn := algebraMap V0 Wn π with hπ'def
  set gpoly : Polynomial Wn := Polynomial.X ^ n - Polynomial.C π' with hgdef
  have hirr : Irreducible gpoly := ABC3.Found.Falt1.eisenstein_X_pow_sub_C π' n hnpos hprime' hnotsq'
  have hB : RingHom.ker (algebraMap (Polynomial Wn) (AdjoinRoot gpoly)) = Ideal.span ({gpoly} : Set (Polynomial Wn)) :=
    ABC3.Found.Falt1.ker_adjoinRoot_mk gpoly
  have hsurj : Function.Surjective (algebraMap (Polynomial Wn) (AdjoinRoot gpoly)) := AdjoinRoot.mk_surjective
  haveI : IsDomain (AdjoinRoot gpoly) := AdjoinRoot.isDomain_of_prime hirr.prime
  have hnWn : (n : Wn) ≠ 0 := by
    intro h
    apply hn'
    have : (algebraMap Wn (FractionRing Wn)) (n:Wn) = (n : FractionRing Wn) := by push_cast; ring
    rw [← this, h, map_zero]
  have hmonic : gpoly.Monic := Polynomial.monic_X_pow_sub_C π' hnpos.ne'
  have hdeg : gpoly.natDegree = n := by
    rw [hgdef, Polynomial.natDegree_sub_C]
    exact Polynomial.natDegree_X_pow n
  have hderivne : algebraMap (Polynomial Wn) (AdjoinRoot gpoly) (Polynomial.derivative gpoly) ≠ 0 := by
    have hderiv : Polynomial.derivative gpoly = Polynomial.C (n : Wn) * Polynomial.X ^ (n-1) := by
      rw [hgdef]
      rw [Polynomial.derivative_sub, Polynomial.derivative_X_pow, Polynomial.derivative_C, sub_zero]
    rw [hderiv]
    show (AdjoinRoot.mk gpoly) (Polynomial.C (n:Wn) * Polynomial.X ^ (n-1)) ≠ 0
    apply AdjoinRoot.mk_ne_zero_of_natDegree_lt hmonic
    · exact mul_ne_zero (Polynomial.C_ne_zero.mpr hnWn) (pow_ne_zero _ Polynomial.X_ne_zero)
    · rw [Polynomial.natDegree_mul (Polynomial.C_ne_zero.mpr hnWn) (pow_ne_zero _ Polynomial.X_ne_zero)]
      rw [Polynomial.natDegree_C, Polynomial.natDegree_X_pow, zero_add, hdeg]
      omega
  haveI : IsScalarTower V0 (Polynomial Wn) (AdjoinRoot gpoly) := by
    apply IsScalarTower.of_algebraMap_eq
    intro x
    rfl
  exact ABC3.Found.Falt1.mapBaseChange_injective_of_nzd gpoly hB hsurj
    (mem_nonZeroDivisors_of_ne_zero hderivne)

/-!
## Exercise 13.7.4 (4)への足場・続き: `Wₙ₊₁` 自身への輸送(2026-09-04)

`falt1_mapBaseChange_injective_adjoinRoot`(`AdjoinRoot g` 版)を、
既存の `mapBaseChange_injective_transport`(代数同型に沿った単射性の
transport、`falt1MapBaseChangeInjective` で既に使ったのと同じ道具)
と `falt1AdjoinRootEquivIntegralClosure`(`AdjoinRoot g ≃ₐ[Wₙ] Wₙ₊₁`)
を組み合わせるだけで `Wₙ₊₁` 自身の単射性が得られた——`omegaCongr` を
経由した自前の自然性証明は不要だった(`mapBaseChange_injective_transport`
の中で既に `omegaCongr` が使われている)。 -/

set_option maxHeartbeats 800000 in
/-- **第二基本完全列・左側の単射性、`Wₙ₊₁` 自身に対して(完成)**:
`e2 : AdjoinRoot g ≃ₐ[Wₙ] Wₙ₊₁` があれば、`falt1_mapBaseChange_injective_adjoinRoot`
の結果を `mapBaseChange_injective_transport` で `Wₙ₊₁` に輸送できる。
`Wₙ₊₁` は抽象的に(`e2` の存在だけを仮定して)扱う——具体的な
`integralClosure` 版は `falt1_mapBaseChange_injective_wn1_concrete`。 -/
theorem falt1_mapBaseChange_injective_wn1
    {V0 Wn Wn1 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    [CommRing Wn1] [Algebra Wn Wn1] [Algebra V0 Wn1] [IsScalarTower V0 Wn Wn1]
    (π : V0) (n : ℕ)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n)
    (e2 : AdjoinRoot ((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)) ≃ₐ[Wn] Wn1) :
    Function.Injective (KaehlerDifferential.mapBaseChange V0 Wn Wn1) := by
  have hA := ABC3.Found.Falt1.falt1_mapBaseChange_injective_adjoinRoot π n hn' hprime' hnotsq' hnpos
  haveI hSTA : IsScalarTower V0 Wn (AdjoinRoot ((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn))) := by
    apply IsScalarTower.of_algebraMap_eq
    intro x
    rfl
  exact ABC3.Found.Falt1.mapBaseChange_injective_transport e2 hA

set_option maxHeartbeats 800000 in
/-- `falt1_mapBaseChange_injective_wn1` の**具体版**: `Wₙ₊₁` を
`integralClosure Wn (AdjoinRoot(baseChange g))`(Falt1 の実際の構成)
として、`falt1AdjoinRootEquivIntegralClosure` で `e2` を構成してから
適用する。これで Exercise 13.7.4 step (4) の完全列の**左側**
(単射性)が Falt1 の具体的な `Wₙ₊₁` に対して完成した。 -/
theorem falt1_mapBaseChange_injective_wn1_concrete
    {V0 Wn : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    (π : V0) (n : ℕ)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hπne0' : algebraMap Wn (FractionRing Wn) (algebraMap V0 Wn π) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n)
    [IsDedekindDomain (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))]
    [Module.IsTorsionFree Wn (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))] :
    Function.Injective (KaehlerDifferential.mapBaseChange V0 Wn
      (integralClosure Wn (AdjoinRoot (((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)).map
        (algebraMap Wn (FractionRing Wn)))))) := by
  set e2 := ABC3.Found.Falt1.falt1AdjoinRootEquivIntegralClosure (algebraMap V0 Wn π) n hn' hπne0' hprime' hnotsq' hnpos
  exact falt1_mapBaseChange_injective_wn1 π n hn' hprime' hnotsq' hnpos e2

/-!
## Exercise 13.7.4 (4)完全列そのもの・長さの加法性(2026-09-04)

第二基本完全列 `0 → Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0} → Ω¹_{Wₙ₊₁/V0} → Ω¹_{Wₙ₊₁/Wₙ} → 0`
の**残り2つの成分**(中央の完全性・右側の全射性)は、実は
**mathlib に無条件で(=分岐や分離性の仮定なしに)既に存在する**と
判明した: `KaehlerDifferential.exact_mapBaseChange_map`(中央の完全性、
任意の塔 `R→A→B` で成立)・`KaehlerDifferential.map_surjective`
(`map R S B B : Ω[B⁄R]→Ω[B⁄S]` の全射性、任意の塔で成立)。左側の
単射性(`falt1_mapBaseChange_injective_wn1`、上で確立済み)と合わせ、
`Module.length_eq_add_of_exact` を適用するだけで**長さの加法公式**
(Brinon-Conrad step (4)-(6) の議論の土台)が閉じる。 -/

set_option maxHeartbeats 800000 in
/-- **`Ω¹` の第二基本完全列から従う長さの加法公式(完成)**:
`length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/V0}) = length_{Wₙ₊₁}(Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}) +
length_{Wₙ₊₁}(Ω¹_{Wₙ₊₁/Wₙ})`。Exercise 13.7.4 step (4)-(6) の
kernel/cokernel 評価(Nakayama+elementary divisors・discriminant の塔)
はこの等式の右辺の2項をそれぞれ評価することに相当する——次回はそちら
へ進む。 -/
theorem falt1_kaehler_length_exact_wn1
    {V0 Wn Wn1 : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    [CommRing Wn1] [Algebra Wn Wn1] [Algebra V0 Wn1] [IsScalarTower V0 Wn Wn1]
    (π : V0) (n : ℕ)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n)
    (e2 : AdjoinRoot ((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)) ≃ₐ[Wn] Wn1) :
    Module.length Wn1 (Ω[Wn1⁄V0]) =
      Module.length Wn1 (TensorProduct Wn Wn1 Ω[Wn⁄V0]) + Module.length Wn1 (Ω[Wn1⁄Wn]) := by
  have hinj := ABC3.Found.Falt1.falt1_mapBaseChange_injective_wn1 π n hn' hprime' hnotsq' hnpos e2
  have hsurj := KaehlerDifferential.map_surjective V0 Wn Wn1
  have hexact := KaehlerDifferential.exact_mapBaseChange_map V0 Wn Wn1
  exact Module.length_eq_add_of_exact (KaehlerDifferential.mapBaseChange V0 Wn Wn1)
    (KaehlerDifferential.map V0 Wn Wn1 Wn1) hinj hsurj hexact

set_option maxHeartbeats 800000 in
/-- **`falt1_kaehler_length_exact_wn1` の右辺第2項を Lemma 1.1
(`falt1CokernelLengthEq`)で `differentIdeal` の言葉に置き換えた版**:
`Ω¹_{Wₙ₊₁/Wₙ}` の長さは Lemma 1.1 により `Wₙ₊₁ ⧸ differentIdeal Wₙ Wₙ₊₁`
の長さに等しい(`x` が `Wₙ₊₁` を(環・体どちらのレベルでも)生成する
という Lemma 1.1 の前提を additional に要求する)。これで
`cancel_conductor_delta` 経由の既存の等式(`differentIdeal` の言葉で
書かれた `hlen_eq`)と、この Kähler 微分の完全列から来る等式とが
**同じ `differentIdeal Wₙ Wₙ₊₁` の長さを共有する**ことが分かり、
2つの証明経路が独立ではなく同じ量を経由して繋がっていることが
確認できた。残るは第1項(`Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}` の長さ、kernel 側
=Exercise 13.7.4 step (4)後半の Nakayama+elementary divisors・
step (5) の discriminant の塔に相当)の評価。 -/
theorem falt1_kaehler_length_exact_wn1_cokernel
    {V0 Wn Wn1 K L : Type*} [CommRing V0] [IsDomain V0] [IsDiscreteValuationRing V0]
    [CommRing Wn] [IsDomain Wn] [IsDiscreteValuationRing Wn] [Algebra V0 Wn]
    [CommRing Wn1] [Algebra Wn Wn1] [Algebra V0 Wn1] [IsScalarTower V0 Wn Wn1]
    [Field K] [Algebra Wn K] [IsFractionRing Wn K] [Field L] [Algebra K L]
    [FiniteDimensional K L] [Algebra.IsSeparable K L]
    [Algebra Wn1 L] [Algebra Wn L] [IsScalarTower Wn K L] [IsScalarTower Wn Wn1 L]
    [IsIntegralClosure Wn1 Wn L] [IsDedekindDomain Wn1] [Module.IsTorsionFree Wn Wn1]
    (π : V0) (n : ℕ)
    (hn' : (n : FractionRing Wn) ≠ 0)
    (hprime' : (Ideal.span ({algebraMap V0 Wn π} : Set Wn)).IsPrime)
    (hnotsq' : algebraMap V0 Wn π ∉ (Ideal.span ({algebraMap V0 Wn π} : Set Wn)) ^ 2)
    (hnpos : 0 < n)
    (e2 : AdjoinRoot ((Polynomial.X ^ n - Polynomial.C (algebraMap V0 Wn π) : Polynomial Wn)) ≃ₐ[Wn] Wn1)
    (x : Wn1) (hint : IsIntegral Wn x) (hadjoin : Algebra.adjoin Wn ({x} : Set Wn1) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap Wn1 L) x} : Set L) = ⊤) :
    Module.length Wn1 (Ω[Wn1⁄V0]) =
      Module.length Wn1 (TensorProduct Wn Wn1 Ω[Wn⁄V0]) +
        Module.length Wn1 (Wn1 ⧸ differentIdeal Wn Wn1) := by
  have hbase := ABC3.Found.Falt1.falt1_kaehler_length_exact_wn1 π n hn' hprime' hnotsq' hnpos e2
  rw [ABC3.Found.Falt1.falt1CokernelLengthEq x hint hadjoin hw] at hbase
  exact hbase

/-!
## 長さの加法公式・右辺第1項(kernel側)の評価(2026-09-04)

`Wₙ₊₁⊗_{Wₙ}Ω¹_{Wₙ/V0}`(kernel 側)は、`Wₙ ⊗_{Wₙ}` ではなく素朴な
「`Ω¹_{Wₙ/V0}` を `Wₙ₊₁` へ base change しただけ」の項なので、Lemma 1.1
(`falt1CokernelIsoLinear`、`Ω¹_{Wₙ/V0} ≃ₗ[Wₙ] Wₙ⧸differentIdeal V0 Wₙ`)を
`V0→Wₙ` にそのまま適用し、`LinearEquiv.baseChange` で `Wₙ₊₁` へ base
change するだけで評価できる——ただし `Wₙ₊₁⊗_{Wₙ}(Wₙ⧸I)` を
`Wₙ₊₁⧸I.map(algebraMap)` に単純化する一般補題(`Algebra.TensorProduct.
tensorQuotientEquiv` と `rid` の合成)がmathlibに直接無かったので
`falt1_tensorQuotientEquiv_algebraMap` として新規に用意した。 -/

set_option maxHeartbeats 1000000 in
/-- **`A ⊗[R] (R ⧸ I) ≃ₐ[A] A ⧸ I.map(algebraMap R A)`(一般補題、新規)**:
`Algebra.TensorProduct.tensorQuotientEquiv`(`A⊗[R](T⧸I)≃(A⊗[R]T)⧸...`)を
`T:=R`の場合に適用し、`Algebra.TensorProduct.rid`(`A⊗[R]R≃A`)で
右辺のテンソル積を消す。 -/
noncomputable def falt1_tensorQuotientEquiv_algebraMap
    {Wn Wn1 : Type*} [CommRing Wn] [CommRing Wn1] [Algebra Wn Wn1] (I : Ideal Wn) :
    TensorProduct Wn Wn1 (Wn ⧸ I) ≃ₐ[Wn1] Wn1 ⧸ (I.map (algebraMap Wn Wn1)) := by
  set e1 := Algebra.TensorProduct.tensorQuotientEquiv (R := Wn) Wn1 Wn Wn1 I with he1
  set e2 := Algebra.TensorProduct.rid Wn Wn1 Wn1 with he2
  set f1 : Wn →+* TensorProduct Wn Wn1 Wn :=
    (Algebra.TensorProduct.includeRight : Wn →ₐ[Wn] TensorProduct Wn Wn1 Wn).toRingHom with hf1
  set g1 : TensorProduct Wn Wn1 Wn →+* Wn1 := e2.toRingEquiv.toRingHom with hg1
  have hcomp : g1.comp f1 = algebraMap Wn Wn1 := by
    apply RingHom.ext
    intro w
    show (e2 : TensorProduct Wn Wn1 Wn →+* Wn1) ((Algebra.TensorProduct.includeRight : Wn →ₐ[Wn] TensorProduct Wn Wn1 Wn) w) = algebraMap Wn Wn1 w
    rw [Algebra.TensorProduct.includeRight_apply]
    show e2 ((1:Wn1) ⊗ₜ[Wn] w) = algebraMap Wn Wn1 w
    rw [he2, Algebra.TensorProduct.rid_tmul, ← Algebra.algebraMap_eq_smul_one]
  have hmapeq : (I.map (algebraMap Wn Wn1) : Ideal Wn1) = Ideal.map g1 (Ideal.map f1 I) := by
    rw [Ideal.map_map, hcomp]
  set e3 := Ideal.quotientEquivAlg (Ideal.map f1 I) (I.map (algebraMap Wn Wn1)) e2 hmapeq with he3
  exact e1.trans e3

set_option maxHeartbeats 1000000 in
/-- **長さの加法公式の右辺第1項(kernel側)を`differentIdeal V0 Wₙ`の
言葉に変換(完成)**: Lemma 1.1(`falt1CokernelIsoLinear`、`V0→Wₙ`)を
`Wₙ₊₁`へ base change し、`falt1_tensorQuotientEquiv_algebraMap`で
単純化するだけ。これで `falt1_kaehler_length_exact_wn1` の右辺2項が
どちらも`differentIdeal`の言葉に変換できた——次はこの結果と
`cancel_conductor_delta`経由の`hlen_eq`(`falt1_cancelConductorDelta_
assembled`)を組み合わせ、`Ω¹_{Wₙ₊₁/V0}`の長さを`differentIdeal V1 Wₙ₊₁`・
`differentIdeal Wₙ Wₙ₊₁`という2つの「小さい」different ideal の長さの
和として書き切ることに進む(Theorem 1.2 の再帰の核心になる可能性)。 -/
theorem falt1_kaehler_length_exact_wn1_kernel
    {V0 K L Wn Wn1 : Type*} [CommRing V0] [IsDedekindDomain V0]
    [Field K] [Algebra V0 K] [IsFractionRing V0 K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing Wn] [Algebra Wn L] [Algebra V0 Wn] [Algebra V0 L]
    [IsScalarTower V0 K L] [IsScalarTower V0 Wn L] [IsIntegralClosure Wn V0 L]
    [IsDedekindDomain Wn] [Module.IsTorsionFree V0 Wn]
    [CommRing Wn1] [Algebra Wn Wn1]
    (w : Wn) (hint : IsIntegral V0 w) (hadjoin : Algebra.adjoin V0 ({w} : Set Wn) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap Wn L) w} : Set L) = ⊤) :
    Module.length Wn1 (TensorProduct Wn Wn1 Ω[Wn⁄V0]) =
      Module.length Wn1 (Wn1 ⧸ (differentIdeal V0 Wn).map (algebraMap Wn Wn1)) := by
  set fbase := ABC3.Found.Falt1.falt1CokernelIsoLinear w hint hadjoin hw with hfbase
  set ftensor := LinearEquiv.baseChange Wn Wn1 Ω[Wn⁄V0] (Wn ⧸ differentIdeal V0 Wn) fbase with hftensor
  set fquot := falt1_tensorQuotientEquiv_algebraMap (Wn := Wn) (Wn1 := Wn1) (differentIdeal V0 Wn) with hfquot
  exact LinearEquiv.length_eq (ftensor.trans fquot.toLinearEquiv)

end ABC3.Found.Falt1
