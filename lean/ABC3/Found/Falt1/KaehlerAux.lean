import Mathlib.RingTheory.Kaehler.Basic
import Mathlib.RingTheory.Kaehler.Polynomial
import Mathlib.LinearAlgebra.Span.Basic
import Mathlib.RingTheory.AdjoinRoot

/-!
# [Falt1] Lemma 1.1 に向けた補助補題(`Found`、sorry 無し)

★★これは `/goal Falt1 Chapter I Found` の一部として、Lemma 1.1
(`Ω_V ⊗_V W → Ω_W` の単射性・余核の長さ)を実際に **証明** する試みの
第一歩。**Lemma 1.1 そのものはまだ証明できていない**——ここにあるのは
その途中で使う再利用可能な補題だけである(正直な進捗として記録する)。

`nzdCotangentEquivQuot`: 非零因子 `f` が生成する単項イデアル `I = (f)`
に対し `I.Cotangent ≃ₗ[A'] A' ⧸ I` が成り立つ(`I.Cotangent = I/I²` の
標準的な計算)。これは Falt1 Lemma 1.1 の証明の核——
`W = V[T]/(f(T))` のとき `Ω_{W/V} ≅ W/(f'(w))`(mathlib の
`KaehlerDifferential.exact_kerCotangentToTensor_mapBaseChange` と
組み合わせて使う)——の1コンポーネント。

**進捗**(2026-09-04):
1. ✅ `RingHom.ker (AdjoinRoot.mk f) = Ideal.span {f}`(`AdjoinRoot f`は
   定義から`Polynomial R ⧸ Ideal.span {f}`なので`Ideal.mk_ker`から従う)。
2. ✅ `kerQuotEquivOmega`: 全射`A→B`に対し
   `(B⊗_A Ω_{A/R}) ⧸ ker(mapBaseChange) ≃ₗ[B] Ω_{B/R}`
   (`LinearMap.quotKerEquivOfSurjective`と`mapBaseChange_surjective`から)。
3. ✅ `ker(mapBaseChange)`の台集合が`kerCotangentToTensor`の値域に一致する
   こと(`exact_kerCotangentToTensor_mapBaseChange`の集合レベルの言い換え)。

**確認できた証明の道筋**(2026-09-04、`lean_check`で各断片を個別に確認
済み——貼り合わせた1つの定理としてはまだ未完成):
- `TensorProduct.AlgebraTensorModule.rid` と `LinearEquiv.baseChange`
  を`KaehlerDifferential.polynomialEquiv`と組ませると
  `TensorProduct A B Ω[A⁄R] ≃ₗ[B] B`(A=Polynomial R の場合)が出る。
- `kerCotangentToTensor_toCotangent`(mathlib既存)より、生成元
  `[f] ∈ I.Cotangent` の像は `1 ⊗ₜ D(f)`、これは上の同一視の下で
  `algebraMap A B (Polynomial.derivative f) =: δ`(= Faltings の
  f'(w))に対応する。
- `I.Cotangent`(`I=(f)`)は`nzdCotangentEquivQuot`により`B`に同型で、
  生成元`[f]`は`1∈B`に対応する。`I`自身は`toSpanSingleton A A f`の
  値域として「Aから生成される巡回加群」であり、
  `range(kerCotangentToTensor) = A∙δ = B∙δ = Ideal.span{δ}`(A→Bが
  全射なのでA-span=B-spanは集合として一致、線形性の張り替えを経ずに
  Submodule.extで示せる)。
- これと`kerQuotEquivOmega`(既存)を合わせれば
  `Ω_{B/R} ≅ B ⧸ Ideal.span{δ}`(δ = f'(root))が得られる見通し。
  ★1つの定理として組み立てる作業はまだ残っている。

**残っている作業**(正直な記録):
4. 上の道筋を実際に1つの定理として組み立てる(`Ω_{B/R} ≅ B/(f'(root))`)。
5. **さらに大きな残作業**: 4 は "第二完全列"(Ω_{V[T]/V}を経由する
   商の計算)であって、Lemma 1.1 が実際に要求している "第一完全列"
   (`Ω_V ⊗_V W → Ω_W → Ω_{W/V} → 0`、V→W の塔に対するもの)とは別の
   完全列である。Lemma 1.1 の**単射性**の主張(第一完全列の最初の射が
   単射)は、この2つの完全列を貼り合わせる追加の議論(Faltings の原文
   だと「f'(w)dT の係数が非零因子だから」という具体的な計算)が必要で、
   まだ手を付けていない。
6. `different(W/V) = (f'(w))` という事実(mathlib の
   `aeval_derivative_mem_differentIdeal`・`conductor_mul_differentIdeal`
   と接続)、および「W/p^δW の長さ」という古典的な「長さ」概念を
   mathlib のどの道具(`Module.length`? 合成列?)で表現するかも未調査。

★見積り: Lemma 1.1 だけでもこの規模の補題を10-20個組み合わせる必要が
あり、数日〜数週間規模の作業になる可能性が高い。§2-4 は Faltings の
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

end ABC3.Found.Falt1
