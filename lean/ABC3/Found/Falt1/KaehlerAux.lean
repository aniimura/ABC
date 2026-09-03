import Mathlib.RingTheory.Kaehler.Basic
import Mathlib.RingTheory.Kaehler.Polynomial
import Mathlib.LinearAlgebra.Span.Basic

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

**残っている作業**(2026-09-04 時点、正直な記録):
1. `RingHom.ker (algebraMap V[T] W) = Ideal.span {f}`(W = AdjoinRoot f
   等の具体形で)を示す。
2. `KaehlerDifferential.polynomialEquiv` と組み合わせて
   `Ω_{W/V} ≅ W ⧸ (f'(w))` を得る(モノジェニックな商の計算)。
3. **さらに大きな残作業**: これは "第二完全列"(Ω_{V[T]/V}を経由する
   商の計算)であって、Lemma 1.1 が実際に要求している "第一完全列"
   (`Ω_V ⊗_V W → Ω_W → Ω_{W/V} → 0`、V→W の塔に対するもの)とは別の
   完全列である。Lemma 1.1 の**単射性**の主張(第一完全列の最初の射が
   単射)は、この2つの完全列を貼り合わせる追加の議論(Faltings の原文
   だと「f'(w)dT の係数が非零因子だから」という具体的な計算)が必要で、
   まだ手を付けていない。
4. `different(W/V) = (f'(w))` という事実(mathlib の
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

end ABC3.Found.Falt1
