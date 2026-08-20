import ABC3.Found.GaloisRep.RootFromIdeal
import Mathlib.RingTheory.ClassGroup.Basic

/-!
# Galois (G5) 第 140 ブロック —— **★★★★D2 の部品**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★D2(`J` が単項であること)の 2 つの部品

`WeilDivisor.lean` の D2 は「因子の類が 0 だから `J` は単項」である。
★その計算に要る 2 つを先に出す。

### ★★部品 1 —— `Σ_{T ∈ E[n]} T = 0`

`E[n] ≅ (ℤ/n)²` は (G1) で取得済である。★一般に

    B が有限可換群 ⟹ Σ_{p ∈ B × B} p = 0

が成り立つ(第 1 成分は `Σ_x Σ_y x = |B|·Σ_x x = 0`、第 2 成分も同様)。
★★`B × B` の形に同型な群なら、同型で送って同じ結論が出る。

★★★これが `Σ_{T ∈ E[n]} T = 0` を与える。

### ★★★部品 2 —— 類が自明なら単項

mathlib の `ClassGroup.mk_eq_one_iff` と `FractionalIdeal.isPrincipal_iff` を繋ぐ
(2026-08-20 実測、どちらも在庫あり)。

## ★★★★★D2 での使われ方

平行移動 `τ_T`(`T ∈ E[n]`)は `[n]∘τ_T = [n]` により `μ` を保つので、
分岐指数 `e_Q` は各ファイバー上で一定である。共通因子 `e` を括り出すと

    Σ_{Q ∈ [n]⁻¹(P)} Q = n²Q₀ = n·P = 0,      Σ_{T ∈ E[n]} T = 0

となり類は 0。★**`e` の値を知る必要はない。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sum_prod_self_eq_zero` | ★★★`Σ_{p ∈ B × B} p = 0` |
| `sum_univ_eq_zero_of_addEquiv` | ★★★★**`A ≃+ B × B` なら `Σ_{a ∈ A} a = 0`** |
| `isPrincipal_of_classGroup_eq_one` | ★★★類が自明な分数イデアルは単項 |
-/

namespace ABC3.Found.GaloisRep

open nonZeroDivisors

/-! ## ★★★部品 1 —— 有限可換群の直積の全元の和 -/

/-- ★★★**有限可換群の直積では全元の和が 0 になる**。

★第 1 成分は `Σ_x Σ_y x = |B|·Σ_x x = 0`(Lagrange)、第 2 成分も同様。 -/
theorem sum_prod_self_eq_zero {B : Type} [AddCommGroup B] [Fintype B] : ∑ p : B × B, p = 0 := by
  refine Prod.ext ?_ ?_
  · rw [Prod.fst_sum, Fintype.sum_prod_type]
    simp only [Finset.sum_const, Finset.card_univ]
    rw [← Finset.smul_sum]
    exact card_nsmul_eq_zero
  · rw [Prod.snd_sum, Fintype.sum_prod_type]
    simp only [Finset.sum_const, Finset.card_univ]
    exact card_nsmul_eq_zero

/-- ★★★★**`B × B` に同型な有限可換群では全元の和が 0 になる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`E[n] ≅ (ℤ/n)²`((G1) で取得済)に当てると `Σ_{T ∈ E[n]} T = 0` が出る。
★★これが D2(因子の類が自明であること)の計算に要る。 -/
theorem sum_univ_eq_zero_of_addEquiv {A B : Type} [AddCommGroup A] [Fintype A]
    [AddCommGroup B] [Fintype B] (e : A ≃+ B × B) : ∑ a : A, a = 0 := by
  have h : ∑ a : A, e a = 0 := by
    rw [Fintype.sum_equiv e.toEquiv (fun a => e a) (fun p => p) (fun _ => rfl)]
    exact sum_prod_self_eq_zero
  have he : e (∑ a : A, a) = 0 := by rw [map_sum]; exact h
  simpa using congrArg e.symm he

/-! ## ★★★部品 2 —— 類が自明なら単項 -/

/-- ★★★**類群で 1 になる分数イデアルは単項である**——mathlib の 2 つを繋いだだけ。 -/
theorem isPrincipal_of_classGroup_eq_one {R K : Type} [CommRing R] [IsDomain R] [Field K]
    [Algebra R K] [IsFractionRing R K] (J : (FractionalIdeal R⁰ K)ˣ)
    (h : ClassGroup.mk K J = 1) :
    ∃ g : K, (J : FractionalIdeal R⁰ K) = FractionalIdeal.spanSingleton R⁰ g :=
  (FractionalIdeal.isPrincipal_iff _).mp (ClassGroup.mk_eq_one_iff.mp h)

/-! ## ★出典の紐付け(`.src`) -/

def sum_univ_eq_zero_of_addEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——n 等分点の総和が 0 であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
