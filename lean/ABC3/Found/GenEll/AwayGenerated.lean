/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ChartSurjectiveCoords
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★`A⁰_{x_i}` は定数と正規化座標で生成される（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6-7。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★これは何か —— `§9-849` が残した最後の穴

`§9-849` は「座標が点を分ける」を `A⁰_{x_i}` **全体**での一致から導いた。
★しかし `northcott_of_projModel` が渡すのは**正規化座標 `x_j/x_i` での一致**だけである。

★★本ファイルはその差を埋める:

> `A⁰_{x_i}` からの環準同型は、**定数 `C c/x_i^0` と正規化座標 `x_j/x_i` での値**で決まる。

★★★これは `Proj` の標準チャート `D₊(x_i) ≅ 𝔸ᴺ` という古典的事実の、
**消費側が要る形**である（mathlib に見当たらない、2026-08-28 実測）。

## ★★★★機構 —— 単項式に分けるだけ

`A⁰_{x_i}` の元は `a/x_i^k`（`a` は `k` 次斉次）である（`Away.mk_surjective`）。
★`a` を単項式に分けると（`MvPolynomial.as_sum`）、各単項式 `c·∏ x_j^{m_j}` は
`∑ m_j = k` を満たすので

    `(c·∏ x_j^{m_j})/x_i^k = c · ∏ (x_j/x_i)^{m_j}`

★★あとは分子についての加法性（`awayMkAddHom`）で足し合わせるだけである。

## ★測定の記録

★★★等式はすべて `HomogeneousLocalization.val_injective` で
`Localization.mk` の計算に落とす。★`Localization.mk_mul` / `mk_pow` / `add_mk_self` の 3 本で足りる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★座標の冪と積の `val` -/

/-- ★**正規化座標の冪の `val`**。 -/
theorem val_pow_projCoord (N : ℕ) (R : Type) [CommRing R] (i j : Fin (N + 1)) (k : ℕ) :
    ((projCoord N R i j) ^ k).val
      = Localization.mk ((MvPolynomial.X j) ^ k)
          (⟨(MvPolynomial.X i : MvPolynomial (Fin (N+1)) R) ^ k,
            ⟨k, rfl⟩⟩ : Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N+1)) R)) := by
  rw [HomogeneousLocalization.val_pow, projCoord, Away.val_mk, Localization.mk_pow]
  congr 1
  ext
  simp

open scoped Classical in
/-- ★★**正規化座標の積の `val`** —— 分母の指数はちょうど次数の和になる。 -/
theorem val_prod_projCoord (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1))
    (S : Finset (Fin (N + 1))) (m : Fin (N + 1) → ℕ) :
    (∏ j ∈ S, (projCoord N R i j) ^ (m j)).val
      = Localization.mk (∏ j ∈ S, (MvPolynomial.X j : MvPolynomial (Fin (N+1)) R) ^ (m j))
          (⟨(MvPolynomial.X i : MvPolynomial (Fin (N+1)) R) ^ (∑ j ∈ S, m j),
            ⟨_, rfl⟩⟩ : Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N+1)) R)) := by
  induction S using Finset.induction with
  | empty => simp
  | insert a S ha ih =>
      rw [Finset.prod_insert ha, Finset.prod_insert ha, Finset.sum_insert ha,
        HomogeneousLocalization.val_mul, val_pow_projCoord, ih, Localization.mk_mul]
      congr 1
      ext
      simp [pow_add]

/-! ## ★★★★★★★★単項式は定数と座標の積である -/

open scoped Classical in
/-- ★★★★★★★★**単項式は定数と正規化座標の積である**。

    `(c·∏ x_j^{m_j})/x_i^k = c · ∏ (x_j/x_i)^{m_j}`   （`∑ m_j = k`）

★これが「`A⁰_{x_i}` は定数と座標で生成される」ことの中身である。 -/
theorem away_mk_monomial (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N + 1))
    (m : Fin (N + 1) →₀ ℕ) (c : R) (k : ℕ) (hk : (∑ j ∈ m.support, m j) = k)
    (hm : (MvPolynomial.monomial m c : MvPolynomial (Fin (N+1)) R)
      ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
        (MvPolynomial.isHomogeneous_X R i) k (MvPolynomial.monomial m c) hm
      = awayConst N R i c * ∏ j ∈ m.support, (projCoord N R i j) ^ (m j) := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_mul, awayConst, Away.val_mk,
    val_prod_projCoord, Localization.mk_mul]
  subst hk
  rw [MvPolynomial.monomial_eq]
  congr 1
  ext
  simp

/-! ## ★★★分子についての加法性 -/

theorem away_mk_add {N : ℕ} {R : Type} [CommRing R] (i : Fin (N + 1)) (k : ℕ)
    (a b : MvPolynomial (Fin (N+1)) R)
    (ha : a ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1))
    (hb : b ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1))
    (hab : a + b ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
        (MvPolynomial.isHomogeneous_X R i) k (a + b) hab
      = Away.mk _ (MvPolynomial.isHomogeneous_X R i) k a ha
        + Away.mk _ (MvPolynomial.isHomogeneous_X R i) k b hb := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_add, Away.val_mk, Away.val_mk,
    Localization.add_mk_self]

theorem away_mk_zero {N : ℕ} {R : Type} [CommRing R] (i : Fin (N + 1)) (k : ℕ)
    (h0 : (0 : MvPolynomial (Fin (N+1)) R)
      ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
        (MvPolynomial.isHomogeneous_X R i) k 0 h0 = 0 := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_zero, Localization.mk_zero]

/-- ★★★**次数 `k` の斉次成分から `A⁰_{x_i}` への加法準同型**。

★これがあると `MvPolynomial.as_sum` の分解を `map_sum` で流し込める。 -/
noncomputable def awayMkAddHom {N : ℕ} (R : Type) [CommRing R] (i : Fin (N + 1)) (k : ℕ) :
    (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1))
      →+ HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i) where
  toFun := fun a => Away.mk _ (MvPolynomial.isHomogeneous_X R i) k a.1 a.2
  map_zero' := away_mk_zero i k _
  map_add' := fun a b => away_mk_add i k a.1 b.1 a.2 b.2 (a + b).2

/-! ## ★★★★★★★★★生成 —— 環準同型は定数と座標で決まる -/

open scoped Classical in
/-- ★★★★★★★★★**`A⁰_{x_i}` からの環準同型は定数と正規化座標での値で決まる**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★これが `§9-849`（座標が点を分ける）を `northcott_of_projModel` が渡す形
——**正規化座標 `x_j/x_i` での一致**——に繋ぐ最後の 1 本である。

★★機構は「元を `a/x_i^k` と書き（`Away.mk_surjective`）、`a` を単項式に分け
（`MvPolynomial.as_sum`）、各単項式に `away_mk_monomial` を当てる」だけである。 -/
theorem ext_of_projCoord {N : ℕ} {R : Type} [CommRing R] [Nontrivial R] (i : Fin (N + 1))
    {S : Type} [CommRing S]
    (ψ χ : HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.X i) →+* S)
    (hc : ∀ c : R, ψ (awayConst N R i c) = χ (awayConst N R i c))
    (hx : ∀ j, ψ (projCoord N R i j) = χ (projCoord N R i j)) :
    ψ = χ := by
  refine RingHom.ext (fun z => ?_)
  obtain ⟨k, a, ha, rfl⟩ := Away.mk_surjective
    (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R) (MvPolynomial.isHomogeneous_X R i) z
  have hdeg : ∀ m ∈ a.support, (∑ j ∈ m.support, m j) = k := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule] at ha
    have h := ha (MvPolynomial.mem_support_iff.1 hm)
    rw [Finsupp.weight_apply] at h
    simpa [Finsupp.sum] using h
  have hmono : ∀ m ∈ a.support,
      (MvPolynomial.monomial m (MvPolynomial.coeff m a) : MvPolynomial (Fin (N+1)) R)
        ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1) := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule]
    have h := MvPolynomial.isHomogeneous_monomial (d := m) (MvPolynomial.coeff m a)
      (n := k) (hdeg m hm)
    simpa using h
  have h1 : (⟨a, ha⟩ : MvPolynomial.homogeneousSubmodule (Fin (N+1)) R (k • 1))
      = ∑ m ∈ a.support.attach,
        ⟨MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a), hmono m.1 m.2⟩ := by
    apply Subtype.ext
    rw [Submodule.coe_sum]
    show a = ∑ m ∈ a.support.attach, MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a)
    rw [Finset.sum_attach a.support (fun m => MvPolynomial.monomial m (MvPolynomial.coeff m a))]
    exact MvPolynomial.as_sum a
  have hdec : Away.mk (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)
        (MvPolynomial.isHomogeneous_X R i) k a ha
      = ∑ m ∈ a.support.attach, (awayConst N R i (MvPolynomial.coeff m.1 a)
          * ∏ j ∈ (m.1).support, (projCoord N R i j) ^ (m.1 j)) := by
    show awayMkAddHom R i k ⟨a, ha⟩ = _
    rw [h1, map_sum]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    exact away_mk_monomial N R i m.1 (MvPolynomial.coeff m.1 a) k (hdeg m.1 m.2) (hmono m.1 m.2)
  rw [hdec, map_sum, map_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [map_mul, map_mul, hc, map_prod, map_prod]
  congr 1
  exact Finset.prod_congr rfl (fun j _ => by rw [map_pow, map_pow, hx])

/-! ## ★出典の紐付け(`.src`) -/

def away_mk_monomial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(単項式は定数と正規化座標の積である)",
    sectionId := "genell-prop-1-4" }

def awayMkAddHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数 k の斉次成分から A⁰_{x_i} への加法準同型)",
    sectionId := "genell-prop-1-4" }

def ext_of_projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(A⁰_{x_i} からの環準同型は定数と正規化座標で決まる)",
    sectionId := "genell-prop-1-4" }

def ext_of_projCoord.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "HomogeneousLocalization.Away.mk_surjective / val_injective"
      (.inMathlib "HomogeneousLocalization.Away.mk_surjective") 2,
    .citation "[mathlib]" "MvPolynomial.as_sum / isHomogeneous_monomial"
      (.inMathlib "MvPolynomial.as_sum") 2,
    .citation "[ABC3]" "ext_of_chart(座標が点を分ける、§9-849)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ext_of_chart") 2,
    .implicitStep
      ("★これは Proj の標準チャート D₊(x_i) ≅ 𝔸ᴺ という古典的事実の" ++
       "**消費側が要る形**である(mathlib に見当たらない、2026-08-28 実測)。" ++
       "★★等式はすべて val_injective で Localization.mk の計算に落とす" ++
       "——Localization.mk_mul / mk_pow / add_mk_self の 3 本で足りる") 4 ]

end ABC3.Found.GenEll
