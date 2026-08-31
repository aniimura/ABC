/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayEvalGen
import ABC3.Found.GenEll.AwayGenerated
import ABC3.Meta.Claim

/-!
# ★★★★★★★★一般の 1 次形式でも「斉次式の評価は `a/f^k`」（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは何か

`§9-861b`（`AwayEvalGen`）で `awayEvalOf`（**任意の 1 次斉次形式 `f`** による評価）を
定義した。★本ファイルはその**中心的な性質**を取る:

    `awayEvalOf f (a) = a/f^k`   （`a` は次数 `k` の斉次式）

★★これは `§9-859` の `awayEval_mk` の一般化であり、証明は `§9-850` の単項式分解
（`away_mk_monomial` ＋ `awayMkAddHom`）を `x_i` から `f` へ**そのまま被せた**ものである。

★★★これがあると

* `awayEvalOf` は**全射**である
* 段 C2c-1 の (a)（可換な四角 `Away.map g f ∘ awayEvalOf f = awayEvalOf (g f) ∘ g`）

がどちらも出る。

## ★測定の記録

★`f` を一般にしても証明は**一字も変わらない**（`isHomogeneous_X R i` を `hf` に
置き換えるだけ）——`x_i` であることを使っていたのは**分母が単項式であること**ではなく
**次数が 1 であること**だけだったからである（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {σ : Type}

/-! ## ★座標の冪と積の `val`（一般形） -/

theorem val_pow_awayCoordOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (j : σ) (k : ℕ) :
    ((awayCoordOf R f hf j) ^ k).val
      = Localization.mk ((MvPolynomial.X j) ^ k)
          (⟨f ^ k, ⟨k, rfl⟩⟩ : Submonoid.powers f) := by
  rw [HomogeneousLocalization.val_pow, awayCoordOf, Away.val_mk, Localization.mk_pow]
  congr 1
  ext
  simp

open scoped Classical in
theorem val_prod_awayCoordOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1)
    (S : Finset σ) (m : σ → ℕ) :
    (∏ j ∈ S, (awayCoordOf R f hf j) ^ (m j)).val
      = Localization.mk (∏ j ∈ S, (MvPolynomial.X j : MvPolynomial σ R) ^ (m j))
          (⟨f ^ (∑ j ∈ S, m j), ⟨_, rfl⟩⟩ : Submonoid.powers f) := by
  induction S using Finset.induction with
  | empty => simp
  | insert a S ha ih =>
      rw [Finset.prod_insert ha, Finset.prod_insert ha, Finset.sum_insert ha,
        HomogeneousLocalization.val_mul, val_pow_awayCoordOf, ih, Localization.mk_mul]
      congr 1
      ext
      simp [pow_add]

/-! ## ★★★★★★★★単項式は定数と座標の積である（一般形） -/

open scoped Classical in
/-- ★★★★★★★★**単項式は定数と正規化座標の積である**（任意の 1 次形式 `f`）。

    `(c·∏ x_j^{m_j})/f^k = c · ∏ (x_j/f)^{m_j}`   （`∑ m_j = k`） -/
theorem away_mk_monomial_of (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1)
    (m : σ →₀ ℕ) (c : R) (k : ℕ) (hk : (∑ j ∈ m.support, m j) = k)
    (hm : (MvPolynomial.monomial m c : MvPolynomial σ R)
      ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule σ R) hf k (MvPolynomial.monomial m c) hm
      = awayConstOf R f hf c * ∏ j ∈ m.support, (awayCoordOf R f hf j) ^ (m j) := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_mul, awayConstOf, Away.val_mk,
    val_prod_awayCoordOf, Localization.mk_mul]
  subst hk
  rw [MvPolynomial.monomial_eq]
  congr 1
  ext
  simp

/-! ## ★★★分子についての加法性（一般形） -/

theorem away_mk_add_of {R : Type} [CommRing R] {f : MvPolynomial σ R}
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (k : ℕ)
    (a b : MvPolynomial σ R)
    (ha : a ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1))
    (hb : b ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1))
    (hab : a + b ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule σ R) hf k (a + b) hab
      = Away.mk _ hf k a ha + Away.mk _ hf k b hb := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_add, Away.val_mk, Away.val_mk,
    Localization.add_mk_self]

theorem away_mk_zero_of {R : Type} [CommRing R] {f : MvPolynomial σ R}
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (k : ℕ)
    (h0 : (0 : MvPolynomial σ R) ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1)) :
    Away.mk (MvPolynomial.homogeneousSubmodule σ R) hf k 0 h0 = 0 := by
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_zero, Localization.mk_zero]

/-- ★★★**次数 `k` の斉次成分から `A⁰_f` への加法準同型**（一般形）。 -/
noncomputable def awayMkAddHomOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (k : ℕ) :
    (MvPolynomial.homogeneousSubmodule σ R (k • 1))
      →+ HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f where
  toFun := fun a => Away.mk _ hf k a.1 a.2
  map_zero' := away_mk_zero_of hf k _
  map_add' := fun a b => away_mk_add_of hf k a.1 b.1 a.2 b.2 (a + b).2

/-! ## ★★★★★★★★斉次式の評価はちょうど `a/f^k` である -/

/-- ★★★★★★★★**斉次式の評価はちょうど `a/f^k` である**（任意の 1 次形式 `f`）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-859` の `awayEval_mk` の一般化。★★証明は一字も変わらない
——`x_i` であることを使っていたのは**分母が単項式であること**ではなく
**次数が 1 であること**だけだったからである。 -/
theorem awayEvalOf_mk (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1)
    (k : ℕ) (a : MvPolynomial σ R)
    (ha : a ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1)) :
    awayEvalOf R f hf a
      = Away.mk (MvPolynomial.homogeneousSubmodule σ R) hf k a ha := by
  classical
  have hdeg : ∀ m ∈ a.support, (∑ j ∈ m.support, m j) = k := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule] at ha
    have h := ha (MvPolynomial.mem_support_iff.1 hm)
    rw [Finsupp.weight_apply] at h
    simpa [Finsupp.sum] using h
  have hmono : ∀ m ∈ a.support,
      (MvPolynomial.monomial m (MvPolynomial.coeff m a) : MvPolynomial σ R)
        ∈ MvPolynomial.homogeneousSubmodule σ R (k • 1) := by
    intro m hm
    rw [MvPolynomial.mem_homogeneousSubmodule]
    have h := MvPolynomial.isHomogeneous_monomial (d := m) (MvPolynomial.coeff m a)
      (n := k) (hdeg m hm)
    simpa using h
  have h1 : (⟨a, ha⟩ : MvPolynomial.homogeneousSubmodule σ R (k • 1))
      = ∑ m ∈ a.support.attach,
        ⟨MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a), hmono m.1 m.2⟩ := by
    apply Subtype.ext
    rw [Submodule.coe_sum]
    show a = ∑ m ∈ a.support.attach, MvPolynomial.monomial m.1 (MvPolynomial.coeff m.1 a)
    rw [Finset.sum_attach a.support (fun m => MvPolynomial.monomial m (MvPolynomial.coeff m a))]
    exact MvPolynomial.as_sum a
  have hR : Away.mk (MvPolynomial.homogeneousSubmodule σ R) hf k a ha
      = ∑ m ∈ a.support.attach, (awayConstOf R f hf (MvPolynomial.coeff m.1 a)
          * ∏ j ∈ (m.1).support, (awayCoordOf R f hf j) ^ (m.1 j)) := by
    show awayMkAddHomOf R f hf k ⟨a, ha⟩ = _
    rw [h1, map_sum]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    exact away_mk_monomial_of R f hf m.1 (MvPolynomial.coeff m.1 a) k (hdeg m.1 m.2)
      (hmono m.1 m.2)
  rw [hR]
  conv_lhs => rw [MvPolynomial.as_sum a,
    ← Finset.sum_attach a.support (fun m => MvPolynomial.monomial m (MvPolynomial.coeff m a))]
  rw [map_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  rw [awayEvalOf, MvPolynomial.eval₂Hom_monomial]
  rfl

/-- ★★**`awayEvalOf` は全射である**。 -/
theorem awayEvalOf_surjective (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    Function.Surjective (awayEvalOf R f hf) := by
  intro z
  obtain ⟨k, a, ha, rfl⟩ := Away.mk_surjective
    (MvPolynomial.homogeneousSubmodule σ R) hf z
  exact ⟨a, awayEvalOf_mk R f hf k a ha⟩

/-! ## ★出典の紐付け(`.src`) -/

def away_mk_monomial_of.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(単項式は定数と正規化座標の積である——一般の 1 次形式)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf_mk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(斉次式の評価は a/f^k である——一般の 1 次形式)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(awayEvalOf は全射である)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf_mk.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayEval_mk(x_i の場合、§9-859)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayEval_mk") 2,
    .citation "[ABC3]" "away_mk_monomial(単項式分解、§9-850)"
      (.inProject "ABC3" "ABC3.Found.GenEll.away_mk_monomial") 2,
    .implicitStep
      ("★f を一般にしても証明は**一字も変わらない**(isHomogeneous_X R i を hf に" ++
       "置き換えるだけ)——x_i であることを使っていたのは**分母が単項式であること**ではなく" ++
       "**次数が 1 であること**だけだったからである(2026-08-28 実測)") 3,
    .implicitStep
      ("★★これで段 C2c-1 の (a)(可換な四角 " ++
       "Away.map g f ∘ awayEvalOf f = awayEvalOf (g f) ∘ g)が書ける" ++
       "——両辺は環準同型なので C と X で見ればよく、" ++
       "X の側は awayEvalOf_mk(k = 1)で潰れる") 4 ]

end ABC3.Found.GenEll
