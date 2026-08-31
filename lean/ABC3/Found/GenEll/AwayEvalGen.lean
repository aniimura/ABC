/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayPolynomial
import ABC3.Meta.Claim

/-!
# ★★★★★★★★一般の 1 次形式による正規化 —— 段 C2c-1 の (a) の準備（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは何か —— なぜ一般化が要るのか

`§9-859` は `A⁰_{x_i}`（**変数**で割った形）について
`awayEval`・`ker_awayEval` を取った。★しかし段 C2c-1 の (a)（可換な四角）では
終域が `Away ℬ (hyperplaneHom (x_i))` であり、
**`hyperplaneHom (x_{i+1}) = y_i` は `rfl` ではない**（`§9-860` の測定）ので、
`Away ℬ f` の `f` が型に現れる以上、変数に特殊化した形では書けない。

★★本ファイルは `f` を**任意の 1 次斉次形式**に一般化する:

| 一般形 | 特殊形（`§9-859`） |
|---|---|
| `awayCoordOf R f hf j` | `projCoord N R i j`（★**定義的に等しい**） |
| `awayConstOf R f hf r` | `awayConst N R i r` |
| `awayEvalOf R f hf` | `awayEval N R i` |

★★★特殊形が一般形の**定義的な**場合であること（`awayCoordOf_X`）を測ってあるので、
`§9-859` の結果はそのまま一般形に持ち上がる。

## ★測定の記録

★`awayCoordOf R (X i) _ j = projCoord N R i j` は **`rfl`** である（2026-08-28 実測）
——`projCoord` が `Away.mk _ (isHomogeneous_X …) 1 (X j) _` だからである。
★★したがって一般化は**書き直しではなく被せ**である。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {σ : Type}

/-! ## ★一般の 1 次形式で割った座標と定数 -/

/-- ★**一般の 1 次形式 `f` による正規化座標** `x_j/f`。 -/
noncomputable def awayCoordOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (j : σ) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f :=
  HomogeneousLocalization.Away.mk _ hf 1 (MvPolynomial.X j)
    (by simpa using MvPolynomial.isHomogeneous_X R j)

/-- ★**一般の 1 次形式 `f` による定数** `C r / f^0`。 -/
noncomputable def awayConstOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (r : R) :
    HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f :=
  HomogeneousLocalization.Away.mk _ hf 0 (MvPolynomial.C r) (by simp)

/-! ## ★★特殊形は一般形の定義的な場合である -/

/-- ★★**`projCoord` は `awayCoordOf` の `f = x_i` の場合である**（`rfl`）。 -/
theorem awayCoordOf_X (N : ℕ) (R : Type) [CommRing R] (i j : Fin (N+1)) :
    awayCoordOf R (MvPolynomial.X i)
      (by simpa using MvPolynomial.isHomogeneous_X R i) j = projCoord N R i j := rfl

/-- ★★**`awayConst` は `awayConstOf` の `f = x_i` の場合である**（`rfl`）。 -/
theorem awayConstOf_X (N : ℕ) (R : Type) [CommRing R] [Nontrivial R] (i : Fin (N+1)) (r : R) :
    awayConstOf R (MvPolynomial.X i)
      (by simpa using MvPolynomial.isHomogeneous_X R i) r = awayConst N R i r := rfl

/-! ## ★★★定数は環準同型である -/

/-- ★★★**一般の `f` についても定数は環準同型である**。

★証明は `§9-859` の `awayConstHom` と同型である——`val_injective` で
`Localization.mk` の計算に落とすだけ。 -/
noncomputable def awayConstHomOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    R →+* HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f where
  toFun := awayConstOf R f hf
  map_one' := by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConstOf, Away.val_mk, HomogeneousLocalization.val_one]
    show Localization.mk (MvPolynomial.C (1 : R)) _ = 1
    rw [map_one]
    rw [show (⟨f ^ 0, ⟨0, rfl⟩⟩ : Submonoid.powers f) = 1 from by ext; simp]
    exact Localization.mk_one
  map_zero' := by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConstOf, Away.val_mk, HomogeneousLocalization.val_zero]
    show Localization.mk (MvPolynomial.C (0 : R)) _ = 0
    rw [map_zero, Localization.mk_zero]
  map_mul' := fun a b => by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConstOf, Away.val_mk, HomogeneousLocalization.val_mul, awayConstOf, awayConstOf,
      Away.val_mk, Away.val_mk, Localization.mk_mul, map_mul]
    congr 1
    ext
    simp
  map_add' := fun a b => by
    refine HomogeneousLocalization.val_injective _ ?_
    rw [awayConstOf, Away.val_mk, HomogeneousLocalization.val_add, awayConstOf, awayConstOf,
      Away.val_mk, Away.val_mk, map_add, Localization.add_mk_self]

/-! ## ★★★★★一般の `f` による評価 -/

/-- ★★★★★**一般の `f` による評価** `x_j ↦ x_j/f`。 -/
noncomputable def awayEvalOf (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    MvPolynomial σ R →+*
      HomogeneousLocalization.Away (MvPolynomial.homogeneousSubmodule σ R) f :=
  MvPolynomial.eval₂Hom (awayConstHomOf R f hf) (fun j => awayCoordOf R f hf j)

theorem awayEvalOf_C (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (r : R) :
    awayEvalOf R f hf (MvPolynomial.C r) = awayConstOf R f hf r := by
  rw [awayEvalOf, MvPolynomial.eval₂Hom_C]
  rfl

theorem awayEvalOf_X (R : Type) [CommRing R] (f : MvPolynomial σ R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (j : σ) :
    awayEvalOf R f hf (MvPolynomial.X j) = awayCoordOf R f hf j := by
  rw [awayEvalOf, MvPolynomial.eval₂Hom_X']

/-- ★★★**特殊形との一致**（`f = x_i`）。 -/
theorem awayEvalOf_eq_awayEval (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N+1)) :
    awayEvalOf R (MvPolynomial.X i)
      (by simpa using MvPolynomial.isHomogeneous_X R i) = awayEval N R i := by
  refine MvPolynomial.ringHom_ext (fun r => ?_) (fun j => ?_)
  · rw [awayEvalOf_C]
    show awayConstOf R (MvPolynomial.X i) _ r = awayEval N R i (MvPolynomial.C r)
    rw [awayEval, MvPolynomial.eval₂Hom_C]
    rfl
  · rw [awayEvalOf_X]
    show awayCoordOf R (MvPolynomial.X i) _ j = awayEval N R i (MvPolynomial.X j)
    rw [awayEval, MvPolynomial.eval₂Hom_X']
    rfl

/-! ## ★出典の紐付け(`.src`) -/

def awayCoordOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(一般の 1 次形式による正規化座標)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(一般の 1 次形式による評価)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf_eq_awayEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(一般形と特殊形の一致)",
    sectionId := "genell-prop-1-4" }

def awayEvalOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayEval / ker_awayEval(段 C2a-2、§9-859)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_awayEval") 2,
    .implicitStep
      ("★なぜ一般化が要るか: 段 C2c-1 の (a)(可換な四角)では終域が " ++
       "Away ℬ (hyperplaneHom (x_i)) であり、**hyperplaneHom (x_{i+1}) = y_i は rfl ではない**" ++
       "(§9-860 の測定)。Away ℬ f の f は**型に現れる**ので、" ++
       "変数に特殊化した形では書けない") 4,
    .implicitStep
      ("★★awayCoordOf R (X i) _ j = projCoord N R i j は **rfl** である(2026-08-28 実測)" ++
       "——projCoord が Away.mk _ (isHomogeneous_X …) 1 (X j) _ だからである。" ++
       "★したがって一般化は**書き直しではなく被せ**である") 3,
    .implicitStep
      ("★★★残るのは awayEvalOf_mk(斉次式の評価は a/f^k である)を一般形で取り、" ++
       "そこから可換な四角 Away.map g f ∘ awayEvalOf f = awayEvalOf (g f) ∘ g を出す段である") 5 ]

end ABC3.Found.GenEll
