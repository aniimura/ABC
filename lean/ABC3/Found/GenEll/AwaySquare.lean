/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwayEvalGenMk
import ABC3.Found.GenEll.Hyperplane
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★可換な四角 —— 段 C2c-1 の (a)（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★これは何か —— 段 C2c-1 の (a)

`§9-860` で段 C2c-1 の逆向きを 3 段に分けた:

| 段 | 内容 | 状態 |
|---|---|---|
| **(a)** | **可換な四角 `Φ ∘ Ψ = Ψ′ ∘ σ`** | ★★**本ファイル** |
| (b) | `ker σ = (x₀)` | ★済（`§9-861`、`KillX0`） |
| (c) | (a)(b) から `ker Φ = (x₀/x_i)` | ★開 |

★本ファイルが取るのは

```
  R[x] --awayEvalOf f--> A⁰_f
   |                      |
   g                  Away.map g f
   v                      v
  R[y] --awayEvalOf (g f)--> B⁰_{g f}
```

が**任意の次数付き環準同型 `g` と任意の 1 次形式 `f`** について可換である、ということである。

## ★★★機構 —— 生成元で見るだけ

★両辺は環準同型なので `MvPolynomial.ringHom_ext` で `C r` と `X j` だけ見ればよい。

| 生成元 | 左辺 | 右辺 |
|---|---|---|
| `C r` | `Away.map g f (C r / f⁰) = g(C r) / (g f)⁰` | `awayEvalOf (g f) (g (C r))` |
| `X j` | `Away.map g f (x_j / f) = g(x_j) / (g f)` | `awayEvalOf (g f) (g (X j))` |

★★右辺を潰すのが `awayEvalOf_mk`（`§9-861c`）である
——`g (C r)` は次数 `0`、`g (X j)` は次数 `1` の斉次式だから。

## ★★★★★測定の記録 —— `g (C r) = C r` は要らなかった

★最初は「`g` が定数を定数に送る」という仮定が要ると読んでいたが、**要らない**。
★★`awayEvalOf_mk` が**次数だけ**で右辺を `Away.mk` に潰すので、
`g (C r)` が何であっても（次数 `0` でありさえすれば）両辺は同じ `Away.mk` になる
（2026-08-28 実測）。
★★★これは「一般化した方が証明が短くなる」場合である。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

variable {σ τ : Type}

/-! ## ★生成元での可換性 -/

/-- ★**定数の側** —— 次数 `0` なので `awayEvalOf_mk` が潰す。 -/
theorem awayMap_awayConstOf (R : Type) [CommRing R]
    (g : (MvPolynomial.homogeneousSubmodule σ R) →+*ᵍ (MvPolynomial.homogeneousSubmodule τ R))
    (f : MvPolynomial σ R) (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (r : R) :
    Away.map g f (awayConstOf R f hf r)
      = awayEvalOf R (g f) (g.map_mem hf) (g (MvPolynomial.C r)) := by
  rw [awayConstOf, Away.map_mk,
    awayEvalOf_mk R (g f) (g.map_mem hf) 0 (g (MvPolynomial.C r))
      (g.map_mem (by simp : (MvPolynomial.C r : MvPolynomial σ R)
        ∈ MvPolynomial.homogeneousSubmodule σ R (0 • 1)))]

/-- ★**座標の側** —— 次数 `1` なので `awayEvalOf_mk` が潰す。 -/
theorem awayMap_awayCoordOf (R : Type) [CommRing R]
    (g : (MvPolynomial.homogeneousSubmodule σ R) →+*ᵍ (MvPolynomial.homogeneousSubmodule τ R))
    (f : MvPolynomial σ R) (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) (j : σ) :
    Away.map g f (awayCoordOf R f hf j)
      = awayEvalOf R (g f) (g.map_mem hf) (g (MvPolynomial.X j)) := by
  rw [awayCoordOf, Away.map_mk,
    awayEvalOf_mk R (g f) (g.map_mem hf) 1 (g (MvPolynomial.X j))
      (g.map_mem (by simpa using MvPolynomial.isHomogeneous_X R j))]

/-! ## ★★★★★★★★★段 C2c-1 の (a) -/

/-- ★★★★★★★★★**可換な四角** —— 段 C2c-1 の (a)。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `Away.map g f ∘ awayEvalOf f = awayEvalOf (g f) ∘ g`

★これで「チャートの上での引き戻し」を**多項式環の側の代入**に翻訳できる。
★★機構は生成元（`C r` と `X j`）で見るだけで、
右辺は `awayEvalOf_mk` が次数 `0`・`1` で潰す。 -/
theorem awayMap_comp_awayEvalOf (R : Type) [CommRing R]
    (g : (MvPolynomial.homogeneousSubmodule σ R) →+*ᵍ (MvPolynomial.homogeneousSubmodule τ R))
    (f : MvPolynomial σ R) (hf : f ∈ MvPolynomial.homogeneousSubmodule σ R 1) :
    (Away.map g f).comp (awayEvalOf R f hf)
      = (awayEvalOf R (g f) (g.map_mem hf)).comp (g.toRingHom) := by
  refine MvPolynomial.ringHom_ext (fun r => ?_) (fun j => ?_)
  · show Away.map g f (awayEvalOf R f hf (MvPolynomial.C r)) = _
    rw [awayEvalOf_C]
    exact awayMap_awayConstOf R g f hf r
  · show Away.map g f (awayEvalOf R f hf (MvPolynomial.X j)) = _
    rw [awayEvalOf_X]
    exact awayMap_awayCoordOf R g f hf j

/-! ## ★★★★★超平面の場合 -/

/-- ★★★★★**超平面の場合** —— `g = hyperplaneHom`（`x₀ ↦ 0`）。 -/
theorem awayMap_comp_awayEvalOf_hyperplane (N : ℕ) (R : Type) [CommRing R]
    (f : MvPolynomial (Fin (N+1)) R)
    (hf : f ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R 1) :
    (Away.map (hyperplaneHom N R) f).comp (awayEvalOf R f hf)
      = (awayEvalOf R (hyperplaneHom N R f)
          ((hyperplaneHom N R).map_mem hf)).comp ((hyperplaneHom N R).toRingHom) :=
  awayMap_comp_awayEvalOf R (hyperplaneHom N R) f hf

/-! ## ★出典の紐付け(`.src`) -/

def awayMap_awayCoordOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(座標の側の可換性)",
    sectionId := "genell-prop-1-4" }

def awayMap_comp_awayEvalOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1 の (a)——可換な四角)",
    sectionId := "genell-prop-1-4" }

def awayMap_comp_awayEvalOf_hyperplane.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(可換な四角——超平面の場合)",
    sectionId := "genell-prop-1-4" }

def awayMap_comp_awayEvalOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayEvalOf_mk(斉次式の評価は a/f^k、§9-861c)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayEvalOf_mk") 2,
    .citation "[mathlib]" "HomogeneousLocalization.Away.map_mk"
      (.inMathlib "HomogeneousLocalization.Away.map_mk") 2,
    .implicitStep
      ("★★★★★最初は「g が定数を定数に送る」という仮定が要ると読んでいたが、**要らない**。" ++
       "awayEvalOf_mk が**次数だけ**で右辺を Away.mk に潰すので、" ++
       "g (C r) が何であっても(次数 0 でありさえすれば)両辺は同じ Away.mk になる" ++
       "(2026-08-28 実測)。★これは「一般化した方が証明が短くなる」場合である") 3,
    .implicitStep
      ("★★段 C2c-1 に残るのは (c) である: (a)(b) から " ++
       "ker (Away.map (hyperplaneHom) (x_i)) = span {x_0/x_i} を出す。" ++
       "★awayEvalOf が全射(§9-861c)なので ker は像で測れ、" ++
       "ker σ = (x₀)(§9-861)と ker (awayEvalOf) = (f − 1) を合わせる") 5 ]

end ABC3.Found.GenEll
