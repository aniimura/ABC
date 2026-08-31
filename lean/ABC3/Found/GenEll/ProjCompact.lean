/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Hyperplane
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim

/-!
# ★★★★★★★`ℙᴺ` は擬コンパクト —— 段 C2c-1 の入口（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★これは何か —— 段 C2c-1 の入口

段 C2c-1 は「超平面因子の引き戻しが座標の生成するイデアルであること」である。
★そのためには**チャートの上でイデアル層を読む**必要があり、mathlib の

    `Scheme.Hom.ker_apply (f) [QuasiCompact f] (U : Y.affineOpens)`
      : `f.ker.ideal U = RingHom.ker (f.app U)`

が使える。★★しかし `QuasiCompact (hyperplaneι N R)` の instance が無かった
（2026-08-28 実測）。

★★★本ファイルはそこを埋める:

| 補題 | 内容 |
|---|---|
| `iSup_projBasicOpen_X_eq_top'` | `D₊(x_j)` は `Proj 𝒜(Fin N)` を覆う（任意の `N`） |
| `projCompactSpace` | ★**`Proj 𝒜(Fin N)` はコンパクト**（有限個のアフィン開で覆えるから） |
| `quasiCompact_hyperplaneι` | ★★したがって `hyperplaneι` は擬コンパクト |
| `hyperplaneIdeal_apply` | ★★★**チャートの上で超平面イデアルを読む** |

## ★測定の記録

★`CompactSpace (Proj 𝒜)` は mathlib に**無い**（2026-08-28 実測）。
`Proj` が擬コンパクトなのは無関係イデアルが有限生成のときだけなので、
一般には成り立たないからである。
★★多項式環の場合は `D₊(x_j)`（`j : Fin N`、**有限個**）で覆えるので出る
——`§9-C2a` の `exists_X_notMem` がそのまま効く。

## ★残っている段（明示）

★★★`RingHom.ker ((hyperplaneι).app (D₊(x_i)))` が
**`x_0/x_i` で生成される**ことを見るのが段 C2c-1 の本体である。
★`Proj.awayι` で `Γ(D₊(x_i)) ≅ A⁰_{x_i}` と読み替え、
`hyperplaneHom`（`x₀ ↦ 0`）の効果を追う。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory ProjectiveSpectrum MvPolynomial

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★標準チャートは `Proj 𝒜(Fin N)` を覆う -/

/-- ★**`D₊(x_j)` は `Proj 𝒜(Fin N)` を覆う**（任意の `N`）。

★`§9-834` の `iSup_projBasicOpen_X_eq_top` は `Fin (N+1)` の形だったが、
本ファイルは `Fin N` の形が要る（`hyperplaneι` の**始域**が `Fin N` だから）。 -/
theorem iSup_projBasicOpen_X_eq_top' (N : ℕ) (R : Type) [CommRing R] :
    (⨆ j : Fin N, Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin N) R)
      (MvPolynomial.X j)) = ⊤ := by
  refine eq_top_iff.2 (fun z _ => TopologicalSpace.Opens.mem_iSup.2 ?_)
  obtain ⟨i, hi⟩ := exists_X_notMem (Fin N) R z
  exact ⟨i, hi⟩

/-! ## ★★★★★★`Proj 𝒜(Fin N)` はコンパクトである -/

/-- ★★★★★★**`Proj 𝒜(Fin N)` はコンパクトである**。

★`D₊(x_j)`（`j : Fin N`、**有限個**）で覆え、各々がアフィン開だから。

## ★測定の記録

`CompactSpace (Proj 𝒜)` は mathlib に**無い**（2026-08-28 実測）
——`Proj` が擬コンパクトなのは無関係イデアルが有限生成のときだけで、一般には成り立たない。
★多項式環の場合は変数が有限個なので出る。 -/
instance projCompactSpace (N : ℕ) (R : Type) [CommRing R] :
    CompactSpace (Proj (MvPolynomial.homogeneousSubmodule (Fin N) R)) := by
  constructor
  have h : (Set.univ : Set (Proj (MvPolynomial.homogeneousSubmodule (Fin N) R)))
      = ⋃ j : Fin N, (Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin N) R)
        (MvPolynomial.X j) : Set _) := by
    have := iSup_projBasicOpen_X_eq_top' N R
    rw [← TopologicalSpace.Opens.coe_top, ← this]
    simp
  rw [h]
  exact isCompact_iUnion (fun j =>
    (Proj.isAffineOpen_basicOpen _ (MvPolynomial.X j)
      (MvPolynomial.isHomogeneous_X R j) one_pos).isCompact)

/-- ★★**したがって `hyperplaneι` は擬コンパクトである**。

★始域がコンパクトで終域が擬分離だから（mathlib の instance）。 -/
instance quasiCompact_hyperplaneι (N : ℕ) (R : Type) [CommRing R] :
    QuasiCompact (hyperplaneι N R) := by infer_instance

/-! ## ★★★★★★★チャートの上で超平面イデアルを読む -/

/-- ★★★★★★★**チャートの上で超平面イデアルを読む** —— 段 C2c-1 の入口。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `hyperplaneIdeal.ideal U = ker ((hyperplaneι).app U)`   （`U` アフィン開）

★★これで「超平面が `x_0/x_i` で切られる」ことを**環の言葉で**言えるようになった。
★★★機構は mathlib の `Scheme.Hom.ker_apply` に `QuasiCompact` を渡すだけである
——その instance が無かったので、`Proj 𝒜(Fin N)` のコンパクト性から作った。 -/
theorem hyperplaneIdeal_apply (N : ℕ) (R : Type) [CommRing R]
    (U : (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)).affineOpens) :
    (hyperplaneIdeal N R).ideal U
      = RingHom.ker ((hyperplaneι N R).app (U : (Proj
          (MvPolynomial.homogeneousSubmodule (Fin (N+1)) R)).Opens)).hom :=
  Scheme.Hom.ker_apply (hyperplaneι N R) U

/-! ## ★出典の紐付け(`.src`) -/

def iSup_projBasicOpen_X_eq_top'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(D₊(x_j) は Proj 𝒜(Fin N) を覆う)",
    sectionId := "genell-prop-1-4" }

def projCompactSpace.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Proj 𝒜(Fin N) はコンパクトである)",
    sectionId := "genell-prop-1-4" }

def hyperplaneIdeal_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートの上で超平面イデアルを読む)",
    sectionId := "genell-prop-1-4" }

def hyperplaneIdeal_apply.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.ker_apply(アフィン開でイデアル層を読む)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.ker_apply") 2,
    .citation "[ABC3]" "exists_X_notMem(無関係イデアルは変数で生成される、段 C2a)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_X_notMem") 2,
    .implicitStep
      ("★CompactSpace (Proj 𝒜) は mathlib に**無い**(2026-08-28 実測)" ++
       "——Proj が擬コンパクトなのは無関係イデアルが有限生成のときだけで、一般には成り立たない。" ++
       "★★多項式環の場合は D₊(x_j)(j : Fin N、有限個)で覆えるので出る") 3,
    .implicitStep
      ("★★★RingHom.ker ((hyperplaneι).app (D₊(x_i))) が **x_0/x_i で生成される**ことを" ++
       "見るのが段 C2c-1 の本体である。" ++
       "★Proj.awayι で Γ(D₊(x_i)) ≅ A⁰_{x_i} と読み替え、hyperplaneHom(x₀ ↦ 0)の効果を追う") 6 ]

end ABC3.Found.GenEll
