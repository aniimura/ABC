/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjectiveSpace
import ABC3.Meta.Claim

/-!
# `ℙⁿ` の**標準アフィン被覆**と**座標比 `x_j/x_i`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★消費側が要るのは `O(1)` ではなく**座標**である

`Found/GenEll/NorthcottCoord.lean` の `northcott_of_projModel` が要求するのは

| 受けるもの | 中身 |
|---|---|
| `crd p : ι → fld p` | ★**同次座標** |
| `idx : ι` | 割る成分 |
| `hinj` | 正規化座標 `crd p j / crd p idx` が単射 |

であって、`O(1)` そのものではない。★★**正規化座標 `x_j/x_i` は `ℙⁿ` の標準アフィン
チャートから直に取れる**——`D₊(x_i) ≅ Spec (A⁰_{x_i})` であり、
`A⁰_{x_i}` は `x_j/x_i` たちで生成される。

★★★本ファイルはその 2 本を作る:

* **`iSup_basicOpen_X_eq_top`** —— `D₊(x_i)` が `ℙⁿ` を覆うこと
* **`projCoord`** —— `x_j/x_i ∈ A⁰_{x_i}`

## ★★★★★被覆の中身は「無関係イデアルは変数で生成される」

`z : ProjectiveSpectrum 𝒜` は**無関係イデアルを含まない**斉次素イデアルである。
★もし全ての `x_i` が `z` に入れば `span{x_i} ⊆ z` となり、
無関係イデアル `⨁_{n>0} 𝒜ₙ` も `z` に入ってしまう。

★★その包含 `irrelevant ≤ span{x_i}` は

    `𝒜ₙ = (𝒜₁)ⁿ`（`homogeneousSubmodule_one_pow`）、
    `𝒜₁ = span_R{x_i}`（`homogeneousSubmodule_one_eq_span_X`）

から `n ≥ 1` で出る。★★★`R`-span と `Ideal.span` は係数環が違うので
`Submodule.restrictScalars` を挟む。

## ★残っている段（明示）

★★**チャートの同型** `Proj| D₊(f) ≅ Spec (A⁰_f)` は mathlib にある
（`AlgebraicGeometry.projIsoSpec`、ただし `LocallyRingedSpace` の水準）。
★★★残るのは「点 `Spec F ⟶ ℙⁿ` をチャートに落として環準同型 `A⁰_{x_i} → F` を取る」
配管と、very ample からの閉埋め込み（段 E）である。
-/

namespace ABC3.Found.GenEll

open MvPolynomial ProjectiveSpectrum AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★正の次数の斉次多項式は変数で生成される -/

/-- ★**正の次数の斉次多項式は変数の生成するイデアルに入る**。 -/
theorem homogeneous_mem_span_X {σ R : Type} [CommRing R] {n : ℕ} (hn : 0 < n)
    (p : MvPolynomial σ R) (hp : p ∈ MvPolynomial.homogeneousSubmodule σ R n) :
    p ∈ Ideal.span (Set.range (MvPolynomial.X : σ → MvPolynomial σ R)) := by
  have hspan : MvPolynomial.homogeneousSubmodule σ R 1
      ≤ (Ideal.span (Set.range (MvPolynomial.X : σ → MvPolynomial σ R))).restrictScalars R := by
    rw [MvPolynomial.homogeneousSubmodule_one_eq_span_X]
    exact Submodule.span_le.2 (fun x hx => Ideal.subset_span hx)
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, (Nat.succ_pred_eq_of_pos hn).symm⟩
  rw [← MvPolynomial.homogeneousSubmodule_one_pow, pow_succ] at hp
  refine Submodule.mul_induction_on hp ?_ ?_
  · intro a _ b hb
    exact Ideal.mul_mem_left _ _ (hspan hb)
  · intro x y hx hy
    exact Ideal.add_mem _ hx hy

/-- ★★**無関係イデアルは変数で生成されるイデアルに含まれる**。 -/
theorem irrelevant_le_span_X (σ R : Type) [CommRing R] :
    (HomogeneousIdeal.irrelevant (MvPolynomial.homogeneousSubmodule σ R)).toIdeal
      ≤ Ideal.span (Set.range (MvPolynomial.X : σ → MvPolynomial σ R)) := by
  classical
  intro x hx
  have hx0 : (DirectSum.decompose (MvPolynomial.homogeneousSubmodule σ R) x 0 :
      MvPolynomial σ R) = 0 :=
    (HomogeneousIdeal.mem_irrelevant_iff (MvPolynomial.homogeneousSubmodule σ R) x).1 hx
  rw [← DirectSum.sum_support_decompose (MvPolynomial.homogeneousSubmodule σ R) x]
  refine Ideal.sum_mem _ (fun i _ => ?_)
  rcases Nat.eq_zero_or_pos i with rfl | hpos
  · rw [hx0]
    exact Ideal.zero_mem _
  · exact homogeneous_mem_span_X hpos _
      (DirectSum.decompose (MvPolynomial.homogeneousSubmodule σ R) x i).2

/-! ## ★★★★★★★★標準アフィン被覆 -/

/-- ★★★**どの点でもある変数が消えない**。 -/
theorem exists_X_notMem (σ R : Type) [CommRing R]
    (z : ProjectiveSpectrum (MvPolynomial.homogeneousSubmodule σ R)) :
    ∃ i : σ, (MvPolynomial.X i : MvPolynomial σ R) ∉ z.asHomogeneousIdeal := by
  by_contra h
  push Not at h
  refine z.not_irrelevant_le ?_
  intro x hx
  have h1 : Ideal.span (Set.range (MvPolynomial.X : σ → MvPolynomial σ R))
      ≤ z.asHomogeneousIdeal.toIdeal := by
    rw [Ideal.span_le]
    rintro _ ⟨i, rfl⟩
    exact h i
  exact h1 (irrelevant_le_span_X σ R hx)

/-- ★★★★★★★★**`ℙⁿ` の標準アフィン被覆** —— `D₊(x_i)` が全体を覆う。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが「どの点にも `x_i ≠ 0` な成分がある」＝**正規化座標が取れる**ことの根拠である。 -/
theorem iSup_basicOpen_X_eq_top (σ R : Type) [CommRing R] :
    (⨆ i : σ, ProjectiveSpectrum.basicOpen (MvPolynomial.homogeneousSubmodule σ R)
      (MvPolynomial.X i)) = ⊤ := by
  refine eq_top_iff.2 (fun z _ => ?_)
  obtain ⟨i, hi⟩ := exists_X_notMem σ R z
  exact TopologicalSpace.Opens.mem_iSup.2 ⟨i, hi⟩

/-! ## ★★★★★★★座標比 `x_j/x_i` -/

/-- ★★★★★★★**正規化座標 `x_j/x_i`** —— `D₊(x_i)` のアフィン環 `A⁰_{x_i}` の元。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`HomogeneousLocalization.Away.mk 𝒜 hf n x hx` は分数 `x / f ^ n` である。
`n = 1`、`f = x_i`、`x = x_j` と取ればよい。 -/
noncomputable def projCoord (N : ℕ) (R : Type) [CommRing R] (i j : Fin (N + 1)) :
    HomogeneousLocalization.Away
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i) :=
  HomogeneousLocalization.Away.mk _ (MvPolynomial.isHomogeneous_X R i) 1 (MvPolynomial.X j)
    (by simpa using MvPolynomial.isHomogeneous_X R j)

/-- ★★**割る成分は `1` である**——`x_i/x_i = 1`。 -/
theorem projCoord_self (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) :
    projCoord N R i i = 1 := by
  apply HomogeneousLocalization.val_injective
  rw [projCoord, HomogeneousLocalization.Away.val_mk]
  simp

/-! ## ★★★★★★★★点をチャートへ落とす（段 C2b） -/

/-- ★★★**体値の点はどれかのチャートに入る**。

★`Spec F` は 1 点なので、`exists_X_notMem` がそのまま効く。 -/
theorem exists_chart_range (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)) :
    ∃ i : Fin (N + 1), Set.range x.base ⊆
      Set.range (Proj.awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base := by
  haveI : Nonempty (Spec (CommRingCat.of F)) :=
    inferInstanceAs (Nonempty (PrimeSpectrum F))
  haveI : Subsingleton (Spec (CommRingCat.of F)) :=
    inferInstanceAs (Subsingleton (PrimeSpectrum F))
  obtain ⟨pt⟩ := ‹Nonempty (Spec (CommRingCat.of F))›
  obtain ⟨i, hi⟩ := exists_X_notMem (Fin (N + 1)) R (x.base pt)
  refine ⟨i, ?_⟩
  rintro _ ⟨q, rfl⟩
  have hq : q = pt := Subsingleton.elim _ _
  subst hq
  have hr := Proj.opensRange_awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
    (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos
  have h2 : x.base q ∈ Proj.basicOpen (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.X i) := hi
  rw [← hr] at h2
  exact h2

/-- ★★★★★★**点をチャートに落として得られる環準同型** `A⁰_{x_i} → F`。

★`Proj.awayι` は開埋め込みなので `IsOpenImmersion.lift` が使え、
`Spec.preimage` で環準同型に直せる。 -/
noncomputable def projChartHom (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (i : Fin (N + 1))
    (hx : Set.range x.base ⊆
      Set.range (Proj.awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base) :
    CommRingCat.of (HomogeneousLocalization.Away
        (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) (MvPolynomial.X i))
      ⟶ CommRingCat.of F :=
  Spec.preimage (IsOpenImmersion.lift
    (Proj.awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos) x hx)

/-- ★★★★★★★★**点の正規化座標** `x_j/x_i ∈ F`。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが `northcott_of_projModel` の `crd p j / crd p idx` そのものである
——`idx = i` と取れば分母は `1` になる（`projPointCoord_self`）。 -/
noncomputable def projPointCoord (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (i : Fin (N + 1))
    (hx : Set.range x.base ⊆
      Set.range (Proj.awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base)
    (j : Fin (N + 1)) : F :=
  (projChartHom N R F x i hx).hom (projCoord N R i j)

/-- ★★★**割る成分の座標は `1`** —— 正規化されている。 -/
theorem projPointCoord_self (N : ℕ) (R : Type) [CommRing R] (F : Type) [Field F]
    (x : Spec (CommRingCat.of F) ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R))
    (i : Fin (N + 1))
    (hx : Set.range x.base ⊆
      Set.range (Proj.awayι (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
        (MvPolynomial.X i) (MvPolynomial.isHomogeneous_X R i) one_pos).base) :
    projPointCoord N R F x i hx i = 1 := by
  show (projChartHom N R F x i hx).hom (projCoord N R i i) = 1
  rw [projCoord_self, map_one]

/-! ### ★出典の紐付け(`.src`) -/

def iSup_basicOpen_X_eq_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——ℙⁿ の標準アフィン被覆 D₊(x_i))",
    sectionId := "genell-prop-1-4" }

def projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——正規化座標 x_j/x_i)",
    sectionId := "genell-prop-1-4" }

def projPointCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(点の正規化座標 x_j/x_i ∈ F)",
    sectionId := "genell-prop-1-4" }

def iSup_basicOpen_X_eq_top.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "MvPolynomial.homogeneousSubmodule_one_eq_span_X(𝒜₁ = span_R{x_i})"
      (.inMathlib "MvPolynomial.homogeneousSubmodule_one_eq_span_X") 6,
    .citation "[mathlib]" "MvPolynomial.homogeneousSubmodule_one_pow(𝒜ₙ = (𝒜₁)ⁿ)"
      (.inMathlib "MvPolynomial.homogeneousSubmodule_one_pow") 6,
    .citation "[mathlib]" "ProjectiveSpectrum.not_irrelevant_le(点は無関係イデアルを含まない)"
      (.inMathlib "ProjectiveSpectrum") 6,
    .implicitStep
      ("★消費側(northcott_of_projModel)が要るのは **同次座標** であって O(1) ではない。" ++
       "★★正規化座標 x_j/x_i は ℙⁿ の標準アフィンチャートから直に取れる") 6,
    .implicitStep
      ("★★★残るのは「点 Spec F ⟶ ℙⁿ をチャートに落として環準同型 A⁰_{x_i} → F を取る」" ++
       "配管と、very ample からの閉埋め込み(段 E)である。" ++
       "チャートの同型 Proj| D₊(f) ≅ Spec (A⁰_f) は mathlib の projIsoSpec にある" ++
       "(ただし LocallyRingedSpace の水準)") 6 ]

end ABC3.Found.GenEll
