/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjSpaceCover
import Mathlib.AlgebraicGeometry.ProjectiveSpectrum.Functor
import ABC3.Meta.Claim

/-!
# `ℙⁿ` の**超平面**とそのイデアル層（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★因子表示なら `O(1)` は**超平面因子**である

台帳の段 C2c は「`O(1)`（twisting sheaf）」であり、mathlib には次数 `n` の
斉次分数の層が無い（次数 0 の `HomogeneousLocalization` だけ）。
★しかし**本プロジェクトは §1 を通して因子表示（`ArithCartier`）で作業している**ので、
そこでは `O(1)` は**超平面因子の類**であり、層そのものは要らない。

★★本ファイルは超平面 `ℙ^{N-1} ↪ ℙ^N` を作り、そのイデアル層を取る:

    `hyperplaneι  : Proj 𝒜(Fin N) ⟶ Proj 𝒜(Fin (N+1))`
    `hyperplaneIdeal ≝ (hyperplaneι).ker`

## ★★★★★★機構は mathlib の `Proj.map` 一本

`AlgebraicGeometry.Proj.map (f : 𝒜 →+*ᵍ ℬ) (hf : ℬ₊ ≤ 𝒜₊.map f) : Proj ℬ ⟶ Proj 𝒜`
が**ある**（2026-08-28 実測、`ProjectiveSpectrum/Functor.lean`）。

★与える次数付き環準同型は `x₀ ↦ 0`、`x_{j+1} ↦ y_j`——すなわち
`MvPolynomial.aeval (Fin.cases 0 X)` である。

### ★★★★★次数付きであることの中身

`aeval g` が斉次性を保つことは mathlib に**無い**（`IsHomogeneous.bind₁` は無い、実測）。
★単項式に分けて示す:

* `aeval g (monomial d r) = C r · ∏ᵢ (g i)^{dᵢ}`
* 各 `g i` は次数 1 の斉次（`0` も次数 1 の斉次である）
* よって積は次数 `Σ dᵢ = deg d` の斉次

### ★★★無関係イデアルの条件は**切断**で出る

`rename Fin.succ` が `aeval g` の切断である（`aeval_hyperGen_rename`）。
★したがって各次数 `n > 0` で `ℬₙ` は `𝒜ₙ` の像に含まれ、
`ℬ₊ ≤ 𝒜₊.map f` が出る。

## ★残っている段（明示）

★★**超平面因子の高さが素朴高さであること**——
`ht_{超平面}(x) = log max|x_i|` ——は別に要る。
★★★そして **very ample**（`n·D` が閉埋め込みによる超平面の引き戻し）と
**Serre の定理**（段 E）が残る。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★★斉次性の道具 -/

/-- ★台の元の次数は斉次次数に等しい。 -/
theorem degree_of_mem_support {σ R : Type} [CommRing R] {n : ℕ} {p : MvPolynomial σ R}
    (hp : p.IsHomogeneous n) {d : σ →₀ ℕ} (hd : d ∈ p.support) : d.degree = n := by
  have h := hp (by simpa using hd)
  rw [Finsupp.degree_eq_weight_one]
  exact h

/-- ★★★**次数 1 の斉次元への代入は斉次性を保つ**。

★mathlib に `IsHomogeneous.bind₁` は無い（2026-08-28 実測）ので、単項式に分けて示す。 -/
theorem aeval_isHomogeneous {σ τ R : Type} [CommRing R] (g : σ → MvPolynomial τ R)
    (hg : ∀ i, (g i).IsHomogeneous 1) {n : ℕ} {p : MvPolynomial σ R}
    (hp : p.IsHomogeneous n) : (MvPolynomial.aeval g p).IsHomogeneous n := by
  classical
  rw [← MvPolynomial.mem_homogeneousSubmodule]
  rw [p.as_sum, map_sum]
  refine Submodule.sum_mem _ (fun d hd => ?_)
  rw [MvPolynomial.mem_homogeneousSubmodule, MvPolynomial.aeval_monomial]
  have hdeg : d.degree = n := degree_of_mem_support hp hd
  have hprod : (d.prod fun i k => g i ^ k).IsHomogeneous d.degree := by
    rw [Finsupp.prod, Finsupp.degree]
    exact MvPolynomial.IsHomogeneous.prod _ _ _
      (fun i _ => by simpa using (hg i).pow (d i))
  have hC : ((algebraMap R (MvPolynomial τ R)) (MvPolynomial.coeff d p)).IsHomogeneous 0 := by
    rw [MvPolynomial.algebraMap_eq]
    exact MvPolynomial.isHomogeneous_C _ _
  have h := hC.mul hprod
  rw [zero_add, hdeg] at h
  exact h

/-- ★★正の次数の斉次元は無関係イデアルに入る。 -/
theorem homogeneous_mem_irrelevant {σ R : Type} [CommRing R] {n : ℕ} (hn : 0 < n)
    (p : MvPolynomial σ R) (hp : p ∈ MvPolynomial.homogeneousSubmodule σ R n) :
    p ∈ HomogeneousIdeal.irrelevant (MvPolynomial.homogeneousSubmodule σ R) := by
  rw [HomogeneousIdeal.mem_irrelevant_iff]
  show (DirectSum.decompose (MvPolynomial.homogeneousSubmodule σ R) p 0 :
    MvPolynomial σ R) = 0
  have h := DirectSum.decompose_of_mem (MvPolynomial.homogeneousSubmodule σ R) hp
  rw [h]
  have h2 := DirectSum.of_eq_of_ne
    (β := fun i => ↥(MvPolynomial.homogeneousSubmodule σ R i)) n 0 ⟨p, hp⟩ (Nat.ne_of_lt hn)
  rw [h2]
  rfl

/-! ## ★★★★★超平面を与える次数付き環準同型 -/

/-- ★`x₀ ↦ 0`、`x_{j+1} ↦ y_j`。 -/
noncomputable def hyperGen (N : ℕ) (R : Type) [CommRing R] :
    Fin (N + 1) → MvPolynomial (Fin N) R :=
  Fin.cases 0 MvPolynomial.X

theorem hyperGen_isHomogeneous (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) :
    (hyperGen N R i).IsHomogeneous 1 := by
  refine Fin.cases ?_ ?_ i
  · show (0 : MvPolynomial (Fin N) R).IsHomogeneous 1
    exact MvPolynomial.isHomogeneous_zero _ _ _
  · intro j
    show (MvPolynomial.X j : MvPolynomial (Fin N) R).IsHomogeneous 1
    exact MvPolynomial.isHomogeneous_X R j

/-- ★★**`rename Fin.succ` は切断である**——これが全射性（無関係イデアルの条件）を出す。 -/
theorem aeval_hyperGen_rename (N : ℕ) (R : Type) [CommRing R] (q : MvPolynomial (Fin N) R) :
    MvPolynomial.aeval (hyperGen N R) (MvPolynomial.rename Fin.succ q) = q := by
  induction q using MvPolynomial.induction_on with
  | C r => simp only [MvPolynomial.rename_C, MvPolynomial.aeval_C]; rfl
  | add p q hp hq => simp only [map_add, hp, hq]
  | mul_X p j hp =>
      simp only [map_mul, MvPolynomial.rename_X, MvPolynomial.aeval_X, hp]
      show p * hyperGen N R (Fin.succ j) = p * MvPolynomial.X j
      rfl

/-- ★★★★★**超平面を与える次数付き環準同型**。 -/
noncomputable def hyperplaneHom (N : ℕ) (R : Type) [CommRing R] :
    (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) →+*ᵍ
      (MvPolynomial.homogeneousSubmodule (Fin N) R) where
  toRingHom := (MvPolynomial.aeval (hyperGen N R)).toRingHom
  map_mem {i x} hx := by
    rw [MvPolynomial.mem_homogeneousSubmodule] at hx ⊢
    exact aeval_isHomogeneous _ (hyperGen_isHomogeneous N R) hx

theorem hyperplaneHom_apply (N : ℕ) (R : Type) [CommRing R] (p : MvPolynomial (Fin (N + 1)) R) :
    hyperplaneHom N R p = MvPolynomial.aeval (hyperGen N R) p := rfl

/-- ★★★**無関係イデアルの条件**——切断 `rename Fin.succ` から出る。 -/
theorem hyperplane_irrelevant_le (N : ℕ) (R : Type) [CommRing R] :
    HomogeneousIdeal.irrelevant (MvPolynomial.homogeneousSubmodule (Fin N) R)
      ≤ HomogeneousIdeal.map (hyperplaneHom N R)
        (HomogeneousIdeal.irrelevant (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)) := by
  classical
  intro q hq
  have hq0 : (DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q 0 :
      MvPolynomial (Fin N) R) = 0 :=
    (HomogeneousIdeal.mem_irrelevant_iff _ q).1 hq
  rw [← DirectSum.sum_support_decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q]
  refine Ideal.sum_mem _ (fun n _ => ?_)
  rcases Nat.eq_zero_or_pos n with rfl | hpos
  · rw [hq0]; exact Ideal.zero_mem _
  · have hc : ((DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q n :
        ↥(MvPolynomial.homogeneousSubmodule (Fin N) R n)) : MvPolynomial (Fin N) R)
        ∈ MvPolynomial.homogeneousSubmodule (Fin N) R n :=
      (DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q n).2
    have hp : MvPolynomial.rename Fin.succ
        ((DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q n :
          ↥(MvPolynomial.homogeneousSubmodule (Fin N) R n)) : MvPolynomial (Fin N) R)
        ∈ MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R n := by
      rw [MvPolynomial.mem_homogeneousSubmodule] at hc ⊢
      exact MvPolynomial.IsHomogeneous.rename_isHomogeneous hc
    have hirr := homogeneous_mem_irrelevant hpos _ hp
    have hmap := Ideal.mem_map_of_mem (hyperplaneHom N R).toRingHom hirr
    rw [show ((hyperplaneHom N R).toRingHom :
          MvPolynomial (Fin (N + 1)) R → MvPolynomial (Fin N) R)
        (MvPolynomial.rename Fin.succ
          ((DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q n :
            ↥(MvPolynomial.homogeneousSubmodule (Fin N) R n)) : MvPolynomial (Fin N) R))
        = ((DirectSum.decompose (MvPolynomial.homogeneousSubmodule (Fin N) R) q n :
            ↥(MvPolynomial.homogeneousSubmodule (Fin N) R n)) : MvPolynomial (Fin N) R)
      from aeval_hyperGen_rename N R _] at hmap
    exact hmap

/-! ## ★★★★★★★★超平面とそのイデアル層 -/

/-- ★★★★★★★★**超平面 `ℙ^{N-1} ↪ ℙ^N`**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は mathlib の `Proj.map` 一本である。 -/
noncomputable def hyperplaneι (N : ℕ) (R : Type) [CommRing R] :
    Proj (MvPolynomial.homogeneousSubmodule (Fin N) R)
      ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R) :=
  Proj.map (hyperplaneHom N R) (hyperplane_irrelevant_le N R)

/-- ★★★★★★★★★**超平面のイデアル層**——因子表示での `O(-1)`。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが「因子表示なら `O(1)` は超平面因子である」の中身である。 -/
noncomputable def hyperplaneIdeal (N : ℕ) (R : Type) [CommRing R] :
    (Proj (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)).IdealSheafData :=
  Scheme.Hom.ker (hyperplaneι N R)

/-! ### ★出典の紐付け(`.src`) -/

def aeval_isHomogeneous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数 1 の斉次元への代入は斉次性を保つ)",
    sectionId := "genell-prop-1-4" }

def hyperplaneι.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——超平面 ℙ^{N-1} ↪ ℙ^N)",
    sectionId := "genell-prop-1-4" }

def hyperplaneIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(語彙——超平面のイデアル層。因子表示での O(-1))",
    sectionId := "genell-prop-1-4" }

def hyperplaneι.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgebraicGeometry.Proj.map(次数付き環準同型から Proj の射)"
      (.inMathlib "AlgebraicGeometry.Proj.map") 6,
    .citation "[mathlib]" "MvPolynomial.IsHomogeneous.rename_isHomogeneous"
      (.inMathlib "MvPolynomial.IsHomogeneous.rename_isHomogeneous") 6,
    .implicitStep
      ("★台帳の段 C2c は O(1)(twisting sheaf)だが、mathlib には次数 n の斉次分数の層が" ++
       "無い(次数 0 の HomogeneousLocalization だけ、2026-08-28 実測)。" ++
       "★★本プロジェクトは §1 を通して因子表示で作業しているので、" ++
       "そこでは O(1) は超平面因子の類であり層そのものは要らない") 6,
    .implicitStep
      ("★★★aeval が斉次性を保つことは mathlib に無い(IsHomogeneous.bind₁ は無い)ので" ++
       "単項式に分けて示した。無関係イデアルの条件は " ++
       "rename Fin.succ が切断であることから出る") 6,
    .implicitStep
      ("★★★★残るのは「超平面因子の高さが素朴高さであること」と、" ++
       "very ample(n·D が閉埋め込みによる超平面の引き戻し)・Serre の定理(段 E)である") 6 ]

end ABC3.Found.GenEll
