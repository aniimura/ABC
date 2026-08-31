/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Hyperplane
import ABC3.Meta.Claim

/-!
# ★★★★★★★★`x₀ ↦ 0` の核は `(x₀)` である —— 段 C2c-1 の (b)（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは何か —— 段 C2c-1 の (b)

`§9-860` で段 C2c-1 の逆向きを 3 段に分けた:

| 段 | 内容 | 状態 |
|---|---|---|
| (a) | 可換な四角 `Φ ∘ Ψ = Ψ′ ∘ σ` | ★開 |
| **(b)** | **`ker σ = (x₀)`**（`σ = aeval hyperGen`） | ★★**本ファイル** |
| (c) | (a)(b) から `ker Φ = (x₀/x_i)` | ★開 |

## ★★★機構 —— 自己代入で差を測る

★`σ` は `x₀ ↦ 0`、`x_{j+1} ↦ y_j` だが、これに `rename Fin.succ` を後置すると
**同じ環の自己代入** `killX0`（`x₀ ↦ 0`、`x_{j+1} ↦ x_{j+1}`）になる。

★★そして `p − killX0 p ∈ (x₀)` が**単項式の帰納法**で出る:

* `C r`: 差は `0`
* `p + q`: 和は閉じている
* `p · x₀`: `killX0 (p·x₀) = killX0 p · 0 = 0` なので差は `p·x₀ ∈ (x₀)`
* `p · x_{j+1}`: 差は `(p − killX0 p)·x_{j+1}`（帰納法の仮定）

★★★`σ p = 0` なら `killX0 p = rename Fin.succ 0 = 0` なので `p = p − killX0 p ∈ (x₀)`。

## ★測定の記録

★`MvPolynomial.aeval (hyperGen N R)` は係数環 `R` が**推論されない**
（`MvPolynomial (Fin N) R` も候補になる）ので `(R := R)` を明示する必要がある
（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory

/-! ## ★自己代入 `x₀ ↦ 0` -/

/-- ★**`x₀ ↦ 0`、`x_{j+1} ↦ x_{j+1}` の自己代入**（`σ` に `rename Fin.succ` を後置したもの）。 -/
noncomputable def killX0 (N : ℕ) (R : Type) [CommRing R] :
    MvPolynomial (Fin (N + 1)) R →+* MvPolynomial (Fin (N + 1)) R :=
  ((MvPolynomial.rename (Fin.succ : Fin N → Fin (N+1))).toRingHom).comp
    ((MvPolynomial.aeval (hyperGen N R)).toRingHom)

theorem killX0_X_zero (N : ℕ) (R : Type) [CommRing R] :
    killX0 N R (MvPolynomial.X (0 : Fin (N+1))) = 0 := by
  show MvPolynomial.rename Fin.succ (MvPolynomial.aeval (hyperGen N R)
    (MvPolynomial.X (0 : Fin (N+1)))) = 0
  rw [MvPolynomial.aeval_X]
  show MvPolynomial.rename Fin.succ (hyperGen N R 0) = 0
  rw [hyperGen]
  simp

theorem killX0_X_succ (N : ℕ) (R : Type) [CommRing R] (j : Fin N) :
    killX0 N R (MvPolynomial.X (Fin.succ j)) = MvPolynomial.X (Fin.succ j) := by
  show MvPolynomial.rename Fin.succ (MvPolynomial.aeval (hyperGen N R)
    (MvPolynomial.X (Fin.succ j))) = _
  rw [MvPolynomial.aeval_X]
  show MvPolynomial.rename Fin.succ (hyperGen N R (Fin.succ j)) = _
  show MvPolynomial.rename Fin.succ (MvPolynomial.X j) = _
  rw [MvPolynomial.rename_X]

theorem killX0_C (N : ℕ) (R : Type) [CommRing R] (r : R) :
    killX0 N R (MvPolynomial.C r) = MvPolynomial.C r := by
  show MvPolynomial.rename Fin.succ (MvPolynomial.aeval (hyperGen N R)
    (MvPolynomial.C r)) = _
  rw [MvPolynomial.aeval_C, MvPolynomial.algebraMap_eq, MvPolynomial.rename_C]

/-! ## ★★★★★差は `(x₀)` に入る -/

/-- ★★★★★**`p − killX0 p ∈ (x₀)`** —— 単項式の帰納法だけである。 -/
theorem sub_killX0_mem (N : ℕ) (R : Type) [CommRing R] (p : MvPolynomial (Fin (N+1)) R) :
    p - killX0 N R p
      ∈ Ideal.span {(MvPolynomial.X (0 : Fin (N+1)) : MvPolynomial (Fin (N+1)) R)} := by
  induction p using MvPolynomial.induction_on with
  | C r => rw [killX0_C, sub_self]; exact Ideal.zero_mem _
  | add p q hp hq =>
      rw [map_add, show p + q - (killX0 N R p + killX0 N R q)
        = (p - killX0 N R p) + (q - killX0 N R q) from by ring]
      exact Ideal.add_mem _ hp hq
  | mul_X p j hp =>
      rw [map_mul]
      refine Fin.cases ?_ ?_ j
      · rw [killX0_X_zero, mul_zero, sub_zero]
        exact Ideal.mul_mem_left _ _ (Ideal.subset_span rfl)
      · intro k
        rw [killX0_X_succ, show p * MvPolynomial.X (Fin.succ k)
          - killX0 N R p * MvPolynomial.X (Fin.succ k)
          = (p - killX0 N R p) * MvPolynomial.X (Fin.succ k) from by ring]
        exact Ideal.mul_mem_right _ _ hp

/-! ## ★★★★★★★★段 C2c-1 の (b) -/

theorem aeval_hyperGen_X_zero (N : ℕ) (R : Type) [CommRing R] :
    MvPolynomial.aeval (R := R) (hyperGen N R) (MvPolynomial.X (0 : Fin (N+1))) = 0 := by
  rw [MvPolynomial.aeval_X]
  show hyperGen N R 0 = 0
  rw [hyperGen]
  simp

/-- ★★★★★★★★**`ker (aeval hyperGen) = (x₀)`** —— 段 C2c-1 の (b)。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★これが「超平面を切る式は `x₀` だけである」ことの多項式環の側の言い方である。
★★機構は `sub_killX0_mem`（差は `(x₀)` に入る）だけである
——`σ p = 0` なら `killX0 p = 0` なので `p = p − killX0 p`。 -/
theorem ker_aeval_hyperGen (N : ℕ) (R : Type) [CommRing R] :
    RingHom.ker ((MvPolynomial.aeval (R := R) (hyperGen N R)).toRingHom)
      = Ideal.span {(MvPolynomial.X (0 : Fin (N+1)) : MvPolynomial (Fin (N+1)) R)} := by
  apply le_antisymm
  · intro p hp
    have hp0 : MvPolynomial.aeval (R := R) (hyperGen N R) p = 0 := hp
    have hk : killX0 N R p = 0 := by
      show MvPolynomial.rename Fin.succ
        (MvPolynomial.aeval (R := R) (hyperGen N R) p) = 0
      rw [hp0, map_zero]
    have h := sub_killX0_mem N R p
    rwa [hk, sub_zero] at h
  · rw [Ideal.span_le]
    rintro z rfl
    rw [SetLike.mem_coe, RingHom.mem_ker]
    exact aeval_hyperGen_X_zero N R

/-! ## ★出典の紐付け(`.src`) -/

def killX0.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(x₀ ↦ 0 の自己代入)",
    sectionId := "genell-prop-1-4" }

def sub_killX0_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(p − killX0 p は (x₀) に入る)",
    sectionId := "genell-prop-1-4" }

def ker_aeval_hyperGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1 の (b)——ker (aeval hyperGen) = (x₀))",
    sectionId := "genell-prop-1-4" }

def ker_aeval_hyperGen.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "hyperGen(x₀ ↦ 0、x_{j+1} ↦ y_j、§9-Hyperplane)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hyperGen") 2,
    .citation "[mathlib]" "MvPolynomial.induction_on(単項式の帰納法)"
      (.inMathlib "MvPolynomial.induction_on") 2,
    .implicitStep
      ("★MvPolynomial.aeval (hyperGen N R) は係数環 R が**推論されない**" ++
       "(MvPolynomial (Fin N) R も候補になる)ので (R := R) を明示する必要がある" ++
       "(2026-08-28 実測)") 2,
    .implicitStep
      ("★★段 C2c-1 の残りは (a)可換な四角 Φ ∘ Ψ = Ψ′ ∘ σ と " ++
       "(c)そこから ker Φ = (x₀/x_i) を出す段である(§9-860 の測定)") 5 ]

end ABC3.Found.GenEll
