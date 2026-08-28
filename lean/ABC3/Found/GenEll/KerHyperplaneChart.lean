/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.AwaySquare
import ABC3.Found.GenEll.AwayKerGen
import ABC3.Found.GenEll.KillX0
import ABC3.Found.GenEll.HyperplaneChart
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 C2c-1 —— 超平面はチャートの上でちょうど `x₀/x_i` で切られる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★これは何か —— 段 C2c-1 が閉じる

`§9-860` で段 C2c-1 の逆向きを 3 段に分けた。3 段とも揃った:

| 段 | 内容 | 場所 |
|---|---|---|
| (a) | 可換な四角 `Φ ∘ Ψ = Ψ′ ∘ σ` | `§9-862`（`AwaySquare`） |
| (b) | `ker σ = (x₀)` | `§9-861`（`KillX0`） |
| **(c)** | **(a)(b) から `ker Φ = (x₀/x_i)`** | ★★**本ファイル** |

★片側（`⊇`）は `§9-858` で既に取ってあった（`span_projCoord_le_ker`）。

## ★★★★★機構 —— 全射を 2 回使う

★`Ψ = awayEvalOf (x_i)` は**全射**（`§9-861c`）なので

    `ker Φ = map Ψ (comap Ψ (ker Φ))`   （`Ideal.map_comap_of_surjective`）

★★`comap Ψ (ker Φ) = ker (Φ ∘ Ψ) = ker (Ψ′ ∘ σ) = comap σ (ker Ψ′)`
であり、`ker Ψ′ = (g(x_i) − 1)`（`§9-862b`）。

★★★`σ` も**全射**（切断 `rename Fin.succ`）なので

    `comap σ (map σ (x_i − 1)) = (x_i − 1) ⊔ ker σ = (x_i − 1) ⊔ (x₀)`
      （`Ideal.comap_map_of_surjective` ＋ `§9-861`）

★★★★最後に `Ψ` で送ると `Ψ(x_i − 1) = x_i/x_i − 1 = 0` が消えて `(x₀/x_i)` だけが残る。

## ★測定の記録

★★★★★「`x_i − 1` の像が消える」ことが、この計算で**唯一の非自明な打ち消し**である
——`A⁰_{x_i}` は `R[x]/(x_i − 1)` なので、`x_i − 1` は既に `0` になっている。
★したがって「超平面 `{x₀ = 0}` のチャート上の方程式が `x₀/x_i = 0` である」ことは、
`ker σ = (x₀)`（`§9-861`）を `Ψ` で押し出しただけである。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★`σ` は全射である -/

/-- ★**`aeval hyperGen` は全射である** —— `rename Fin.succ` が切断だから。 -/
theorem surjective_hyperplaneHom (N : ℕ) (R : Type) [CommRing R] :
    Function.Surjective ((hyperplaneHom N R).toRingHom) := fun q =>
  ⟨MvPolynomial.rename Fin.succ q, aeval_hyperGen_rename N R q⟩

/-- ★**`ker σ = (x₀)`**（`§9-861` を `hyperplaneHom` の言葉に直したもの）。 -/
theorem ker_hyperplaneHom (N : ℕ) (R : Type) [CommRing R] :
    RingHom.ker ((hyperplaneHom N R).toRingHom)
      = Ideal.span {(MvPolynomial.X (0 : Fin (N+1)) : MvPolynomial (Fin (N+1)) R)} :=
  ker_aeval_hyperGen N R

/-! ## ★★★★★★★★★★段 C2c-1 の本体 -/

/-- ★★★★★★★★★★**段 C2c-1** —— チャートの上で超平面はちょうど `x₀/x_i` で切られる。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `ker (Away.map (hyperplaneHom) (x_i)) = (x₀/x_i)`

★`§9-858` は `⊇` しか取っていなかった。★★本補題が `⊆` を埋め、段 C2c-1 が閉じる。 -/
theorem ker_awayMap_hyperplane (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) :
    RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i))
      = Ideal.span {projCoord N R i 0} := by
  have hf : (MvPolynomial.X i : MvPolynomial (Fin (N+1)) R)
      ∈ MvPolynomial.homogeneousSubmodule (Fin (N+1)) R 1 := by
    simpa using MvPolynomial.isHomogeneous_X R i
  have hΨsurj : Function.Surjective (awayEvalOf R (MvPolynomial.X i) hf) :=
    awayEvalOf_surjective R _ hf
  -- `comap Ψ (ker Φ) = (x_i − 1) ⊔ (x₀)`
  have hcomap :
      Ideal.comap (awayEvalOf R (MvPolynomial.X i) hf)
          (RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i)))
        = Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)}
          ⊔ Ideal.span {(MvPolynomial.X (0 : Fin (N+1)) : MvPolynomial (Fin (N+1)) R)} := by
    have h1 : Ideal.comap (awayEvalOf R (MvPolynomial.X i) hf)
        (RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i)))
        = RingHom.ker ((Away.map (hyperplaneHom N R) (MvPolynomial.X i)).comp
            (awayEvalOf R (MvPolynomial.X i) hf)) := by
      rw [RingHom.ker, RingHom.ker, Ideal.comap_comap]
    rw [h1, awayMap_comp_awayEvalOf R (hyperplaneHom N R) (MvPolynomial.X i) hf]
    have h2 : RingHom.ker (((awayEvalOf R (hyperplaneHom N R (MvPolynomial.X i))
          ((hyperplaneHom N R).map_mem hf))).comp ((hyperplaneHom N R).toRingHom))
        = Ideal.comap ((hyperplaneHom N R).toRingHom)
            (RingHom.ker (awayEvalOf R (hyperplaneHom N R (MvPolynomial.X i))
              ((hyperplaneHom N R).map_mem hf))) := by
      rw [RingHom.ker, RingHom.ker, Ideal.comap_comap]
    rw [h2, ker_awayEvalOf]
    have h3 : Ideal.span {(hyperplaneHom N R (MvPolynomial.X i) - 1 :
          MvPolynomial (Fin N) R)}
        = Ideal.map ((hyperplaneHom N R).toRingHom)
            (Ideal.span {(MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R)}) := by
      rw [Ideal.map_span, Set.image_singleton, map_sub, map_one]
      rfl
    rw [h3, Ideal.comap_map_of_surjective _ (surjective_hyperplaneHom N R)]
    congr 1
    exact ker_hyperplaneHom N R
  -- `ker Φ = map Ψ (comap Ψ (ker Φ))`
  rw [← Ideal.map_comap_of_surjective (awayEvalOf R (MvPolynomial.X i) hf) hΨsurj
    (RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i))), hcomap,
    Ideal.map_sup, Ideal.map_span, Ideal.map_span, Set.image_singleton, Set.image_singleton]
  have hzero : (awayEvalOf R (MvPolynomial.X i) hf)
      (MvPolynomial.X i - 1 : MvPolynomial (Fin (N+1)) R) = 0 := by
    rw [map_sub, map_one, awayEvalOf_self, sub_self]
  have hx0 : (awayEvalOf R (MvPolynomial.X i) hf)
      (MvPolynomial.X (0 : Fin (N+1)) : MvPolynomial (Fin (N+1)) R)
      = projCoord N R i 0 := by
    rw [awayEvalOf_X]
    rfl
  rw [hzero, hx0, Ideal.span_singleton_eq_bot.2 rfl, bot_sup_eq]

/-- ★★**`§9-858` の片側は本補題の系である** —— 記録として残す。 -/
theorem span_projCoord_le_ker' (N : ℕ) (R : Type) [CommRing R] (i : Fin (N + 1)) :
    Ideal.span {projCoord N R i 0}
      ≤ RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i)) :=
  (ker_awayMap_hyperplane N R i).ge

/-! ## ★出典の紐付け(`.src`) -/

def surjective_hyperplaneHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(aeval hyperGen は全射である)",
    sectionId := "genell-prop-1-4" }

def ker_awayMap_hyperplane.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1——超平面はチャート上で x₀/x_i で切られる)",
    sectionId := "genell-prop-1-4" }

def ker_awayMap_hyperplane.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "awayMap_comp_awayEvalOf(可換な四角、段 C2c-1 の (a)、§9-862)"
      (.inProject "ABC3" "ABC3.Found.GenEll.awayMap_comp_awayEvalOf") 2,
    .citation "[ABC3]" "ker_aeval_hyperGen(ker σ = (x₀)、段 C2c-1 の (b)、§9-861)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_aeval_hyperGen") 2,
    .citation "[ABC3]" "ker_awayEvalOf(ker (awayEvalOf f) = (f − 1)、§9-862b)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_awayEvalOf") 2,
    .citation "[mathlib]" "Ideal.map_comap_of_surjective / Ideal.comap_map_of_surjective"
      (.inMathlib "Ideal.comap_map_of_surjective") 2,
    .implicitStep
      ("★★★★★「x_i − 1 の像が消える」ことが、この計算で**唯一の非自明な打ち消し**である" ++
       "——A⁰_{x_i} は R[x]/(x_i − 1) なので x_i − 1 は既に 0 になっている。" ++
       "★したがって「超平面 {x₀ = 0} のチャート上の方程式が x₀/x_i = 0 である」ことは、" ++
       "ker σ = (x₀)(§9-861)を Ψ で押し出しただけである") 3,
    .implicitStep
      ("★★段 C2c-1 が閉じたので、残るのは段 C2c-3" ++
       "(Fubini–Study 型の計量が Definition 1.1 を満たすこと)と " ++
       "E0 / E1d / F2c である") 4 ]

end ABC3.Found.GenEll
