/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArithCartierNpow
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★Northcott は射影埋め込みに沿って移送される（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★これは何か —— 一般の `X` への移送

`§9-878`（段 C2e）は `ℙᴺ_ℤ` の超平面因子についての Northcott だった。
★`§9-879`（`ht_{ψ^*E}(x) = ht_E(ψ ∘ x)`）を使うと、
**射 `ψ : X ⟶ ℙᴺ` に沿ってそのまま移送できる**。

    `{p | ht_{ψ^*E}(x_p) ≤ C}` は有限

★★これが原文の「`L_ℚ` が豊富 ⟹ 射影埋め込み ⟹ Northcott」の**移送の段**である。

## ★★★仮定の読み方

| 仮定 | 意味 |
|---|---|
| `hdiv`・`hcont` | `E` は `ℙᴺ` の超平面因子で、計量は Fubini–Study に対して連続 |
| `hdeg` | 次数が `d` 以下 |
| `hcomp` | 各点を `ψ` で送ったものが**整な同次座標**を持つ |
| `hinj` | その座標が点を分ける |

★`hinj` が「`X` の点を分ける」条件になっているのが要点で、
これは `ψ` が**埋め込み**であることの帰結である（`§9-849`・`§9-851`）。

## ★残っている段（明示）

★★一般の `X` への還元に残るのは **`ψ^*(超平面) = n·D`**（段 E3、`L^n` が非常に豊富）だけである。
★`§9-880` の `ht_{D^n} = n·ht_D` と本ファイルを繋げば、`ht_D` についての Northcott になる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-- ★★★★★★★★★★★★★★**Northcott は射に沿って移送される**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-878`（`ℙᴺ` の Northcott）＋ `§9-879`（高さの関手性）である。 -/
theorem northcott_comap (N d : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {X : Scheme.{0}} (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
    (x : ∀ p, Fin (N+1) → NumberField.RingOfIntegers (fld p))
    (hcomp : ∀ p, haveI := hnf p; ∃ (i₀ : Fin (N+1))
      (φ : CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
        ⟶ CommRingCat.of (NumberField.RingOfIntegers (fld p))),
      xF p ≫ ψ = Spec.map φ ≫ chartA N ℤ i₀ ∧ x p ≠ 0 ∧
        (∀ k, x p k = φ.hom (projCoord N ℤ i₀ k) * x p i₀) ∧ x p 0 ≠ 0)
    (idx : Fin (N+1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N+1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | haveI := hnf p; htArith (fld p) (E.comap ψ) (xF p) ≤ C}.Finite := by
  refine northcott_hyperplane' N d E hdiv hcont fld hnf hdeg x
    (fun p => haveI := hnf p; htArith (fld p) (E.comap ψ) (xF p)) (fun p => ?_) idx hinj C
  haveI := hnf p
  obtain ⟨i₀, φ, hfac, hx, hw, h0⟩ := hcomp p
  exact ⟨i₀, φ, by rw [htArith_comap, hfac], hx, hw, h0⟩

/-! ## ★出典の紐付け(`.src`) -/

def northcott_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Northcott は射影埋め込みに沿って移送される)",
    sectionId := "genell-prop-1-4" }

def northcott_comap.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_hyperplane(ℙᴺ の Northcott、段 C2e、§9-878)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_hyperplane") 2,
    .citation "[ABC3]" "htArith_comap(高さの関手性、段 C2f、§9-879)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_comap") 2,
    .implicitStep
      ("★hinj が「X の点を分ける」条件になっているのが要点で、" ++
       "これは ψ が**埋め込み**であることの帰結である(§9-849・§9-851)") 3,
    .implicitStep
      ("★★一般の X への還元に残るのは **ψ^*(超平面) = n·D**" ++
       "(段 E3、L^n が非常に豊富)だけである。" ++
       "§9-880 の ht_{D^n} = n·ht_D と本ファイルを繋げば ht_D についての Northcott になる") 5 ]

end ABC3.Found.GenEll
