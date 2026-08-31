/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottComap
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★[GenEll] Proposition 1.4, (iv) —— 非常に豊富な因子について（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★これは何か —— 鎖の閉じ方

★`§9-881`（段 C2h）: 射 `ψ : X ⟶ ℙᴺ` に沿って Northcott が移送される
★★`§9-880`（段 C2g）: `ht_{D^n} = n·ht_D`

★★★この 2 本を繋ぐと、**`D` の `n` 倍が（計量込みで）超平面の引き戻しであるような
`D` については、`ht_D` で有界な点の集合が有限である**。

    `D^{⊗n} = ψ^*Ē`  ⟹  `{p | ht_D(x_p) ≤ C}` は有限

★★★★これが原文の
「`L_ℚ` が豊富 ⟹ ある冪が非常に豊富 ⟹ 射影埋め込み ⟹ Northcott」
の**最後の 2 段**（射影埋め込みから先）である。

## ★★★仮定 `hDn` の読み方

★`hDn : D.npow n = E.comap ψ` は「`D` の `n` 倍が、超平面の算術因子 `Ē` の
`ψ` に沿った引き戻しに（**計量込みで**）等しい」という条件である。
★★因子だけの条件（`IsVeryAmpleDiv`）より強いが、
計量の差は `Proposition 1.4, (iii)` で吸収できる（`§9-872`）ので、
BD-類の水準では同じことである。

## ★残っている段（明示）

★★★**`hDn` を満たす `n` の存在**（＝ Serre の定理「豊富なら或る冪は非常に豊富」）が残る。
★これは段 E3（`Found/GenEll/GlobalChartToProj.lean` の
`isImmersion_globalChartToProj` まで来ている）の締めくくりである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-- ★★★★★★★★★★★★★★★**[GenEll] Proposition 1.4, (iv)**——
非常に豊富な因子について（射影埋め込みから先）。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-880`（`ht_{D^n} = n·ht_D`）と `§9-881`（移送）を繋いだものである。 -/
theorem northcott_of_veryAmple (N d n : ℕ) (hn : 0 < n)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {X : Scheme.{0}} (ψ : X ⟶ Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ))
    (D : ArithCartier X) (hDn : D.npow n = E.comap ψ)
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (xF : ∀ p, haveI := hnf p; specRingOfIntegers (fld p) ⟶ X)
    (hDne : ∀ p, haveI := hnf p; pullbackIdeal (fld p) D.divisor (xF p) ≠ 0)
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
    {p : P | haveI := hnf p; htArith (fld p) D (xF p) ≤ C}.Finite := by
  refine Set.Finite.subset
    (northcott_comap N d E hdiv hcont ψ fld hnf hdeg xF x hcomp idx hinj ((n : ℝ) * C))
    (fun p hp => ?_)
  haveI := hnf p
  have hp' : htArith (fld p) D (xF p) ≤ C := hp
  have hkey : htArith (fld p) (E.comap ψ) (xF p) = (n : ℝ) * htArith (fld p) D (xF p) := by
    rw [← hDn, htArith_npow (fld p) D (xF p) (hDne p) n]
  show htArith (fld p) (E.comap ψ) (xF p) ≤ (n : ℝ) * C
  rw [hkey]
  exact mul_le_mul_of_nonneg_left hp' (Nat.cast_nonneg n)

/-! ## ★出典の紐付け(`.src`) -/

def northcott_of_veryAmple.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(非常に豊富な因子についての Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_of_veryAmple.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_comap(射に沿った移送、段 C2h、§9-881)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_comap") 2,
    .citation "[ABC3]" "htArith_npow(ht_{D^n} = n·ht_D、段 C2g、§9-880)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_npow") 2,
    .implicitStep
      ("★仮定 hDn : D.npow n = E.comap ψ は「D の n 倍が、超平面の算術因子 Ē の " ++
       "ψ に沿った引き戻しに(**計量込みで**)等しい」という条件である。" ++
       "★★因子だけの条件(IsVeryAmpleDiv)より強いが、計量の差は Proposition 1.4, (iii) で" ++
       "吸収できる(§9-872)ので、BD-類の水準では同じことである") 4,
    .implicitStep
      ("★★★**hDn を満たす n の存在**(＝ Serre の定理「豊富なら或る冪は非常に豊富」)が残る。" ++
       "これは段 E3(GlobalChartToProj.lean の isImmersion_globalChartToProj まで来ている)の" ++
       "締めくくりである") 6 ]

end ABC3.Found.GenEll
