/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightNormalizationBridge
import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★[GenEll] Proposition 1.4, (iv) —— `ℙᴺ` の超平面因子について（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★これは何か —— Serre の道の到着点

★`§9-876`（段 C2c）：`|ht_E(x) − log H(x)/[F:ℚ]| ≤ C`
★★`§9-877`（段 C2d）：正規化の橋 ⟹ `mulHeight ≤ exp(d·ht + const)`（`hcmp` の形）
★★★`NorthcottCoord.lean` の `northcott_of_projModel`：`hcmp` ＋ `hinj` ⟹ 有限

★★★★これを繋ぐと、**`ℙᴺ_ℤ` の超平面因子については Northcott 性が
（計量の連続性以外の仮定なしで）出る**。

## ★★★仮定の読み方

| 仮定 | 意味 |
|---|---|
| `hdiv`・`hcont` | `E` は超平面因子で、計量は Fubini–Study に対して連続（`Definition 1.1`） |
| `hdeg` | 次数が `d` 以下（原文の `X(ℚ)_{≤d}`） |
| `hht` | 各点が**整な同次座標**を持ち、その高さが `ht` である |
| `hinj` | 正規化座標が点を分ける（`§9-849`・`§9-851`） |

★`hht` は `§9-854`（`exists_integral_repr`）と `§9-C2b`（点はチャートに入る）から作る。

## ★測定の記録

★★★★★結論の高さが **`d·ht`** になっているのは `§9-877` の測定
（`Height.mulHeight` は相対、`htArith` は絶対）による。
★`{p | ht p ≤ C} ⊆ {p | d·ht p ≤ d·C}`（`d ≥ 0`）なので、
**絶対高さで有界な点の集合の有限性はここから出る**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField MvPolynomial HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-- ★★★★★★★★★★★★★★★★**[GenEll] Proposition 1.4, (iv)**——`ℙᴺ` の超平面因子について。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★Serre の道（段 A〜段 F）の到着点である。 -/
theorem northcott_hyperplane (N d : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (x : ∀ p, Fin (N+1) → NumberField.RingOfIntegers (fld p))
    (ht : P → ℝ)
    (hht : ∀ p, haveI := hnf p; ∃ (i₀ : Fin (N+1))
      (ψ : CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
        ⟶ CommRingCat.of (NumberField.RingOfIntegers (fld p))),
      ht p = htArith (fld p) E (Spec.map ψ ≫ chartA N ℤ i₀) ∧ x p ≠ 0 ∧
        (∀ k, x p k = ψ.hom (projCoord N ℤ i₀ k) * x p i₀) ∧ x p 0 ≠ 0)
    (idx : Fin (N+1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N+1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | (d : ℝ) * ht p ≤ C}.Finite := by
  obtain ⟨const, hconst, hbound⟩ := mulHeight_le_exp_htArith N d E hdiv hcont
  refine northcott_of_projModel d (fun p => (d : ℝ) * ht p) fld hnf hdeg
    (fun p k => ((x p k : fld p))) idx const (fun p => ?_) hinj C
  haveI := hnf p
  obtain ⟨i₀, ψ, hteq, hx, hw, h0⟩ := hht p
  have hd : (Module.finrank ℚ (fld p) : ℝ) ≤ (d : ℝ) := by exact_mod_cast hdeg p
  have := hbound (fld p) hd i₀ ψ (x p) hx hw h0
  rw [hteq]
  exact this

/-- ★★**絶対高さで有界な点も有限である** —— `d ≥ 0` なので包含から出る。 -/
theorem northcott_hyperplane' (N d : ℕ)
    (E : ArithCartier (Proj (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ)))
    (hdiv : E.divisor = hyperplaneIdeal N ℤ)
    (hcont : @Continuous _ _ (projArcModel N).topology _
      (fun p => E.green p - greenFS N p))
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (x : ∀ p, Fin (N+1) → NumberField.RingOfIntegers (fld p))
    (ht : P → ℝ)
    (hht : ∀ p, haveI := hnf p; ∃ (i₀ : Fin (N+1))
      (ψ : CommRingCat.of (HomogeneousLocalization.Away
          (MvPolynomial.homogeneousSubmodule (Fin (N+1)) ℤ) (MvPolynomial.X i₀))
        ⟶ CommRingCat.of (NumberField.RingOfIntegers (fld p))),
      ht p = htArith (fld p) E (Spec.map ψ ≫ chartA N ℤ i₀) ∧ x p ≠ 0 ∧
        (∀ k, x p k = ψ.hom (projCoord N ℤ i₀ k) * x p i₀) ∧ x p 0 ≠ 0)
    (idx : Fin (N+1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N+1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | ht p ≤ C}.Finite := by
  refine Set.Finite.subset
    (northcott_hyperplane N d E hdiv hcont fld hnf hdeg x ht hht idx hinj ((d : ℝ) * C))
    (fun p hp => ?_)
  have hp' : ht p ≤ C := hp
  exact mul_le_mul_of_nonneg_left hp' (Nat.cast_nonneg d)

/-! ## ★出典の紐付け(`.src`) -/

def northcott_hyperplane.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ℙᴺ の超平面因子についての Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_hyperplane'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(絶対高さで有界な点も有限)",
    sectionId := "genell-prop-1-4" }

def northcott_hyperplane.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_of_projModel(hcmp ＋ hinj ⟹ 有限)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_of_projModel") 2,
    .citation "[ABC3]" "mulHeight_le_exp_htArith(hcmp の形、段 C2d、§9-877)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mulHeight_le_exp_htArith") 2,
    .citation "[ABC3]" "abs_htArith_sub_log_mulHeight_le(段 C2c、§9-876)"
      (.inProject "ABC3" "ABC3.Found.GenEll.abs_htArith_sub_log_mulHeight_le") 2,
    .implicitStep
      ("★仮定 hht(各点が**整な同次座標**を持ち、その高さが ht である)は " ++
       "§9-854(exists_integral_repr)と §9-C2b(点はチャートに入る)から作る") 4,
    .implicitStep
      ("★★仮定 hinj(正規化座標が点を分ける)は §9-849・§9-851 から出る" ++
       "——ただし本補題では消費側のデータとして受けている") 4,
    .implicitStep
      ("★★★これは **ℙᴺ の超平面因子**についての Northcott であり、" ++
       "原文の「X が ℤ-固有で L_ℚ が豊富」の場合そのものではない。" ++
       "★一般の X への還元は段 E3(ample ⟹ 射影埋め込み)が担う") 5 ]

end ABC3.Found.GenEll
