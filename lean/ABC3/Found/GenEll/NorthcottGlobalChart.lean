/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightGlobalChart
import ABC3.Found.GenEll.NorthcottComap
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★[GenEll] Proposition 1.4, (iv) —— `X` のチャートの上で（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★これは何か —— `X` の点についての Northcott

★`§9-881`（段 C2h）: Northcott は射 `ψ : X ⟶ ℙᴺ` に沿って移送される
★★`§9-887`（段 E3g）: チャート射に沿った点は `Spec.map ρ ≫ chartA i` に分解する

★★★これを繋ぐと、**切断の非零開 `X_{s_i}` の点について Northcott が出る**:

    `{ p | ht( (超平面)^*(globalChartToProj i) )(x_p) ≤ C }` は有限

★★★★計量については**仮定が要らない**——`hyperplaneArith N` の Green 関数は
`greenFS` そのものなので、差は `0`（定数関数）で連続だからである。

## ★★★残っている仮定（明示）

| 仮定 | 出どころ |
|---|---|
| `hdeg` | 原文の `X(ℚ)_{≤d}` |
| `hw`・`hx`・`h0` | `§9-854`（`exists_integral_repr`）——点の同次座標の整数表示 |
| `hinj` | `§9-849`・`§9-851`（正規化座標が点を分ける） |

★★★★★これが `Proposition 1.4, (iv)` の**チャート版**である。
`X` 全体へ広げるには `X = ⋃ X_{s_i}`（段 E2）で有限和を取ればよい。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★計量の仮定は要らない -/

/-- ★**`hyperplaneArith N` と Fubini–Study の差は `0` である**（したがって連続）。 -/
theorem continuous_green_hyperplaneArith_sub (N : ℕ) :
    @Continuous _ _ (projArcModel N).topology _
      (fun p => (hyperplaneArith N).green p - greenFS N p) := by
  show @Continuous _ _ (projArcModel N).topology _ (fun p => greenFS N p - greenFS N p)
  simpa using (@continuous_const _ _ (projArcModel N).topology _ (0:ℝ))

/-! ## ★★★★★★★★★★★★★★★★チャートの上の Northcott -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★★★**[GenEll] Proposition 1.4, (iv)**——`X` のチャートの上で。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`§9-881`（移送）と `§9-887`（チャート射の分解）の合成である。
★★計量については**仮定が要らない**——差が `0`（定数関数）だからである。 -/
theorem northcott_globalChart (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) (d : ℕ)
    {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (xF : ∀ p, haveI := hnf p;
      specRingOfIntegers (fld p) ⟶ (nonVanishing M (s i)).toScheme)
    (x : ∀ p, Fin (N+1) → NumberField.RingOfIntegers (fld p))
    (hx : ∀ p, x p ≠ 0)
    (hw : ∀ p, haveI := hnf p; ∀ k, x p k
      = (Spec.preimage (xF p ≫ globalChartMorphism M hM φ s i)).hom
          (projCoord N ℤ i k) * x p i)
    (h0 : ∀ p, x p 0 ≠ 0)
    (idx : Fin (N+1))
    (hinj : Function.Injective (fun (p : P) (k : Fin (N+1)) =>
      ((((x p k : fld p)) / ((x p idx : fld p)) : fld p) : ℂ)))
    (C : ℝ) :
    {p : P | haveI := hnf p;
      htArith (fld p) ((hyperplaneArith N).comap (globalChartToProj M hM φ s i))
        (xF p) ≤ C}.Finite :=
  northcott_comap N d (hyperplaneArith N) rfl (continuous_green_hyperplaneArith_sub N)
    (globalChartToProj M hM φ s i) fld hnf hdeg xF x
    (fun p => haveI := hnf p;
      ⟨i, Spec.preimage (xF p ≫ globalChartMorphism M hM φ s i),
        globalChart_factor M hM φ s i (fld p) (xF p), hx p, hw p, h0 p⟩)
    idx hinj C

/-! ## ★出典の紐付け(`.src`) -/

def continuous_green_hyperplaneArith_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(hyperplaneArith と Fubini–Study の差は 0)",
    sectionId := "genell-prop-1-4" }

def northcott_globalChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(X のチャートの上での Northcott)",
    sectionId := "genell-prop-1-4" }

def northcott_globalChart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "northcott_comap(射に沿った移送、段 C2h、§9-881)"
      (.inProject "ABC3" "ABC3.Found.GenEll.northcott_comap") 2,
    .citation "[ABC3]" "globalChart_factor(チャート射に沿った点の分解、段 E3g、§9-887)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalChart_factor") 2,
    .implicitStep
      ("★計量については**仮定が要らない**——hyperplaneArith N の Green 関数は " ++
       "greenFS そのものなので、差は 0(定数関数)で連続だからである") 2,
    .implicitStep
      ("★★残っている仮定は hdeg(原文の X(ℚ)_{≤d})・hw/hx/h0(§9-854 の整数表示)・" ++
       "hinj(§9-849・§9-851 の点の分離)の 3 つで、いずれも出どころが判っている") 3,
    .implicitStep
      ("★★★これが Proposition 1.4, (iv) の**チャート版**である。" ++
       "X 全体へ広げるには X = ⋃ X_{s_i}(段 E2)で有限和を取ればよい") 4 ]

end ABC3.Found.GenEll
