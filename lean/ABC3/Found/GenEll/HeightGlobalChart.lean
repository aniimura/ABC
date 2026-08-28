/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HyperplanePullbackGlobal
import ABC3.Found.GenEll.GreenGlobal
import ABC3.Found.GenEll.ArithCartierComap
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★`X` の点の高さは比の同次座標の素朴高さである（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★★これは何か —— 一般の `X` での高さの明示公式

★`§9-871`: `ℙᴺ` では `htArith (hyperplaneArith N) = log H(x)/[F:ℚ]`
★★`§9-879`: 高さは射に沿って関手的（`ht_{ψ^*E}(x) = ht_E(ψ ∘ x)`）
★★★`§9-886`: チャート射に沿った点は `Spec.map ρ ≫ chartA i` に分解する

★★★★これを繋ぐと、**`X` の点の高さが比の同次座標の素朴高さで書ける**:

    `htArith F ( (hyperplaneArith N)^* (globalChartToProj i) ) x
       = log H( (s_k/s_i)(x) )_k / [F : ℚ]`

★これが「豊富な直線束の高さは、切断で作った射影座標の素朴高さである」ことの中身である。

## ★★★機構 —— 3 本の合成だけ

★`htArith_comap`（`§9-879`）で `ℙᴺ` の点の高さに移し、
`globalChart_factor`（本ファイル）でチャート射を `Spec.map ρ ≫ chartA i` に分解し、
`htArith_hyperplaneArith`（`§9-871`）を当てる。

## ★残っている段（明示）

★★仮定 `hw`（点の同次座標が `ρ(x_k/x_i)` の分母を払った整数表示である）は
`§9-854`（`exists_integral_repr`）から作る。
★★★これで `northcott_comap`（`§9-881`）へ渡せば、
**`X` の点についての Northcott** になる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★チャート射に沿った点の分解 -/

/-- ★**チャート射に沿った点は `Spec.map ρ ≫ chartA i` に分解する**。 -/
theorem globalChart_factor (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme) :
    xF ≫ globalChartToProj M hM φ s i
      = Spec.map (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)) ≫ chartA N ℤ i := by
  rw [Spec.map_preimage, globalChartToProj, Category.assoc]
  rfl

/-! ## ★★★★★★★★★★★★★★高さの明示公式 -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★★★★★**`X` の点の高さは比の同次座標の素朴高さである**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `htArith F ( (hyperplaneArith N)^* (globalChartToProj i) ) x = log H(x_k)/[F:ℚ]`

★これが「豊富な直線束の高さは、切断で作った射影座標の素朴高さである」ことの中身である。
★★機構は `§9-879`（関手性）・本ファイルの分解・`§9-871`（`ℙᴺ` での等式）の合成だけ。 -/
theorem htArith_comap_globalChart (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1))
    (F : Type) [Field F] [NumberField F]
    (xF : specRingOfIntegers F ⟶ (nonVanishing M (s i)).toScheme)
    (x : Fin (N+1) → NumberField.RingOfIntegers F) (hx : x ≠ 0)
    (hw : ∀ k, x k = (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)).hom
      (projCoord N ℤ i k) * x i)
    (h0 : x 0 ≠ 0) :
    htArith F ((hyperplaneArith N).comap (globalChartToProj M hM φ s i)) xF
      = Real.log (Height.mulHeight (fun k => ((x k : F)))) / (Module.finrank ℚ F : ℝ) := by
  rw [htArith_comap]
  rw [show xF ≫ globalChartToProj M hM φ s i
      = Spec.map (Spec.preimage (xF ≫ globalChartMorphism M hM φ s i)) ≫ chartA N ℤ i from
    globalChart_factor M hM φ s i F xF]
  exact htArith_hyperplaneArith F N i _ x hx hw h0

/-! ## ★出典の紐付け(`.src`) -/

def globalChart_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャート射に沿った点の分解)",
    sectionId := "genell-prop-1-4" }

def htArith_comap_globalChart.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(X の点の高さは比の同次座標の素朴高さである)",
    sectionId := "genell-prop-1-4" }

def htArith_comap_globalChart.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htArith_hyperplaneArith(ℙᴺ での等式、§9-871)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_hyperplaneArith") 2,
    .citation "[ABC3]" "htArith_comap(高さの関手性、段 C2f、§9-879)"
      (.inProject "ABC3" "ABC3.Found.GenEll.htArith_comap") 2,
    .implicitStep
      ("★これが「豊富な直線束の高さは、切断で作った射影座標の素朴高さである」ことの中身である") 3,
    .implicitStep
      ("★★仮定 hw(点の同次座標が ρ(x_k/x_i) の分母を払った整数表示である)は " ++
       "§9-854(exists_integral_repr)から作る。" ++
       "★★★これで northcott_comap(§9-881)へ渡せば **X の点についての Northcott** になる") 4 ]

end ABC3.Found.GenEll
