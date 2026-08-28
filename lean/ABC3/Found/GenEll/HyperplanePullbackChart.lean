/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalAwayHom
import ABC3.Found.GenEll.KerHyperplaneChart
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★チャートの上で超平面の引き戻しは比の切断である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★これは何か —— 非常に豊富であることの中身

`§9-882`（段 C2i）は仮定 `D^{⊗n} = ψ^*Ē`（`n` 倍が超平面の引き戻し）を置いた。
★その仮定の**中身**は「チャートの上で超平面の引き戻しが**比の切断で生成される**」ことである。

★★本ファイルはそれを環の水準で取る:

    `(globalAwayHom i)_* ( ker (Away.map hyperplaneHom (x_i)) ) = ( s_0/s_i )`

★★★左辺は `§9-863` により `(x_0/x_i)` の像であり、
`§9-842` の `globalAwayHom_projCoord`（`x_j/x_i ↦ s_j/s_i`）で右辺になる。

## ★★★機構 —— 在庫 2 本の合成

| 道具 | 役割 |
|---|---|
| `ker_awayMap_hyperplane`（`§9-863`、段 C2c-1） | `ker (Away.map hyperplaneHom (x_i)) = (x_0/x_i)` |
| `globalAwayHom_projCoord`（`§9-842`） | `x_j/x_i ↦ s_j/s_i` |

★あとは `Ideal.map_span` だけである。

## ★残っている段（明示）

★★★★これを `Scheme.IdealSheafData` の水準へ持ち上げると
`(globalChartToProj i)^* (超平面) = div(s_0)|_{X_{s_i}}` になる
——`§9-884` の `divisorOfSection_eq_ofIdealTop`（`X_{s_i}` の上では `M` は自明）が受ける。
★それが `§9-882` の仮定 `hDn` の `n = 1` の場合（＝**非常に豊富**）である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization
open ABC3.Found.Arakelov

attribute [local instance] MvPolynomial.gradedAlgebra

variable {X : Scheme.{0}}

/-! ## ★正規化座標の像 -/

/-- ★**正規化座標の生成するイデアルは比の切断の生成するイデアルへ行く**。 -/
theorem map_globalAwayHom_span_projCoord (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (φ : R →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i j : Fin (N + 1)) :
    Ideal.map (globalAwayHom M hM φ s i) (Ideal.span {projCoord N R i j})
      = Ideal.span {globalRatio M hM (s j) (s i)} := by
  rw [Ideal.map_span, Set.image_singleton, globalAwayHom_projCoord]

/-! ## ★★★★★★★★★★★★超平面の引き戻し -/

/-- ★★★★★★★★★★★★**チャートの上で超平面の引き戻しは比の切断である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `(globalAwayHom i)_* ( ker (Away.map hyperplaneHom (x_i)) ) = ( s_0/s_i )`

★これが「非常に豊富」であることの**中身**である
——超平面を引き戻すと、切断 `s_0` の零点因子（`s_i` で自明化して測ったもの）になる。
★★機構は `§9-863`（段 C2c-1）と `§9-842` の合成だけである。 -/
theorem map_globalAwayHom_ker_hyperplane (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    {N : ℕ}
    (φ : ℤ →+* (Γ(X, (⊤ : X.Opens)) : Type))
    (s : Fin (N + 1) → (M.obj (op ⊤) : Type)) (i : Fin (N + 1)) :
    Ideal.map (globalAwayHom M hM φ s i)
        (RingHom.ker (Away.map (hyperplaneHom N ℤ) (MvPolynomial.X i)))
      = Ideal.span {globalRatio M hM (s 0) (s i)} := by
  rw [ker_awayMap_hyperplane, map_globalAwayHom_span_projCoord]

/-! ## ★出典の紐付け(`.src`) -/

def map_globalAwayHom_span_projCoord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(正規化座標の生成するイデアルは比の切断へ行く)",
    sectionId := "genell-prop-1-4" }

def map_globalAwayHom_ker_hyperplane.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの上で超平面の引き戻しは比の切断である)",
    sectionId := "genell-prop-1-4" }

def map_globalAwayHom_ker_hyperplane.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "ker_awayMap_hyperplane(段 C2c-1、§9-863)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_awayMap_hyperplane") 2,
    .citation "[ABC3]" "globalAwayHom_projCoord(x_j/x_i ↦ s_j/s_i、§9-842)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalAwayHom_projCoord") 2,
    .implicitStep
      ("★これが「非常に豊富」であることの**中身**である" ++
       "——超平面を引き戻すと、切断 s_0 の零点因子(s_i で自明化して測ったもの)になる") 3,
    .implicitStep
      ("★★これを Scheme.IdealSheafData の水準へ持ち上げると " ++
       "(globalChartToProj i)^* (超平面) = div(s_0)|_{X_{s_i}} になる" ++
       "——§9-884 の divisorOfSection_eq_ofIdealTop(X_{s_i} の上では M は自明)が受ける。" ++
       "★それが §9-882 の仮定 hDn の n = 1 の場合(＝非常に豊富)である") 4 ]

end ABC3.Found.GenEll
