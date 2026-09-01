/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.TateSigmaGalAct
import ABC3.Meta.Claim

/-!
# 第 1288 ブロック —— **`σ` の作用は変数変換と可換**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——節点 8 の要

☆`E ⊗ K` と Tate モデルは**変数変換 `C` で結ばれる**（`exists_tate_model`・
`tateModel_baseChange`、どちらも在庫）。

★`C` の係数は基礎局所体の整数環にあるので `σ` で固定される。
したがって `σ` の作用は変数変換と可換である——本ブロックがそれである。

☆座標では `σ(vcX C x) = vcX C (σ x)`（第 1284）であり、
点はその座標で決まる（第 1286 と同じ道具）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {K : Type} [Field K]

/-- ★★★★**`galAct` は原点でない点を原点でない点に送る**——★**無条件**（第 1288）。 -/
theorem galAct_ne_zero (σK : K →+* K) (W : WeierstrassCurve K) (hW : W.map σK = W)
    {P : W.toAffine.Point} (hP : P ≠ 0) : galAct σK W hW P ≠ 0 := by
  intro hcon
  refine hP ?_
  have h1 : rhPoint σK W P = 0 := by
    have h2 := congrArg (pointEquivOfEq hW).symm hcon
    simpa [galAct, rhPointHom] using h2
  exact (rhPoint_eq_zero_iff σK W P).1 h1

/-- ★★★★**`vcPoint` は原点でない点を原点でない点に送る**——★**無条件**（第 1288）。 -/
theorem vcPoint_ne_zero (C : VariableChange K) (W : WeierstrassCurve K)
    {P : W.toAffine.Point} (hP : P ≠ 0) : vcPoint C W P ≠ 0 := by
  cases P with
  | zero => exact absurd rfl hP
  | some x y h =>
      intro hcon
      rw [vcPoint_some] at hcon
      exact absurd hcon (by simp)

/-- ★★★★**原点でなければ `vcPoint` の座標は `vcX`・`vcY`**——★**無条件**（第 1288）。 -/
theorem pointCoords_vcPoint' (C : VariableChange K) (W : WeierstrassCurve K)
    (P : W.toAffine.Point) (hP : P ≠ 0) :
    pointCoords (vcPoint C W P)
      = (vcX C (pointCoords P).1, vcY C (pointCoords P).1 (pointCoords P).2) := by
  cases P with
  | zero => exact absurd rfl hP
  | some x y h => rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**`σ` の作用は変数変換と可換**——★**無条件**（第 1288）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`σ` が変数変換 `C` の係数を固定し、曲線 `W` も固定するなら、
`vcPoint C W` は `σ` の作用と可換である。 -/
theorem galAct_vcPoint (σK : K →+* K) (W : WeierstrassCurve K) (C : VariableChange K)
    (hC : C.map σK = C) (hW : W.map σK = W)
    (P : W.toAffine.Point) :
    galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW])
        (vcPoint C W P)
      = vcPoint C W (galAct σK W hW P) := by
  by_cases hP : P = 0
  · subst hP
    rw [vcPoint_zero, map_zero, map_zero, vcPoint_zero]
  · refine pointCoords_injective
      (galAct_ne_zero σK (C • W) _ (vcPoint_ne_zero C W hP))
      (vcPoint_ne_zero C W (galAct_ne_zero σK W hW hP)) ?_
    rw [pointCoords_galAct, pointCoords_vcPoint' C W P hP,
      pointCoords_vcPoint' C W _ (galAct_ne_zero σK W hW hP), pointCoords_galAct]
    simp only []
    rw [map_vcX_fixed σK C hC, map_vcY_fixed σK C hC]

/-! ## ★出典の紐付け(`.src`) -/

def galAct_ne_zero.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galAct は原点でない点を原点でない点に送る。★無条件)",
    sectionId := "genell-thm-3-8" }

def vcPoint_ne_zero.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(vcPoint は原点でない点を原点でない点に送る。★無条件)",
    sectionId := "genell-thm-3-8" }

def pointCoords_vcPoint'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(原点でなければ vcPoint の座標は vcX・vcY。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct_vcPoint.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ の作用は変数変換と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct_vcPoint.needs : List ProofObligation :=
  [ .citation "[ABC3]" "map_vcX_fixed・map_vcY_fixed(第 1284、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.map_vcX_fixed") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1288）**——節点 8 の要である。" ++
       "☆`E ⊗ K` と Tate モデルは変数変換 `C` で結ばれ（`tateModel_baseChange`、在庫）、" ++
       "`C` の係数は基礎局所体の整数環にあるので `σ` で固定される。" ++
       "★したがって `σ` の作用は両者の間で対応する。") 2 ]

end ABC3.Found.GenEll
