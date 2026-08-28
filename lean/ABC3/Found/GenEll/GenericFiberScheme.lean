/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GenericFiber
import ABC3.Found.GenEll.NorthcottVeryAmple
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 F2c —— 閉包のイデアル層を生成ファイバーへ落とす（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★★★これは何か —— 段 F2c の締めくくり

`§9-GenericFiber` は**アフィンでの核**を取った:

    `ker (A_M → B) = (ker φ) · A_M`   （`ker_locMap`）

★台帳は「`Scheme.IdealSheafData` の水準へ持ち上げるには
`Scheme.Hom.ker_apply` を各チャートに当てて貼るだけで、**新しい数学は要らない**」
と測っていた。★★本ファイルはその持ち上げである。

    `ker (locMap) = (f.ker.ideal U) · A_M`
    `A_M / (f.ker.ideal U)·A_M ≅ Γ(Y, f⁻¹U)`（`locMap` が全射なら）

★★★すなわち「**閉包を生成ファイバーへ落とすと `Y` そのものに戻る**」ことが、
イデアル層の言葉で書けた。

## ★★★機構 —— 在庫 2 本の合成

| 道具 | 役割 |
|---|---|
| `Scheme.Hom.ker_apply`（mathlib） | `f.ker.ideal U = ker (f.app U)`（`U` アフィン開・`f` 擬コンパクト） |
| `ker_locMap`（`§9-GenericFiber`） | 局所化は核と可換 |

★測定どおり、**新しい数学は要らなかった**（2026-08-28 実測）。

## ★配管の記録

★★`Γ(X, U)` を含む式では `rw` が「型が instances 透明度で正しくない」で落ちるので、
`set_option backward.isDefEq.respectTransparency false` を置く
（`tools/lean-idioms.md` の項）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★**閉包のイデアル層の生成ファイバーへの像**。

    `ker (A_M → Γ(Y, f⁻¹U)) = (f.ker.ideal U) · A_M`

★`Scheme.Hom.ker_apply`（アフィン開で核を読む）と `ker_locMap`（局所化は核と可換）の合成。 -/
theorem ker_ideal_localized {X Y : Scheme.{0}} (f : Y ⟶ X) [QuasiCompact f]
    (U : X.affineOpens) (Aq : Type) [CommRing Aq]
    (M : Submonoid Γ(X, U.1)) [Algebra (Γ(X, U.1)) Aq] [IsLocalization M Aq]
    (hM : ∀ m ∈ M, IsUnit ((f.app U.1).hom m)) :
    RingHom.ker (locMap Aq M (f.app U.1).hom hM)
      = Ideal.map (algebraMap (Γ(X, U.1)) Aq) (f.ker.ideal U) := by
  rw [Scheme.Hom.ker_apply f U, ker_locMap]

set_option backward.isDefEq.respectTransparency false in
/-- ★★★★★★★★★★**段 F2c** —— 閉包を生成ファイバーへ落とすと `Y` に戻る。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

    `A_M / (f.ker.ideal U)·A_M ≅ Γ(Y, f⁻¹U)`

★`locMap` が全射（＝ `Y` が生成ファイバーの側で閉部分スキーム）のときである。 -/
noncomputable def genericFiberEquiv {X Y : Scheme.{0}} (f : Y ⟶ X) [QuasiCompact f]
    (U : X.affineOpens) (Aq : Type) [CommRing Aq]
    (M : Submonoid Γ(X, U.1)) [Algebra (Γ(X, U.1)) Aq] [IsLocalization M Aq]
    (hM : ∀ m ∈ M, IsUnit ((f.app U.1).hom m))
    (hsurj : Function.Surjective (locMap Aq M (f.app U.1).hom hM)) :
    (Aq ⧸ Ideal.map (algebraMap (Γ(X, U.1)) Aq) (f.ker.ideal U))
      ≃+* Γ(Y, f ⁻¹ᵁ U.1) :=
  (Ideal.quotEquivOfEq (by rw [Scheme.Hom.ker_apply f U])).trans
    (quotKerMapEquiv Aq M (f.app U.1).hom hM hsurj)

/-! ## ★出典の紐付け(`.src`) -/

def ker_ideal_localized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(閉包のイデアル層の生成ファイバーへの像)",
    sectionId := "genell-remark-1-4-1" }

def genericFiberEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(段 F2c——閉包を生成ファイバーへ落とすと Y に戻る)",
    sectionId := "genell-remark-1-4-1" }

def genericFiberEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.ker_apply(アフィン開で核を読む)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.ker_apply") 2,
    .citation "[ABC3]" "ker_locMap(局所化は核と可換、§9-GenericFiber)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ker_locMap") 2,
    .implicitStep
      ("★台帳の測定どおり、**新しい数学は要らなかった**——" ++
       "Scheme.Hom.ker_apply と ker_locMap の合成だけである(2026-08-28 実測)") 2,
    .implicitStep
      ("★★配管: Γ(X, U) を含む式では rw が「型が instances 透明度で正しくない」で落ちるので、" ++
       "set_option backward.isDefEq.respectTransparency false を置く" ++
       "(tools/lean-idioms.md の項)") 2 ]

end ABC3.Found.GenEll
