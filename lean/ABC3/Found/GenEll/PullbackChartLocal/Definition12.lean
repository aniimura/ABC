/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottGlobalToProj
import ABC3.Found.GenEll.PullbackLocalization
import ABC3.Meta.Claim

/-!
# PullbackChartLocal —— `[GenEll] Definition 1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MvPolynomial HomogeneousLocalization NumberField
open ABC3.Found.Arakelov
variable {X : Scheme.{0}}

/-! ## ★★★合成に沿った引き戻し -/

/-- ★★★**`comap` に沿った引き戻しは合成の引き戻しである**。

★`Scheme.IdealSheafData.comap_comp` をそのまま `pullbackIdealOf` の定義に流すだけ。 -/
theorem pullbackIdealOf_comap {B : CommRingCat.{0}} {X Y : Scheme.{0}}
    (D : Y.IdealSheafData) (f : X ⟶ Y) (x : Spec B ⟶ X) :
    pullbackIdealOf B (D.comap f) x = pullbackIdealOf B D (x ≫ f) := by
  rw [pullbackIdealOf, pullbackIdealOf, Scheme.IdealSheafData.comap_comp]

/-! ## ★出典の紐付け(`.src`) -/

def pullbackIdealOf_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(comap に沿った引き戻しは合成の引き戻しである)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
