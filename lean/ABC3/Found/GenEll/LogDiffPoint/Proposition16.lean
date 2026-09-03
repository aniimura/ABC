/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogDiffTower
import ABC3.Found.GenEll.MinField
import ABC3.Found.GenEll.UPoint
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.LogDiffPoint.Definition15

/-!
# LogDiffPoint —— `[GenEll] Proposition 1.6` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

@[simp] theorem logCondAt_algPointOff (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} {D : X.IdealSheafData} (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D xF ≠ 0) :
    logCondAt (algPointOff F xF h) = logCond F D xF := rfl

/-- ★★**負の対照** —— `D = ⊤`（空因子）なら `log-cond = 0`。

★引き戻しの向きか被約化を取り違えていれば、ここが `0` にならない。 -/
@[simp] theorem logCondAt_top {X : Scheme.{0}} (p : AlgPointOff X (⊤ : X.IdealSheafData)) :
    logCondAt p = 0 := by
  letI := p.instField
  letI := p.instNF
  exact logCond_top p.fld p.map

/-- ★★**[GenEll] Proposition 1.6 の非アルキメデス側**（点の言葉で）。

> `log-cond_D(x) ≤ deg_F(D_x)`

★被約化で次数は減る（`deg_idealADiv_radical_le`）。 -/
theorem logCondAt_le {X : Scheme.{0}} {D : X.IdealSheafData} (p : AlgPointOff X D) :
    logCondAt p
      ≤ (letI := p.instField; letI := p.instNF;
         degNormalized (idealADiv p.fld (pullbackIdeal p.fld D p.map))) := by
  letI := p.instField
  letI := p.instNF
  exact logCond_le_degNormalized_pullback p.fld D p.map p.off

def logCondAt_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(非アルキメデス側——点の言葉で)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
