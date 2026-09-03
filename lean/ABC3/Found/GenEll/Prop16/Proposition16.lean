/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArcModel
import ABC3.Found.GenEll.UPoint

/-!
# Prop16 —— `[GenEll] Proposition 1.6` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
open scoped LinearAlgebra.Projectivization
variable (F : Type) [Field F] [NumberField F]
variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-! ## ★★★`Proposition 1.6` -/

/-- ★★★**[GenEll] Proposition 1.6**(Conductor Bounded by the Height)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

    `log-cond_D(x) ≤ ht_D(x) + C`   (`C` は `F` にも `x` にも依らない)

★★★**射影モデルと Green 関数の連続性から、定数が自動で出る。**

★機構は 3 段:
1. `X^arc` はコンパクト(`ArcModel.compactSpace`)
2. 連続な Green 関数は下に有界(`ArcModel.exists_bound`)
3. 導手は高さ + 定数で抑えられる(`ArchBound.logCond_le_htArith_add`)

★★定数が `F` に依らないのは**正規化のおかげ**である
(`Σ_v mult v = [F:ℚ]`。`ArchBound.lean` の測定)。 -/
theorem prop_1_6 (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D : ArithCartier X) (hg : @Continuous _ _ M.topology _ D.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        pullbackIdeal F D.divisor xF ≠ 0 →
        logCond F D.divisor xF ≤ htArith F D xF + C := by
  obtain ⟨C, hC, hlo, -⟩ := M.exists_bound D.green hg
  exact ⟨C, hC, fun F _ _ xF h => logCond_le_htArith_add F D xF C hC h hlo⟩

def prop_1_6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(射影モデルを与えられたものとして——ℤ-固有からの構成は含まない)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
