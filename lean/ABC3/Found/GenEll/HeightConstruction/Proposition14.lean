/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArchPoint
import ABC3.Found.GenEll.HeightConstruction.Definition12

/-!
# HeightConstruction —— `[GenEll] Proposition 1.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
variable (F : Type) [Field F] [NumberField F]

/-- ★**空の算術因子の高さは 0**(非空虚性の witness)。 -/
@[simp] theorem htArith_one {X : Scheme.{0}} (xF : specRingOfIntegers F ⟶ X) :
    htArith F (ArithCartier.one X) xF = 0 := by
  rw [htArith, pullbackADiv_one, degNormalized_zero]

/-- ★★★**引き戻しの加法性** —— 残る 1 本を仮説に置いた形。

★★アルキメデス側は**無条件**、有限素点側は `PullbackMul` と非零性のみに依る。 -/
theorem pullbackADiv_tensor {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hmul : PullbackMul F xF)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    pullbackADiv F (D.tensor E) xF = pullbackADiv F D xF + pullbackADiv F E xF := by
  refine Prod.ext ?_ ?_
  · show (pullbackADiv F (D.tensor E) xF).fin
      = (pullbackADiv F D xF + pullbackADiv F E xF).fin
    rw [ADiv.fin_add, pullbackADiv_fin, pullbackADiv_fin, pullbackADiv_fin]
    show (idealADiv F (pullbackIdeal F (D.divisor * E.divisor) xF)).fin = _
    rw [hmul D.divisor E.divisor, idealADiv_mul F _ _ hD hE, ADiv.fin_add]
  · show (pullbackADiv F (D.tensor E) xF).arc
      = (pullbackADiv F D xF + pullbackADiv F E xF).arc
    rw [ADiv.arc_add]
    exact pullbackADiv_arc_tensor F D E xF

/-- ★★★**`Proposition 1.4, (i)`** —— 構成した高さの加法性。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**posit された `HeightTheoryData` の上ではこれは偽だった**
(`Check/GenEll/HeightAxiomGap.lean` に反例)。
構成に置き換えると**真になる**——それが本ファイルの主眼である。 -/
theorem htArith_tensor {X : Scheme.{0}} (D E : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hmul : PullbackMul F xF)
    (hD : pullbackIdeal F D.divisor xF ≠ 0) (hE : pullbackIdeal F E.divisor xF ≠ 0) :
    htArith F (D.tensor E) xF = htArith F D xF + htArith F E xF := by
  rw [htArith, htArith, htArith, pullbackADiv_tensor F D E xF hmul hD hE,
    degNormalized_add]

def htArith_tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (i)(構成した高さの加法性——引き戻しの積保存を仮説に置く)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
