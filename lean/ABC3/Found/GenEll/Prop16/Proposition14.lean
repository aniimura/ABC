/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArcModel
import ABC3.Found.GenEll.UPoint
import ABC3.Found.GenEll.Prop16.Proposition16

/-!
# Prop16 —— `[GenEll] Proposition 1.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField
open scoped LinearAlgebra.Projectivization
variable (F : Type) [Field F] [NumberField F]
variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-- ★★**BD-class の水準**(原文の `≲` そのもの)。

★★★定数が `F` に依らないので、**すべての数体の上で同時に**成り立つ。 -/
theorem prop_1_6_bdge (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D : ArithCartier X) (hg : @Continuous _ _ M.topology _ D.green)
    (F : Type) [Field F] [NumberField F]
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0) :
    BDge (fun xF : specRingOfIntegers F ⟶ X => logCond F D.divisor xF)
      (fun xF => htArith F D xF) := by
  obtain ⟨C, hC, hlo, -⟩ := M.exists_bound D.green hg
  exact logCond_bdge_htArith_of_bddBelow F D C hC hlo h

/-! ## ★★★`Proposition 1.4, (ii)` の一様形 -/

/-- ★★★**[GenEll] Proposition 1.4, (ii)** —— 高さは下に一様に有界。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**定数 `C` は `F` にも `x` にも依らない**——`Prop 1.6` と同じ機構である。

★`HeightNonneg.lean` の `htArith_nonneg` は `g ≥ 0` を仮定していたが、
本定理は**射影モデルと連続性**だけでよい——それが原文の仮定である。 -/
theorem prop_1_4_ii (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D : ArithCartier X) (hg : @Continuous _ _ M.topology _ D.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        -C ≤ htArith F D xF := by
  obtain ⟨C, hC, hlo, -⟩ := M.exists_bound D.green hg
  refine ⟨C, hC, fun F _ _ xF => ?_⟩
  have h1 : (0:ℝ) ≤ degNormalized (idealADiv F (pullbackIdeal F D.divisor xF)) :=
    degNormalized_nonneg_of_isEffective F _ (idealADiv_isEffective F _)
  have h2 : -C ≤ (archADiv F D.green xF).sum (fun _ r => r)
      / (Module.finrank ℚ F : ℝ) :=
    archADiv_sum_div_finrank_ge F D.green xF C hC hlo
  have h3 : htArith F D xF
      = degNormalized (idealADiv F (pullbackIdeal F D.divisor xF))
        + (archADiv F D.green xF).sum (fun _ r => r) / (Module.finrank ℚ F : ℝ) := by
    rw [htArith, degNormalized, degNormalized, deg_pullbackADiv]
    ring
  rw [h3]; linarith

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文は `X` が ℤ-固有であることから射影埋め込みを得るが、
本ファイルは `ArcModel` として**与えられたものとして受けている**。 -/

def prop_1_4_ii.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(射影モデルを与えられたものとして——一様な定数つき)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
