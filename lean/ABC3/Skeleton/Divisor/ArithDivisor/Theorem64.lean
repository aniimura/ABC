/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.ProductFormula
import Mathlib.NumberTheory.NumberField.Units.DirichletTheorem
import ABC3.Skeleton.Divisor.ArithDivisor.Example63

/-!
# ArithDivisor —— `[FrdI] Theorem 6.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Skeleton.Divisor

open ABC3.Meta NumberField
universe u
variable (L : Type u) [Field L] [NumberField L]

/-- ★★★**`δ_A : Pic_Φ(A) ≅ ℝ`** —— `Theorem 6.4, (i)` の最後。

原文 (FrdI p.115):
> support whose image under degarith

★★原文は「an immediate consequence of the well-known **Dirichlet unit theorem**」と書く。
mathlib に在る(`NumberField.Units.dirichletUnitTheorem`)。 -/
theorem degArith_surjective_and_kernel_eq_image :
    Function.Surjective (degArith L) := by
  sorry

def degArith_surjective_and_kernel_eq_image.src : Source :=
  { paper := "FrdI", pdfPage := 114, item := "Theorem 6.4, (i) — δ_A : Pic_Φ(A) ≅ ℝ",
    sectionId := "frdi-thm-6-4" }

/-- ★★原文が「well-known」で畳んだ所。mathlib に在る。 -/
def degArith_surjective_and_kernel_eq_image.needs : List ProofObligation :=
  [ .citation "[mathlib]" "NumberField.Units.dirichletUnitTheorem(Dirichlet 単数定理)"
      (.inMathlib "NumberField.Units.dirichletUnitTheorem") 114,
    .derivation "Φ^birat(L) ⊗ ℝ の像が deg^arith の核に一致する(単数定理の階数の主張)" 114,
    .implicitStep "★原文は「an immediate consequence of the well-known Dirichlet unit theorem」で畳む" 114 ]

end ABC3.Skeleton.Divisor
