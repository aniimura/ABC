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
mathlib に在る(`NumberField.Units.dirichletUnitTheorem`)。

## ★★判断 D16（2026-09-06）—— **畳めなかった。記録に留める**

★**橋は届かない**。原典 `Theorem 6.4, (i)` の結論そのものは
`Found/Divisor/Ex63RlfPic.lean` の `rlfDeltaA`（`Pic_Φ(A) ≃+ ℝ`、sorry 0）として
既に閉じているが、本宣言の `degArith` は `Example63.lean` の **`sorry` 本体の `def`** であり、
`Found` 側の `arithDegreeLin` とは型の上でも項の上でも何の関係も持たない。
薄い橋を架けるには `degArith` に本体を与える（= `Example63.lean` の `sorry` を閉じる）ことが要り、
それは本持ち場の外である。

★★★**測定（2026-09-06）——現在の型のままでは、この `sorry` は原理的に閉じない**。
`degArith` の本体は引数に依存しないので `degArith L x = degArith L y` が `rfl` で通る
（すなわち `degArith L` は定数関数に展開される）。
ゆえに `Function.Surjective (degArith L)` からは `(0 : ℝ) = 1` が出て、
現在の環境では**反証可能**である。
★この反証は `degArith` 経由で `sorryAx` に依存するので、**原典についての主張ではなく、
この place-holder についての主張**である。`Example63.lean` の `sorry` が閉じれば消える。
★statement は変えない。

## ★★★**名前が約束した「核 = 像」を足してはならない**

宣言名は「全射**かつ**核 = 像」だが、型にあるのは全射性——**易しい半分**——だけである。
★素直に核の条を足すと**偽になる**。
`Check/FrdI/Thm64PicDegenerate.lean` の `no_nonzero_arch_kernel`（sorry 0）が、
非実現化の `ArithPhiGp` ではアルキメデス方向の直線が丸ごと核に入るのに
主因子の像は可算なので `ℝ` が可算になって矛盾する、を証明している。
条件は「**無限素点が 2 つ以上**」=「**単数の階数 ≥ 1**」
（`two_infinite_places_iff_units_rank_pos`）であり、
★**原文が Dirichlet 単数定理を引く必要が生じるのと同じ条件**である。
`L = ℚ(√2)` がその最小の実例（`arithPic_ker_not_principal_subgroup_qsqrt2`）。

⇒ 原典に忠実な水準は**実現化した `Φ^rlf`**（`Found/Divisor/Ex63RlfPic.lean`）の側であり、
非実現化の `ArithPhiGp` にこの条を足す修復は `False` を作り込む。 -/
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
    .implicitStep "★原文は「an immediate consequence of the well-known Dirichlet unit theorem」で畳む" 114,
    .citation "[ABC3]" "rlfDeltaA(実現化した Φ^rlf の水準での Theorem 6.4, (i)。Pic_Φ(A) ≃+ ℝ、sorry 0)"
      (.inProject "ABC3" "ABC3.Found.Divisor.rlfDeltaA") 114,
    .citation "[ABC3]" "no_nonzero_arch_kernel(非実現化の ArithPhiGp に核の等式を足すと偽になることの証明、sorry 0)"
      (.inProject "ABC3" "ABC3.Check.FrdI.no_nonzero_arch_kernel") 114,
    .implicitStep
      ("★判断 D16(2026-09-06、記録のみ): 非実現化の ArithPhiGp では「核 = 主因子の像」は偽である"
       ++ "(無限素点が 2 つ以上 = 単数の階数 ≥ 1 のとき)。"
       ++ "原典に忠実な水準は実現化した Φ^rlf の側(Found/Divisor/Ex63RlfPic.lean の rlfDeltaA)であり、"
       ++ "statement には核の条を足さない。"
       ++ "★degArith は sorry 本体の def なので定数関数に展開され、"
       ++ "この全射性は現在の環境では反証可能である(原典ではなく place-holder についての主張)") 114 ]

end ABC3.Skeleton.Divisor
