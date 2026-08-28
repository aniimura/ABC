/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17Bookkeeping
import ABC3.Found.Divisor.ArithFundamental
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★幾何の紐の第 1 本 —— ノルムと剰余次数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`§9-971` で `Proposition 1.7, (i)` の左の `≲` が等式である**算術の中身**を取った。
そこに残っていた「幾何の紐」は 3 本である:

| 本 | 内容 | 状態 |
|---|---|---|
| 1 | `log N(P) = f_P · log N(p)` | ★**本ファイル** |
| 2 | `∑_{P∣p} e_P f_P = [K:F]` | ★mathlib の `Ideal.sum_ramification_inertia`（実測） |
| 3 | 条件 (b) `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red` を導手の台の対応に翻訳 | ☆残っている |

★★★本ファイルは第 1 本を取る——`Found/Divisor/ArithFundamental.lean` の
`absNorm_eq_pow_inertiaDeg_rel`（`N(𝔓) = N(𝔭)^{f}`）に `log` を当てるだけである。

★★第 2 本は mathlib にそのままある（2026-08-29 実測）:

    `Ideal.sum_ramification_inertia S K L (hp : p ≠ ⊥) :
       ∑ P ∈ primesOverFinset p S, e(P∣p) * f(P∣p) = finrank K L`

★したがって `Proposition 1.7` に残るのは**第 3 本（幾何の条件 (b)）だけ**である。
-/

namespace ABC3.Found.GenEll

open NumberField Ideal

/-! ## ★★★★★★★★★★★★★★ノルムの対数は剰余次数倍 -/

/-- ★★★★★★★★★★★★★★**`log N(P) = f_P · log N(p)`**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★`Found/Divisor/ArithFundamental.lean` の `absNorm_eq_pow_inertiaDeg_rel`
（`N(𝔓) = N(𝔭)^{f(𝔓/𝔭)}`）に `log` を当てるだけである。
★★これが `§9-971` の帳簿に入れる第 1 本である。 -/
theorem log_absNorm_eq_inertiaDeg_mul (L M : Type) [Field L] [Field M]
    [NumberField L] [NumberField M] [Algebra L M]
    (P : Ideal (𝓞 M)) [P.IsPrime] (hPne : P ≠ ⊥) :
    Real.log (Ideal.absNorm P)
      = (((P.under (𝓞 L)).inertiaDeg P : ℕ) : ℝ)
        * Real.log (Ideal.absNorm (P.under (𝓞 L))) := by
  rw [ABC3.Found.Divisor.absNorm_eq_pow_inertiaDeg_rel L M P hPne]
  push_cast
  rw [Real.log_pow]

/-! ## ★出典の紐付け(`.src`) -/

def log_absNorm_eq_inertiaDeg_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log N(P) = f_P · log N(p))",
    sectionId := "genell-prop-1-7" }

def log_absNorm_eq_inertiaDeg_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "absNorm_eq_pow_inertiaDeg_rel(N(𝔓) = N(𝔭)^{f})"
      (.inProject "ABC3" "ABC3.Found.Divisor.absNorm_eq_pow_inertiaDeg_rel") 2,
    .citation "[mathlib]" "Ideal.sum_ramification_inertia(基本等式——幾何の紐の第 2 本)"
      (.inMathlib "Ideal.sum_ramification_inertia") 2,
    .implicitStep
      ("★★★★★測定(2026-08-29): §9-971 の帳簿に残っていた「幾何の紐」3 本のうち、" ++
       "第 1 本(log N(P) = f_P·log N(p))は本ファイルで、" ++
       "第 2 本(∑ e_P f_P = [K:F])は **mathlib にそのまま**ある" ++
       "(Ideal.sum_ramification_inertia、実測)。" ++
       "★したがって Proposition 1.7 に残るのは**第 3 本だけ**である" ++
       "——原文の条件 (b)(D_ℚ = φ_ℚ^{-1}(E_ℚ)_red)を導手の台の対応に翻訳する段") 6 ]

end ABC3.Found.GenEll
