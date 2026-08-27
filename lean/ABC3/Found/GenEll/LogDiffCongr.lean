/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogDiffValue
import ABC3.Found.GenEll.MinFieldBaseChange

/-!
# [GenEll] Definition 1.5, (iii) —— **同型な数体は同じ `log-diff` をもつ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-defined log-different function log-diffX on X(Q).

## ★★これで well-defined 性の 2 つの側がそろう

`log-diff` を「点の**最小定義体**で測る」と定めたので、
`X(ℚ̄)` の関数として well-defined であるには 2 つが要る:

| 側 | 主張 | 取った場所 |
|---|---|---|
| **幾何** | 同じ点の最小定義体は互いに対応する | ★`minField_baseChange`（`MinFieldBaseChange.lean`） |
| **代数** | 対応する体は同じ `log-diff` を与える | ★★**本ファイル** |

## ★★★代数の側は 1 行で落ちた

`LogDiffValue.lean` が **`log-diff = log|disc F| / [F:ℚ]`** まで計算していたので、
残るのは

* `disc` が同型で不変 —— mathlib の **`NumberField.discr_eq_discr_of_algEquiv`**
* `[F:ℚ]` が同型で不変 —— `LinearEquiv.finrank_eq`

の 2 本だけである。★環同型が自動的に ℚ-代数同型になる（標数 0）のを挟むだけで繋がる。

★★★**`log-diff` を「型が付くだけの式」で止めず古典的な公式まで落としてあったことが、
ここで効いた。** 判別式まで降りていなければ、同型不変性は差積イデアルの
レベルで戦うことになっていた。
-/

namespace ABC3.Found.GenEll

open NumberField

/-- ★**環同型は自動的に ℚ-代数同型である**（標数 0 だから `ℚ` 上で可換）。 -/
noncomputable def ratAlgEquivOfRingEquiv {F K : Type} [Field F] [NumberField F]
    [Field K] [NumberField K] (e : F ≃+* K) : F ≃ₐ[ℚ] K :=
  { e with commutes' := fun q => by simp }

/-- ★★★★★**同型な数体は同じ `log-diff` をもつ**。

原文 (GenEll p.8):
> determines a well-defined log-different function log-diffX on X(Q).

★★これが `Definition 1.5, (iii)` の well-defined 性の**代数の側**である
（幾何の側は `minField_baseChange`）。

★中身は `log-diff = log|disc F| / [F:ℚ]`（`logDiffOfField_eq`）に
`discr_eq_discr_of_algEquiv` と `LinearEquiv.finrank_eq` を当てるだけである。 -/
theorem logDiffOfField_congr {F K : Type} [Field F] [NumberField F]
    [Field K] [NumberField K] (e : F ≃+* K) :
    logDiffOfField F = logDiffOfField K := by
  have halg : F ≃ₐ[ℚ] K := ratAlgEquivOfRingEquiv e
  have hdisc : NumberField.discr F = NumberField.discr K :=
    NumberField.discr_eq_discr_of_algEquiv F halg
  have hrank : Module.finrank ℚ F = Module.finrank ℚ K :=
    LinearEquiv.finrank_eq halg.toLinearEquiv
  rw [logDiffOfField_eq, logDiffOfField_eq, hdisc, hrank]

/-- ★★同型な数体は同じ次数をもつ（上の証明で使った部分を単独で）。 -/
theorem finrank_congr {F K : Type} [Field F] [NumberField F]
    [Field K] [NumberField K] (e : F ≃+* K) :
    Module.finrank ℚ F = Module.finrank ℚ K :=
  LinearEquiv.finrank_eq (ratAlgEquivOfRingEquiv e).toLinearEquiv

/-- ★★同型な数体は同じ判別式をもつ。 -/
theorem discr_congr {F K : Type} [Field F] [NumberField F]
    [Field K] [NumberField K] (e : F ≃+* K) :
    NumberField.discr F = NumberField.discr K :=
  NumberField.discr_eq_discr_of_algEquiv F (ratAlgEquivOfRingEquiv e)

/-! ### ★出典の紐付け(`.src`) -/

def ratAlgEquivOfRingEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(環同型は自動的に ℚ-代数同型)",
    sectionId := "genell-def-1-5" }

def logDiffOfField_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(同型な数体は同じ log-diff をもつ)",
    sectionId := "genell-def-1-5" }

def logDiffOfField_congr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiffOfField_eq(log-diff = log|disc F| / [F:ℚ])"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_eq") 8,
    .citation "[mathlib]" "NumberField.discr_eq_discr_of_algEquiv(判別式は ℚ-代数同型で不変)"
      (.inMathlib "NumberField.discr_eq_discr_of_algEquiv") 8,
    .implicitStep
      ("★★これで Definition 1.5, (iii) の well-defined 性は 2 つの側とも取れた —— " ++
       "幾何の側は minField_baseChange(最小定義体は底変換で対応する)、" ++
       "代数の側が本定理である") 8 ]

end ABC3.Found.GenEll
