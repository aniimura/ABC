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

/-! ## ★★★★★★2 つの側を繋ぐ -/

/-- ★**部分体の像は元の部分体と同型**（体からの射は単射だから）。 -/
noncomputable def subfieldMapEquiv {F K : Type} [Field F] [Field K]
    (E : Subfield F) (f : F →+* K) : E ≃+* (E.map f) := by
  refine RingEquiv.ofBijective
    ({ toFun := fun x => ⟨f x.1, ⟨x.1, x.2, rfl⟩⟩
       map_one' := by ext; simp
       map_mul' := fun a b => by ext; simp
       map_zero' := by ext; simp
       map_add' := fun a b => by ext; simp } : E →+* (E.map f)) ⟨?_, ?_⟩
  · intro a b hab
    exact Subtype.ext (f.injective (congrArg Subtype.val hab))
  · rintro ⟨_, x, hx, rfl⟩
    exact ⟨⟨x, hx⟩, rfl⟩

open AlgebraicGeometry CategoryTheory in
/-- ★★★★★★**同じ点の最小定義体は同じ `log-diff` を与える**。

原文 (GenEll p.8):
> determines a well-defined log-different function log-diffX on X(Q).

★★★**これが `Definition 1.5, (iii)` の well-defined 性そのもの**である ——
幾何の側（`minField_baseChange`：最小定義体は底変換で対応する）と
代数の側（`logDiffOfField_congr`：同型な数体は同じ `log-diff`）を繋いだ形。

★対比: `logDiffAt`（点の**定義体**で測る）は底変換で**増える**
（`logDiffAt_le_baseChange`）。★★**最小定義体で測ると不変になる** —— それが本定理である。 -/
theorem logDiff_minField_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] {X : Scheme}
    (xF : Spec (CommRingCat.of F) ⟶ X) :
    haveI := NumberField.of_subfield (minField F xF)
    haveI := NumberField.of_subfield
      (minField K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF))
    logDiffOfField (minField K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF))
      = logDiffOfField (minField F xF) := by
  haveI := NumberField.of_subfield (minField F xF)
  haveI := NumberField.of_subfield
    (minField K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF))
  haveI := NumberField.of_subfield ((minField F xF).map (algebraMap F K))
  -- ★`minField K (bc) ≃+* (minField F xF).map alg ≃+* minField F xF`
  have e : (minField K (Spec.map (CommRingCat.ofHom (algebraMap F K)) ≫ xF))
      ≃+* ((minField F xF).map (algebraMap F K)) :=
    RingEquiv.subfieldCongr (minField_baseChange F K xF)
  have e2 : (minField F xF) ≃+* ((minField F xF).map (algebraMap F K)) :=
    subfieldMapEquiv (minField F xF) (algebraMap F K)
  exact (logDiffOfField_congr e).trans (logDiffOfField_congr e2).symm

/-! ### ★出典の紐付け(`.src`) -/

def subfieldMapEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(部分体の像は元の部分体と同型)",
    sectionId := "genell-def-1-5" }

def logDiff_minField_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(同じ点の最小定義体は同じ log-diff を与える)",
    sectionId := "genell-def-1-5" }

def logDiff_minField_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "minField_baseChange(幾何の側——最小定義体は底変換で対応する)"
      (.inProject "ABC3" "ABC3.Found.GenEll.minField_baseChange") 8,
    .citation "[ABC3]" "logDiffOfField_congr(代数の側——同型な数体は同じ log-diff)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_congr") 8,
    .implicitStep
      ("★★★これが Definition 1.5, (iii) の well-defined 性そのものである。" ++
       "★対比: logDiffAt(点の定義体で測る)は底変換で増える(logDiffAt_le_baseChange)。" ++
       "最小定義体で測ると不変になる —— それが本定理") 8 ]

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
