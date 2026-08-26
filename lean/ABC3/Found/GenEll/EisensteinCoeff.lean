import ABC3.Found.GenEll.EisensteinMatch
import Mathlib.NumberTheory.ModularForms.LevelOne.DimensionFormula
import Mathlib.NumberTheory.ModularForms.QExpansion

/-!
# GenEll 第 345 ブロック —— **★★★★★`E₄`・`E₆` の q 展開係数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★(i) の最後——`E₄³ − E₆² = 1728Δ` に向けて

第 344 までで `g₂ = (4π⁴/3)E₄`・`g₃ = (8π⁶/27)E₆` と
`g₂³ − 27g₃² = (64π¹²/27)(E₄³ − E₆²)` が出た。★残るのは **`E₄³ − E₆² = 1728Δ`** である。

★★その筋は **Sturm 境界**:`f := E₄³ − E₆² − 1728Δ` は重さ 12 の modular form で、
q 展開の位数が `12/12 = 1` より大きければ `f = 0`(mathlib の `sturm_bound_levelOne`)。
★★★したがって **0 次と 1 次の係数**を突き合わせればよい。

## ★★★★★2026-08-26 の実測——mathlib の在庫

| 段 | mathlib |
|---|---|
| `E_k` の q 展開係数(`σ` と Bernoulli で) | ✅ `EisensteinSeries.E_qExpansion_coeff` |
| `Δ` の 1 次係数 = `1` | ✅ `discriminant_qExpansion_coeff_one` |
| `Δ ≠ 0` | ✅ `discriminant_ne_zero` |
| Sturm 境界(レベル 1) | ✅ `ModularForm.sturm_bound_levelOne` |
| q 展開の積・和・スカラー倍 | ✅ `qExpansion_mul`・`_add`・`_smul` |
| **`E₄³ − E₆² = 1728Δ`** | ★**無い** |

★★★★係数の公式は `-(2k/B_k)·σ_{k-1}(m)` で、`B₄ = −1/30`・`B₆ = 1/42` を入れると
`240`・`−504` になる。★**Bernoulli 数の値は第 342 で確認したとおり在庫にある**
(`Found/GaloisRep/TateEquationAnalytic.lean`、Tate 曲線の `q` 展開のために積んだもの)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `E4_coeff_zero`・`E4_coeff_one` | ★★★★**`E₄ = 1 + 240q + …`** |
| `E6_coeff_zero`・`E6_coeff_one` | ★★★★**`E₆ = 1 − 504q + …`** |
-/

namespace ABC3.Found.GenEll

open UpperHalfPlane EisensteinSeries ArithmeticFunction ModularForm

/-! ## ★★★★`E₄` の係数 -/

/-- ★★★`E₄` の定数項は `1`。 -/
theorem E4_coeff_zero :
    (qExpansion 1 (ModularForm.E (k := 4) (by norm_num))).coeff 0 = 1 :=
  EisensteinSeries.E_qExpansion_coeff_zero (by norm_num) (by decide)

/-- ★★★★**`E₄` の `q` の係数は `240`**——`B₄ = −1/30` から。 -/
theorem E4_coeff_one :
    (qExpansion 1 (ModularForm.E (k := 4) (by norm_num))).coeff 1 = 240 := by
  rw [EisensteinSeries.E_qExpansion_coeff (by norm_num) (by decide) 1]
  norm_num [ABC3.Found.GaloisRep.bernoulli_four_val]

/-! ## ★★★★`E₆` の係数 -/

/-- ★★★`E₆` の定数項は `1`。 -/
theorem E6_coeff_zero :
    (qExpansion 1 (ModularForm.E (k := 6) (by norm_num))).coeff 0 = 1 :=
  EisensteinSeries.E_qExpansion_coeff_zero (by norm_num) (by decide)

/-- ★★★★**`E₆` の `q` の係数は `−504`**——`B₆ = 1/42` から。 -/
theorem E6_coeff_one :
    (qExpansion 1 (ModularForm.E (k := 6) (by norm_num))).coeff 1 = -504 := by
  rw [EisensteinSeries.E_qExpansion_coeff (by norm_num) (by decide) 1]
  norm_num [ABC3.Found.GaloisRep.bernoulli_six_val]

/-! ## ★出典の紐付け(`.src`) -/

def E4_coeff_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def E6_coeff_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
