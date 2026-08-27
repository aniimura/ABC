/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.FlatModel
import ABC3.Found.GenEll.ProjectiveSpace
import ABC3.Meta.Claim

/-!
# **`ℤ`-射影かつ `ℤ`-平坦なモデル**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

## ★★★★★★★★★台帳 `ample-and-projective-embedding` の段 F2（の 2/3）

段 F は原文第 2 文の括弧「`ℤ`-flat, `ℤ`-projective model」である。

| 段 | 内容 | 状態 |
|---|---|---|
| F1 | 閉包（スキーム論的像）の切断が **`ℤ`-平坦** | ★閉じた（`FlatModel.lean`） |
| F2a | 閉包が **`ℤ`-射影**（`ℙᴺ_ℤ` への閉埋め込み） | ★★**本ファイル**（`subschemeι` がそのまま閉埋め込み） |
| F2b | 閉包が **`ℤ`-固有** | ★★**本ファイル**（`ℙᴺ_ℤ` の固有性と合成） |
| F2c | 生成ファイバーが元の `Y` に戻ること | ★開 |

## ★★★★★★機構は 2 行である

`Found/GenEll/ProjectiveSpace.lean` で `ℙⁿ_R` が `Spec R` 上固有だと示した。
★閉埋め込みは固有であり、固有射の合成は固有だから、

    `(ker f).subschemeι ≫ projSpaceOverSpec N ℤ`

は固有である。★★`IsProper` は合成の instance を持たないので
`exact IsProper.mk` と書く（`tools\lean-idioms.md`）。

★★★**射影性は別に示す必要がない**——`subschemeι` は mathlib の instance で
すでに閉埋め込みである。「`ℤ`-射影」とは `ℙᴺ_ℤ` の閉部分スキームであることだからである。

## ★残っている段（明示）

★★**生成ファイバーが元の `Y` に戻ること**（F2c）だけである。
それは「`Y` が `ℚ`-スキームなら、スキーム論的像の `ℚ` への底変換が `Y` に等しい」
という主張で、`Y ⟶ image` が開埋め込み（かつ全射）であることに帰着する。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-- ★★★★★★**スキーム論的像は `ℤ` 上固有である**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

★機構は「閉埋め込みは固有」＋「`ℙᴺ_ℤ` は `Spec ℤ` 上固有」＋「固有射の合成は固有」。 -/
theorem isProper_subscheme_projSpace (N : ℕ) {Y : Scheme.{0}}
    (f : Y ⟶ projSpace N ℤ) [QuasiCompact f] :
    IsProper ((Scheme.Hom.ker f).subschemeι ≫ projSpaceOverSpec N ℤ) :=
  IsProper.mk

/-- ★★★★★★★★★**`ℤ`-射影かつ `ℤ`-平坦なモデル**。

原文 (GenEll p.8):
> Remark 1.4.1. Observe that it follows immediately from the definitions, together with Proposition 1.4, (iii), that the theory of

`ℚ`-スキーム `Y` を `ℙᴺ_ℤ` の中で閉包すると、

1. ★`ℙᴺ_ℤ` への**閉埋め込み**である（＝`ℤ`-射影）
2. ★★`Spec ℤ` 上**固有**である
3. ★★★切断が **`ℤ`-平坦**である（`FlatModel.lean`）

★★★★これが原文第 2 文の括弧の中身である。残るのは
「生成ファイバーが元の `Y` に戻ること」だけである。 -/
theorem projective_flat_model (N : ℕ) {Y : Scheme.{0}}
    (f : Y ⟶ projSpace N ℤ) [QuasiCompact f]
    (U : (projSpace N ℤ).affineOpens)
    [Algebra ℚ Γ(Y, f ⁻¹ᵁ (U : (projSpace N ℤ).Opens))] :
    IsClosedImmersion (Scheme.Hom.ker f).subschemeι
    ∧ IsProper ((Scheme.Hom.ker f).subschemeι ≫ projSpaceOverSpec N ℤ)
    ∧ Module.Flat ℤ Γ((Scheme.Hom.ker f).subscheme,
        (Scheme.Hom.ker f).subschemeι ⁻¹ᵁ (U : (projSpace N ℤ).Opens)) :=
  ⟨inferInstance, IsProper.mk, flat_int_subschemeObj f U⟩

/-! ### ★出典の紐付け(`.src`) -/

def isProper_subscheme_projSpace.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℚ-スキームの閉包が ℤ-固有であること)",
    sectionId := "genell-rem-1-4-1" }

def projective_flat_model.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Remark 1.4.1(ℤ-射影かつ ℤ-平坦なモデル。生成ファイバーが Y に戻ることは含まない)",
    sectionId := "genell-rem-1-4-1" }

def projective_flat_model.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "projSpaceOverSpec_isProper(ℙⁿ_R は Spec R 上固有)"
      (.inProject "ABC3" "ABC3.Found.GenEll.projSpaceOverSpec_isProper") 8,
    .citation "[ABC3]" "flat_int_subschemeObj(閉包の切断が ℤ-平坦)"
      (.inProject "ABC3" "ABC3.Found.GenEll.flat_int_subschemeObj") 8,
    .citation "[mathlib]" "Scheme.IdealSheafData.subschemeι(閉埋め込み)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.subschemeι") 8,
    .implicitStep
      ("★射影性は別に示す必要がない——subschemeι は mathlib の instance で" ++
       "すでに閉埋め込みであり、「ℤ-射影」とは ℙᴺ_ℤ の閉部分スキームであることだからである") 8,
    .implicitStep
      ("★★残っているのは**生成ファイバーが元の Y に戻ること**(F2c)だけである。" ++
       "それは Y ⟶ image の ℚ への底変換が同型であることに帰着する") 8 ]

end ABC3.Found.GenEll
