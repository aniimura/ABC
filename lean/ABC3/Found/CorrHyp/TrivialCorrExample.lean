/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.Instance4

/-!
# [CorrHyp] `Definition 1.2`(`Corr.IsTrivial`)の非空虚な実現

`ModularExample.lean` の `corr_witness_isTrivial_degenerate` は
`corrHypInstance`(`Space := FuchsianGroup`、`FEt := PLift ∘ IsFiniteIndexIn`)
において `Corr.IsTrivial` が**任意の correspondence について自動的に真**に
なってしまう(`FEt` が `Prop` を包んだだけの型なので証明無関係性により
同じ `X Y` の間の `FEt` が常に唯一)ことを指摘し、`.src` を付けずに
「正直な限界」として記録していた。

`corrHypInstance4`(`Space := QcqsSpace`、`FEt := QcqsFEt`)ではこの問題が
**構造的に起こらない**——`QcqsFEt A B := {f : A.1 ⟶ B.1 // IsFinite f.left ∧
Etale f.left}` は実際のスキーム射の部分型であり `Prop` の詰め替えではない
ので、一般の `A ≠ B` に対しては `FEt A B` が空(有限エタール射が存在しない)
こともあれば複数の相異なる元を持つこともある——`Corr.IsTrivial` は
**c ごとに真偽が変わる非自明な述語**である。ここでは最も単純な(しかし
正真正銘の)例——`X = Y = C := basePt4`、`α = β = γ := QcqsFEt.idFEt`——
で `Definition 1.2` を実際に witness する。 -/

namespace ABC3.Found.CorrHyp

open CategoryTheory

/-- `FEtK` の恒等射どうしの合成は恒等射(圏の`id_comp`則、`Subtype.ext`で
証明成分を消し込む)。 -/
theorem FEtK_comp_idFEt_idFEt (X : Over BaseK) :
    FEtK.comp (FEtK.idFEt X) (FEtK.idFEt X) = FEtK.idFEt X := by
  apply Subtype.ext
  show (𝟙 X ≫ 𝟙 X : X ⟶ X) = 𝟙 X
  simp

open ABC3.Skeleton.CorrHyp in
/-- `corrHypInstance4` における最も単純な correspondence(`C = X = Y = basePt4`、
両脚とも恒等射)。 -/
noncomputable def corr_witness4 : Corr corrHypInstance4 basePt4 basePt4 :=
  ⟨basePt4, QcqsFEt.idFEt basePt4, QcqsFEt.idFEt basePt4⟩

open ABC3.Skeleton.CorrHyp in
/-- **[CorrHyp] `Definition 1.2`(`Corr.IsTrivial`)の `corrHypInstance4` に
おける実現**——`γ := idFEt basePt4` で `α = comp β γ` が(`id ≫ id = id`
から)成り立つ。`FEt` が本物の射の型であるこの instance では、この述語は
`c` ごとに真偽が変わる非自明な内容を持つ(`corrHypInstance` の場合との
違いは本ファイル冒頭のdocstringを参照)。

★**sorry 無し**。標準3公理のみ。 -/
theorem corr_witness4_isTrivial : Corr.IsTrivial corrHypInstance4 corr_witness4 := by
  refine ⟨QcqsFEt.idFEt basePt4, ?_⟩
  show QcqsFEt.idFEt basePt4 =
    corrHypInstance4.comp (QcqsFEt.idFEt basePt4) (QcqsFEt.idFEt basePt4)
  show FEtK.idFEt basePt4.1 = FEtK.comp (FEtK.idFEt basePt4.1) (FEtK.idFEt basePt4.1)
  rw [FEtK_comp_idFEt_idFEt]

def corr_witness4_isTrivial.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 4, item := "Definition 1.2", sectionId := "corrhyp-def-1-2" }

end ABC3.Found.CorrHyp
