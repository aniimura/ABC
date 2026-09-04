/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.HyperbolicCoreExample
import ABC3.Skeleton.CorrHyp.Section5

/-!
# [CorrHyp] `Definition 5.2`(`hyperbolicCoreGeneral`)の `corrHypInstance` における実現

`Definition 3.1`(`HyperbolicCoreExample.lean`)と同じ「命名段」パターン
——`Section5.lean` の `hyperbolicCoreGeneral X h := D.core X` も非
arithmetic性の仮定 `h`(ここでは `D.Ext X` について)を定義の中で使わない。
`corrHypInstance` は `core := id`・`Ext := id`(`Instance.lean` 冒頭で
placeholder と明記済み)なので、`X := FG_SL2Z` に対して
`hyperbolicCoreGeneral_witness = FG_SL2Z` そのものになる——
`HyperbolicCoreExample.lean` と同じ透明性のスタイルで正直に記録する。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Skeleton.CorrHyp

/-- `corrHypInstance.Ext FG_SL2Z = FG_SL2Z`(`Ext := id`)の下でも
`FG_SL2Z` は非arithmetic——`FG_SL2Z_not_arithmetic`(`HyperbolicCoreExample.lean`)
と同じ理由(`MargulisArithmetic`/`ShimuraArithmetic` が恒偽)。 -/
theorem FG_SL2Z_ext_not_arithmetic :
    ¬ Arithmetic corrHypInstance (corrHypInstance.Gamma (corrHypInstance.Ext FG_SL2Z)) := by
  simp [Arithmetic, corrHypInstance]

/-- **[CorrHyp] `Definition 5.2`(`hyperbolicCoreGeneral`)の `corrHypInstance`
における実現**——`X := FG_SL2Z`。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def hyperbolicCoreGeneral_witness : corrHypInstance.Space :=
  hyperbolicCoreGeneral corrHypInstance FG_SL2Z FG_SL2Z_ext_not_arithmetic

def hyperbolicCoreGeneral_witness.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 12, item := "Definition 5.2", sectionId := "corrhyp-def-5-2" }

end ABC3.Found.CorrHyp
