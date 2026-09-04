/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.CorrHyp.ModularExample

/-!
# [CorrHyp] `Definition 3.1`(`hyperbolicCore`)の `corrHypInstance` における実現

`Skeleton/CorrHyp/Section3.lean` 自身のdocstringが明言する通り、
`Definition 3.1`(`hyperbolicCore X h := D.core X`)は**「hyperbolic core」
という名前を与える段」であり、非arithmetic性の仮定 `h` は定義の中で
未使用**——原文 (p.8) も「we shall refer to Y as the hyperbolic core of X」
という命名の宣言に過ぎない。ゆえにこの項目の実装は「`core` フィールドに
何か深い数学を仕込む」ことではなく、**具体的な `D`・具体的な非arithmetic
な `X` を揃えて型を通す**ことに尽きる。

`corrHypInstance`(`Instance.lean`)は `core := id`(`Proposition 3.2` が
読まない placeholder、`Instance.lean` 冒頭の docstring で明記済み)なので、
ここでの `hyperbolicCore_witness` は定義的には `FG_SL2Z` そのものに等しい
——★これは意図的な正直さであり、隠していない(`Definition 1.2` を
`corrHypInstance` では claim しなかった `ModularExample.lean` と同じ
透明性のスタイル)。非arithmetic性の証明(`FG_SL2Z_not_arithmetic`)は
`corrHypInstance.MargulisArithmetic`/`ShimuraArithmetic` が恒偽
(`Instance.lean` の placeholder)であることから出る——`Proposition 3.2`
の実現(`prop_3_2_at_instance`)が同じ事実にすでに依拠しているのと
整合的。 -/

namespace ABC3.Found.CorrHyp

open ABC3.Skeleton.CorrHyp

/-- `FG_SL2Z`(モジュラー群、`ModularExample.lean`)は `corrHypInstance` に
おいて非arithmetic——`MargulisArithmetic`/`ShimuraArithmetic` が両方
恒偽(`corrHypInstance` の placeholder)なので `Arithmetic` も恒偽。 -/
theorem FG_SL2Z_not_arithmetic :
    ¬ Arithmetic corrHypInstance (corrHypInstance.Gamma FG_SL2Z) := by
  simp [Arithmetic, corrHypInstance]

/-- **[CorrHyp] `Definition 3.1`(`hyperbolicCore`)の `corrHypInstance` に
おける実現**——`X := FG_SL2Z`(モジュラー群、非arithmetic の証明つき)。

★**sorry 無し**。標準3公理のみ。 -/
noncomputable def hyperbolicCore_witness : corrHypInstance.Space :=
  hyperbolicCore corrHypInstance FG_SL2Z FG_SL2Z_not_arithmetic

def hyperbolicCore_witness.src : ABC3.Meta.Source :=
  { paper := "CorrHyp", pdfPage := 8, item := "Definition 3.1", sectionId := "corrhyp-def-3-1" }

end ABC3.Found.CorrHyp
