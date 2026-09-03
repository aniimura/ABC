import ABC3.Meta.Claim
import ABC3.Interface.LocProP.EtaleSetup

/-!
# [LocProP] Lemma 0.4 —— 実装(`CohomologyComparisonSetup` 上で)

`Skeleton/LocProP/Section0.lean` から委譲される実装本体。
`Interface/LocProP/EtaleSetup.lean` の `compare_bijective` が posit している
主張そのものを取り出すだけなので証明は空(1 行)——本体は posit 側にある。
-/

namespace ABC3.Found.LocProP

open ABC3.Interface.LocProP

theorem cohomologyComparison_bijective (E : CohomologyComparisonSetup)
    (i r : ℤ) : Function.Bijective (E.compare i r) := E.compare_bijective i r

/-- ★★**非退化性の対照** —— `CohomologyComparisonSetup.nonvacuous` の具体例
(`Hgrp = Het = ℤ`、`compare = id`)で実際に成り立つことを確認する。 -/
theorem cohomologyComparison_bijective.nonvacuous :
    Function.Bijective (CohomologyComparisonSetup.nonvacuous.some.compare 0 0) :=
  cohomologyComparison_bijective CohomologyComparisonSetup.nonvacuous.some 0 0

/-- `Found/` 側でも本項目([LocProP] Lemma 0.4)が受け持たれていることを示す出典
(G8: `Skeleton/LocProP/Section0.lean` の `etaleGroupCohomologyComparison_bijective` から
本ファイルへ委譲していることの対照)。 -/
def cohomologyComparison_bijective.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 15, item := "Lemma 0.4",
    sectionId := "locprop-lemma-0-4" }

end ABC3.Found.LocProP
