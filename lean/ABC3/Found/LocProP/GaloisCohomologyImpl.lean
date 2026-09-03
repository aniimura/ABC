import ABC3.Meta.Claim
import ABC3.Interface.LocProP.GaloisCohomologyReview

/-!
# [LocProP] §2 —— 実装(`GaloisCohomologyReviewSetup` 上で)

`Skeleton/LocProP/Section2.lean` から委譲される実装本体。posit している
主張(同型の**存在**として `Nonempty (_ ≃+ _)`)を束ねるだけ——本体は posit 側にある。
★型注釈を書かない(`E.H0DeltaRprf` 等の instance が外部コンテキストで見えない罠、
`tools/lean-idioms.md` 参照)。
-/

namespace ABC3.Found.LocProP

open ABC3.Interface.LocProP

def lemma21 (E : GaloisCohomologyReviewSetup) :=
  And.intro (fun j => Nonempty.intro (E.lemma_2_1_i j)) (fun j => Nonempty.intro (E.lemma_2_1_ii j))

def lemma22 (E : GaloisCohomologyReviewSetup) :=
  And.intro (Nonempty.intro E.lemma_2_2_i) (And.intro E.lemma_2_2_ii (Nonempty.intro E.lemma_2_2_iii))

def lemma23 (E : GaloisCohomologyReviewSetup) :=
  And.intro (fun n => Nonempty.intro (E.lemma_2_3_i n)) E.lemma_2_3_ii

def prop24 (E : GaloisCohomologyReviewSetup) :=
  And.intro E.prop24_left_injective (And.intro E.prop24_right_surjective E.prop24_exact)

def lemma21.nonvacuous := lemma21 GaloisCohomologyReviewSetup.example
def lemma22.nonvacuous := lemma22 GaloisCohomologyReviewSetup.example
def lemma23.nonvacuous := lemma23 GaloisCohomologyReviewSetup.example
def prop24.nonvacuous := prop24 GaloisCohomologyReviewSetup.example

def lemma21.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.1", sectionId := "locprop-lemma-2-1" }
def lemma22.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.2", sectionId := "locprop-lemma-2-2" }
def lemma23.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.3", sectionId := "locprop-lemma-2-3" }
def prop24.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 21, item := "Proposition 2.4",
    sectionId := "locprop-prop-2-4" }

end ABC3.Found.LocProP
