import ABC3.Meta.Claim
import ABC3.Interface.LocProP.WeightZeroQuotient

/-!
# [LocProP] §3 —— 実装(`WeightZeroQuotientSetup` 上で)

`Skeleton/LocProP/Section3.lean` から委譲される実装本体。posit している
主張そのものを取り出すだけ——本体は posit 側にある。
-/

namespace ABC3.Found.LocProP

open ABC3.Interface.LocProP

def kerCX (E : WeightZeroQuotientSetup) := E.lemma_3_1
def zXParametrized (E : WeightZeroQuotientSetup) := E.ZX
def etaVanishes (E : WeightZeroQuotientSetup) := E.prop_3_3
def zXSingle (E : WeightZeroQuotientSetup) := E.ZX
def zSplits (E : WeightZeroQuotientSetup) := E.prop_3_5

def kerCX.nonvacuous := kerCX WeightZeroQuotientSetup.example
def etaVanishes.nonvacuous := etaVanishes WeightZeroQuotientSetup.example
def zSplits.nonvacuous := zSplits WeightZeroQuotientSetup.example

def kerCX.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 22, item := "Lemma 3.1", sectionId := "locprop-lemma-3-1" }
def zXParametrized.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 24, item := "Definition 3.2",
    sectionId := "locprop-def-3-2" }
def etaVanishes.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 25, item := "Proposition 3.3",
    sectionId := "locprop-prop-3-3" }
def zXSingle.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 27, item := "Definition 3.4",
    sectionId := "locprop-def-3-4" }
def zSplits.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 27, item := "Proposition 3.5",
    sectionId := "locprop-prop-3-5" }

end ABC3.Found.LocProP
