import ABC3.Meta.Claim
import ABC3.Interface.LocProP.JGeometric

/-!
# [LocProP] §4 —— 実装(`JGeometricSetup` 上で)

`Skeleton/LocProP/Section4.lean` から委譲される実装本体。posit している
主張そのものを取り出すだけ——本体は posit 側にある。
-/

namespace ABC3.Found.LocProP

open ABC3.Interface.LocProP

def lPointExists (E : JGeometricSetup) := E.lemma_4_1
def conjugationFormula (E : JGeometricSetup) := E.lemma_4_2
def isJGeometricDef (E : JGeometricSetup) := E.isJGeometric
def splIffJGeometric (E : JGeometricSetup) := E.prop_4_4

def lPointExists.nonvacuous := lPointExists JGeometricSetup.example
def conjugationFormula.nonvacuous := conjugationFormula JGeometricSetup.example
def splIffJGeometric.nonvacuous := splIffJGeometric JGeometricSetup.example

def lPointExists.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 30, item := "Lemma 4.1", sectionId := "locprop-lemma-4-1" }
def conjugationFormula.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Lemma 4.2", sectionId := "locprop-lemma-4-2" }
def isJGeometricDef.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Definition 4.3",
    sectionId := "locprop-def-4-3" }
def splIffJGeometric.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Proposition 4.4",
    sectionId := "locprop-prop-4-4" }

end ABC3.Found.LocProP
