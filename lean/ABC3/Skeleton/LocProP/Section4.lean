import ABC3.Meta.Claim
import ABC3.Interface.LocProP.JGeometric
import ABC3.Found.LocProP.JGeometricImpl

/-!
# [LocProP] §4 —— J-Geometric Sections(4 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.30-31。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

Lemma 4.1・4.2・Definition 4.3・Proposition 4.4 はすべて
`Interface/LocProP/JGeometric.lean` の `JGeometricSetup` の上で述べる。
★§4 は 4/4(2026-09-04 時点)。これで LocProP §0-§4 の 21 項目すべてが揃った。
-/

namespace ABC3.Skeleton.LocProP

open ABC3.Interface.LocProP

def geometricPointCriterion (E : JGeometricSetup) := ABC3.Found.LocProP.lPointExists E
def geometricPointCriterion.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 30, item := "Lemma 4.1", sectionId := "locprop-lemma-4-1" }
def geometricPointCriterion.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[BK]" "Lemma 3.8.1" 23, .otherPaper "[BK]" "Example 3.11" 29 ]

def actionConjugationFormula (E : JGeometricSetup) := ABC3.Found.LocProP.conjugationFormula E
def actionConjugationFormula.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Lemma 4.2", sectionId := "locprop-lemma-4-2" }

def isJGeometric (E : JGeometricSetup) := ABC3.Found.LocProP.isJGeometricDef E
def isJGeometric.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Definition 4.3",
    sectionId := "locprop-def-4-3" }

def splCriterion (E : JGeometricSetup) := ABC3.Found.LocProP.splIffJGeometric E
def splCriterion.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 31, item := "Proposition 4.4",
    sectionId := "locprop-prop-4-4" }

end ABC3.Skeleton.LocProP
