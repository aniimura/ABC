import ABC3.Meta.Claim
import ABC3.Interface.LocProP.WeightZeroQuotient
import ABC3.Found.LocProP.WeightZeroQuotientImpl

/-!
# [LocProP] §3 —— The Weight Zero Quotient(5 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.22-27。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

Lemma 3.1・Definition 3.2/3.4・Proposition 3.3/3.5 はすべて
`Interface/LocProP/WeightZeroQuotient.lean` の `WeightZeroQuotientSetup`
(Malčev completion・Tannakian category を posit したもの)の上で述べる。
★§3 は 5/5(2026-09-04 時点)。
-/

namespace ABC3.Skeleton.LocProP

open ABC3.Interface.LocProP

def kernelOfCX (E : WeightZeroQuotientSetup) := ABC3.Found.LocProP.kerCX E
def kernelOfCX.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 22, item := "Lemma 3.1", sectionId := "locprop-lemma-3-1" }

def weightZeroQuotientParam (E : WeightZeroQuotientSetup) := ABC3.Found.LocProP.zXParametrized E
def weightZeroQuotientParam.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 24, item := "Definition 3.2",
    sectionId := "locprop-def-3-2" }

def etaAlphaVanishes (E : WeightZeroQuotientSetup) := ABC3.Found.LocProP.etaVanishes E
def etaAlphaVanishes.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 25, item := "Proposition 3.3",
    sectionId := "locprop-prop-3-3" }

def weightZeroQuotientSingle (E : WeightZeroQuotientSetup) := ABC3.Found.LocProP.zXSingle E
def weightZeroQuotientSingle.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 27, item := "Definition 3.4",
    sectionId := "locprop-def-3-4" }

def weightZeroSequenceSplits (E : WeightZeroQuotientSetup) := ABC3.Found.LocProP.zSplits E
def weightZeroSequenceSplits.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 27, item := "Proposition 3.5",
    sectionId := "locprop-prop-3-5" }

end ABC3.Skeleton.LocProP
