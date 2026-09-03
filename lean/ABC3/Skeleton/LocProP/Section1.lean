import ABC3.Meta.Claim
import ABC3.Interface.LocProP.OrdinaryCase
import ABC3.Found.LocProP.OrdinaryCaseImpl

/-!
# [LocProP] §1 —— The Ordinary Case(4 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.17-19。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

Definition 1.1・Lemma 1.2・1.3・1.4 はすべて `Interface/LocProP/OrdinaryCase.lean` の
`OrdinaryCaseSetup`(Jacobian・ordinary abelian variety・étale cohomology 経由の構成を
posit したもの)の上で述べる。★§1 は 4/4(2026-09-04 時点)。
-/

namespace ABC3.Skeleton.LocProP

open ABC3.Interface.LocProP

/-- **[LocProP] Definition 1.1** —— `E.isOrdinary` の posit そのもの。 -/
def isOrdinary (E : OrdinaryCaseSetup) : Prop := E.isOrdinary

def isOrdinary.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 17, item := "Definition 1.1",
    sectionId := "locprop-def-1-1" }

def isOrdinary.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "Jacobian(曲線の)とordinary abelian variety(p-rank = 次元)はmathlibに0件。
       Interface/LocProP/OrdinaryCase.lean のisOrdinary : Propとしてposit した
       (証明していない)。" 17 ]

/-- **[LocProP] Lemma 1.2**。実装は `Found/LocProP/OrdinaryCaseImpl.lean` へ委譲。 -/
def thetaXUnique (E : OrdinaryCaseSetup) := ABC3.Found.LocProP.thetaUnique E

def thetaXUnique.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 18, item := "Lemma 1.2", sectionId := "locprop-lemma-1-2" }

/-- **[LocProP] Lemma 1.3**。実装は `Found/LocProP/OrdinaryCaseImpl.lean` へ委譲。 -/
def deltaXEtKernel (E : OrdinaryCaseSetup) := ABC3.Found.LocProP.kerEqIInter E

def deltaXEtKernel.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 18, item := "Lemma 1.3", sectionId := "locprop-lemma-1-3" }

/-- **[LocProP] Lemma 1.4**。実装は `Found/LocProP/OrdinaryCaseImpl.lean` へ委譲。 -/
def alphaXEtEqTheta (E : OrdinaryCaseSetup) := ABC3.Found.LocProP.alphaEtEqTheta E

def alphaXEtEqTheta.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 19, item := "Lemma 1.4", sectionId := "locprop-lemma-1-4" }

end ABC3.Skeleton.LocProP
