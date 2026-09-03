import ABC3.Meta.Claim
import ABC3.Interface.LocProP.GaloisCohomologyReview
import ABC3.Found.LocProP.GaloisCohomologyImpl

/-!
# [LocProP] §2 —— Review of Galois Cohomology(4 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.20-21。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

Lemma 2.1-2.3・Proposition 2.4 はすべて `Interface/LocProP/GaloisCohomologyReview.lean`
の `GaloisCohomologyReviewSetup`(Faltings の almost étale extension 理論・Hodge-Tate
分解を posit したもの)の上で述べる。★§2 は 4/4(2026-09-04 時点)。
-/

namespace ABC3.Skeleton.LocProP

open ABC3.Interface.LocProP

def galoisCohomology21 (E : GaloisCohomologyReviewSetup) := ABC3.Found.LocProP.lemma21 E
def galoisCohomology21.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.1", sectionId := "locprop-lemma-2-1" }
def galoisCohomology21.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[Falt1]" "Theorem 4.4, (i)(iv)" 17 ]

def galoisCohomology22 (E : GaloisCohomologyReviewSetup) := ABC3.Found.LocProP.lemma22 E
def galoisCohomology22.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.2", sectionId := "locprop-lemma-2-2" }
def galoisCohomology22.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[Tate]" "Theorem 1(i)(iii)" 10, .otherPaper "[Tate]" "Theorem 2(ii)" 10 ]

def galoisCohomology23 (E : GaloisCohomologyReviewSetup) := ABC3.Found.LocProP.lemma23 E
def galoisCohomology23.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 20, item := "Lemma 2.3", sectionId := "locprop-lemma-2-3" }

def hodgeTateExactSeq (E : GaloisCohomologyReviewSetup) := ABC3.Found.LocProP.prop24 E
def hodgeTateExactSeq.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 21, item := "Proposition 2.4",
    sectionId := "locprop-prop-2-4" }
def hodgeTateExactSeq.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "原文は[Falt2](Crystalline Cohomology and p-adic Galois Representations、
       JHU Press 1989)のTheorem 6.2(Comparison Theorem)を引くが、[Falt2]は
       paywallで未取得。代わりにBrinon-Conrad(papers.json短縮タグBC)の
       Theorem 2.2.3(同じ内容の comparison isomorphism)を等価な代替典拠
       として登記した(2026-09-04)。" 21,
    .otherPaper "[BC]" "Theorem 2.2.3" 13 ]

end ABC3.Skeleton.LocProP
