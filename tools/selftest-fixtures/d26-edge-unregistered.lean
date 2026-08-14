-- D26 fixture: otherPaper が未登記の論文を指す(落ちるべき)
namespace Fixture
theorem badedge : True := trivial
def badedge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def badedge.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[NoSuchPaper]" "Theorem 1.1" 3 ]
end Fixture
