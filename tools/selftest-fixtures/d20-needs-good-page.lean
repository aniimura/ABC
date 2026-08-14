-- D20 fixture: .needs のページが範囲内(通るべき)
namespace Fixture
theorem goodpage : True := trivial
def goodpage.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def goodpage.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "範囲内のページ" 3 ]
end Fixture
