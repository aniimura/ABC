-- D19 fixture: .needs が範囲外の物理ページを指す(落ちるべき)
namespace Fixture
theorem badpage : True := trivial
def badpage.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def badpage.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "存在しないページを指す" 999 ]
end Fixture
