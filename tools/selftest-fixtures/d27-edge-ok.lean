-- D27 fixture: otherPaper が登記済みの論文の範囲内のページを指す(通るべき)
-- 「全部落とす器具は全部通す器具と同じくらい無情報」なので、正しい辺の対照を置く。
namespace Fixture
theorem goodedge : True := trivial
def goodedge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def goodedge.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[pGC]" "Proposition 1.2" 5 ]
end Fixture
