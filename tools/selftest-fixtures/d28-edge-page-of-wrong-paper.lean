-- D28 fixture: otherPaper のページが「引用元の論文」では範囲内だが
--              「辺の先の論文」では範囲外(落ちるべき)。
-- ★これが 2026-08-15 に見つかった実際の穴の形である。
--   以前の器具は `.needs` のページを **所有論文** の pdfPages と比べていたので、
--   別論文を指す辺は事実上検査されていなかった。
--   ここでは所有論文 IUTchIII(199 ページ)では 150 は範囲内、
--   辺の先 pGC(9 ページ)では範囲外。
namespace Fixture
theorem wrongpaperpage : True := trivial
def wrongpaperpage.src : ABC3.Meta.Source :=
  { paper := "IUTchIII", pdfPage := 174, item := "Corollary 3.12", sectionId := "cor-3-12-theta" }
def wrongpaperpage.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[pGC]" "Proposition 1.1" 150 ]
end Fixture
