-- D35 fixture: Found/ の `.src` が正しい locator を指す(通るべき)。
-- 「全部落とす器具は全部通す器具と同じくらい無情報」なので、正しい対照を置く。
namespace Fixture
def implWithGoodSrc : True := trivial
def implWithGoodSrc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
end Fixture
