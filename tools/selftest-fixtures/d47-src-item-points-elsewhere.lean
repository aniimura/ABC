-- D47 fixture: `.src` の `item` と `sectionId` が**別の項目**を指している(落とすべき)。
-- ★これが素通りすると locator は静かにずれる(backlog M31。実木で 3 件見つかった)。
--   `prop-1-1` は実在し、pdfPage も範囲内なので、旧 G1 は 3 つとも通していた。
namespace Fixture
def implWithDriftedSrc : True := trivial
def implWithDriftedSrc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.2", sectionId := "prop-1-1" }
end Fixture
