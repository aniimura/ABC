-- D33 fixture: Found/ の `.src` が 1_Structured に無い sectionId を指す(落とすべき)。
-- ★`Found/` の `.src` は**任意**だが、書いたら `Skeleton/` と同じ厳しさで検証する。
--   これが無いと「実装の件数」を嘘の出典で水増しできてしまう。
namespace Fixture
def implWithBadSrc : True := trivial
def implWithBadSrc.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1",
    sectionId := "this-section-does-not-exist" }
end Fixture
