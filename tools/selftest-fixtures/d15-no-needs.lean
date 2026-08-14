-- D15 fixture: Skeleton の theorem に .needs が無い(落ちるべき)
namespace Fixture
theorem unscoped : True := trivial
def unscoped.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
end Fixture
