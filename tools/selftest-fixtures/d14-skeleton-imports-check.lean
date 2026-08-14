-- D14 fixture: Skeleton から Check を import している(落ちるべき)
import ABC3.Check.PGC.Section1Discriminating
namespace Fixture
theorem laundered : True := trivial
def laundered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def laundered.needs : List ABC3.Meta.ProofObligation := []
end Fixture
