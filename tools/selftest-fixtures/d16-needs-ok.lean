-- D16 fixture: .needs があり、空リストも明示的な主張として認める(通るべき)
namespace Fixture
theorem scoped : True := trivial
def scoped.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def scoped.needs : List ABC3.Meta.ProofObligation := []
end Fixture
