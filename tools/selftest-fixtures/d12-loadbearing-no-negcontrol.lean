-- D12 fixture: load-bearing の印はあるが負の対照が無い(落ちるべき)
namespace Fixture
theorem central : True := trivial
def central.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def central.loadBearing : ABC3.Meta.LoadBearing := { consumer := "Theorem 4.2" }
def central.needs : List ABC3.Meta.ProofObligation := []
end Fixture
