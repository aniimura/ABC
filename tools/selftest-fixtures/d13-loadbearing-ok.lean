-- D13 fixture: load-bearing + 負の対照あり(通るべき)
namespace Fixture
theorem central : True := trivial
def central.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def central.loadBearing : ABC3.Meta.LoadBearing := { consumer := "Theorem 4.2" }
def central.negControl : ABC3.Meta.NegControl :=
  { dropped := "フィルトレーションの保存", witness := "Fixture.central_fails_without_filt" }
def central.needs : List ABC3.Meta.ProofObligation := []
end Fixture
