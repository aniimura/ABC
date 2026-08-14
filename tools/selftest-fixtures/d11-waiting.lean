-- D11 fixture: まだ witness を作れないが、何を待っているか書いてある(通るべき)
namespace Fixture
structure Awaited where
  f : Nat → Nat
def Awaited.waiting : ABC3.Meta.WaitingFor :=
  { what := "局所類体論の相互法則", trackB := "Found/LocalCFT" }
end Fixture
