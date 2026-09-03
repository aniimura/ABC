-- D40 fixture: sorry 0 の条なし `.src` が Skeleton にある(落ちるべき)
-- ★非空虚性の対照はあるので G9 は鳴らない —— G8 だけを較正する。
namespace Fixture
theorem bare_done : True := trivial
def bare_done.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def bare_done.needs : List ABC3.Meta.ProofObligation := []
theorem bare_done_nonvacuous : True := trivial
def bare_done_nonvacuous.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1(非空虚性の対照)", sectionId := "prop-1-1" }
def bare_done_nonvacuous.needs : List ABC3.Meta.ProofObligation := []
end Fixture
