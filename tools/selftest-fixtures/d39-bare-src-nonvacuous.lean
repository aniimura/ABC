-- D39 fixture: 条なし `.src` + 非空虚性の対照あり(通るべき)
-- ★本体に sorry があるので G8 も鳴らない。
namespace Fixture
theorem bare_ok : True := sorry
def bare_ok.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def bare_ok.needs : List ABC3.Meta.ProofObligation := []
theorem bare_ok_nonvacuous : True := trivial
def bare_ok_nonvacuous.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1(非空虚性の対照)", sectionId := "prop-1-1" }
def bare_ok_nonvacuous.needs : List ABC3.Meta.ProofObligation := []
end Fixture
