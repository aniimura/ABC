-- D38 fixture: 条なし `.src` なのに非空虚性の対照が無い(落ちるべき)
-- ★本体に sorry があるので G8 は鳴らない —— G9 だけを較正する。
namespace Fixture
theorem bare_no_witness : True := sorry
def bare_no_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def bare_no_witness.needs : List ABC3.Meta.ProofObligation := []
end Fixture
