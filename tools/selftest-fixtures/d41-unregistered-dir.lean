-- D41 fixture: 論文でも登記済みの理論でもないディレクトリに置かれている(落ちるべき)
-- ★sorry があるので G8 は鳴らず、nonvacuous があるので G9 も鳴らない —— G10 だけを較正する。
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
