-- D30 fixture: **同一論文内**の辺も、辺として範囲検査される(落ちるべき)
-- ★2026-08-15 の測定: `otherPaper` という名前に反して、tag は同じ論文でよい。
--   実際、現在の 11 本の辺のうち 6 本が同一論文内である。
--   名前が誤解を招くので、ここで「同一論文内でも検査される」を固定する。
namespace Fixture
theorem sameedge : True := trivial
def sameedge.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def sameedge.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[pGC]" "Proposition 1.2" 999 ]
end Fixture
