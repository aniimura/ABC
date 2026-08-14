-- D17 fixture: コメント内の引用が該当ページに無い(落ちるべき)
namespace Fixture
/-- 原文 (pGC p.3):
> This sentence certainly does not appear on page three of that paper.
-/
theorem quoted : True := trivial
def quoted.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def quoted.needs : List ABC3.Meta.ProofObligation := []
end Fixture
