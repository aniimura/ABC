-- D18 fixture: 装飾記法つきの正しい引用(通るべき)
namespace Fixture
/-- 原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.
-/
theorem quotedOk : True := trivial
def quotedOk.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def quotedOk.needs : List ABC3.Meta.ProofObligation := []
end Fixture
