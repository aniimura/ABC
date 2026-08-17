-- D37 fixture: 引用行の末尾に docstring 終端 `-/` が乗っていても照合できる(通るべき)
-- ★2026-08-17 の回帰: `-/` が引用本文に混ざり、照合が末尾 2 文字で落ちた
--   (「61/63 文字まで一致 / 次に来るはず: "-/"」)。
namespace Fixture
/-- 原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K. -/
theorem quotedEndsDocstring : True := trivial
def quotedEndsDocstring.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def quotedEndsDocstring.needs : List ABC3.Meta.ProofObligation := []
end Fixture
