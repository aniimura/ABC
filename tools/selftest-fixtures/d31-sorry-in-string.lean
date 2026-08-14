-- D31 fixture: 文字列リテラル内の "sorry" を台帳に数えない(通るべき)
-- ★2026-08-15 実測の誤検出。`.needs` の `LeanStatus.inProject` は
--   「あちらの宣言は sorry 無しである」と**文字列で**書く場所なので、
--   潰さないと台帳が +1 される(Skeleton/AbsTopIII/LogShell.lean で実際に起きた)。
namespace Fixture
theorem mentions : True := trivial
def mentions.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def mentions.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "あちらの補題は sorry 無しで通っている" 3 ]
end Fixture
