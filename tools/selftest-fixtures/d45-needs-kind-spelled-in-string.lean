-- D45 fixture: `.needs` の**文字列の中**に obligation の綴り(先頭ドットつき)がある。
-- 区切り検出が文字列を潰さずに走ると、ここで obligation が 1 件 → 2 件に割れ、
-- 前半の断片が引用符の対応を失って **本文中の数値 1036 を物理ページとして拾う**。
-- pGC は 9 頁なので「p.1036 は範囲外」という**偽の NG** が出る(第 1036 で実際に出た)。
-- 正しくは 1 件・p.3 と読めて **通る** べき。
namespace Fixture
theorem spelledinstring : True := trivial
def spelledinstring.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def spelledinstring.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "第 1036 で .implicitStep という綴りを本文から外して回避した" 3 ]
end Fixture
