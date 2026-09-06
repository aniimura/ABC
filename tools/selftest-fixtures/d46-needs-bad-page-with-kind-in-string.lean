-- D46 fixture: D45 の対。**文字列を潰しても本物の範囲外は落とせる**ことの確認。
-- 文字列の中に `.otherPaper` と綴ってあるが、obligation は 1 件であり、
-- その物理ページ 999 は pGC(9 頁)の範囲外なので **落ちる** べき。
-- (D45 だけだと「マスクを入れて全部素通りにした」でも PASS してしまう。)
namespace Fixture
theorem badpagewithspelling : True := trivial
def badpagewithspelling.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def badpagewithspelling.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "本文に .otherPaper と綴ってあるが、これは記録であって区切りではない" 999 ]
end Fixture
