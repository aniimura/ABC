-- D29 fixture: docstring が `foo.needs` に言及していても、本体を取り違えない(落ちるべき)
-- ★2026-08-15 実測の失敗形。以前の器具は `src.indexOf("foo.needs")` で
--   **docstring 中の散文**を先に拾い、そこから最初の `[` までを本体と誤認していた。
--   その結果 `.needs` は「件数」には数えられるのに、中身は 1 件も集計されず、
--   ページ範囲も辺も検査されないまま素通りした。
namespace Fixture
/-- 探索範囲は `mention.needs` に書いた [これは散文であって本体ではない] -/
theorem mention : True := trivial
def mention.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }
def mention.needs : List ABC3.Meta.ProofObligation :=
  [ .folklore "本体はこちら。範囲外のページ" 999 ]
end Fixture
