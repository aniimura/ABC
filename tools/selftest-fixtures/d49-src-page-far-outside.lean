-- D49 fixture: `.src` の `pdfPage` が項目の頁の範囲から**大きく外れ**、
-- しかもそのファイルは**その頁の逐語を 1 行も引いていない**(落とすべき)。
-- ★規則 R′(backlog M43、メタ第 13〜14 回)の退行の見張り。
--   `prop-1-1` は p.3 の項目で、範囲は [3..4](次に頁が進む項目は p.4 の `remark-1-end`)。
--   `pdfPage := 7` は `lo−1 = 2` より大きく `hi = 4` より大きいので範囲外。
--   ★`sectionId` は実在し、`item` も `data-item` と一致し、頁も 1..9 の範囲内なので、
--     ★**規則 R′ が無いと 3 つの既存 G1 条項は全部これを通す**(実際、実木で
--     `FrdI#frdi-def-2-4` の `pdfPage := 51` が 4 件、こうして静かに残っていた)。
namespace Fixture
def implWithPageFarOutside : True := trivial
def implWithPageFarOutside.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Proposition 1.1", sectionId := "prop-1-1" }
end Fixture
