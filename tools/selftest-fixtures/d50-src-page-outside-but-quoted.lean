-- D50 fixture: `.src` の `pdfPage` は項目の頁の範囲の外だが、
-- **そのファイルがその頁の逐語を引いている**(通るべき)。★D49 と**対**である。
-- ★対で置く理由(実測。メタ第 13 回):
--   「範囲の外なら落とす」だけを入れると、実木で **6 件中 3 件(50%)が誤報**になる。
--   誤報の中身は「§1 の設定の頁を引いて Proposition 1.4 の足場にする」
--   「1 つの実装が複数頁を引く」という**正当な運用**である。
--   `原文 (タグ p.N):` を第 2 条件に足すと、その 3 件が全部消えて誤報 0 になった。
--   ★この対を外すと、正当な引用に**嘘の頁**を書かせる逆インセンティブが生まれる
--     (G1 の `Found/` 非対称・D47/D48 と同じ理由)。
namespace Fixture
/-- 原文 (pGC p.7):
> where the vertical morphisms are induced by field inclusions
-/
def implPageOutsideButQuoted : True := trivial
def implPageOutsideButQuoted.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 7, item := "Proposition 1.1", sectionId := "prop-1-1" }
end Fixture
