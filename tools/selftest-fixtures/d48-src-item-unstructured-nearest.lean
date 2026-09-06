-- D48 fixture: 構造化されていない項目を**最寄りの section** にぶら下げる(通るべき)。
-- ★D47 と対にすること。「item と data-item が違えば落とす」だけを入れると、
--   まだ構造化していない細かい Remark に**嘘の id** を書かせる逆インセンティブになる。
--   落とすのは「主張している項目の section が実在する」ときだけである。
--   pGC に `Remark 1.1.1` の section は無い(2026-09-06 実測: data-item は 33 種)。
namespace Fixture
def implOnNearestSection : True := trivial
def implOnNearestSection.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Remark 1.1.1(Proposition 1.1 に付随する注)",
    sectionId := "prop-1-1" }
end Fixture
