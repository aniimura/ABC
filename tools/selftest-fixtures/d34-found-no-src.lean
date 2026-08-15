-- D34 fixture: Found/ の宣言に `.src` が無い(通るべき)。
-- ★`Found/` には補助補題や witness など原典項目に対応しない宣言が多数ある。
--   `.src` を必須にすると、それらに嘘の出典を書かせる逆インセンティブが生まれる。
namespace Fixture
theorem helperLemma : True := trivial
end Fixture
