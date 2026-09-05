-- D42 fixture: `.absent` に再実行できる検索パターンが無い(落とすべき)。
-- ★2026-08-14 の実失敗(charP_iff_prime_eq_zero を「無い」と書いた)と、
--   2026-09-05 の 4 件(ULift.field / continuousCohomology / ProfiniteCompletion /
--   CompactSpace Gal)が根拠。探索範囲を書かない「無い」は再現も反証もできない。
namespace Fixture

def someResult.status : String :=
  toString (ABC3.Meta.LeanStatus.absent "mathlib に相当する宣言は無い(実測)")

end Fixture
