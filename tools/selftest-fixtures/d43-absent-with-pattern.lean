-- D43 fixture: `.absent` に再実行できるパターンと測定日がある(通るべき)。
-- ★この形なら `node tools/absent-recheck.mjs` が索引に当て直せる——
--   「無い」が覆ったときに、書き直しを始める前に鳴らせる。
namespace Fixture

def someResult.status : String :=
  toString (ABC3.Meta.LeanStatus.absent
    "mathlib 全体を re:`ZzzNoSuchDeclZzz`→0 で検索して 0 件(2026-09-05 実測)")

end Fixture
