-- D36 fixture: 宇宙注釈 `.{u, v}` を付けた structure でも witness を見つけられる(通るべき)
-- ★2026-08-17 の回帰: 宣言名の正規表現が `.{u, v}` の `.` を名前に含めてしまい、
--   G2 が `Universed..waiting` を探して偽陽性を出した。
namespace Fixture
universe u v
structure Universed.{u, v} where
  Point : Type u
  Bundle : Type v
def Universed.waiting : ABC3.Meta.WaitingFor :=
  { what := "宇宙多相な台", trackB := "Found/Fixture" }
end Fixture
