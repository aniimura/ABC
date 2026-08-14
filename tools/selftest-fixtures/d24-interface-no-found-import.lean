-- D24 fixture: Interface が Found を import していない(通るべき)
--   witness を持つこと自体は禁じていない——禁じているのは import の向きだけ。
--   実装を要する witness は Found 側のファイルから Interface 名前空間へ足す。
import ABC3.Skeleton.PGC.Setup
namespace Fixture
structure Decoupled where
  f : Nat → Nat
theorem Decoupled.nonvacuous : Nonempty Decoupled := ⟨{ f := id }⟩
end Fixture
