-- D23 fixture: Interface から Found を import している(落ちるべき)
--   これを許すと Skeleton が Found を推移的に引き、
--   「実装が無くても statement を書ける」という2トラック構成の要点が壊れる。
import ABC3.Found.PGC.LocalFieldNorm
namespace Fixture
structure Coupled where
  f : Nat → Nat
theorem Coupled.nonvacuous : Nonempty Coupled := ⟨{ f := id }⟩
end Fixture
