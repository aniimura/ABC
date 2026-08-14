-- D7 fixture: Interface/ の structure に witness も waiting も無い(落ちるべき)
namespace Fixture
structure NoWitness where
  f : Nat → Nat
  bad : ∀ n, f n < f (n + 1) ∧ f n < 10
end Fixture
