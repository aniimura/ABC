-- D10 fixture: 非空虚 witness を持つ Interface(通るべき)
namespace Fixture
structure Good where
  f : Nat → Nat
  mono : ∀ n, f n < f (n + 1)
theorem Good.nonvacuous : Nonempty Good :=
  ⟨{ f := id, mono := fun n => Nat.lt_succ_self n }⟩
end Fixture
