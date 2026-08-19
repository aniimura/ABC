import ABC3.Found.GaloisRep.UniversalF2
import Mathlib.Algebra.Field.ZMod

namespace ABC3.Skeleton.ProbeDom
example : IsDomain (MvPolynomial (Fin 5) (ZMod 2)) := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  infer_instance
example : IsDomain (Polynomial (Polynomial (MvPolynomial (Fin 5) (ZMod 2)))) := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  infer_instance
end ABC3.Skeleton.ProbeDom
