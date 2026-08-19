import ABC3.Found.GaloisRep.UniversalF2

namespace ABC3.Skeleton.ProbeDom
example : IsDomain (MvPolynomial (Fin 5) (ZMod 2)) := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  haveI : Field (ZMod 2) := inferInstance
  haveI : IsDomain (ZMod 2) := inferInstance
  infer_instance
end ABC3.Skeleton.ProbeDom
