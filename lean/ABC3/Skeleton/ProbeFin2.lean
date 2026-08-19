import ABC3.Found.GaloisRep.TwoTorsion

namespace ABC3.Skeleton.ProbeFin2

open WeierstrassCurve Polynomial WeierstrassCurve.Affine ABC3.Found.GaloisRep

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- 2-捩れでは `y` は `x` から決まる。 -/
theorem two_torsion_y_unique (h2 : (2 : F) ≠ 0) {x y : F}
    (hy : y = W.toAffine.negY x y) : y = -(W.a₁ * x + W.a₃) / 2 := by
  rw [WeierstrassCurve.Affine.negY] at hy
  field_simp
  linear_combination -hy

example (h4 : (4 : F) ≠ 0) : {x : F | W.Ψ₂Sq.IsRoot x}.Finite :=
  Polynomial.finite_setOf_isRoot (W.Ψ₂Sq_ne_zero h4)

end ABC3.Skeleton.ProbeFin2
