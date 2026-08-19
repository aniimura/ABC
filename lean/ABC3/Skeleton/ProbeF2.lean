import ABC3.Found.GaloisRep.PsiRec

universe u

namespace ABC3.Skeleton.ProbeF2

open WeierstrassCurve Polynomial Polynomial.Bivariate ABC3.Found.GaloisRep

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★`psiComp 1 = ψ₂`。 -/
theorem psiComp_one : psiComp W 1 = W.ψ₂ := by
  rw [psiComp]
  exact complEDS₂_one _ _ _

/-- ★`ψ 1 = 1`。 -/
theorem psi_one : W.ψ 1 = 1 := by
  rw [WeierstrassCurve.ψ]
  exact normEDS_one _ _ _

end ABC3.Skeleton.ProbeF2

namespace ABC3.Skeleton.ProbeF2
open WeierstrassCurve Polynomial Polynomial.Bivariate ABC3.Found.GaloisRep

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★★標数 2 では `omegaNum 1 = 0`。 -/
theorem omegaNum_one_char_two : omegaNum W 1 = 0 := by
  rw [omegaNum_one, psi2_char_two, phi_one]
  simp
  ring

end ABC3.Skeleton.ProbeF2
