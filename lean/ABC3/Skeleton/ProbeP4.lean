import ABC3.Found.GaloisRep.UniversalF2
import ABC3.Found.GaloisRep.OmegaThree

namespace ABC3.Skeleton.ProbeP4

open WeierstrassCurve Polynomial Polynomial.Bivariate ABC3.Found.GaloisRep

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

theorem key4 :
    psiComp W 4 * W.ψ₂ ^ 2
      = (W.ψ 4 * (W.ψ₂ * W.ψ 4 ^ 2
          + C (C W.a₁) * (W.ψ 5 * W.ψ 3))) * W.ψ₂ ^ 2 := by
  have hR : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2 : (2 : R[X][Y]) = 0 := by
    rw [show (2 : R[X][Y]) = C (C (2 : R)) by rw [map_ofNat, map_ofNat], hR, map_zero, map_zero]
  have hd : C W.preΨ₄ = W.ψ₂ ^ 4 + C (C W.a₁) * C W.Ψ₃ * W.ψ₂ := by
    rw [psi2_char_two, ← map_pow, ← map_mul, ← map_mul, ← map_add]
    exact congrArg C (preP4_frob_char_two W)
  have h5 : W.ψ 5 = C W.preΨ₄ * W.ψ₂ * W.ψ₂ ^ 3 - C W.Ψ₃ ^ 3 := by
    have h := psi_odd W 2
    norm_num at h
    exact h
  have h6 : W.ψ 6 * W.ψ₂ = W.ψ₂ ^ 2 * C W.Ψ₃ * W.ψ 5 - C W.Ψ₃ * (C W.preΨ₄ * W.ψ₂) ^ 2 := by
    have h := psi_even W 3
    norm_num at h
    rw [psi_four] at h
    linarith_or_polyrith
  sorry

end ABC3.Skeleton.ProbeP4
