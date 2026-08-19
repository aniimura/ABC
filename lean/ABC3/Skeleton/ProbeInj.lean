import ABC3.Found.Arakelov.ADegEmb

namespace ABC3.Skeleton.ProbeInj

open NumberField

variable {F : Type} [Field F] [NumberField F]
variable {B : Type} [CommRing B] [IsDomain B] [CharZero B]

theorem injective_of_charZero (f : (𝓞 F) →+* B) : Function.Injective f := by
  rw [injective_iff_map_eq_zero]
  intro x hx
  by_contra hne
  have hint : IsIntegral ℤ x := IsIntegral.of_finite ℤ x
  have hmin := minpoly.aeval ℤ x
  have hc0 : (minpoly ℤ x).coeff 0 ≠ 0 := by
    intro h0
    exact hne (minpoly.eq_zero_iff_coeff_zero_eq_zero (R := ℤ) hint |>.mp h0)
  sorry

end ABC3.Skeleton.ProbeInj
