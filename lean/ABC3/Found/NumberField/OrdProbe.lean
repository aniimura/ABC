import ABC3.Found.Divisor.ArithDivisor
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace

open NumberField IsDedekindDomain
open scoped NumberField Classical

noncomputable def ordAt {L : Type} [Field L] [NumberField L]
    (v : HeightOneSpectrum (𝓞 L)) (x : L) : ℤ :=
  if h : v.valuation L x = 0 then 0 else -(WithZero.unzero h).toAdd

example {L : Type} [Field L] [NumberField L] (v : HeightOneSpectrum (𝓞 L)) (x : L) (hx : x ≠ 0) :
    ((FinitePlace.mk v) x : ℝ) = (Ideal.absNorm v.asIdeal : ℝ) ^ (-(ordAt v x)) := by
  sorry
