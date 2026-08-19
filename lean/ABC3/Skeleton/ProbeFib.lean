import ABC3.Found.Arakelov.ADegEmb
import Mathlib.FieldTheory.PrimitiveElement

namespace ABC3.Skeleton.ProbeFib

open NumberField

open scoped Classical

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

theorem card_fiber (σ : F →+* ℂ) :
    (Finset.univ.filter (fun τ : K →+* ℂ => τ.comp (algebraMap F K) = σ)).card
      = Module.finrank F K := by
  classical
  letI : Algebra F ℂ := σ.toAlgebra
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  have hequiv : {τ : K →+* ℂ // τ.comp (algebraMap F K) = σ} ≃ (K →ₐ[F] ℂ) :=
    ⟨fun t => AlgHom.mk t.1 (fun r => congrArg (fun h : F →+* ℂ => h r) t.2),
      fun g => ⟨g.toRingHom, by ext r; exact g.commutes r⟩,
      fun _ => rfl, fun _ => rfl⟩
  rw [← Fintype.card_subtype, Fintype.card_congr hequiv, AlgHom.card]

end ABC3.Skeleton.ProbeFib
