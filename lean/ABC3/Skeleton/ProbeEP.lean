import ABC3.Found.Arakelov.ADegEmb

namespace ABC3.Skeleton.ProbeEP

open AlgebraicGeometry CategoryTheory NumberField ABC3.Interface.Arakelov

theorem embPoint_comp (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
    [Algebra F K] [IsScalarTower (𝓞 F) (𝓞 K) K] (τ : K →+* ℂ) :
    embPoint K τ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
      = embPoint F (τ.comp (algebraMap F K)) := by
  show Spec.map (CommRingCat.ofHom (τ.comp (algebraMap (𝓞 K) K)))
      ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
    = Spec.map (CommRingCat.ofHom ((τ.comp (algebraMap F K)).comp (algebraMap (𝓞 F) F)))
  rw [← Spec.map_comp]
  congr 1

end ABC3.Skeleton.ProbeEP
