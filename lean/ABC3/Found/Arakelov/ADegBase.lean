import ABC3.Found.Arakelov.ADegEmb
import Mathlib.FieldTheory.PrimitiveElement

/-!
# Arakelov (D2) 第 301 ブロック —— **★★★★★★★D2 —— `APicSpecData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★★★★★★底変換不変性——**埋め込みの数え上げ**で出る

    Σ_{τ : K →+* ℂ} f(τ|_F) = [K:F] · Σ_{σ : F →+* ℂ} f(σ)

★各 `σ` の上のファイバーは `K →ₐ[F] ℂ`(ℂ を `σ` で `F`-代数と見る)であり、
その濃度は mathlib の **`AlgHom.card`** で `[K:F]` である。

★★`[K:ℚ] = [F:ℚ]·[K:F]` で正規化の分母がちょうど約分される。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `embPoint_comp` | ★複素点の底変換 |
| `card_fiber` | ★★★ファイバーの濃度は `[K:F]`(`AlgHom.card`) |
| `sum_emb_baseChange` | ★★★★埋め込みにわたる和の底変換 |
| `degFOf_baseChange` | ★★★★★**次数の底変換不変性** |
| `APicSpecData.nonvacuous` | ★★★★★★★**D2 達成** |
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

open scoped Classical

attribute [local instance] aPicGroup

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-- ★**複素点の底変換**——`𝓞_F → 𝓞_K` に沿って合成すると埋め込みの制限になる。 -/
theorem embPoint_comp [IsScalarTower (𝓞 F) (𝓞 K) K] (τ : K →+* ℂ) :
    embPoint K τ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
      = embPoint F (τ.comp (algebraMap F K)) := by
  show Spec.map (CommRingCat.ofHom (τ.comp (algebraMap (𝓞 K) K)))
      ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))
    = Spec.map (CommRingCat.ofHom ((τ.comp (algebraMap F K)).comp (algebraMap (𝓞 F) F)))
  rw [← Spec.map_comp]
  congr 1

/-- ★★★**ファイバーの濃度は `[K:F]`**——mathlib の `AlgHom.card` そのもの。 -/
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

/-- ★★★★**埋め込みにわたる和の底変換**。 -/
theorem sum_emb_baseChange (f : (F →+* ℂ) → ℝ) :
    ∑ τ : K →+* ℂ, f (τ.comp (algebraMap F K))
      = (Module.finrank F K : ℝ) * ∑ σ : F →+* ℂ, f σ := by
  classical
  have h := Finset.sum_fiberwise_of_maps_to
      (s := (Finset.univ : Finset (K →+* ℂ))) (t := (Finset.univ : Finset (F →+* ℂ)))
      (g := fun τ : K →+* ℂ => τ.comp (algebraMap F K))
      (f := fun τ : K →+* ℂ => f (τ.comp (algebraMap F K)))
      (fun τ _ => Finset.mem_univ _)
  rw [← h, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun σ _ => ?_)
  have hin : ∀ τ ∈ Finset.univ.filter (fun τ : K →+* ℂ => τ.comp (algebraMap F K) = σ),
      f (τ.comp (algebraMap F K)) = f σ := by
    intro τ hτ
    rw [(Finset.mem_filter.1 hτ).2]
  rw [Finset.sum_congr rfl hin, Finset.sum_const, card_fiber F K σ, nsmul_eq_mul]

/-- ★★★★★**次数の底変換不変性**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★`[K:ℚ] = [F:ℚ]·[K:F]` が正規化の分母をちょうど約分する。 -/
theorem degFOf_baseChange [IsScalarTower (𝓞 F) (𝓞 K) K]
    (x : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F)))
      × Multiplicative (arcCM (Spec (CommRingCat.of (𝓞 F))))) :
    degFOf K (aPicDataWitness.pullback
        (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) x) = degFOf F x := by
  haveI : FiniteDimensional ℚ F := inferInstance
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  have htower : Module.finrank ℚ F * Module.finrank F K = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  have hF : (Module.finrank ℚ F : ℝ) ≠ 0 := by
    have := Module.finrank_pos (R := ℚ) (M := F); positivity
  have hFK : (Module.finrank F K : ℝ) ≠ 0 := by
    have := Module.finrank_pos (R := F) (M := K); positivity
  have hsum : (∑ τ : K →+* ℂ,
        (x.2 : arcCM (Spec (CommRingCat.of (𝓞 F)))).fn
          (embPoint K τ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))))
      = (Module.finrank F K : ℝ)
        * ∑ σ : F →+* ℂ, (x.2 : arcCM _).fn (embPoint F σ) := by
    rw [← sum_emb_baseChange F K (fun σ => (x.2 : arcCM _).fn (embPoint F σ))]
    exact Finset.sum_congr rfl (fun τ _ => congrArg _ (embPoint_comp F K τ))
  show -(∑ τ : K →+* ℂ, (x.2 : arcCM _).fn
      (embPoint K τ ≫ Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))))
      / (Module.finrank ℚ K : ℝ) = _
  rw [hsum, degFOf, ← htower]
  push_cast
  field_simp

/-- ★★★★★★★**`APicSpecData` の実装**。 -/
noncomputable def aPicSpecDataWitness : APicSpecData where
  toAPicData := aPicDataWitness
  degF := fun F _ _ x => degFOf F x
  degF_mul := fun F _ _ x y => degFOf_mul F x y
  degF_baseChange := fun F K _ _ _ _ _ _ x => degFOf_baseChange F K x
  degF_scale := fun F _ _ L m c => degFOf_scale F L m c

/-- ★★★★★★★**`APicSpecData` は非空虚である**——D2 達成。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★これが Arakelov 理論の 9 件のうち **D2** である。 -/
theorem APicSpecData.nonvacuous : Nonempty APicSpecData :=
  ⟨aPicSpecDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def aPicSpecDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(正規化次数 deg_F とその底変換不変性)",
    sectionId := "genell-def-1-1-ii" }

end ABC3.Interface.Arakelov
