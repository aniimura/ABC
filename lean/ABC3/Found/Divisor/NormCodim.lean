/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormBFunctor
import Mathlib.RingTheory.Ideal.GoingDown
import Mathlib.RingTheory.Ideal.Height
import Mathlib.RingTheory.IntegralClosure.GoingDown
import Mathlib.AlgebraicGeometry.Morphisms.SchemeTheoreticallyDominant

/-!
# 支配射で余次元は上がらない(鎖 `normalize` の `codim-preserved`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> that map into [possibly subvarieties of codimension ≥1 of] prime divisors of DK.

## ★★★going-down 1 本に帰着する

`g : V[M] → V[L]` は**整射**であり、`V[L]` は**正規**(`normObj_isNormalScheme`)、
`g` は**支配的**なので切断への射は**単射**である。
★したがって mathlib の
`instance [IsDomain S] [FaithfulSMul R S] [Algebra.IsIntegral R S] [IsIntegrallyClosed R] :
Algebra.HasGoingDown R S` がアフィン開の切断に当たり、
★★素イデアルの**高さが上がらない**:

  `(Q.under R).height ≤ Q.height`

これは `Ideal.exists_ltSeries_of_hasGoingDown`(鎖の持ち上げ)そのものである。

★★★茎は切断の局所化(`IsAffineOpen.isLocalization_stalk`)なので、
`IsLocalization.AtPrime.ringKrullDim_eq_height` で
`ringKrullDim (𝒪_{V[L],g(x)}) ≤ ringKrullDim (𝒪_{V[M],x})` になる。

★素イデアルの対応は `IsAffineOpen.comap_primeIdealOf_appLE`(mathlib)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `height_under_le_of_hasGoingDown` | ★★**going-down ⟹ 高さは上がらない** |
| `ringKrullDim_stalk_normMap_le` | ★★★★**`dim 𝒪_{V[L],g(x)} ≤ dim 𝒪_{V[M],x}`** |
| `normFFUnits_mem_BSubgroup'` | ★★★★★★**`B(L) → B(M)`(仮定なし)** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

/-! ## ★1. going-down ⟹ 高さは上がらない -/

/-- ★★**going-down があれば素イデアルの高さは上がらない**。

★中身は `Ideal.exists_ltSeries_of_hasGoingDown`(鎖の持ち上げ)1 本。 -/
theorem height_under_le_of_hasGoingDown {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    [Algebra.HasGoingDown R S] (Q : Ideal S) [hQ : Q.IsPrime] :
    (Q.under R).height ≤ Q.height := by
  haveI : (Q.under R).IsPrime := Ideal.IsPrime.under R Q
  rw [PrimeSpectrum.height_eq_orderHeight (⟨Q.under R, ‹(Q.under R).IsPrime›⟩ : PrimeSpectrum R),
    PrimeSpectrum.height_eq_orderHeight (⟨Q, hQ⟩ : PrimeSpectrum S), Order.height_le_iff']
  intro l hl
  haveI : Q.LiesOver l.last.asIdeal := ⟨by rw [hl]⟩
  obtain ⟨L, hlen, hlast, -⟩ := Ideal.exists_ltSeries_of_hasGoingDown l Q
  rw [← hlen, ← hlast]
  exact Order.length_le_height_last

/-- ★★**整拡大では素イデアルの高さは下がらない**(incomparability)。

★中身は `Ideal.comap_lt_comap_of_integral_mem_sdiff`(比較不能性)1 本 ——
`S` の鎖を縮約すると `R` の**狭義**の鎖になる。 -/
theorem height_le_height_under_of_isIntegral {R S : Type*} [CommRing R] [CommRing S]
    [Algebra R S] [Algebra.IsIntegral R S] (Q : Ideal S) [hQ : Q.IsPrime] :
    Q.height ≤ (Q.under R).height := by
  haveI : (Q.under R).IsPrime := Ideal.IsPrime.under R Q
  rw [PrimeSpectrum.height_eq_orderHeight (⟨Q, hQ⟩ : PrimeSpectrum S),
    PrimeSpectrum.height_eq_orderHeight (⟨Q.under R, ‹(Q.under R).IsPrime›⟩ : PrimeSpectrum R),
    Order.height_le_iff']
  intro l hl
  have hmono : StrictMono (fun i => PrimeSpectrum.comap (algebraMap R S) (l i)) := by
    intro i j hij
    have hlt : l i < l j := l.strictMono hij
    haveI : (l i).asIdeal.IsPrime := (l i).isPrime
    have hnle : ¬ ((l j).asIdeal ≤ (l i).asIdeal) := by
      intro h
      exact (ne_of_lt hlt) (le_antisymm hlt.le h)
    obtain ⟨x, hxJ, hxI⟩ := SetLike.not_le_iff_exists.mp hnle
    exact Ideal.comap_lt_comap_of_integral_mem_sdiff hlt.le ⟨hxJ, hxI⟩
      (Algebra.IsIntegral.isIntegral x)
  have hlast : RelSeries.last (LTSeries.mk l.length _ hmono)
      = (⟨Q.under R, ‹(Q.under R).IsPrime›⟩ : PrimeSpectrum R) := by
    show PrimeSpectrum.comap (algebraMap R S) l.last = _
    rw [hl]
    rfl
  rw [← hlast]
  exact Order.length_le_height_last

/-! ## ★2. 茎の次元は上がらない -/

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

set_option maxHeartbeats 1000000 in
/-- ★★★★**支配射で茎の次元は上がらない**。 -/
theorem ringKrullDim_stalk_normMap_le {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (x : normObj V M) :
    ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base x))
      ≤ ringKrullDim ((normObj V M).presheaf.stalk x) := by
  obtain ⟨_, ⟨U, hUaff, rfl⟩, hxU, -⟩ := V.isBasis_affineOpens.exists_subset_of_mem_open
    (Set.mem_univ ((normDown V M).base x)) isOpen_univ
  have hUne : (U : Set V).Nonempty := ⟨_, hxU⟩
  have hULaff : IsAffineOpen (normDown V L ⁻¹ᵁ U) := hUaff.preimage (normDown V L)
  have hUMeq : (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U) = normDown V M ⁻¹ᵁ U := by
    rw [← Scheme.Hom.comp_preimage, normMap_normDown]
  have hUMaff : IsAffineOpen ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) := by
    rw [hUMeq]; exact hUaff.preimage (normDown V M)
  have hxUM : x ∈ (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U) := by rw [hUMeq]; exact hxU
  have hgxUL : (normMap V f).base x ∈ normDown V L ⁻¹ᵁ U := hxUM
  haveI : Nonempty (normDown V L ⁻¹ᵁ U) := ⟨⟨(normMap V f).base x, hgxUL⟩⟩
  haveI : Nonempty ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) := ⟨⟨x, hxUM⟩⟩
  letI : Algebra Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) :=
    ((normMap V f).app (normDown V L ⁻¹ᵁ U)).hom.toAlgebra
  haveI : IsIntegrallyClosed Γ(normObj V L, normDown V L ⁻¹ᵁ U) :=
    isIntegrallyClosed_normSections V L hUaff hUne
  haveI : Algebra.IsIntegral Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) :=
    ⟨(normMap V f).isIntegral_app _ hULaff⟩
  haveI : FaithfulSMul Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) := by
    rw [faithfulSMul_iff_algebraMap_injective]
    haveI := IsSchemeTheoreticallyDominant.of_isDominant (normMap V f)
    exact (normMap V f).app_injective _
  letI : Algebra Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      ((normObj V L).presheaf.stalk
        ((⟨(normMap V f).base x, hgxUL⟩ : normDown V L ⁻¹ᵁ U) : normObj V L)) :=
    TopCat.Presheaf.algebra_section_stalk (normObj V L).presheaf ⟨_, hgxUL⟩
  haveI := hULaff.isLocalization_stalk ⟨(normMap V f).base x, hgxUL⟩
  letI : Algebra Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U))
      ((normObj V M).presheaf.stalk
        ((⟨x, hxUM⟩ : (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) : normObj V M)) :=
    TopCat.Presheaf.algebra_section_stalk (normObj V M).presheaf ⟨x, hxUM⟩
  haveI := hUMaff.isLocalization_stalk ⟨x, hxUM⟩
  rw [IsLocalization.AtPrime.ringKrullDim_eq_height
      (hULaff.primeIdealOf ⟨(normMap V f).base x, hgxUL⟩).asIdeal _,
    IsLocalization.AtPrime.ringKrullDim_eq_height
      (hUMaff.primeIdealOf ⟨x, hxUM⟩).asIdeal _]
  have hunder : (hUMaff.primeIdealOf ⟨x, hxUM⟩).asIdeal.under
        Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      = (hULaff.primeIdealOf ⟨(normMap V f).base x, hgxUL⟩).asIdeal := by
    have h := IsAffineOpen.comap_primeIdealOf_appLE (f := normMap V f) (x := x)
      (normDown V L ⁻¹ᵁ U) hULaff ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) hUMaff le_rfl hxUM
    have happ : ((normMap V f).appLE (normDown V L ⁻¹ᵁ U)
        ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) le_rfl)
        = (normMap V f).app (normDown V L ⁻¹ᵁ U) := (Scheme.Hom.app_eq_appLE _).symm
    rw [happ] at h
    exact congrArg PrimeSpectrum.asIdeal h
  rw [← hunder]
  exact_mod_cast height_under_le_of_hasGoingDown _

set_option maxHeartbeats 1000000 in
/-- ★★★★**整射で茎の次元は下がらない**(比較不能性)。 -/
theorem ringKrullDim_stalk_le_normMap {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (x : normObj V M) :
    ringKrullDim ((normObj V M).presheaf.stalk x)
      ≤ ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base x)) := by
  obtain ⟨_, ⟨U, hUaff, rfl⟩, hxU, -⟩ := V.isBasis_affineOpens.exists_subset_of_mem_open
    (Set.mem_univ ((normDown V M).base x)) isOpen_univ
  have hUne : (U : Set V).Nonempty := ⟨_, hxU⟩
  have hULaff : IsAffineOpen (normDown V L ⁻¹ᵁ U) := hUaff.preimage (normDown V L)
  have hUMeq : (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U) = normDown V M ⁻¹ᵁ U := by
    rw [← Scheme.Hom.comp_preimage, normMap_normDown]
  have hUMaff : IsAffineOpen ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) := by
    rw [hUMeq]; exact hUaff.preimage (normDown V M)
  have hxUM : x ∈ (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U) := by rw [hUMeq]; exact hxU
  have hgxUL : (normMap V f).base x ∈ normDown V L ⁻¹ᵁ U := hxUM
  haveI : Nonempty (normDown V L ⁻¹ᵁ U) := ⟨⟨(normMap V f).base x, hgxUL⟩⟩
  haveI : Nonempty ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) := ⟨⟨x, hxUM⟩⟩
  letI : Algebra Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) :=
    ((normMap V f).app (normDown V L ⁻¹ᵁ U)).hom.toAlgebra
  haveI : Algebra.IsIntegral Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) :=
    ⟨(normMap V f).isIntegral_app _ hULaff⟩
  letI : Algebra Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      ((normObj V L).presheaf.stalk
        ((⟨(normMap V f).base x, hgxUL⟩ : normDown V L ⁻¹ᵁ U) : normObj V L)) :=
    TopCat.Presheaf.algebra_section_stalk (normObj V L).presheaf ⟨_, hgxUL⟩
  haveI := hULaff.isLocalization_stalk ⟨(normMap V f).base x, hgxUL⟩
  letI : Algebra Γ(normObj V M, (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U))
      ((normObj V M).presheaf.stalk
        ((⟨x, hxUM⟩ : (normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) : normObj V M)) :=
    TopCat.Presheaf.algebra_section_stalk (normObj V M).presheaf ⟨x, hxUM⟩
  haveI := hUMaff.isLocalization_stalk ⟨x, hxUM⟩
  rw [IsLocalization.AtPrime.ringKrullDim_eq_height
      (hUMaff.primeIdealOf ⟨x, hxUM⟩).asIdeal _,
    IsLocalization.AtPrime.ringKrullDim_eq_height
      (hULaff.primeIdealOf ⟨(normMap V f).base x, hgxUL⟩).asIdeal _]
  have hunder : (hUMaff.primeIdealOf ⟨x, hxUM⟩).asIdeal.under
        Γ(normObj V L, normDown V L ⁻¹ᵁ U)
      = (hULaff.primeIdealOf ⟨(normMap V f).base x, hgxUL⟩).asIdeal := by
    have h := IsAffineOpen.comap_primeIdealOf_appLE (f := normMap V f) (x := x)
      (normDown V L ⁻¹ᵁ U) hULaff ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) hUMaff le_rfl hxUM
    have happ : ((normMap V f).appLE (normDown V L ⁻¹ᵁ U)
        ((normMap V f) ⁻¹ᵁ (normDown V L ⁻¹ᵁ U)) le_rfl)
        = (normMap V f).app (normDown V L ⁻¹ᵁ U) := (Scheme.Hom.app_eq_appLE _).symm
    rw [happ] at h
    exact congrArg PrimeSpectrum.asIdeal h
  rw [← hunder]
  exact_mod_cast height_le_height_under_of_isIntegral _

/-! ## ★2.5. 素因子の上には素因子が在る -/

/-- ★**`V[M] → V[L]` は全射** —— 整射は閉写像であり、支配的だから。 -/
instance normMap_surjective {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Surjective (normMap V f) :=
  surjective_of_isDominant_of_isClosed_range (normMap V f)
    (by simpa using (normMap V f).isClosedMap _ isClosed_univ)

/-- ★★★**素因子の上には素因子が在る** —— 次元が上がりも下がりもしないから。 -/
theorem exists_primeDivisorPt_over {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (s : PrimeDivisorPt (normObj V L)) :
    ∃ w : PrimeDivisorPt (normObj V M), (normMap V f).base w.1 = s.1 := by
  obtain ⟨w0, hw0⟩ := (normMap V f).surjective s.1
  refine ⟨⟨w0, ?_⟩, hw0⟩
  have heq : ringKrullDim ((normObj V M).presheaf.stalk w0)
      = ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base w0)) :=
    le_antisymm (ringKrullDim_stalk_le_normMap V f w0) (ringKrullDim_stalk_normMap_le V f w0)
  show ringKrullDim ((normObj V M).presheaf.stalk w0) = 1
  rw [heq, hw0]
  exact s.2

/-! ## ★3. `B(L) → B(M)`(仮定なし) -/

/-- ★★★★★★**`B(L) → B(M)`** —— 仮定なし。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported -/
theorem normFFUnits_mem_BSubgroup' (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    (u : ((normObj V L).functionField)ˣ)
    (hu : u ∈ BSubgroup V DK L (normObj_isNormalScheme V L)) :
    normFFUnits V f u ∈ BSubgroup V DK M (normObj_isNormalScheme V M) :=
  normFFUnits_mem_BSubgroup V DK f
    (fun y => y.2 ▸ ringKrullDim_stalk_normMap_le V f y.1) u hu

/-- ★★★★★★**`B(L) → B(M)` を群の射として**。 -/
noncomputable def BSubgroupMap (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)] :
    BSubgroup V DK L (normObj_isNormalScheme V L)
      →* BSubgroup V DK M (normObj_isNormalScheme V M) where
  toFun u := ⟨normFFUnits V f (u : ((normObj V L).functionField)ˣ),
    normFFUnits_mem_BSubgroup' V DK f _ u.2⟩
  map_one' := by
    refine Subtype.ext ?_
    show normFFUnits V f 1 = 1
    exact map_one _
  map_mul' a b := by
    refine Subtype.ext ?_
    show normFFUnits V f ((a : ((normObj V L).functionField)ˣ) * b)
      = normFFUnits V f (a : ((normObj V L).functionField)ˣ)
        * normFFUnits V f (b : ((normObj V L).functionField)ˣ)
    exact map_mul _ _ _

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Example 6.1` の `B(L)` の関手性(仮定なし)。 -/
def BSubgroupMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — B(L) の関手性(支配射で余次元は上がらない)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
