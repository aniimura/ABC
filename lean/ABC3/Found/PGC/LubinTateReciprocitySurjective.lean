import ABC3.Found.PGC.LubinTateActionBijective

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
`reciprocityMap` の全射性(`sorry` 無し)

`Found/PGC/LubinTateActionBijective.lean::reciprocityMap`が**全射**
であること——`Gal(K.closure/K.carrier)`が`ψ_nの根`に**推移的**に
作用することの言い換え——を示す。

## 鍵となる事実

`ψ_n`は`𝒪_K`上既約(`irreducible_iteratedLubinTatePsi`、
`LubinTateActionPsi.lean`、既出)。Gaussの補題で`K.carrier`上へ
持ち上げると、`ψ_n`の根`x,y`はどちらも`minpoly K.carrier x =
minpoly K.carrier y`(ともに`ψ_n`自身、`minpoly.eq_of_irreducible_
of_monic`)を満たす——**同じ最小多項式を持つ2元**であり、mathlibの
標準的なGalois理論の補題`minpoly.exists_algEquiv_of_root'`
(`K.closure`が`K.carrier`の代数閉包であること、`IsAlgClosure.normal`
から`Normal K.carrier K.closure`が言えることを使う)により、
`σ(x)=y`となる`σ∈Gal(K.closure/K.carrier)`が直接存在する——
`AdjoinRoot`・`IsAlgClosed.lift`を手作業で組み立てる必要は無かった
(instance diamondを一切経由しない、はるかに軽い経路)。

これと`unitActionQuotientLift`の全射性(既出、`unitActionQuotient
BijOn_bijective`)を組み合わせると、任意の単数類`U`に対し
`y:=U·x`は`ψ_nの根`なので、上の事実から`σ(x)=y`となる`σ`が存在し、
一意性(`existsUnique_unitActionQuotient_eq_algEquiv`)から
`reciprocityMap σ = U`が出る——`reciprocityMap`の全射性。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **`ψ_n`の根は同じ最小多項式を持つ2元なので、Galois共役**——
`x,y`がともに`ψ_nの根`なら`σ(x)=y`となる`σ∈Gal(K.closure/K.carrier)`
が存在する。`irreducible_iteratedLubinTatePsi`(既出)をGaussの
補題で`K.carrier`上へ持ち上げ、`minpoly.exists_algEquiv_of_root'`
(mathlibの標準的なGalois理論)を適用するだけ。 -/
theorem exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x y : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hyψ : y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    ∃ σ : K.closure ≃ₐ[K.carrier] K.closure, σ x = y := by
  haveI := valuationRing_isDVR K
  haveI := IsAlgClosure.normal K.carrier K.closure
  set psiK : Polynomial K.carrier :=
    (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).map (algebraMap (𝒪[K.carrier]) K.carrier)
    with hpsiK
  have hmap : algebraMap (𝒪[K.carrier]) K.closure =
      (algebraMap K.carrier K.closure).comp (algebraMap (𝒪[K.carrier]) K.carrier) :=
    IsScalarTower.algebraMap_eq (𝒪[K.carrier]) K.carrier K.closure
  have hrw : (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) =
      psiK.map (algebraMap K.carrier K.closure) := by
    rw [hmap, ← Polynomial.map_map, hpsiK]
  have hget : ∀ z : K.closure, z ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn →
      Polynomial.aeval z psiK = 0 := by
    intro z hzψ
    rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, hrw, Polynomial.mem_roots'] at hzψ
    have hthis := hzψ.2
    rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] at hthis
    exact hthis
  have hxroot := hget x hxψ
  have hyroot := hget y hyψ
  have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
  have hirr : Irreducible psiK :=
    (Polynomial.IsPrimitive.irreducible_iff_irreducible_map_fraction_map hmonic.isPrimitive).mp
      (irreducible_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)
  have hpsiKmonic : psiK.Monic := hmonic.map _
  have hminpolyx : psiK = minpoly K.carrier x := minpoly.eq_of_irreducible_of_monic hirr hxroot hpsiKmonic
  have hy' : Polynomial.aeval y (minpoly K.carrier x) = 0 := hminpolyx ▸ hyroot
  exact minpoly.exists_algEquiv_of_root' (Algebra.IsAlgebraic.isAlgebraic x) hy'

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMap`は全射**: `Gal(K.closure/K.carrier)→(𝒪_K)^×⧸
principalUnits K π n`という写像が全射であること——`Gal(K.closure/
K.carrier)`が`ψ_nの根`に推移的に作用することの言い換え。任意の
単数類`U`に対し`y:=U·x`(`ψ_nの根`、`unitActionQuotientLift`の
全射性から)へ上の共役の存在定理を適用し、一意性から
`reciprocityMap σ = U`を得る。 -/
theorem reciprocityMap_surjective
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Function.Surjective (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro U
  set y : K.closure := (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) with hy_def
  have hyψ : y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn :=
    unitActionQuotientLift_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn
      hmem U
  obtain ⟨σ, hσ⟩ := exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn
    x y hxψ hyψ
  refine ⟨σ, (existsUnique_unitActionQuotient_eq_algEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
    σ).unique (reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ) ?_⟩
  rw [hσ]

end ABC3.Found.PGC
