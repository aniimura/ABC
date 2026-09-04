import ABC3.Found.PGC.LubinTateReciprocitySurjective
import ABC3.Found.PGC.LubinTateActionEquivariance

/-!
# `K.carrier⟮x⟯` の自己同型はすべて大域的な `σ` の制限として実現される(`sorry` 無し)

`Found/PGC/LubinTateActionEquivariance.lean::algEquivRestrictSelf`
(`σ:K.closure≃ₐ[K.carrier]K.closure`を`K.carrier⟮x⟯`自身への
自己同型に制限したもの)が**全射**であること——`K.carrier⟮x⟯`の
任意の自己同型`τ`に対し、`algEquivRestrictSelf σ = τ`となる大域的な
`σ`が存在すること——を示す。

## 証明の構造

`τ`の`x`での値`τ(x)`(`K.carrier⟮x⟯`の元、`K.closure`へ座標を
戻したもの)が`ψ_nの根`に留まることを、`Found/PGC/LubinTateReciprocity
Surjective.lean`と同じ手筋(`ψ_n`の既約性→`minpoly K.carrier x=ψ_n`
→`τ`が`K.carrier`-代数準同型であることから`Polynomial.aeval_algHom_
apply`で`τ(x)`も同じ多項式の根)で示す——ただし今回は`τ`がすでに
`K.carrier⟮x⟯`の中の自己同型なので、大域的な`Gal(K.closure/
K.carrier)`のGalois同変性を経由する必要は無い。これで`τ(x)∈ψ_nの根`
が分かれば、`exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints`
(既出)で`σ(x)=τ(x)`となる大域的な`σ`が存在し、あとは`K.carrier⟮x⟯`
が`x`で生成されること(`IntermediateField.algHom_ext_of_eq_adjoin`)
から`algEquivRestrictSelf σ=τ`が従う。

## 使い道

これで`res:=algEquivRestrictSelf`の全射性が確立され、
`reciprocityMap`が(`res`のファイバー上で一定であることと合わせて)
`Gal(K.carrier⟮x⟯/K.carrier)`上へ全単射的に降りることを示す
準備が整った——`Gal(K.carrier⟮x⟯/K.carrier)≅(𝒪_K/π^n)^×`という
古典的なLubin-Tate理論の主定理の、最後の仕上げに向けた土台。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`algEquivRestrictSelf`は全射**: `K.carrier⟮x⟯`の任意の自己同型`τ`
に対し、`algEquivRestrictSelf σ = τ`となる大域的な
`σ∈Gal(K.closure/K.carrier)`が存在する。 -/
theorem algEquivRestrictSelf_surjective
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
    Function.Surjective (algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro τ
  have hτxψ : (τ ⟨x, hmem⟩ : K.closure) ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
    haveI := valuationRing_isDVR K
    haveI := IsAlgClosure.normal K.carrier K.closure
    have hmonic := (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
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
    have hxroot : Polynomial.aeval x psiK = 0 := by
      rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, hrw, Polynomial.mem_roots'] at hxψ
      have hthis := hxψ.2
      rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] at hthis
      exact hthis
    have hinj : Function.Injective (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val :=
      RingHom.injective _
    have hLval : Polynomial.aeval (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) psiK
        = 0 := by
      apply hinj
      rw [map_zero, ← hxroot]
      exact (Polynomial.aeval_algHom_apply (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val
        (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) psiK).symm
    set φ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) →ₐ[K.carrier] K.closure :=
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).val.comp τ.toAlgHom with hφ
    have hval_eq2 := Polynomial.aeval_algHom_apply φ
      (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) psiK
    have hφapp : φ (⟨x, hmem⟩ : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) =
        (τ ⟨x, hmem⟩ : K.closure) := rfl
    rw [hφapp] at hval_eq2
    have hτxroot : Polynomial.aeval (τ ⟨x, hmem⟩ : K.closure) psiK = 0 := by
      rw [hval_eq2, hLval, map_zero]
    rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, hrw, Polynomial.mem_roots']
    refine ⟨((hmonic.map _).map _).ne_zero, ?_⟩
    show Polynomial.IsRoot _ (τ ⟨x, hmem⟩ : K.closure)
    rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]
    exact hτxroot
  obtain ⟨σ, hσ⟩ := exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn
    x (τ ⟨x, hmem⟩ : K.closure) hxψ hτxψ
  refine ⟨σ, ?_⟩
  apply AlgEquiv.coe_algHom_injective
  apply IntermediateField.algHom_ext_of_eq_adjoin K.carrier rfl
  intro z hz
  simp only [Set.mem_singleton_iff] at hz
  subst hz
  ext
  show σ z = (τ ⟨z, hmem⟩ : K.closure)
  exact hσ

end ABC3.Found.PGC
