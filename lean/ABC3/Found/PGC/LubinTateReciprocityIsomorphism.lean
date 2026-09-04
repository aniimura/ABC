import ABC3.Found.PGC.LubinTateReciprocityGaloisSurjective
import ABC3.Found.PGC.LubinTateReciprocityHomomorphism

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
Lubin-Tate 相互律の主定理: `Gal(K.carrier⟮x⟯/K.carrier) ≅ (𝒪_K/π^n)^×`(`sorry` 無し)

古典的なLubin-Tate理論の主定理——原始的な`π^n`-捩れ点`x`を1つ添加した
拡大`K.carrier⟮x⟯/K.carrier`のGalois群が`(𝒪_K/π^n𝒪_K)^×`と自然に
同型であること——を、本セッションで積み上げた全ての結果を組み合わせ
て完成させる。

## 構成

`Found/PGC/LubinTateActionBijective.lean::reciprocityMap`
(`Gal(K.closure/K.carrier)→(𝒪_K)^×⧸principalUnits K π n`)は
①well-defined(存在・一意性)②全単射(`reciprocityMap_surjective`)
③群準同型(`reciprocityMap_mul`)のすべてを満たすが、**定義域が
大きすぎる**(`K.closure`全体の絶対Galois群であり、単射ではない)。

これを`K.carrier⟮x⟯`の自己同型群まで**降ろす**ため:

1. **`res:=algEquivRestrictSelf`**(`Found/PGC/LubinTateAction
   Equivariance.lean`)は`Gal(K.closure/K.carrier)→Gal(K.carrier⟮x⟯/
   K.carrier)`の制限写像で、**全射**
   (`algEquivRestrictSelf_surjective`)かつ**群準同型**
   (`algEquivRestrictSelf_mul`)。
2. **`galoisUnitReciprocityMap`**: `reciprocityMap`と全く同じ構成
   (`existsUnique`+`Classical.choose`)を`τ∈Gal(K.carrier⟮x⟯/
   K.carrier)`に対して直接行い、`τ↦u_τ`(`u_τ·x=τ(x)`となる一意な
   単数類)を得る。`res`を経由して`reciprocityMap`と一致すること
   (`galoisUnitReciprocityMap_eq_reciprocityMap`)から、
   全射性(`reciprocityMap_surjective`+`res`の全射性)・
   単射性(`reciprocityMap_spec`+`algEquivRestrictSelf_eq_of_eq`)・
   群準同型性(`reciprocityMap_mul`+`res`の群準同型性)のすべてが
   `reciprocityMap`側から直接移植される。
3. これらを`MulEquiv.ofBijective`で束ね、**`galoisUnitReciprocity
   Equiv:Gal(K.carrier⟮x⟯/K.carrier)≃*(𝒪_K)^×⧸principalUnits K π n`**
   という全単射的群準同型(=群同型)を得る。
4. 既存の**`principalUnitsQuotientEquiv`**
   (`(𝒪_K)^×⧸principalUnits K π n≃*(𝒪_K/π^n𝒪_K)^×`、
   `Found/PGC/AdjoinIntegers.lean`)と合成し、**`galoisReciprocity
   Equiv:Gal(K.carrier⟮x⟯/K.carrier)≃*(𝒪_K/π^n𝒪_K)^×`**——
   古典的なLubin-Tate理論の主定理そのものを得る。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `σ(x)=τ(x)`ならば`algEquivRestrictSelf σ=algEquivRestrictSelf τ`
——`K.carrier⟮x⟯=K.carrier[x]`が`x`で生成されること
(`IntermediateField.algHom_ext_of_eq_adjoin`)から。 -/
theorem algEquivRestrictSelf_eq_of_eq
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ τ : K.closure ≃ₐ[K.carrier] K.closure) (h : σ x = τ x) :
    algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ =
      algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ := by
  apply AlgEquiv.coe_algHom_injective
  apply IntermediateField.algHom_ext_of_eq_adjoin K.carrier rfl
  intro z hz
  simp only [Set.mem_singleton_iff] at hz
  subst hz
  ext
  show σ z = τ z
  exact h

/-- `ψ_nの根`の任意の元`x`は、`K.carrier⟮x⟯`の任意の自己同型`τ`の
下でも`ψ_nの根`に留まる——`Found/PGC/LubinTateReciprocitySurjective.
lean::exists_algEquiv_of_mem_iteratedLubinTatePsiTorsionPoints`の
証明と同じ手筋(`ψ_n`の既約性→`minpoly K.carrier x=ψ_n`→`τ`が
`K.carrier`-代数準同型であることから同じ多項式の根)を、大域的な
`Gal(K.closure/K.carrier)`を経由せず直接`τ:K.carrier⟮x⟯≃ₐ[K.carrier]
K.carrier⟮x⟯`について行う。 -/
theorem algEquivL_mem_iteratedLubinTatePsiTorsionPoints
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
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (τ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (τ ⟨x, hmem⟩ : K.closure) ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
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

/-- `K.carrier⟮x⟯`の自己同型`τ`に対し、`U·x=τ(x)`となる単数類`U`が
ただ一つ存在する——`reciprocityMap`と全く同じ構成
(`unitActionQuotientBijOn`の全単射性)を、大域的な`σ`ではなく
`τ`について直接行う。 -/
theorem existsUnique_unitActionQuotient_eq_algEquivL
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ∃! U : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n,
      (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = (τ ⟨x, hmem⟩ : K.closure) := by
  have hτxψ : (τ ⟨x, hmem⟩ : K.closure) ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn :=
    algEquivL_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hmem τ
  obtain ⟨U, hU⟩ := (unitActionQuotientBijOn_bijective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem).2
    ⟨(τ ⟨x, hmem⟩ : K.closure), hτxψ⟩
  have hU0 := Subtype.ext_iff.mp hU
  simp only [unitActionQuotientBijOn] at hU0
  refine ⟨U, hU0, ?_⟩
  intro U'' hU''
  apply unitActionQuotientBijOn_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  apply Subtype.ext
  show (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U'') :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) =
    (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)
  rw [hU'']
  exact hU0.symm

/-- **相互写像`τ↦u_τ`**(`K.carrier⟮x⟯`の自己同型群版)——
`reciprocityMap`の`K.carrier⟮x⟯`版。 -/
noncomputable def galoisUnitReciprocityMap
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (𝒪[K.carrier])ˣ ⧸ principalUnits K π n :=
  (existsUnique_unitActionQuotient_eq_algEquivL K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).choose

theorem galoisUnitReciprocityMap_spec
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
        (galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ)) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = (τ ⟨x, hmem⟩ : K.closure) :=
  (existsUnique_unitActionQuotient_eq_algEquivL K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ).choose_spec.1

/-- `algEquivRestrictSelf`は群準同型。 -/
theorem algEquivRestrictSelf_mul
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ σ' : K.closure ≃ₐ[K.carrier] K.closure) :
    algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (σ * σ') =
      algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ *
        algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ' := by
  apply AlgEquiv.ext
  intro z
  apply Subtype.ext
  rw [coe_algEquivRestrictSelf]
  show σ (σ' (z : K.closure)) =
    ((algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ *
      algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ') z : K.closure)
  rw [AlgEquiv.mul_apply, coe_algEquivRestrictSelf, coe_algEquivRestrictSelf]

/-- `galoisUnitReciprocityMap`は`res:=algEquivRestrictSelf`を経由
すれば`reciprocityMap`に一致する——両者とも「`x`での値」だけで
一意に定まることから。 -/
theorem galoisUnitReciprocityMap_eq_reciprocityMap
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
      (algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ) =
      reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ := by
  apply (existsUnique_unitActionQuotient_eq_algEquivL K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
    (algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ)).unique
    (galoisUnitReciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
      (algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ))
  rw [coe_algEquivRestrictSelf]
  exact reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ

/-- **`galoisUnitReciprocityMap`は群準同型**——`res`の群準同型性+
`reciprocityMap`の群準同型性から。 -/
theorem galoisUnitReciprocityMap_mul
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (τ ρ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (τ * ρ) =
      galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ *
        galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem ρ := by
  obtain ⟨σ, hσ⟩ := algEquivRestrictSelf_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ
  obtain ⟨σ', hσ'⟩ := algEquivRestrictSelf_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem ρ
  have hτρ : τ * ρ = algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (σ * σ') := by
    rw [algEquivRestrictSelf_mul, hσ, hσ']
  rw [hτρ, ← hσ, ← hσ', galoisUnitReciprocityMap_eq_reciprocityMap, galoisUnitReciprocityMap_eq_reciprocityMap,
    galoisUnitReciprocityMap_eq_reciprocityMap, reciprocityMap_mul]

/-- **`galoisUnitReciprocityMap`は単射**——`res`の全射性+
`algEquivRestrictSelf_eq_of_eq`から。 -/
theorem galoisUnitReciprocityMap_injective
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
    Function.Injective (galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro τ ρ heq
  obtain ⟨σ, hσ⟩ := algEquivRestrictSelf_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ
  obtain ⟨σ', hσ'⟩ := algEquivRestrictSelf_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem ρ
  rw [← hσ, ← hσ'] at heq ⊢
  rw [galoisUnitReciprocityMap_eq_reciprocityMap, galoisUnitReciprocityMap_eq_reciprocityMap] at heq
  have hσx : σ x = σ' x := by
    have h1 := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
    have h2 := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ'
    rw [heq] at h1
    rw [← h1, h2]
  rw [algEquivRestrictSelf_eq_of_eq K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ σ' hσx]

/-- **`galoisUnitReciprocityMap`は全射**——`reciprocityMap`の全射性
から直接。 -/
theorem galoisUnitReciprocityMap_surjective
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
    Function.Surjective (galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro U
  obtain ⟨σ, hσ⟩ := reciprocityMap_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem U
  exact ⟨algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ,
    (galoisUnitReciprocityMap_eq_reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ).trans hσ⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K.carrier⟮x⟯/K.carrier) ≃* (𝒪_K)^×⧸principalUnits K π n`**
——`galoisUnitReciprocityMap`の単射性・全射性・群準同型性を
`MulEquiv.ofBijective`で束ねただけ。 -/
noncomputable def galoisUnitReciprocityEquiv
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
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≃*
      (𝒪[K.carrier])ˣ ⧸ principalUnits K π n :=
  MulEquiv.ofBijective
    (MonoidHom.mk' (galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem)
      (galoisUnitReciprocityMap_mul K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem))
    ⟨galoisUnitReciprocityMap_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem,
      galoisUnitReciprocityMap_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**古典的なLubin-Tate理論の主定理**:
**`Gal(K.carrier⟮x⟯/K.carrier) ≃* (𝒪_K/π^n𝒪_K)^×`**——
`galoisUnitReciprocityEquiv`と既存の`principalUnitsQuotientEquiv`
(`Found/PGC/AdjoinIntegers.lean`)を合成するだけ。 -/
noncomputable def galoisReciprocityEquiv
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
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≃*
      (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])))ˣ :=
  (galoisUnitReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem).trans
    (principalUnitsQuotientEquiv K hπmax n hn)

end ABC3.Found.PGC
