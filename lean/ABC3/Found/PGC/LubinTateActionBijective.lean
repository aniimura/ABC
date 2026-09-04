import ABC3.Found.PGC.LubinTateActionInjective
import ABC3.Found.PGC.QuotientCardinality

/-!
# `(𝒪_K)^×⧸principalUnits K π n` と `ψ_n` の根の全単射(`sorry` 無し)

★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★ **本セッションの
到達点**: `Found/PGC/LubinTateActionInjective.lean::
unitActionQuotientLift_injective`(単射性)と `Found/PGC/
QuotientCardinality.lean::card_principalUnitsQuotient`/`Found/PGC/
LubinTateDistinguishedSeparable.lean::card_iteratedLubinTatePsi
TorsionPoints`(両者とも `q^n-q^{n-1}` に一致する濃度、既出)を
組み合わせ、

  `(𝒪_K)^×⧸principalUnits K π n ≃ ψ_n の根`

という**古典的なLubin-Tate理論の核心的な事実**(原始的な
`π^n`-捩れ点の全体が `(𝒪_K/π^n)^×` と自然に一対一対応する)を、
有限集合の鳩の巣原理(単射+濃度一致⟹全単射)だけで確立する。

## 構成

`unitActionQuotientBijOn`: `unitActionQuotientLift` を(既存の
`unit_action_mem_iteratedLubinTatePsiTorsionPoints` により)
`ψ_n の根`という部分集合へ**値域を制限**した写像
(`↥(iteratedLubinTatePsiTorsionPoints ...)` への写像、`Finset`の
`Sort`強制)。単射性は `unitActionQuotientLift_injective` と
「`adjoinIntegers K x → K.carrier⟮x⟯ → K.closure` の二重coeが
単射(部分環・部分体の埋め込みの合成)」を組み合わせるだけ。

全射性は `Function.Injective.surjective_of_finite`(有限型の間の
単射で、`Nonempty (α≃β)` が分かっていれば全射になる、という
mathlibの一般論)を、`Fintype.card_eq.mp`(濃度が等しければ
`Nonempty (α≃β)`)で得た同型の証拠に適用するだけ——`Nat.card`
(`card_principalUnitsQuotient`)と`Finset.card`
(`card_iteratedLubinTatePsiTorsionPoints`、`Fintype.card_coe`
経由)を `Nat.card_eq_fintype_card` で橋渡しする。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `unitActionQuotientLift` の値域を `ψ_n の根` へ制限した写像
——`unit_action_mem_iteratedLubinTatePsiTorsionPoints` による
値域の制限を `QuotientGroup.induction_on` で一般の `U` へ持ち上げる
だけ。 -/
theorem unitActionQuotientLift_mem_iteratedLubinTatePsiTorsionPoints
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
    (U : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n) :
    (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn := by
  induction U using QuotientGroup.induction_on with
  | H u =>
    rw [unitActionQuotientLift_mk]
    exact unit_action_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x
      hxψ hmem hxn u

/-- **`unitActionQuotientLift` の値域を `ψ_n の根` へ制限した写像**。 -/
noncomputable def unitActionQuotientBijOn
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
    ((𝒪[K.carrier])ˣ ⧸ principalUnits K π n) →
      (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure) :=
  fun U => ⟨_, unitActionQuotientLift_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf
    n hn x hxψ hxn hmem U⟩

/-- `unitActionQuotientBijOn` の単射性——`unitActionQuotientLift_
injective` と、`adjoinIntegers K x → K.carrier⟮x⟯ → K.closure` の
二重coeが単射(`Subtype.coe_injective` を2回)であることから。 -/
theorem unitActionQuotientBijOn_injective
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
    Function.Injective (unitActionQuotientBijOn K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  intro U1 U2 heq
  apply unitActionQuotientLift_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  have h0 := Subtype.ext_iff.mp heq
  simp only [unitActionQuotientBijOn] at h0
  have h1 : (↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U1) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) =
      (↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U2) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    Subtype.coe_injective h0
  exact Subtype.coe_injective h1

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**古典的な
Lubin-Tate理論の核心的な全単射**: `(𝒪_K)^×⧸principalUnits K π n ≃
ψ_n の根`。単射性(`unitActionQuotientBijOn_injective`)+濃度の
一致(`card_principalUnitsQuotient`=`card_iteratedLubinTatePsiTorsion
Points`、ともに`q^n-q^{n-1}`)から、有限集合の間の単射は自動的に
全単射(`Function.Injective.surjective_of_finite`)。 -/
theorem unitActionQuotientBijOn_bijective
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
    Function.Bijective (unitActionQuotientBijOn K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) := by
  have hinj := unitActionQuotientBijOn_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  haveI := finite_principalUnitsQuotient K hπmax n hn
  haveI : Fintype ((𝒪[K.carrier])ˣ ⧸ principalUnits K π n) := Fintype.ofFinite _
  have hcard : Fintype.card ((𝒪[K.carrier])ˣ ⧸ principalUnits K π n) =
      Fintype.card (iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn : Finset K.closure) := by
    rw [Fintype.card_coe, ← Nat.card_eq_fintype_card,
      card_principalUnitsQuotient K hq hπmax hπne0 n hn,
      card_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn]
  obtain ⟨e⟩ := Fintype.card_eq.mp hcard
  exact ⟨hinj, Function.Injective.surjective_of_finite e hinj⟩

end ABC3.Found.PGC
