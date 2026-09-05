import ABC3.Found.PGC.LubinTateReciprocityMapLimitKernel
import ABC3.Found.PGC.UnramifiedExtension

/-!
# `Γ_K` から `K^×` の二つの因子への全射——相互律の二つの半分

`Found/PGC/UnitsSplit.lean` で `K^× ≅ ℤ × 𝒪_K^×` と分解した。古典的局所
類体論の相互律 `Γ_K^ab ≅ (K^×)^∧` は、この二つの因子それぞれが `Γ_K` の
商として現れることを言う。本ファイルはその**二つの全射**を並べる:

* **完全分岐の半分**(節目(5)、Lubin-Tate):
  `Γ_K ↠ 𝒪_K^×`——`reciprocityMapLimit`(`CompatibleUnits` への全射)と
  `unitReductionHom` が同型(単射は `IsHausdorff`、全射は `IsAdicComplete`)
  であることを合わせる。
* **不分岐の半分**(本日の不分岐拡大理論):
  `Γ_K ↠ ℤ/n`(任意の `n ≥ 1`)——`exists_surjective_absGal_to_zmod`。

`Ẑ = lim ℤ/n` を作れば後者は `Γ_K ↠ Ẑ` になる(`Ẑ` は mathlib に不在、
2026-09-05 実測)。両者を「合わせて全体」にするには Lubin-Tate の主定理
`K^ab = K^ur · K_π` が要る——そこが次の山。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- **`𝒪_K^× ≃* CompatibleUnits`**——`unitReductionHom` は単射
(`IsHausdorff`)かつ全射(`IsAdicComplete`)。 -/
noncomputable def unitsEquivCompatibleUnits (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    (𝒪[K.carrier])ˣ ≃* CompatibleUnits K hπmax :=
  MulEquiv.ofBijective (unitReductionHom K hπmax)
    ⟨unitReductionHom_injective K hπmax, unitReductionHom_surjective K hπmax⟩

/-- **★★完全分岐の半分: `Γ_K ↠ 𝒪_K^×`**——Lubin-Tate 相互律の射影極限
(節目(5))を、実際の単数群の言葉に翻訳したもの。 -/
theorem exists_surjective_absGal_to_units (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ pp ^ ff) :
    ∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →* (𝒪[K.carrier])ˣ, Function.Surjective φ := by
  refine ⟨((unitsEquivCompatibleUnits K hπmax).symm :
      CompatibleUnits K hπmax ≃* (𝒪[K.carrier])ˣ).toMonoidHom.comp
    (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf), ?_⟩
  exact (unitsEquivCompatibleUnits K hπmax).symm.surjective.comp
    (reciprocityMapLimit_surjective K hq hπmax hπne0 f hf0 hf1 hf)

/-- **★★★相互律の二つの半分**——`K^× ≅ ℤ × 𝒪_K^×` の各因子が `Γ_K` の
商として現れる。左は完全分岐(Lubin-Tate、節目(5))、右は不分岐
(不分岐拡大理論)。 -/
theorem exists_surjective_absGal_both_halves (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ pp ^ ff) :
    (∃ φ : (K.closure ≃ₐ[K.carrier] K.closure) →* (𝒪[K.carrier])ˣ, Function.Surjective φ)
      ∧ (∀ n : ℕ, n ≠ 0 → ∃ ψ : (K.closure ≃ₐ[K.carrier] K.closure)
          →* Multiplicative (ZMod n), Function.Surjective ψ) :=
  ⟨exists_surjective_absGal_to_units K hq hπmax hπne0 f hf0 hf1 hf,
    fun n hn => exists_surjective_absGal_to_zmod K n hn⟩

end ABC3.Found.PGC
