import ABC3.Found.PGC.LubinTateReciprocityIsomorphism
import ABC3.Found.PGC.LubinTateNormalAdjoin

/-!
# `[K(Λ_n) : K] = q^n − q^{n−1}`

Lubin-Tate 理論の主定理 `galoisReciprocityEquiv`
(`Gal(K(x)/K) ≃* (𝒪_K/π^n)^×`、`Found/PGC/LubinTateReciprocityIsomorphism.lean`)
から**次数**を読む。

三つを繋ぐだけ:

1. `normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints`(`K(x)/K` は正規)
   ——標数 0 なので分離的、したがって Galois
2. `IsGalois.card_aut_eq_finrank`(`|Gal| = [K(x):K]`)
3. `galoisReciprocityEquiv` + `card_units_quotient_span_pi_pow`
   (`|(𝒪_K/π^n)^×| = q^n − q^{n−1}`、`Found/PGC/QuotientCardinality.lean`)

★これで `K_π,n/K` が「次数 `q^{n−1}(q−1)` の Galois 拡大で Galois 群が
`(𝒪_K/π^n)^×`」という古典的な形で揃った。`n = 1` の場合は
`Found/PGC/LubinTateTotallyRamified.lean::exists_isTotallyRamifiedAdjoin_lubinTate`
が**完全分岐**であることまで示している(`φ_1` が `𝒪_K` 上 Eisenstein だから)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **★★★★★★★★★★★★★★★★★★★★★★`[K(Λ_n) : K] = q^n − q^{n−1}`**。 -/
theorem finrank_adjoin_iteratedLubinTatePsi
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x
      hxψ hxn hmem
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  have hcard : Nat.card ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsGalois.card_aut_eq_finrank K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  rw [← hcard, Nat.card_congr (galoisReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x
    hxψ hxn hmem).toEquiv]
  exact card_units_quotient_span_pi_pow K hq hπmax hπne0 n hn

end ABC3.Found.PGC
