import ABC3.Found.PGC.TotallyRamified
import ABC3.Found.PGC.LubinTateActionWeierstrass

/-!
# Lubin-Tate の第一段 `K(Λ_1)/K` は完全分岐

`φ_1`(原始 `π`-捩れ点の最小多項式候補)は **`𝒪_K` 上の Eisenstein**
(`LubinTateActionWeierstrass.lean::isEisensteinAt_iteratedLubinTatePrimitive_one`)
なので、`TotallyRamified.lean::exists_isTotallyRamifiedAdjoin_of_eisenstein` に
そのまま渡せる:

**`K(Λ_1)/K` は次数 `q-1` の完全分岐拡大**

`n ≥ 2` の `ψ_n` は `𝒪_{K(Λ_{n-1})}` 上の Eisenstein なので、塔の合成
(慣性次数の乗法性)が別途要る——それは次の課題。

## 何に効くか

`TotallyRamified.lean::finrank_eq_one_of_mem_unramifiedClosure_of_le`
(「完全分岐 ∩ K^ur = K」)と合わせると、`K(Λ_1)` は `K^ur` と
**線型無関係**である。相互律 `Γ_K^ab ≅ (K^×)^∧` が要求する
`Γ_K ↠ 𝒪_K^× × Ẑ` の全射性の、最初の一段。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- **★★★★★★★★`K(Λ_1)/K` は次数 `q-1` の完全分岐拡大**。 -/
theorem exists_isTotallyRamifiedAdjoin_lubinTate (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff)) :
    ∃ α : K.closure, IsTotallyRamifiedAdjoin K α
      ∧ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        = pp ^ ff - 1 := by
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπne0 (Subtype.ext h)
  have hdegpos : 0 < (iteratedLubinTatePrimitive hq hπmax hπne0 f hf0 hf1 hf 1).natDegree := by
    rw [natDegree_iteratedLubinTatePrimitive hq hπmax hπne0 f hf0 hf1 hf 1, pow_one]
    have h2 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
    omega
  obtain ⟨α, hα, hrank⟩ := exists_isTotallyRamifiedAdjoin_of_eisenstein K hπ0 hπmax _
    (monic_iteratedLubinTatePrimitive hq hπmax hπne0 f hf0 hf1 hf 1)
    (isEisensteinAt_iteratedLubinTatePrimitive_one hq hπmax hπne0 f hf0 hf1 hf) hdegpos
  refine ⟨α, hα, ?_⟩
  rw [hrank, natDegree_iteratedLubinTatePrimitive hq hπmax hπne0 f hf0 hf1 hf 1, pow_one]

end ABC3.Found.PGC
