import Mathlib.FieldTheory.LinearDisjoint
import ABC3.Found.PGC.LubinTateTotallyRamified

/-!
# 不分岐拡大と完全分岐拡大は線型無関係

`TotallyRamified.lean::totallyRamifiedAdjoin_inf_unramifiedClosure`
(`K(α) ⊓ K^ur = K`)を、mathlib の

* `IntermediateField.LinearDisjoint.of_inf_eq_bot`
  (`[IsGalois F A]` + 有限次 + `A ⊓ B = ⊥` ⟹ 線型無関係)
* `IntermediateField.LinearDisjoint.finrank_sup`
  (`[A ⊔ B : F] = [A:F] · [B:F]`)

に渡すだけで、**合成体の次数が積になる**。

★`IsGalois` を要求されるのは**不分岐側**でよい——不分岐拡大は常に Galois
(`UnramifiedExtension.lean::exists_isGalois_isUnramifiedAdjoin`)なので、
完全分岐側には正規性を要求しない。これが向きの選択の要点。

## 何に効くか

`Γ_K ↠ 𝒪_K^× × Ẑ`(相互律の全射性)の要。
`[K(Λ_n) · K_{ur,m} : K] = (q^n − q^{n−1}) · m` が言えれば、
`Gal(K(Λ_n)·K_{ur,m}/K) ≅ (𝒪_K/π^n)^× × ℤ/m` へ進める。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- 不分岐な単項拡大は `K^ur` に含まれる。 -/
theorem adjoin_le_unramifiedClosure_of_isUnramifiedAdjoin (K : PAdicLocalField p)
    {y : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hy : IsUnramifiedAdjoin K y) :
    IntermediateField.adjoin K.carrier ({y} : Set K.closure) ≤ unramifiedClosure K := by
  rw [IntermediateField.adjoin_le_iff]
  rintro w rfl
  exact (mem_unramifiedClosure_iff_isUnramified K w).mpr hy

/-- **★★★★★★★★★★不分岐 `K(y)` と完全分岐 `K(α)` は交わらない**。 -/
theorem inf_eq_bot_of_isUnramified_of_isTotallyRamified (K : PAdicLocalField p)
    {y α : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure))]
    (hy : IsUnramifiedAdjoin K y) (hα : IsTotallyRamifiedAdjoin K α) :
    IntermediateField.adjoin K.carrier ({y} : Set K.closure)
      ⊓ IntermediateField.adjoin K.carrier ({α} : Set K.closure) = ⊥ := by
  refine le_antisymm ?_ bot_le
  calc IntermediateField.adjoin K.carrier ({y} : Set K.closure)
        ⊓ IntermediateField.adjoin K.carrier ({α} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier ({α} : Set K.closure) ⊓ unramifiedClosure K := by
        rw [inf_comm]
        exact inf_le_inf_left _ (adjoin_le_unramifiedClosure_of_isUnramifiedAdjoin K hy)
    _ = ⊥ := totallyRamifiedAdjoin_inf_unramifiedClosure K hα

/-- **★★★★★★★★★★★★★★合成体の次数は積**——
`[K(y)·K(α) : K] = [K(y):K] · [K(α):K]`(`y` 不分岐、`α` 完全分岐)。 -/
theorem finrank_sup_of_isUnramified_of_isTotallyRamified (K : PAdicLocalField p)
    {y α : K.closure}
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure))]
    [IsGalois K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))]
    (hy : IsUnramifiedAdjoin K y) (hα : IsTotallyRamifiedAdjoin K α) :
    Module.finrank K.carrier
        ((IntermediateField.adjoin K.carrier ({y} : Set K.closure))
          ⊔ (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) :
          IntermediateField K.carrier K.closure)
      = Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure))
        * Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) :=
  (IntermediateField.LinearDisjoint.of_inf_eq_bot
    (inf_eq_bot_of_isUnramified_of_isTotallyRamified K hy hα)).finrank_sup

/-- **★★★★★★★★★★★★★★★★★★`[K(Λ_n) · K_{ur,m} : K] = (q^n − q^{n−1}) · m`**。

`exists_isGalois_isUnramifiedAdjoin`(次数 `m` の不分岐拡大——Galois)と
`exists_isTotallyRamifiedAdjoin_lubinTate_psi`(次数 `q^n − q^{n−1}` の
完全分岐拡大)を `finrank_sup_of_isUnramified_of_isTotallyRamified` に渡す。

★相互律 `Γ_K^ab ≅ 𝒪_K^× × Ẑ` が要求する「`K_π` と `K^ur` が独立」の、
有限段での形。 -/
theorem exists_finrank_sup_lubinTate_unramified (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (m : ℕ) (hm : m ≠ 0) :
    ∃ y α : K.closure,
      Module.finrank K.carrier
        ((IntermediateField.adjoin K.carrier ({y} : Set K.closure))
          ⊔ (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) :
          IntermediateField K.carrier K.closure)
        = m * ((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1)) := by
  obtain ⟨y, hym, hyu, hygal⟩ := exists_isGalois_isUnramifiedAdjoin K m hm
  obtain ⟨α, hα, hαd⟩ :=
    exists_isTotallyRamifiedAdjoin_lubinTate_psi K hq hπmax hπne0 f hf0 hf1 hf n hn
  haveI := hygal
  refine ⟨y, α, ?_⟩
  rw [finrank_sup_of_isUnramified_of_isTotallyRamified K hyu hα, hym, hαd]

end ABC3.Found.PGC
