import ABC3.Found.PGC.TotallyRamified
import ABC3.Found.PGC.LubinTateActionWeierstrass
import ABC3.Found.PGC.LubinTateDegree

/-!
# Lubin-Tate の第一段 `K(Λ_1)/K` は完全分岐

`φ_1`(原始 `π`-捩れ点の最小多項式候補)は **`𝒪_K` 上の Eisenstein**
(`LubinTateActionWeierstrass.lean::isEisensteinAt_iteratedLubinTatePrimitive_one`)
なので、`TotallyRamified.lean::exists_isTotallyRamifiedAdjoin_of_eisenstein` に
そのまま渡せる:

**`K(Λ_1)/K` は次数 `q-1` の完全分岐拡大**

## ★★★2026-09-05: `n ≥ 2` も塔なしで出た(上の見立ては誤りだった)

ここには「`n ≥ 2` の `ψ_n` は `𝒪_{K(Λ_{n-1})}` 上の Eisenstein なので、
塔の合成(慣性次数の乗法性)が別途要る」と書いていた。**誤りだった**——
`LubinTateActionPsi.lean::isEisensteinAt_iteratedLubinTatePsi` は
**`𝒪_K`(= 一般の `A`)上**の Eisenstein 性を主張している。塔は要らない。

理由は `ψ_n(X) = ψ_1(φ_{n-1}(X))` の形を見れば明らかで、`𝔪` を法として
`φ_{n-1}(X) ≡ X^{q^{n-1}}`・`ψ_1(Y) ≡ Y^{q-1}` だから
`ψ_n(X) ≡ X^{q^{n-1}(q-1)}`、定数項は `ψ_n(0) = ψ_1(0) = π ∉ 𝔪²`。

したがって `exists_isTotallyRamifiedAdjoin_lubinTate_psi` が
**すべての `n ≥ 1`** について直ちに出る。

## 何に効くか

`TotallyRamified.lean::finrank_eq_one_of_mem_unramifiedClosure_of_le`
(「完全分岐 ∩ K^ur = K」)と合わせると、`K(Λ_1)` は `K^ur` と
**線型無関係**である。相互律 `Γ_K^ab ≅ (K^×)^∧` が要求する
`Γ_K ↠ 𝒪_K^× × Ẑ` の全射性の、最初の一段。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued Classical

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

/-- **★★★★★★★★★★★★★★`K(Λ_n)/K` は次数 `q^n − q^{n−1}` の完全分岐拡大**
(すべての `n ≥ 1`)。

`ψ_n` が `𝒪_K` 上の Eisenstein であること
(`isEisensteinAt_iteratedLubinTatePsi`)を
`exists_isTotallyRamifiedAdjoin_of_eisenstein` に渡すだけ——**塔は要らない**。

`Found/PGC/LubinTateDegree.lean::finrank_adjoin_iteratedLubinTatePsi` が
同じ次数を Galois 群の側(`Gal ≅ (𝒪_K/π^n)^×`)から出しているので、
両者は一致する。 -/
theorem exists_isTotallyRamifiedAdjoin_lubinTate_psi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) :
    ∃ α : K.closure, IsTotallyRamifiedAdjoin K α
      ∧ Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure))
        = (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) := by
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπne0 (Subtype.ext h)
  have hq1 : 1 < pp ^ ff := hq ▸ Fintype.one_lt_card
  have hdegpos : 0 < (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).natDegree := by
    rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]
    have h1 : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := by
      refine Nat.pow_lt_pow_right hq1 ?_
      omega
    omega
  obtain ⟨α, hα, hrank⟩ := exists_isTotallyRamifiedAdjoin_of_eisenstein K hπ0 hπmax _
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
    (isEisensteinAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) hdegpos
  exact ⟨α, hα, hrank.trans (natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)⟩

/-! ## 与えられた捩れ点 `x` そのものについての完全分岐性

`exists_isTotallyRamifiedAdjoin_lubinTate_psi` は「そういう `α` が**ある**」
という形で、`galoisReciprocityEquiv`(`Gal(K(x)/K) ≃* (𝒪_K/π^n)^×`)や
`normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints`(正規性)が語る
`x` と同じ元だとは言えていない。**同じ `x` について三つが揃う**形にする。 -/

/-- **★★★★★★★★★★★★★★★★`ψ_n` の根 `x` について `K(x)/K` は完全分岐**。

これで同じ `x` について
* `Gal(K(x)/K) ≃* (𝒪_K/π^n)^×`(`galoisReciprocityEquiv`)
* `[K(x):K] = q^n − q^{n−1}`(`finrank_adjoin_iteratedLubinTatePsi`)
* `K(x)/K` は正規(`normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints`)
* `K(x)/K` は完全分岐(本定理)

が揃う。 -/
theorem isTotallyRamifiedAdjoin_iteratedLubinTatePsi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])] {ff : ℕ}
    (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    IsTotallyRamifiedAdjoin K x := by
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπne0 (Subtype.ext h)
  have hroot0 : Polynomial.eval x (Polynomial.map (algebraMap 𝒪[K.carrier] K.closure)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0 := by
    classical
    have hmem' : x ∈ (Polynomial.map (algebraMap 𝒪[K.carrier] K.closure)
        (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)).roots :=
      Multiset.mem_toFinset.mp hxψ
    exact (Polynomial.mem_roots'.mp hmem').2
  refine isTotallyRamifiedAdjoin_of_eisenstein K x hπ0 hπmax _
    (isDistinguishedAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn).monic
    (isEisensteinAt_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn) ?_ ?_
  · rw [natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn,
      finrank_adjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem]
  · exact hroot0

end ABC3.Found.PGC
