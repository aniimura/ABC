import ABC3.Found.PGC.LubinTateReciprocityLimitCompat

/-!
# 隣接するレベルで両立する `reciprocityMap` の値を持つ2つの `σ` は、
生成元の上で一致する(`sorry` 無し)

節目(5)のさらに先(古典的Lubin-Tate理論が実際に主張する
`Gal(K_π/K)≅𝒪_K^×`、`K_π:=K(Λ_∞)`)へ向けた布石。`reciprocityMapLimit`
の全射性を示すには、`v∈CompatibleUnits`を1つの大域的な`τ:K_π≃ₐ[K.
carrier]K_π`へ持ち上げる必要があり、そのために「各レベル`n`で
`reciprocityMap_surjective`から得た`σ_n`が、隣接するレベル間で
`(psiGenSeq n).pt`の上で一致する」ことを示す必要がある——本ファイルは
その核心部分を確立する。

## 鍵となる観察

`reciprocityMap_pred_eq_map_succ`(既出、`LubinTateTowerCompatible.
lean`)は「同じ`σ`について」のn跨ぎ両立性を述べている。ここでは
**異なる**`σ,σ'`について、それぞれの`reciprocityMap`の値が
(`QuotientGroup.map`で)両立する**ならば**、`σ,σ'`は生成元の上で
**一致する**ことを示す——`reciprocityMap_pred_eq_map_succ`(同じσの
n跨ぎ)・`reciprocityMap_congr`(既出、値の一致だけで決まる)・
`reciprocityMap_spec`(一意性の源、既出)を2回組み合わせるだけで、
新しい罠は無かった。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**隣接レベルで両立する`reciprocityMap`の値を持つ`σ,σ'`は、下の
生成元の上で一致する**: `σ`をレベル`m+2`(生成元`(psiGenStepResult
m (psiGenSeq m)).pt`)で評価した値を`QuotientGroup.map`でレベル
`m+1`へ落としたものが、`σ'`をレベル`m+1`(生成元`(psiGenSeq m).pt`)
で評価した値に一致するなら、`σ((psiGenSeq m).pt)=σ'((psiGenSeq
m).pt)`。`reciprocityMap_pred_eq_map_succ`で`σ`自身のn跨ぎ両立性を
経由し、`reciprocityMap_congr`で表現を`psiGenSeq`基準へ揃えてから、
`reciprocityMap_spec`(一意性の源)を2回使うだけ。 -/
theorem algEquiv_eq_of_reciprocityMap_eq_of_map_eq
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (σ σ' : K.closure ≃ₐ[K.carrier] K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
          Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure))]
    (hUU' : QuotientGroup.map (principalUnits K π (m + 1 + 1)) (principalUnits K π (m + 1)) (MonoidHom.id _)
        (principalUnits_succ_le K π (m + 1))
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) (by omega)
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hψ
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hn
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hmem σ) =
      reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ') :
    σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
      σ' ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) := by
  set g := psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m) with hg_def
  have hyψ : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
  have hyn : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
  have hymem : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
          IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)} : Set K.closure) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem
  haveI hyfd : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
      ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
        IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)} : Set K.closure)) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  have hbridge := reciprocityMap_pred_eq_map_succ K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    g.pt g.hψ g.hn g.hmem hyψ hyn hymem σ
  have hcongr := reciprocityMap_congr K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt g.hcompat hyψ hyn hymem
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
  have hUeq := hcongr.symm.trans (hbridge.symm.trans hUU')
  have h1 := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
  have h2 := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ'
  rw [hUeq] at h1
  rw [← h1, h2]

end ABC3.Found.PGC
