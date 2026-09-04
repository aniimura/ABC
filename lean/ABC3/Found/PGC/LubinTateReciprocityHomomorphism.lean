import ABC3.Found.PGC.LubinTateActionEquivariance

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
`reciprocityMap` の準同型性(`sorry` 無し)——本セッションの集大成

`Found/PGC/LubinTateActionBijective.lean::reciprocityMap`(各`σ∈
Gal(K.closure/K.carrier)`に`u_σ·x=σ(x)`を満たす単数の類`u_σ`を
対応させる写像)が**群準同型**であること:

  `reciprocityMap (σ*τ) = reciprocityMap σ * reciprocityMap τ`

を示す。これで`Gal(K.closure/K.carrier)→(𝒪_K)^×⧸principalUnits
K π n`という写像が、集合論的な全単射(実質的には`K(Λ_n)`上への
制限を経由して)だけでなく**準同型**でもあることが確立され、
古典的なLubin-Tate理論の主定理`Gal(K(Λ_n)/K)≅(𝒪_K/π^n)^×`
(`principalUnitsQuotientEquiv`と合わせれば右辺も得られる)へ向けた
数学的な道具がすべて揃った。

## 証明の構造

`u_σ,u_τ`(`reciprocityMap σ,τ`の単数代表元)を取り、
`(σ*τ)(x)=σ(τ(x))=σ(u_τ·x)`(`reciprocityMap_spec`)に
**Galois同変性**(`algEquiv_lubinTateActionAtTorsionPoint_comm`、
`Found/PGC/LubinTateActionEquivariance.lean`)を`a:=u_τ`で適用すると
`σ(u_τ·x)="u_τ·σ(x)"`(`lubinTateActionAtAlgEquivPoint`、`x`自身の
座標系での評価)。ここで`σ(x)=u_σ·x`(`reciprocityMap_spec`)に対応
する`adjoinIntegersRestrictSelfAlgHom σ`の値と`u_σ·x`
(`lubinTateActionAtTorsionPoint`)が(座標が一致するので)**同じ
`adjoinIntegers K x`の元**であることを`lubinTateEvalAtPoint_congr`
(既出、まさにこの「同じ値・異なる証明項」を吸収するための補題)で
処理し、`lubinTateAction_mul`(乗法性、既出)で
`"u_τ·σ(x)"=(u_τ*u_σ)·x`まで計算する。最後に一意性
(`existsUnique_unitActionQuotient_eq_algEquiv`)と`𝒪_K`の可換性
(`u_τ*u_σ=u_σ*u_τ`)から結論する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMap` は群準同型**: `reciprocityMap (σ*τ) =
reciprocityMap σ * reciprocityMap τ`。 -/
theorem reciprocityMap_mul
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
    (σ τ : K.closure ≃ₐ[K.carrier] K.closure) :
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (σ * τ) =
      reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ *
        reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ := by
  obtain ⟨uσ, huσ⟩ := QuotientGroup.mk_surjective
    (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ)
  obtain ⟨uτ, huτ⟩ := QuotientGroup.mk_surjective
    (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ)
  have hspecσ := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
  have hspecτ := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem τ
  rw [← huσ, unitActionQuotientLift_mk] at hspecσ
  rw [← huτ, unitActionQuotientLift_mk] at hspecτ
  -- `σ(x)` を表す `adjoinIntegersRestrictSelfAlgHom σ` の値は、`u_σ·x` と(座標系
  -- の中の元として)一致する——両者の `K.closure` への座標がともに `σ(x)` だから。
  have hw : adjoinIntegersRestrictSelfAlgHom K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
        (⟨⟨x, hmem⟩, mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
          K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem⟩ : adjoinIntegers K x) =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (uσ : 𝒪[K.carrier]) := by
    show adjoinIntegersRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ _ = _
    apply Subtype.ext; apply Subtype.ext
    rw [coe_adjoinIntegersRestrictSelf, hspecσ]
  -- `u_τ` を `σ(x)` で評価したものは、`hw` による橋渡し(`lubinTateEvalAtPoint_congr`)
  -- を経由して `lubinTateAction_mul` により `(u_τ*u_σ)·x` に一致する。
  have hval : lubinTateActionAtAlgEquivPoint K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ uτ =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
        ((uτ : 𝒪[K.carrier]) * (uσ : 𝒪[K.carrier])) := by
    show lubinTateEvalAtPoint K x _
      (hasEval_adjoinIntegersRestrictSelfAlgHom_mk K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ)
      (LubinTateAction hq hπmax f hf0 hf1 hf uτ) = _
    rw [lubinTateEvalAtPoint_congr K x
      (hasEval_adjoinIntegersRestrictSelfAlgHom_mk K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ)
      hw (LubinTateAction hq hπmax f hf0 hf1 hf uτ)]
    exact lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem uτ uσ
  have hcomm := algEquiv_lubinTateActionAtTorsionPoint_comm K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
    σ uτ
  rw [hspecτ, hval] at hcomm
  -- `mk(u_τ*u_σ)` は `σ*τ` の定義性質を満たす。
  have hspecmul : (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      (QuotientGroup.mk (s := principalUnits K π n) (uτ * uσ))) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = (σ * τ) x := by
    rw [unitActionQuotientLift_mk]
    show (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
      ((uτ * uσ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) = _
    push_cast
    exact hcomm.symm
  -- 一意性から `reciprocityMap (σ*τ) = mk(u_τ*u_σ) = mk(u_σ)*mk(u_τ)`(`𝒪_K` の可換性)。
  have huniq := (existsUnique_unitActionQuotient_eq_algEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
    (σ * τ)).unique
    (reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (σ * τ)) hspecmul
  rw [huniq, ← huσ, ← huτ, ← QuotientGroup.mk_mul, mul_comm uτ uσ]

end ABC3.Found.PGC
