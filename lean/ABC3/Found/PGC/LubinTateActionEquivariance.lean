import ABC3.Found.PGC.LubinTateActionBijective
import ABC3.Found.PGC.PowerSeriesAevalComm

/-!
# Galois同変性への足場: `σ` は `K.carrier⟮x⟯` を自己同型として保つ(`sorry` 無し)

`Gal(K.closure/K.carrier)`の各`σ`に対する`σ(a·x)=a·σ(x)`
(Galois同変性そのもの)を目指す最初の段階として、これまで
「異なる`adjoinIntegers K x`・`adjoinIntegers K (σx)`の間を
橋渡しする必要がある(cross-point instance bridging)」と警戒
されてきた壁を、**`σ(x)`自身が`K.carrier⟮x⟯`の中に留まる**という
事実(単数の作用の全単射性、`Found/PGC/LubinTateActionBijective.lean`
の直接の系)を使って回避する。

## 鍵となる観察

`σ(x)`は`ψ_nの根`に留まる(`algEquiv_mem_iteratedLubinTatePsi
TorsionPoints_of_mem`、既出)。ところが`ψ_nの根`は`{u·x:u∈(𝒪_K)^×}`
そのもの(`unitActionQuotientBijOn`の**全射性**、既出)であり、
`u·x`は定義から常に`adjoinIntegers K x⊆K.carrier⟮x⟯`に入る——
よって`σ(x)∈K.carrier⟮x⟯`が、新しい構成を一切経由せず既存の結果
だけから出る。

これを`IntermediateField.adjoin_map`(`σ`による像の adjoin は
`σ(x)`の adjoin)と組み合わせると、`σ`は`K.carrier⟮x⟯`を
**自分自身の中へ**写す(`map σ.toAlgHom L≤L`)。単射性(体の
埋め込みは常に単射)+有限次元性から、次元の等式
(`IntermediateField.eq_of_le_of_finrank_eq`)で**全射性**まで
持ち上がり、`σ`は`K.carrier⟮x⟯`の**自己同型**に制限される
(`algEquivRestrictSelf`)——`adjoinIntegers K x`と`adjoinIntegers
K (σx)`という異なる2つの環を橋渡しする必要が、そもそも生じない。

さらに`norm_algEquiv_eq`(既出、`σ`はspectralNormの等長写像)から
この制限もノルムを保つ(`norm_algEquivRestrictSelf`)ので、
`adjoinIntegers K x`(ノルム`≤1`の部分環)へもそのまま制限できる
(`adjoinIntegersRestrictSelf`)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- `ψ_n の根` は `K.carrier⟮x⟯` に含まれる——単数の作用の全射性
(`unitActionQuotientBijOn`、既出)から、任意の`ψ_nの根`の元は
`u·x`(`u∈(𝒪_K)^×`)の形に書け、これは常に`adjoinIntegers K x⊆
K.carrier⟮x⟯`に入る。 -/
theorem iteratedLubinTatePsiTorsionPoints_subset_adjoin
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
    (y : K.closure) (hy : y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    y ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  obtain ⟨U, hU⟩ := (unitActionQuotientBijOn_bijective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem).2
    ⟨y, hy⟩
  have : y = (↑(↑(unitActionQuotientLift K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem U) :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) :=
    (congrArg Subtype.val hU).symm
  rw [this]
  exact SetLike.coe_mem _

/-- **`σ`は`K.carrier⟮x⟯`を自分自身の中へ写す(実は上へも)**:
`IntermediateField.map σ.toAlgHom K.carrier⟮x⟯ = K.carrier⟮x⟯`。
`σ(x)∈K.carrier⟮x⟯`(上の補題)から`≤`が出て、単射性+有限次元の
次元の等式から`=`まで持ち上がる。 -/
theorem algEquiv_map_adjoin_eq
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
    IntermediateField.map σ.toAlgHom (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) =
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
  have hσxψ : σ x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn :=
    algEquiv_mem_iteratedLubinTatePsiTorsionPoints_of_mem K hq hπmax hπne0 f hf0 hf1 hf n hn σ x hxψ
  have hσxL : σ x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
    iteratedLubinTatePsiTorsionPoints_subset_adjoin K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
      (σ x) hσxψ
  have hle : IntermediateField.map σ.toAlgHom (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) ≤
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) := by
    rw [IntermediateField.adjoin_map, Set.image_singleton]
    exact IntermediateField.adjoin_simple_le_iff.mpr hσxL
  have hdim : Module.finrank K.carrier
      (IntermediateField.map σ.toAlgHom (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) =
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    (LinearEquiv.finrank_eq
      (IntermediateField.equivMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        σ.toAlgHom).toLinearEquiv).symm
  exact IntermediateField.eq_of_le_of_finrank_eq hle hdim

/-- **`σ`を`K.carrier⟮x⟯`自身への自己同型に制限したもの**——
`algEquiv_map_adjoin_eq`(像が`K.carrier⟮x⟯`自身に一致)を経由して
`IntermediateField.equivMap`(単射性からの像への同型)を
`IntermediateField.equivOfEq`(像=`K.carrier⟮x⟯`という等式からの
同一視)と合成するだけ。`adjoinIntegers K x`と`adjoinIntegers
K (σx)`という異なる2つの環を橋渡しする必要が、これで無くなる。 -/
noncomputable def algEquivRestrictSelf
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
    IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure) :=
  (IntermediateField.equivMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) σ.toAlgHom).trans
    (IntermediateField.equivOfEq
      (algEquiv_map_adjoin_eq K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ))

/-- `algEquivRestrictSelf`の座標——`K.closure`へ戻せば`σ`そのもの
(`rfl`)。 -/
theorem coe_algEquivRestrictSelf
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
    (σ : K.closure ≃ₐ[K.carrier] K.closure)
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    (↑(algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ z) : K.closure) =
      σ (z : K.closure) := by
  rfl

/-- `algEquivRestrictSelf`はノルムを保つ——`norm_algEquiv_eq`
(既出)からの直接の帰結。 -/
theorem norm_algEquivRestrictSelf
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
    (σ : K.closure ≃ₐ[K.carrier] K.closure)
    (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :
    ‖algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ z‖ = ‖z‖ := by
  show spectralNorm K.carrier K.closure
      (↑(algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ z) : K.closure) =
    spectralNorm K.carrier K.closure (z : K.closure)
  rw [coe_algEquivRestrictSelf]
  exact norm_algEquiv_eq K σ (z : K.closure)

/-- **`σ`を`adjoinIntegers K x`自身への自己写像に制限したもの**——
ノルムを保つ(上の補題)ので、ノルム`≤1`の部分環をそのまま保つ。 -/
noncomputable def adjoinIntegersRestrictSelf
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
    (σ : K.closure ≃ₐ[K.carrier] K.closure) (z : adjoinIntegers K x) : adjoinIntegers K x :=
  ⟨algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
      (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)),
    by
      show ‖algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
        (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ ≤ 1
      rw [norm_algEquivRestrictSelf]
      exact z.2⟩

/-- `adjoinIntegersRestrictSelf`の座標——`K.closure`へ戻せば`σ`
そのもの。 -/
theorem coe_adjoinIntegersRestrictSelf
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
    (σ : K.closure ≃ₐ[K.carrier] K.closure) (z : adjoinIntegers K x) :
    (↑(↑(adjoinIntegersRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ z) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) =
      σ (↑(↑z : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) := by
  show (↑(algEquivRestrictSelf K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
      (z : IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure) = _
  rw [coe_algEquivRestrictSelf]

end ABC3.Found.PGC
