import ABC3.Found.PGC.RamifiedUnramifiedDisjoint

/-!
# 合成体への同時全射 `Γ_K ↠ Gal(E₁/K) × Gal(E₂/K)`

第 999 で `[E₁ ⊔ E₂ : K] = [E₁:K]·[E₂:K]`(`E₁` 不分岐・`E₂` 完全分岐)が
出たので、あとは**純粋な群論**で同時全射が従う。

## 体の塔を作らない

`AlgEquiv.restrictNormalHom` で `Gal(E₁⊔E₂/K) → Gal(E₁/K)` を作る道は、
`↥E₁ → ↥(E₁⊔E₂)` の代数構造を手で入れることになり、中間体の 2 層をまたぐ
`rfl` が kernel を止める(`tools/lean-idioms.md` #59)。

代わりに **`Γ_K` の側で閉じる**:

* `fixingSubgroup (E₁ ⊔ E₂) = fixingSubgroup E₁ ⊓ fixingSubgroup E₂`
* `[Γ_K : fixingSubgroup E] = [E : K]`
  (mathlib `IntermediateField.finrank_eq_fixingSubgroup_index`)
* `[Γ : H₁ ⊓ H₂] = [Γ:H₁]·[Γ:H₂]` なら `Γ → Γ/H₁ × Γ/H₂` は全射(群論)

これで体の塔は一度も要らない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

/-! ## 群論の部分 -/

/-- **★★★★★★★★指数が積なら商への同時写像は全射**。

`Γ → Γ/H₁ × Γ/H₂` の核は `H₁ ⊓ H₂` なので `Γ/(H₁⊓H₂) ↪ Γ/H₁ × Γ/H₂`。
指数が積なら濃度が一致するので単射から全単射が出る。 -/
theorem surjective_prod_quotient_of_index_mul {G : Type*} [Group G]
    (H₁ H₂ : Subgroup G) [H₁.Normal] [H₂.Normal]
    (hfin₁ : Finite (G ⧸ H₁)) (hfin₂ : Finite (G ⧸ H₂))
    (h : (H₁ ⊓ H₂).index = H₁.index * H₂.index) :
    Function.Surjective
      (fun g : G => ((QuotientGroup.mk g : G ⧸ H₁), (QuotientGroup.mk g : G ⧸ H₂))) := by
  haveI := hfin₁
  haveI := hfin₂
  set φ : G →* (G ⧸ H₁) × (G ⧸ H₂) :=
    (QuotientGroup.mk' H₁).prod (QuotientGroup.mk' H₂) with hφ
  have hker : φ.ker = H₁ ⊓ H₂ := by
    ext g
    simp [hφ, MonoidHom.mem_ker, Prod.ext_iff, QuotientGroup.eq_one_iff]
  have hinj : Function.Injective (QuotientGroup.kerLift φ) := QuotientGroup.kerLift_injective φ
  have hcard : Nat.card (G ⧸ φ.ker) = Nat.card ((G ⧸ H₁) × (G ⧸ H₂)) := by
    rw [hker, Nat.card_prod]
    simpa [Subgroup.index] using h
  have hbij : Function.Bijective (QuotientGroup.kerLift φ) :=
    (Nat.bijective_iff_injective_and_card _).mpr ⟨hinj, hcard⟩
  intro y
  obtain ⟨q, hq⟩ := hbij.2 y
  obtain ⟨g, rfl⟩ := QuotientGroup.mk_surjective q
  exact ⟨g, hq⟩

/-! ## 体論の部分 -/

/-- 固定部分群は sup を inf に移す——`fixingSubgroup`/`fixedField` の
Galois 接続(`IntermediateField.le_iff_le`)から。 -/
theorem fixingSubgroup_sup {F E : Type*} [Field F] [Field E] [Algebra F E]
    (A B : IntermediateField F E) :
    (A ⊔ B).fixingSubgroup = A.fixingSubgroup ⊓ B.fixingSubgroup := by
  have hanti : ∀ (C D : IntermediateField F E), C ≤ D → D.fixingSubgroup ≤ C.fixingSubgroup := by
    intro C D hCD σ hσ
    rw [IntermediateField.mem_fixingSubgroup_iff] at hσ ⊢
    exact fun z hz => hσ z (hCD hz)
  refine le_antisymm (le_inf (hanti _ _ le_sup_left) (hanti _ _ le_sup_right)) ?_
  rw [← IntermediateField.le_iff_le]
  exact sup_le
    ((IntermediateField.le_iff_le (A.fixingSubgroup ⊓ B.fixingSubgroup) A).mpr inf_le_left)
    ((IntermediateField.le_iff_le (A.fixingSubgroup ⊓ B.fixingSubgroup) B).mpr inf_le_right)

variable {p : ℕ} [Fact p.Prime]

/-- Galois な中間体の固定部分群は正規。 -/
theorem normal_fixingSubgroup_of_isGalois (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier E]
    [IsGalois K.carrier E] : (E.fixingSubgroup).Normal := by
  haveI := isGalois_closure K
  exact (InfiniteGalois.normal_iff_isGalois E).mpr inferInstance

/-- **★★★★★★★★★★★★★★★★`Γ_K ↠ Γ_K/H₁ × Γ_K/H₂`**——
`E₁`・`E₂` が Galois で、合成体の次数が積なら。 -/
theorem surjective_absGal_prod_quotient (K : PAdicLocalField p)
    (E₁ E₂ : IntermediateField K.carrier K.closure)
    [FiniteDimensional K.carrier E₁] [FiniteDimensional K.carrier E₂]
    [IsGalois K.carrier E₁] [IsGalois K.carrier E₂]
    (hmul : Module.finrank K.carrier (E₁ ⊔ E₂ : IntermediateField K.carrier K.closure)
      = Module.finrank K.carrier E₁ * Module.finrank K.carrier E₂) :
    Function.Surjective (fun g : K.absGal =>
      ((QuotientGroup.mk g : K.absGal ⧸ E₁.fixingSubgroup),
       (QuotientGroup.mk g : K.absGal ⧸ E₂.fixingSubgroup))) := by
  haveI := normal_fixingSubgroup_of_isGalois K E₁
  haveI := normal_fixingSubgroup_of_isGalois K E₂
  haveI : FiniteDimensional K.carrier (E₁ ⊔ E₂ : IntermediateField K.carrier K.closure) :=
    IntermediateField.finiteDimensional_sup E₁ E₂
  have hi₁ : (E₁.fixingSubgroup).index = Module.finrank K.carrier E₁ :=
    (IntermediateField.finrank_eq_fixingSubgroup_index E₁).symm
  have hi₂ : (E₂.fixingSubgroup).index = Module.finrank K.carrier E₂ :=
    (IntermediateField.finrank_eq_fixingSubgroup_index E₂).symm
  have hinf : (E₁.fixingSubgroup ⊓ E₂.fixingSubgroup).index
      = Module.finrank K.carrier (E₁ ⊔ E₂ : IntermediateField K.carrier K.closure) := by
    rw [← fixingSubgroup_sup]
    exact (IntermediateField.finrank_eq_fixingSubgroup_index _).symm
  have hfin₁ : Finite (K.absGal ⧸ E₁.fixingSubgroup) := by
    haveI : (E₁.fixingSubgroup).FiniteIndex := ⟨by rw [hi₁]; exact Module.finrank_pos.ne'⟩
    exact Subgroup.finite_quotient_of_finiteIndex
  have hfin₂ : Finite (K.absGal ⧸ E₂.fixingSubgroup) := by
    haveI : (E₂.fixingSubgroup).FiniteIndex := ⟨by rw [hi₂]; exact Module.finrank_pos.ne'⟩
    exact Subgroup.finite_quotient_of_finiteIndex
  exact surjective_prod_quotient_of_index_mul _ _ hfin₁ hfin₂ (by rw [hinf, hi₁, hi₂, hmul])

/-- **★★★★★★★★★★★★★★★★★★★★`Γ_K ↠ Γ_K/H_{K(y)} × Γ_K/H_{K(Λ_n)}`**

`y` は次数 `m` の不分岐拡大の生成元、`x` は `ψ_n` の捩れ点。
第 999 の次数の積と、上の群論を繋いだもの。

★`Γ_K/H_{K(x)} ≅ Gal(K(x)/K) ≃* (𝒪_K/π^n)^×`(`galoisReciprocityEquiv`)、
`Γ_K/H_{K(y)} ≅ Gal(K(y)/K) ≃* ℤ/m`(`exists_gal_mulEquiv_zmod`)なので、
これは相互律 `Γ_K ↠ 𝒪_K^× × Ẑ` の**有限段での形**である。 -/
theorem surjective_absGal_prod_lubinTate_unramified (K : PAdicLocalField p)
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (m : ℕ) (hm : m ≠ 0) :
    ∃ y : K.closure,
      Module.finrank K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) = m
      ∧ Function.Surjective (fun g : K.absGal =>
        ((QuotientGroup.mk g :
            K.absGal ⧸ (IntermediateField.adjoin K.carrier ({y} : Set K.closure)).fixingSubgroup),
         (QuotientGroup.mk g :
            K.absGal ⧸ (IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup)))
      := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  obtain ⟨y, hym, hyu, hygal⟩ := exists_isGalois_isUnramifiedAdjoin K m hm
  haveI := hygal
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x
      hxψ hxn hmem
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  have hx := isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x
    hxψ hxn hmem
  exact ⟨y, hym,
    surjective_absGal_prod_quotient K _ _
      (finrank_sup_of_isUnramified_of_isTotallyRamified K hyu hx)⟩

/-- **`Γ_K / Gal(K̄/E) ≃* Gal(E/K)`**——`E` が `K` 上正規なら。
`IntermediateField.restrictNormalHom_ker` と
`AlgEquiv.restrictNormalHom_surjective` を繋いだだけ
(`InertiaIdentification.lean::absGalQuotKerEquivUnramifiedGal` の一般化)。 -/
noncomputable def absGalQuotFixingSubgroupEquiv (K : PAdicLocalField p)
    (E : IntermediateField K.carrier K.closure) [Normal K.carrier E]
    [E.fixingSubgroup.Normal] :
    (K.absGal ⧸ E.fixingSubgroup) ≃* (E ≃ₐ[K.carrier] E) := by
  haveI := IsAlgClosure.normal K.carrier K.closure
  exact (QuotientGroup.quotientMulEquivOfEq
      (IntermediateField.restrictNormalHom_ker E).symm).trans
    (QuotientGroup.quotientKerEquivOfSurjective _
      (AlgEquiv.restrictNormalHom_surjective (F := K.carrier) (K₁ := E) (E := K.closure)))

/-- **★★★★★★★★★★★★★★★★★★★★★★★★★★`Γ_K ↠ ℤ/m × (𝒪_K/π^n)^×`**

**相互律 `Γ_K ↠ 𝒪_K^× × Ẑ` の、有限段での完全な形**。

* 第一成分: `Γ_K ↠ Gal(K(y)/K) ≃* ℤ/m`(`y` は次数 `m` の不分岐拡大の生成元、
  `exists_gal_mulEquiv_zmod`)
* 第二成分: `Γ_K ↠ Gal(K(Λ_n)/K) ≃* (𝒪_K/π^n)^×`(`galoisReciprocityEquiv`)
* 同時全射性: 第 999(次数が積)+ 第 1001(群論) -/
theorem exists_surjective_absGal_zmod_prod_units (K : PAdicLocalField p)
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (m : ℕ) (hm : m ≠ 0) :
    ∃ F : K.absGal →* (Multiplicative (ZMod m)
        × (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set 𝒪[K.carrier]))ˣ),
      Function.Surjective F := by
  classical
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  obtain ⟨y, hym, hyu, ⟨ezm⟩⟩ := exists_gal_mulEquiv_zmod K m hm
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    normal_of_isUnramifiedAdjoin K y hyu
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({y} : Set K.closure)) := ⟨⟩
  haveI : Normal K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn x
      hxψ hxn hmem
  haveI : Algebra.IsSeparable K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IntermediateField.isSeparable_tower_bot K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  haveI : IsGalois K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := ⟨⟩
  have hx := isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x
    hxψ hxn hmem
  have hsurj := surjective_absGal_prod_quotient K
    (IntermediateField.adjoin K.carrier ({y} : Set K.closure))
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (finrank_sup_of_isUnramified_of_isTotallyRamified K hyu hx)
  haveI := normal_fixingSubgroup_of_isGalois K
    (IntermediateField.adjoin K.carrier ({y} : Set K.closure))
  haveI := normal_fixingSubgroup_of_isGalois K
    (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
  set e₁ : (K.absGal ⧸ (IntermediateField.adjoin K.carrier
      ({y} : Set K.closure)).fixingSubgroup) ≃* Multiplicative (ZMod m) :=
    (absGalQuotFixingSubgroupEquiv K
      (IntermediateField.adjoin K.carrier ({y} : Set K.closure))).trans ezm with he₁
  set e₂ : (K.absGal ⧸ (IntermediateField.adjoin K.carrier
      ({x} : Set K.closure)).fixingSubgroup)
      ≃* (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set 𝒪[K.carrier]))ˣ :=
    (absGalQuotFixingSubgroupEquiv K
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))).trans
      (galoisReciprocityEquiv K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem) with he₂
  refine ⟨(e₁.toMonoidHom.comp (QuotientGroup.mk' _)).prod
    (e₂.toMonoidHom.comp (QuotientGroup.mk' _)), fun z => ?_⟩
  obtain ⟨g, hg⟩ := hsurj (e₁.symm z.1, e₂.symm z.2)
  refine ⟨g, ?_⟩
  have h1 : (QuotientGroup.mk g :
      K.absGal ⧸ (IntermediateField.adjoin K.carrier
        ({y} : Set K.closure)).fixingSubgroup) = e₁.symm z.1 := congrArg Prod.fst hg
  have h2 : (QuotientGroup.mk g :
      K.absGal ⧸ (IntermediateField.adjoin K.carrier
        ({x} : Set K.closure)).fixingSubgroup) = e₂.symm z.2 := congrArg Prod.snd hg
  show (e₁ (QuotientGroup.mk g), e₂ (QuotientGroup.mk g)) = z
  rw [h1, h2, MulEquiv.apply_symm_apply, MulEquiv.apply_symm_apply]

/-- **★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Γ_K^{ab} ↠ ℤ/m × (𝒪_K/π^n)^×`**

Prop 1.2 が数えるのは `Γ_K` そのものではなく**アーベル化** `Γ_K^{ab}` の不変量
(捩れの prime-to-p 部分が `q−1` 個、pro-p 部分の階数が `[K:ℚ_p]+1`)なので、
第 1002 をアーベル化に持ち上げておく。行き先が可換なので
`Abelianization.lift` を当てるだけ。 -/
theorem exists_surjective_abelianization_zmod_prod_units (K : PAdicLocalField p)
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
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (m : ℕ) (hm : m ≠ 0) :
    ∃ F : Abelianization K.absGal →* (Multiplicative (ZMod m)
        × (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set 𝒪[K.carrier]))ˣ),
      Function.Surjective F := by
  obtain ⟨F, hF⟩ := exists_surjective_absGal_zmod_prod_units K hq hπmax hπne0 f hf0 hf1 hf
    n hn x hxψ hxn hmem m hm
  refine ⟨Abelianization.lift (G := K.absGal)
    (A := Multiplicative (ZMod m)
      × (𝒪[K.carrier] ⧸ Ideal.span ({π ^ n} : Set 𝒪[K.carrier]))ˣ) F, fun z => ?_⟩
  obtain ⟨g, hg⟩ := hF z
  refine ⟨Abelianization.of g, ?_⟩
  rw [Abelianization.lift_apply_of]
  exact hg

end ABC3.Found.PGC
