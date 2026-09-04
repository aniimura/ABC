import ABC3.Found.PGC.LubinTateReciprocityCompactness
import ABC3.Found.PGC.LubinTateReciprocityAlgEquivCompat
import ABC3.Found.PGC.LubinTateReciprocityMapLimit

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
`reciprocityMapLimit` は全射(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の最終組み立て、部品(iii)の
前半——`reciprocityMapLimit:Gal(K.closure/K.carrier)→*Compatible
Units K hπmax`の**全射性**——を、`Found/PGC/LubinTateReciprocity
Compactness.lean`のコンパクト性論法で完成させる。

## 証明の骨格

`v∈CompatibleUnits`を任意に取る。
1. `unitReductionHom_surjective`(既出)で`v`を単一の大域単数`u`
   へ持ち上げる。
2. 各`m:ℕ`について`C m:={σ:σ((psiGenSeq m).pt)=u·(psiGenSeq
   m).pt}`(`exists_algEquiv_apply_eq`——`reciprocityMap_surjective`
   +`reciprocityMap_spec`+`unitActionQuotientLift_mk`の組み合わせ)
   を定義する。
3. `reciprocityMap_eq_mk_of_apply_eq`(`reciprocityMap_spec`の逆
   方向、`unitActionQuotientLift`の単射性から)で「`σ(x)=u·x`」と
   「`reciprocityMap(x)(σ)=mk u`」の言い換えを確立する。
4. `algEquiv_apply_eq_of_succ_apply_eq`(上の言い換え+既出の
   `algEquiv_eq_of_reciprocityMap_eq_of_map_eq`+`QuotientGroup.
   map_mk`)で`C(m+1)⊆C m`(入れ子)を示す。
5. `isClosed_algEquiv_eq`(既出)で各`C m`が閉、`compactSpace_
   algEquiv`(既出)+閉部分集合のコンパクト性で`C 0`がコンパクト、
   `IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_
   isClosed`(mathlib)で交わり`⋂ C m`が空でないことを結論する。
6. その交わりの元`σ`が**すべての`m`で同時に**`σ((psiGenSeq
   m).pt)=u·(psiGenSeq m).pt`を満たすので、
   `reciprocityMap_eq_mk_of_apply_eq`+`principalUnitsQuotient
   Equiv_apply_mk`(既出、(i))で`reciprocityMapLimitFamily σ(m+1)
   =unitReductionQuotientMap(m+1)(u)=v(m+1)`(`u`の構成から)。
   `n=0`成分は`(𝒪_K/π^0)^×`が自明群なので自動的に一致する。

新しい数学的な難所は無かった——既に確立済みの部品を「配管」で
組み合わせるだけで届いた。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **`reciprocityMap`の値の特徴づけ**(`reciprocityMap_spec`の逆
方向): `σ(x)=u·x`ならば`reciprocityMap(x)(σ)=QuotientGroup.mk u`
——`unitActionQuotientLift`の単射性(`unitActionQuotientLift_
injective`、既出)と`unitActionQuotientLift_mk`(既出)から。 -/
theorem reciprocityMap_eq_mk_of_apply_eq
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
    (u : (𝒪[K.carrier])ˣ) (σ : K.closure ≃ₐ[K.carrier] K.closure)
    (hσ : σ x = (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) :
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ = QuotientGroup.mk u := by
  apply unitActionQuotientLift_injective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  have hspec := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
  have hmk := unitActionQuotientLift_mk K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem u
  apply Subtype.ext; apply Subtype.ext
  rw [hmk, ← hσ, ← hspec]

/-- **`C m`の非空性**: `u`に対応する`σ`が`reciprocityMap_surjective`
から存在する。 -/
theorem exists_algEquiv_apply_eq
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
    (u : (𝒪[K.carrier])ˣ) :
    ∃ σ : K.closure ≃ₐ[K.carrier] K.closure, σ x =
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) := by
  obtain ⟨σ, hσ⟩ := reciprocityMap_surjective K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem (QuotientGroup.mk u)
  refine ⟨σ, ?_⟩
  have hspec := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
  rw [hσ, unitActionQuotientLift_mk] at hspec
  exact hspec.symm

/-- **`C(m+1)⊆C m`**(入れ子性): `σ`が`(m+1)`段目で`u`の作用を実現
するなら、`m`段目で`u`の作用を実現する`σ'`(与えられたもの)と
`(psiGenSeq m).pt`の上で一致する——`reciprocityMap_eq_mk_of_apply_eq`
で両者の`reciprocityMap`の値を`QuotientGroup.mk u`に翻訳してから
既出の`algEquiv_eq_of_reciprocityMap_eq_of_map_eq`を適用するだけ
(`QuotientGroup.map`が`mk u`を`mk u`自身へ送ることは`QuotientGroup.
map_mk`+`MonoidHom.id`から`rfl`)。 -/
theorem algEquiv_apply_eq_of_succ_apply_eq
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ) (u : (𝒪[K.carrier])ˣ) (σ σ' : K.closure ≃ₐ[K.carrier] K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
          Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure))]
    (hσ : σ ((psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt) =
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1)
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hn
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hmem
          (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier
          ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
            Set K.closure)) : K.closure))
    (hσ' : σ' ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier
          ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) : K.closure)) :
    σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) = σ' ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) := by
  apply algEquiv_eq_of_reciprocityMap_eq_of_map_eq K hq hπmax hπne0 f hf0 hf1 hf m σ σ'
  rw [reciprocityMap_eq_mk_of_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) (by omega) _
      (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hψ _ _ u σ hσ,
    reciprocityMap_eq_mk_of_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) _
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ _ _ u σ' hσ',
    QuotientGroup.map_mk]
  rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMapLimit`は全射**。節目(5)の部品(iii)の前半、完成。 -/
theorem reciprocityMapLimit_surjective
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff)) :
    Function.Surjective (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf) := by
  intro v
  obtain ⟨u, hu⟩ := unitReductionHom_surjective K hπmax v
  haveI := isGalois_carrier_closure K
  haveI := compactSpace_algEquiv K
  haveI : ∀ m, FiniteDimensional K.carrier
      (IntermediateField.adjoin K.carrier ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) :=
    fun m => (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  haveI : ∀ m, FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
      ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
        Set K.closure)) :=
    fun m => (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hfd
  set C : ℕ → Set (K.closure ≃ₐ[K.carrier] K.closure) := fun m =>
    {σ | σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier
          ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) : K.closure)} with hC_def
  have htn : ∀ m, (C m).Nonempty := fun m =>
    exists_algEquiv_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem u
  have htcl : ∀ m, IsClosed (C m) := fun m => by
    obtain ⟨σ₀, hσ₀⟩ := htn m
    exact isClosed_algEquiv_eq K _ _ σ₀ hσ₀
  have htd : ∀ m, C (m + 1) ⊆ C m := by
    intro m σ hσ
    obtain ⟨σ', hσ'⟩ := htn m
    have hσ2 : σ ((psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt) =
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1)
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hn
          (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hmem
          (u : 𝒪[K.carrier])) :
        IntermediateField.adjoin K.carrier
          ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
            Set K.closure)) : K.closure) := hσ
    show σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) = _
    rw [algEquiv_apply_eq_of_succ_apply_eq K hq hπmax hπne0 f hf0 hf1 hf m u σ σ' hσ2 hσ']
    exact hσ'
  have h0 : IsCompact (C 0) := (htcl 0).isCompact
  obtain ⟨σ, hσ⟩ := IsCompact.nonempty_iInter_of_sequence_nonempty_isCompact_isClosed C htd htn h0 htcl
  refine ⟨σ, ?_⟩
  apply Subtype.ext
  funext n
  match n with
  | 0 =>
    have h0sub : Subsingleton (𝒪[K.carrier] ⧸ Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])))ˣ := by
      have h00 : Ideal.span ({π ^ 0} : Set (𝒪[K.carrier])) = ⊤ := by
        rw [pow_zero]; exact Ideal.span_singleton_one
      rw [h00]
      constructor
      intro a b
      apply Units.ext
      exact Subsingleton.elim _ _
    exact h0sub.elim _ _
  | (m + 1) =>
    have hCm : σ ∈ C m := Set.mem_iInter.mp hσ m
    show reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1) = _
    show principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) = _
    rw [reciprocityMap_eq_mk_of_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) _
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ _ _ u σ hCm,
      principalUnitsQuotientEquiv_apply_mk]
    exact congrFun (congrArg Subtype.val hu) (m + 1)

end ABC3.Found.PGC
