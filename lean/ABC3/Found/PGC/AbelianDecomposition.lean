import ABC3.Found.PGC.LubinTateClosure
import ABC3.Found.PGC.LubinTateTotallyRamified

/-!
# `K_π ⊓ K^ur = K` と `Gal(K_π · K^ur / K) ≅ 𝒪_K^× × Gal(K^ur/K)`(`sorry` 無し)

経路 Λ の節点 Λ4。Λ3(`Found/PGC/LubinTateClosure.lean`)で作った
`K_π`(全捩れ点で生成した中間体、`Gal(K_π/K) ≅ 𝒪_K^×`)と、
不分岐閉包 `K^ur` が **交わらない**ことを示し、合成体の Galois 群が
直積に分解することを出す。

## 在庫にあったもの / 無かったもの

**在庫にあった**(有限段):

* `TotallyRamified.lean::totallyRamifiedAdjoin_inf_unramifiedClosure`
  ——`K(α) ⊓ K^ur = ⊥`(`α` は完全分岐な単項生成元)
* `TotallyRamified.lean::finrank_eq_one_of_mem_unramifiedClosure_of_le`
  ——`K(x) ≤ K(y)`・`x ∈ K^ur`・`K(y)/K` 完全分岐 ⟹ `[K(x):K] = 1`
* `LubinTateTotallyRamified.lean::isTotallyRamifiedAdjoin_iteratedLubinTatePsi`
  ——`ψ_n` の根 `x` について `K(x)/K` は完全分岐
* `RamifiedUnramifiedDisjoint.lean`——線型無関係と次数の積(有限段)

**在庫に無かった**(本ファイルで足した):

* **無限次での** `K_π ⊓ K^ur = ⊥`
  (`lubinTateClosure_inf_unramifiedClosure`)。
  `K_π` の元はどれか 1 つの `K(x_{m+1})` に入る
  (Λ3 の `exists_mem_adjoin_psiGenSeq_of_mem_lubinTateClosure`、
  塔が増大列であることから)ので、有限段の結果に帰着する。
* **無限次 Galois での直積分解**
  (`galSupEquivProd`: `A ⊓ B = ⊥` ⟹ `Gal(A·B/k) ≅ Gal(A/k) × Gal(B/k)`)。
  副有限位相を使う: `Gal(k̄/A)`・`Gal(k̄/B)` は閉(mathlib の
  `InfiniteGalois.fixingSubgroup_isClosed`)、`Gal(k̄/k)` はコンパクト
  かつ T2 なので、片方が正規なら 2 つの積 `N_A · N_B` はコンパクト、
  従って閉。無限次 Galois 対応(`InfiniteGalois.fixingSubgroup_
  fixedField`)により `N_A ⊔ N_B = ⊤`、これが全射性を与える。

## ★やっていないこと(Λ5 の担当)

`Gal(K^ur/K) ≅ Ẑ`(`ProfiniteGrp.ProfiniteCompletion.completion`)は
**本ファイルでは主張していない**——これは Frobenius 元の同定
(経路 Λ の Λ5)そのものなので、そちらに回す。したがって本ファイルの
結論は `Gal(K_π · K^ur / K) ≅ 𝒪_K^× × Gal(K^ur/K)` の形で止めてある。
`Gal(K^ur/K) ≅ Ẑ` が入れば、そのまま `𝒪_K^× × Ẑ` になる。

## ★設計上の注意(守ったこと)

* 代数的な `Abelianization` は**使っていない**(副有限群では有限指数
  部分群が開とは限らないため)。本ファイルは商群を一切作らず、
  `QuotientGroup.quotientKerEquivOfSurjective` を核の一致に対して
  使うだけ。
* 惰性群 `inertia` を経由していない——不分岐側は
  `unramifiedClosure K` という**体**として直接扱っている。
  `Found/PGC/InertiaIdentification.lean` も import していない
  (`Normal K.carrier K^ur` は `UnramifiedExtension.lean::
  normal_unramifiedClosure` から直接取る)。
  ★ただし `Found/PGC/TotallyRamified.lean`(有限段の非交叉を供給する)
  の import 閉包が既に `UnramifiedCriterion` 経由で
  `SubgroupCorrespondenceConstruction` と `Skeleton/PGC/Section1Cor13.lean`
  に届いている。本ファイルの主張・証明はそれらを一切参照しておらず
  (`#print axioms` に `sorryAx` は出ない)、
  `inertia`/`residueCardinality`/`subgroupCorrespondence` という語も
  出てこないが、ファイル依存としては既存構造の上に乗っている。
* 部分群や `n` を自由なパラメータとして結論に出していない——
  仮説なしの形(`exists_lubinTateUnramified_decomposition`)では
  すべて `∃` の内側に閉じ込めてある。

## 逸脱(記録)

Λ3 と同じく、`K_π` は Lubin-Tate 級数 `f` にも依存する形になっている
(`f` 非依存性は Λ6 の担当)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 無限次 Galois 拡大での直積分解(一般論) -/

section General

variable {k E : Type*} [Field k] [Field E] [Algebra k E]

/-- `fixedField` は部分群の `⊔` を中間体の `⊓` に移す。 -/
theorem fixedField_sup (H H' : Subgroup (E ≃ₐ[k] E)) :
    IntermediateField.fixedField (H ⊔ H') =
      IntermediateField.fixedField H ⊓ IntermediateField.fixedField H' := by
  ext x
  simp only [IntermediateField.mem_inf, IntermediateField.mem_fixedField_iff]
  constructor
  · intro h
    exact ⟨fun g hg => h g (le_sup_left (a := H) (b := H') hg),
      fun g hg => h g (le_sup_right (a := H) (b := H') hg)⟩
  · rintro ⟨h1, h2⟩ g hg
    have hle : H ⊔ H' ≤ MulAction.stabilizer (E ≃ₐ[k] E) x := by
      rw [sup_le_iff]
      exact ⟨fun a ha => h1 a ha, fun a ha => h2 a ha⟩
    exact hle hg

/-- 正規な中間体の固定部分群は正規部分群——制限写像の核だから。 -/
theorem normal_fixingSubgroup (B : IntermediateField k E) [Normal k B] :
    B.fixingSubgroup.Normal := by
  rw [← ker_restrictNormalHom_eq_fixingSubgroup B]
  infer_instance

/-- 2 つの固定部分群が生成する部分群は**閉**。

片方が正規なら `↑(N_A ⊔ N_B) = ↑N_A * ↑N_B`(`Subgroup.mul_normal`)で、
`N_A`・`N_B` は閉(`InfiniteGalois.fixingSubgroup_isClosed`)、
`Gal(E/k)` はコンパクトなので両者はコンパクト、積もコンパクト、
`Gal(E/k)` は T2 なので閉。 -/
theorem isClosed_sup_fixingSubgroup [IsGalois k E]
    (A B : IntermediateField k E) [Normal k B] :
    IsClosed ((A.fixingSubgroup ⊔ B.fixingSubgroup : Subgroup (E ≃ₐ[k] E)) :
      Set (E ≃ₐ[k] E)) := by
  haveI := normal_fixingSubgroup B
  rw [Subgroup.mul_normal]
  exact (((InfiniteGalois.fixingSubgroup_isClosed A).isCompact).mul
    ((InfiniteGalois.fixingSubgroup_isClosed B).isCompact)).isClosed

/-- ★**`A ⊓ B = k` なら `Gal(E/A)` と `Gal(E/B)` は `Gal(E/k)` 全体を
生成する**——無限次 Galois 対応(`InfiniteGalois.fixingSubgroup_
fixedField`)を、上で示した閉性に対して適用する。 -/
theorem sup_fixingSubgroup_eq_top [IsGalois k E]
    (A B : IntermediateField k E) [Normal k B] (h : A ⊓ B = ⊥) :
    A.fixingSubgroup ⊔ B.fixingSubgroup = ⊤ := by
  have hcl := InfiniteGalois.fixingSubgroup_fixedField
    (⟨A.fixingSubgroup ⊔ B.fixingSubgroup, isClosed_sup_fixingSubgroup A B⟩ :
      ClosedSubgroup (E ≃ₐ[k] E))
  rw [show ((⟨A.fixingSubgroup ⊔ B.fixingSubgroup, isClosed_sup_fixingSubgroup A B⟩ :
      ClosedSubgroup (E ≃ₐ[k] E)) : Subgroup (E ≃ₐ[k] E)) =
      A.fixingSubgroup ⊔ B.fixingSubgroup from rfl] at hcl
  rw [fixedField_sup, InfiniteGalois.fixedField_fixingSubgroup,
    InfiniteGalois.fixedField_fixingSubgroup, h, IntermediateField.fixingSubgroup_bot] at hcl
  exact hcl.symm

/-- `A ⊔ B` は `A ∪ B` で生成される中間体。 -/
theorem sup_eq_adjoin_union (A B : IntermediateField k E) :
    A ⊔ B = IntermediateField.adjoin k ((A : Set E) ∪ (B : Set E)) := by
  refine le_antisymm (sup_le ?_ ?_) ?_
  · intro x hx
    exact IntermediateField.subset_adjoin k _ (Or.inl hx)
  · intro x hx
    exact IntermediateField.subset_adjoin k _ (Or.inr hx)
  · rw [IntermediateField.adjoin_le_iff]
    rintro x (hx | hx)
    · exact le_sup_left (a := A) (b := B) hx
    · exact le_sup_right (a := A) (b := B) hx

/-- ★`Gal(E/A·B) = Gal(E/A) ∩ Gal(E/B)`。 -/
theorem fixingSubgroup_sup (A B : IntermediateField k E) :
    (A ⊔ B).fixingSubgroup = A.fixingSubgroup ⊓ B.fixingSubgroup := by
  ext σ
  rw [sup_eq_adjoin_union, mem_fixingSubgroup_adjoin_iff]
  simp only [Subgroup.mem_inf, IntermediateField.mem_fixingSubgroup_iff, Set.mem_union]
  constructor
  · intro h
    exact ⟨fun x hx => h x (Or.inl hx), fun x hx => h x (Or.inr hx)⟩
  · rintro ⟨h1, h2⟩ x (hx | hx)
    · exact h1 x hx
    · exact h2 x hx

/-- `σ ↦ (σ|_A, σ|_B)`。 -/
noncomputable def restrictPairHom (A B : IntermediateField k E) [Normal k A] [Normal k B] :
    (E ≃ₐ[k] E) →* ((A : Type _) ≃ₐ[k] (A : Type _)) × ((B : Type _) ≃ₐ[k] (B : Type _)) :=
  (AlgEquiv.restrictNormalHom (A : Type _)).prod (AlgEquiv.restrictNormalHom (B : Type _))

theorem ker_restrictPairHom (A B : IntermediateField k E) [Normal k A] [Normal k B] :
    (restrictPairHom A B).ker = (A ⊔ B).fixingSubgroup := by
  rw [restrictPairHom, MonoidHom.ker_prod, ker_restrictNormalHom_eq_fixingSubgroup,
    ker_restrictNormalHom_eq_fixingSubgroup, fixingSubgroup_sup]

/-- ★**`A ⊓ B = k` なら `σ ↦ (σ|_A, σ|_B)` は全射**。

`σ|_A = a` なる `σ` と `τ|_B = b` なる `τ` を取り、`σ⁻¹τ = n_A n_B`
(`sup_fixingSubgroup_eq_top` から)と分解して `ρ := σ n_A = τ n_B⁻¹`
と置けばよい。 -/
theorem restrictPairHom_surjective [IsGalois k E]
    (A B : IntermediateField k E) [Normal k A] [Normal k B] (h : A ⊓ B = ⊥) :
    Function.Surjective (restrictPairHom A B) := by
  rintro ⟨a, b⟩
  obtain ⟨σ, hσ⟩ := AlgEquiv.restrictNormalHom_surjective (F := k) (K₁ := (A : Type _)) E a
  obtain ⟨τ, hτ⟩ := AlgEquiv.restrictNormalHom_surjective (F := k) (K₁ := (B : Type _)) E b
  haveI := normal_fixingSubgroup B
  have hmem : σ⁻¹ * τ ∈ ((A.fixingSubgroup ⊔ B.fixingSubgroup : Subgroup (E ≃ₐ[k] E)) :
      Set (E ≃ₐ[k] E)) := by
    rw [sup_fixingSubgroup_eq_top A B h]
    trivial
  rw [Subgroup.mul_normal] at hmem
  obtain ⟨na, hna, nb, hnb, hab⟩ := hmem
  have hna1 : AlgEquiv.restrictNormalHom (A : Type _) na = 1 :=
    (restrictNormalHom_eq_one_iff' A na).mpr hna
  have hnb1 : AlgEquiv.restrictNormalHom (B : Type _) nb = 1 :=
    (restrictNormalHom_eq_one_iff' B nb).mpr hnb
  have heq : σ * na = τ * nb⁻¹ := by
    have hna' : na = σ⁻¹ * τ * nb⁻¹ := by rw [← hab]; group
    rw [hna']; group
  refine ⟨σ * na, Prod.ext ?_ ?_⟩
  · show AlgEquiv.restrictNormalHom (A : Type _) (σ * na) = a
    rw [map_mul, hna1, mul_one, hσ]
  · show AlgEquiv.restrictNormalHom (B : Type _) (σ * na) = b
    rw [heq, map_mul, map_inv, hnb1, inv_one, mul_one, hτ]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`A ⊓ B = k` なら `Gal(A·B/k) ≅ Gal(A/k) × Gal(B/k)`**(無限次でも)。

`Gal(E/k) ↠ Gal(A·B/k)` と `Gal(E/k) ↠ Gal(A/k) × Gal(B/k)` の核が
どちらも `(A ⊔ B).fixingSubgroup` に一致することから。 -/
noncomputable def galSupEquivProd [IsGalois k E]
    (A B : IntermediateField k E) [Normal k A] [Normal k B] (h : A ⊓ B = ⊥) :
    ((A ⊔ B : IntermediateField k E) ≃ₐ[k] (A ⊔ B : IntermediateField k E)) ≃*
      ((A : Type _) ≃ₐ[k] (A : Type _)) × ((B : Type _) ≃ₐ[k] (B : Type _)) := by
  have hker : MonoidHom.ker (AlgEquiv.restrictNormalHom (F := k) (K₁ := E)
        (E := ((A ⊔ B : IntermediateField k E) : Type _))) =
      MonoidHom.ker (restrictPairHom A B) := by
    rw [ker_restrictNormalHom_eq_fixingSubgroup, ker_restrictPairHom]
  exact
    (QuotientGroup.quotientKerEquivOfSurjective _
      (AlgEquiv.restrictNormalHom_surjective (F := k)
        (K₁ := ((A ⊔ B : IntermediateField k E) : Type _)) E)).symm.trans
      ((QuotientGroup.quotientMulEquivOfEq hker).trans
        (QuotientGroup.quotientKerEquivOfSurjective _ (restrictPairHom_surjective A B h)))

end General

/-! ## 1. `K_π ⊓ K^ur = K` -/

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-- 各 `K(x_{m+1})/K` は完全分岐——`x_{m+1}` は `ψ_{m+1}` の根だから。 -/
theorem isTotallyRamifiedAdjoin_psiGenSeq (m : ℕ) :
    @IsTotallyRamifiedAdjoin p _ K (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd := by
  haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  exact isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
    (by omega) _ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K_π ⊓ K^ur = K`**——完全分岐側と不分岐側は無限次でも交わらない。

`x ∈ K_π` はどれか 1 つの `K(x_{m+1})` に入る(塔が増大列)ので、
有限段の `finrank_eq_one_of_mem_unramifiedClosure_of_le` に帰着する。 -/
theorem lubinTateClosure_inf_unramifiedClosure :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊓ unramifiedClosure K = ⊥ := by
  refine le_antisymm ?_ bot_le
  intro x hx
  obtain ⟨hxπ, hxur⟩ := hx
  obtain ⟨m, hxm⟩ :=
    exists_mem_adjoin_psiGenSeq_of_mem_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf hxπ
  haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  have hle : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) :=
    IntermediateField.adjoin_simple_le_iff.mpr hxm
  have h1 := finrank_eq_one_of_mem_unramifiedClosure_of_le K hle hxur
    (isTotallyRamifiedAdjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)
  have hbot : IntermediateField.adjoin K.carrier ({x} : Set K.closure) = ⊥ :=
    IntermediateField.finrank_eq_one_iff.mp h1
  have hx0 : x ∈ (⊥ : IntermediateField K.carrier K.closure) := by
    rw [← hbot]
    exact IntermediateField.mem_adjoin_simple_self K.carrier x
  exact hx0

/-! ## 2. `Gal(K_π · K^ur / K) ≅ 𝒪_K^× × Gal(K^ur/K)` -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K_π · K^ur / K) ≅ 𝒪_K^× × Gal(K^ur/K)`**。

`galSupEquivProd`(交わりが自明な 2 つの正規中間体の合成)に
`lubinTateClosure_inf_unramifiedClosure` を渡し、第 1 成分を
Λ3 の `lubinTateClosureGalEquivUnits` で `𝒪_K^×` に取り替えるだけ。

★第 2 成分 `Gal(K^ur/K)` を `Ẑ` と同定するのは Λ5(Frobenius)の担当。 -/
noncomputable def lubinTateUnramifiedGalEquivProd :
    ((lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
        IntermediateField K.carrier K.closure) ≃ₐ[K.carrier]
      (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
        IntermediateField K.carrier K.closure)) ≃*
      ((𝒪[K.carrier])ˣ × (unramifiedClosure K ≃ₐ[K.carrier] unramifiedClosure K)) := by
  haveI := isGalois_carrier_closure K
  haveI := normal_unramifiedClosure K
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  exact (galSupEquivProd (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) (unramifiedClosure K)
      (lubinTateClosure_inf_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf)).trans
    ((lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf).prodCongr (MulEquiv.refl _))

/-! ## 3. 仮説なしの形 -/

/-- ★★★★**仮説なしの Λ4**: 任意の p 進局所体 `K` について、
`K̄` の中に `K` 上正規な中間体 `K_π` が在って

* `K_π ⊓ K^ur = K`
* `Gal(K_π/K) ≅ 𝒪_K^×`
* `Gal(K_π · K^ur / K) ≅ 𝒪_K^× × Gal(K^ur/K)`

★すべて `∃` の内側に閉じ込めてある(自由なパラメータを結論に
出さない、という設計上の制約を守るため)。 -/
theorem exists_lubinTateUnramified_decomposition (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧
      E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      Nonempty
        (((E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure) ≃ₐ[K.carrier]
            (E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure)) ≃*
          ((𝒪[K.carrier])ˣ × (unramifiedClosure K ≃ₐ[K.carrier] unramifiedClosure K))) := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨π, hπirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer π).mp hπirr
  have hπne0 : π ≠ 0 := hπirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  exact ⟨lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    lubinTateClosure_inf_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf,
    ⟨lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf⟩,
    ⟨lubinTateUnramifiedGalEquivProd K hq hπmax hπne0 f hf0 hf1 hf⟩⟩

end ABC3.Found.PGC
