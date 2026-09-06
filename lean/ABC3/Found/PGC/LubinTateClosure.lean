import ABC3.Found.PGC.LubinTateReciprocityMapLimitKernel
import ABC3.Found.PGC.LubinTateNormalAdjoin
import ABC3.Found.PGC.AbsGalUnitsSurjective

/-!
# `K_π := ⋃_n K(Λ_n)` と `Gal(K_π/K) ≅ 𝒪_K^×`(`sorry` 無し)

経路 Λ の節点 Λ3。Lubin-Tate 塔の**極限側**を、体のオブジェクト
`K_π`(`K.closure` の中間体)として実際に構成し、その Galois 群を
単数群 `𝒪_K^×` と同一視する。

## 何が在庫にあって、何を足したか

在庫(すべて `sorry` 無し):

| 部品 | 供給元 |
|---|---|
| `reciprocityMapLimit : Gal(K̄/K) →* CompatibleUnits` | `LubinTateReciprocityMapLimit.lean` |
| その全射性 `reciprocityMapLimit_surjective` | `LubinTateReciprocityMapLimitSurjective.lean` |
| その核の特徴づけ `mem_ker_reciprocityMapLimit_iff` | `LubinTateReciprocityMapLimitKernel.lean` |
| `(𝒪_K)^× ≃* CompatibleUnits`(`unitsEquivCompatibleUnits`) | `AbsGalSurjections.lean` |
| compatible な生成元列 `psiGenSeq` | `LubinTateGeneratorSequence.lean` |
| `K(x)/K` の正規性 | `LubinTateNormalAdjoin.lean` |
| `Λ_n` の原始的な点は生成元 1 つで生成される | `LubinTateActionEquivariance.lean` |

★核の特徴づけ `mem_ker_reciprocityMapLimit_iff` は「`σ` が生成元列
`(psiGenSeq m).pt` を**すべて固定する**」という形で、体 `K_π` を
構成せずに `Gal(K̄/K_π)` を述べていた。本ファイルはその `K_π` を
実際に作り、**核 = `K_π` の固定部分群**であることを示して、
`Gal(K_π/K) ≅ 𝒪_K^×` を出す。

足したのは:

1. `mem_fixingSubgroup_adjoin_iff`(一般の体論): 生成集合の上で
   固定すれば `adjoin` 全体を固定する。
2. `lubinTateClosure`(= `K_π`)の定義と、`K_π = K(x_1, x_2, …)`
   (`lubinTateClosure_eq_adjoin_genSet`)。
3. `normal_lubinTateClosure`: `K_π/K` は正規。
4. `ker_reciprocityMapLimit_eq_fixingSubgroup`: 核 = `Gal(K̄/K_π)`。
5. `lubinTateClosureGalEquivUnits`: `Gal(K_π/K) ≃* 𝒪_K^×`。

## 定義の読み方(原典との対応)

原典の `K_π := ⋃_n K(Λ_n)`(体の増大列の合併)を、ここでは
`K(⋃_n Λ_n)`(全捩れ点で生成した体)として定義する
(`lubinTateClosure`)。合併される体の族は増大列なので両者は一致する
——形式化では後者のほうが `IntermediateField.adjoin` 1 つで書けて
扱いやすい。増大性そのものは `monotone_adjoin_psiGenSeq` と
`exists_mem_adjoin_psiGenSeq_of_mem_lubinTateClosure` で明示する。

## 逸脱(記録)

* **`f` への依存**: `lubinTateClosure` は一意化元 `π` だけでなく
  Lubin-Tate 級数 `f` にも依存する形で定義されている。古典論では
  `K_π` は `f` に依らないが、その独立性は Λ6(Lubin-Tate-Rosen)の
  担当なので、ここでは主張しない。
* **`Gal(K_π/K) ≅ 𝒪_K^×` は群同型としてのみ**。位相群としての同型
★★**2026-09-06 解消**: 位相版は `Found/PGC/LubinTateClosureTopology.lean` の
`lubinTateClosureGalContinuousMulEquivUnits` で閉じた(sorry 0)。
★`toMulEquiv` は本ファイルの `lubinTateClosureGalEquivUnits` そのものである。
  (両辺の副有限位相/`π` 進位相が対応すること)は本ファイルでは
  主張していない——`𝒪_K^×` の位相を `lim (𝒪_K/π^n)^×` と同一視する
  部品がまだ在庫に無い。残った節点として報告する。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 一般の体論の補題 -/

/-- **生成集合の上で固定すれば `adjoin` 全体を固定する**。

`σ` が固定する元全体は中間体(`IntermediateField.fixedField
(Subgroup.closure {σ})`)なので、`adjoin_le_iff` がそのまま効く。 -/
theorem mem_fixingSubgroup_adjoin_iff {F E : Type*} [Field F] [Field E] [Algebra F E]
    (S : Set E) (σ : E ≃ₐ[F] E) :
    σ ∈ (IntermediateField.adjoin F S).fixingSubgroup ↔ ∀ x ∈ S, σ x = x := by
  constructor
  · intro hσ x hx
    exact (IntermediateField.mem_fixingSubgroup_iff _ σ).mp hσ x
      (IntermediateField.subset_adjoin F S hx)
  · intro hS
    rw [IntermediateField.mem_fixingSubgroup_iff]
    have hle : IntermediateField.adjoin F S ≤
        IntermediateField.fixedField (Subgroup.closure ({σ} : Set (E ≃ₐ[F] E))) := by
      rw [IntermediateField.adjoin_le_iff]
      intro x hx
      simp only [SetLike.mem_coe, IntermediateField.mem_fixedField_iff]
      intro g hg
      induction hg using Subgroup.closure_induction with
      | mem y hy => rw [Set.mem_singleton_iff] at hy; subst hy; exact hS x hx
      | one => rfl
      | mul a b _ _ ha hb => show a (b x) = x; rw [hb, ha]
      | inv a _ ha =>
          show a⁻¹ x = x
          conv_lhs => rw [← ha]
          exact a.symm_apply_apply x
    intro x hx
    have hx' := hle hx
    simp only [IntermediateField.mem_fixedField_iff] at hx'
    exact hx' σ (Subgroup.subset_closure rfl)

/-- **`E` への制限が自明 ⟺ `E` を各点固定する**(一般の体の塔で)。
`Found/NumberField/SplCompositum.lean::restrictNormalHom_eq_one_iff`
(`ℚ` 固定版)の一般化。 -/
theorem restrictNormalHom_eq_one_iff' {F Ω : Type*} [Field F] [Field Ω] [Algebra F Ω]
    (E : IntermediateField F Ω) [Normal F E] (σ : Ω ≃ₐ[F] Ω) :
    AlgEquiv.restrictNormalHom (E : Type _) σ = 1 ↔ σ ∈ E.fixingSubgroup := by
  rw [IntermediateField.mem_fixingSubgroup_iff]
  constructor
  · intro h x hx
    have hh := AlgEquiv.restrictNormal_commutes σ (E : Type _) ⟨x, hx⟩
    have h2 : (σ.restrictNormal (E : Type _)) ⟨x, hx⟩ = ⟨x, hx⟩ := DFunLike.congr_fun h ⟨x, hx⟩
    rw [h2] at hh
    exact hh.symm
  · intro h
    refine AlgEquiv.ext fun y => ?_
    refine FaithfulSMul.algebraMap_injective (E : Type _) Ω ?_
    have hh := AlgEquiv.restrictNormal_commutes σ (E : Type _) y
    show algebraMap (E : Type _) Ω ((σ.restrictNormal (E : Type _)) y) = algebraMap (E : Type _) Ω y
    rw [hh]
    exact h (algebraMap (E : Type _) Ω y) y.2

/-- `ker(restrictNormalHom E) = E.fixingSubgroup`。 -/
theorem ker_restrictNormalHom_eq_fixingSubgroup {F Ω : Type*} [Field F] [Field Ω] [Algebra F Ω]
    (E : IntermediateField F Ω) [Normal F E] :
    MonoidHom.ker (AlgEquiv.restrictNormalHom (F := F) (K₁ := Ω) (E := (E : Type _))) =
      E.fixingSubgroup := by
  ext σ
  exact restrictNormalHom_eq_one_iff' E σ

/-! ## 1. `K_π` の定義 -/

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

/-- 全段の捩れ点の合併 `Λ_∞ = ⋃_n Λ_n`。 -/
noncomputable def lubinTateTorsionSet : Set K.closure :=
  ⋃ n : ℕ,
    ((iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n : Finset K.closure) :
      Set K.closure)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`K_π := K(Λ_∞) = ⋃_n K(Λ_n)`**——Lubin-Tate 塔の合併。 -/
noncomputable def lubinTateClosure : IntermediateField K.carrier K.closure :=
  IntermediateField.adjoin K.carrier (lubinTateTorsionSet K hq hπmax hπne0 f hf0 hf1 hf)

/-- compatible な生成元列がなす集合 `{x_1, x_2, …}`。 -/
noncomputable def lubinTateGenSet : Set K.closure :=
  Set.range (fun m : ℕ => (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt)

/-! ## 2. `K_π = K(x_1, x_2, …)` -/

/-- `Λ_0 = {0}`——`|Λ_0| = q^0 = 1` と `0 ∈ Λ_0` から。 -/
theorem iteratedLubinTateTorsionPoints_zero :
    iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf 0 = {0} := by
  have hcard := card_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf 0
  rw [pow_zero] at hcard
  exact (Finset.card_eq_one.mp hcard).elim (fun a ha => by
    have h0 := zero_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf 0
    rw [ha, Finset.mem_singleton] at h0
    rw [ha, h0])

/-- `K(x_{m+1}) ≤ K(x_1, x_2, …)`。 -/
theorem adjoin_psiGenSeq_le_adjoin_genSet (m : ℕ) :
    IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) ≤
      IntermediateField.adjoin K.carrier (lubinTateGenSet K hq hπmax hπne0 f hf0 hf1 hf) :=
  IntermediateField.adjoin_simple_le_iff.mpr
    (IntermediateField.subset_adjoin K.carrier _ ⟨m, rfl⟩)

/-- **`Λ_n ⊆ K(x_1, x_2, …)`**——`n` についての帰納法。
`Λ_{n+1} = Λ_n ∪ (Λ_{n+1} の原始的な部分)`
(`iteratedLubinTateTorsionPoints_eq_union`)で分け、原始的な側は
`iteratedLubinTatePsiTorsionPoints_subset_adjoin`(同じ段の原始的な
点は生成元 1 つで生成される)に渡す。 -/
theorem iteratedLubinTateTorsionPoints_subset_adjoin_genSet (n : ℕ) :
    ∀ x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
      x ∈ IntermediateField.adjoin K.carrier (lubinTateGenSet K hq hπmax hπne0 f hf0 hf1 hf) := by
  induction n with
  | zero =>
    intro x hx
    rw [iteratedLubinTateTorsionPoints_zero, Finset.mem_singleton] at hx
    subst hx
    exact zero_mem _
  | succ m ih =>
    intro x hx
    rw [iteratedLubinTateTorsionPoints_eq_union K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega),
      Finset.mem_union] at hx
    rcases hx with hx | hx
    · exact ih x (by simpa using hx)
    · haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
      exact adjoin_psiGenSeq_le_adjoin_genSet K hq hπmax hπne0 f hf0 hf1 hf m
        (iteratedLubinTatePsiTorsionPoints_subset_adjoin K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
          (by omega) _ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem x hx)

/-- 生成元は捩れ点。 -/
theorem lubinTateGenSet_subset_torsionSet :
    lubinTateGenSet K hq hπmax hπne0 f hf0 hf1 hf ⊆
      lubinTateTorsionSet K hq hπmax hπne0 f hf0 hf1 hf := by
  rintro _ ⟨m, rfl⟩
  exact Set.mem_iUnion.mpr ⟨m + 1, (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn⟩

/-- ★**`K_π = K(x_1, x_2, …)`**——全捩れ点で生成しても、compatible な
生成元列だけで生成しても同じ体。 -/
theorem lubinTateClosure_eq_adjoin_genSet :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf =
      IntermediateField.adjoin K.carrier (lubinTateGenSet K hq hπmax hπne0 f hf0 hf1 hf) := by
  refine le_antisymm ?_ ?_
  · rw [lubinTateClosure, IntermediateField.adjoin_le_iff]
    intro x hx
    obtain ⟨n, hxn⟩ := Set.mem_iUnion.mp hx
    exact iteratedLubinTateTorsionPoints_subset_adjoin_genSet K hq hπmax hπne0 f hf0 hf1 hf n x hxn
  · exact IntermediateField.adjoin.mono _ _ _
      (lubinTateGenSet_subset_torsionSet K hq hπmax hπne0 f hf0 hf1 hf)

/-- ★`K_π = ⨆_m K(x_{m+1})`——増大する有限次拡大の上限としての表示。 -/
theorem lubinTateClosure_eq_iSup :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf =
      ⨆ m : ℕ, IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) := by
  rw [lubinTateClosure_eq_adjoin_genSet, lubinTateGenSet, Set.range_eq_iUnion,
    IntermediateField.adjoin_iUnion]

/-! ## 3. `K_π/K` は正規 -/

/-- ★**`K_π/K` は正規**——各 `K(x_{m+1})` が正規
(`normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints`)で、
`K_π` はそれらの上限(`IntermediateField.normal_iSup`)。 -/
theorem normal_lubinTateClosure :
    Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) := by
  rw [lubinTateClosure_eq_iSup]
  haveI : ∀ m : ℕ, Normal K.carrier
      (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) := by
    intro m
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
    exact normal_adjoin_of_mem_iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf
      (m + 1) (by omega) _ (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem
  exact IntermediateField.normal_iSup K.carrier K.closure _

/-! ## 4. 核 = `Gal(K̄/K_π)`、そして `Gal(K_π/K) ≅ 𝒪_K^×` -/

/-- ★★**`ker(reciprocityMapLimit) = Gal(K̄/K_π)`**——
`mem_ker_reciprocityMapLimit_iff`(生成元をすべて固定する)を、
`K_π` の固定部分群の言葉に翻訳しただけ。 -/
theorem ker_reciprocityMapLimit_eq_fixingSubgroup :
    MonoidHom.ker (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf) =
      (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup := by
  ext σ
  rw [mem_ker_reciprocityMapLimit_iff, lubinTateClosure_eq_adjoin_genSet,
    mem_fixingSubgroup_adjoin_iff]
  constructor
  · rintro h _ ⟨m, rfl⟩
    exact h m
  · intro h m
    exact h _ ⟨m, rfl⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K_π/K) ≅ 𝒪_K^×`**(群同型)。

`Gal(K̄/K) ↠ Gal(K_π/K)`(`restrictNormalHom`、`K_π/K` 正規)と
`Gal(K̄/K) ↠ CompatibleUnits`(`reciprocityMapLimit`)の**核が一致**
することから、第一同型定理を 2 回使って得る。最後に
`unitsEquivCompatibleUnits`(`𝒪_K^× ≃* CompatibleUnits`)で
単数群へ戻す。 -/
noncomputable def lubinTateClosureGalEquivUnits :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
      lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) ≃* (𝒪[K.carrier])ˣ := by
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  have hker : MonoidHom.ker (AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (E := (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _))) =
      MonoidHom.ker (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf) := by
    rw [ker_restrictNormalHom_eq_fixingSubgroup, ker_reciprocityMapLimit_eq_fixingSubgroup]
  exact
    (QuotientGroup.quotientKerEquivOfSurjective _
      (AlgEquiv.restrictNormalHom_surjective (F := K.carrier)
        (K₁ := (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)) K.closure)).symm.trans
      ((QuotientGroup.quotientMulEquivOfEq hker).trans
        ((QuotientGroup.quotientKerEquivOfSurjective _
          (reciprocityMapLimit_surjective K hq hπmax hπne0 f hf0 hf1 hf)).trans
          (unitsEquivCompatibleUnits K hπmax).symm))

/-! ## 5. 塔の増大性(Λ4 で使う) -/

/-- `x_{m+1} ∈ K(x_{m+2})`——`π·x_{m+2} = x_{m+1}`(`psiGenSeq_compat`)の
右辺が `K(x_{m+2})` の元の像だから。 -/
theorem psiGenSeq_pt_mem_adjoin_succ (m : ℕ) :
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt ∈
      IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf (m + 1)).pt} : Set K.closure) := by
  have hc := psiGenSeq_compat K hq hπmax hπne0 f hf0 hf1 hf m
  rw [← hc]
  exact SetLike.coe_mem _

/-- ★塔 `K(x_1) ⊆ K(x_2) ⊆ …` は増大列。 -/
theorem monotone_adjoin_psiGenSeq :
    Monotone (fun m : ℕ => IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) :=
  monotone_nat_of_le_succ fun m =>
    IntermediateField.adjoin_simple_le_iff.mpr
      (psiGenSeq_pt_mem_adjoin_succ K hq hπmax hπne0 f hf0 hf1 hf m)

/-- ★`K_π` の元はどれか 1 つの `K(x_{m+1})` に入る(塔が増大列だから)。 -/
theorem exists_mem_adjoin_psiGenSeq_of_mem_lubinTateClosure
    {x : K.closure} (hx : x ∈ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :
    ∃ m : ℕ, x ∈ IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) := by
  rw [lubinTateClosure_eq_iSup] at hx
  have hdir := (monotone_adjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf).directed_le
  have hcoe := IntermediateField.coe_iSup_of_directed hdir
  have hx' : x ∈ (↑(⨆ m : ℕ, IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) : Set K.closure) := hx
  rw [hcoe] at hx'
  exact Set.mem_iUnion.mp hx'

/-! ## 6. 仮説なしの形 -/

/-- ★★★★**仮説なしの Λ3**: 任意の p 進局所体 `K` について、
`K̄` の中に `K` 上正規な中間体 `K_π` が在って `Gal(K_π/K) ≅ 𝒪_K^×`。

仮説(adic 完備性・剰余体の有限性と位数・一意化元・Lubin-Tate 級数)は
すべて本リポジトリで既に構築済みなので、`AbsGalUnitsSurjective.lean`
と同じ組み上げで消える。 -/
theorem exists_lubinTateClosure_galEquiv_units (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧ Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) := by
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
    ⟨lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf⟩⟩

end ABC3.Found.PGC
