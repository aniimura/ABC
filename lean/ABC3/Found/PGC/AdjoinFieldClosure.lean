import ABC3.Found.PGC.SubgroupCorrespondenceConstruction

/-!
# 代数閉包の同一視 `(K(x))‾ ≃ K‾`

`adjoinField K x` を `PAdicLocalField` として見ると、その代数閉包は
`AlgebraicClosure ↥K⟮x⟯` という**別の型**になる。しかし `K.closure` 自身が
`K(x)` 上の代数閉包でもあるので、両者は同型である。

これが要る理由: Lubin-Tate 塔は `K(Λ_n)/K(Λ_{n-1})` という**相対的な**拡大の
連なりとして与えられる(`ψ_n` は `𝒪_{K(Λ_{n-1})}` 上の Eisenstein)。
一方 `Found/PGC/TotallyRamified.lean` の完全分岐の道具は
`K.closure` の中の中間体について書かれている。両者を繋ぐには
`(K(Λ_{n-1}))‾ ≃ K‾` が要る。

## 配管の注意

`Algebra ↥K⟮x⟯ K.closure` などのインスタンスは、**中間体を直接書けば**
見つかるが `(adjoinField K x).carrier` と書くと見つからない
(`adjoinField` は `def` なので instance 探索が展開しない)。
そこで明示的に `instance` として与え直す。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-- `K.closure` は `K(x)`-代数(中間体を直接書けば見つかるインスタンスの言い換え)。 -/
noncomputable instance algebraAdjoinFieldClosure (K : PAdicLocalField p) (x : K.closure) :
    Algebra ((adjoinField K x).carrier) K.closure :=
  (inferInstance : Algebra (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure)

instance isAlgebraicAdjoinFieldClosure (K : PAdicLocalField p) (x : K.closure) :
    Algebra.IsAlgebraic ((adjoinField K x).carrier) K.closure :=
  (inferInstance :
    Algebra.IsAlgebraic (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure)

/-- `K.closure` は `K(x)` 上の**代数閉包**でもある。 -/
instance isAlgClosureAdjoinField (K : PAdicLocalField p) (x : K.closure) :
    IsAlgClosure ((adjoinField K x).carrier) K.closure := by
  constructor
  · infer_instance
  · infer_instance

/-- **★★★★★代数閉包の同一視** `(K(x))‾ ≃ K‾`(`K(x)`-代数として)。

これで「`K(x)` 上の拡大」を `K.closure` の中の中間体として捉え直せる
——Lubin-Tate 塔の相対的な段を `K` からの絶対的な塔に繋ぐ鍵。 -/
noncomputable def closureEquivAdjoinField (K : PAdicLocalField p) (x : K.closure) :
    (adjoinField K x).closure ≃ₐ[(adjoinField K x).carrier] K.closure :=
  IsAlgClosure.equiv _ _ _

/-! ## `Γ_{L}` は `Γ_K` の開部分群

`Interface/PGC/LocalFieldData.lean::SubgroupCorrespondence.waiting` は
待っているものを二つ挙げていた:

> 開部分群 H ⊆ Γ_K に対応する中間体が p進局所体であり、
> **その絶対 Galois 群が H であること**

前者は `Found/PGC/SubgroupCorrespondenceConstruction.lean::fixedFieldLocalField`
で解消済み。**後者を本節で解消する**——代数閉包の同一視
`closureEquivAdjoinField` と mathlib の
`IntermediateField.fixingSubgroupEquiv` を繋ぐだけ。 -/

/-- **★★★★★★`Γ_{K(x)}` は `K(x)` の固定部分群**。 -/
noncomputable def absGalAdjoinFieldEquiv (K : PAdicLocalField p) (x : K.closure) :
    (adjoinField K x).absGal ≃*
      ((IntermediateField.adjoin K.carrier ({x} : Set K.closure)).fixingSubgroup) :=
  (AlgEquiv.autCongr (closureEquivAdjoinField K x)).trans
    (IntermediateField.fixingSubgroupEquiv
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure))).symm

noncomputable instance algebraFixedFieldClosure (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    Algebra ((fixedFieldLocalField K H hH).carrier) K.closure :=
  (inferInstance : Algebra (IntermediateField.fixedField H) K.closure)

instance isAlgebraicFixedFieldClosure (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    Algebra.IsAlgebraic ((fixedFieldLocalField K H hH).carrier) K.closure :=
  (inferInstance : Algebra.IsAlgebraic (IntermediateField.fixedField H) K.closure)

instance isAlgClosureFixedField (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    IsAlgClosure ((fixedFieldLocalField K H hH).carrier) K.closure := by
  constructor
  · infer_instance
  · infer_instance

noncomputable def closureEquivFixedField (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    (fixedFieldLocalField K H hH).closure
      ≃ₐ[(fixedFieldLocalField K H hH).carrier] K.closure :=
  IsAlgClosure.equiv _ _ _

/-- **★★★★★★★開部分群 `H` に対応する体 `L_H` の絶対 Galois 群は `H` そのもの**。

`Interface/PGC/LocalFieldData.lean::SubgroupCorrespondence.waiting` が
待っていたもう半分——これで `SubgroupCorrespondence` の待ちは解消した。 -/
noncomputable def absGalFixedFieldEquiv (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    (fixedFieldLocalField K H hH).absGal ≃* H := by
  haveI := isGalois_closure K
  have hclosed : IsClosed (H : Set K.absGal) := Subgroup.isClosed_of_isOpen H hH
  have h2 : (IntermediateField.fixedField H).fixingSubgroup = H :=
    InfiniteGalois.fixingSubgroup_fixedField (⟨H, hclosed⟩ : ClosedSubgroup K.absGal)
  exact ((AlgEquiv.autCongr (closureEquivFixedField K H hH)).trans
    (IntermediateField.fixingSubgroupEquiv (IntermediateField.fixedField H)).symm).trans
    (MulEquiv.subgroupCongr h2)

/-! ## `AlgEquiv.autCongr` は Krull 位相で連続

`Found/PGC/GaloisTransferContinuous.lean::continuous_galMulEquivOf`(第 969)は
「二つの p進局所体の**台**の同型から誘導される `Γ_K ≃ Γ_{K'}` は連続」を
言っていた。ここで要るのはその双対——**台は同じで代数閉包が違う**場合。

一般に、`F`-代数の同型 `e : Ω₁ ≃ₐ[F] Ω₂` から誘導される
`AlgEquiv.autCongr e : Gal(Ω₁/F) ≃* Gal(Ω₂/F)` は連続である。
証明は第 969 と同じ形だが、原始元定理の代わりに
`IntermediateField.map` で有限次中間体を引き戻せばよく、**より短い**。 -/

/-- **★★★★★`AlgEquiv.autCongr e` は連続**——`F`-代数の同型 `e : Ω₁ ≃ₐ[F] Ω₂`
から誘導される `Gal(Ω₁/F) → Gal(Ω₂/F)`。 -/
theorem continuous_autCongr {F Ω₁ Ω₂ : Type*} [Field F] [Field Ω₁] [Field Ω₂]
    [Algebra F Ω₁] [Algebra F Ω₂] (e : Ω₁ ≃ₐ[F] Ω₂) :
    Continuous (AlgEquiv.autCongr e) := by
  refine continuous_of_continuousAt_one (AlgEquiv.autCongr e).toMonoidHom ?_
  rw [ContinuousAt, map_one, Filter.tendsto_def]
  intro s hs
  obtain ⟨E₂, hE₂fin, hE₂sub⟩ := (krullTopology_mem_nhds_one_iff F Ω₂ s).mp hs
  haveI := hE₂fin
  rw [krullTopology_mem_nhds_one_iff]
  refine ⟨IntermediateField.map (e.symm : Ω₂ →ₐ[F] Ω₁) E₂, ?_, ?_⟩
  · exact (IntermediateField.equivMap E₂ (e.symm : Ω₂ →ₐ[F] Ω₁)).toLinearEquiv.finiteDimensional
  · intro g hg
    apply hE₂sub
    rw [SetLike.mem_coe] at hg ⊢
    rw [IntermediateField.mem_fixingSubgroup_iff] at hg ⊢
    intro z hz
    have hgz : g ((e.symm : Ω₂ →ₐ[F] Ω₁) z) = (e.symm : Ω₂ →ₐ[F] Ω₁) z := hg _ ⟨z, hz, rfl⟩
    show e (g (e.symm z)) = z
    rw [show g (e.symm z) = e.symm z from hgz, AlgEquiv.apply_symm_apply]

/-- **★★★★★`autCongr` は位相群の同型**。 -/
noncomputable def autCongrContinuousMulEquiv {F Ω₁ Ω₂ : Type*} [Field F] [Field Ω₁] [Field Ω₂]
    [Algebra F Ω₁] [Algebra F Ω₂] (e : Ω₁ ≃ₐ[F] Ω₂) :
    ContinuousMulEquiv (Ω₁ ≃ₐ[F] Ω₁) (Ω₂ ≃ₐ[F] Ω₂) where
  toMulEquiv := AlgEquiv.autCongr e
  continuous_toFun := continuous_autCongr e
  continuous_invFun := by
    show Continuous (⇑(AlgEquiv.autCongr e).symm)
    rw [AlgEquiv.autCongr_symm]
    exact continuous_autCongr e.symm

/-- 代数閉包の同一視は**位相群**の同型を与える:
`Γ_{K(x)} = Gal((K(x))‾/K(x)) ≃ₜ* Gal(K‾/K(x))`。 -/
noncomputable def absGalAdjoinFieldContinuousEquiv (K : PAdicLocalField p) (x : K.closure) :
    ContinuousMulEquiv ((adjoinField K x).closure ≃ₐ[(adjoinField K x).carrier]
        (adjoinField K x).closure)
      (K.closure ≃ₐ[(adjoinField K x).carrier] K.closure) :=
  autCongrContinuousMulEquiv (closureEquivAdjoinField K x)

/-- 同じく固定体の側。 -/
noncomputable def absGalFixedFieldContinuousEquiv (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    ContinuousMulEquiv ((fixedFieldLocalField K H hH).closure
        ≃ₐ[(fixedFieldLocalField K H hH).carrier] (fixedFieldLocalField K H hH).closure)
      (K.closure ≃ₐ[(fixedFieldLocalField K H hH).carrier] K.closure) :=
  autCongrContinuousMulEquiv (closureEquivFixedField K H hH)

/-! ## `fixingSubgroupEquiv` は連続(有限次中間体の場合)

`IntermediateField.fixingSubgroupEquiv E : E.fixingSubgroup ≃* Gal(Ω/E)` は
集合としては「同じ自己同型を見る向きを変える」だけだが、位相は
**部分空間位相(`Gal(Ω/F)` から)** と **Krull 位相(`Gal(Ω/E)` の)** で
別物なので、同相であることは示す必要がある。

`E/F` が有限次のときは順方向が短く済む: `E'/E` 有限次なら
`E'.restrictScalars F` は `F` 上も有限次(塔)なので、その固定部分群は
`Gal(Ω/F)` で**開**であり、その引き戻しが求める近傍になる。 -/

/-- **★★★★★`fixingSubgroupEquiv` は連続**(`E/F` 有限次)。 -/
theorem continuous_fixingSubgroupEquiv {F Ω : Type*} [Field F] [Field Ω] [Algebra F Ω]
    (E : IntermediateField F Ω) [FiniteDimensional F E] :
    Continuous (E.fixingSubgroupEquiv) := by
  refine continuous_of_continuousAt_one (E.fixingSubgroupEquiv).toMonoidHom ?_
  rw [ContinuousAt, map_one, Filter.tendsto_def]
  intro s hs
  obtain ⟨E', hE'fin, hE'sub⟩ := (krullTopology_mem_nhds_one_iff (↥E) Ω s).mp hs
  haveI := hE'fin
  haveI : FiniteDimensional F (E'.restrictScalars F) :=
    (Module.Finite.trans (R := F) (↥E) (↥E') : FiniteDimensional F ↥E')
  have hopen : IsOpen ((E'.restrictScalars F).fixingSubgroup : Set (Ω ≃ₐ[F] Ω)) :=
    IntermediateField.fixingSubgroup_isOpen _
  have hnhd : (Subtype.val ⁻¹' ((E'.restrictScalars F).fixingSubgroup : Set (Ω ≃ₐ[F] Ω)))
      ∈ nhds (1 : E.fixingSubgroup) := by
    refine IsOpen.mem_nhds (hopen.preimage continuous_subtype_val) ?_
    show (1 : Ω ≃ₐ[F] Ω) ∈ (E'.restrictScalars F).fixingSubgroup
    exact one_mem _
  refine Filter.mem_of_superset hnhd ?_
  intro σ hσ
  apply hE'sub
  rw [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff]
  intro z hz
  exact (IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ (z : Ω) hz

/-- **★★★★★逆向きも連続**。

`F` 上有限次の `E''` に対して `E` 上有限次の `E'` を取ればよい。原始元を
取る必要はない——**合成体 `E ⊔ E''`** を `extendScalars` で `E` 上の中間体と
見れば、`E/F`・`E''/F` がともに有限次なので `E ⊔ E''` も有限次
(`IntermediateField.finiteDimensional_sup`)、したがって `E` 上も有限次。 -/
theorem continuous_fixingSubgroupEquiv_symm {F Ω : Type*} [Field F] [Field Ω] [Algebra F Ω]
    (E : IntermediateField F Ω) [FiniteDimensional F E] :
    Continuous (E.fixingSubgroupEquiv.symm) := by
  refine continuous_of_continuousAt_one (E.fixingSubgroupEquiv.symm).toMonoidHom ?_
  rw [ContinuousAt, map_one, Filter.tendsto_def]
  intro s hs
  rw [nhds_subtype, Filter.mem_comap] at hs
  obtain ⟨t, ht, hts⟩ := hs
  obtain ⟨E'', hfin, hsub⟩ := (krullTopology_mem_nhds_one_iff F Ω t).mp ht
  haveI := hfin
  set E' : IntermediateField (↥E) Ω :=
    IntermediateField.extendScalars (le_sup_left : E ≤ E ⊔ E'') with hE'
  haveI hF : FiniteDimensional F ↥E' := by
    show FiniteDimensional F ↥(IntermediateField.extendScalars (le_sup_left : E ≤ E ⊔ E''))
    exact (inferInstance : FiniteDimensional F ↥(E ⊔ E''))
  haveI : FiniteDimensional (↥E) ↥E' := Module.Finite.of_restrictScalars_finite F (↥E) ↥E'
  rw [krullTopology_mem_nhds_one_iff]
  refine ⟨E', inferInstance, ?_⟩
  intro σ hσ
  apply hts
  show (E.fixingSubgroupEquiv.symm σ : Ω ≃ₐ[F] Ω) ∈ t
  apply hsub
  rw [SetLike.mem_coe, IntermediateField.mem_fixingSubgroup_iff]
  intro z hz
  have hmem : z ∈ E' := by
    rw [hE', IntermediateField.mem_extendScalars]
    exact (le_sup_right : E'' ≤ E ⊔ E'') hz
  exact (IntermediateField.mem_fixingSubgroup_iff _ _).mp hσ z hmem

/-- **`E.fixingSubgroup ≃ₜ* Gal(Ω/E)`**——部分空間位相と Krull 位相の一致。 -/
noncomputable def fixingSubgroupContinuousMulEquiv {F Ω : Type*} [Field F] [Field Ω] [Algebra F Ω]
    (E : IntermediateField F Ω) [FiniteDimensional F E] :
    ContinuousMulEquiv (E.fixingSubgroup) (Ω ≃ₐ[E] Ω) where
  toMulEquiv := E.fixingSubgroupEquiv
  continuous_toFun := continuous_fixingSubgroupEquiv E
  continuous_invFun := continuous_fixingSubgroupEquiv_symm E

/-- 等しい部分群の間の同一視は同相。 -/
noncomputable def subgroupCongrContinuousMulEquiv {G : Type*} [Group G] [TopologicalSpace G]
    {H₁ H₂ : Subgroup G} (h : H₁ = H₂) : ContinuousMulEquiv H₁ H₂ where
  toMulEquiv := MulEquiv.subgroupCongr h
  continuous_toFun := continuous_induced_rng.2 continuous_subtype_val
  continuous_invFun := continuous_induced_rng.2 continuous_subtype_val

/-- **★★★★★★★★★★`Γ_{L_H} ≃ₜ* H`——位相群の同型として**。

`absGalFixedFieldEquiv`(群同型)の位相版。第 992 の
`inertia_recoverable_of_residueCard_transport` が第二仮定として要求していた形。

三段の合成:
1. `absGalFixedFieldContinuousEquiv`(代数閉包の同一視、第 993)
2. `fixingSubgroupContinuousMulEquiv` の逆(部分空間位相 = Krull 位相)
3. `InfiniteGalois.fixingSubgroup_fixedField`(閉部分群の Galois 対応) -/
noncomputable def absGalFixedFieldCME (K : PAdicLocalField p) (H : Subgroup K.absGal)
    (hH : IsOpen (H : Set K.absGal)) :
    ContinuousMulEquiv (fixedFieldLocalField K H hH).absGal H := by
  haveI := isGalois_closure K
  haveI := finiteDimensional_fixedField_of_isOpen K H hH
  have hclosed : IsClosed (H : Set K.absGal) := Subgroup.isClosed_of_isOpen H hH
  have h2 : (IntermediateField.fixedField H).fixingSubgroup = H :=
    InfiniteGalois.fixingSubgroup_fixedField (⟨H, hclosed⟩ : ClosedSubgroup K.absGal)
  exact ((absGalFixedFieldContinuousEquiv K H hH).trans
    (fixingSubgroupContinuousMulEquiv (IntermediateField.fixedField H)).symm).trans
    (subgroupCongrContinuousMulEquiv h2)

end ABC3.Found.PGC
