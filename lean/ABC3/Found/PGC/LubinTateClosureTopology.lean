import ABC3.Found.PGC.LubinTateClosure
import Mathlib.FieldTheory.Galois.Profinite

/-!
# `Gal(K_π/K) ≃ₜ* 𝒪_K^×`(位相群としての同型、`sorry` 無し)

経路 Λ の節点 Λ3 の**未達分**。`Found/PGC/LubinTateClosure.lean` は
`lubinTateClosureGalEquivUnits : Gal(K_π/K) ≃* 𝒪_K^×` を**群同型として
のみ**出していた(同ファイルの「逸脱(記録)」に明記されている)。本
ファイルはそれを**位相群の同型 `≃ₜ*`(`ContinuousMulEquiv`)へ持ち上げる**。

## 何が在庫にあって、何を足したか

在庫(すべて `sorry` 無し):

| 部品 | 供給元 |
|---|---|
| `lubinTateClosureGalEquivUnits : Gal(K_π/K) ≃* 𝒪_K^×` | `LubinTateClosure.lean` |
| `mem_fixingSubgroup_adjoin_iff` | `LubinTateClosure.lean` |
| `principalUnits`(主単数群)・`principalUnits_eq_ker` | `AdjoinIntegers.lean` |
| `principalUnits_succ_le`(単調減少) | `LubinTateTowerCompatible.lean` |
| `unitsEquivCompatibleUnits : 𝒪_K^× ≃* CompatibleUnits` | `AbsGalSurjections.lean` |
| `reciprocityMapLimit` とその全射性 | `LubinTateReciprocityMapLimit*.lean` |
| `compactSpace_integer : CompactSpace 𝒪[K]` | `Found/ResidueFieldFinite.lean` |

★**mathlib に在った(見積を大きく下回った理由)**:

* `instance [IsGalois k K] : CompactSpace Gal(K/k)`
  (`Mathlib/FieldTheory/Galois/Profinite.lean`)——`Gal(K_π/K)` のコンパクト性。
* `Continuous.homeoOfEquivCompactToT2`——コンパクト → T2 の連続全単射は同相。
* `instance [T1Space α] [ContinuousMul α] [CompactSpace α] : CompactSpace αˣ`
  (`Mathlib/Topology/Algebra/Group/Basic.lean`)——使わなかったが、逆向きの
  経路(`𝒪_K^× → Gal(K_π/K)` の連続性)を選んでも通る。
* `Units.isEmbedding_embedProduct`——単数群の位相の扱い。
* `IntermediateField.fixingSubgroup_isOpen`——有限次中間体の固定部分群は開。

足したのは:

1. `𝒪_K` の位相と `π` 進フィルターの対応(§1):
   `norm_le_of_mem_span_pow` / `mem_span_pow_of_norm_le`(`x ∈ π^N𝒪_K ↔
   ‖x‖ ≤ ‖π‖^N`)、`norm_pi_lt_one`、
   **`isOpen_span_pow`**(`π^N𝒪_K` は開)、
   **`isOpen_principalUnits`**(主単数群 `1+π^N𝒪_K` は `𝒪_K^×` で開)、
   **`exists_principalUnits_subset_of_mem_nhds`**(主単数群は `𝒪_K^×` の
   `1` の**基本近傍系**をなす)。
2. `reciprocityMapLimitFamily_succ_eq_one_iff`(§2): 相互律の**層ごと**の核。
   在庫の `mem_ker_reciprocityMapLimit_iff` は「全層で固定する」という形
   だったので、位相を論じるために 1 層分に切り出した。
3. `preimage_principalUnits_eq_fixingSubgroup`(§3):
   `Φ⁻¹(1+π^{m+1}𝒪_K) = Gal(K_π/K(x_m))`。両辺とも「`x_m` を固定する」。
4. `continuous_lubinTateClosureGalEquivUnits` と
   **`lubinTateClosureGalContinuousMulEquivUnits : Gal(K_π/K) ≃ₜ* 𝒪_K^×`**。

## 証明の骨格

`Φ := lubinTateClosureGalEquivUnits` の連続性は `1` での連続性に帰着する
(`continuous_of_continuousAt_one`)。`𝒪_K^×` の `1` の近傍 `V` に対し
主単数群 `1+π^N𝒪_K ⊆ V` を取り(§1)、その `Φ` による引き戻しが
`Gal(K_π/K(x_N))`(§3)——`K(x_N)/K` は有限次なので**開**——であること
から `Φ⁻¹(V)` は `1` の近傍。あとは `Gal(K_π/K)` がコンパクト
(`K_π/K` は Galois)・`𝒪_K^×` が Hausdorff なので
`Continuous.homeoOfEquivCompactToT2` が逆写像の連続性を与える。

## 逸脱(記録)

* §3 の補題群は `[Normal K.carrier ↥K_π]` を**インスタンス変数として**
  受け取る形にしてある。`AlgEquiv.restrictNormalHom` が主張の**型**に
  正規性を要求するためで、`normal_lubinTateClosure`(在庫、`sorry` 無し)
  が実際にそれを供給する(§4 の仮説なしの形で消える)。`Normal` は
  `Prop` クラスなので、どのインスタンスを渡しても命題は同じ。
* `f` への依存は `lubinTateClosure` と同じ扱い(引数)のまま。新しい
  パラメータは増やしていない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical
open Filter Topology

/-! ## 1. `𝒪_K` の位相と `π` 進フィルター -/

section Integers

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)

/-- `x ∈ π^N 𝒪_K` なら `‖x‖ ≤ ‖π‖^N`。 -/
theorem norm_le_of_mem_span_pow {π : 𝒪[K.carrier]} (N : ℕ) (x : 𝒪[K.carrier])
    (hx : x ∈ Ideal.span ({π ^ N} : Set (𝒪[K.carrier]))) :
    ‖(x : K.carrier)‖ ≤ ‖(π : K.carrier)‖ ^ N := by
  obtain ⟨c, hc⟩ := Ideal.mem_span_singleton.mp hx
  have hcle : ‖(c : K.carrier)‖ ≤ 1 := by
    have h := c.2; rw [Valued.integer.mem_iff] at h; exact h
  have hxc : (x : K.carrier) = (π : K.carrier) ^ N * (c : K.carrier) := by
    rw [hc]; push_cast; ring
  rw [hxc, norm_mul, norm_pow]
  nlinarith [pow_nonneg (norm_nonneg (π : K.carrier)) N, norm_nonneg (c : K.carrier)]

/-- 逆向き: `‖x‖ ≤ ‖π‖^N` なら `x ∈ π^N 𝒪_K`(`x/π^N` のノルムが `1` 以下
なので `𝒪_K` の元)。 -/
theorem mem_span_pow_of_norm_le {π : 𝒪[K.carrier]} (hπne0 : π ≠ 0) (N : ℕ)
    (x : 𝒪[K.carrier]) (hx : ‖(x : K.carrier)‖ ≤ ‖(π : K.carrier)‖ ^ N) :
    x ∈ Ideal.span ({π ^ N} : Set (𝒪[K.carrier])) := by
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπne0 (Subtype.ext h)
  have hpos : (0 : ℝ) < ‖(π : K.carrier)‖ ^ N :=
    pow_pos (norm_pos_iff.mpr hπ0) N
  have hcmem : (x : K.carrier) / (π : K.carrier) ^ N ∈ 𝒪[K.carrier] := by
    rw [Valued.integer.mem_iff, norm_div, norm_pow]
    rw [div_le_one hpos]
    exact hx
  refine Ideal.mem_span_singleton.mpr ⟨⟨_, hcmem⟩, ?_⟩
  apply Subtype.ext
  push_cast
  field_simp

/-- 一意化元のノルムは `1` 未満。 -/
theorem norm_pi_lt_one {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π}) :
    ‖(π : K.carrier)‖ < 1 := by
  have hmem : π ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := by
    rw [hπmax]; exact Ideal.mem_span_singleton_self π
  have hnu : ¬ IsUnit π := fun h => (IsLocalRing.mem_maximalIdeal π).mp hmem h
  have hle : ‖(π : K.carrier)‖ ≤ 1 := by
    have h := π.2; rw [Valued.integer.mem_iff] at h; exact h
  exact lt_of_le_of_ne hle (fun h => hnu (Valued.integer.isUnit_iff_norm_eq_one.mpr h))

/-- ★**`π^n 𝒪_K` は `𝒪_K` の `0` の基本近傍系をなす**——`‖π‖ < 1` なので
`‖π‖^N` はいくらでも小さくでき、`π^N 𝒪_K` は半径 `‖π‖^N` の閉球に入る。 -/
theorem exists_span_pow_subset_of_mem_nhds {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    {U : Set (𝒪[K.carrier])} (hU : U ∈ 𝓝 1) :
    ∃ N : ℕ, {x : 𝒪[K.carrier] | x - 1 ∈ Ideal.span ({π ^ N} : Set (𝒪[K.carrier]))} ⊆ U := by
  obtain ⟨ε, hε, hball⟩ := Metric.mem_nhds_iff.mp hU
  obtain ⟨N, hN⟩ := exists_pow_lt_of_lt_one hε (norm_pi_lt_one K hπmax)
  refine ⟨N, fun x hx => hball ?_⟩
  have h1 : ‖((x - 1 : 𝒪[K.carrier]) : K.carrier)‖ ≤ ‖(π : K.carrier)‖ ^ N :=
    norm_le_of_mem_span_pow K N _ hx
  have h2 : dist x (1 : 𝒪[K.carrier]) = ‖((x - 1 : 𝒪[K.carrier]) : K.carrier)‖ := by
    show dist (x : K.carrier) ((1 : 𝒪[K.carrier]) : K.carrier) = _
    push_cast
    rw [dist_eq_norm]
  simp only [Metric.mem_ball, h2]
  linarith

/-- ★**`π^N 𝒪_K` は `𝒪_K` で開**——加法部分群であって `0` の近傍
(半径 `‖π‖^N` の開球)を含むから。 -/
theorem isOpen_span_pow {π : 𝒪[K.carrier]} (hπne0 : π ≠ 0) (N : ℕ) :
    IsOpen ((Ideal.span ({π ^ N} : Set (𝒪[K.carrier])) : Set (𝒪[K.carrier]))) := by
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hπne0 (Subtype.ext h)
  have hpos : (0 : ℝ) < ‖(π : K.carrier)‖ ^ N := pow_pos (norm_pos_iff.mpr hπ0) N
  have hsub : Metric.ball (0 : 𝒪[K.carrier]) (‖(π : K.carrier)‖ ^ N) ⊆
      (Ideal.span ({π ^ N} : Set (𝒪[K.carrier])) : Set (𝒪[K.carrier])) := by
    intro x hx
    refine mem_span_pow_of_norm_le K hπne0 N x (le_of_lt ?_)
    have hd : dist x (0 : 𝒪[K.carrier]) = ‖(x : K.carrier)‖ := by
      show dist (x : K.carrier) ((0 : 𝒪[K.carrier]) : K.carrier) = _
      push_cast
      rw [dist_eq_norm, sub_zero]
    simpa [hd] using hx
  exact (Ideal.span ({π ^ N} : Set (𝒪[K.carrier]))).toAddSubgroup.isOpen_of_mem_nhds
    (g := 0) (Filter.mem_of_superset (Metric.ball_mem_nhds 0 hpos) hsub)

/-- ★★**主単数群 `1 + π^N 𝒪_K` は `𝒪_K^×` で開**——`Units.val` は連続で、
`principalUnits` はその下での開集合 `1 + π^N 𝒪_K` の逆像。 -/
theorem isOpen_principalUnits {π : 𝒪[K.carrier]} (hπne0 : π ≠ 0) (N : ℕ) :
    IsOpen ((principalUnits K π N : Set (𝒪[K.carrier])ˣ)) := by
  have hopen : IsOpen {x : 𝒪[K.carrier] |
      x - 1 ∈ Ideal.span ({π ^ N} : Set (𝒪[K.carrier]))} := by
    have hcont : Continuous fun x : 𝒪[K.carrier] => x - 1 := continuous_id.sub continuous_const
    exact (isOpen_span_pow K hπne0 N).preimage hcont
  exact hopen.preimage Units.continuous_val

/-- ★★**主単数群は `𝒪_K^×` の `1` の基本近傍系をなす**——単数群の位相は
`Units.embedProduct`(`u ↦ (u, u⁻¹)`)による誘導位相なので、`u` と `u⁻¹`
の両方を `1 + π^N 𝒪_K` に押し込めればよい。`principalUnits` は部分群なので
`u ∈ principalUnits` から `u⁻¹ ∈ principalUnits` が出る。 -/
theorem exists_principalUnits_subset_of_mem_nhds {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    {V : Set (𝒪[K.carrier])ˣ} (hV : V ∈ 𝓝 1) :
    ∃ N : ℕ, (principalUnits K π N : Set (𝒪[K.carrier])ˣ) ⊆ V := by
  rw [(Units.isEmbedding_embedProduct (M := 𝒪[K.carrier])).toIsInducing.nhds_eq_comap,
    show Units.embedProduct (𝒪[K.carrier]) 1 = ((1 : 𝒪[K.carrier]), (1 : (𝒪[K.carrier])ᵐᵒᵖ)) by
      simp,
    Filter.mem_comap] at hV
  obtain ⟨W, hW, hWV⟩ := hV
  obtain ⟨U₁, hU₁, U₂, hU₂, hprod⟩ := mem_nhds_prod_iff.mp hW
  have hU₂' : MulOpposite.op ⁻¹' U₂ ∈ 𝓝 (1 : 𝒪[K.carrier]) := by
    have hc : Continuous (MulOpposite.op : 𝒪[K.carrier] → (𝒪[K.carrier])ᵐᵒᵖ) :=
      MulOpposite.continuous_op
    exact hc.continuousAt (x := (1 : 𝒪[K.carrier])) (by simpa using hU₂)
  obtain ⟨N, hN⟩ := exists_span_pow_subset_of_mem_nhds K hπmax (Filter.inter_mem hU₁ hU₂')
  refine ⟨N, fun u hu => hWV ?_⟩
  have h1 : (u : 𝒪[K.carrier]) - 1 ∈ Ideal.span ({π ^ N} : Set (𝒪[K.carrier])) := hu
  have h2 : ((u⁻¹ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) - 1 ∈
      Ideal.span ({π ^ N} : Set (𝒪[K.carrier])) := (principalUnits K π N).inv_mem hu
  refine hprod ?_
  rw [Units.embedProduct_apply]
  exact ⟨(hN h1).1, (hN h2).2⟩

end Integers

/-! ## 2. `reciprocityMapLimit` の層ごとの核 -/

section Recip

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

/-- ★**層ごとの核**: `reciprocityMapLimit σ` の第 `m+1` 成分が `1` である
ことと、`σ` が level `m+1` の生成元 `x_m` を固定することは同値。
`Found/PGC/LubinTateReciprocityMapLimitKernel.lean` の
`mem_ker_reciprocityMapLimit_iff` は「全 `m` で固定する」という形で
述べられていたが、位相を論じるには**各層ごと**の形が要る——証明は
同ファイルの各方向をそのまま 1 層に切り出したもの。 -/
theorem reciprocityMapLimitFamily_succ_eq_one_iff
    (σ : K.closure ≃ₐ[K.carrier] K.closure) (m : ℕ) :
    reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1) = 1 ↔
      σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt := by
  haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  constructor
  · intro h1
    change principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) = 1 at h1
    rw [MulEquiv.map_eq_one_iff] at h1
    have hspec := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
    rw [h1, ← QuotientGroup.mk_one, unitActionQuotientLift_mk] at hspec
    rw [← hspec, show ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) = 1 from rfl,
      lubinTateActionAtTorsionPoint_one]
  · intro hσm
    show principalUnitsQuotientEquiv K hπmax (m + 1) (by omega)
        (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ) = 1
    have h1 : σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
        (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1)
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
          (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])) :
          IntermediateField.adjoin K.carrier
            ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) : K.closure) := by
      rw [show ((1 : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) = 1 from rfl,
        lubinTateActionAtTorsionPoint_one]
      exact hσm
    rw [reciprocityMap_eq_mk_of_apply_eq K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) _
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ _ _ 1 σ h1,
      QuotientGroup.mk_one]
    exact map_one _

end Recip

/-! ## 3. `Gal(K_π/K) ≃ₜ* 𝒪_K^×` -/

/-- `φ` が全射のとき、`(G ⧸ ker φ ≃* H)` の逆写像は `φ g` を `mk g` に戻す。 -/
theorem quotientKerEquivOfSurjective_symm_apply {G H : Type*} [Group G] [Group H] (φ : G →* H)
    (hφ : Function.Surjective φ) (g : G) :
    (QuotientGroup.quotientKerEquivOfSurjective φ hφ).symm (φ g) = QuotientGroup.mk g := by
  rw [MulEquiv.symm_apply_eq]
  rfl

section Main

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

/-- `Gal(K̄/K) →* 𝒪_K^×`(相互律の単数群版)——`reciprocityMapLimit` の後に
`unitsEquivCompatibleUnits.symm` を合成しただけ。 -/
noncomputable def reciprocityUnits :
    (K.closure ≃ₐ[K.carrier] K.closure) →* (𝒪[K.carrier])ˣ :=
  ((unitsEquivCompatibleUnits K hπmax).symm :
      CompatibleUnits K hπmax ≃* (𝒪[K.carrier])ˣ).toMonoidHom.comp
    (reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf)

/-- ★`reciprocityUnits σ` が level `m+1` の主単数であることと、`σ` が
生成元 `x_m` を固定することは同値。 -/
theorem mem_principalUnits_reciprocityUnits_iff
    (σ : K.closure ≃ₐ[K.carrier] K.closure) (m : ℕ) :
    reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ ∈ principalUnits K π (m + 1) ↔
      σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt := by
  rw [principalUnits_eq_ker, MonoidHom.mem_ker]
  have h : unitsEquivCompatibleUnits K hπmax
      (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ)
      = reciprocityMapLimit K hq hπmax hπne0 f hf0 hf1 hf σ := MulEquiv.apply_symm_apply _ _
  have hkey : (Units.map ((Ideal.Quotient.mk
        (Ideal.span ({π ^ (m + 1)} : Set (𝒪[K.carrier])))).toMonoidHom))
      (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ)
      = reciprocityMapLimitFamily K hq hπmax hπne0 f hf0 hf1 hf σ (m + 1) :=
    congrFun (congrArg Subtype.val h) (m + 1)
  rw [hkey]
  exact reciprocityMapLimitFamily_succ_eq_one_iff K hq hπmax hπne0 f hf0 hf1 hf σ m

/-- 生成元 `x_m` を `K_π` の元とみたもの。 -/
noncomputable def lubinTateClosureGen (m : ℕ) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf :=
  ⟨(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt, by
    rw [lubinTateClosure_eq_adjoin_genSet]
    exact IntermediateField.subset_adjoin K.carrier _ ⟨m, rfl⟩⟩

/-- `K_π` は `K` 上整——`K̄` が `K` 上代数的だから。 -/
theorem isIntegral_lubinTateClosureGen (m : ℕ) :
    IsIntegral K.carrier (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m) := by
  have hclos : IsIntegral K.carrier
      ((lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m : K.closure)) :=
    (Algebra.IsAlgebraic.isAlgebraic
      (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m : K.closure)).isIntegral
  exact (isIntegral_algebraMap_iff
    (algebraMap (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) K.closure).injective).mp hclos

/-- level `m` の中間体 `K(x_m) ⊆ K_π`。 -/
noncomputable def lubinTateClosureLevel (m : ℕ) :
    IntermediateField K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :=
  IntermediateField.adjoin K.carrier
    ({lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m} :
      Set (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf))

instance finiteDimensional_lubinTateClosureLevel (m : ℕ) :
    FiniteDimensional K.carrier (lubinTateClosureLevel K hq hπmax hπne0 f hf0 hf1 hf m) :=
  IntermediateField.adjoin.finiteDimensional
    (isIntegral_lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m)

variable [Normal K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)]

/-- ★`Φ ∘ (K̄ への制限) = reciprocityUnits`——`lubinTateClosureGalEquivUnits`
の定義(2 回の第一同型定理)を実際に 1 段ずつほどく。 -/
theorem lubinTateClosureGalEquivUnits_restrictNormalHom
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf
        (AlgEquiv.restrictNormalHom
          (E := (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)) σ)
      = reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ := by
  simp only [lubinTateClosureGalEquivUnits, MulEquiv.trans_apply,
    quotientKerEquivOfSurjective_symm_apply, QuotientGroup.quotientMulEquivOfEq_mk]
  rfl

/-- ★★**`Φ⁻¹(1 + π^{m+1}𝒪_K) = Gal(K_π/K(x_m))`**——両辺とも「`x_m` を
固定する」で言い換えられる。 -/
theorem preimage_principalUnits_eq_fixingSubgroup (m : ℕ) :
    {τ : lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
        lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf |
      lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf τ ∈
        principalUnits K π (m + 1)} =
      ((lubinTateClosureLevel K hq hπmax hπne0 f hf0 hf1 hf m).fixingSubgroup :
        Set (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
          lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)) := by
  ext τ
  obtain ⟨σ, rfl⟩ := AlgEquiv.restrictNormalHom_surjective
    (F := K.carrier)
    (K₁ := (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)) K.closure τ
  have hcom := AlgEquiv.restrictNormal_commutes σ
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)
    (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m)
  have hfix : (AlgEquiv.restrictNormalHom
        (E := (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)) σ)
        (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m) =
      lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m ↔
      σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) =
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt := by
    rw [Subtype.ext_iff]
    show ((σ.restrictNormal (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)
        (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m) : K.closure) = _) ↔ _
    rw [show ((σ.restrictNormal (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf : Type _)
        (lubinTateClosureGen K hq hπmax hπne0 f hf0 hf1 hf m) : K.closure))
      = σ ((psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt) from hcom]
    exact Iff.rfl
  simp only [Set.mem_setOf_eq, SetLike.mem_coe,
    lubinTateClosureGalEquivUnits_restrictNormalHom,
    mem_principalUnits_reciprocityUnits_iff, lubinTateClosureLevel,
    mem_fixingSubgroup_adjoin_iff, Set.mem_singleton_iff, forall_eq]
  exact hfix.symm


/-- `K_π/K` は Galois——正規(`normal_lubinTateClosure`)かつ、標数 `0` なので
分離的。 -/
theorem isGalois_lubinTateClosure :
    IsGalois K.carrier (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) := by
  haveI := charZero_carrier K
  haveI : PerfectField K.carrier := PerfectField.ofCharZero
  haveI : Algebra.IsSeparable K.carrier
      (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :=
    Algebra.IsAlgebraic.isSeparable_of_perfectField
  exact ⟨⟩

/-- ★★★**`Gal(K_π/K) → 𝒪_K^×` は連続**——`1` での連続性だけ見ればよい
(`continuous_of_continuousAt_one`)。`𝒪_K^×` の `1` の近傍 `V` を取ると、
主単数群が基本近傍系をなす(`exists_principalUnits_subset_of_mem_nhds`)
ので `1 + π^N𝒪_K ⊆ V`、その引き戻しは `Gal(K_π/K(x_N))`
(`preimage_principalUnits_eq_fixingSubgroup`)で、`K(x_N)/K` は有限次
なので開(`IntermediateField.fixingSubgroup_isOpen`)。 -/
theorem continuous_lubinTateClosureGalEquivUnits :
    Continuous (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf) := by
  apply continuous_of_continuousAt_one
  simp only [ContinuousAt, map_one, Filter.tendsto_def]
  intro V hV
  obtain ⟨N, hN⟩ := exists_principalUnits_subset_of_mem_nhds K hπmax hV
  have hsub : ((lubinTateClosureLevel K hq hπmax hπne0 f hf0 hf1 hf N).fixingSubgroup :
        Set (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
          lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)) ⊆
      (lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf) ⁻¹' V := by
    rw [← preimage_principalUnits_eq_fixingSubgroup K hq hπmax hπne0 f hf0 hf1 hf N]
    exact fun τ hτ => hN (principalUnits_succ_le K π N hτ)
  exact Filter.mem_of_superset
    (IsOpen.mem_nhds (IntermediateField.fixingSubgroup_isOpen _)
      (lubinTateClosureLevel K hq hπmax hπne0 f hf0 hf1 hf N).fixingSubgroup.one_mem) hsub

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Gal(K_π/K) ≃ₜ* 𝒪_K^×`(位相群としての同型)**。

`lubinTateClosureGalEquivUnits`(群同型)が連続であることと、
`Gal(K_π/K)` がコンパクト(`K_π/K` が Galois、mathlib の `InfiniteGalois`)・
`𝒪_K^×` が Hausdorff であることから、`Continuous.homeoOfEquivCompactToT2`
で逆写像の連続性が自動的に出る。 -/
noncomputable def lubinTateClosureGalContinuousMulEquivUnits :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ≃ₐ[K.carrier]
      lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) ≃ₜ* (𝒪[K.carrier])ˣ := by
  haveI := isGalois_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  have hcont : Continuous
      ⇑(lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf).toEquiv :=
    continuous_lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf
  exact ⟨lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf,
    hcont.homeoOfEquivCompactToT2.continuous_toFun,
    hcont.homeoOfEquivCompactToT2.continuous_invFun⟩

end Main

/-! ## 4. 仮説なしの形 -/

/-- ★★★★**仮説なしの Λ3(位相版)**: 任意の p 進局所体 `K` について、
`K̄` の中に `K` 上正規な中間体 `K_π` が在って
`Gal(K_π/K) ≃ₜ* 𝒪_K^×`(**位相群としての同型**)。

`Found/PGC/LubinTateClosure.lean::exists_lubinTateClosure_galEquiv_units`
(群同型版)の位相版。仮説の組み上げ方はそちらと同一。 -/
theorem exists_lubinTateClosure_continuousMulEquiv_units {p : ℕ} [Fact p.Prime]
    (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧
        Nonempty ((E ≃ₐ[K.carrier] E) ≃ₜ* (𝒪[K.carrier])ˣ) := by
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
  haveI := normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  exact ⟨lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    ⟨lubinTateClosureGalContinuousMulEquivUnits K hq hπmax hπne0 f hf0 hf1 hf⟩⟩


end ABC3.Found.PGC
