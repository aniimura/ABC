import ABC3.Found.PGC.UnramifiedResidueField
import ABC3.Found.PGC.TotallyRamified

/-!
# `𝒪_{K̂^{ur}}` の極大イデアルは `π` で生成される(`sorry` 無し)

経路 Λ の節点 **Λ6 の入口 M1**。Λ6(Dwork の補題 = `Art_π` の π 非依存性)は
`𝒪_{K̂^{ur}}` 上の逐次近似で、**各段で「`π` で割る」**ができなければ収束を
管理できない。本ファイルはその土台を置く:

```
𝔪_{𝒪_{K̂^ur}} = (π)      ‖K̂^{ur}‖ = ‖K‖      𝒪_{K̂^ur} は離散付値環
```

`Found/PGC/UnramifiedCompletion.lean`(Λ5)は `K̂^{ur}` を
**ノルム(スペクトルノルム)の完備化**として作り、「`𝒪_{K^{ur}}` の `𝔪`-進完備化の
商体と同じものだが一致は未証明」と逸脱を記録した。本ファイルはその穴のうち
**下流(Λ6)が実際に使う部分**——値群が拡大でも完備化でも増えないこと——を埋める。

## 到達点

| | 主張 |
|---|---|
| 1 | `norm_le_norm_uniformizer_of_isUnramifiedAdjoin`:`K(x)/K` が不分岐なら `K(x)` の値に `(‖π‖, 1)` の隙間は無い |
| 2 | `norm_le_norm_uniformizer_unramifiedCompletion`:同じことが `K̂^{ur}` でも成り立つ(稠密性) |
| 3 | **`maximalIdeal_unramifiedCompletionInt_eq_span`**:`𝔪_{𝒪_{K̂^ur}} = (π)` |
| 4 | `exists_norm_eq_zpow_unramifiedCompletion` / `norm_range_unramifiedCompletion`:**`‖K̂^{ur}‖ = ‖K‖`** |
| 5 | `isPrincipalIdealRing_unramifiedCompletionInt` / `isDiscreteValuationRing_unramifiedCompletionInt`:**`𝒪_{K̂^{ur}}` は DVR** |
| 6 | `exists_eq_uniformizer_mul`:`‖z‖ < 1` なら `z = π·w`(Λ6 が使う「π で割る」) |

## 証明の骨(段取り測定の見立てとの差)

段取り測定は「**稠密性＋値群の離散性**で済むはず」と見立てた。**半分当たっている**。

* 稠密性(2)はその通り。`exists_norm_sub_lt` で `‖z - x‖ < ‖z‖` を取れば
  超距離から `‖x‖ = ‖z‖` になり、値は `K^{ur}` 側の値に一致する。
* しかし「**有限段では `e = 1` なので値群が `‖K‖` と同じ**」という部分は
  在庫に**無かった**。`Found/PGC/UnramifiedExtension.lean` の
  `ramificationIndex`(= `Ideal.ramificationIdx`)とノルムを繋ぐ補題が
  片方向(`TotallyRamified.lean` の「`‖x‖^n = ‖π‖` ⟹ 完全分岐」)しか無い。
* そこで(1)は**逆向きに**組んだ。`‖π‖ < ‖y‖ < 1` を仮定すると
  `π = (π/y)·y` で `π/y` も `y` も `𝔪_L` に入るから `π ∈ 𝔪_L²`、すなわち
  `map 𝔪_K ≤ 𝔪_L²`、よって `2 ≤ e`(`le_ramificationIdx_of_map_le_pow`)。
  不分岐 `e = 1` に矛盾。**DVR の一意化元 `ϖ` を経由しない**のが要点で、
  「`𝒪_L` が DVR」も「`(π) = 𝔪_L^e`」も要らない。
* 「値群の離散性」は(1)(2)から**出てくる**もので、前提として使ってはいない。
  (4)は `exists_mem_Ioc_zpow` で `‖π‖^{n+1} < ‖z‖ ≤ ‖π‖^n` を取り、
  `z/π^n` に(2)を当てて `‖z/π^n‖ = 1` を得る形。

## ★設計上の注意(守ったこと)

* **`inertia` を経由していない**。`Interface` の `SubgroupCorrespondence` /
  `ResidueCardinality` は本ファイルの主張にも証明にも現れない
  (Corollary 1.3 との循環を避けるため)。不分岐側は `unramifiedClosure K` と
  `IsUnramifiedAdjoin` で直接扱っている。
* **`Valued` を `K^{ur}` の instance にしていない**。Λ5 の設計判断をそのまま
  引き継ぎ、`K̂^{ur}` 側の `NormedField.toValued` だけを使う。
* **自由なパラメータを結論に出していない**。`π` は `hπmax`(極大イデアルの
  生成元)で釘づけされており、非空虚性は
  `exists_maximalIdeal_unramifiedCompletionInt_eq_span`(`π` を `∃` の内側に
  閉じ込めた形)で保証してある。(4)(5)は `π` を含まない。
* **`sorry` 本体の `def` を作っていない**(D21)。
* `Abelianization` を使っていない。

## 逸脱(記録)

* Λ5 の逸脱(「ノルムの完備化 vs `𝔪`-進完備化」)は**まだ完全には閉じていない**。
  本ファイルが閉じたのは「値群が一致する」「整数環が DVR」という
  **下流が使う部分**であって、「`𝒪_{K̂^{ur}}` が `𝒪_{K^{ur}}` の `𝔪`-進完備化と
  環同型」という主張そのものは証明していない。必要になったら
  (4)(5)+ `Found/PGC/UnramifiedCompletion.lean` の `residueFieldEquiv` から
  組めるはずなので、そのときに節点を立てること。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued NNReal

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. 有限段——不分岐なら値に `(‖π‖, 1)` の隙間は無い -/

/-- **★★★★★★★★不分岐拡大の値群は `K` の値群を超えない**(隙間の形)。

`K(x)/K` が不分岐で `y ∈ K(x)` が `‖y‖ < 1` なら `‖y‖ ≤ ‖π‖`。

証明は背理法。`‖π‖ < ‖y‖ < 1` とすると `y` も `π/y` も `𝔪_L` に入り、
`π = (π/y)·y ∈ 𝔪_L²`。ゆえに `map 𝔪_K ≤ 𝔪_L²` で `2 ≤ e`
(`le_ramificationIdx_of_map_le_pow`)、不分岐 `e = 1` と矛盾する。

★DVR の一意化元を経由しない——`𝒪_L` が DVR であることも
`(π) = 𝔪_L^e` も使わない。

退化の自己検査:`hu`(不分岐)を落とすと**偽**。Eisenstein の根 `α`
(`‖α‖^n = ‖π‖`、`n ≥ 2`)は `‖π‖ < ‖α‖ < 1` を満たす。 -/
theorem norm_le_norm_uniformizer_of_isUnramifiedAdjoin (K : PAdicLocalField p) (x : K.closure)
    (hu : IsUnramifiedAdjoin K x)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (y : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (hy : ‖y‖ < 1) :
    ‖y‖ ≤ ‖(π : K.carrier)‖ := by
  haveI := isLocalRing_adjoinIntegers K x
  by_contra hcon
  rw [not_le] at hcon
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hypos : 0 < ‖y‖ := lt_trans hπpos hcon
  have hy0 : y ≠ 0 := norm_pos_iff.mp hypos
  have hπL : ‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier))‖
      = ‖(π : K.carrier)‖ := norm_algebraMap' _ _
  have hcnorm : ‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier)) / y‖ < 1 := by
    rw [norm_div, hπL, div_lt_one hypos]
    exact hcon
  have hcle : ‖(algebraMap K.carrier
      (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) (π : K.carrier)) / y‖ ≤ 1 :=
    le_of_lt hcnorm
  have hamem : (⟨y, le_of_lt hy⟩ : adjoinIntegers K x)
      ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) :=
    (mem_maximalIdeal_adjoinIntegers K x _).mpr hy
  have hbmem : (⟨_, hcle⟩ : adjoinIntegers K x)
      ∈ IsLocalRing.maximalIdeal (adjoinIntegers K x) :=
    (mem_maximalIdeal_adjoinIntegers K x _).mpr hcnorm
  have heq : (algebraMap 𝒪[K.carrier] (adjoinIntegers K x)) π
      = (⟨_, hcle⟩ : adjoinIntegers K x) * (⟨y, le_of_lt hy⟩ : adjoinIntegers K x) := by
    apply Subtype.ext
    show algebraMap K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      (π : K.carrier) = _
    push_cast
    rw [div_mul_cancel₀ _ hy0]
  have hmap : Ideal.map (algebraMap 𝒪[K.carrier] (adjoinIntegers K x))
      (IsLocalRing.maximalIdeal 𝒪[K.carrier])
      ≤ (IsLocalRing.maximalIdeal (adjoinIntegers K x)) ^ 2 := by
    rw [hπmax, Ideal.map_span, Ideal.span_le, Set.image_singleton, Set.singleton_subset_iff,
      SetLike.mem_coe, heq, sq]
    exact Ideal.mul_mem_mul hbmem hamem
  have h2 : 2 ≤ ramificationIndex K x :=
    le_ramificationIdx_of_map_le_pow (ramificationIndex_ne_zero K x) hmap
  rw [show ramificationIndex K x = 1 from hu] at h2
  omega

/-- `K^{ur}` の元は `K(z)/K` が不分岐なので上の補題がそのまま効く。 -/
theorem norm_le_norm_uniformizer_of_mem_unramifiedClosure (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    {z : K.closure} (hz : z ∈ unramifiedClosure K) (hlt : ‖z‖ < 1) :
    ‖z‖ ≤ ‖(π : K.carrier)‖ :=
  norm_le_norm_uniformizer_of_isUnramifiedAdjoin K z
    ((mem_unramifiedClosure_iff_isUnramified K z).mp hz) hπ0 hπmax
    ⟨z, IntermediateField.mem_adjoin_simple_self K.carrier z⟩ hlt

/-! ## 2. 完備化へ持ち上げる(稠密性) -/

/-- **★★★★★★★★★★`K̂^{ur}` の値にも `(‖π‖, 1)` の隙間は無い**。

`z ≠ 0` なら `‖z - x‖ < ‖z‖` なる `x ∈ K^{ur}` が取れ(`exists_norm_sub_lt`)、
超距離から `‖x‖ = ‖z‖`。あとは有限段の補題。

退化の自己検査:`z = 0` は別扱い(`0 ≤ ‖π‖`)。`hlt` を落とすと偽
(`‖z‖ = 1` の元はいくらでもある)。 -/
theorem norm_le_norm_uniformizer_unramifiedCompletion (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (z : unramifiedCompletion K) (hlt : ‖z‖ < 1) :
    ‖z‖ ≤ ‖(π : K.carrier)‖ := by
  rcases eq_or_ne z 0 with rfl | hz0
  · simp
  have hzpos : 0 < ‖z‖ := norm_pos_iff.mpr hz0
  obtain ⟨x, hx⟩ := exists_norm_sub_lt K z hzpos
  have hxle : ‖(x : unramifiedCompletion K)‖ ≤ ‖z‖ := by
    have h := IsUltrametricDist.norm_add_le_max z (-(z - (x : unramifiedCompletion K)))
    simp only [neg_sub, add_sub_cancel] at h
    refine h.trans (max_le le_rfl ?_)
    rw [norm_sub_rev]
    exact hx.le
  have hxge : ‖z‖ ≤ ‖(x : unramifiedCompletion K)‖ := by
    have h := IsUltrametricDist.norm_add_le_max (x : unramifiedCompletion K)
      (z - (x : unramifiedCompletion K))
    rw [add_sub_cancel] at h
    rcases max_cases ‖(x : unramifiedCompletion K)‖ ‖z - (x : unramifiedCompletion K)‖ with
      ⟨he, -⟩ | ⟨he, -⟩
    · rwa [he] at h
    · rw [he] at h; linarith
  have hxeq : ‖x‖ = ‖z‖ := by
    rw [← norm_coe_unramifiedCompletion]
    exact le_antisymm hxle hxge
  rw [← hxeq]
  exact norm_le_norm_uniformizer_of_mem_unramifiedClosure K hπ0 hπmax x.2 (hxeq ▸ hlt)

/-! ## 3. `𝔪_{𝒪_{K̂^ur}} = (π)` -/

/-- 極大イデアルの生成元はノルムが `1` 未満。 -/
theorem norm_lt_one_of_maximalIdeal_eq_span (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    ‖(π : K.carrier)‖ < 1 := by
  have hπle : ‖(π : K.carrier)‖ ≤ 1 := by
    have h := π.2; rw [Valued.integer.mem_iff] at h; exact h
  have hπmem : π ∈ IsLocalRing.maximalIdeal 𝒪[K.carrier] := by
    rw [hπmax]; exact Ideal.mem_span_singleton_self π
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff,
    Valued.integer.isUnit_iff_norm_eq_one] at hπmem
  exact lt_of_le_of_ne hπle hπmem

/-- `π ∈ 𝒪_K` の `𝒪_{K̂^{ur}}` における像。 -/
noncomputable def uniformizerCompletionInt (K : PAdicLocalField p) (π : 𝒪[K.carrier]) :
    ↥(unramifiedCompletionInt K) :=
  ⟨algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier), by
    rw [mem_unramifiedCompletionInt, norm_algebraMap_unramifiedCompletion]
    have h := π.2; rw [Valued.integer.mem_iff] at h; exact h⟩

@[simp] theorem uniformizerCompletionInt_coe (K : PAdicLocalField p) (π : 𝒪[K.carrier]) :
    ((uniformizerCompletionInt K π : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      = algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier) := rfl

/-- **★★★★★★★★★★★★★★★★(Λ6-M1)`𝒪_{K̂^{ur}}` の極大イデアルは `π` で生成される**。

`e(K^{ur}/K) = 1` の完備化版。`⊇` は `‖π‖ < 1` から、`⊆` は
「`‖w‖ < 1` なら `‖w‖ ≤ ‖π‖`」(§2)から `w/π ∈ 𝒪_{K̂^{ur}}`。

★Λ6(Dwork の補題)の逐次近似が各段で「`π` で割る」ために使う
(`exists_eq_uniformizer_mul`)。

退化の自己検査。

* `hπmax` を落とすと `π` が自由なパラメータになり主張は偽
  (`π` を `π²` に替えれば右辺は小さくなる)。
* `hπ0` を落とすと `π = 0` で右辺が `⊥`、左辺は `⊥` でないので偽
  (`isDiscreteValuationRing_unramifiedCompletionInt` の証明を参照)。
* 非空虚性は `exists_maximalIdeal_unramifiedCompletionInt_eq_span` で保証。 -/
theorem maximalIdeal_unramifiedCompletionInt_eq_span (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
      = Ideal.span {uniformizerCompletionInt K π} := by
  have hπlt := norm_lt_one_of_maximalIdeal_eq_span K hπmax
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hπC : algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier) ≠ 0 := by
    rw [← norm_ne_zero_iff, norm_algebraMap_unramifiedCompletion]
    exact ne_of_gt hπpos
  apply le_antisymm
  · intro w hw
    rw [mem_maximalIdeal_unramifiedCompletionInt] at hw
    have hle := norm_le_norm_uniformizer_unramifiedCompletion K hπ0 hπmax
      (w : unramifiedCompletion K) hw
    rw [Ideal.mem_span_singleton']
    have hc : ‖(w : unramifiedCompletion K)
        / algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)‖ ≤ 1 := by
      rw [norm_div, norm_algebraMap_unramifiedCompletion, div_le_one hπpos]
      exact hle
    refine ⟨⟨_, (mem_unramifiedCompletionInt K _).mpr hc⟩, ?_⟩
    apply Subtype.ext
    show (w : unramifiedCompletion K)
        / algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)
        * algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)
      = (w : unramifiedCompletion K)
    exact div_mul_cancel₀ _ hπC
  · rw [Ideal.span_le, Set.singleton_subset_iff, SetLike.mem_coe,
      mem_maximalIdeal_unramifiedCompletionInt, uniformizerCompletionInt_coe,
      norm_algebraMap_unramifiedCompletion]
    exact hπlt

/-- `𝒪_K` の一意化元(`𝒪_K` は DVR)。非空虚性の材料。 -/
theorem exists_uniformizer_carrierIntegers (K : PAdicLocalField p) :
    ∃ π : 𝒪[K.carrier], (π : K.carrier) ≠ 0 ∧
      IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π} := by
  haveI := isDiscreteValuationRing_carrierIntegers K
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  exact ⟨ϖ, fun h => hϖ.ne_zero (Subtype.ext h), hϖ.maximalIdeal_eq⟩

/-- **★★★★★★★★★★★★(Λ6-M1)極大イデアルは `π` で生成される**(`π` を `∃` の
内側に閉じ込めた非空虚形)。★結論に自由なパラメータが出ない。 -/
theorem exists_maximalIdeal_unramifiedCompletionInt_eq_span (K : PAdicLocalField p) :
    ∃ π : 𝒪[K.carrier], (π : K.carrier) ≠ 0 ∧
      IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π} ∧
      IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
        = Ideal.span {uniformizerCompletionInt K π} := by
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  exact ⟨π, hπ0, hπmax, maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax⟩

/-- **`π` で割る**——`‖z‖ < 1` なら `z = π·w`。★Λ6 の逐次近似が各段で使う形。 -/
theorem exists_eq_uniformizer_mul (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    {z : ↥(unramifiedCompletionInt K)} (hz : ‖(z : unramifiedCompletion K)‖ < 1) :
    ∃ w : ↥(unramifiedCompletionInt K), z = uniformizerCompletionInt K π * w := by
  have hmem : z ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) :=
    (mem_maximalIdeal_unramifiedCompletionInt K z).mpr hz
  rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax,
    Ideal.mem_span_singleton'] at hmem
  obtain ⟨w, hw⟩ := hmem
  exact ⟨w, by rw [← hw, mul_comm]⟩

/-! ## 4. `‖K̂^{ur}‖ = ‖K‖` -/

/-- **★★★★★★★★★★`K̂^{ur}` の値は `‖π‖` の整数冪**。

`exists_mem_Ioc_zpow` で `‖π‖^{n+1} < ‖z‖ ≤ ‖π‖^n` を取り、`z/π^n` に §2 を
当てる:`‖z/π^n‖ < 1` なら `‖z‖ ≤ ‖π‖^{n+1}` で矛盾、ゆえに `‖z/π^n‖ = 1`。

★「値群の離散性」は前提ではなく**ここで出てくる**もの。 -/
theorem exists_norm_eq_zpow_unramifiedCompletion (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    {z : unramifiedCompletion K} (hz : z ≠ 0) :
    ∃ n : ℤ, ‖z‖ = ‖(π : K.carrier)‖ ^ n := by
  have hπlt := norm_lt_one_of_maximalIdeal_eq_span K hπmax
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hzpos : 0 < ‖z‖ := norm_pos_iff.mpr hz
  obtain ⟨n, hn⟩ := exists_mem_Ioc_zpow hzpos ((one_lt_inv₀ hπpos).mpr hπlt)
  rw [Set.mem_Ioc, inv_zpow, ← zpow_neg, inv_zpow, ← zpow_neg] at hn
  obtain ⟨hn1, hn2⟩ := hn
  refine ⟨-(n + 1), ?_⟩
  have hmpos : 0 < ‖(π : K.carrier)‖ ^ (-(n + 1)) := zpow_pos hπpos _
  have hPnorm : ‖algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier)‖
      = ‖(π : K.carrier)‖ := norm_algebraMap_unramifiedCompletion K _
  have hwnorm : ‖z / (algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier))
      ^ (-(n + 1))‖ = ‖z‖ / ‖(π : K.carrier)‖ ^ (-(n + 1)) := by
    rw [norm_div, norm_zpow, hPnorm]
  have hle1 : ‖z / (algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier))
      ^ (-(n + 1))‖ ≤ 1 := by
    rw [hwnorm, div_le_one hmpos]; exact hn2
  rcases lt_or_eq_of_le hle1 with hlt1 | heq1
  · exfalso
    have hkey := norm_le_norm_uniformizer_unramifiedCompletion K hπ0 hπmax _ hlt1
    rw [hwnorm, div_le_iff₀ hmpos] at hkey
    have hexp : ‖(π : K.carrier)‖ ^ (-(n + 1)) * ‖(π : K.carrier)‖
        = ‖(π : K.carrier)‖ ^ (-n) := by
      rw [← zpow_add_one₀ (ne_of_gt hπpos)]; congr 1; ring
    rw [mul_comm, hexp] at hkey
    linarith
  · rw [hwnorm, div_eq_one_iff_eq (ne_of_gt hmpos)] at heq1
    exact heq1

/-- `K̂^{ur}` の各元のノルムは `K` のある元のノルムに一致する。 -/
theorem exists_norm_eq_norm_carrier (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    ∃ a : K.carrier, ‖z‖ = ‖a‖ := by
  rcases eq_or_ne z 0 with rfl | hz
  · exact ⟨0, by simp⟩
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  obtain ⟨n, hn⟩ := exists_norm_eq_zpow_unramifiedCompletion K hπ0 hπmax hz
  exact ⟨(π : K.carrier) ^ n, by rw [norm_zpow]; exact hn⟩

/-- **★★★★★★★★★★★★(Λ6-M1)`‖K̂^{ur}‖ = ‖K‖`**——値群は不分岐拡大でも
完備化でも増えない。★結論に自由なパラメータが出ない形。 -/
theorem norm_range_unramifiedCompletion (K : PAdicLocalField p) :
    (Set.range fun z : unramifiedCompletion K => ‖z‖)
      = (Set.range fun a : K.carrier => ‖a‖) := by
  apply Set.eq_of_subset_of_subset
  · rintro - ⟨z, rfl⟩
    obtain ⟨a, ha⟩ := exists_norm_eq_norm_carrier K z
    exact ⟨a, ha.symm⟩
  · rintro - ⟨a, rfl⟩
    exact ⟨algebraMap K.carrier (unramifiedCompletion K) a,
      norm_algebraMap_unramifiedCompletion K a⟩

/-! ## 5. `𝒪_{K̂^{ur}}` は離散付値環 -/

/-- `𝒪_{K̂^{ur}}` の `0` でない元のノルムは `‖π‖` の**自然数**冪。 -/
theorem exists_norm_eq_pow_unramifiedCompletionInt (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπ0 : (π : K.carrier) ≠ 0)
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    {w : ↥(unramifiedCompletionInt K)} (hw : w ≠ 0) :
    ∃ n : ℕ, ‖(w : unramifiedCompletion K)‖ = ‖(π : K.carrier)‖ ^ n := by
  have hπlt := norm_lt_one_of_maximalIdeal_eq_span K hπmax
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hw0 : (w : unramifiedCompletion K) ≠ 0 := fun h => hw (Subtype.ext h)
  obtain ⟨m, hm⟩ := exists_norm_eq_zpow_unramifiedCompletion K hπ0 hπmax hw0
  have hwle : ‖(w : unramifiedCompletion K)‖ ≤ 1 := (mem_unramifiedCompletionInt K _).mp w.2
  have hm0 : 0 ≤ m := by
    have h1 : (1 : ℝ) < ‖(π : K.carrier)‖⁻¹ := (one_lt_inv₀ hπpos).mpr hπlt
    have h2 : ‖(π : K.carrier)‖⁻¹ ^ (-m) ≤ 1 := by
      rw [inv_zpow, zpow_neg, inv_inv, ← hm]; exact hwle
    have h3 := (zpow_le_one_iff_right₀ h1).mp h2
    omega
  refine ⟨m.toNat, ?_⟩
  rw [hm, ← zpow_natCast, Int.toNat_of_nonneg hm0]

/-- **`𝒪_{K̂^{ur}}` の任意のイデアルは単項**。

`I ≠ ⊥` なら `{n | ∃ a ∈ I, a ≠ 0 ∧ ‖a‖ = ‖π‖^n}` の最小元 `n₀` を取り、
それを実現する `a` が生成元。`b ∈ I` は `‖b‖ = ‖π‖^k`(`k ≥ n₀`)なので
`‖b/a‖ ≤ 1`、すなわち `b = (b/a)·a`。 -/
theorem isPrincipalIdealRing_unramifiedCompletionInt (K : PAdicLocalField p) :
    IsPrincipalIdealRing ↥(unramifiedCompletionInt K) := by
  classical
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  have hπpos : 0 < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hπle : ‖(π : K.carrier)‖ ≤ 1 := le_of_lt (norm_lt_one_of_maximalIdeal_eq_span K hπmax)
  constructor
  intro I
  rcases eq_or_ne I ⊥ with rfl | hI
  · exact ⟨⟨0, by simp⟩⟩
  have hex : ∃ n : ℕ, ∃ a : ↥(unramifiedCompletionInt K), a ∈ I ∧ a ≠ 0 ∧
      ‖(a : unramifiedCompletion K)‖ = ‖(π : K.carrier)‖ ^ n := by
    obtain ⟨a, haI, ha0⟩ := (Submodule.ne_bot_iff I).mp hI
    obtain ⟨n, hn⟩ := exists_norm_eq_pow_unramifiedCompletionInt K hπ0 hπmax ha0
    exact ⟨n, a, haI, ha0, hn⟩
  obtain ⟨a, haI, ha0, han⟩ := Nat.find_spec hex
  have ha0' : (a : unramifiedCompletion K) ≠ 0 := fun h => ha0 (Subtype.ext h)
  refine ⟨⟨a, le_antisymm ?_ ?_⟩⟩
  · intro b hbI
    rcases eq_or_ne b 0 with rfl | hb0
    · exact Submodule.zero_mem _
    obtain ⟨k, hk⟩ := exists_norm_eq_pow_unramifiedCompletionInt K hπ0 hπmax hb0
    have hkge : Nat.find hex ≤ k := by
      by_contra hc
      exact Nat.find_min hex (not_le.mp hc) ⟨b, hbI, hb0, hk⟩
    have hnorm_le : ‖(b : unramifiedCompletion K)‖ ≤ ‖(a : unramifiedCompletion K)‖ := by
      rw [hk, han]
      exact pow_le_pow_of_le_one (le_of_lt hπpos) hπle hkge
    rw [Ideal.mem_span_singleton']
    have hc : ‖(b : unramifiedCompletion K) / (a : unramifiedCompletion K)‖ ≤ 1 := by
      rw [norm_div, div_le_one (norm_pos_iff.mpr ha0')]
      exact hnorm_le
    refine ⟨⟨_, (mem_unramifiedCompletionInt K _).mpr hc⟩, ?_⟩
    apply Subtype.ext
    show (b : unramifiedCompletion K) / (a : unramifiedCompletion K)
        * (a : unramifiedCompletion K) = (b : unramifiedCompletion K)
    exact div_mul_cancel₀ _ ha0'
  · rw [Ideal.span_le, Set.singleton_subset_iff]
    exact haI

/-- **★★★★★★★★★★★★★★(Λ6-M1)`𝒪_{K̂^{ur}}` は離散付値環**。

単項イデアル環(上)＋局所環(`ValuationSubring`)＋極大イデアルが `⊥` でない
(`π ≠ 0`)。★結論に自由なパラメータが出ない形。

退化の自己検査:`K^{ur}` を `K.closure` に替えると**偽**——`𝒪_{ℂ_p}` の値群は
`ℚ` で稠密なので DVR ではない。効いているのは不分岐性(§1)である。 -/
theorem isDiscreteValuationRing_unramifiedCompletionInt (K : PAdicLocalField p) :
    IsDiscreteValuationRing ↥(unramifiedCompletionInt K) := by
  haveI := isPrincipalIdealRing_unramifiedCompletionInt K
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  refine IsDiscreteValuationRing.mk ?_
  intro hbot
  have hmem : uniformizerCompletionInt K π
      ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
    rw [maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax]
    exact Ideal.mem_span_singleton_self _
  rw [hbot, Ideal.mem_bot] at hmem
  have hzero : algebraMap K.carrier (unramifiedCompletion K) (π : K.carrier) = 0 :=
    congrArg Subtype.val hmem
  rw [← norm_eq_zero, norm_algebraMap_unramifiedCompletion] at hzero
  exact hπ0 (norm_eq_zero.mp hzero)

end ABC3.Found.PGC
