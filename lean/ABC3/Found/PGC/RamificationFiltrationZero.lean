import ABC3.Found.PGC.UnramifiedBaseChangeInvariance
import ABC3.Found.PGC.InertiaTransport

/-!
# `Γ_K^0 = I_K` —— ★★`Skeleton/PGC/Section2.lean` の `prop_2_2` が要求する最後の 1 本

`Skeleton/PGC/Section2.lean` の `prop_2_2` は「依拠する境界外の結果」として 3 つを挙げる:

1. `RamificationFiltration p`(Herbrand の定理) —— Y19e の
   `Found/PGC/UnramifiedBaseChangeInvariance.lean::ramificationFiltration` が
   **仮定ゼロ**で構成した。
2. 上付き↔下付き番号付けの変換 —— Y18 の
   `Found/PGC/UpperRamificationGroup.lean::upperRamificationGroup` ほか。
3. ★**`Γ_K^0 = I_K`(Corollary 1.3 の系)** —— ★★**本ファイルがこれを埋める。**

到達点は

```
ramificationFiltration_Gv_zero (K : PAdicLocalField p) :
    (ramificationFiltration p).Gv K 0 = absInertia K
```

であり、**仮定は 1 つも無い**。★`ramificationFiltration p` は Y19e が構成した
**本物**であって、退化 witness(`trivialStageFiltration` / `inertiaStageFiltration`)
ではない(§4 の自己検査を見よ)。

## 何が入ったか

### §1 抽象核(純位相群論) —— 分岐・付値・Galois の語彙が 1 つも出てこない

* `Subgroup.iInf_sup_eq_topologicalClosure` :
  `𝒩` が **開正規**部分群からなる `1` の近傍基のとき
  `⨅ N ∈ 𝒩, (H ⊔ N) = H` の閉包。★`H` について**何も仮定しない**一般形。
* `Subgroup.iInf_sup_eq_self_of_isClosed` : `H` が閉なら `⨅ N ∈ 𝒩, (H ⊔ N) = H`。
  ★こちらは `N` の開性すら要らない(`ContinuousMul` + 正規性 + 近傍基だけ)。
* `Subgroup.isClosed_iInf_sup` : ★**`H` が閉という仮定は落とせない**ことの表明
  ——`⨅ N ∈ 𝒩, (H ⊔ N)` は常に閉部分群だから、`H` が閉でなければ等式は偽である。

### §2 `StageFiltration` への橋(まだ抽象。位相群 `Γ` の話しかしない)

* `StageFiltration.limit_eq_of_S_eq_sup` :
  各段が `S N v = H ⊔ N`(`H` 閉)の形なら、逆極限は `H` そのもの。

### §3 具体層

* `absGalStageFiltration_limit_of_nonpos` : `v ≤ 0` で `Γ_K^v = I_K`。
* `ramificationFiltration_Gv_of_nonpos` / ★★`ramificationFiltration_Gv_zero`。
* `ramificationFiltration_Gv_zero_eq_inertia` :
  `Skeleton/PGC/Section1Defs.lean` の `inertia` との一致。
* `ramificationFiltration_Gv_zero_recoverable` :
  ★**`Γ_K^0` は `Γ_K` から群論的に復元できる**(Corollary 1.3 + 上の一致)。
  ★これが原文の「`Γ_K^0 = I_K`(Corollary 1.3 の系)」の中身である。
* `ramificationFiltration_Gv_le_absInertia` : `v ≥ 0` で `Γ_K^v ⊆ I_K`
  ——原典 §2 が使う `v > 0` の族に効く形。

### §4 退化の自己検査

* `absInertia_ne_top` : ★**`I_K ≠ Γ_K`**(`Γ_K/I_K` は任意の `ℤ/n` へ全射する)。
* `trivialStageFiltration_limit_zero` : Y19 の退化 witness では `Γ^0 = ⊤` であり、
  ★★**本定理は偽になる**(`trivialStageFiltration_limit_zero_ne_absInertia`、
  `absInertia_ne_top` により**無条件**)。
  すなわち本定理は「どんな段データでも成り立つ空虚な主張」ではない。
* `exists_coe_ramificationFiltration_mul_coe_ne_inertiaStage` :
  ★★**構成した `Γ_K^v` は Y19b+c の上からの近似 `inertiaStageFiltration` とも
  値が食い違う**(`v > 0` の側で)。すなわち `Γ_K^0 = I_K` を
  `inertiaStageFiltration`(そこでは自明に成り立つ)で「出した」ことにしていない。

## 数学(なぜ `v ≤ 0` で `I_K` なのか)

有限 Galois 拡大 `L/K` に対する上付き番号付けの分岐群は `Gal(L/K)^v = I(L/K)`
(`v ≤ 0`)である(Y19d `stageUpperRamification_of_nonpos`)。`Γ_K` への引き戻しは
`I_K ⊔ N`(`N = Gal(K̄/L)`、Y19d `absGalStage_of_nonpos`)。したがって

  `Γ_K^v = ⋂_{N 開正規} (I_K ⊔ N)`  (`v ≤ 0`)

であり、`I_K` は閉部分群(`isClosed_absInertia`)で、開正規部分群は `Γ_K` の `1` の
近傍基をなす(`absGal_exists_mem_openNormalBase_subset`)から、§1 の抽象核により
右辺は `I_K` に等しい。

★★**`I_K` が閉であることを落とすと偽になる**: §1 の `Subgroup.isClosed_iInf_sup` が
示すとおり、左辺は常に閉部分群であり、一般には `H` の閉包が出る。

## 逸脱の記録

* 原典 [pGC] は `Γ_K^v` を `v > 0` でしか使わない。`Interface/PGC/LocalFieldData.lean` の
  `RamificationFiltration` が(`Antitone` を書きやすくするため)全実数上で定義している
  ——これは Y19 の設計判断であり、本ファイルはそれに乗る。`v = 0` は
  原文の「`Γ_K^0 = I_K`」がまさに指す点である。
* import は指示された `ABC3.Found.PGC.UnramifiedBaseChangeInvariance` に加えて
  `ABC3.Found.PGC.InertiaTransport` を 1 本足した(§3 の
  `ramificationFiltration_Gv_zero_recoverable` に Corollary 1.3 の本体
  `inertia_recoverable_real` が要るため)。★他ファイルは変更していない。
* ★`Skeleton/**` `Interface/**` `1_Structured/**` `CorrHyp/**` には手を出していない。

原文 (pGC p.3):
> The inertia subgroup I_K ⊆ Γ_K can be determined group-theoretically from Γ_K.
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC
open scoped Pointwise

/-! ## §1 ★★抽象核 —— 純位相群論

★分岐・付値・Galois の語彙は 1 つも出てこない。位相群 `Γ`、その部分群 `H`、
`1` の近傍基をなす正規部分群の族 `𝒩` だけの話である。 -/

section AbstractCore

variable {Γ : Type*} [Group Γ] [TopologicalSpace Γ]

/-- ★★**抽象核(閉部分群版)**: `𝒩` が正規部分群からなる `1` の近傍基で `H` が閉なら

  `⨅ N ∈ 𝒩, (H ⊔ N) = H`.

★仮定は `ContinuousMul`(左移動の連続性)だけ。`𝒩` の元の**開性は要らない**。

★★**`H` が閉であることは落とせない**(`Subgroup.isClosed_iInf_sup` を見よ)。

証明: `⊇` は `le_sup_left`。`⊆` は、`x ∈ H ⊔ N = H · N`(`N` 正規)を
`x = h · n` と書けば `h = x n⁻¹ ∈ xN` であり、`N` が `1` の近傍基を走るので
`x` の任意の近傍が `H` と交わる——すなわち `x ∈ closure H = H`。 -/
theorem Subgroup.iInf_sup_eq_self_of_isClosed [ContinuousMul Γ] {H : Subgroup Γ}
    (hH : IsClosed (H : Set Γ)) {𝒩 : Set (Subgroup Γ)}
    (hnormal : ∀ ⦃N⦄, N ∈ 𝒩 → N.Normal)
    (hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ 𝒩, (N : Set Γ) ⊆ U) :
    ⨅ N : 𝒩, (H ⊔ (N : Subgroup Γ)) = H := by
  refine le_antisymm ?_ (le_iInf fun _ => le_sup_left)
  intro x hx
  rw [Subgroup.mem_iInf] at hx
  show x ∈ (H : Set Γ)
  rw [← hH.closure_eq, mem_closure_iff_nhds]
  intro U hU
  have hcont : ContinuousAt (fun g : Γ => x * g) 1 :=
    (continuous_const.mul continuous_id).continuousAt
  have hU1 : (fun g : Γ => x * g) ⁻¹' U ∈ nhds (1 : Γ) := by
    apply hcont.preimage_mem_nhds
    simpa using hU
  obtain ⟨N, hN, hNU⟩ := hbasis _ hU1
  haveI := hnormal hN
  have hmem : x ∈ ((H ⊔ N : Subgroup Γ) : Set Γ) := hx ⟨N, hN⟩
  -- ★`lean-idioms.md` #147: mathlib の `Subgroup.mul_normal` は `↑(H ⊔ N) = ↑H * ↑N` の向き。
  rw [Subgroup.mul_normal] at hmem
  obtain ⟨h, hh, n, hn, hhn⟩ := hmem
  refine ⟨h, ?_, hh⟩
  have hx' : h = x * n⁻¹ := by rw [← hhn]; group
  rw [hx']
  exact hNU (inv_mem hn)

/-- ★**抽象核(一般形)**: `𝒩` が**開**正規部分群からなる `1` の近傍基なら

  `⨅ N ∈ 𝒩, (H ⊔ N) = H` の位相閉包

——`H` について何も仮定しない。★上の閉部分群版はこの系である。 -/
theorem Subgroup.iInf_sup_eq_topologicalClosure [IsTopologicalGroup Γ] (H : Subgroup Γ)
    {𝒩 : Set (Subgroup Γ)} (hopen : ∀ ⦃N⦄, N ∈ 𝒩 → IsOpen (N : Set Γ))
    (hnormal : ∀ ⦃N⦄, N ∈ 𝒩 → N.Normal)
    (hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ 𝒩, (N : Set Γ) ⊆ U) :
    ⨅ N : 𝒩, (H ⊔ (N : Subgroup Γ)) = H.topologicalClosure := by
  refine le_antisymm ?_ ?_
  · intro x hx
    rw [Subgroup.mem_iInf] at hx
    show x ∈ ((H.topologicalClosure : Subgroup Γ) : Set Γ)
    rw [Subgroup.topologicalClosure_coe, mem_closure_iff_nhds]
    intro U hU
    have hcont : ContinuousAt (fun g : Γ => x * g) 1 :=
      (continuous_const.mul continuous_id).continuousAt
    have hU1 : (fun g : Γ => x * g) ⁻¹' U ∈ nhds (1 : Γ) := by
      apply hcont.preimage_mem_nhds
      simpa using hU
    obtain ⟨N, hN, hNU⟩ := hbasis _ hU1
    haveI := hnormal hN
    have hmem : x ∈ ((H ⊔ N : Subgroup Γ) : Set Γ) := hx ⟨N, hN⟩
    rw [Subgroup.mul_normal] at hmem
    obtain ⟨h, hh, n, hn, hhn⟩ := hmem
    refine ⟨h, ?_, hh⟩
    have hx' : h = x * n⁻¹ := by rw [← hhn]; group
    rw [hx']
    exact hNU (inv_mem hn)
  · exact le_iInf fun N => H.topologicalClosure_minimal le_sup_left
      (Subgroup.isClosed_of_isOpen _ (Subgroup.isOpen_mono le_sup_right (hopen N.2)))

/-- ★★**`H` が閉であるという仮定は落とせない**: `𝒩` が開部分群からなるなら
`⨅ N ∈ 𝒩, (H ⊔ N)` は**常に閉部分群**である。

したがって `H` が閉でなければ `⨅ N ∈ 𝒩, (H ⊔ N) = H` は偽になる
(実際に出るのは `H` の閉包である——`Subgroup.iInf_sup_eq_topologicalClosure`)。 -/
theorem Subgroup.isClosed_iInf_sup [IsTopologicalGroup Γ] (H : Subgroup Γ)
    {𝒩 : Set (Subgroup Γ)} (hopen : ∀ ⦃N⦄, N ∈ 𝒩 → IsOpen (N : Set Γ)) :
    IsClosed (((⨅ N : 𝒩, (H ⊔ (N : Subgroup Γ))) : Subgroup Γ) : Set Γ) := by
  rw [Subgroup.coe_iInf]
  exact isClosed_iInter fun N =>
    Subgroup.isClosed_of_isOpen _ (Subgroup.isOpen_mono le_sup_right (hopen N.2))

end AbstractCore

/-! ## §2 `StageFiltration` への橋(まだ抽象)

★ここも位相群 `Γ` の話しかしない。 -/

/-- ★★**段データの逆極限が `H` になる十分条件**:
`base` が `1` の近傍基で、各段が `S N v = H ⊔ N`(`H` 閉)の形をしているなら
`Γ^v = H`。

★`StageFiltration.limit` は `⨅ N : base, S N v` なので、§1 の抽象核をそのまま
当てるだけである。 -/
theorem StageFiltration.limit_eq_of_S_eq_sup {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [ContinuousMul Γ] (F : StageFiltration Γ)
    (hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ F.base, (N : Set Γ) ⊆ U)
    {H : Subgroup Γ} (hH : IsClosed (H : Set Γ)) {v : ℝ}
    (hS : ∀ ⦃N⦄, N ∈ F.base → F.S N v = H ⊔ N) :
    F.limit v = H := by
  have hcongr : F.limit v = ⨅ N : F.base, (H ⊔ (N : Subgroup Γ)) := by
    rw [StageFiltration.limit]
    exact iInf_congr fun N => hS N.2
  rw [hcongr]
  exact Subgroup.iInf_sup_eq_self_of_isClosed hH (fun _ hN => F.normal_base hN) hbasis

/-! ## §3 ★★★具体層 —— `Γ_K^0 = I_K` -/

variable {p : ℕ} [Fact p.Prime]

/-- ★`Γ_K` では開正規部分群が段データの `base` をなし、それが `1` の近傍基である。 -/
theorem absGalStageFiltration_basis (K : PAdicLocalField p) (hcompat)
    (U : Set K.absGal) (hU : U ∈ nhds (1 : K.absGal)) :
    ∃ N ∈ (absGalStageFiltration K hcompat).base, (N : Set K.absGal) ⊆ U := by
  obtain ⟨N, hN, hNU⟩ := absGal_exists_mem_openNormalBase_subset K hU
  exact ⟨N, hN, hNU⟩

/-- ★★**`v ≤ 0` では `Γ_K^v = I_K`** —— Y19d の `absGalStage_of_nonpos`
(各段が `I_K ⊔ N`)を §1 の抽象核に通しただけ。 -/
theorem absGalStageFiltration_limit_of_nonpos (K : PAdicLocalField p) (hcompat) {v : ℝ}
    (hv : v ≤ 0) : (absGalStageFiltration K hcompat).limit v = absInertia K := by
  refine (absGalStageFiltration K hcompat).limit_eq_of_S_eq_sup
    (absGalStageFiltration_basis K hcompat) (isClosed_absInertia K) ?_
  intro N hN
  exact absGalStage_of_nonpos K hN hv

/-- ★★**Y19e の本物のフィルトレーションについて `v ≤ 0` で `Γ_K^v = I_K`**。 -/
theorem ramificationFiltration_Gv_of_nonpos (K : PAdicLocalField p) {v : ℝ} (hv : v ≤ 0) :
    (ramificationFiltration p).Gv K v = absInertia K := by
  show (absGalStageFiltration K (absGalStage_compat unramifiedBaseChange K)).limit v = _
  exact absGalStageFiltration_limit_of_nonpos K _ hv

/-- ★★★★★★**`Γ_K^0 = I_K`** —— `Skeleton/PGC/Section2.lean` の `prop_2_2` が
「依拠する境界外の結果」として挙げた 3 つのうちの最後の 1 本。

★★**仮定は 1 つも無い。** `ramificationFiltration p` は Y19e が
`Found/PGC/UnramifiedBaseChangeInvariance.lean` で構成した**本物**である
(退化 witness ではない——§4 の自己検査)。 -/
theorem ramificationFiltration_Gv_zero (K : PAdicLocalField p) :
    (ramificationFiltration p).Gv K 0 = absInertia K :=
  ramificationFiltration_Gv_of_nonpos K le_rfl

def ramificationFiltration_Gv_zero.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- ★`Skeleton/PGC/Section1Defs.lean` の `inertia`(Corollary 1.3 の主語)との一致。 -/
theorem ramificationFiltration_Gv_zero_eq_inertia (K : PAdicLocalField p) :
    (ramificationFiltration p).Gv K 0
      = inertia (residueCardinality p) (subgroupCorrespondence p) K := by
  rw [ramificationFiltration_Gv_zero, absInertia_eq_inertia]

/-- ★★**原文の「`Γ_K^0 = I_K`(Corollary 1.3 の系)」の中身**:
`Γ_K^0` は `Γ_K` から群論的に復元できる。

`Γ_K^0 = I_K`(上)と Corollary 1.3(`inertia_recoverable_real`)の合成である。
★`prop_2_2` は「`v > 0` の族が与えられている」という設定なので、`v = 0` の項が
群論的に復元できることを別途言う必要がある——それがこれである。 -/
theorem ramificationFiltration_Gv_zero_recoverable :
    (inertiaObject (residueCardinality p) (subgroupCorrespondence p)).RecoverableFromAbsGal :=
  inertia_recoverable_real

def ramificationFiltration_Gv_zero_recoverable.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Corollary 1.3", sectionId := "cor-1-3" }

/-- ★`v ≥ 0` では `Γ_K^v ⊆ I_K` —— 反単調性と `Γ_K^0 = I_K` から即座に出る。
★原典 §2 が使うのは `v > 0` の族なので、実際に効くのはこの形である。 -/
theorem ramificationFiltration_Gv_le_absInertia (K : PAdicLocalField p) {v : ℝ} (hv : 0 ≤ v) :
    (ramificationFiltration p).Gv K v ≤ absInertia K := by
  rw [← ramificationFiltration_Gv_zero K]
  exact (ramificationFiltration p).antitone K hv

/-! ## §4 ★★退化の自己検査

★★「どんな段データでも `Γ^0 = I_K` になる」のではないこと、および
「`inertiaStageFiltration`(上からの近似、そこでは自明)で出したのではない」ことを
型で検査する。 -/

/-- Y19 の退化 witness `trivialStageFiltration` では `Γ^0 = ⊤`。 -/
theorem trivialStageFiltration_limit_zero (Γ : Type*) [Group Γ] [TopologicalSpace Γ] :
    (trivialStageFiltration Γ).limit 0 = ⊤ := by
  rw [eq_top_iff]
  intro x _
  rw [StageFiltration.mem_limit_iff]
  intro N _
  show x ∈ (if (0 : ℝ) ≤ 0 then (⊤ : Subgroup Γ) else N)
  simp

/-- ★★**`I_K ≠ Γ_K`** —— `K` には必ず非自明な不分岐拡大がある。

`Γ_K/I_K ≅ Gal(K^ur/K)` は任意の `n ≥ 1` に対して `ℤ/n` へ全射する
(`exists_surjective_quotKer_to_zmod`)ので、`I_K = ⊤` なら `ℤ/2` が自明になってしまう。

★これで下の退化検査が**無条件**になる。 -/
theorem absInertia_ne_top (K : PAdicLocalField p) : absInertia K ≠ ⊤ := by
  intro htop
  have hker : ∀ σ : K.absGal, σ ∈ (AlgEquiv.restrictNormalHom
      (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K)).ker := by
    intro σ
    rw [ker_restrictNormalHom_unramifiedClosure K, ← absInertia_eq_inertia K, htop]
    trivial
  obtain ⟨φ, hφ⟩ := exists_surjective_quotKer_to_zmod K 2 two_ne_zero
  obtain ⟨q, hq⟩ := hφ (Multiplicative.ofAdd (1 : ZMod 2))
  obtain ⟨σ, rfl⟩ := QuotientGroup.mk_surjective q
  have h1 : (QuotientGroup.mk σ : K.absGal ⧸ (AlgEquiv.restrictNormalHom
      (F := K.carrier) (K₁ := K.closure) (unramifiedClosure K)).ker) = 1 :=
    (QuotientGroup.eq_one_iff σ).2 (hker σ)
  rw [h1, map_one] at hq
  have h2 : (0 : ZMod 2) = 1 := by simpa using congrArg Multiplicative.toAdd hq
  exact absurd h2 (by decide)

/-- ★★**本定理は空虚ではない**: 退化 witness `trivialStageFiltration` に対して
`Γ^0 = I_K` は**偽**である(そこでは `Γ^0 = ⊤ ≠ I_K`)。 -/
theorem trivialStageFiltration_limit_zero_ne_absInertia (K : PAdicLocalField p) :
    (trivialStageFiltration K.absGal).limit 0 ≠ absInertia K := by
  rw [trivialStageFiltration_limit_zero]
  exact fun h => absInertia_ne_top K h.symm

/-- ★★★**構成した `Γ_K^v` は Y19b+c の上からの近似 `inertiaStageFiltration` と
値が食い違う**。

`I_K ≰ N` なる開正規部分群 `N` に対し、ある `v` で
`Γ_K^v · N ≠ (I_K ⊔ N)` である(`inertiaStageFiltration` の第 `N` 段は
`v ≥ 0` で常に `I_K ⊔ N`)。★すなわち §3 は `inertiaStageFiltration`(そこでは
`Γ^0 = I_K` が自明に成り立つ)で「出した」ものではない。 -/
theorem exists_coe_ramificationFiltration_mul_coe_ne_inertiaStage (K : PAdicLocalField p)
    {N : Subgroup K.absGal} (hN : N ∈ openNormalBase K.absGal) (hne : ¬ absInertia K ≤ N) :
    ∃ v : ℝ, (((ramificationFiltration p).Gv K v : Subgroup K.absGal) : Set K.absGal)
        * (N : Set K.absGal)
      ≠ (((inertiaStageFiltration K).S N v : Subgroup K.absGal) : Set K.absGal) := by
  obtain ⟨v, hv⟩ := exists_absGalStage_ne_inertiaStageFiltration K hN hne
  refine ⟨v, ?_⟩
  rw [coe_ramificationFiltration_mul_coe K hN v]
  exact fun h => hv (SetLike.coe_injective h)

#print axioms ramificationFiltration_Gv_zero
#print axioms absInertia_ne_top
#print axioms trivialStageFiltration_limit_zero_ne_absInertia
#print axioms ramificationFiltration_Gv_zero_eq_inertia
#print axioms ramificationFiltration_Gv_zero_recoverable
#print axioms Subgroup.iInf_sup_eq_self_of_isClosed
#print axioms Subgroup.iInf_sup_eq_topologicalClosure

end ABC3.Found.PGC
