import ABC3.Found.PGC.AbsGalRamificationFiltration
import ABC3.Found.PGC.LubinTateClosureTopology
import ABC3.Found.PGC.LubinTateReciprocityMapLimitSurjective

/-!
# `Γ^v_K` の `Γ^ab_K` における像 = `U^v_K` —— pGC §2 の分岐入力

典拠: S. Mochizuki, *A Version of the Grothendieck Conjecture for p-adic Local Fields*
(pGC) 物理 p.4(Proposition 2.1 の直後の段落)。原典が
"it is well-known (Theorem 1 of [3], p. 155)" と 1 語で畳んでいる箇所であり、
[3] は J.-P. Serre, *Corps Locaux* (Local Fields) 第 XV 章 §2 Theorem 1 である。

原文 (pGC p.4):
> Then we shall denote by Γ^v_K ⊆ Γ_K the higher ramification group associated to the number v in the "upper numbering" (see, e.g., [3], p. 155).

原文 (pGC p.4):
> Let us denote by U^v_K ⊆ U_K the subgroup 1 + m^n_K, where m_K ⊆ O_K is the maximal ideal, and n is the unique integer for which n − 1 < v ≤ n.

原文 (pGC p.4):
> Then it is well-known (Theorem 1 of [3], p. 155) that the image of Γ^v_K in Γ^ab_K is equal to U^v_K ⊆ U_K.

## ★★添字のずれ —— 最初に確かめたこと

持ち場が名指しした 2 つの「ずれ」は **同じ**であり、打ち消し合う。

| どこ | ずれ |
|---|---|
| 原典 | `U^v_K = 1 + m^n_K`、`n` は `n − 1 < v ≤ n` を満たす唯一の整数。`v = n : ℕ` なら `n = n` |
| 在庫 `mem_principalUnits_reciprocityUnits_iff` | `reciprocityUnits σ ∈ principalUnits K π (m+1) ↔ σ x_m = x_m` |

在庫の `psiGenSeq ... m` は **`ψ_{m+1}` の根**(原始的な `π^{m+1}`-捩れ点)であり、
`K(x_m) = K_{f,m+1}` は Lubin-Tate 塔の **第 `m+1` 段**である。すなわち在庫の
`(m+1)` は「レベル `m+1` の主単数群 ↔ 第 `m+1` 段」という**ずれの無い**対応であって、
`psiGenSeq` の添字が 0 始まりなだけである。したがって

    `Art^{-1}(Γ^v_K) ⊆ U^n_K`  (`n = ⌈v⌉`)

を出すのに **1 のずれは入らない**。本ファイルは `⌈v⌉₊`(`Nat.ceil`)で原典の丸めを実装し、
`natCeil_eq_of_sub_one_lt_le` で「`⌈v⌉₊` が原典の `n` に一致する(`v ≥ 0`)」ことを、
`eq_of_sub_one_lt_le` で「原典の `n` が一意である」ことを証明した。
★`v = 0` でも一致する: 原典の条件 `n − 1 < 0 ≤ n` は `n = 0` を与え、
`U^0_K = 1 + m^0_K = O_K^×`。本ファイルの `principalUnits_zero_eq_top` が `= ⊤` を言う。

## 本ファイルが埋めたもの / 埋めていないもの(★正直な線引き)

**埋めたもの**

* §1–§4 **抽象核**(★分岐・付値・Galois の語彙が 1 語も出てこない):
  - `coe_mul_coe_iInter_of_directed` —— コンパクト位相群で
    `S · ⋂ T_i = ⋂ (S · T_i)`(`S` 閉、`T_i` 閉で下降有向)。
  - `map_eq_of_coe_mul_coe_ker_eq` —— ★**純群論・選択公理を使わない**:
    全射 `ρ` について `S · ker ρ = ρ^{-1}(U)` なら `ρ(S) = U`。
  - `map_eq_of_coe_mul_coe_eq_comap` —— 上の 2 つの合成。
  - `StageFiltration.map_limit_le_of_stage_le_comap` /
    `StageFiltration.map_limit_eq_of_stage_eq_comap` —— Y19 の段データへの橋。
  - `StageFiltration.ofAntitoneNormal`(+ `_limit`)—— 反単調な正規部分群の族から
    段データを作る構成子。
  - `eq_of_sub_one_lt_le` —— 「`n − 1 < v ≤ n` を満たす整数は高々 1 つ」。
* §5–§6 **具体層**:
  - `iInf_principalUnits_eq_bot` —— `⋂_n (1 + m^n) = 1`(`IsHausdorff`)。
  - `reciprocityUnits_surjective`。
  - `absGalPrincipalLevel n := Art^{-1}(U^n_K) ⊆ Γ_K` と、それが
    **`Gal(K̄/K_{f,n})`**(`= K(x_{n-1})` の固定化部分群)に等しいこと
    (`absGalPrincipalLevel_eq_fixingSubgroup`)、開かつ正規であること。
* §7 **主定理**:
  - `map_limit_reciprocityUnits_le_principalUnits`(★包含 `⊆`。コンパクト性は不要)、
  - `map_limit_reciprocityUnits_eq_principalUnits` /
    `map_limit_reciprocityUnits_eq_upperPrincipalUnits`(★等号)。
* §8 **仮定の非空虚性**: `exists_stageFiltration_map_limit_eq_upperPrincipalUnits`
  —— §7 の仮定 2 本を**同時に満たす段データが存在する**ことを構成で示す。

**埋めていないもの(★次のノード)**

★★**「有限段の上付き分岐群」を与える段データ `F` 自身は本ファイルでは作らない。**
`Found/PGC/RamificationFiltrationBuild.lean` が `absGalStage` を作っているが
`compat` が残っており、`StageFiltration K.absGal` はまだ立っていない。本ファイルは
**その `F` を仮定として受け取り、`F` に課される条件を 2 本だけに絞った**:

1. `absGalPrincipalLevel (n+i) ∈ F.base`(★本ファイルが証明済み ——
   `absGalPrincipalLevel_mem_openNormalBase`。`F.base = openNormalBase` なら自動)、
2. `F.S (absGalPrincipalLevel (n+i)) v = absGalPrincipalLevel n`
   —— これが原典の "well-known" の**有限段での中身**、すなわち
   `Gal(K_{f,n+i}/K)^n = U^n_K/U^{n+i}_K` である。

★★2 の**特殊形 `i = 0`(`Gal(K_{f,n}/K)^n = {id}`)だけで包含 `⊆` は出る**
(`map_limit_reciprocityUnits_le_principalUnits`)。木の在庫
`upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap`(Yoshida Prop 6.14 の
`n = 1` 形)がまさにその形であり、そこに残っている仮定は原典の段 1(`hρ`)1 本である。
★**等号を出すには 2 を `i` すべてで要る**(＝ Prop 6.14 の**鋭い形**)。
木にはまだ「消える」側しか無い。

## ★原典より短い道(名指し)

原典(および Serre XV §2)は上付き番号付けの Herbrand 関数を経由して
`Art(U^n) = Γ^n` を出すが、本ファイルの等号は
**「コンパクト群で閉部分群と下降有向な閉部分群族の積は交換する」**という
1 本の位相補題(`coe_mul_coe_iInter_of_directed`)しか使わない。
Herbrand 関数も微分も出てこない。核 `ker(Art)` を
`⋂_i Art^{-1}(U^{n+i})` として捉え直した(`iInf_absGalPrincipalLevel_eq_ker`)のが要点で、
これは `⋂_n (1 + m^n) = 1` すなわち `𝒪_K` の分離性そのものである。

## ★見立てが外れた点(隠さない)

当初は「在庫の `mem_principalUnits_reciprocityUnits_iff` の `(m+1)` と原典の `⌈v⌉` が
1 ずれる」と見込んでいた。★**外れていた**(上の「添字のずれ」の節)。
`psiGenSeq` が 0 始まりなだけで、レベルと段は一致している。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★**`Γ^ab_K` ではなく `𝒪_K^×` への写像 `reciprocityUnits` で述べた。**
   原典は「`Γ^v_K` の `Γ^ab_K` における像が `U^v_K`」と言い、`U_K = 𝒪_K^×` を
   `Γ^ab_K` の部分群と見ている。木にあるのは Lubin-Tate 相互律
   `reciprocityUnits : Γ_K →* 𝒪_K^×`(`Γ_K → Γ^ab_K` の `𝒪_K^×` 成分)なので、
   本ファイルはその成分での等式を述べる。両者が一致するには
   「`v ≥ 0` では `Γ^v_K ⊆ Γ^0_K =` 惰性群」(不分岐成分が消えること)が要る。
   ★**その包含は本ファイルでは証明していない。** 消費側が `Γ^ab_K` の言葉を要るなら
   1 ノード足すこと。
2. ★**`v` は実数全体で定義した**(原典は `v ≥ 0`)。`upperPrincipalUnits K π v :=
   principalUnits K π ⌈v⌉₊` なので `v < 0` では `⊤` になる。`v ≥ 0` では原典と一致する
   (`natCeil_eq_of_sub_one_lt_le`)。Y19 の逸脱 1(`StageFiltration` を `ℝ` 全体で見る)
   を引き継いだだけである。
3. ★**主定理は自然数 `v = n` の場合である。** 実数 `v` の場合は
   「段データ `F` の側で `F.S N v` が `⌈v⌉` にしか依らない」ことが要り、
   それは `F` の中身(＝まだ立っていない有限段の上付き分岐群)に踏み込む。
   ★**別ノードにした。**
4. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ)。
   `Found/PGC/LubinTate*.lean` は読んだだけである(D28)。
5. ★本ファイルは原典が名前を付けていない段落なので `.src` を持たない
   (`AbsGalRamificationFiltration.lean` / `RamificationFiltrationBuild.lean` と同じ理由)。
-/

namespace ABC3.Found.PGC

open scoped Pointwise

/-! ## §1 抽象核 —— コンパクト位相群での積と下降交わりの交換

★分岐・付値・Galois の語彙が 1 語も出てこない。 -/

/-- ★★**抽象核**: コンパクト位相群において、**閉**部分群 `S` と、**下降有向**な
**閉**部分群の族 `T` について

    `S · (⋂ᵢ Tᵢ) = ⋂ᵢ (S · Tᵢ)`.

`⊆` は自明。`⊇` は **コンパクト性**(下降有向な空でない閉集合族の交わりは空でない)で、
`x ∈ ⋂ᵢ (S · Tᵢ)` に対して `(x • Tᵢ) ∩ S` たちの交わりから元を取る。

★`S` が閉であることは落とせない(閉包を取ると等式が壊れる)。
★分岐・付値・Galois の語彙は 1 つも出てこない。 -/
theorem coe_mul_coe_iInter_of_directed {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] [CompactSpace Γ] {ι : Type*} [Nonempty ι]
    {S : Subgroup Γ} (hS : IsClosed (S : Set Γ))
    {T : ι → Subgroup Γ} (hT : ∀ i, IsClosed ((T i : Set Γ)))
    (hdir : Directed (· ≥ ·) T) :
    (S : Set Γ) * ((⨅ i, T i : Subgroup Γ) : Set Γ) = ⋂ i, ((S : Set Γ) * (T i : Set Γ)) := by
  refine Set.Subset.antisymm ?_ ?_
  · rintro _ ⟨a, ha, b, hb, rfl⟩
    exact Set.mem_iInter.2 fun i => ⟨a, ha, b, (Subgroup.mem_iInf.1 hb) i, rfl⟩
  · intro x hx
    simp only [Set.mem_iInter] at hx
    set t : ι → Set Γ := fun i => (x • (T i : Set Γ)) ∩ (S : Set Γ) with ht
    have htn : ∀ i, (t i).Nonempty := by
      intro i
      obtain ⟨s, hs, b, hb, hsb⟩ := hx i
      refine ⟨s, ?_, hs⟩
      rw [Set.mem_smul_set_iff_inv_smul_mem, smul_eq_mul]
      have hb2 : x⁻¹ * s = b⁻¹ := by rw [← hsb]; group
      rw [hb2]; exact inv_mem hb
    have htcl : ∀ i, IsClosed (t i) := fun i => ((hT i).smul x).inter hS
    have htd : Directed (· ⊇ ·) t := by
      intro i j
      obtain ⟨k, hki, hkj⟩ := hdir i j
      exact ⟨k, Set.inter_subset_inter_left _ (Set.smul_set_mono hki),
        Set.inter_subset_inter_left _ (Set.smul_set_mono hkj)⟩
    obtain ⟨y, hy⟩ := IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed
      t htd htn (fun i => (htcl i).isCompact) htcl
    simp only [Set.mem_iInter, ht, Set.mem_inter_iff] at hy
    refine ⟨y, (hy (Classical.arbitrary ι)).2, y⁻¹ * x, ?_, by group⟩
    refine Subgroup.mem_iInf.2 fun i => ?_
    have hyi := (hy i).1
    rw [Set.mem_smul_set_iff_inv_smul_mem, smul_eq_mul] at hyi
    simpa [mul_inv_rev] using inv_mem hyi

/-! ## §2 抽象核 —— 「像が `U`」を「`S · ker = ρ⁻¹(U)`」に翻訳する -/

/-- ★★**抽象核(純群論)**: `ρ` が全射で `S · ker ρ = ρ^{-1}(U)` なら `ρ(S) = U`。

★`#print axioms` が `[propext, Quot.sound]` ——**選択公理を使わない**。
★位相も分岐も出てこない。 -/
theorem map_eq_of_coe_mul_coe_ker_eq {Γ Q : Type*} [Group Γ] [Group Q] {ρ : Γ →* Q}
    (hρ : Function.Surjective ρ) {S : Subgroup Γ} {U : Subgroup Q}
    (h : (S : Set Γ) * ((ρ.ker : Subgroup Γ) : Set Γ)
      = ((Subgroup.comap ρ U : Subgroup Γ) : Set Γ)) :
    S.map ρ = U := by
  have hle : S ≤ Subgroup.comap ρ U := by
    intro s hs
    have hmem : s ∈ ((Subgroup.comap ρ U : Subgroup Γ) : Set Γ) := by
      rw [← h]; exact ⟨s, hs, 1, one_mem _, mul_one s⟩
    exact hmem
  refine le_antisymm (Subgroup.map_le_iff_le_comap.2 hle) ?_
  intro u hu
  obtain ⟨x, rfl⟩ := hρ u
  have hx : x ∈ ((Subgroup.comap ρ U : Subgroup Γ) : Set Γ) := hu
  rw [← h] at hx
  obtain ⟨s, hs, k, hk, rfl⟩ := hx
  exact ⟨s, hs, by rw [map_mul, MonoidHom.mem_ker.1 hk, mul_one]⟩

/-- ★★**抽象核**(§1 と §2 の合成): コンパクト位相群 `Γ`、全射 `ρ : Γ →* Q`、
**閉**部分群 `S`、下降有向な**閉**部分群の族 `T` で `⋂ᵢ Tᵢ = ker ρ` のとき、
`S · Tᵢ = ρ^{-1}(U)` がすべての `i` で成り立てば `ρ(S) = U`。 -/
theorem map_eq_of_coe_mul_coe_eq_comap {Γ Q : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] [CompactSpace Γ] [Group Q] {ρ : Γ →* Q}
    (hρ : Function.Surjective ρ) {S : Subgroup Γ} (hS : IsClosed (S : Set Γ))
    {ι : Type*} [Nonempty ι] {T : ι → Subgroup Γ} (hTcl : ∀ i, IsClosed ((T i : Set Γ)))
    (hdir : Directed (· ≥ ·) T) (hker : (⨅ i, T i) = ρ.ker)
    {U : Subgroup Q}
    (h : ∀ i, (S : Set Γ) * ((T i : Subgroup Γ) : Set Γ)
      = ((Subgroup.comap ρ U : Subgroup Γ) : Set Γ)) :
    S.map ρ = U := by
  refine map_eq_of_coe_mul_coe_ker_eq hρ ?_
  rw [← hker, coe_mul_coe_iInter_of_directed hS hTcl hdir]
  simp only [h]
  exact Set.iInter_const _

/-! ## §3 段フィルトレーション(Y19)への橋 -/

namespace StageFiltration

variable {Γ : Type*} [Group Γ] [TopologicalSpace Γ] (F : StageFiltration Γ)

/-- ★**包含 `⊆` の側**: ある段 `N ∈ base` で `S N v ≤ ρ^{-1}(U)` なら
`ρ(Γ^v) ≤ U`。★コンパクト性も全射性も要らない。 -/
theorem map_limit_le_of_stage_le_comap {Q : Type*} [Group Q] (ρ : Γ →* Q) {N : Subgroup Γ}
    (hN : N ∈ F.base) {v : ℝ} {U : Subgroup Q} (hle : F.S N v ≤ Subgroup.comap ρ U) :
    (F.limit v).map ρ ≤ U :=
  Subgroup.map_le_iff_le_comap.2 ((F.limit_le hN v).trans hle)

/-- ★★**等号**: 段の族 `T` が `base` に属し、下降有向で `⋂ᵢ Tᵢ = ker ρ` を満たし、
どの段でも `S Tᵢ v = ρ^{-1}(U)` なら `ρ(Γ^v) = U`。

★Y19 の `coe_limit_mul_coe`(極限が各段へ全射する)を §2 の抽象核に流し込むだけ。 -/
theorem map_limit_eq_of_stage_eq_comap [IsTopologicalGroup Γ] [CompactSpace Γ]
    {Q : Type*} [Group Q] {ρ : Γ →* Q}
    (hρ : Function.Surjective ρ) {ι : Type*} [Nonempty ι] {T : ι → Subgroup Γ}
    (hbase : ∀ i, T i ∈ F.base) (hdir : Directed (· ≥ ·) T) (hker : (⨅ i, T i) = ρ.ker)
    {v : ℝ} {U : Subgroup Q} (h : ∀ i, F.S (T i) v = Subgroup.comap ρ U) :
    (F.limit v).map ρ = U :=
  map_eq_of_coe_mul_coe_eq_comap (T := T) hρ (F.isClosed_limit v)
    (fun i => Subgroup.isClosed_of_isOpen _ (F.isOpen_base (hbase i))) hdir hker
    (fun i => by rw [F.coe_limit_mul_coe (hbase i) v, h i])

/-- ★**抽象核**: 反単調な正規部分群の族 `A : ℝ → Subgroup Γ` から段データ
`S N v := A v ⊔ N` を作る。`compat` は `N ≤ M` のとき
`(A v ⊔ N) · M = A v ⊔ (N ⊔ M) = A v ⊔ M` という束の計算だけである。 -/
def ofAntitoneNormal (A : ℝ → Subgroup Γ) (hnormal : ∀ v, (A v).Normal)
    (hanti : Antitone A) : StageFiltration Γ :=
  StageFiltration.ofOpenNormal (fun N v => A v ⊔ N)
    (fun _ _ _ => le_sup_right)
    (fun N hN v => by haveI := hnormal v; haveI := hN.2; infer_instance)
    (fun _ _ _ _ hab => sup_le_sup_right (hanti hab) _)
    (fun M hM N hN hNM v => by
      haveI : M.Normal := hM.2
      rw [← Subgroup.mul_normal, sup_assoc, sup_eq_right.2 hNM])

@[simp] theorem ofAntitoneNormal_S (A : ℝ → Subgroup Γ) (h1 : ∀ v, (A v).Normal)
    (h2 : Antitone A) (N : Subgroup Γ) (v : ℝ) :
    (ofAntitoneNormal A h1 h2).S N v = A v ⊔ N := rfl

@[simp] theorem ofAntitoneNormal_base (A : ℝ → Subgroup Γ) (h1 : ∀ v, (A v).Normal)
    (h2 : Antitone A) : (ofAntitoneNormal A h1 h2).base = openNormalBase Γ := rfl

/-- `ofAntitoneNormal` の極限は `A v` そのもの(`A v` が閉で、`base` が `1` の近傍基のとき)。
★Y19 の一意性 `eq_limit_of_isClosed` の直接の帰結。コンパクト性は要らない。 -/
theorem ofAntitoneNormal_limit [IsTopologicalGroup Γ] (A : ℝ → Subgroup Γ)
    (h1 : ∀ v, (A v).Normal) (h2 : Antitone A)
    (hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ openNormalBase Γ, (N : Set Γ) ⊆ U)
    (hcl : ∀ v, IsClosed ((A v : Subgroup Γ) : Set Γ)) (v : ℝ) :
    (ofAntitoneNormal A h1 h2).limit v = A v :=
  ((ofAntitoneNormal A h1 h2).eq_limit_of_isClosed hbasis (hcl v)
    (fun M hM => by haveI : M.Normal := hM.2; rw [← Subgroup.mul_normal]; rfl)).symm

end StageFiltration

/-! ## §4 抽象核 —— 原典の丸め `n − 1 < v ≤ n` -/

/-- ★**抽象核(実数と整数)**: 原典の「`n − 1 < v ≤ n` を満たす**唯一の**整数 `n`」の
**一意性**。★分岐も付値も出てこない。 -/
theorem eq_of_sub_one_lt_le {v : ℝ} {m n : ℤ} (hm1 : (m : ℝ) - 1 < v) (hm2 : v ≤ m)
    (hn1 : (n : ℝ) - 1 < v) (hn2 : v ≤ n) : m = n := by
  have h1 : (m : ℝ) < n + 1 := by linarith
  have h2 : (n : ℝ) < m + 1 := by linarith
  have h1' : m < n + 1 := by exact_mod_cast h1
  have h2' : n < m + 1 := by exact_mod_cast h2
  omega

/-! ## §5 単数の側 —— `U^v_K` -/

open ABC3.Skeleton.PGC
open IsLocalRing
open scoped NormedField Valued

section Units

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p) (π : 𝒪[K.carrier])

/-- `U^0_K = 1 + m^0_K = 𝒪_K^×`。★原典の丸めは `v = 0` で `n = 0` を与えるので、
`U^0_K = U_K` である。 -/
theorem principalUnits_zero_eq_top : principalUnits K π 0 = ⊤ := by
  ext u
  simp only [Subgroup.mem_top, iff_true, mem_principalUnits_iff]
  exact ⟨(u : 𝒪[K.carrier]) - 1, by ring⟩

/-- `n ↦ U^n_K` は反単調(在庫 `principalUnits_succ_le` から)。 -/
theorem antitone_principalUnits : Antitone (principalUnits K π) :=
  antitone_nat_of_succ_le (principalUnits_succ_le K π)

/-- ★**原典の `U^v_K`**(`v : ℝ`): `1 + m^n_K`、`n` は `n − 1 < v ≤ n` を満たす唯一の整数。
`v ≥ 0` では `n = ⌈v⌉₊`(`natCeil_eq_of_sub_one_lt_le`)。

★逸脱: `v < 0` でも定義される(`⌈v⌉₊ = 0` なので `⊤`)。冒頭「逸脱の記録 2」。 -/
noncomputable def upperPrincipalUnits (v : ℝ) : Subgroup (𝒪[K.carrier])ˣ :=
  principalUnits K π ⌈v⌉₊

@[simp] theorem upperPrincipalUnits_natCast (n : ℕ) :
    upperPrincipalUnits K π (n : ℝ) = principalUnits K π n := by
  rw [upperPrincipalUnits, Nat.ceil_natCast]

end Units

/-- ★**`⌈v⌉₊` が原典の丸めと一致する**(`v ≥ 0`)。`v = 0` の場合も含む
(そこでは `n = 0`、`U^0_K = 𝒪_K^×`)。 -/
theorem natCeil_eq_of_sub_one_lt_le {v : ℝ} (hv : 0 ≤ v) {n : ℕ}
    (h1 : (n : ℝ) - 1 < v) (h2 : v ≤ n) : ⌈v⌉₊ = n := by
  rcases Nat.eq_zero_or_pos n with hn | hn
  · subst hn
    have hz : v = 0 := le_antisymm (by simpa using h2) hv
    simp [hz]
  · refine (Nat.ceil_eq_iff (by omega)).2 ⟨?_, h2⟩
    rwa [Nat.cast_sub hn, Nat.cast_one]

/-- ★★**`⋂_n (1 + m^n_K) = 1`** —— `𝒪_K` の `m`-進分離性(`IsAdicComplete` から出る
`IsHausdorff`)そのもの。これが §7 の「核の同定」の中身である。 -/
theorem iInf_principalUnits_eq_bot {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    (⨅ n : ℕ, principalUnits K π n) = ⊥ := by
  refine le_antisymm ?_ bot_le
  intro u hu
  rw [Subgroup.mem_iInf] at hu
  have hz : (u : 𝒪[K.carrier]) - 1 = 0 := by
    refine IsHausdorff.haus (I := IsLocalRing.maximalIdeal 𝒪[K.carrier])
      inferInstance _ (fun n => ?_)
    rw [SModEq.sub_mem]
    have h2 : (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ^ n
        = Ideal.span ({π ^ n} : Set 𝒪[K.carrier]) := by
      rw [hπmax, Ideal.span_singleton_pow]
    have h1 : ((u : 𝒪[K.carrier]) - 1) - 0 ∈ (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ^ n := by
      rw [h2, sub_zero]; exact hu n
    simpa using h1
  have hu1 : u = 1 := Units.ext (by rw [Units.val_one]; linear_combination hz)
  exact Subgroup.mem_bot.2 hu1

/-! ## §6 具体層 —— `Γ_K` の中の `Art^{-1}(U^n_K)` -/

section Concrete

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))

/-- Lubin-Tate 相互律 `Γ_K →* 𝒪_K^×` は**全射**
(在庫 `reciprocityMapLimit_surjective` に同型を 1 つ合成するだけ)。 -/
theorem reciprocityUnits_surjective :
    Function.Surjective (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf) := by
  intro u
  obtain ⟨σ, hσ⟩ := reciprocityMapLimit_surjective K hq hπmax hπne0 f hf0 hf1 hf
    (unitsEquivCompatibleUnits K hπmax u)
  exact ⟨σ, by simp [reciprocityUnits, hσ]⟩

/-- ★`Γ_K` の中の**第 `n` レベル** `Art^{-1}(U^n_K)`。
`absGalPrincipalLevel_eq_fixingSubgroup` により、これは
`Gal(K̄/K_{f,n})`(Lubin-Tate 塔の第 `n` 段の固定化部分群)に等しい。 -/
noncomputable def absGalPrincipalLevel (n : ℕ) : Subgroup K.absGal :=
  Subgroup.comap (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf) (principalUnits K π n)

/-- ★★**`Art^{-1}(U^{m+1}_K) = Gal(K̄/K(x_m))`** —— 在庫
`mem_principalUnits_reciprocityUnits_iff` の言い換え。★ここで添字の対応が決まる:
`psiGenSeq ... m` は `ψ_{m+1}` の根なので `K(x_m) = K_{f,m+1}` は塔の第 `m+1` 段であり、
レベル `m+1` と段 `m+1` は**ずれない**(冒頭「添字のずれ」の節)。 -/
theorem absGalPrincipalLevel_eq_fixingSubgroup (m : ℕ) :
    absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1) =
      (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)).fixingSubgroup := by
  ext σ
  simp only [absGalPrincipalLevel, Subgroup.mem_comap,
    mem_principalUnits_reciprocityUnits_iff, mem_fixingSubgroup_adjoin_iff,
    Set.mem_singleton_iff, forall_eq]

/-- 可換群への準同型の引き戻しなので**正規**。 -/
instance normal_absGalPrincipalLevel (n : ℕ) :
    (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n).Normal := by
  unfold absGalPrincipalLevel
  infer_instance

/-- **開**部分群。`n = 0` では `⊤`、`n = m+1` では有限次拡大 `K(x_m)/K` の
固定化部分群(`IntermediateField.fixingSubgroup_isOpen`)。 -/
theorem isOpen_absGalPrincipalLevel (n : ℕ) :
    IsOpen ((absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n : Subgroup K.absGal) :
      Set K.absGal) := by
  cases n with
  | zero =>
    have h : absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf 0 = ⊤ := by
      ext σ
      simp [absGalPrincipalLevel, principalUnits_zero_eq_top K]
    rw [h]; simp
  | succ m =>
    rw [absGalPrincipalLevel_eq_fixingSubgroup]
    haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
    exact IntermediateField.fixingSubgroup_isOpen _

/-- ★したがって射影系の底(開正規部分群の全体)に属する。
★**主定理の仮定 1 は、`F.base = openNormalBase` なら自動で満たされる。** -/
theorem absGalPrincipalLevel_mem_openNormalBase (n : ℕ) :
    absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n ∈ openNormalBase K.absGal :=
  ⟨isOpen_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n,
    normal_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n⟩

theorem absGalPrincipalLevel_antitone {i m : ℕ} (h : i ≤ m) :
    absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf m
      ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf i :=
  Subgroup.comap_mono (antitone_principalUnits K π h)

theorem directed_absGalPrincipalLevel (n : ℕ) :
    Directed (· ≥ ·) (fun i : ℕ => absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i)) :=
  fun i j => ⟨max i j,
    absGalPrincipalLevel_antitone K hq hπmax hπne0 f hf0 hf1 hf
      (Nat.add_le_add_left (le_max_left i j) n),
    absGalPrincipalLevel_antitone K hq hπmax hπne0 f hf0 hf1 hf
      (Nat.add_le_add_left (le_max_right i j) n)⟩

/-- ★★**核の同定**: `⋂_{i} Art^{-1}(U^{n+i}_K) = ker(Art)`。
これは `⋂_n U^n_K = 1`(`iInf_principalUnits_eq_bot`)の引き戻しであり、
★**原典より短い道**の要点である(Herbrand 関数を使わずに済む)。 -/
theorem iInf_absGalPrincipalLevel_eq_ker (n : ℕ) :
    (⨅ i : ℕ, absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i))
      = (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf).ker := by
  refine le_antisymm ?_ ?_
  · intro σ hσ
    rw [Subgroup.mem_iInf] at hσ
    have hall : reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ
        ∈ (⨅ m : ℕ, principalUnits K π m) := by
      refine Subgroup.mem_iInf.2 fun m => ?_
      exact antitone_principalUnits K π (Nat.le_add_left m n) (hσ m)
    rw [iInf_principalUnits_eq_bot K hπmax] at hall
    exact MonoidHom.mem_ker.2 (Subgroup.mem_bot.1 hall)
  · intro σ hσ
    refine Subgroup.mem_iInf.2 fun i => ?_
    show reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ ∈ principalUnits K π (n + i)
    rw [MonoidHom.mem_ker.1 hσ]
    exact one_mem _

/-! ## §7 主定理 -/

/-- ★★**包含 `⊆`**: `Γ_K^v` の像は `U^n_K` に含まれる。

仮定 `hstage` は「**第 `n` 段での上付き分岐群が消える**」——
`Gal(K_{f,n}/K)^v = {id}`、すなわち `S N v = N`(`le_S` と合わせて等号)——の形であり、
これは Yoshida Proposition 6.14(木の
`upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap`)そのものである。

★コンパクト性も相互律の全射性も使わない。★`v` は実数のままでよい。 -/
theorem map_limit_reciprocityUnits_le_principalUnits (F : StageFiltration K.absGal)
    (n : ℕ) (hN : absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n ∈ F.base) {v : ℝ}
    (hstage : F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) v
      ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) :
    (F.limit v).map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
      ≤ principalUnits K π n :=
  F.map_limit_le_of_stage_le_comap _ hN hstage

/-- 同じ仮定から `Γ_K^v ⊆ Art^{-1}(U^n_K)` そのもの。 -/
theorem limit_le_absGalPrincipalLevel (F : StageFiltration K.absGal) (n : ℕ)
    (hN : absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n ∈ F.base) {v : ℝ}
    (hstage : F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) v
      ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) :
    F.limit v ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n :=
  (F.limit_le hN v).trans hstage

/-- ★体の言葉での言い換え: `Γ_K^v` は `K_{f,m+1} = K(x_m)` を**各点固定する**。 -/
theorem limit_le_fixingSubgroup (F : StageFiltration K.absGal) (m : ℕ)
    (hN : absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1) ∈ F.base) {v : ℝ}
    (hstage : F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1)) v
      ≤ absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (m + 1)) :
    F.limit v ≤ (IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)).fixingSubgroup := by
  rw [← absGalPrincipalLevel_eq_fixingSubgroup]
  exact limit_le_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf F (m + 1) hN hstage

/-- ★★★★**主定理(等号)** —— pGC p.4 の "the image of `Γ^v_K` in `Γ^ab_K` is equal to
`U^v_K`" の、`v = n : ℕ` の場合。

仮定 `hstage` は原典の "well-known" の**有限段での中身**
`Gal(K_{f,n+i}/K)^n = U^n_K/U^{n+i}_K` である(★`i` すべてで要る)。
`hbase` は `F.base = openNormalBase K.absGal` なら
`absGalPrincipalLevel_mem_openNormalBase` で自動。

★★仮定が空虚でないことは §8 で構成的に示す。 -/
theorem map_limit_reciprocityUnits_eq_principalUnits (F : StageFiltration K.absGal)
    (n : ℕ) {v : ℝ}
    (hbase : ∀ i : ℕ, absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i) ∈ F.base)
    (hstage : ∀ i : ℕ, F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i)) v
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) :
    (F.limit v).map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
      = principalUnits K π n := by
  haveI := compactSpace_absGal K
  exact F.map_limit_eq_of_stage_eq_comap
    (reciprocityUnits_surjective K hq hπmax hπne0 f hf0 hf1 hf) hbase
    (directed_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n)
    (iInf_absGalPrincipalLevel_eq_ker K hq hπmax hπne0 f hf0 hf1 hf n) hstage

/-- ★★★★★**主定理を原典の記号で**: `Art(Γ_K^v) = U^v_K`(`v = n : ℕ`)。

原文 (pGC p.4):
> Then it is well-known (Theorem 1 of [3], p. 155) that the image of Γ^v_K in Γ^ab_K is equal to U^v_K ⊆ U_K.

★逸脱: `Γ^ab_K` ではなく `𝒪_K^×` 成分で述べている(冒頭「逸脱の記録 1」)。
★逸脱: `v` は自然数(冒頭「逸脱の記録 3」)。 -/
theorem map_limit_reciprocityUnits_eq_upperPrincipalUnits (F : StageFiltration K.absGal)
    (n : ℕ)
    (hbase : ∀ i : ℕ, absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i) ∈ F.base)
    (hstage : ∀ i : ℕ, F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i)) (n : ℝ)
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) :
    (F.limit (n : ℝ)).map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
      = upperPrincipalUnits K π (n : ℝ) := by
  rw [upperPrincipalUnits_natCast]
  exact map_limit_reciprocityUnits_eq_principalUnits K hq hπmax hπne0 f hf0 hf1 hf F n hbase hstage

/-! ## §8 仮定の非空虚性(★退化 witness。原典の `Γ_K^v` ではない) -/

/-- **退化 witness**: `S N v := Art^{-1}(U^{⌈v⌉₊}_K) ⊔ N`。

★★これは §7 の仮定 2 本を**同時に満たす**段データであり、主定理が空虚でないことを示す。
★★**原典の `Γ_K^v` ではない**(この `F` の極限は `Art^{-1}(U^{⌈v⌉₊}_K)` そのもので、
本物の `Γ_K^v` はそれよりずっと小さい)。Y19 の `trivialStageFiltration` と同じ扱いをすること
——G2 の非空虚 witness として使ってはならない。 -/
noncomputable def lubinTateStageFiltration : StageFiltration K.absGal :=
  StageFiltration.ofAntitoneNormal
    (fun v => absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf ⌈v⌉₊)
    (fun v => normal_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf ⌈v⌉₊)
    (fun _ _ hab => absGalPrincipalLevel_antitone K hq hπmax hπne0 f hf0 hf1 hf
      (Nat.ceil_mono hab))

theorem lubinTateStageFiltration_base :
    (lubinTateStageFiltration K hq hπmax hπne0 f hf0 hf1 hf).base = openNormalBase K.absGal := rfl

theorem lubinTateStageFiltration_S (N : Subgroup K.absGal) (v : ℝ) :
    (lubinTateStageFiltration K hq hπmax hπne0 f hf0 hf1 hf).S N v
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf ⌈v⌉₊ ⊔ N := rfl

theorem lubinTateStageFiltration_limit (v : ℝ) :
    (lubinTateStageFiltration K hq hπmax hπne0 f hf0 hf1 hf).limit v
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf ⌈v⌉₊ :=
  StageFiltration.ofAntitoneNormal_limit _ _ _
    (fun _ hU => absGal_exists_mem_openNormalBase_subset K hU)
    (fun w => Subgroup.isClosed_of_isOpen _
      (isOpen_absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf ⌈w⌉₊)) v

theorem lubinTateStageFiltration_stage (n i : ℕ) :
    (lubinTateStageFiltration K hq hπmax hπne0 f hf0 hf1 hf).S
        (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i)) (n : ℝ)
      = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n := by
  rw [lubinTateStageFiltration_S, Nat.ceil_natCast]
  exact sup_eq_left.2 (absGalPrincipalLevel_antitone K hq hπmax hπne0 f hf0 hf1 hf
    (Nat.le_add_right n i))

/-- ★★★**主定理の仮定は空虚でない**: §7 の 2 本の仮定を同時に満たす段データが存在し、
そのとき結論 `Art(Γ^v) = U^v_K` が実際に成り立つ。 -/
theorem exists_stageFiltration_map_limit_eq_upperPrincipalUnits (n : ℕ) :
    ∃ F : StageFiltration K.absGal,
      (∀ i : ℕ, absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i) ∈ F.base) ∧
      (∀ i : ℕ, F.S (absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf (n + i)) (n : ℝ)
        = absGalPrincipalLevel K hq hπmax hπne0 f hf0 hf1 hf n) ∧
      (F.limit (n : ℝ)).map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
        = upperPrincipalUnits K π (n : ℝ) := by
  refine ⟨lubinTateStageFiltration K hq hπmax hπne0 f hf0 hf1 hf,
    fun i => absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf (n + i),
    fun i => lubinTateStageFiltration_stage K hq hπmax hπne0 f hf0 hf1 hf n i, ?_⟩
  exact map_limit_reciprocityUnits_eq_upperPrincipalUnits K hq hπmax hπne0 f hf0 hf1 hf _ n
    (fun i => absGalPrincipalLevel_mem_openNormalBase K hq hπmax hπne0 f hf0 hf1 hf (n + i))
    (fun i => lubinTateStageFiltration_stage K hq hπmax hπne0 f hf0 hf1 hf n i)

end Concrete

end ABC3.Found.PGC
