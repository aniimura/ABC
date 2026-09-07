import ABC3.Found.PGC.UpperRamificationGroup
import ABC3.Found.PGC.FilteredGroup
import ABC3.Found.PGC.LubinTateReciprocityCompactness
import ABC3.Interface.PGC.LocalFieldData
import Mathlib.Topology.Algebra.ClopenNhdofOne

/-!
# 絶対 Galois 群の上の分岐フィルトレーション —— 有限段の族から逆極限を作る

典拠: S. Mochizuki, *A Version of the Grothendieck Conjecture for p-adic Local Fields*
(pGC) の **Definition 2.3**(物理 p.5)と、そこに入力として現れる高次分岐群の族
`{Γ_K^v}`(物理 p.4)。有限段の上付き番号付けは T. Yoshida,
*Local Class Field Theory via Lubin-Tate Theory* の **Definition 6.12 / Corollary 6.13(i)**
(物理 p.17)を使う(`Found/PGC/UpperRamificationGroup.lean`)。

原文 (pGC p.4):
> Then we shall denote by Γ^v_K ⊆ Γ_K the higher ramification group associated to the
> number v in the "upper numbering" (see, e.g., [3], p. 155).

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0.

## 本ファイルが埋めるもの / 埋めないもの(★正直な線引き)

`Interface/PGC/LocalFieldData.lean` の `RamificationFiltration p` は
`Gv` / `isClosed` / `isNormal` / `antitone` の 4 つ組を要求する。本ファイルは

* **抽象核**(§1–§3、位相群の言葉だけ。分岐・付値・Galois の語彙が 1 つも出てこない):
  「開正規部分群の下降有向族 `base` と、各段の(引き戻し形の)フィルトレーション `S N v`
  が **両立している**(`S N v · M = S M v`)」というデータから、極限
  `Γ^v := ⋂_{N ∈ base} S N v` を作り、**閉・正規・反単調**を出す。さらに
  - `Γ^v · M = S M v`(★**極限が各段へ全射する**。`CompactSpace` を使う)、
  - その性質を持つ**閉**部分群は `Γ^v` に限る(★一意性。`base` が 1 の近傍基のとき)
  を証明する。
* **具体層**(§6): `Γ := K.absGal` へ代入して `RamificationFiltration p` の
  **4 フィールドすべて**を、各 `K` ごとの段データ `StageFiltration K.absGal` から作る
  (`ramificationFiltrationOfStages`)。`K.absGal` が
  コンパクト(`compactSpace_algEquiv`)・完全不連結であることから、
  §1 の全射性・一意性がそのまま具体層で使える。

**埋めていないもの(★次のノード)**: 各 `K` に対する段データ
`StageFiltration K.absGal` の**中身**、すなわち有限 Galois 部分拡大 `L/K` ごとの
`Gal(L/K)^v` の Γ_K への引き戻し。これは Y18 の `upperRamificationGroup` を
`Gal(L/K)` に当てる仕事であり、下の「完全分岐の穴」を通る必要がある。

## ★★完全分岐の穴(測定して選んだ道と、その理由)

Yoshida §6.1 の設定(したがって `upperRamificationGroup` / `herbrandPhiGroup`)は
`Algebra.adjoin A {α} = ⊤`(`α` は一意化元)を要求する——すなわち**完全分岐**である。
一方 `RamificationFiltration` は**任意の**有限 Galois 部分拡大にわたる逆極限を要求し、
一般の有限 Galois 拡大は完全分岐ではない。3 つの道を測った:

* **道 C(完全分岐な有限 Galois 部分拡大だけで逆極限を張る)は使えない。**
  その族は**有向でない**——完全分岐拡大の合成は完全分岐とは限らない。
  反例: `K = ℚ_3` で `K(√3)` と `K(√-3)` はどちらも完全分岐(次数 2 で分岐)だが、
  合成体は `√3 · √-3 = 3√-1` を含むので `√-1` を含み、`ℚ_3(√-1)` は
  **不分岐**二次拡大である(`-1` は `mod 3` で平方非剰余)。ゆえに合成体は
  完全分岐でない。共終どころか有向ですらないので、逆極限の底として使えない。
* **道 A / B(惰性部分に寄せる)が正しい道である。** `v > 0` では
  `Γ_K^v ⊆ I_K` であり、有限段では `L₀ := L ∩ K^ur` と置くと `L/L₀` は完全分岐で、
  `Gal(L/K)^v = Gal(L/L₀)^v`(`v > 0`)。Yoshida の設定はこの `L/L₀` にそのまま当たる。
  ★この道で残る数学は 2 つ: (a) 下付き分岐群 `G_i`(`i ≥ 0`)が
  **不分岐底変換 `L/K ↝ L/L₀` で変わらない**こと(付値が同じだから、が理由)、
  (b) `L ⊆ L′` のとき惰性群の間の写像 `Gal(L′/L′₀) ↠ Gal(L/L₀)` が全射で、
  その核による商が Yoshida Corollary 6.13(i) の `G/H` に一致すること。
  ★★どちらも本ファイルの抽象核の**外側**にあり、新しいノードである。
* 本ファイルは**どちらの道でも共通に効く部分**(逆極限の構成と、その 4 性質、
  全射性、一意性、そして Y18 の結論を Γ レベルへ移す `comap_compat_of_coe_mul_coe_eq`)を
  切り出した。★**穴そのもの(有限段の `S N v` の構成)は埋まっていない**——
  それは「不分岐部分の除去 + 惰性群の全射性」という別の数学であり、新ノードとして報告する。

## ★Y18 の `upperRamification_coe_mul_coe_eq` はどう入るか

Y18 の結論は**集合の積の等式**
`↑(G^m) * ↑H = ↑((G/H)^m の引き戻し)`(`G` は有限段の Galois 群)である。
本ファイルの `StageFiltration.compat` は**まさにこの形**をしている。
有限段 `G = Γ_K/N` から `Γ_K` へ移すのは §4 の
`comap_compat_of_coe_mul_coe_eq`(全射準同型による引き戻しは集合の積を保つ)で、
これは純群論である。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 原文の `v` の範囲は `v > 0` だが、`Interface` の `RamificationFiltration` に合わせて
   **全実数**で定義する(`Antitone` を全実数上で書くための設計判断。
   `Interface/PGC/LocalFieldData.lean` の同じ逸脱を引き継いだだけで、新しい逸脱ではない)。
2. 段データの添字を「開正規部分群 `N`」で取り、`S N v` を**商群 `Γ/N` の部分群ではなく
   `Γ` への引き戻し**として扱う。これは Y18 が Corollary 6.13(i) を
   `G^m H/H` ではなく `G^m · H`(`G` の中の集合)で書いているのに合わせたもので、
   商群の型を作らないぶん配管が軽い。商群の言葉での形は `map_limit_eq` に置いた。
3. 段データ `S` は `base` の外でも値を持つ**全域関数**にしてある(条件は `base` の上だけ)。
   `base` の外の値は `limit` にも定理にも一切効かない。

## 退化の自己検査

* **両立性(`compat`)を落とすと `coe_limit_mul_coe` が壊れる**——極限が各段へ全射しなくなり、
  「`Γ^v` の `Gal(L/K)` への像が `Gal(L/K)^v`」という、逆極限で定義する意味そのものが消える。
  ★本ファイルはその全射性を**定理として**出しているので、この退化は自動的に排除される。
* **`isClosed` は落とせない**(構造体が要求している)。★測った結果:
  各段 `S N v` は開部分群 `N` を含むので**開**、位相群の開部分群は**閉**
  (`Subgroup.isClosed_of_isOpen`)。その交わりだから閉——見立てのとおりだった。
* **`antitone` の向き**: `Antitone` の定義 `a ≤ b → f b ≤ f a` に `a := v2, b := v1` を
  代入すると原文の「`v1 ≥ v2 ⟹ Γ^{v1} ⊆ Γ^{v2}`」になる。★上付きは番号が大きいほど小さい。
* **`v ≤ 0` の扱い**: 抽象核は `v ≤ 0` に何も課さない(段データが決める)。
  Y18 は `m ≤ 0` で `G^m = ⊤` としており、その段データを入れれば極限も `v ≤ 0` で `⊤` になる。
* §8 の `trivialStageFiltration` は**退化した witness** である(`v > 0` で `Γ^v = ⊥`)。
  構造体が空虚でないことしか言っていない。★**原典の `Γ_K^v` ではない**ので、
  G2 の非空虚 witness として使ってはならない。
-/

namespace ABC3.Found.PGC

open scoped Pointwise

/-! ## §1 抽象核 —— 両立する段フィルトレーション -/

/-- **抽象核**: 位相群 `Γ` の上の「射影系の各段のフィルトレーション」。

`base` は開正規部分群の下降有向族(`Γ` の有限商の射影系の添字)であり、
`S N v` は第 `N` 段の商 `Γ/N` 上のフィルトレーションの `Γ` への**引き戻し**である。
`compat` が Herbrand の定理(Yoshida Corollary 6.13(i))に対応する。

★分岐・付値・Galois の語彙は 1 つも出てこない。 -/
structure StageFiltration (Γ : Type*) [Group Γ] [TopologicalSpace Γ] where
  /-- 射影系の添字となる開正規部分群の族 -/
  base : Set (Subgroup Γ)
  /-- `base` の元は開部分群 -/
  isOpen_base : ∀ ⦃N⦄, N ∈ base → IsOpen (N : Set Γ)
  /-- `base` の元は正規部分群 -/
  normal_base : ∀ ⦃N⦄, N ∈ base → N.Normal
  /-- `base` は下降有向(射影系であること) -/
  directed_base : ∀ ⦃M⦄, M ∈ base → ∀ ⦃N⦄, N ∈ base → ∃ P ∈ base, P ≤ M ⊓ N
  /-- 第 `N` 段のフィルトレーションの `Γ` への引き戻し -/
  S : Subgroup Γ → ℝ → Subgroup Γ
  /-- 引き戻しなので `N` を含む -/
  le_S : ∀ ⦃N⦄, N ∈ base → ∀ v, N ≤ S N v
  /-- 各段は正規 -/
  normal_S : ∀ ⦃N⦄, N ∈ base → ∀ v, (S N v).Normal
  /-- 各段は反単調 -/
  antitone_S : ∀ ⦃N⦄, N ∈ base → Antitone (S N)
  /-- ★★**両立性**(Yoshida Corollary 6.13(i) の形): `N ≤ M` のとき `S N v · M = S M v`。
  商の言葉では「`(Γ/N)^v` の `Γ/M` への像が `(Γ/M)^v`」。 -/
  compat : ∀ ⦃M⦄, M ∈ base → ∀ ⦃N⦄, N ∈ base → N ≤ M → ∀ v,
    (S N v : Set Γ) * (M : Set Γ) = (S M v : Set Γ)

namespace StageFiltration

variable {Γ : Type*} [Group Γ] [TopologicalSpace Γ] (F : StageFiltration Γ)

/-- **極限フィルトレーション** `Γ^v := ⋂_{N ∈ base} S N v`。 -/
def limit (v : ℝ) : Subgroup Γ := ⨅ N : F.base, F.S (N : Subgroup Γ) v

theorem limit_le {N : Subgroup Γ} (hN : N ∈ F.base) (v : ℝ) : F.limit v ≤ F.S N v :=
  iInf_le (fun N : F.base => F.S (N : Subgroup Γ) v) ⟨N, hN⟩

theorem mem_limit_iff {v : ℝ} {x : Γ} : x ∈ F.limit v ↔ ∀ N ∈ F.base, x ∈ F.S N v := by
  simp only [limit, Subgroup.mem_iInf, Subtype.forall]

/-- 段は添字の部分群について単調(`compat` から出る。`1 ∈ M` を使うだけ)。 -/
theorem S_mono {M N : Subgroup Γ} (hM : M ∈ F.base) (hN : N ∈ F.base) (hNM : N ≤ M) (v : ℝ) :
    F.S N v ≤ F.S M v := by
  intro x hx
  have h : x ∈ (F.S N v : Set Γ) * (M : Set Γ) := ⟨x, hx, 1, one_mem _, mul_one x⟩
  rw [F.compat hM hN hNM v] at h
  exact h

/-- 各段は**開**部分群(開部分群 `N` を含むから)。 -/
theorem isOpen_S [IsTopologicalGroup Γ] {N : Subgroup Γ} (hN : N ∈ F.base) (v : ℝ) :
    IsOpen (F.S N v : Set Γ) :=
  Subgroup.isOpen_mono (F.le_S hN v) (F.isOpen_base hN)

/-- 各段は**閉**部分群(位相群の開部分群は閉)。 -/
theorem isClosed_S [IsTopologicalGroup Γ] {N : Subgroup Γ} (hN : N ∈ F.base) (v : ℝ) :
    IsClosed (F.S N v : Set Γ) :=
  Subgroup.isClosed_of_isOpen _ (F.isOpen_S hN v)

/-- ★`RamificationFiltration.isClosed` に対応: 極限は閉部分群。 -/
theorem isClosed_limit [IsTopologicalGroup Γ] (v : ℝ) : IsClosed (F.limit v : Set Γ) := by
  rw [limit, Subgroup.coe_iInf]
  exact isClosed_iInter fun N => F.isClosed_S N.2 v

/-- ★`RamificationFiltration.isNormal` に対応: 極限は正規部分群。 -/
theorem normal_limit (v : ℝ) : (F.limit v).Normal :=
  Subgroup.normal_iInf_normal fun N => F.normal_S N.2 v

/-- ★`RamificationFiltration.antitone` に対応: 極限は反単調
(`v1 ≥ v2 → Γ^{v1} ⊆ Γ^{v2}`)。 -/
theorem antitone_limit : Antitone F.limit := by
  intro a b hab
  simp only [limit, le_iInf_iff]
  exact fun N => le_trans (iInf_le _ N) (F.antitone_S N.2 hab)

/-! ## §2 ★★★極限は各段へ全射する(両立性が効くところ) -/

/-- ★★★**極限フィルトレーションは各段へ全射する**: `Γ^v · M = S M v`。

これが「有限段の族の逆極限として `Γ^v` を定義してよい」ことの中身である。
`⊆` は各段が部分群であることから、`⊇` は**コンパクト性**(下降有向な空でない閉集合族の
交わりは空でない)から出る。

★★両立性 `compat` を落とすとこの定理は壊れる(`S N v` たちが噛み合わず、
`x ∈ S M v` に対して `xM` と `S N v` の交わりが空になりうる)。 -/
theorem coe_limit_mul_coe [IsTopologicalGroup Γ] [CompactSpace Γ] {M : Subgroup Γ}
    (hM : M ∈ F.base) (v : ℝ) :
    (F.limit v : Set Γ) * (M : Set Γ) = (F.S M v : Set Γ) := by
  refine Set.Subset.antisymm ?_ ?_
  · rintro _ ⟨a, ha, b, hb, rfl⟩
    exact mul_mem (F.limit_le hM v ha) (F.le_S hM v hb)
  · intro x hx
    have hclM : IsClosed (x • (M : Set Γ)) :=
      (Subgroup.isClosed_of_isOpen _ (F.isOpen_base hM)).smul x
    haveI : Nonempty {N : Subgroup Γ // N ∈ F.base ∧ N ≤ M} := ⟨⟨M, hM, le_rfl⟩⟩
    set t : {N : Subgroup Γ // N ∈ F.base ∧ N ≤ M} → Set Γ :=
      fun N => (x • (M : Set Γ)) ∩ (F.S N.1 v : Set Γ) with ht
    have htn : ∀ N, (t N).Nonempty := by
      intro N
      have hxm : x ∈ (F.S N.1 v : Set Γ) * (M : Set Γ) := by
        rw [F.compat hM N.2.1 N.2.2 v]; exact hx
      obtain ⟨s, hs, m, hm, hsm⟩ := hxm
      refine ⟨s, ?_, hs⟩
      rw [Set.mem_smul_set_iff_inv_smul_mem, smul_eq_mul]
      have hxs : x⁻¹ * s = m⁻¹ := by rw [← hsm]; group
      rw [hxs]
      exact inv_mem hm
    have htcl : ∀ N, IsClosed (t N) := fun N => hclM.inter (F.isClosed_S N.2.1 v)
    have htd : Directed (· ⊇ ·) t := by
      rintro ⟨N₁, hN₁, hN₁M⟩ ⟨N₂, hN₂, hN₂M⟩
      obtain ⟨P, hP, hPle⟩ := F.directed_base hN₁ hN₂
      refine ⟨⟨P, hP, le_trans hPle (le_trans inf_le_left hN₁M)⟩, ?_, ?_⟩
      · exact Set.inter_subset_inter_right _
          (F.S_mono hN₁ hP (le_trans hPle inf_le_left) v)
      · exact Set.inter_subset_inter_right _
          (F.S_mono hN₂ hP (le_trans hPle inf_le_right) v)
    obtain ⟨y, hy⟩ := IsCompact.nonempty_iInter_of_directed_nonempty_isCompact_isClosed
      t htd htn (fun N => (htcl N).isCompact) htcl
    simp only [Set.mem_iInter, ht, Set.mem_inter_iff] at hy
    have hyM : y ∈ x • (M : Set Γ) := (hy ⟨M, hM, le_rfl⟩).1
    have hylim : y ∈ F.limit v := by
      rw [F.mem_limit_iff]
      intro N hN
      obtain ⟨P, hP, hPle⟩ := F.directed_base hN hM
      exact F.S_mono hN hP (le_trans hPle inf_le_left) v
        (hy ⟨P, hP, le_trans hPle inf_le_right⟩).2
    rw [Set.mem_smul_set_iff_inv_smul_mem, smul_eq_mul] at hyM
    exact ⟨y, hylim, y⁻¹ * x, by simpa [mul_inv_rev] using inv_mem hyM, by group⟩

/-- 上の全射性を部分群の束の言葉で書いた形: `Γ^v ⊔ M = S M v`。 -/
theorem limit_sup_eq [IsTopologicalGroup Γ] [CompactSpace Γ] {M : Subgroup Γ}
    (hM : M ∈ F.base) (v : ℝ) : F.limit v ⊔ M = F.S M v := by
  haveI := F.normal_base hM
  apply SetLike.coe_injective
  rw [Subgroup.mul_normal]
  exact F.coe_limit_mul_coe hM v

/-- ★**商群の言葉での全射性**: `Γ^v` の `Γ/M` への像は `S M v` の像
(すなわち第 `M` 段のフィルトレーション `(Γ/M)^v`)にちょうど一致する。

★これが「`Γ^v` を有限段の逆極限として定義してよい」ことの、原文に一番近い言い方である。 -/
theorem map_limit_eq [IsTopologicalGroup Γ] [CompactSpace Γ] {M : Subgroup Γ} [M.Normal]
    (hM : M ∈ F.base) (v : ℝ) :
    (F.limit v).map (QuotientGroup.mk' M) = (F.S M v).map (QuotientGroup.mk' M) := by
  have hb : M.map (QuotientGroup.mk' M) = ⊥ := by
    rw [Subgroup.map_eq_bot_iff, QuotientGroup.ker_mk']
  rw [← F.limit_sup_eq hM v, Subgroup.map_sup, hb, sup_bot_eq]

/-! ## §3 一意性 —— この定義以外にありえない -/

/-- ★★**一意性**: `base` が `1` の近傍基であるとき、
「各段へ全射する**閉**部分群」は極限フィルトレーションに限る。

★コンパクト性は使わない(存在の側だけがコンパクト性を要る)。
★`T` が閉であることは落とせない: `T` の閉包は `T` と同じ段を持つ。 -/
theorem eq_limit_of_isClosed [IsTopologicalGroup Γ]
    (hbasis : ∀ U ∈ nhds (1 : Γ), ∃ N ∈ F.base, (N : Set Γ) ⊆ U)
    {T : Subgroup Γ} (hTcl : IsClosed (T : Set Γ)) {v : ℝ}
    (hT : ∀ M ∈ F.base, (T : Set Γ) * (M : Set Γ) = (F.S M v : Set Γ)) :
    T = F.limit v := by
  refine le_antisymm ?_ ?_
  · intro x hx
    rw [F.mem_limit_iff]
    intro N hN
    have h : x ∈ (T : Set Γ) * (N : Set Γ) := ⟨x, hx, 1, one_mem _, mul_one x⟩
    rw [hT N hN] at h
    exact h
  · intro x hx
    have hmem : x ∈ (T : Set Γ) := by
      rw [← hTcl.closure_eq, mem_closure_iff_nhds]
      intro U hU
      have hcont : ContinuousAt (fun g : Γ => x * g) 1 :=
        (continuous_const.mul continuous_id).continuousAt
      have hU1 : (fun g : Γ => x * g) ⁻¹' U ∈ nhds (1 : Γ) := by
        apply hcont.preimage_mem_nhds
        simpa using hU
      obtain ⟨N, hN, hNU⟩ := hbasis _ hU1
      have hxS : x ∈ (T : Set Γ) * (N : Set Γ) := by
        rw [hT N hN]; exact F.limit_le hN v hx
      obtain ⟨t, htT, n, hn, htn⟩ := hxS
      refine ⟨t, ?_, htT⟩
      have hteq : t = x * n⁻¹ := by rw [← htn]; group
      rw [hteq]
      exact hNU (inv_mem hn)
    exact hmem

/-- ★**pGC Definition 2.3 への橋**: 極限フィルトレーションは
`FilteredGroup`(原文の filtered group)そのものである。 -/
def toFilteredGroup [IsTopologicalGroup Γ] : FilteredGroup where
  G := Γ
  Gv := F.limit
  isClosed := F.isClosed_limit
  isNormal := F.normal_limit
  antitone := F.antitone_limit

@[simp] theorem toFilteredGroup_Gv [IsTopologicalGroup Γ] (v : ℝ) :
    F.toFilteredGroup.Gv v = F.limit v := rfl

end StageFiltration

/-! ## §4 Y18 の結論を `Γ` へ移す(純群論) -/

/-- ★★**全射準同型による引き戻しは「集合の積の等式」を保つ**。

Yoshida Corollary 6.13(i)(`Found/PGC/UpperRamificationGroup.lean` の
`upperRamification_coe_mul_coe_eq`)は有限段の Galois 群 `G` の中の等式
`↑(G^m) * ↑H = ↑((G/H)^m)` である。`f : Γ ↠ G` を全射準同型(`Γ_K → Gal(L/K)`)とすると、
本補題によりそれがそのまま `StageFiltration.compat` の形になる。

★分岐も付値も出てこない——`f` の全射性と `1 ∈ B` だけを使う。 -/
theorem comap_compat_of_coe_mul_coe_eq {Γ G : Type*} [Group Γ] [Group G] (f : Γ →* G)
    (hf : Function.Surjective f) {A B C : Subgroup G}
    (h : (A : Set G) * (B : Set G) = (C : Set G)) :
    ((A.comap f : Subgroup Γ) : Set Γ) * ((B.comap f : Subgroup Γ) : Set Γ)
      = ((C.comap f : Subgroup Γ) : Set Γ) := by
  ext x
  constructor
  · rintro ⟨a, ha, b, hb, rfl⟩
    show f (a * b) ∈ C
    rw [map_mul, ← SetLike.mem_coe, ← h]
    exact ⟨f a, ha, f b, hb, rfl⟩
  · intro hx
    have hfx : f x ∈ (A : Set G) * (B : Set G) := by rw [h]; exact hx
    obtain ⟨a, ha, b, hb, hab⟩ := hfx
    obtain ⟨α, hα⟩ := hf a
    refine ⟨α, ?_, α⁻¹ * x, ?_, by group⟩
    · show f α ∈ A; rw [hα]; exact ha
    · show f (α⁻¹ * x) ∈ B
      rw [map_mul, map_inv, hα, ← hab]
      simpa using hb

/-! ## §5 開正規部分群の族 —— 射影系の底 -/

/-- 開正規部分群の全体。副有限群ではこれが射影系の底になる。 -/
def openNormalBase (Γ : Type*) [Group Γ] [TopologicalSpace Γ] : Set (Subgroup Γ) :=
  {N : Subgroup Γ | IsOpen (N : Set Γ) ∧ N.Normal}

theorem isOpen_of_mem_openNormalBase {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    {N : Subgroup Γ} (hN : N ∈ openNormalBase Γ) : IsOpen (N : Set Γ) := hN.1

theorem normal_of_mem_openNormalBase {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    {N : Subgroup Γ} (hN : N ∈ openNormalBase Γ) : N.Normal := hN.2

/-- 開正規部分群の族は**下降有向**(共通部分を取ればよい)。★射影系であることの中身。 -/
theorem directed_openNormalBase {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    {M N : Subgroup Γ} (hM : M ∈ openNormalBase Γ) (hN : N ∈ openNormalBase Γ) :
    ∃ P ∈ openNormalBase Γ, P ≤ M ⊓ N := by
  haveI := hM.2; haveI := hN.2
  refine ⟨M ⊓ N, ⟨?_, inferInstance⟩, le_rfl⟩
  rw [Subgroup.coe_inf]
  exact hM.1.inter hN.1

/-- 副有限群(コンパクト + 完全不連結な位相群)では、開正規部分群は `1` の近傍基をなす。
mathlib の `ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one` の言い換え。 -/
theorem exists_mem_openNormalBase_subset {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] [CompactSpace Γ] [TotallyDisconnectedSpace Γ]
    {U : Set Γ} (hU : U ∈ nhds (1 : Γ)) :
    ∃ N ∈ openNormalBase Γ, (N : Set Γ) ⊆ U := by
  obtain ⟨V, hVU, hVopen, hV1⟩ := mem_nhds_iff.1 hU
  obtain ⟨H, hH⟩ := ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one hVopen hV1
  exact ⟨H.toOpenSubgroup.toSubgroup, ⟨H.toOpenSubgroup.isOpen', H.isNormal'⟩, hH.trans hVU⟩

/-- 段データを「開正規部分群の族」の上で作るための構成子
(`base` まわりの 3 条件は自動で埋まる)。 -/
def StageFiltration.ofOpenNormal {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    (S : Subgroup Γ → ℝ → Subgroup Γ)
    (le_S : ∀ ⦃N⦄, N ∈ openNormalBase Γ → ∀ v, N ≤ S N v)
    (normal_S : ∀ ⦃N⦄, N ∈ openNormalBase Γ → ∀ v, (S N v).Normal)
    (antitone_S : ∀ ⦃N⦄, N ∈ openNormalBase Γ → Antitone (S N))
    (compat : ∀ ⦃M⦄, M ∈ openNormalBase Γ → ∀ ⦃N⦄, N ∈ openNormalBase Γ → N ≤ M → ∀ v,
      (S N v : Set Γ) * (M : Set Γ) = (S M v : Set Γ)) :
    StageFiltration Γ where
  base := openNormalBase Γ
  isOpen_base _ hN := hN.1
  normal_base _ hN := hN.2
  directed_base _ hM _ hN := directed_openNormalBase hM hN
  S := S
  le_S := le_S
  normal_S := normal_S
  antitone_S := antitone_S
  compat := compat

@[simp] theorem StageFiltration.ofOpenNormal_base {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    (S : Subgroup Γ → ℝ → Subgroup Γ) (h1 h2 h3 h4) :
    (StageFiltration.ofOpenNormal S h1 h2 h3 h4).base = openNormalBase Γ := rfl

@[simp] theorem StageFiltration.ofOpenNormal_S {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    (S : Subgroup Γ → ℝ → Subgroup Γ) (h1 h2 h3 h4) :
    (StageFiltration.ofOpenNormal S h1 h2 h3 h4).S = S := rfl

/-! ## §6 具体層 —— `RamificationFiltration p` の 4 フィールドを埋める -/

open ABC3.Skeleton.PGC ABC3.Interface.PGC

variable {p : ℕ} [Fact p.Prime]

/-- `Γ_K` はコンパクト(`compactSpace_algEquiv` の言い換え)。 -/
theorem compactSpace_absGal (K : PAdicLocalField p) : CompactSpace K.absGal :=
  compactSpace_algEquiv K

/-- ★`Γ_K` では開正規部分群が `1` の近傍基をなす
(コンパクト + 完全不連結。`Γ_K` は副有限群である)。 -/
theorem absGal_exists_mem_openNormalBase_subset (K : PAdicLocalField p)
    {U : Set K.absGal} (hU : U ∈ nhds (1 : K.absGal)) :
    ∃ N ∈ openNormalBase K.absGal, (N : Set K.absGal) ⊆ U := by
  haveI := compactSpace_absGal K
  exact exists_mem_openNormalBase_subset hU

/-- ★★★**`RamificationFiltration p` の 4 フィールドすべて**を、各 `K` ごとの段データから作る。

`Gv K v` は「有限 Galois 部分拡大 `L/K` ごとの `Gal(L/K)^v` の引き戻し」の共通部分
——すなわち**逆極限**である。

★入力の `F` (各 `K` の段データ)を作るのが残った仕事であり、それが
冒頭「完全分岐の穴」の節の内容である。本定義は**その 1 点に問題を還元した**。 -/
noncomputable def ramificationFiltrationOfStages (F : ∀ K : PAdicLocalField p, StageFiltration K.absGal) :
    RamificationFiltration p where
  Gv K v := (F K).limit v
  isClosed K v := (F K).isClosed_limit v
  isNormal K v := (F K).normal_limit v
  antitone K := (F K).antitone_limit

@[simp] theorem ramificationFiltrationOfStages_Gv
    (F : ∀ K : PAdicLocalField p, StageFiltration K.absGal) (K : PAdicLocalField p) (v : ℝ) :
    (ramificationFiltrationOfStages F).Gv K v = (F K).limit v := rfl

/-- ★★**構成した `Γ_K^v` は各有限段へ全射する**(`Γ_K` のコンパクト性を使う)。

`M` が開正規部分群、`L := (K̄)^M` が対応する有限 Galois 拡大のとき、これは
「`Γ_K^v` の `Gal(L/K)` への像が第 `M` 段のフィルトレーションに一致する」ことである。 -/
theorem ramificationFiltrationOfStages_coe_mul_coe
    (F : ∀ K : PAdicLocalField p, StageFiltration K.absGal) (K : PAdicLocalField p)
    {M : Subgroup K.absGal} (hM : M ∈ (F K).base) (v : ℝ) :
    (((ramificationFiltrationOfStages F).Gv K v : Subgroup K.absGal) : Set K.absGal)
        * (M : Set K.absGal)
      = ((F K).S M v : Set K.absGal) := by
  haveI := compactSpace_absGal K
  exact (F K).coe_limit_mul_coe hM v

/-- ★★**一意性の具体層**: 各段へ全射する閉部分群は `Γ_K^v` に限る
(段データの底が開正規部分群の全体であるとき)。 -/
theorem ramificationFiltrationOfStages_unique
    (F : ∀ K : PAdicLocalField p, StageFiltration K.absGal) (K : PAdicLocalField p)
    (hbase : (F K).base = openNormalBase K.absGal)
    {T : Subgroup K.absGal} (hTcl : IsClosed (T : Set K.absGal)) {v : ℝ}
    (hT : ∀ M ∈ (F K).base, (T : Set K.absGal) * (M : Set K.absGal) = ((F K).S M v : Set K.absGal)) :
    T = (ramificationFiltrationOfStages F).Gv K v := by
  refine (F K).eq_limit_of_isClosed ?_ hTcl hT
  intro U hU
  obtain ⟨N, hN, hNU⟩ := absGal_exists_mem_openNormalBase_subset K hU
  exact ⟨N, by rw [hbase]; exact hN, hNU⟩

/-! ## §7 次のノードへの橋 —— 有限段への全射と、その核

段データ `S N v` を作る側が最初に要るのは
「有限 Galois 部分拡大 `L/K` に対する全射 `Γ_K ↠ Gal(L/K)` と、その核が開正規部分群であること」
である。そこだけ先に置いておく(★ここも純粋な体論で、分岐は出てこない)。

段データを作る手順は:
1. `L` を有限次正規中間体とし、`N := L.fixingSubgroup ∈ openNormalBase K.absGal`
   (`fixingSubgroup_mem_openNormalBase`)、
2. `f := AlgEquiv.restrictNormalHom ↥L : Γ_K ↠ Gal(L/K)`(`surjective_restrictNormalHom`)、
3. `S N v := Subgroup.comap f (Gal(L/K)^v)`(Y18 の `upperRamificationGroup`)、
4. `compat` は Y18 の `upperRamification_coe_mul_coe_eq` に §4 の
   `comap_compat_of_coe_mul_coe_eq` を合わせる。

★4 の段階で「完全分岐の穴」(冒頭)を通る必要がある。 -/

/-- 有限段への制限 `Γ = Gal(E/F) ↠ Gal(L/F)` は**全射**(`L` は正規中間体)。 -/
theorem surjective_restrictNormalHom {F E : Type*} [Field F] [Field E] [Algebra F E]
    [Normal F E] (L : IntermediateField F E) [Normal F L] :
    Function.Surjective (AlgEquiv.restrictNormalHom (F := F) (K₁ := E) (L : Type _)) :=
  AlgEquiv.restrictNormalHom_surjective (F := F) (K₁ := (L : Type _)) E

/-- その**核は `L` の固定化部分群**。 -/
theorem ker_restrictNormalHom_eq_fixingSubgroup {F E : Type*} [Field F] [Field E] [Algebra F E]
    (L : IntermediateField F E) [Normal F L] :
    (AlgEquiv.restrictNormalHom (F := F) (K₁ := E) (L : Type _)).ker = L.fixingSubgroup := by
  ext σ
  rw [MonoidHom.mem_ker, IntermediateField.mem_fixingSubgroup_iff]
  constructor
  · intro h x hx
    have h2 := congrArg (fun (τ : L ≃ₐ[F] L) => (τ ⟨x, hx⟩ : E)) h
    have hc := AlgEquiv.restrictNormal_commutes σ (L : Type _) ⟨x, hx⟩
    simp only [AlgEquiv.restrictNormalHom, MonoidHom.mk'_apply, AlgEquiv.one_apply] at h2
    exact hc.symm.trans h2
  · intro h
    ext x
    simp only [AlgEquiv.restrictNormalHom, MonoidHom.mk'_apply, AlgEquiv.one_apply]
    exact (AlgEquiv.restrictNormal_commutes σ (L : Type _) x).trans (h x x.2)

/-- ★**有限次正規中間体の固定化部分群は射影系の底に属する**
(開: `IntermediateField.fixingSubgroup_isOpen`、正規: 準同型の核)。 -/
theorem fixingSubgroup_mem_openNormalBase (K : PAdicLocalField p)
    (L : IntermediateField K.carrier K.closure) [Normal K.carrier L]
    [FiniteDimensional K.carrier L] :
    L.fixingSubgroup ∈ openNormalBase K.absGal := by
  refine ⟨L.fixingSubgroup_isOpen, ?_⟩
  rw [← ker_restrictNormalHom_eq_fixingSubgroup (E := K.closure) L]
  infer_instance

/-! ## §8 退化 witness(★原典の `Γ_K^v` ではない) -/

/-- **退化した段データ**: `v ≤ 0` で `⊤`、`v > 0` で `N` 自身。

★これは `StageFiltration` の 9 条件が空虚でないことを示すだけのものであり、
**原典の高次分岐群ではない**(副有限群では極限が `v > 0` で `⊥` になる)。
G2 の非空虚 witness として使ってはならない。 -/
noncomputable def trivialStageFiltration (Γ : Type*) [Group Γ] [TopologicalSpace Γ] : StageFiltration Γ :=
  StageFiltration.ofOpenNormal (fun N v => if v ≤ 0 then ⊤ else N)
    (fun _ _ v => by by_cases hv : v ≤ 0 <;> simp [hv])
    (fun _ hN v => by
      by_cases hv : v ≤ 0
      · simp only [hv, if_true]; infer_instance
      · simp only [hv, if_false]; exact hN.2)
    (fun _ _ a b hab => by
      by_cases hb : b ≤ 0
      · have ha : a ≤ 0 := le_trans hab hb
        simp [ha, hb]
      · by_cases ha : a ≤ 0 <;> simp [ha, hb])
    (fun M hM _ _ hNM v => by
      haveI := hM.2
      by_cases hv : v ≤ 0
      · simp only [hv, if_true]
        rw [← Subgroup.mul_normal, top_sup_eq]
      · simp only [hv, if_false]
        rw [← Subgroup.mul_normal, sup_eq_right.2 hNM])

/-- ★退化の確認: 副有限(かつ `T1`)な `Γ` では、上の witness の極限は `v > 0` で `⊥`。
★**原典の `Γ_K^v` はこうならない**(`Γ_K^v` は `v` が小さいところで非自明)。 -/
theorem trivialStageFiltration_limit_eq_bot {Γ : Type*} [Group Γ] [TopologicalSpace Γ]
    [IsTopologicalGroup Γ] [CompactSpace Γ] [TotallyDisconnectedSpace Γ] [T1Space Γ]
    {v : ℝ} (hv : ¬ v ≤ 0) : (trivialStageFiltration Γ).limit v = ⊥ := by
  ext x
  simp only [Subgroup.mem_bot]
  constructor
  · intro hx
    rw [StageFiltration.mem_limit_iff] at hx
    by_contra hne
    have hopen : IsOpen ({x}ᶜ : Set Γ) := isOpen_compl_singleton
    have hmem : ({x}ᶜ : Set Γ) ∈ nhds (1 : Γ) :=
      hopen.mem_nhds (by simpa [eq_comm] using hne)
    obtain ⟨N, hN, hNU⟩ := exists_mem_openNormalBase_subset hmem
    have hxN : x ∈ N := by
      have := hx N hN
      simpa [trivialStageFiltration, hv] using this
    exact hNU hxN rfl
  · rintro rfl
    exact one_mem _

/-- ★**`RamificationFiltration p` が空虚でないことの確認**(組み立てが本当に
`Interface` の構造体の 4 フィールドに嵌まることを型で検査するためのもの)。

★★**これは原典の `{Γ_K^v}` ではない**——`v > 0` で `Γ_K^v = ⊥` になる退化した族である
(`trivialStageFiltration_limit_eq_bot`)。`Skeleton/PGC/Section2.lean` の
`prop_2_1` / `prop_2_2` にこれを渡してはならないし、G2 の非空虚 witness にもならない。
本物を作るには冒頭「完全分岐の穴」の節が要求する段データが要る。 -/
noncomputable def trivialRamificationFiltration : RamificationFiltration p :=
  ramificationFiltrationOfStages (fun K => trivialStageFiltration K.absGal)

end ABC3.Found.PGC
