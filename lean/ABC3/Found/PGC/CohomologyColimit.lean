import ABC3.Found.PGC.Transgression
import Mathlib.Algebra.Colimit.Module
import Mathlib.Topology.Algebra.OpenSubgroup
import Mathlib.Topology.Algebra.ClopenNhdofOne

/-!
# [pGC] Proposition 1.1 —— 開正規部分群の塔についての colimit(離散側)

`Found/PGC/Transgression.lean` は 5 項完全列の右 2 箇所

  `H¹(S,A)^{G/S} --tg--> H²(G ⧸ S, Aˢ) --inf--> H²(G,A)`,  `range(tg) = ker(inf)`

を**仮定なしで**証明した。本ファイルはその **colimit 版**である。すなわち
`G` の正規部分群の塔 `S i ↓ 1` について

  `colim_i H¹(S i, A)`  と  `colim_i H²(G ⧸ S i, A^{S i}) → H²(G,A)`

を離散側で構成し、★**source の colimit が消えるなら比較射は極限で単射**を証明する。

## ★★★本ファイルの到達点(射程)

1. ★**有向系と colimit を構成した**:
   `h1Sys` / `H1Colim`(`H¹` の制限系)、`h2Sys` / `H2Colim`(`H²` の inflation 系)、
   `invH1Sys`(`H¹(S,A)^{G/S}` の系)、比較射 `H2ColimToH2`。
2. ★★**`colim_i H¹_sm(S i, A) = 0`(★仮定なし)** —— `smoothH1Colim_eq_zero`。
   ここで `H¹_sm` は「`1` の近傍で恒等的に消える 1-コサイクル」が定める部分加群であり、
   ★`↥A` が離散なら「`1` で連続な 1-コサイクル」と同値(`vanishesNearOne_of_continuousAt`)。
   仮定つきの全体版は `H1Colim_eq_zero`(仮定 `IsSmoothTowerH1`)。
3. ★★**比較射は極限で単射**(`H2ColimToH2_injective`)。
   ★colimit を経由しない使いやすい形は `exists_inflH2Step_eq_zero`。
4. **(c)-3(連続コホモロジーとの比較関手)には触れていない。** どこで止まったかは下の
   「★★(c)-3 との境界」を参照。

## ★★★★正直な注意 —— 離散側では `colim_i H¹(S i, A) = 0` は**偽**である

★★mathlib の `groupCohomology` は**離散**群コホモロジーであり、`H¹(S,A)` は
**連続とは限らない**コサイクルの類を含む。したがって

> `S` が `1` の近傍基底をなす開正規部分群の塔であっても
> `colim_i H¹(S i, A) = 0` は**一般には成り立たない**。

★反例(数学の水準。Lean では書いていない): `G = Ẑ`、`A = ℚ/ℤ`(自明作用)。
`H¹(S,A) = Hom(S, ℚ/ℤ)`(★抽象群の準同型)。`φ : Ẑ → ℚ/ℤ` の類が `nẐ` で消えるのは
`φ` が `Ẑ/nẐ = ℤ/n` を経由するときに限るが、連続でない `φ` はどの `nẐ` の上でも消えない。

★★**そこで本ファイルは 2 通りに分けた**。

* **仮定なしの版**: 「`1` の近傍で消える」1-コサイクルだけを見る
  (`smoothCocycles₁` / `smoothH1` / `smoothH1Colim_eq_zero`)。
  ★これは**定理であって仮定ではない**。★`smoothH1_eq_top_of_continuousAt` により、
  「どの 1-コサイクルも `1` で連続」なら `H¹_sm = H¹` に潰れる。
* **仮定つきの版**: `IsSmoothTowerH1`(どの段のどの 1-コサイクルも塔の先で消える)を
  仮定して `H1Colim_eq_zero` / `H2ColimToH2_injective` を得る。
  ★`IsSmoothTowerH1` は**空虚ではない**: `isSmoothTowerH1_of_continuousAt`(連続 + 離散)と
  `isSmoothTowerH1_of_exists_bot`(塔が `⊥` に達する場合)を置いてある。

★★**「通すための偽の仮説」は作っていない。** `IsSmoothTowerH1` は
連続コチェイン(= (c)-3)の設定でちょうど成り立つ性質を代数的に切り出したものである。

## ★★★測定 —— mathlib に何があり何が無いか(2026-09-07、コマンドと出力)

```
awk -F'\t' '$2 ~ /^Module\.DirectLimit/ {print $2"\t"$3}' .cache/mathlib-index.txt
```
→ **32 件**(`Module.DirectLimit` / `.of` / `.lift` / `.exists_of` / `.of_f` / `.lift_of` …)。
★**colimit の在庫は在る。** ★`Module.DirectLimit` は `[DecidableEq ι]` を要求する。
★★**`DirectedSystem` は `Module.DirectLimit` の定義にも `exists_of` にも要らない**
(要るのは `[Nonempty ι]` と `[IsDirectedOrder ι]` だけ)。
それでも本ファイルは有向系の公理を別に証明してある
(`resH1Step_self` / `resH1Step_trans` / `inflH2Step_self` / `inflH2Step_trans`)。

```
awk -F'\t' '$2 ~ /DirectedSystem|IsDirected|DirectLimit/ && $3 ~ /Colimit|Directed/ ...'
```
→ `DirectedSystem`(`Order/DirectedInverseSystem.lean:69`)、
`IsDirectedOrder`(`Order/Directed.lean:150`)、`DirectLimit`(同 `:113`)ほか。

```
grep -inE "OpenSubgroup|OpenNormalSubgroup" .cache/mathlib-index.txt | grep -iE "nhds|basis"
```
→ `ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one`
(`Topology/Algebra/ClopenNhdofOne.lean:48`)。★**副有限群で開正規部分群が `1` の近傍基底を
なすことは mathlib に在る**(束ねられていない形で使える)。これを
`isNhdsOneBasis_profiniteTower` に使った。

★★**無いもの**(直前の波の測定と同じ): `H^n(S,A)` への `G ⧸ S` の共役作用。
本ファイルも `Transgression.lean` の `IsGInvariantH1` / `invariantsH1` をそのまま使う。

★重複検査(#158): 本ファイルで新しく付ける 37 個の名前について
`grep -rl --include=*.lean` を `lean/ABC3` 全体に、`grep -c` を
`.cache/mathlib-index.txt` に走らせ、★**ABC3 側 0 件 / mathlib 側 0 件**を確認した。

## ★★設計 —— 抽象核と具体層

★★**§1 と §2 と §3 は「コホモロジー・分岐・付値・Galois の語彙が 1 語も出ない核」**である。

* **§1(colimit の核)** —— `Module.DirectLimit` と `Preorder` だけ。
  - `directLimit_eq_zero_of_eventually` —— 各元が先で消えるなら colimit は 0。
  - `lift_injective_of_eventually_zero` —— 単射性の素形。
  - ★★`lift_injective_of_ker_eq_range` —— **「source の colimit が消えるなら比較射は単射」**。
    ★これが本ファイルの抽象核で、`ker(g i) = range(t i)` と可換な四角形と
    source の消滅だけから結論する。★**コホモロジーを一切知らない。**
* **§2(離散性の核)** —— 位相群と `Zero` だけ。
  - ★`vanishesNearOne_of_continuousAt` —— **`DiscreteTopology` を使う唯一の補題**。
  - `exists_le_forall_eq_zero` —— 近傍基底 + 有向性で「塔の先で恒等的に 0」に変える。
* **§3(塔の核)** —— `[Group G] [AddCommGroup A]` と `ρ : G → A →+ A` だけ(#209)。
  - ★★`IsTgLift.of_le` —— **持ち上げは部分群を小さくしても持ち上げのまま**。
    ★これがそのまま「inflation と transgression が可換」(`inflH2Step_tgClass`)の中身。
    ★**証明は 1 行**(3 つのフィールドを `hle` で押し込むだけ)。

具体層(§4 以降)は核に代入するだけである。

## ★正規性・可判定性・選択公理をどこで使ったか(★直前の波の作法)

★**正規性を使う宣言**: `quotStep` / `invStep` / `inflH2Step` とその系(`G ⧸ S` を作るため)、
`invH1Step` / `inflH2Step_transgressionLin` / `h2Sys` / `H2Colim` / `H2ColimToH2` /
`invH1Sys` / `ker_inflH2Lin` / `exists_inflH2Step_eq_zero`。

★**正規性を使わない宣言**(★11 本): §1 の 4 本、§2 の 4 本、§3 の 3 本。
さらに `resH1Step` / `resH1Step_H1π` / `coe_mapCocycles₁_inclusion` / `resH1Step_self` /
`resH1Step_trans` / `resH1Step_H1π_eq_zero` / `h1Sys` / `H1Colim` / `IsSmoothTowerH1` /
`h1Sys_eventually_zero` / `H1Colim_eq_zero` / `subsingleton_H1Colim` /
`smoothCocycles₁` / `smoothH1` / `resH1Step_mem_smoothH1` / `smoothH1Sys` /
`SmoothH1Colim` / `smoothH1Colim_eq_zero` も**正規性を要しない**
(★`H¹` の制限は正規性と無関係だから)。

★**可判定性**: `Module.DirectLimit` が `[DecidableEq ι]` を要求するので、
colimit を作る宣言はすべて `[DecidableEq ι]` を**引数で受ける**。
★**`Classical.decEq` を勝手に差し込まない**(型が実例に依存してしまうため)。

★**選択公理**: `#print axioms` は次のとおり。
`IsTgLift.of_le` / `IsCocycleOn.of_le` / `IsInvarianceWitness.of_le` /
`exists_le_forall_eq_zero` は ★**`[propext, Quot.sound]` のみ**。
それ以外(colimit を触るもの・`Rep` を触るもの)は
`[propext, Classical.choice, Quot.sound]`。
★注意: 直前の波と同じく、`#print axioms` では「選択を使っていない」ことは確認できない
(`Module.DirectLimit` の構成そのものが `Classical.choice` を引く)。

## ★★(c)-3 との境界 —— どこで止まったか

★**触れていない。** 理由は次の 1 点である。

> mathlib には**連続コホモロジー**が次数 0 以外に無い(直前の波の測定)。
> したがって「`colim_i H^n(G ⧸ S i, A^{S i}) ≅ H^n_cont(G, A)`」を述べる**相手側の対象が無い**。

★本ファイルが (c)-3 に**寄せた**のは次の 2 点だけである。

* `VanishesNearOne` / `smoothCocycles₁` —— 連続コチェインの**次数 1 の影**を、
  位相を使わない代数的な述語として切り出した。
* `vanishesNearOne_of_continuousAt` —— その影が本当に「連続」から出ることを証明した。

★**残っているのは「連続コチェイン複体を次数 `n` で定義すること」**であり、
それは本持ち場の射程外である((c)-3)。

## ★退化の自己検査

* `resH1Step_self` / `resH1Step_trans` / `inflH2Step_self` / `inflH2Step_trans` ——
  有向系の公理(★`Module.DirectLimit` には要らないが、系が本当に関手的であることの確認)。
* `smoothCocycles₁_bot` —— `S = ⊥` ではすべての 1-コサイクルがなめらか。
  ★`Transgression.lean` の `invariantsH1_bot` と整合する。
* `H2ColimToH2_injective_of_exists_bot` —— 塔が `⊥` に達すれば
  比較射の単射性は**仮定なし**で出る。★`inflation₂_bijective_bot` の塔版。
* `smoothH1_eq_top_of_continuousAt` —— `H¹_sm` が `H¹` 全体になる場合。
  ★「なめらかな部分」が痩せすぎていないことの確認。
* `isNhdsOneBasis_profiniteTower` / `smoothH1Colim_profinite_eq_zero` ——
  ★**副有限群の開正規部分群の塔という本物の実例**で全部が動くことの確認。
* ★**`axiom` も `sorry` も置いていない。**

## 逸脱の記録

* **逸脱 1**: 原典 (pGC) は `Γ_K` の**連続**コホモロジーを使うが、本ファイルは
  mathlib の**離散**群コホモロジー `groupCohomology` で書いている。
  ★`Transgression.lean` / `InflationRestrictionH2.lean` / `GroupCohomologyFinite.lean` の
  逸脱 1 と同じ理由・同じ範囲。★上の「正直な注意」で、この逸脱が
  「`colim H¹ = 0`」の真偽を変えることを明示した。
* **逸脱 2**: 塔の添字を「開正規部分群の全体」に固定せず、
  **反単調な部分群の族 `S : ι → Subgroup G`** で受けた。
  ★原典は `Γ_K` の開部分群を走らせるが、本ファイルの主張はどれも
  「有向な塔」しか使わないためである。★副有限の実例は §13 に置いた。
* **逸脱 3**: 「`G ⧸ S`-不変部分」は `Transgression.lean` の `invariantsH1`
  (コチェインの述語)をそのまま使う。★mathlib に `H^n(S,A)` への共役作用が無いため。
-/

namespace ABC3.Found.PGC

open CategoryTheory groupCohomology

universe u v

/-- 原典の位置(★`Transgression.lean` の `transgression.src` と同じ項目)。

原文 (pGC p.3):
> Proposition 1.1: The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

★原典は 5 項完全列も transgression も colimit も名指ししていない(「well-known」で畳んでいる)。
本ファイルは mathlib 側の欠落(連続コホモロジーが次数 0 しか無い)を、
離散側の colimit で可能なところまで埋める位置づけである。 -/
def cohomologyColimit.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-! ## 1. ★★★抽象核 I —— 有向系の colimit(語彙ゼロ)

★`Module.DirectLimit` と `Preorder` しか出てこない。
★コホモロジー・分岐・付値・Galois の語彙は 1 語も無い。 -/

section AbstractCoreColimit

variable {R : Type*} [Ring R] {ι : Type*} [Preorder ι] [DecidableEq ι]
  {M : ι → Type*} [∀ i, AddCommGroup (M i)] [∀ i, Module R (M i)]

/-- 有向系のすべての元が「先で消える」なら colimit は 0。 -/
theorem directLimit_eq_zero_of_eventually (f : ∀ i j, i ≤ j → M i →ₗ[R] M j)
    [Nonempty ι] [IsDirectedOrder ι]
    (h : ∀ (i : ι) (x : M i), ∃ j, ∃ hij : i ≤ j, f i j hij x = 0)
    (z : Module.DirectLimit M f) : z = 0 := by
  obtain ⟨i, x, rfl⟩ := Module.DirectLimit.exists_of z
  obtain ⟨j, hij, hx⟩ := h i x
  rw [← Module.DirectLimit.of_f (f := f) (i := i) (j := j) (hij := hij) (x := x), hx, map_zero]

/-- `directLimit_eq_zero_of_eventually` の `Subsingleton` 版。 -/
theorem subsingleton_directLimit_of_eventually (f : ∀ i j, i ≤ j → M i →ₗ[R] M j)
    [Nonempty ι] [IsDirectedOrder ι]
    (h : ∀ (i : ι) (x : M i), ∃ j, ∃ hij : i ≤ j, f i j hij x = 0) :
    Subsingleton (Module.DirectLimit M f) :=
  ⟨fun z w => by
    rw [directLimit_eq_zero_of_eventually f h z, directLimit_eq_zero_of_eventually f h w]⟩

/-- 極限での単射性(素形) —— 各段で「`g i` の核に入る元は先で消える」なら
`lift g` は単射。 -/
theorem lift_injective_of_eventually_zero
    (f : ∀ i j, i ≤ j → M i →ₗ[R] M j) [Nonempty ι] [IsDirectedOrder ι]
    {P : Type*} [AddCommGroup P] [Module R P]
    (g : ∀ i, M i →ₗ[R] P) (Hg : ∀ i j hij x, g j (f i j hij x) = g i x)
    (h : ∀ (i : ι) (x : M i), g i x = 0 → ∃ j, ∃ hij : i ≤ j, f i j hij x = 0) :
    Function.Injective (Module.DirectLimit.lift R ι M f g Hg) := by
  rw [injective_iff_map_eq_zero]
  intro z hz
  obtain ⟨i, x, rfl⟩ := Module.DirectLimit.exists_of z
  rw [Module.DirectLimit.lift_of] at hz
  obtain ⟨j, hij, hx⟩ := h i x hz
  rw [← Module.DirectLimit.of_f (f := f) (i := i) (j := j) (hij := hij) (x := x), hx, map_zero]

/-- **★★★★本ファイルの抽象核** —— 「source の colimit が消えるなら比較射は極限で単射」。

各段で `ker(g i) = range(t i)` が成り立ち、`t` が系の射と可換で、
source の系 `N` のすべての元が先で消えるなら、`lift g` は単射である。

★★**コホモロジーを一切知らない**(`Module` と `Preorder` だけ)。 -/
theorem lift_injective_of_ker_eq_range
    {N : ι → Type*} [∀ i, AddCommGroup (N i)] [∀ i, Module R (N i)]
    (f : ∀ i j, i ≤ j → M i →ₗ[R] M j) (u : ∀ i j, i ≤ j → N i →ₗ[R] N j)
    (t : ∀ i, N i →ₗ[R] M i) [Nonempty ι] [IsDirectedOrder ι]
    {P : Type*} [AddCommGroup P] [Module R P]
    (g : ∀ i, M i →ₗ[R] P) (Hg : ∀ i j hij x, g j (f i j hij x) = g i x)
    (hsq : ∀ i j (hij : i ≤ j) (y : N i), f i j hij (t i y) = t j (u i j hij y))
    (hker : ∀ i, LinearMap.ker (g i) = LinearMap.range (t i))
    (hvan : ∀ (i : ι) (y : N i), ∃ j, ∃ hij : i ≤ j, u i j hij y = 0) :
    Function.Injective (Module.DirectLimit.lift R ι M f g Hg) := by
  refine lift_injective_of_eventually_zero f g Hg ?_
  intro i x hx
  have hmem : x ∈ LinearMap.range (t i) := by rw [← hker i]; exact hx
  obtain ⟨y, rfl⟩ := hmem
  obtain ⟨j, hij, hy⟩ := hvan i y
  exact ⟨j, hij, by rw [hsq i j hij y, hy, map_zero]⟩

end AbstractCoreColimit

/-! ## 2. ★★★抽象核 II —— 離散性(語彙ゼロ)

★位相群 `G` と `Zero` を持つ空間 `X` しか出てこない。
★**`DiscreteTopology` を使うのは `vanishesNearOne_of_continuousAt` の 1 本だけ**である。 -/

section AbstractCoreSmooth

open Topology

variable {G : Type*} [Group G] [TopologicalSpace G] {ι : Type*}

/-- 部分群の族 `S` が `1 ∈ G` の近傍基底をなす。 -/
def IsNhdsOneBasis (S : ι → Subgroup G) : Prop :=
  ∀ V ∈ 𝓝 (1 : G), ∃ i, ((S i : Subgroup G) : Set G) ⊆ V

/-- `T` 上の関数 `z` が `1` の近傍で恒等的に消える。 -/
def VanishesNearOne {X : Type*} [Zero X] (T : Subgroup G) (z : ↥T → X) : Prop :=
  ∃ V ∈ 𝓝 (1 : G), ∀ t : ↥T, (t : G) ∈ V → z t = 0

/-- ★★**`X` の離散性を使う唯一の場所**。`z` が `1` で連続で `z 1 = 0` なら、
`z` は `1` の近傍で恒等的に `0` である。 -/
theorem vanishesNearOne_of_continuousAt {X : Type*} [Zero X] [TopologicalSpace X]
    [DiscreteTopology X] (T : Subgroup G) (z : ↥T → X) (h1 : z 1 = 0)
    (hz : ContinuousAt z 1) : VanishesNearOne T z := by
  have hz0 : ({0} : Set X) ∈ 𝓝 (z (1 : ↥T)) := by
    rw [h1]; exact (isOpen_discrete ({0} : Set X)).mem_nhds rfl
  have hmem : z ⁻¹' ({0} : Set X) ∈ 𝓝 (1 : ↥T) := hz hz0
  rw [mem_nhds_subtype (T : Set G) (1 : ↥T) _] at hmem
  obtain ⟨V, hV, hsub⟩ := hmem
  refine ⟨V, by simpa using hV, fun t ht => ?_⟩
  exact hsub (by simpa using ht)

/-- 近傍基底と有向性から、`1` の近傍で消える関数は塔の**ある段の上で恒等的に** `0`。 -/
theorem exists_le_forall_eq_zero {X : Type*} [Zero X] [Preorder ι] [IsDirectedOrder ι]
    (S : ι → Subgroup G) (hS : Antitone S) (hb : IsNhdsOneBasis S) (i : ι)
    (z : ↥(S i) → X) (hz : VanishesNearOne (S i) z) :
    ∃ j, ∃ _ : i ≤ j, ∀ t : ↥(S i), (t : G) ∈ S j → z t = 0 := by
  obtain ⟨V, hV, hzV⟩ := hz
  obtain ⟨i', hi'⟩ := hb V hV
  obtain ⟨j, hij, hi'j⟩ := exists_ge_ge i i'
  exact ⟨j, hij, fun t ht => hzV t (hi' (hS hi'j ht))⟩

end AbstractCoreSmooth

/-! ## 3. ★★★抽象核 III —— 塔に沿った持ち上げ(語彙ゼロ)

★`[Group G] [AddCommGroup A]` と `ρ : G → A →+ A` だけ(#209)。
★★**`IsTgLift.of_le` が「inflation と transgression が可換」の中身のすべて**である。 -/

section AbstractCoreTower

variable {G : Type u} {A : Type v} [Group G] [AddCommGroup A]

/-- ★★**抽象核** —— transgression の持ち上げは部分群を小さくしても持ち上げのままである。

★`IsTgLift` の 3 つのフィールドはどれも「`s ∈ S`」の下でしか要求されないので、
`T ≤ S` に押し込むだけでよい。★**証明は 1 行**。 -/
theorem IsTgLift.of_le {ρ : G → A →+ A} {S T : Subgroup G} {f F : G → A}
    (hF : IsTgLift ρ S f F) (hle : T ≤ S) : IsTgLift ρ T f F :=
  ⟨hF.one, fun g s hs => hF.right g s (hle hs), fun g s hs => hF.conj g s (hle hs)⟩

/-- ★抽象核 —— `S` 上の 1-コサイクルは `T ≤ S` 上でも 1-コサイクル。 -/
theorem IsCocycleOn.of_le {ρ : G → A →+ A} {S T : Subgroup G} {f : G → A}
    (hf : IsCocycleOn ρ S f) (hle : T ≤ S) : IsCocycleOn ρ T f :=
  fun s hs t ht => hf s (hle hs) t (hle ht)

/-- ★抽象核 —— 不変性の証人も `T ≤ S` に落ちる。 -/
theorem IsInvarianceWitness.of_le {ρ : G → A →+ A} {S T : Subgroup G} {f a : G → A}
    (ha : IsInvarianceWitness ρ S f a) (hle : T ≤ S) : IsInvarianceWitness ρ T f a :=
  fun g s hs => ha g s (hle hs)

end AbstractCoreTower

/-! ## 4. 具体層 I —— `H¹` の制限(★正規性を要しない) -/

section ConcreteRestrictionH1

variable {k G : Type} [CommRing k] [Group G]

/-- **`H¹(S,A) → H¹(T,A)`**(`T ≤ S` への制限、`k`-線形写像)。

★塔 `S ↓ 1` に沿った `H¹` の有向系の 1 段。★`S` の正規性は要らない。 -/
noncomputable def resH1Step (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S) :
    (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k) →ₗ[k]
      (groupCohomology (Rep.res T.subtype A) 1 : ModuleCat k) :=
  ModuleCat.Hom.hom (groupCohomology.map (A := Rep.res S.subtype A) (B := Rep.res T.subtype A)
    (Subgroup.inclusion hle) (𝟙 _) 1)

/-- `resH1Step` のコサイクル水準での記述。 -/
theorem resH1Step_H1π (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S)
    (z : cocycles₁ (Rep.res S.subtype A)) :
    resH1Step A S T hle (ConcreteCategory.hom (H1π (Rep.res S.subtype A)) z)
      = ConcreteCategory.hom (H1π (Rep.res T.subtype A))
          (ConcreteCategory.hom (mapCocycles₁ (A := Rep.res S.subtype A)
            (B := Rep.res T.subtype A) (Subgroup.inclusion hle) (𝟙 _)) z) :=
  H1π_comp_map_apply (A := Rep.res S.subtype A) (B := Rep.res T.subtype A)
    (Subgroup.inclusion hle) (𝟙 _) z

/-- 制限されたコサイクルの値は元のコサイクルの値。 -/
theorem coe_mapCocycles₁_inclusion (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S)
    (z : cocycles₁ (Rep.res S.subtype A)) (t : ↥T) :
    ((ConcreteCategory.hom (mapCocycles₁ (A := Rep.res S.subtype A)
      (B := Rep.res T.subtype A) (Subgroup.inclusion hle) (𝟙 _)) z : ↥T → A)) t
      = (z : ↥S → A) (Subgroup.inclusion hle t) := by
  rw [coe_mapCocycles₁]
  rfl

/-- ★有向系の公理その 1 —— 恒等段は恒等射。 -/
theorem resH1Step_self (A : Rep k G) (S : Subgroup G)
    (x : (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k)) :
    resH1Step A S S le_rfl x = x := by
  induction x using groupCohomology.H1_induction_on with
  | h z => rw [resH1Step_H1π]; congr 1

/-- ★有向系の公理その 2 —— 段の合成則。 -/
theorem resH1Step_trans (A : Rep k G) (S T U : Subgroup G) (hTS : T ≤ S) (hUT : U ≤ T)
    (x : (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k)) :
    resH1Step A T U hUT (resH1Step A S T hTS x) = resH1Step A S U (hUT.trans hTS) x := by
  induction x using groupCohomology.H1_induction_on with
  | h z => rw [resH1Step_H1π, resH1Step_H1π, resH1Step_H1π]; congr 1

/-- ★★**`T` 上で恒等的に消えるコサイクルの類は `H¹(T,A)` で 0**。

★これが「`colim_i H¹(S i,A) = 0`」の唯一のコホモロジー的な中身である。 -/
theorem resH1Step_H1π_eq_zero (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S)
    (z : cocycles₁ (Rep.res S.subtype A))
    (hz : ∀ t : ↥S, (t : G) ∈ T → (z : ↥S → A) t = 0) :
    resH1Step A S T hle (ConcreteCategory.hom (H1π (Rep.res S.subtype A)) z) = 0 := by
  rw [resH1Step_H1π]
  have h0 : (ConcreteCategory.hom (mapCocycles₁ (A := Rep.res S.subtype A)
      (B := Rep.res T.subtype A) (Subgroup.inclusion hle) (𝟙 _)) z) = 0 := by
    ext t
    rw [coe_mapCocycles₁]
    show (Rep.Hom.hom (𝟙 (Rep.res (Subgroup.inclusion hle) (Rep.res S.subtype A))))
        ((z : ↥S → A) (Subgroup.inclusion hle t)) = (0 : ↥T → A) t
    rw [hz (Subgroup.inclusion hle t) t.2]
    simp
  rw [h0, map_zero]

end ConcreteRestrictionH1

/-! ## 5. ★★★`colim_i H¹(S i, A)`(仮定つきの版) -/

section H1Colimit

variable {k G : Type} [CommRing k] [Group G] {ι : Type} [Preorder ι]

/-- 塔 `S`(反単調な部分群の族)に沿った `H¹` の有向系。 -/
noncomputable def h1Sys (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) (i j : ι)
    (hij : i ≤ j) :
    (groupCohomology (Rep.res (S i).subtype A) 1 : ModuleCat k) →ₗ[k]
      (groupCohomology (Rep.res (S j).subtype A) 1 : ModuleCat k) :=
  resH1Step A (S i) (S j) (hS hij)

/-- ★★★**塔 `S` に沿った `colim_i H¹(S i, A)`**(離散側)。 -/
abbrev H1Colim (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) [DecidableEq ι] : Type :=
  Module.DirectLimit (fun i => (groupCohomology (Rep.res (S i).subtype A) 1 : ModuleCat k))
    (h1Sys A S hS)

/-- ★★塔 `S` に関する **1-コサイクルのなめらかさ** ——
どの段のどの 1-コサイクルも、塔のもっと先の段の上では恒等的に消える。

★★**離散側ではこれは一般には成り立たない**(ファイル冒頭の反例)。
成り立つ十分条件は `isSmoothTowerH1_of_continuousAt`(連続 + 離散)と
`isSmoothTowerH1_of_exists_bot`(塔が `⊥` に達する)である。 -/
def IsSmoothTowerH1 (A : Rep k G) (S : ι → Subgroup G) : Prop :=
  ∀ (i : ι) (z : cocycles₁ (Rep.res (S i).subtype A)),
    ∃ j, ∃ _ : i ≤ j, ∀ t : ↥(S i), (t : G) ∈ S j → (z : ↥(S i) → A) t = 0

/-- なめらかなら、`H¹` の有向系のすべての元は先で消える。 -/
theorem h1Sys_eventually_zero (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S)
    (hsm : IsSmoothTowerH1 A S) (i : ι)
    (x : (groupCohomology (Rep.res (S i).subtype A) 1 : ModuleCat k)) :
    ∃ j, ∃ hij : i ≤ j, h1Sys A S hS i j hij x = 0 := by
  induction x using groupCohomology.H1_induction_on with
  | h z =>
    obtain ⟨j, hij, hz⟩ := hsm i z
    exact ⟨j, hij, resH1Step_H1π_eq_zero A (S i) (S j) (hS hij) z hz⟩

/-- ★★★★**`colim_i H¹(S i, A) = 0`**(なめらかさを仮定)。
★transgression の source が極限で消える。 -/
theorem H1Colim_eq_zero (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) [DecidableEq ι]
    [Nonempty ι] [IsDirectedOrder ι] (hsm : IsSmoothTowerH1 A S) (z : H1Colim A S hS) : z = 0 :=
  directLimit_eq_zero_of_eventually (h1Sys A S hS) (h1Sys_eventually_zero A S hS hsm) z

/-- ★★★★`colim_i H¹(S i, A) = 0`(`Subsingleton` 版)。 -/
theorem subsingleton_H1Colim (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) [DecidableEq ι]
    [Nonempty ι] [IsDirectedOrder ι] (hsm : IsSmoothTowerH1 A S) :
    Subsingleton (H1Colim A S hS) :=
  subsingleton_directLimit_of_eventually (h1Sys A S hS) (h1Sys_eventually_zero A S hS hsm)

end H1Colimit

/-! ## 6. ★`IsSmoothTowerH1` が空虚でないこと -/

section H1SmoothNonvacuous

variable {k G : Type} [CommRing k] [Group G] {ι : Type} [Preorder ι]

/-- ★★★**`A` が離散であることを使う場所** —— どの段のどの 1-コサイクルも `1` で連続なら、
塔はなめらかである。

★`vanishesNearOne_of_continuousAt` が `DiscreteTopology ↥A` を使う唯一の補題。 -/
theorem isSmoothTowerH1_of_continuousAt [TopologicalSpace G] [IsDirectedOrder ι]
    (A : Rep k G) [TopologicalSpace ↥A] [DiscreteTopology ↥A] (S : ι → Subgroup G)
    (hS : Antitone S) (hb : IsNhdsOneBasis S)
    (hcont : ∀ (i : ι) (z : cocycles₁ (Rep.res (S i).subtype A)),
      ContinuousAt (z : ↥(S i) → ↥A) 1) :
    IsSmoothTowerH1 A S := fun i z =>
  exists_le_forall_eq_zero S hS hb i (z : ↥(S i) → ↥A)
    (vanishesNearOne_of_continuousAt (S i) (z : ↥(S i) → ↥A) (cocycles₁_map_one z) (hcont i z))

/-- ★退化の自己検査 —— 塔が `⊥` に達するなら、なめらかさは自動的に成り立つ。 -/
theorem isSmoothTowerH1_of_exists_bot (A : Rep k G) (S : ι → Subgroup G)
    (h : ∀ i, ∃ j, i ≤ j ∧ S j = ⊥) : IsSmoothTowerH1 A S := by
  intro i z
  obtain ⟨j, hij, hj⟩ := h i
  refine ⟨j, hij, fun t ht => ?_⟩
  rw [hj] at ht
  have h1 : t = 1 := Subtype.ext (by simpa using ht)
  rw [h1]
  exact cocycles₁_map_one z

end H1SmoothNonvacuous

/-! ## 7. 具体層 II —— `H²` の inflation(★正規性を使う) -/

section ConcreteInflationH2

variable {k G : Type} [CommRing k] [Group G]

/-- 塔の 1 段の商群の写像 `G ⧸ T →* G ⧸ S`(`T ≤ S`)。 -/
def quotStep (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S) : G ⧸ T →* G ⧸ S :=
  QuotientGroup.map T S (MonoidHom.id G) (by simpa using hle)

@[simp] theorem quotStep_mk (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S) (g : G) :
    quotStep S T hle (QuotientGroup.mk' T g) = QuotientGroup.mk' S g := rfl

/-- 部分群が小さくなれば不変部分は大きくなる。 -/
theorem invariants_mono (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S) :
    Representation.invariants (A.ρ.comp S.subtype) ≤
      Representation.invariants (A.ρ.comp T.subtype) := fun _ hx t => hx ⟨(t : G), hle t.2⟩

/-- 塔の 1 段の係数の写像 `Aˢ ↪ Aᵀ`(表現の射)。 -/
noncomputable def invStep (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S) :
    Rep.res (quotStep S T hle) (A.quotientToInvariants S) ⟶ A.quotientToInvariants T :=
  Rep.ofHom ⟨Submodule.inclusion (invariants_mono A S T hle), by
    intro q
    induction q using QuotientGroup.induction_on with | @H g => ext _; rfl⟩

/-- **`H²(G ⧸ S, Aˢ) → H²(G ⧸ T, Aᵀ)`**(inflation の 1 段、`k`-線形写像)。 -/
noncomputable def inflH2Step (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) :
    (groupCohomology (A.quotientToInvariants S) 2 : ModuleCat k) →ₗ[k]
      (groupCohomology (A.quotientToInvariants T) 2 : ModuleCat k) :=
  ModuleCat.Hom.hom (groupCohomology.map (A := A.quotientToInvariants S)
    (B := A.quotientToInvariants T) (quotStep S T hle) (invStep A S T hle) 2)

/-- `inflH2Step` のコサイクル水準での記述。 -/
theorem inflH2Step_H2π (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S)
    (x : cocycles₂ (A.quotientToInvariants S)) :
    inflH2Step A S T hle (ConcreteCategory.hom (H2π (A.quotientToInvariants S)) x)
      = ConcreteCategory.hom (H2π (A.quotientToInvariants T))
          (ConcreteCategory.hom (mapCocycles₂ (A := A.quotientToInvariants S)
            (B := A.quotientToInvariants T) (quotStep S T hle) (invStep A S T hle)) x) :=
  H2π_comp_map_apply (A := A.quotientToInvariants S) (B := A.quotientToInvariants T)
    (quotStep S T hle) (invStep A S T hle) x

/-- **`H²(G ⧸ S, Aˢ) → H²(G, A)`**(inflation、`k`-線形写像)。
★`Transgression.lean` の `range_transgressionLin` が使っているものと同じ射。 -/
noncomputable def inflH2Lin (A : Rep k G) (S : Subgroup G) [S.Normal] :
    (groupCohomology (A.quotientToInvariants S) 2 : ModuleCat k) →ₗ[k]
      (groupCohomology A 2 : ModuleCat k) :=
  ModuleCat.Hom.hom (groupCohomology.map (A := A.quotientToInvariants S) (B := A)
    (QuotientGroup.mk' S) (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)

/-- ★有向系の公理その 1 —— 恒等段は恒等射。 -/
theorem inflH2Step_self (A : Rep k G) (S : Subgroup G) [S.Normal]
    (x : (groupCohomology (A.quotientToInvariants S) 2 : ModuleCat k)) :
    inflH2Step A S S le_rfl x = x := by
  induction x using groupCohomology.H2_induction_on with
  | h y =>
    rw [inflH2Step_H2π]
    congr 1
    ext p q
    induction p using QuotientGroup.induction_on with | @H g =>
    induction q using QuotientGroup.induction_on with | @H h =>
    rfl

/-- ★有向系の公理その 2 —— 段の合成則。 -/
theorem inflH2Step_trans (A : Rep k G) (S T U : Subgroup G) [S.Normal] [T.Normal] [U.Normal]
    (hTS : T ≤ S) (hUT : U ≤ T)
    (x : (groupCohomology (A.quotientToInvariants S) 2 : ModuleCat k)) :
    inflH2Step A T U hUT (inflH2Step A S T hTS x) = inflH2Step A S U (hUT.trans hTS) x := by
  induction x using groupCohomology.H2_induction_on with
  | h y =>
    rw [inflH2Step_H2π, inflH2Step_H2π, inflH2Step_H2π]
    congr 1
    ext p q
    induction p using QuotientGroup.induction_on with | @H g =>
    induction q using QuotientGroup.induction_on with | @H h =>
    rfl

/-- `inflH2Lin` のコサイクル水準での記述。 -/
theorem inflH2Lin_H2π (A : Rep k G) (S : Subgroup G) [S.Normal]
    (y : cocycles₂ (A.quotientToInvariants S)) :
    inflH2Lin A S (ConcreteCategory.hom (H2π (A.quotientToInvariants S)) y)
      = ConcreteCategory.hom (H2π A) (ConcreteCategory.hom
          (mapCocycles₂ (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
            (Rep.ofHom (A.ρ.quotientToInvariants_lift S))) y) :=
  H2π_comp_map_apply (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
    (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) y

/-- ★★**inflation の推移性** —— 塔の段の inflation を合成すると全体の inflation になる。
★これが colimit から `H²(G,A)` への比較射が定義できる理由である。 -/
theorem inflH2Lin_inflH2Step (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S)
    (x : (groupCohomology (A.quotientToInvariants S) 2 : ModuleCat k)) :
    inflH2Lin A T (inflH2Step A S T hle x) = inflH2Lin A S x := by
  induction x using groupCohomology.H2_induction_on with
  | h y =>
    rw [inflH2Step_H2π, inflH2Lin_H2π, inflH2Lin_H2π]
    congr 1

end ConcreteInflationH2

/-! ## 8. ★★★★inflation と transgression の可換性 -/

section TransgressionTower

variable {k G : Type} [CommRing k] [Group G]

/-- ★★**inflation と transgression の可換性(類の水準)**。

★中身は抽象核 `IsTgLift.of_le` —— 同じ持ち上げ `F` が `T ≤ S` の持ち上げにもなるので、
`dF` を `G ⧸ T` に降ろしたものは `G ⧸ S` に降ろしたものの inflation にほかならない。 -/
theorem inflH2Step_tgClass (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S)
    {f F : G → A} (hF : IsTgLift (repAddHom A) S f F) :
    inflH2Step A S T hle (tgClass A S hF) = tgClass A T (hF.of_le hle) := by
  show inflH2Step A S T hle (ConcreteCategory.hom (H2π (A.quotientToInvariants S))
      (tgCocycle A S hF))
    = ConcreteCategory.hom (H2π (A.quotientToInvariants T)) (tgCocycle A T (hF.of_le hle))
  rw [inflH2Step_H2π]
  congr 1
  ext p q
  induction p using QuotientGroup.induction_on with | @H g =>
  induction q using QuotientGroup.induction_on with | @H h =>
  rfl

/-- ★制限は `G ⧸ S`-不変性を保つ。 -/
theorem isGInvariantH1_res (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S)
    {c : (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k)} (hc : IsGInvariantH1 A S c) :
    IsGInvariantH1 A T (resH1Step A S T hle c) := by
  obtain ⟨f, F, hF, e⟩ := hc
  refine ⟨f, F, hF.of_le hle, ?_⟩
  rw [← e, resH1Step_H1π]
  congr 1

/-- **`H¹(S,A)^{G/S} → H¹(T,A)^{G/T}`**(transgression の source の有向系の 1 段)。 -/
noncomputable def invH1Step (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal] (hle : T ≤ S) :
    invariantsH1 A S →ₗ[k] invariantsH1 A T where
  toFun c := ⟨resH1Step A S T hle (c : _), isGInvariantH1_res A S T hle c.2⟩
  map_add' _ _ := Subtype.ext (map_add _ _ _)
  map_smul' _ _ := Subtype.ext (map_smul _ _ _)

@[simp] theorem invH1Step_coe (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (c : invariantsH1 A S) :
    ((invH1Step A S T hle c : invariantsH1 A T) : _) = resH1Step A S T hle (c : _) := rfl

/-- ★★★★**可換な四角形** —— `inf ∘ tg = tg ∘ res`。

★これが「比較射が極限で単射になる」ことの配管である。 -/
theorem inflH2Step_transgressionLin (A : Rep k G) (S T : Subgroup G) [S.Normal] [T.Normal]
    (hle : T ≤ S) (c : invariantsH1 A S) :
    inflH2Step A S T hle (transgressionLin A S c)
      = transgressionLin A T (invH1Step A S T hle c) := by
  obtain ⟨f, F, hF, e⟩ := id c.2
  have hfc : (ConcreteCategory.hom (H1π (Rep.res T.subtype A)))
      ⟨fun t : T => f (t : G),
        mem_cocycles₁_restrict A T ((hF.of_le hle).isCocycleOn (repAddHom_one A))⟩
      = ((invH1Step A S T hle c : invariantsH1 A T) : _) := by
    show _ = resH1Step A S T hle (c : _)
    rw [← e, resH1Step_H1π]
    congr 1
  rw [transgressionLin_apply, transgression_eq A S c.2 hF e, inflH2Step_tgClass,
    transgressionLin_apply, transgression_eq A T (invH1Step A S T hle c).2 (hF.of_le hle) hfc]

end TransgressionTower

/-! ## 9. ★★★★★`colim_i H²(G ⧸ S i, A^{S i})` と比較射 -/

section H2Colimit

variable {k G : Type} [CommRing k] [Group G] {ι : Type} [Preorder ι]

/-- 塔 `S` に沿った `H²(G ⧸ S, Aˢ)` の有向系(inflation)。 -/
noncomputable def h2Sys (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal] (hS : Antitone S)
    (i j : ι) (hij : i ≤ j) :
    (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k) →ₗ[k]
      (groupCohomology (A.quotientToInvariants (S j)) 2 : ModuleCat k) :=
  inflH2Step A (S i) (S j) (hS hij)

/-- ★★★★**塔 `S` に沿った `colim_i H²(G ⧸ S i, A^{S i})`**(離散側)。 -/
abbrev H2Colim (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal] (hS : Antitone S)
    [DecidableEq ι] : Type :=
  Module.DirectLimit
    (fun i => (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k)) (h2Sys A S hS)

/-- ★★★★**比較射** `colim_i H²(G ⧸ S i, A^{S i}) → H²(G, A)`。 -/
noncomputable def H2ColimToH2 (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal]
    (hS : Antitone S) [DecidableEq ι] :
    H2Colim A S hS →ₗ[k] (groupCohomology A 2 : ModuleCat k) :=
  Module.DirectLimit.lift k ι _ (h2Sys A S hS) (fun i => inflH2Lin A (S i))
    (fun i j hij x => inflH2Lin_inflH2Step A (S i) (S j) (hS hij) x)

/-- 比較射の各段での値は inflation そのもの。 -/
theorem H2ColimToH2_of (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal] (hS : Antitone S)
    [DecidableEq ι] (i : ι)
    (x : (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k)) :
    H2ColimToH2 A S hS
        (Module.DirectLimit.of k ι
          (fun i => (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k))
          (h2Sys A S hS) i x)
      = inflH2Lin A (S i) x :=
  Module.DirectLimit.lift_of (fun i => inflH2Lin A (S i))
    (fun i j hij x => inflH2Lin_inflH2Step A (S i) (S j) (hS hij) x) x

/-- 塔 `S` に沿った `H¹(S i, A)^{G/S i}` の有向系。 -/
noncomputable def invH1Sys (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal]
    (hS : Antitone S) (i j : ι) (hij : i ≤ j) :
    invariantsH1 A (S i) →ₗ[k] invariantsH1 A (S j) :=
  invH1Step A (S i) (S j) (hS hij)

/-- ★`ker(inf) = range(tg)`(`Transgression.lean` の `range_transgressionLin` の言い換え)。 -/
theorem ker_inflH2Lin (A : Rep k G) (S : Subgroup G) [S.Normal] :
    LinearMap.ker (inflH2Lin A S) = LinearMap.range (transgressionLin A S) :=
  (range_transgressionLin A S).symm

/-- なめらかなら、transgression の source の有向系も先で消える。 -/
theorem invH1Sys_eventually_zero (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal]
    (hS : Antitone S) (hsm : IsSmoothTowerH1 A S) (i : ι) (y : invariantsH1 A (S i)) :
    ∃ j, ∃ hij : i ≤ j, invH1Sys A S hS i j hij y = 0 := by
  obtain ⟨j, hij, h⟩ := h1Sys_eventually_zero A S hS hsm i (y : _)
  exact ⟨j, hij, Subtype.ext h⟩

/-- ★★★★★**比較射は極限で単射である**。

★`ker(inf) = range(tg)`(★仮定なし、`Transgression.lean` の `range_transgressionLin`)と
`colim_i H¹(S i,A) = 0`(なめらかさ)を合わせたもの。
★★中身はすべて抽象核 `lift_injective_of_ker_eq_range` に入っている。 -/
theorem H2ColimToH2_injective (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal]
    (hS : Antitone S) [DecidableEq ι] [Nonempty ι] [IsDirectedOrder ι]
    (hsm : IsSmoothTowerH1 A S) : Function.Injective (H2ColimToH2 A S hS) :=
  lift_injective_of_ker_eq_range (h2Sys A S hS) (invH1Sys A S hS)
    (fun i => transgressionLin A (S i)) (fun i => inflH2Lin A (S i)) _
    (fun i j hij y => inflH2Step_transgressionLin A (S i) (S j) (hS hij) y)
    (fun i => ker_inflH2Lin A (S i)) (invH1Sys_eventually_zero A S hS hsm)

end H2Colimit

/-! ## 10. ★★★★★仮定なしの `colim = 0` —— なめらかな部分 -/

section SmoothH1Colimit

open Topology

variable {k G : Type} [CommRing k] [Group G] [TopologicalSpace G] {ι : Type} [Preorder ι]

/-- **`1` の近傍で恒等的に消える 1-コサイクル**のなす部分加群。

★`↥A` が離散なら「`1` で連続な 1-コサイクル」と同値(`vanishesNearOne_of_continuousAt`)。 -/
def smoothCocycles₁ (A : Rep k G) (S : Subgroup G) :
    Submodule k (cocycles₁ (Rep.res S.subtype A)) where
  carrier := {z | VanishesNearOne S ((z : ↥S → ↥A))}
  add_mem' := by
    rintro z w ⟨V, hV, hz⟩ ⟨W, hW, hw⟩
    refine ⟨V ∩ W, Filter.inter_mem hV hW, fun t ht => ?_⟩
    show (z : ↥S → ↥A) t + (w : ↥S → ↥A) t = 0
    rw [hz t ht.1, hw t ht.2, add_zero]
  zero_mem' := ⟨Set.univ, Filter.univ_mem, fun _ _ => rfl⟩
  smul_mem' a z := by
    rintro ⟨V, hV, hz⟩
    refine ⟨V, hV, fun t ht => ?_⟩
    show a • (z : ↥S → ↥A) t = 0
    rw [hz t ht, smul_zero]

/-- **なめらかな類のなす `H¹(S,A)` の部分加群** `H¹_sm(S,A)`。 -/
noncomputable def smoothH1 (A : Rep k G) (S : Subgroup G) :
    Submodule k (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k) :=
  Submodule.map (ModuleCat.Hom.hom (H1π (Rep.res S.subtype A))) (smoothCocycles₁ A S)

theorem mem_smoothH1 (A : Rep k G) (S : Subgroup G)
    {c : (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k)} :
    c ∈ smoothH1 A S ↔ ∃ z : cocycles₁ (Rep.res S.subtype A),
      VanishesNearOne S ((z : ↥S → ↥A)) ∧
        ConcreteCategory.hom (H1π (Rep.res S.subtype A)) z = c :=
  Submodule.mem_map

/-- なめらかさは制限で保たれる。 -/
theorem resH1Step_mem_smoothH1 (A : Rep k G) (S T : Subgroup G) (hle : T ≤ S)
    {c : (groupCohomology (Rep.res S.subtype A) 1 : ModuleCat k)} (hc : c ∈ smoothH1 A S) :
    resH1Step A S T hle c ∈ smoothH1 A T := by
  obtain ⟨z, ⟨V, hV, hz⟩, rfl⟩ := (mem_smoothH1 A S).1 hc
  refine (mem_smoothH1 A T).2 ⟨ConcreteCategory.hom (mapCocycles₁ (A := Rep.res S.subtype A)
    (B := Rep.res T.subtype A) (Subgroup.inclusion hle) (𝟙 _)) z, ⟨V, hV, fun t ht => ?_⟩, ?_⟩
  · rw [coe_mapCocycles₁_inclusion]
    exact hz (Subgroup.inclusion hle t) ht
  · exact (resH1Step_H1π A S T hle z).symm

/-- 塔 `S` に沿った `H¹_sm` の有向系。 -/
noncomputable def smoothH1Sys (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) (i j : ι)
    (hij : i ≤ j) : smoothH1 A (S i) →ₗ[k] smoothH1 A (S j) :=
  LinearMap.restrict (resH1Step A (S i) (S j) (hS hij))
    (fun _ hx => resH1Step_mem_smoothH1 A (S i) (S j) (hS hij) hx)

/-- ★★★★**なめらかな部分の colimit** `colim_i H¹_sm(S i, A)`。 -/
abbrev SmoothH1Colim (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S) [DecidableEq ι] : Type :=
  Module.DirectLimit (fun i => (smoothH1 A (S i) : Type)) (smoothH1Sys A S hS)

/-- ★★★★★**仮定なしの `colim_i H¹_sm(S i, A) = 0`**。

★`S` が `1` の近傍基底で有向であれば、なめらかな類は**必ず**塔の先で消える。
★★`↥A` の離散性はここには現れない —— 離散性は
「連続コサイクル ⟹ なめらか」(`vanishesNearOne_of_continuousAt`)の側で使われる。 -/
theorem smoothH1Colim_eq_zero (A : Rep k G) (S : ι → Subgroup G) (hS : Antitone S)
    (hb : IsNhdsOneBasis S) [DecidableEq ι] [Nonempty ι] [IsDirectedOrder ι]
    (z : SmoothH1Colim A S hS) : z = 0 := by
  refine directLimit_eq_zero_of_eventually (smoothH1Sys A S hS) (fun i y => ?_) z
  obtain ⟨w, hw, hwy⟩ := (mem_smoothH1 A (S i)).1 y.2
  obtain ⟨j, hij, hz⟩ := exists_le_forall_eq_zero S hS hb i ((w : ↥(S i) → ↥A)) hw
  refine ⟨j, hij, Subtype.ext ?_⟩
  show resH1Step A (S i) (S j) (hS hij) (y : _) = 0
  rw [← hwy]
  exact resH1Step_H1π_eq_zero A (S i) (S j) (hS hij) w hz

end SmoothH1Colimit

/-! ## 11. ★退化の自己検査と、colimit を経由しない形 -/

section Degeneracy

open Topology

variable {k G : Type} [CommRing k] [Group G] {ι : Type} [Preorder ι]

/-- `H1π` は全射(`H1_induction_on` の言い換え)。 -/
theorem H1π_surjective (B : Rep k G) :
    Function.Surjective (ConcreteCategory.hom (H1π B)) := by
  intro c
  induction c using groupCohomology.H1_induction_on with | h z => exact ⟨z, rfl⟩

/-- ★★**colimit を経由しない使いやすい形** —— `inf` の核は塔の先で消える。 -/
theorem exists_inflH2Step_eq_zero (A : Rep k G) (S : ι → Subgroup G) [∀ i, (S i).Normal]
    (hS : Antitone S) (hsm : IsSmoothTowerH1 A S) (i : ι)
    (x : (groupCohomology (A.quotientToInvariants (S i)) 2 : ModuleCat k))
    (h : inflH2Lin A (S i) x = 0) :
    ∃ j, ∃ hij : i ≤ j, inflH2Step A (S i) (S j) (hS hij) x = 0 := by
  have hmem : x ∈ LinearMap.range (transgressionLin A (S i)) := by
    rw [← ker_inflH2Lin]; exact h
  obtain ⟨y, rfl⟩ := hmem
  obtain ⟨j, hij, hy⟩ := invH1Sys_eventually_zero A S hS hsm i y
  refine ⟨j, hij, ?_⟩
  rw [inflH2Step_transgressionLin A (S i) (S j) (hS hij) y]
  show transgressionLin A (S j) (invH1Sys A S hS i j hij y) = 0
  rw [hy, map_zero]

/-- ★退化の自己検査 —— 塔が `⊥` に達すれば、比較射の単射性は**仮定なし**で出る。
★`InflationRestrictionH2.lean` の `inflation₂_bijective_bot` の塔版。 -/
theorem H2ColimToH2_injective_of_exists_bot (A : Rep k G) (S : ι → Subgroup G)
    [∀ i, (S i).Normal] (hS : Antitone S) [DecidableEq ι] [Nonempty ι] [IsDirectedOrder ι]
    (hbot : ∀ i, ∃ j, i ≤ j ∧ S j = ⊥) : Function.Injective (H2ColimToH2 A S hS) :=
  H2ColimToH2_injective A S hS (isSmoothTowerH1_of_exists_bot A S hbot)

/-- ★退化の自己検査 —— `S = ⊥` ではすべての 1-コサイクルがなめらか。
★`Transgression.lean` の `invariantsH1_bot` と整合する。 -/
theorem smoothCocycles₁_bot [TopologicalSpace G] (A : Rep k G) :
    smoothCocycles₁ A (⊥ : Subgroup G) = ⊤ := by
  refine eq_top_iff.2 (fun z _ => ⟨Set.univ, Filter.univ_mem, fun t _ => ?_⟩)
  have h1 : t = 1 := Subtype.ext (Subgroup.mem_bot.1 t.2)
  rw [h1]
  exact cocycles₁_map_one z

/-- ★★**`H¹_sm = H¹` になる場合** —— `↥A` が離散でどの 1-コサイクルも `1` で連続なら、
なめらかな部分は `H¹(S,A)` 全体である。
★このとき `smoothH1Colim_eq_zero` は「`colim_i H¹(S i,A) = 0`」そのものになる。 -/
theorem smoothH1_eq_top_of_continuousAt [TopologicalSpace G] (A : Rep k G) [TopologicalSpace ↥A]
    [DiscreteTopology ↥A] (S : Subgroup G)
    (hcont : ∀ z : cocycles₁ (Rep.res S.subtype A), ContinuousAt ((z : ↥S → ↥A)) 1) :
    smoothH1 A S = ⊤ := by
  refine eq_top_iff.2 (fun c _ => ?_)
  obtain ⟨z, rfl⟩ := H1π_surjective (Rep.res S.subtype A) c
  exact (mem_smoothH1 A S).2 ⟨z, vanishesNearOne_of_continuousAt S ((z : ↥S → ↥A))
    (cocycles₁_map_one z) (hcont z), rfl⟩

end Degeneracy

/-! ## 12. ★★★★副有限群の実例 —— 開正規部分群の塔

★★**ここが「空虚でない」ことの最終確認である。** 副有限群 `G` の開正規部分群の全体
(包含の逆順)は有向で `1` の近傍基底をなし、上の機械がそのまま動く。 -/

section ProfiniteTower

open Topology

variable (G : Type) [Group G] [TopologicalSpace G]

/-- ★副有限群の**開正規部分群の全体**(包含の逆順)を添字とする塔。 -/
def profiniteTower : (OpenNormalSubgroup G)ᵒᵈ → Subgroup G :=
  fun H => ((OrderDual.ofDual H : OpenNormalSubgroup G) : Subgroup G)

theorem nonempty_openNormalSubgroup : Nonempty (OpenNormalSubgroup G) :=
  ⟨⟨(⊤ : OpenSubgroup G), by
    show ((⊤ : OpenSubgroup G) : Subgroup G).Normal
    exact Subgroup.normal_top⟩⟩

theorem antitone_profiniteTower : Antitone (profiniteTower G) := fun _ _ h => h

instance instNormalProfiniteTower (H : (OpenNormalSubgroup G)ᵒᵈ) : (profiniteTower G H).Normal :=
  (OrderDual.ofDual H : OpenNormalSubgroup G).isNormal'

variable [IsTopologicalGroup G] [CompactSpace G] [TotallyDisconnectedSpace G]

/-- ★副有限群では開正規部分群の全体が `1` の近傍基底をなす
(★mathlib の `ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one`)。 -/
theorem isNhdsOneBasis_profiniteTower : IsNhdsOneBasis (profiniteTower G) := by
  intro V hV
  obtain ⟨U, hUV, hU, hU1⟩ := mem_nhds_iff.1 hV
  obtain ⟨H, hH⟩ := ProfiniteGrp.exist_openNormalSubgroup_sub_open_nhds_of_one hU hU1
  exact ⟨OrderDual.toDual H, hH.trans hUV⟩

variable {k : Type} [CommRing k] {G}

/-- ★★★★**副有限群の開正規部分群の塔について、なめらかな部分の colimit は 0**
(★仮定なし)。 -/
theorem smoothH1Colim_profinite_eq_zero (A : Rep k G)
    [DecidableEq ((OpenNormalSubgroup G)ᵒᵈ)]
    (z : SmoothH1Colim A (profiniteTower G) (antitone_profiniteTower G)) : z = 0 := by
  haveI : Nonempty ((OpenNormalSubgroup G)ᵒᵈ) := nonempty_openNormalSubgroup G
  exact smoothH1Colim_eq_zero A (profiniteTower G) (antitone_profiniteTower G)
    (isNhdsOneBasis_profiniteTower G) z

end ProfiniteTower

end ABC3.Found.PGC
