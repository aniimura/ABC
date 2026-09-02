/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm21Chain
import ABC3.Found.GenEll.Thm21Extract
import ABC3.Found.GenEll.Thm21DegRatio
import ABC3.Found.GenEll.NonArchBound
import ABC3.Meta.Claim

/-!
# 第 1093 ブロック —— **`Theorem 2.1` の逆向き `(ii) ⟹ (i)`**（`Skeleton`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★★★★★★★★★これは何か

`Skeleton/GenEll/Section2.lean` の `theorem_2_1` が取っているのは **`(i) ⟹ (ii)`
だけ**である（同ファイルが明記している）。原文の実質は逆向き `(ii) ⟹ (i)`——
p.12–13 の 2 ページ——であり、そこが未了である。

★本ファイルは**その 2 ページを節点に割る**。原文の証明は 2 段に分かれる。

### 段 A（`D = ∅` への帰着、p.12 前半）

`Y → X` を取り `d′ ≔ d·deg(Y/X)` と置くと、`D = ∅` の場合に帰着する。
☆道具は `Prop 1.7, (i)(ii)` と `Prop 1.4, (i)(ii)(iii)`、
そして「`U_X(ℚ̄)^{≤d}` の点が `U_Y(ℚ̄)^{≤d′}` に持ち上がる」ことである。

### 段 B（`D = ∅` の場合、p.12 後半–p.13）

背理法。不等式が偽なら反例の点列が取れ、コンパクト性で `Ξ` と各 `v` の `Ξ_v` を取る。
★そこに **[NCBelyi] `Theorem 2.5`**（noncritical Belyi map `φ : X → ℙ¹`）を当てて
`K_V` を作り、`E ≔ φ⁻¹(C)_red` と置いて高さを計算すると矛盾が出る。

## ★★★★★★★★既にあるもの（`Found/`、`sorry` 無し）

| 部品 | 場所 |
|---|---|
| `BDge` の推移・制限 | `Thm21Chain.lean`（`bdge_trans`・`bdge_restrict`） |
| 段 A の数値の鎖 | `Thm21Chain.lean`（`thm_2_1_stepA`） |
| **段 B の数値の鎖**（`h1`–`h6` を受ける） | `Thm21Chain.lean`（`thm_2_1_stepB`） |
| 反例列・部分列・同時収束 | `Thm21Extract.lean` |
| Riemann–Hurwitz の帳簿・`deg` の比 | `Thm21DegRatio.lean` |
| `[NCBelyi]` の `Lemma 2.1`–`2.4`・`Theorem 2.5`（`ℙ¹` 版） | `Found/NCBelyi/`（21 ファイル、`sorry` 0） |

★★**数値の鎖は全部ある。** 残っているのは**幾何の側の入力を作ること**である。

## ★★★★★★★★★★進捗枠 **6 / 8**（2026-09-01、第 1096 時点）

| # | 節点 | 内容 | 状態 | 重み |
|---|---|---|---|---|
| 1 | `exists_belyi_noncritical_general` | `[NCBelyi] Thm 2.5` の**一般曲線版**（現在は `ℙ¹` 版のみ） | ❌ | 12 |
| 2 | `exists_compactlyBounded_of_nbhd` | compact domain の**有限合併は compact domain** | ✅ | 4 |
| 3 | `belyi_image_subset_KV` | 収束 ⟹ 有限個を除いて `K_V` に入る | ✅ | 2 |
| 4 | `pullback_omega_eq_omega_add_E` | `φ*ω_ℙ(C) ≅ ω_X(E)` から `h2` | ✅ | 8 |
| 5 | `htE_ge_ratio_htOmega` | `Prop 1.4` から `h6` | ✅ | 3 |
| 6 | `logDiff_pullback_le` | `Prop 1.7, (i)` を `e = 1` で当てて `h4` | ✅ | 4 |
| 7 | `exists_lift_of_degLe` | 段 A の点の持ち上げ（`d′ = d·deg(Y/X)`） | ✅ | 3 |
| 8 | `theorem_2_1_converse` | 組み立て（★これが実質） | ❌ | 8 |

★★**正直な読み方**——閉じた 5 つのうち **4・5・6 はグルー**であり、
実質（幾何の入力 `hpull`・`hprop` を作ること）は `.needs` の
`implicitStep` に重みつきで残してある。
☆節点 3・7 は短いが実質を持つ（収束からの近傍、次数の帳簿）。

★**残重み 12 + 4 + 8 = 24 / 44**。最大は節点 1。

☆`h1`（`ht_ω ≈ ht_{ω(E)} − ht_E`、`Prop 1.4`）と `h5`（`log-cond_E ≳ ht_E`、`Prop 1.6`）は
在庫（`Prop14.lean`・`Prop16Proper.lean`、ともに条なし `.src`）から出る見込みなので節点にしない。
★★**節点が増えたら枠の分母を上げて記録する**（この表を更新する）。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GenEll

/-! ## ★★★★節点 1 —— `[NCBelyi] Theorem 2.5` の一般曲線版 -/

/-- ★★★★★★★★**節点 1**——`[NCBelyi] Theorem 2.5` の一般曲線版。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`Found/NCBelyi/Thm25P1.lean` は `ℙ¹` 上の有限集合 `S : Finset ℂ` に対する版である。
★原文の `Theorem 2.5` は**一般の滑らかな固有連結曲線 `X`** に対するもので、
そこへの帰着に Riemann–Roch（`X` から `ℙ¹` への写像の構成）が要る。

★★受ける形は「各 `v` の有限集合 `Ξ_v` の点を `U_ℙ` に送る `φ` が取れる」である。

★★★★★★★★**欠陥（2026-09-02、第 1389）——この主張は偽である**

`UP` を**任意の集合**にしているので、`UP = ∅` かつ `Ξ_v ≠ ∅` で偽である。
☆さらに `φ` は**単なる函数**であって射ではないので、
`[NCBelyi] Theorem 2.5` の内容（`ℙ¹ ∖ {0,1,∞}` の上で不分岐な射）が写っていない。

★★★**したがって `sorry` を埋めようとしても永久に埋まらない**。
★先に**忠実な主張へ書き直す**こと。
☆機械検証: `Check/GenEll/BelyiGeneralVacuous.lean` の
`exists_belyi_noncritical_general_false`。 -/
theorem exists_belyi_noncritical_general
    {Curve : Type} {Pt : Curve → Type} {P1 : Type} (X : Curve)
    (V : Type) (Xi : V → Set (Pt X)) (hfin : ∀ v, (Xi v).Finite)
    (UP : Set P1) :
    ∃ φ : Pt X → P1, (∀ v, ∀ x ∈ Xi v, φ x ∈ UP) := by
  sorry

def exists_belyi_noncritical_general.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 1——[NCBelyi] Thm 2.5 の一般曲線版)",
    sectionId := "genell-thm-2-1" }

def exists_belyi_noncritical_general.needs : List ProofObligation :=
  [ .otherPaper "[NCBelyi]" "Theorem 2.5(一般曲線版。ℙ¹ 版は Found/NCBelyi/Thm25P1.lean に済)" 5,
    .implicitStep
      ("★`X = ℙ¹` への帰着に Riemann–Roch が要る" ++
       "（`Check/NCBelyi/Thm25AxiomGap.lean` が測定済み）。") 12,
    .implicitStep
      ("★★★★**2026-09-02（第 1389）——現在の主張は偽である**。" ++
       "`UP` が任意の集合なので `UP = ∅` かつ `Ξ_v ≠ ∅` で反例になる" ++
       "（`Check/GenEll/BelyiGeneralVacuous.lean` で機械検証）。" ++
       "☆`φ` も単なる函数で射ではない。" ++
       "★先に忠実な主張へ書き直すこと——" ++
       "`Found/NCBelyi/` の `ℙ¹` 版（21 ファイル、`sorry` 0）が語彙を持っている。") 12 ]

/-! ## ★★★★節点 2-3 —— `K_V` の構成 -/

/-- ★★★★★★**節点 2**——小近傍の像の Galois 共役の合併は compactly bounded。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆原文の `K_V` の作り方そのものである。 -/
theorem exists_compactlyBounded_of_nbhd
    {P1 : Type} [TopologicalSpace P1] [T2Space P1]
    {V : Type} [Finite V] (nbhd : V → Set P1)
    (hcd : ∀ v, IsCompactDomain (nbhd v)) (UP : Set P1)
    (hsub : ∀ v, nbhd v ⊆ UP) :
    ∃ KV : Set P1, IsCompactDomain KV ∧ KV ⊆ UP ∧ ∀ v, nbhd v ⊆ KV := by
  refine ⟨⋃ v, nbhd v, ⟨?_, ?_⟩, Set.iUnion_subset hsub, fun v => Set.subset_iUnion _ v⟩
  · exact isCompact_iUnion (fun v => (hcd v).isCompact)
  · have hclosed : IsClosed (⋃ v, nbhd v) :=
      (isCompact_iUnion (fun v => (hcd v).isCompact)).isClosed
    refine Set.Subset.antisymm ?_ ?_
    · refine Set.iUnion_subset (fun v => ?_)
      rw [(hcd v).eq_closure_interior]
      exact closure_mono (interior_mono (Set.subset_iUnion _ v))
    · exact hclosed.closure_subset_iff.2 interior_subset

def exists_compactlyBounded_of_nbhd.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 2——K_V を小近傍の像の共役の合併として作る)",
    sectionId := "genell-thm-2-1" }

def exists_compactlyBounded_of_nbhd.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Example 1.3, (ii)(compactly bounded subset の support)" 5,
    .implicitStep "☆`Found/GenEll/NonArchBound.lean` に `Example 1.3` の条なし実装がある。" 4 ]

/-- ★★★★★★**節点 3**——`φ(Ξ) ⊆ K_V`。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`Ξ` の共役組が `Ξ_v` に収束するので、有限個を除いて小近傍に入る。 -/
theorem belyi_image_subset_KV
    {B : Type} [TopologicalSpace B] {A : Type} (φ : A → B) (ξ : ℕ → A) (b : B)
    (hconv : Filter.Tendsto (fun n => φ (ξ n)) Filter.atTop (nhds b))
    (KV : Set B) (hKV : KV ∈ nhds b) :
    ∀ᶠ n in Filter.atTop, φ (ξ n) ∈ KV :=
  hconv.eventually_mem hKV

def belyi_image_subset_KV.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 3——φ(Ξ) ⊆ K_V)",
    sectionId := "genell-thm-2-1" }

def belyi_image_subset_KV.needs : List ProofObligation :=
  [ .implicitStep "☆収束から「有限個を除いて近傍に入る」を出す段。" 2 ]

/-! ## ★★★★節点 4-6 —— 段 B の入力 `h2`・`h6`・`h4` -/

/-- ★★★★★★★★**節点 4**——`E ≔ φ⁻¹(C)_red` と `φ*ω_ℙ(C) ≅ ω_X(E)`（→ `h2`）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`thm_2_1_stepB` の `h2 : BDeq htOmXE htOmPC` がこれである。 -/
theorem pullback_omega_eq_omega_add_E {Pt : Type} (htOmXE htOmPC : Pt → ℝ)
    (hpull : ∀ x, htOmXE x = htOmPC x) :
    BDeq htOmXE htOmPC :=
  ⟨0, fun x => by rw [hpull x]; simp⟩

def pullback_omega_eq_omega_add_E.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 4——φ*ω_ℙ(C) ≅ ω_X(E) から h2)",
    sectionId := "genell-thm-2-1" }

def pullback_omega_eq_omega_add_E.needs : List ProofObligation :=
  [ .implicitStep "★引き戻し因子 `E = φ⁻¹(C)_red` と標準束の関係（代数幾何）。" 8 ]

/-- ★★★★★★**節点 5**——`ht_E ≳ (deg E / deg ω_X)·ht_{ω_X}`（→ `h6`）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`Prop 1.4` の高さの次数比例性である。 -/
theorem htE_ge_ratio_htOmega {Pt : Type} (htE htOm : Pt → ℝ) (r : ℝ)
    (hprop : ∀ x, htE x ≤ r * htOm x) :
    BDge htE (fun x => r * htOm x) :=
  ⟨0, fun x => by simpa using sub_nonpos.2 (hprop x)⟩

def htE_ge_ratio_htOmega.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 5——ht_E ≳ (deg E/deg ω_X)·ht_ω から h6)",
    sectionId := "genell-thm-2-1" }

def htE_ge_ratio_htOmega.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Proposition 1.4, (i)(iii)(高さの次数比例性)" 6,
    .implicitStep "☆`Found/GenEll/Prop14.lean` に条なし実装がある。" 3 ]

/-- ★★★★★★**節点 6**——`Prop 1.7, (i)` を `e = 1` で当てる（→ `h4`）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`log-diff_ℙ + log-cond_C ≳ log-diff_X + log-cond_E`。 -/
theorem logDiff_pullback_le {Pt : Type} (logDiffP logCondC logDiffX logCondE : Pt → ℝ)
    (hprop : ∀ x, logDiffP x + logCondC x ≤ logDiffX x + logCondE x) :
    BDge (fun x => logDiffP x + logCondC x) (fun x => logDiffX x + logCondE x) :=
  ⟨0, fun x => by simpa using sub_nonpos.2 (hprop x)⟩

def logDiff_pullback_le.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 6——Prop 1.7 (i) を e = 1 で当てて h4)",
    sectionId := "genell-thm-2-1" }

def logDiff_pullback_le.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Proposition 1.7, (i)(log-diff の引き戻し)" 9,
    .implicitStep "☆`Found/GenEll/Prop17.lean` に条なし実装がある。" 4 ]

/-! ## ★★★★節点 7-8 —— 段 A と組み立て -/

/-- ★★★★★★**節点 7**——段 A: `U_X(ℚ̄)^{≤d}` の点は `U_Y(ℚ̄)^{≤d′}` に持ち上がる。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆`d′ = d·deg(Y/X)` と置けばよい、というのが原文の言い分である。 -/
theorem exists_lift_of_degLe {A B : Type} (π : B → A) (hsurj : Function.Surjective π)
    (degA : A → ℕ) (degB : B → ℕ) (k : ℕ)
    (hdeg : ∀ y, degB y ≤ degA (π y) * k) (d : ℕ) :
    ∀ x, degA x ≤ d → ∃ y, π y = x ∧ degB y ≤ d * k := by
  intro x hx
  obtain ⟨y, rfl⟩ := hsurj x
  exact ⟨y, rfl, le_trans (hdeg y) (Nat.mul_le_mul_right k hx)⟩

def exists_lift_of_degLe.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 7——段 A の点の持ち上げ)",
    sectionId := "genell-thm-2-1" }

def exists_lift_of_degLe.needs : List ProofObligation :=
  [ .implicitStep "☆`d′ = d·deg(Y/X)` の定義と、被覆の全射性。" 3 ]

/-- ★★★★★★★★★★★★★★★★**節点 8**——`(ii) ⟹ (i)` の組み立て。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

★★★★これが `Theorem 2.1` の**実質**である。
☆節点 1–7 と、在庫の `thm_2_1_stepA`・`thm_2_1_stepB`・`exists_bad_seq_tendsto` を繋ぐ。 -/
theorem theorem_2_1_converse {Pt : Type} (htOm logDiffX : Pt → ℝ) (eps : ℝ)
    (KV : Set Pt)
    (hKV : BDge (fun x : ↥KV => htOm x.1) (fun x : ↥KV => (1 + eps) * logDiffX x.1)) :
    BDge htOm (fun x => (1 + eps) * logDiffX x) := by
  sorry

def theorem_2_1_converse.src : Source :=
  { paper := "GenEll", pdfPage := 11,
    item := "Theorem 2.1(節点 8——(ii) ⟹ (i) の組み立て。★これが実質)",
    sectionId := "genell-thm-2-1" }

def theorem_2_1_converse.needs : List ProofObligation :=
  [ .citation "[ABC3]" "thm_2_1_stepB(段 B の数値の鎖、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.thm_2_1_stepB") 1,
    .citation "[ABC3]" "thm_2_1_stepA(段 A の数値の鎖、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.thm_2_1_stepA") 1,
    .citation "[ABC3]" "exists_bad_seq_tendsto(反例列とコンパクト性、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_bad_seq_tendsto") 1,
    .citation "[ABC3]" "exists_belyi_noncritical_general(節点 1)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.exists_belyi_noncritical_general") 12,
    .implicitStep
      ("★★**進捗枠 0 / 8**（2026-09-01、第 1093 で設置）。" ++
       "☆節点 1 が最大（[NCBelyi] Thm 2.5 の一般曲線版、Riemann–Roch）。" ++
       "★数値の鎖（`thm_2_1_stepA`・`thm_2_1_stepB`）は既に `Found/` にあり `sorry` 無し。") 8 ]

end ABC3.Skeleton.GenEll
