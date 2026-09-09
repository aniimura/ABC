import ABC3.Found.PGC.TraceGainKernel

/-!
# [pGC] ★★測定: `Gained*` 族は出口に**一度も繋がっていない**（0/26, 0/23）＋ 繋ぎ金具

## 持ち場（前波で「次に測るべきこと」とした点）

前波の私の言葉（逐語）:

> `Gained*` 族 26 件のどれが `AxLemmaGraded` / `AxSenTate` まで繋がっているか。
> ★私の `AxWildDescent` 側の 3 本とどちらが安いかも**測っていません**。

## ★★★測定の結果（`tools/gained-family-reach2.mjs`、推移閉包・両方向）

| 測ったもの | 値 |
|---|---|
| `Gained*` 族（`^theorem exists_norm_sub_algebraMap_le_prod_axDecay*`） | ★**26 件 / 14 ファイル** |
| 結論に `AxSenTate` / `AxLemmaGraded` を持つ theorem | ★**23 件** |
| 出口から**推移的に** `Gained*` 族へ届くもの | ★★**0 / 23** |
| `Gained*` 族から**推移的に**出口へ届くもの | ★★**0 / 26** |

★★**`Gained*` 族は木の出口と完全に切り離されている。** 両方向とも 0 である。

★健全性検査（script が 0 を返すのがバグでないことの確認、同スクリプト末尾）:
既知の道 4 本（`axSenTate_of_deep_bound ⇝ axLemma_of_wildDescent_Icc` など）は
すべて `YES` を返す。⇒ ★0 は script の欠陥ではない。

### ★ファイル数の訂正（自己訂正 17 度目）

前波で私は「26 件 / **16 ファイル**」と書いた。★**14 が正しい**（本体の測定と一致）。
★私は `grep` の 26 行の出力を**目で見て**ファイル名を数えていた。
★★**予防**: 数は必ずコマンドの出力そのものを使う（`| wc -l`）。目視で数えない。
本波は `grep -rln ... | wc -l` → `14` を実行してから書いた。

## ★測定の解釈 —— 2 つは競合ではなく**相補**

* `Gained*` 族の結論は `∃ y : F, ‖x − algebraMap F M y‖ ≤ (∏ axDecay)·‖τ x − x‖`
  ——★**1 つの層 `M/F`（`[M:F] = p`）の中の `x`** についての主張。
* 出口 `AxWildDescent` が要るのは `x : K.closure`（**任意の深さ**）から `x'` を作ること。

⇒ ★**どちらが安いかという問いは立たない。** 欠けているのは**繋ぎ金具**である。

仮説の数（同スクリプト、★section variable は数に入らないので**両側とも下限**）:

| 宣言 | 仮説（下限） |
|---|---|
| `DeepLossExit.axSenTate_of_cyclotomic_loss`（私） | **5** |
| `BaseLayerConstant.axSenTate_of_deep_bound`（私） | **5** |
| `DeepDescentPairDirect.axSenTate_of_axLemmaGraded`（木） | **2** |
| `Gained*` 族のうち最小（`PureStepSetup.lean:422`） | 2 |
| `Gained*` 族の中央値あたり（`GainedDescentBridge.lean:527`） | **26** |
| `Gained*` 族の最大（`TotallyRamifiedValueGroup.lean:430`） | **33** |

## ★本ファイルが足す繋ぎ金具

★★**「深さ `k` の元が `K` の元で `c k · ε` 以内に近似できる」なら `AxWildDescent K c`**（§1）。

`Gained*` 族の結論はまさにこの形（`F = K.carrier` と取れば）なので、
★**族を出口に繋ぐのに要るのは「`F` を `K` に取り、深さと層の次数を合わせる」ことだけ**である。

§2 は「深さ 1 だけ改善する」版。★`Gained*` 族は `[M:F] = p` の**1 層**しか扱わないので、
実際に繋げられるのは深さ 1 の段である。他の深さは在庫の `axWildDescent_prime`（定数 `p`）で埋める。

## 逸脱の記録

- §1 は `AxWildDescent` の定義を展開して `x' := algebraMap a` を入れるだけで、
  ★新しい数学は 0。★価値は「族の結論の形が出口の形に**そのまま嵌まる**」ことを型にした点。
- ★測定スクリプトの既知の限界: 同名の宣言が 2 つあると先に見つけた方だけを見る
  （`exists_norm_sub_algebraMap_le_prod_axDecay` は `PureStepSetup.lean:422` と
  `GainedDescentBridge.lean:527` の 2 つがある）。★0/26 の結論は族の**全 26 名**について
  取っているのでこの限界の影響を受けない。
-/

namespace ABC3.Found.PGC

namespace GainedToWildStep

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 ★繋ぎ金具 —— 「各深さで `K` の元に近づける」⇒ `AxWildDescent` -/

section Connector

/-- ★★**深さごとの `K`-近似から `AxWildDescent` が出る**。

`happ` : 深さ `k ≥ 1` の元 `x` が `∀σ ‖σx − x‖ ≤ ε` を満たすなら
`‖x − a‖ ≤ c k · ε` なる `a ∈ K` が取れる。

⇒ `AxWildDescent K c`。★`x' := a` は `K` の元なので

* `wildDepth K x' = 0 < wildDepth K x`（`wildDepth_algebraMap`）、
* `∀σ ‖σx' − x'‖ = 0`（★`ε` の伸びは**ゼロ**）

が自動である。

★★`Gained*` 族（26 件）の結論 `∃ y : F, ‖x − algebraMap F M y‖ ≤ C·‖τ x − x‖` は
`F = K.carrier` と取ればちょうどこの形である。 -/
theorem axWildDescent_of_depth_approx (K : PAdicLocalField p) {c : ℕ → ℝ}
    (hc : ∀ k, 0 ≤ c k)
    (happ : ∀ (k : ℕ), 1 ≤ k → ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure,
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) → wildDepth K x = k →
      ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖ ≤ c k * ε) :
    AxWildDescent K c := by
  intro ε hε x hx hdvd
  have hk1 : 1 ≤ wildDepth K x :=
    Nat.one_le_iff_ne_zero.mpr (fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd)
  obtain ⟨a, ha⟩ := happ (wildDepth K x) hk1 ε hε x hx rfl
  refine ⟨algebraMap K.carrier K.closure a, ?_, ha, ?_⟩
  · rw [wildDepth_algebraMap]
    omega
  · intro σ
    rw [smul_closure_def, AlgEquiv.commutes, sub_self, norm_zero]
    exact mul_nonneg (hc _) hε

end Connector

/-! ## §2 深さ 1 だけ改善する版（`Gained*` 族は 1 層しか扱わない） -/

section DepthOne

/-- ★★**深さ 1 だけ `C` に改善し、他の深さは在庫の `p` で埋める**。

★`Gained*` 族が扱うのは `[M:F] = p` の **1 層**なので、実際に繋げられるのはここである。
他の深さは `WildDepthFieldDescent.axWildDescent_prime`（定数 `p`、無条件）で埋まる。 -/
theorem axWildDescent_of_depth_one_approx (K : PAdicLocalField p) {C : ℝ} (hC : 0 ≤ C)
    (happ : ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      wildDepth K x = 1 → ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖ ≤ C * ε) :
    AxWildDescent K (fun k => if k = 1 then C else (p : ℝ)) := by
  intro ε hε x hx hdvd
  by_cases hd : wildDepth K x = 1
  · obtain ⟨a, ha⟩ := happ ε hε x hx hd
    refine ⟨algebraMap K.carrier K.closure a, ?_, ?_, ?_⟩
    · rw [wildDepth_algebraMap, hd]
      omega
    · simpa [hd] using ha
    · intro σ
      rw [smul_closure_def, AlgEquiv.commutes, sub_self, norm_zero]
      simp only [hd, if_pos]
      exact mul_nonneg hC hε
  · obtain ⟨x', hlt, h1, h2⟩ := axWildDescent_prime K ε hε x hx hdvd
    exact ⟨x', hlt, by simpa [hd] using h1, fun σ => by simpa [hd] using h2 σ⟩

end DepthOne

/-! ## §3 使っている公理の一覧 -/

#print axioms axWildDescent_of_depth_approx
#print axioms axWildDescent_of_depth_one_approx

end GainedToWildStep

end ABC3.Found.PGC
