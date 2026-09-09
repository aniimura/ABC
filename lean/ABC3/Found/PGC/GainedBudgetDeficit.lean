import ABC3.Found.PGC.GainedPriceIsEssential
import ABC3.Found.PGC.WildDescentMultiStep

/-!
# [pGC] ★★★ギャップは「値段」ではなく **着地の深さ** —— 過払いの因子は `∏_{Icc 1 d'}` ちょうど

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> 「1 層下りても `ε` を増やさない」形の帰納が木に在るか。…
> ★`axLemmaGraded_of_multi` が **`∏` を増やさずに底まで降ろしている**はずなので、
> ★**その証明が `ε` をどう扱っているかを読むこと**。★まだ読んでいません。

## ★★読んだ結果（3 点、いずれも逐語）

### (1) `AxWildDescentMulti` は**望遠鏡的**（`WildDescentMultiStep.lean:233`）

```lean
def AxWildDescentMulti (K : PAdicLocalField p) (c : ℕ → ℝ) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    p ∣ (minpoly K.carrier x).natDegree →
      ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧
        ‖x - x'‖ ≤ (∏ k ∈ Finset.Icc (wildDepth K x' + 1) (wildDepth K x), c k) * ε
```

★1 歩の値段は ★**飛ばした段だけの積** `∏_{Icc (d'+1) d}` である。
★互いに素な区間なので `d → 0` まで合成すると `∏_{Icc 1 d}` に**畳まれる**（積の積にならない）。

### (2) `ε` の伸びは木でも**無料**（同ファイル `:233` の docstring、逐語）

> * `ε` の伸びの条件を**書かない**(`UltraCore` の抽象核 1 でただで出るため)。

★実体は `UltraCore.norm_smul_sub_self_le_of_norm_sub`（`WildDescentMultiStep.lean:142`）。

### (3) 三者は同値（`DeepDescentPairDirect.lean:295/:301`）

`AxDeepDescentPair` ⟺ `AxWildDescentMulti K (axDecay p)` ⟺ `AxLemmaGraded`、
そして `axSenTate_of_axLemmaGraded`（`:306`）で `AxSenTate` に出る。

## ★★★本波の答え —— ギャップは「値段」ではなく **着地の深さ**

* `AxWildDescentMulti` が許すのは、深さ `d` から `d'` へ下りて `∏_{Icc (d'+1) d} axDecay p k`。
* `Gained*` が払うのは `∏_{Icc 1 (k+1)} axDecay p j`（`k+1 = d`）で、
  ★**着地は `F = twr g p k`、すなわち深さ `d' ≤ d − 1`（`0` ではない）**。

⇒ ★★`Gained*` の値段は「**深さ 0 まで下りる**ときの値段」である。
★★★**過払いの因子はちょうど `∏_{Icc 1 d'} axDecay p k`**（§2 `gained_overpay_factor`、等式）。
⇒ `d' = 0` なら過払いは無い。★`d' ≥ 1` なら値段は予算を**真に超える**（§2）。

★★**ゆえに、`Gained*` を `AxWildDescentMulti` に嵌めるのに要るのは「値段を下げる」ことではなく
「着地を `K`（深さ 0）まで持っていく」こと**である。★そしてそれは `AxLemmaGraded` そのもの。

## ★★★★申告 —— 私の 3 波が木の重複だった（3 件目）

| 波 | 私のファイル | 木の既存 |
|---|---|---|
| 第 1128 | `TraceGainKernel` §1–§3 | `CyclicLayerDescent.lean:412/:461`（★木の方が良い評価） |
| 第 1131 | `DepthApproxToAxLemma` | `DeepDescentPairDirect.lean:223 AxLemmaGraded` ＋ `:306` |
| 第 1132 | `StepwiseApproxFree` §1・§3 | ★**`UltraCore.norm_smul_sub_self_le_of_norm_sub`**（`WildDescentMultiStep.lean:142`） |

★第 1132 で私は「木の `WildDepthDescent:517` は第 3 条件を 6 行かけて証明している。
本補題があればそれは要らない」と書いたが、★**木は既に `UltraCore` に一般形を持っていた**。
★私の形（`f : M → M` が加法的かつ等長、群作用を要らない）は木の形（モノイド作用 `G`）と
**互いに一般化の向きが違う**が、★**新規性は無い**。★新規結果として数えないこと。

### ★根本原因を測った（予防 4 本目）

木には ★**144 個の名前空間**がある（`grep -rhn "^namespace " lean/ABC3/Found/PGC/*.lean |
sed 's/.*namespace //' | sort -u | wc -l`）。そのうち `UltraCore` / `ProdCore` / `PairLedger` /
`JumpArith` / `GainedDescent` などは ★**抽象核だけを置く名前空間**である。
★私は毎回「自分の言葉で書いた結論の形」で grep しており、★**抽象核の置き場を見ていなかった**。

★★**予防（4 本目）: 抽象核を書く前に、`grep -rn '抽象核' lean/ABC3/Found/PGC/*.lean | head -40`
で「抽象核」と自称している節を先に見る。**（既存の 3 本: 数は `wc -l`／grep は実行してから書く／
名前を引くときは定義も開く。）★本波はこれを**書く前に**実行した。

## 逸脱の記録

- §1・§2 は `Finset` の積の分割だけで、★体も分岐も `axDecay` の中身も使わない
  （`1 ≤ axDecay` のみ）。★過払いの**等式**なので、どちらが大きいかの議論に依らない。
- ★本ファイルは `Gained*` の着地を深くする方法を**与えていない**。★与えられるとも書いていない。
-/

namespace ABC3.Found.PGC

namespace GainedBudgetDeficit

/-! ## §1 ★抽象核 —— `Icc 1 d` を `d'` で切る -/

section Split

variable {R : Type*} [CommMonoid R]

/-- ★`∏_{Icc 1 d} f = (∏_{Icc 1 d'} f) · (∏_{Icc (d'+1) d} f)`（`d' ≤ d`）。

★`Finset.Icc 1 n = Finset.Ioc 0 n` に直してから `Finset.prod_Ioc_consecutive`
（`lean-idioms` #302 の最後の段と同じ直し方）。 -/
theorem prod_Icc_split_at (f : ℕ → R) {d' d : ℕ} (h : d' ≤ d) :
    (∏ k ∈ Finset.Icc 1 d', f k) * (∏ k ∈ Finset.Icc (d' + 1) d, f k)
      = ∏ k ∈ Finset.Icc 1 d, f k := by
  have e1 : Finset.Icc 1 d' = Finset.Ioc 0 d' := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have e2 : Finset.Icc (d' + 1) d = Finset.Ioc d' d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  have e3 : Finset.Icc 1 d = Finset.Ioc 0 d := by
    ext k; simp only [Finset.mem_Icc, Finset.mem_Ioc]; omega
  rw [e1, e2, e3]
  exact Finset.prod_Ioc_consecutive f (Nat.zero_le d') h

end Split

/-! ## §2 ★★過払いの因子は `∏_{Icc 1 d'}` ちょうど -/

section Deficit

variable {p : ℕ} [Fact p.Prime]

omit [Fact p.Prime] in
/-- ★★★**過払いの因子の等式**。

`Gained*` が払う `∏_{Icc 1 d}` と `AxWildDescentMulti` が許す `∏_{Icc (d'+1) d}` の比は
★ちょうど `∏_{Icc 1 d'} axDecay p k` である。

⇒ ★`d' = 0` なら比は `1`（過払い無し）。★`d' ≥ 1` なら比は `> 1`。 -/
theorem gained_overpay_factor {d' d : ℕ} (h : d' ≤ d) :
    (∏ k ∈ Finset.Icc 1 d', axDecay p k) * (∏ k ∈ Finset.Icc (d' + 1) d, axDecay p k)
      = ∏ k ∈ Finset.Icc 1 d, axDecay p k :=
  prod_Icc_split_at (axDecay p) h

/-- ★★**`d' ≥ 1` なら `Gained*` の値段は `AxWildDescentMulti` の予算を真に超える**。

⇒ ★★★`Gained*`（着地は `F`、深さ `d' = d − 1 ≥ 1`）は、`d ≥ 2` では
`AxWildDescentMulti K (axDecay p)` を**満たさない**。
★足りないのは「値段を下げる」ことではなく「**着地を深さ 0（`K`）まで持っていく**」ことである。 -/
theorem gained_price_gt_budget {d' d : ℕ} (hd' : 1 ≤ d') (h : d' ≤ d) :
    (∏ k ∈ Finset.Icc (d' + 1) d, axDecay p k) < ∏ k ∈ Finset.Icc 1 d, axDecay p k := by
  have hsplit := gained_overpay_factor (p := p) h
  have hA : (1 : ℝ) < ∏ k ∈ Finset.Icc 1 d', axDecay p k := by
    have hmem : (1 : ℕ) ∈ Finset.Icc 1 d' := by
      simp only [Finset.mem_Icc]; omega
    have hrest : (1 : ℝ) ≤ ∏ k ∈ (Finset.Icc 1 d').erase 1, axDecay p k :=
      Finset.one_le_prod (fun i _ => one_le_axDecay p i)
    have h1 : (1 : ℝ) < axDecay p 1 := WildStepFieldSupply.one_lt_axDecay_one
    have := Finset.mul_prod_erase (Finset.Icc 1 d') (fun k => axDecay p k) hmem
    nlinarith [this]
  have hB : (1 : ℝ) ≤ ∏ k ∈ Finset.Icc (d' + 1) d, axDecay p k :=
    Finset.one_le_prod (fun i _ => one_le_axDecay p i)
  nlinarith [hsplit]

omit [Fact p.Prime] in
/-- ★**着地が深さ 0 のときだけ値段と予算が一致する**（等式の側）。 -/
theorem budget_eq_of_landing_zero (d : ℕ) :
    (∏ k ∈ Finset.Icc (0 + 1) d, axDecay p k) = ∏ k ∈ Finset.Icc 1 d, axDecay p k := by
  norm_num

end Deficit

/-! ## §3 使っている公理の一覧 -/

#print axioms prod_Icc_split_at
#print axioms gained_overpay_factor
#print axioms gained_price_gt_budget
#print axioms budget_eq_of_landing_zero

end GainedBudgetDeficit

end ABC3.Found.PGC
