import ABC3.Found.PGC.WildStepFieldSupply

/-!
# [pGC] ★★★`k ≥ 1` は在った ——「深さごとの `K`-近似」から `AxSenTate` が**直接**出る

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> `k ≥ 1` の層（wild 深さ 2 以上）を扱う `Gained*` 族の members が**存在するか**を測ること。
> ★**それらの `k` が wild 深さと一致するかは測っていません**。

## ★★★測定の結果（読んで確かめた。見込みではない）

**答: 存在する。そして `k + 1` が wild 深さである。**

`GainedTowerStep.lean:404 exists_norm_sub_algebraMap_le_prod_axDecay_of_tower_jumps` の
docstring（逐語）:

> ★入力は実際の跳びの列 `u` と塔のデータだけで、`τ : M →+ M`、**`1 ≤ k` 可**。

`k` の意味は、同ファイル `:366 _of_tower` の仮説を読んで確定した:

* `hdegj : ∀ j, j < k → (minpoly (E j) π).natDegree = p ^ (k + 1 - j)`
* `heM   : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e)`
* `hvalj : ∀ j, j < k → ∀ d : E j, d ≠ 0 → ∃ m, ‖algebraMap (E j) M d‖ = ‖π‖ ^ (p^{k+1-j} * m)`

⇒ `[M : F] = p^{k+1}` で全分岐（`e_M = p^{k+1}·e_F`）。
⇒ ★★`M = K(x)`, `F = K` と取れば **`k + 1 = v_p([K(x):K]) = wildDepth K x`**。

★★★**前波の「深さ 1 だけでは足りない」（`prod_depth_one_only_unbounded`）で塞がった道は、
`PureStepSetup`（`k = 0` 専用）の道だけだった。** `GainedTowerStep` / `GainedTowerModel` の
`_of_tower*` 系は `k` が動く。

仮説の数（`sed -n ... | grep -oE '^\s+[({\[]' | wc -l`、★目視で数えていない）:
`_of_tower` = **19**、`_of_tower_jumps` = **19**。

## ★★★本ファイルの結論 —— `AxWildDescent` を**通らない**出口

`Gained*` 族の結論は

  `∃ y : F, ‖x − algebraMap F M y‖ ≤ (∏ j ∈ Finset.Icc 1 (k+1), axDecay p j) · ‖τ x − x‖`

で、★積は**その層 1 つ**についてすでに閉じた形（`∏_{j=1}^{深さ}`）になっている。
⇒ ★★**段を積む必要が無い。** `F = K.carrier` と取り、`‖τ x − x‖ ≤ ε` を使えば

  `∃ y : K, ‖x − y‖ ≤ (∏_{j=1}^{wildDepth K x} axDecay p j)·ε ≤ axConstant p · ε`

で、これは **`AxLemma K (axConstant p)` そのもの**である（§2）。

⇒ ★★★**`AxWildDescent` も `axLemma_of_wildDescent_Icc` も基底段の別扱い（`K₀`）も要らない。**
★私が第 1125〜1130 で積んだ `AxWildDescent` 側の 5 本
（`FiniteExceptionsIcc` / `BaseLayerConstant` / `DeepLossExit` / `GainedToWildStep` /
`WildStepFieldSupply`）は、★**この道では使わない**。
★それらは「1 段ずつ降りる」道のためのもので、`Gained*` は「一気に降りる」道である。

★上界 `∏_{Icc 1 k} axDecay p j ≤ axConstant p` は私の
`BaseLayerConstant.prod_axDecay_le`（第 1126 波）を使う。★あの波の切り出しがここで効いた。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `axLemma_of_depth_bounds` | ★**抽象核**（深さごとの上界に一様な蓋があれば `AxLemma`） |
| §2 | `axLemma_axConstant_of_prod_axDecay_approx` | ★★`d k = ∏_{Icc 1 k} axDecay p j` を代入 |
| §2 | `axSenTate_of_prod_axDecay_approx` | ★★★出口 |
| §3 | `depth_bound_le_axConstant` | 蓋が `axConstant p` であることの明示 |

## ★残っているちょうど 1 点

`happ` —— すなわち「`x : K.closure` の wild 深さ `k` に対し、
`(∏_{j=1}^{k} axDecay p j)·ε` 以内で `K` の元に近づけること」。
★`GainedTowerStep._of_tower_jumps` の 19 仮説を `F := K.carrier`, `M := K(x)` で供給すれば
これになる。★19 のうち何本が `PAdicLocalField` から自動かは**まだ測っていない**。

## 逸脱の記録

- §1 は `AxLemma` の定義（`AxLemma.lean:248`）を展開して深さで場合分けするだけで、
  ★新しい数学は 0。★価値は「`Gained*` 族の結論の形が出口に**直接**嵌まる」ことを型にした点。
- ★本ファイルは `happ` を**供給していない**。★供給できるとも書いていない。★形を測っただけである。
- `AxLemma K C` の `C` は `axConstant p` ちょうどで、★**`p^{K₀}` の余分が付かない**
  （`AxWildDescent` 側の道では付いていた。`BaseLayerConstant.axLemma_of_deep_bound` 参照）。
-/

namespace ABC3.Found.PGC

namespace DepthApproxToAxLemma

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 ★抽象核 —— 深さごとの上界に一様な蓋があれば `AxLemma` -/

section Kernel

/-- ★★**深さごとの `K`-近似から `AxLemma` が直接出る**。

`happ` : wild 深さ `k` の元 `x` が `∀σ ‖σx − x‖ ≤ ε` を満たすなら
`‖x − a‖ ≤ d k · ε` なる `a ∈ K` が取れる。
`hd` : `d k ≤ C`（★深さについて**一様な蓋**）。

⇒ `AxLemma K C`。

★★`AxWildDescent` を**経由しない**。段を積まないので基底段の別扱いも要らない。
★`Gained*` 族（`GainedTowerStep.lean:366/:404` ほか）の結論はちょうど `happ` の形である。 -/
theorem axLemma_of_depth_bounds (K : PAdicLocalField p) {C : ℝ} {d : ℕ → ℝ}
    (hd : ∀ k, d k ≤ C)
    (happ : ∀ (k : ℕ) (ε : ℝ), 0 ≤ ε → ∀ x : K.closure, wildDepth K x = k →
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖ ≤ d k * ε) :
    AxLemma K C := by
  intro ε hε x hx
  obtain ⟨a, ha⟩ := happ (wildDepth K x) ε hε x rfl hx
  exact ⟨a, ha.trans (mul_le_mul_of_nonneg_right (hd _) hε)⟩

end Kernel

/-! ## §2 ★★★出口 —— `Gained*` 族の結論の形をそのまま入れる -/

section Exit

/-- 蓋は `axConstant p`（★私の `BaseLayerConstant.prod_axDecay_le`、第 1126 波）。 -/
theorem depth_bound_le_axConstant (k : ℕ) :
    (∏ j ∈ Finset.Icc 1 k, axDecay p j) ≤ axConstant p :=
  BaseLayerConstant.prod_axDecay_le p k

/-- ★★★**`Gained*` 族の結論を深さごとに供給すると `AxLemma K (axConstant p)`**。

★`p^{K₀}` の余分が付かない（`AxWildDescent` 側の道では付いていた）。 -/
theorem axLemma_axConstant_of_prod_axDecay_approx (K : PAdicLocalField p)
    (happ : ∀ (k : ℕ) (ε : ℝ), 0 ≤ ε → ∀ x : K.closure, wildDepth K x = k →
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖
        ≤ (∏ j ∈ Finset.Icc 1 k, axDecay p j) * ε) :
    AxLemma K (axConstant p) :=
  axLemma_of_depth_bounds K depth_bound_le_axConstant happ

/-- ★★★★**そこから `AxSenTate K`**。

★`AxWildDescent` も `axLemma_of_wildDescent_Icc` も基底段の別扱いも通らない。 -/
theorem axSenTate_of_prod_axDecay_approx (K : PAdicLocalField p)
    (happ : ∀ (k : ℕ) (ε : ℝ), 0 ≤ ε → ∀ x : K.closure, wildDepth K x = k →
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
      ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖
        ≤ (∏ j ∈ Finset.Icc 1 k, axDecay p j) * ε) :
    AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant p))
    (axLemma_axConstant_of_prod_axDecay_approx K happ)

end Exit

/-! ## §3 使っている公理の一覧 -/

#print axioms axLemma_of_depth_bounds
#print axioms depth_bound_le_axConstant
#print axioms axLemma_axConstant_of_prod_axDecay_approx
#print axioms axSenTate_of_prod_axDecay_approx

end DepthApproxToAxLemma

end ABC3.Found.PGC
