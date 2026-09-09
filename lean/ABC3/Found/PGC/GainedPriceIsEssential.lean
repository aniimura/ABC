import ABC3.Found.PGC.StepwiseApproxFree
import ABC3.Found.PGC.GainedTowerDescent

/-!
# [pGC] ★★★`Gained*` の `∏` は**緩いのではなく本質** ＋ 私の前々波は木の重複だった

## 持ち場（前波で「次の 1 点」とした点）

前波の私の言葉（逐語）:

> `Gained*` の 1 層の値段が `∏_{j=1}^{k+1} axDecay p j` なのは**上界として緩いだけ**か、
> それとも本質かを測ること。★**`∏` ではなく `axDecay p (k+1)` 単項で押さえられるか**。

## ★★★答: **本質。単項には落ちない**（実データ 7 例すべてで破れる）

`GainedTowerDescent.lean:435 towerBudget_iff` が与える**積の条件**（逐語の形）:

  `p^{J/(p^k e)} ≤ ∏_{i ∈ Icc 1 k} axDecay p i  ↔  (p−1)²·J ≤ (p^{k+1} − p)·e`

**単項の条件**（本波で導出、§1 で `ℤ` の形を型にした）:

  `p^{J/(p^k e)} ≤ axDecay p k  ↔  (p−1)·J ≤ p·e`

（`axDecay p k = p^{(1/(p−1))p^{1−k}}` なので
`J/(p^k e) ≤ p^{1−k}/(p−1) ⟺ (p−1)J ≤ p^k·e·p^{1−k} = p·e`。）

★**2 つの閾値の比は `(p^k − 1)/(p−1) = 1 + p + ⋯ + p^{k−1}`**。
⇒ `k = 1` で一致、★`k ≥ 2` では積の条件が真に緩い。

### 測定（`tools/gained-single-vs-prod.py`、厳密整数、`gainedLoss` を定義どおり計算）

| 例 | `n` | `Λ = gainedLoss` | 積 `(p−1)²Λ ≤ (p^{n+1}−p)e` | 単項 `(p−1)Λ ≤ p·e` |
|---|---|---|---|---|
| `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` `e=2`, `t=(2,8)` | 2 | **8** | `32 ≤ 48` OK | `16 ≤ 6` ★**NG** |
| `ℚ₃(ζ₈₁)/ℚ₃(ζ₉)` `e=6`, `t=(8,26)` | 2 | 26 | `104 ≤ 144` OK | `52 ≤ 18` ★NG |
| `ℚ₃(ζ₈₁)/ℚ₃(ζ₃)` `e=2`, `t=(2,8,26)` | 3 | 28 | `112 ≤ 156` OK | `56 ≤ 6` ★NG |
| `p=5`, `n=2` / `n=3` | 2/3 | 24 / 128 | OK / OK | ★NG / ★NG |
| `p=2`, `n=2` / `n=3` | 2/3 | 3 / 8 | OK / OK | ★NG / ★NG |

★★**7 例すべてで単項は破れる**（`p ∈ {2,3,5}`、`n ≥ 2`）。★`t` の値は木の
`GainedTowerModel.lean` §10 `ChainCheck`（`u 0 = 2`, `u 1 = 8` ほか）から取った。

⇒ ★★★**前波の負の結果（`StepwiseApproxFree.prod_of_prod_unbounded`）と合わせると、
`Gained*` の 1 層あたりの値段を素朴に積む道は閉じない。**

## ★★在庫の測定 —— 私の前々波は木の重複だった（申告）

`DeepDescentPairDirect.lean:223`（逐語）:

```lean
def AxLemmaGraded (K : PAdicLocalField p) : Prop :=
  ∀ ε : ℝ, 0 ≤ ε → ∀ x : K.closure, (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
    ∃ a : K.carrier, ‖x - algebraMap K.carrier K.closure a‖
      ≤ (∏ i ∈ Finset.Icc 1 (wildDepth K x), axDecay p i) * ε
```

★★これは私の `DepthApproxToAxLemma.axLemma_axConstant_of_prod_axDecay_approx`（第 1131 波）の
**仮説そのもの**であり、出口 `axSenTate_of_axLemmaGraded`（`:306`）も既に在る。
⇒ ★**あのファイルは木の重複だった**（`prod_axDecay_le` で `axConstant p` に丸める点だけが差）。
★新規結果として数えないこと。★第 1128 波の `TraceGainKernel` に続き 2 度目である。

★★**測っていれば防げた。** `grep -rn "def AxLemmaGraded" lean/ABC3/Found/PGC/*.lean` は 0.3 秒。
★私は `axSenTate_of_axLemmaGraded` の**名前**は第 1131 で引用しておきながら、
★**定義を開かなかった**。★名前を引くときは定義も開く。

## ★これで分かった「本当のギャップ」

* `Gained*` : `M` → `F`（**1 層**、`[M:F] = p`）を値段 `∏_{j≤k+1} axDecay p j` で。
* `AxLemmaGraded` : `M` → **`K`（底まで）** を値段 `∏_{j≤wildDepth} axDecay p j` で。

★**値段は同じ形。違うのは「1 層」か「底まで」か**だけである。
★前波の `norm_sub_le_of_approx` で次の `ε` は `∏·ε` にしかならないので、
★素朴な帰納は積の積になって閉じない（前波 §2）。
⇒ ★★**残っているのは「1 層下りても `ε` を増やさない」形の帰納**である。

## 逸脱の記録

- §1 は `ℤ` の不等式だけで、★体も分岐も出ない。★`rpow` の同値（`towerBudget_iff`）は
  木に在るのでそちらを使う前提で、本ファイルは**閾値の比較だけ**を型にした。
- §2 の `gainedLoss` の値（`Λ = 8`）は ★**定義を展開して `norm_num` が計算した**もので、
  Python の計算とは独立である（`norm_num [gainedLoss, jumpSum, tZeta27]` の 1 行）。
- ★本ファイルは「単項に落ちない」ことを 7 例で示しただけで、★**すべての塔について
  落ちないことは証明していない**（`t` の一般形を仮定していない）。
-/

namespace ABC3.Found.PGC

namespace GainedPriceIsEssential

/-! ## §1 ★抽象核 —— 2 つの予算条件の比較（`ℤ` だけ） -/

section Budget

/-- ★**単項の条件は積の条件より真に強い**（`k ≥ 1`）。

積: `(p−1)²·J ≤ (p^{k+1} − p)·e`（木の `GainedDescent.towerBudget_iff`）
単項: `(p−1)·J ≤ p·e`（`p^{J/(p^k e)} ≤ axDecay p k` と同値）

★逆は `k = 1` を除いて成り立たない（§2 の反例）。 -/
theorem single_le_prod_budget {p J e : ℤ} {k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k)
    (he : 0 ≤ e) (h : (p - 1) * J ≤ p * e) :
    (p - 1) ^ 2 * J ≤ (p ^ (k + 1) - p) * e := by
  have hp1 : (0 : ℤ) ≤ p - 1 := by linarith
  have hstep : (p - 1) ^ 2 * J ≤ (p - 1) * (p * e) := by
    have := mul_le_mul_of_nonneg_left h hp1
    nlinarith
  have hpow : p * p ≤ p ^ (k + 1) := by
    calc p * p = p ^ 2 := by ring
      _ ≤ p ^ (k + 1) := pow_le_pow_right₀ (by linarith) (by omega)
  nlinarith

/-- ★2 つの閾値の比が `1 + p + ⋯ + p^{k−1}` であることの `k = 2` の姿。

`(p³ − p) = (p + 1)·((p − 1)·p)` である。★積の条件の右辺 `(p³−p)·e` を `(p−1)` で割ると
`(p+1)·p·e` で、単項の条件の右辺 `p·e` の ★**`p + 1` 倍**（`= 1 + p`）。 -/
theorem budget_ratio_two (p : ℤ) : (p ^ 3 - p) = (p + 1) * ((p - 1) * p) := by ring

end Budget

/-! ## §2 ★★反例 —— `ℚ₃(ζ₂₇)` の実データで単項は破れる -/

section CounterExample

open GainedDescent

/-- `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の跳び `t = (t₁,t₂) = (2,8)`（木の
`GainedTowerModel.lean` §10 `ChainCheck`: `u 0 = 2`, `u 1 = 8`）。 -/
def tZeta27 : ℕ → ℤ := fun j => if j = 1 then 2 else if j = 2 then 8 else 0

/-- ★`gainedLoss 3 tZeta27 2 = 8`（★定義を展開して機械が計算した。
`tools/gained-single-vs-prod.py` の値と独立に一致）。 -/
theorem gainedLoss_zeta27 : gainedLoss (3 : ℤ) tZeta27 2 = 8 := by
  norm_num [gainedLoss, jumpSum, tZeta27]

/-- ★**積の条件は成り立つ**: `(3−1)²·8 = 32 ≤ (3³−3)·2 = 48`。 -/
theorem prod_budget_holds_zeta27 :
    ((3 : ℤ) - 1) ^ 2 * 8 ≤ ((3 : ℤ) ^ (2 + 1) - 3) * 2 := by norm_num

/-- ★★★**単項の条件は破れる**: `(3−1)·8 = 16 > 3·2 = 6`。

⇒ ★`p^{Λ/(p²·e)} ≤ axDecay 3 2` は**偽**。
★`Gained*` の 1 層の値段は `axDecay p k` 単項には落ちない。 -/
theorem single_budget_fails_zeta27 : ¬ (((3 : ℤ) - 1) * 8 ≤ 3 * 2) := by norm_num

/-- ★★`k = 2` で「単項 ⇒ 積」は使えても逆は使えないことの証拠
（`J = 8`, `e = 2`, `p = 3` は積を満たし単項を満たさない）。 -/
theorem prod_not_imp_single :
    (((3 : ℤ) - 1) ^ 2 * 8 ≤ ((3 : ℤ) ^ (2 + 1) - 3) * 2) ∧ ¬ (((3 : ℤ) - 1) * 8 ≤ 3 * 2) :=
  ⟨prod_budget_holds_zeta27, single_budget_fails_zeta27⟩

end CounterExample

/-! ## §3 使っている公理の一覧 -/

#print axioms single_le_prod_budget
#print axioms budget_ratio_two
#print axioms gainedLoss_zeta27
#print axioms prod_budget_holds_zeta27
#print axioms single_budget_fails_zeta27
#print axioms prod_not_imp_single

end GainedPriceIsEssential

end ABC3.Found.PGC
