import ABC3.Found.PGC.AxDecayKnobCircular
import ABC3.Found.PGC.WildDescentDistanceOnly

/-!
# [pGC] 予算関数 `F` への載せ替えは**得をしない** —— 縛りは `c` ではなく `g` にある

## ★配られた問いへの答え（3 つとも測った）

**(1) `exists_mem_of_descent_budget`（`AxTowerDecay.lean:327`）の字面を読んだ。**

```
hF1 : ∀ d, 1 ≤ F d          hg0 : ∀ d, 0 ≤ g d
hFc : ∀ d, d ≠ 0 → c d ≤ F d
hFg : ∀ d d', d' < d → g d * F d' ≤ F d
⇒ ∀ x ε, 0 ≤ ε → P ε x → ∃ y ∈ S, ‖x − y‖ ≤ F (deg x) * ε
```

★**`c` にかかる条件は `c d ≤ F d` だけ**である。★`F` を束ねているのは `hFg`、すなわち
★**`g`（`ε` の伸び）だけ**である。

**(2) ★★★4 つの穴のうち②（一様定数）の射程を測り直した —— ★自分の結果の訂正。**

* ★`c ≡ C`（**一様定数の損失**）は予算形では**許される**:
  `g d ≤ 1` なら `F ≡ C` が 3 条件を全部満たす（§1 `budget_const_of_growth_le_one`）。
  ⇒ ★★**前波までの②（`const_descent_no_go`）は「損失 `c` が一様定数では届かない」
  という主張として読んではいけない。** それは `g ≡ c` と置いた形
  （`exists_mem_of_descent_prod` ＝ `axLemma_of_wildDescent`）にしか当たらない。
* ★逆に、★**②の内容は `g` に移る**: `C ≤ g d`（`1 < C`）が一様に成り立つなら
  `C^d ≤ F d`（§1 `pow_le_budget`）で `F` は非有界（`not_bddAbove_budget`）。
  ⇒ ★★**縛りは損失ではなく `ε` の伸びである。**

**(3) ★★★しかし載せ替えても得をしない —— `g ≤ 1` の予算形は `AxLemma` と同値（循環）。**

§2 で予算形の `hstep`（`c ≡ C`, `g ≡ θ`）が `AxWildDescentDecay K C θ` そのもの
であることを示し（`axWildDescentDecay_iff_budgetStep`。違いは `wildDepth ≠ 0` と
`p ∣ deg` の言い換えだけ）、前波の `AxKnob.axWildDescentDecay_iff_axLemma` に繋いだ:

★★`budgetStep_iff_axLemma` —— `1 ≤ C`・`0 ≤ θ ≤ 1` で **予算形 ⟺ `AxLemma K C`**。
★極端な場合 `θ = 0`（`ε` が 1 段で **0 に落ちる**）でも同値である
（`budgetStep_zero_iff_axLemma`）。

⇒ ★★★**「`c` と `g` を分けて `F` で束ねる」という載せ替えは、`g ≤ 1` を狙う限り
前波で塞いだ循環に戻る。** ★得をするとすれば「`g d > 1` の段が有限個だけ」という
中間の形だが、それも `F` が有界なら同じ議論で `AxLemma K (sup F)` と同値になる
（★本波では `sup F` の形は形式化していない）。

## ★`g d < 1`（深い段で `ε` が減る）を測った

`AxTowerDecay.lean:80-84` は `g d < 1` の理由として
「激しく分岐した拡大では**跡が小さい**ので `σ^p x − x = Tr(σx − x)` が `σx − x` より
真に小さくなる」と書いている。★測った結果は次のとおり。

* ★★**別の量である。** そこで小さくなるのは `σ^p x − x`（**同じ `x`**、`σ` の**冪**）だが、
  予算形が要求するのは `∀ σ ∈ Γ_K` についての `‖σ • x′ − x′‖`（**新しい `x′`**、
  **すべての** `σ`）である。★`Γ_K` は `σ` 自身も含むので、`σ` の冪に限った評価では届かない。
* ★超距離からの評価は `g ≤ max(c, 1)` しか与えない（§2 `growth_le_max_of_dist`、
  中身は `WildDescentDistanceOnly.lean:161 norm_smul_sub_self_le_max`）。
  ⇒ ★**`g d < 1` を得るには `c d < 1`（`ε` より真に良い近似）が要る。**
  ★これは「1 段で `ε` より良くなる」という主張で、`AxLemma K C` の `C < 1` に近い。

## ★証拠 5 の「食い違い」を測った —— ★食い違っていない

* `NormalizedTraceDescent.lean:262 traceLoss_eq_sharpLoss_rpow`: 跡は指数を `(p−1)` 倍にする
  ⇒ **損失** `c` にとっては**悪い**。
* `AxTowerDecay.lean:83`: 跡は小さい ⇒ **伸び** `g` にとっては**良い**。

★**同じ `(p−1)` 倍が、損失では悪く、伸びでは良い。向きが逆なので矛盾しない。**
⇒ ★木の 2 か所の記述は両立する（★本波の測定。どちらの docstring も書き換えていない）。

## ★4 つの穴の現状（本波の再測定）

| 穴 | 予算形の下で |
|---|---|
| ①不分岐（`TotallyRamifiedLayer.lean:342`） | ★生きている（分岐の話で、予算の形と無関係） |
| ②一様定数（`PGroupDescentToAxWild.const_descent_no_go`） | ★★**射程を訂正**: `c` ではなく `g` についての主張になる（§1） |
| ③幾何減衰（`JumpGeometricDecay.no_uniform_geometric_of_harith`） | ★生きている（跳びの変位の比の話で、`c`/`g` の分離と無関係） |
| ④測定点（`FirstJumpRouteEquiv.not_exists_firstJump_data_cyclotomic`） | ★★`c k ≤ axDecay p k` を前提にしている。予算形では `c` は `F` で抑えればよいので★**当たらない可能性がある**（★本波では確定していない） |

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ項目（pGC 物理 p.6 Corollary 3.1）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 は `ℝ` と `ℕ` だけで、分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。
4. ★`sup F` を使う「`g d > 1` が有限個」の形は**形式化していない**（上に明記）。
-/

namespace ABC3.Found.PGC

namespace BudgetRecast

/-! ## §1 抽象核 —— 予算 `F` を縛るのは `g` であって `c` ではない -/

section Core

/-- ★★抽象核: `g` が一様に `C ≥ 1` 以上なら予算は `C^d` 以上に膨らむ。

`exists_mem_of_descent_budget`（`AxTowerDecay.lean:327`）の `hFg`/`hF1` だけを使う。
★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。 -/
theorem pow_le_budget {g F : ℕ → ℝ} {C : ℝ} (hC : 1 ≤ C) (hF1 : ∀ d, 1 ≤ F d)
    (hg : ∀ d, C ≤ g d) (hFg : ∀ d d', d' < d → g d * F d' ≤ F d) :
    ∀ d : ℕ, C ^ d ≤ F d := by
  intro d
  induction d with
  | zero => simpa using hF1 0
  | succ m ih =>
      have h1 : g (m + 1) * F m ≤ F (m + 1) := hFg (m + 1) m (by omega)
      have h2 : C * F m ≤ g (m + 1) * F m :=
        mul_le_mul_of_nonneg_right (hg (m + 1))
          (le_trans zero_le_one (hF1 m))
      have h3 : C * C ^ m ≤ C * F m :=
        mul_le_mul_of_nonneg_left ih (le_trans zero_le_one hC)
      calc C ^ (m + 1) = C * C ^ m := by ring
        _ ≤ C * F m := h3
        _ ≤ g (m + 1) * F m := h2
        _ ≤ F (m + 1) := h1

/-- ★★`1 < C ≤ g d` なら予算 `F` は有界にできない。
⇒ ★**②（一様定数の no-go）の内容は損失 `c` ではなく伸び `g` についての主張である。** -/
theorem not_bddAbove_budget {g F : ℕ → ℝ} {C B : ℝ} (hC : 1 < C) (hF1 : ∀ d, 1 ≤ F d)
    (hg : ∀ d, C ≤ g d) (hFg : ∀ d d', d' < d → g d * F d' ≤ F d) :
    ¬ (∀ d : ℕ, F d ≤ B) := by
  intro hB
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt B hC
  have := pow_le_budget (le_of_lt hC) hF1 hg hFg n
  linarith [hB n]

/-- ★★★**損失は一様定数でよい** —— `g d ≤ 1` なら `F ≡ C` が
`exists_mem_of_descent_budget` の 3 条件を全部満たす。

⇒ ★`PGroupDescentToAxWild.const_descent_no_go`（②）を
「損失 `c` が一様定数では届かない」と読むのは**誤り**である。
★あれは `g ≡ c` と置いた形（`exists_mem_of_descent_prod`）にしか当たらない。 -/
theorem budget_const_of_growth_le_one {c g : ℕ → ℝ} {C : ℝ} (hC : 1 ≤ C)
    (hc : ∀ d, c d ≤ C) (hg : ∀ d, g d ≤ 1) :
    (∀ d : ℕ, (1:ℝ) ≤ (fun _ : ℕ => C) d)
      ∧ (∀ d : ℕ, d ≠ 0 → c d ≤ (fun _ : ℕ => C) d)
      ∧ (∀ d d' : ℕ, d' < d → g d * (fun _ : ℕ => C) d' ≤ (fun _ : ℕ => C) d) := by
  refine ⟨fun _ => hC, fun d _ => hc d, fun d _ _ => ?_⟩
  have : g d * C ≤ 1 * C := mul_le_mul_of_nonneg_right (hg d) (le_trans zero_le_one hC)
  simpa using this

end Core

/-! ## §2 具体層 —— 予算形の `hstep` は `AxWildDescentDecay` そのもの -/

section Concrete

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- 予算形の `hstep`（`c ≡ C`, `g ≡ θ`）は `AxWildDescentDecay K C θ` そのもの。
★違いは `wildDepth K x ≠ 0` と `p ∣ deg minpoly x` の言い換えだけである。 -/
theorem axWildDescentDecay_iff_budgetStep (K : PAdicLocalField p) {C θ : ℝ} :
    AxWildDescentDecay K C θ ↔
      (∀ (ε : ℝ) (x : K.closure), 0 ≤ ε → (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
        wildDepth K x ≠ 0 →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧ ‖x - x'‖ ≤ C * ε ∧
          ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ θ * ε) := by
  constructor
  · intro h ε x hε hx hd
    refine h ε hε x hx ?_
    by_contra hdvd
    exact hd ((wildDepth_eq_zero_iff K x).mpr hdvd)
  · intro h ε hε x hx hdvd
    exact h ε x hε hx (fun h0 => (wildDepth_eq_zero_iff K x).mp h0 hdvd)

/-- ★★★★**本ファイルの主結果** —— `1 ≤ C`・`0 ≤ θ ≤ 1` のとき
**予算形の `hstep` ⟺ `AxLemma K C`**。

⇒ ★★「`c` と `g` を分けて `F` で束ねる」載せ替えは、`g ≤ 1` を狙う限り
前波で塞いだ循環（`AxDecayKnobCircular.axWildDescentDecay_iff_axLemma`）に戻る。 -/
theorem budgetStep_iff_axLemma (K : PAdicLocalField p) {C : ℝ} (hC : 1 ≤ C) {θ : ℝ}
    (h0 : 0 ≤ θ) (h1 : θ ≤ 1) :
    (∀ (ε : ℝ) (x : K.closure), 0 ≤ ε → (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
        wildDepth K x ≠ 0 →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧ ‖x - x'‖ ≤ C * ε ∧
          ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ θ * ε)
      ↔ AxLemma K C :=
  (axWildDescentDecay_iff_budgetStep K).symm.trans
    (AxKnob.axWildDescentDecay_iff_axLemma K hC h0 h1)

/-- ★超距離からの評価は `g ≤ max(c, 1)` しか与えない
（`WildDescentDistanceOnly.lean:161 norm_smul_sub_self_le_max`）。
⇒ ★`g d < 1` を得るには `c d < 1`、すなわち `ε` より真に良い 1 段の近似が要る。 -/
theorem growth_le_max_of_dist (K : PAdicLocalField p) {ε D : ℝ} {x x' : K.closure}
    (hx : ∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) (hd : ‖x - x'‖ ≤ D) (σ : K.absGal) :
    ‖σ • x' - x'‖ ≤ max D ε :=
  norm_smul_sub_self_le_max (G := K.absGal) (M := K.closure)
    (fun g z => norm_smul_closure K g z) hx hd σ

end Concrete

/-! ## §3 ★★★極端な場合（`g ≡ 0`）でも同値 -/

section Summary

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★極端な場合 —— `g ≡ 0`（`ε` が 1 段で **0 に落ちる**）でも
予算形は `AxLemma K C` と**同値**である。★載せ替えの余地はここには無い。 -/
theorem budgetStep_zero_iff_axLemma (K : PAdicLocalField p) {C : ℝ} (hC : 1 ≤ C) :
    (∀ (ε : ℝ) (x : K.closure), 0 ≤ ε → (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) →
        wildDepth K x ≠ 0 →
        ∃ x' : K.closure, wildDepth K x' < wildDepth K x ∧ ‖x - x'‖ ≤ C * ε ∧
          ∀ σ : K.absGal, ‖σ • x' - x'‖ ≤ 0 * ε)
      ↔ AxLemma K C :=
  budgetStep_iff_axLemma K hC le_rfl zero_le_one

end Summary

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def pow_le_budget.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def budget_const_of_growth_le_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def budgetStep_iff_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def growth_le_max_of_dist.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms pow_le_budget
#print axioms not_bddAbove_budget
#print axioms budget_const_of_growth_le_one
#print axioms axWildDescentDecay_iff_budgetStep
#print axioms budgetStep_iff_axLemma
#print axioms growth_le_max_of_dist
#print axioms budgetStep_zero_iff_axLemma

end BudgetRecast

end ABC3.Found.PGC
