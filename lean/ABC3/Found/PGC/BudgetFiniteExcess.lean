import ABC3.Found.PGC.BudgetRecast
import ABC3.Found.PGC.FirstJumpRouteEquiv

/-!
# [pGC] ★★★④ は**記法の問題だった** —— 予算形（積の形）では破れが消える

## ★配られた問いへの答え

**★★★④（`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の破れ）は予算形で消える。**（§3、`cyclotomic_break_vanishes_in_prod`）

`WildDescentDistanceOnly.lean:282 axDecay_three_two_lt_cyclotomic_layer_loss` が測った破れは
「1 段の損失が `3^{4/18}` で、`axDecay 3 2 = 3^{3/18}` を `3^{1/18}` だけ超える」である。
★★しかし `axLemma_of_wildDescent_Icc`（`AxEpsilonDecay.lean:440`）が要求するのは
★**各段の点ごとの比較ではなく、積が `C` 以下であること**だけである。

本ファイルは実際に列を作って測った:

```
cEx 1 = axDecay 3 1 = 3^{9/18},  cEx 2 = 3^{4/18},  cEx d = 1 (d ≥ 3)
```

* ★`axDecay 3 2 < cEx 2`（点ごとには**破れている**。`cEx_two_gt_axDecay`）
* ★★`∀ n, ∏_{j ∈ Icc 1 n} cEx j = 3^{13/18} ≤ 3^{13.5/18} = axConstant 3`
  （`cEx_prod_le`。`axConstant 3 = 3^{3/4}`）

⇒ ★★★**`AxWildDescent K cEx` から `AxLemma K (axConstant 3)` と `AxSenTate K` が出る**
（§4 `axLemma_of_wildDescent_cEx` / `axSenTate_of_wildDescent_cEx`）。
★すなわち★**④ は「段ごとに一様な予算 `axDecay p k`」という記法の問題であって、
本物の穴ではない。** 木が `WildDescentDistanceOnly.lean:94-97` で
「★★**合成(`AxLemma`)は無傷である**」と地の文で書いていたことを、
★本ファイルが**定理にした**（`13/18 ≤ 3/4` の 1 行に落ちる）。

**★`g d > 1` の段が有限個の形（`sup F`）を形式化した。**（§1 `budget_of_eventually_le_one`）

`exists_mem_of_descent_budget`（`AxTowerDecay.lean:327`）の 3 条件は、
`g d ≤ G`（全部）・`g d ≤ 1`（`d > N`）・`c d ≤ C` から

  `F d := C · G^{min d N}`

で**全部満たされる**。★`hFg` は `d' < d` の**すべての対**を要求するが、`min` を使えば通る
（★持ち場が「`hFg` が何を要求するかは読んでいない」と書いた点の答え）。
★`F` は `C·G^N` で**有界**（`budget_le_of_eventually_le_one`）。
⇒ ★★**悪い段が有限個なら予算は有界**である。

★★ただし前波で測ったとおり、★**有界な `F` の予算形は `AxLemma K (sup F)` と同値**
（`BudgetRecast.budgetStep_iff_axLemma`）なので、★これは道を短くしない。
★**効いているのは §3 の方**（点ごとの比較を積の比較に緩める）である。

## ★`c d < 1` は③と同じ不等式か —— ★別の量である（★形式化していない。測定として記録）

前波で「`g d < 1` には `c d < 1` が要る」と測った（`BudgetRecast.growth_le_max_of_dist`）。
③（`JumpGeometricDecay.no_uniform_geometric_of_harith`）は
`‖σ^{p^k}π − π‖ / ‖σπ − π‖ ≥ ‖p‖^{1/(p−1)}` という★**同じ `π` の、`σ` の冪による変位の比**
についての下界である。一方 `c d` は★**`x` と部分体との距離**を `ε` で割ったものである。
⇒ ★**別の量なので、③はそのままでは `c d < 1` を否定しない。**
★（`c d < 1` を否定する測定は本波では行っていない。）

## ★4 つの穴の現状（本波の再測定）

| 穴 | 現状 |
|---|---|
| ①不分岐（`TotallyRamifiedLayer.lean:342`） | ★生きている |
| ②一様定数（`PGroupDescentToAxWild.const_descent_no_go`） | ★射程は `g`（前波で訂正） |
| ③幾何減衰（`JumpGeometricDecay.no_uniform_geometric_of_harith`） | ★生きている（変位の比の話） |
| ④測定点（`FirstJumpRouteEquiv.not_exists_firstJump_data_cyclotomic`） | ★★★**消えた**（本波、§3）。記法の問題だった |

★★**穴は 4 つから 3 つに減った。**

## ★まだ測っていないこと（正直に）

* `WildDescentDistanceOnly.lean:98-99` の「`K = ℚ₃(ζ₃)` の 40 個の巡回 3 次拡大で
  `d(x,N)` を実際に測る必要がある。★測っていない」は★**今も未測定**である
  （本波でも測っていない）。★ただし §3 により、そこを測らなくても
  ★**④ は `AxLemma` の道を塞がない**。
* 本ファイルの `cEx` は★**実在する降下から作った列ではない**。
  「点ごとの比較が破れても積は通る」ことを示す**証拠列**である。
  ★実際の降下がこの `cEx` を実現するかは別問題で、測っていない。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ項目（pGC 物理 p.6 Corollary 3.1）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1・§2 は `ℝ` と `ℕ` と `Finset` だけで、分岐・付値・Galois の語彙が 1 語も出ない。
4. §4 は `Fact (Nat.Prime 3)` を `variable` で受けている（`p = 3` の具体例なので
   global instance にしていない）。
-/

namespace ABC3.Found.PGC

namespace BudgetFinite

/-! ## §1 抽象核 —— 悪い段が有限個なら予算 `F` は有界に取れる -/

section Core

/-- ★★**悪い段が有限個なら予算は有界に取れる** —— `exists_mem_of_descent_budget`
（`AxTowerDecay.lean:327`）の 3 条件を `F d := C · G^{min d N}` が全部満たす。

`hFg` は `d' < d` の**すべての対**を要求するが、`min d N` を使えば通る。
★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。 -/
theorem budget_of_eventually_le_one {c g : ℕ → ℝ} {C G : ℝ} {N : ℕ}
    (hC : 1 ≤ C) (hG : 1 ≤ G)
    (hc : ∀ d, c d ≤ C) (hgG : ∀ d, g d ≤ G) (hg1 : ∀ d, N < d → g d ≤ 1) :
    (∀ d : ℕ, (1:ℝ) ≤ C * G ^ (min d N))
      ∧ (∀ d : ℕ, d ≠ 0 → c d ≤ C * G ^ (min d N))
      ∧ (∀ d d' : ℕ, d' < d → g d * (C * G ^ (min d' N)) ≤ C * G ^ (min d N)) := by
  have hC0 : (0:ℝ) < C := lt_of_lt_of_le zero_lt_one hC
  have hG0 : (0:ℝ) < G := lt_of_lt_of_le zero_lt_one hG
  have hpow1 : ∀ n : ℕ, (1:ℝ) ≤ G ^ n := fun n => one_le_pow₀ hG
  have hbase : ∀ d : ℕ, C ≤ C * G ^ (min d N) := by
    intro d
    nlinarith [hpow1 (min d N)]
  refine ⟨fun d => le_trans hC (hbase d), fun d _ => le_trans (hc d) (hbase d), ?_⟩
  intro d d' hlt
  rcases le_or_gt d N with hdN | hdN
  · have hmin : min d N = d := min_eq_left hdN
    have hmin' : min d' N = d' := min_eq_left (le_of_lt (lt_of_lt_of_le hlt hdN))
    rw [hmin, hmin']
    have hstep : G * G ^ d' ≤ G ^ d := by
      have : G ^ (d' + 1) ≤ G ^ d := pow_le_pow_right₀ hG (by omega)
      calc G * G ^ d' = G ^ (d' + 1) := by ring
        _ ≤ G ^ d := this
    calc g d * (C * G ^ d') ≤ G * (C * G ^ d') :=
          mul_le_mul_of_nonneg_right (hgG d) (by positivity)
      _ = C * (G * G ^ d') := by ring
      _ ≤ C * G ^ d := mul_le_mul_of_nonneg_left hstep (le_of_lt hC0)
  · have hmin : min d N = N := min_eq_right (le_of_lt hdN)
    have hle : G ^ (min d' N) ≤ G ^ N := pow_le_pow_right₀ hG (min_le_right _ _)
    rw [hmin]
    calc g d * (C * G ^ (min d' N)) ≤ 1 * (C * G ^ (min d' N)) :=
          mul_le_mul_of_nonneg_right (hg1 d hdN) (by positivity)
      _ = C * G ^ (min d' N) := by ring
      _ ≤ C * G ^ N := mul_le_mul_of_nonneg_left hle (le_of_lt hC0)

/-- ★その `F` は `C · G^N` で有界。 -/
theorem budget_le_of_eventually_le_one {C G : ℝ} {N : ℕ} (hC : 1 ≤ C) (hG : 1 ≤ G) (d : ℕ) :
    C * G ^ (min d N) ≤ C * G ^ N := by
  have hC0 : (0:ℝ) < C := lt_of_lt_of_le zero_lt_one hC
  exact mul_le_mul_of_nonneg_left (pow_le_pow_right₀ hG (min_le_right _ _)) (le_of_lt hC0)

end Core

/-! ## §2 抽象核 —— 途中から `1` になる列の有限積は最大でも `Icc 1 N` の積 -/

section ProdCore

/-- ★抽象核: `N` より先で `c d = 1` かつ `1 ≤ c d` なら、
`Icc 1 n` 上の積は `Icc 1 N` 上の積を超えない（`n` に依らない上界）。 -/
theorem prod_Icc_le_of_eventually_one {c : ℕ → ℝ} {N : ℕ}
    (hc1 : ∀ d, 1 ≤ c d) (heq : ∀ d, N < d → c d = 1) (n : ℕ) :
    ∏ j ∈ Finset.Icc 1 n, c j ≤ ∏ j ∈ Finset.Icc 1 N, c j := by
  rcases le_or_gt n N with h | h
  · exact Finset.prod_le_prod_of_subset_of_one_le
      (Finset.Icc_subset_Icc le_rfl h)
      (fun i _ => le_trans zero_le_one (hc1 i))
      (fun i _ _ => hc1 i)
  · refine le_of_eq (Finset.prod_subset (Finset.Icc_subset_Icc le_rfl (le_of_lt h)) ?_).symm
    intro x hx hnx
    simp only [Finset.mem_Icc] at hx
    simp only [Finset.mem_Icc, not_and, not_le] at hnx
    exact heq x (hnx hx.1)

end ProdCore

/-! ## §3 ★★★④ は積の形では消える（`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の測定点） -/

section Cyclotomic

open scoped BigOperators

/-- ★★測定点 `ℚ₃(ζ₂₇)/ℚ₃(ζ₃)` の破れを含む**証拠列**。

`cEx 1 = axDecay 3 1 = 3^{9/18}`、`cEx 2 = 3^{4/18}`（★破れ）、`cEx d = 1`（`d ≥ 3`）。
★実在する降下から作った列では**ない**（点ごとの比較が破れても積は通ることの証拠）。 -/
noncomputable def cEx : ℕ → ℝ :=
  fun d => if d = 1 then axDecay 3 1
    else if d = 2 then (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) else 1

theorem one_le_cEx (d : ℕ) : 1 ≤ cEx d := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  unfold cEx
  split_ifs with h1 h2
  · exact one_le_axDecay 3 1
  · exact Real.one_le_rpow (by norm_num) (by positivity)
  · exact le_rfl

theorem cEx_eq_one_of_gt (d : ℕ) (hd : 2 < d) : cEx d = 1 := by
  unfold cEx
  split_ifs with h1 h2
  · omega
  · omega
  · rfl

theorem cEx_prod_Icc_two :
    ∏ j ∈ Finset.Icc 1 2, cEx j
      = axDecay 3 1 * (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) := by
  have hset : Finset.Icc (1:ℕ) 2 = {1, 2} := by decide
  rw [hset, Finset.prod_insert (by decide), Finset.prod_singleton]
  norm_num [cEx]

/-- ★点ごとには破れている（`WildDescentDistanceOnly.lean:282` の測定そのもの）。 -/
theorem cEx_two_gt_axDecay : axDecay 3 2 < cEx 2 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h : cEx 2 = (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) := by norm_num [cEx]
  rw [h]
  exact axDecay_three_two_lt_cyclotomic_layer_loss

/-- ★したがって `axLemma_of_axDecay` の `hdecay`（点ごとの比較）は満たさない。 -/
theorem cEx_not_pointwise_le_axDecay : ¬ (∀ d : ℕ, cEx d ≤ axDecay 3 d) := by
  intro h
  exact absurd (h 2) (not_le.mpr cEx_two_gt_axDecay)

/-- ★★**合成は通る** —— `3^{9/18} · 3^{4/18} = 3^{13/18} ≤ 3^{13.5/18} = axConstant 3`。
★木が `WildDescentDistanceOnly.lean:94-97` で地の文で書いた「合成は無傷」の中身。 -/
theorem composite_bound :
    axDecay 3 1 * (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) ≤ axConstant 3 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  have h3 : (0:ℝ) < 3 := by norm_num
  rw [JumpArith.axDecay_three_one_eq, ← Real.rpow_add h3, axConstant]
  refine Real.rpow_le_rpow_of_exponent_le (by norm_num) ?_
  norm_num

/-- ★★★`axLemma_of_wildDescent_Icc`（`AxEpsilonDecay.lean:440`）の積の仮説は**満たされる**。 -/
theorem cEx_prod_le (n : ℕ) : ∏ j ∈ Finset.Icc 1 n, cEx j ≤ axConstant 3 := by
  refine (prod_Icc_le_of_eventually_one one_le_cEx cEx_eq_one_of_gt n).trans ?_
  rw [cEx_prod_Icc_two]
  exact composite_bound

end Cyclotomic

/-! ## §4 その列から `AxLemma` / `AxSenTate` が出る -/

section Route

open ABC3.Skeleton.PGC

variable [Fact (Nat.Prime 3)]

/-- ★★★★**④ は本物の穴ではない** —— 点ごとの比較が破れている `cEx` でも
`AxWildDescent K cEx` から `AxLemma K (axConstant 3)` が出る。 -/
theorem axLemma_of_wildDescent_cEx (K : PAdicLocalField 3) (h : AxWildDescent K cEx) :
    AxLemma K (axConstant 3) :=
  axLemma_of_wildDescent_Icc K one_le_cEx cEx_prod_le h

theorem axSenTate_of_wildDescent_cEx (K : PAdicLocalField 3) (h : AxWildDescent K cEx) :
    AxSenTate K :=
  axSenTate_of_axLemma (le_trans zero_le_one (one_le_axConstant 3))
    (axLemma_of_wildDescent_cEx K h)

end Route

/-! ## §5 まとめの 1 本 -/

section Summary

/-- ★★★★**本ファイルの主結果を 1 本にした形** ——
点ごとには `axDecay 3 2 < cEx 2`（破れ）だが、積は `axConstant 3` に収まる。
⇒ ★④ は「段ごとに一様な予算」という**記法の問題**であった。 -/
theorem cyclotomic_break_vanishes_in_prod :
    axDecay 3 2 < cEx 2 ∧ (∀ n : ℕ, ∏ j ∈ Finset.Icc 1 n, cEx j ≤ axConstant 3) :=
  ⟨cEx_two_gt_axDecay, cEx_prod_le⟩

end Summary

/-! ## §6 `.src` と 使っている公理の一覧 -/

section Src

def budget_of_eventually_le_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def prod_Icc_le_of_eventually_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def cyclotomic_break_vanishes_in_prod.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axSenTate_of_wildDescent_cEx.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms budget_of_eventually_le_one
#print axioms budget_le_of_eventually_le_one
#print axioms prod_Icc_le_of_eventually_one
#print axioms one_le_cEx
#print axioms cEx_prod_Icc_two
#print axioms cEx_two_gt_axDecay
#print axioms cEx_not_pointwise_le_axDecay
#print axioms composite_bound
#print axioms cEx_prod_le
#print axioms axLemma_of_wildDescent_cEx
#print axioms axSenTate_of_wildDescent_cEx
#print axioms cyclotomic_break_vanishes_in_prod

end BudgetFinite

end ABC3.Found.PGC
