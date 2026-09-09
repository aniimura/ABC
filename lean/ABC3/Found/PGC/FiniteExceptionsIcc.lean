import ABC3.Found.PGC.FiniteExceptions

/-!
# [pGC] ★★前波の形は `axDecay` に**使えなかった** —— `Icc 1 n` 版に直す

## ★★自己訂正（13 度目）

前波の `FiniteExceptions.axLemma_of_eventually` は
`hD : ∀ s : Finset ℕ, ∏ k ∈ s, d k ≤ C`（**任意の `Finset`**）を要求している。
★**この形は `d = axDecay p` では満たされない。**

理由: `axDecay p k` の指数は `(1/(p−1))·(1/p)^{k−1}` で、★`ℕ` の引き算なので `0 − 1 = 0`。
⇒ ★**`axDecay p 0 = axDecay p 1`**（§3 で `rw [axDecay, axDecay]` だけで証明した）。
`{0,1}` を含む集合では同じ因子が 2 回入り、指数の総和が `axConstant` の分を超える。

★★**木も同じことを言っていた**（`AxEpsilonDecay.lean:439`、本波で読んだ）:

> ★`AxTowerDecay.axLemma_of_wildDescent` の、積の条件を **`Finset.Icc 1 n` だけ**に
> 弱めた版。★**ずらした減衰では `c 0` に条件が付かないので、この形が要る**。

★前波の私は「`Icc` と繋ぐには単調性が 1 本要る」と書いたが、★**「あると便利」ではなく
「そうでないと使えない」**だった。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `prod_Icc_le_of_eventually` | ★`Icc 1 n` 版の抽象核（部分集合への単調性を内側で処理） |
| §2 | `axLemma_of_eventually_Icc` | ★★出口（`axLemma_of_wildDescent_Icc` に接続） |
| §3 | `axDecay_zero_eq_one` | ★訂正の証拠 |

★単調性は外の補題を探さず、★**`prod_filter_mul_prod_filter_not` ＋ `le_mul_of_one_le_left`**
で内側に埋めた（`ℝ` は乗法の順序モノイドではないので
`Finset.prod_le_prod_of_subset_of_one_le'` がそのままでは当たらない）。

## ★在庫の測定（本体が渡した行番号を使った）

- `axLemma_of_wildDescent_Icc` は ★**`AxEpsilonDecay.lean:440`** に在り、名前空間は
  `ABC3.Found.PGC` **のみ**（`grep -n '^namespace' lean/ABC3/Found/PGC/AxEpsilonDecay.lean`
  → `:197` の 1 行だけ）。★#359 の確認を先にした結果、修飾なしで引けた。

## 逸脱の記録

- ★前波の `FiniteExceptions.axLemma_of_eventually`（任意 `Finset` 版）は**削除していない**
  （規約: 既存ファイルは触らない）。★`d` が `k = 0` でも条件を持つ場合には使える。
  ★`axDecay` には**使えない**ことを本ファイルに記録した。
- §3 は `axDecay` の定義を 2 回展開するだけで閉じる。★`ℕ` の引き算が原因だと分かる形にした。
-/

namespace ABC3.Found.PGC

namespace FiniteExceptionsIcc

open Finset

/-! ## §1 `Icc 1 n` 版の抽象核 -/

section Kernel

/-- ★★前波の `FiniteExceptions.prod_le_of_eventually` の **`Icc 1 n` 版**。

★前波の形（任意の `Finset`）は ★**`d = axDecay p` では使えない** ——
`axDecay p 0 = axDecay p 1`（`ℕ` の引き算で `0 − 1 = 0`）なので `{0,1}` を含む集合では
指数が `axConstant` の分を超える。★木も同じことを言っている
（`AxEpsilonDecay.lean:439`「ずらした減衰では `c 0` に条件が付かないので、この形が要る」）。 -/
theorem prod_Icc_le_of_eventually {c d : ℕ → ℝ} {K₀ : ℕ} {B C : ℝ}
    (hc1 : ∀ k, 1 ≤ c k) (hd1 : ∀ k, 1 ≤ d k) (hB : 1 ≤ B)
    (hsmall : ∀ k, k < K₀ → c k ≤ B)
    (hbig : ∀ k, K₀ ≤ k → c k ≤ d k)
    (hD : ∀ n : ℕ, ∏ k ∈ Finset.Icc 1 n, d k ≤ C)
    (n : ℕ) : ∏ k ∈ Finset.Icc 1 n, c k ≤ B ^ K₀ * C := by
  classical
  have hc0 : ∀ k, (0 : ℝ) ≤ c k := fun k => le_trans zero_le_one (hc1 k)
  have hd0 : ∀ k, (0 : ℝ) ≤ d k := fun k => le_trans zero_le_one (hd1 k)
  have hsplit : (∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), c k)
      * ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), c k
      = ∏ k ∈ Finset.Icc 1 n, c k :=
    Finset.prod_filter_mul_prod_filter_not _ _ c
  have h1 : ∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), c k ≤ B ^ K₀ := by
    have hle : ∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), c k
        ≤ ∏ _k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), B :=
      Finset.prod_le_prod (fun i _ => hc0 i)
        (fun i hi => hsmall i (Finset.mem_filter.mp hi).2)
    have hcard : ((Finset.Icc 1 n).filter (fun k => k < K₀)).card ≤ K₀ := by
      have hsub : (Finset.Icc 1 n).filter (fun k => k < K₀) ⊆ range K₀ := by
        intro i hi
        exact Finset.mem_range.mpr (Finset.mem_filter.mp hi).2
      simpa using Finset.card_le_card hsub
    calc ∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), c k
        ≤ ∏ _k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), B := hle
      _ = B ^ ((Finset.Icc 1 n).filter (fun k => k < K₀)).card := by rw [Finset.prod_const]
      _ ≤ B ^ K₀ := pow_le_pow_right₀ hB hcard
  have h2 : ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), c k ≤ C := by
    have hstep : ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), c k
        ≤ ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), d k :=
      Finset.prod_le_prod (fun i _ => hc0 i)
        (fun i hi => hbig i (not_lt.mp (Finset.mem_filter.mp hi).2))
    have hdsplit : (∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), d k)
        * ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), d k
        = ∏ k ∈ Finset.Icc 1 n, d k :=
      Finset.prod_filter_mul_prod_filter_not _ _ d
    have hone : (1 : ℝ) ≤ ∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), d k :=
      Finset.one_le_prod (fun i _ => hd1 i)
    have hmono : ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), d k
        ≤ ∏ k ∈ Finset.Icc 1 n, d k := by
      calc ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), d k
          ≤ (∏ k ∈ (Finset.Icc 1 n).filter (fun k => k < K₀), d k)
            * ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), d k :=
            le_mul_of_one_le_left (Finset.prod_nonneg fun i _ => hd0 i) hone
        _ = ∏ k ∈ Finset.Icc 1 n, d k := hdsplit
    exact le_trans hstep (le_trans hmono (hD n))
  have hp2 : (0 : ℝ) ≤ ∏ k ∈ (Finset.Icc 1 n).filter (fun k => ¬ k < K₀), c k :=
    Finset.prod_nonneg fun i _ => hc0 i
  rw [← hsplit]
  exact mul_le_mul h1 h2 hp2 (pow_nonneg (le_trans zero_le_one hB) K₀)

end Kernel

/-! ## §2 ★出口（`Icc` 版） -/

section Exit

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**基底段を別扱いした出口（`Icc` 版）**。★こちらが `axDecay` に使える形である。 -/
theorem axLemma_of_eventually_Icc (K : PAdicLocalField p) {c d : ℕ → ℝ} {K₀ : ℕ} {B C : ℝ}
    (hc : ∀ k, 1 ≤ c k) (hd1 : ∀ k, 1 ≤ d k) (hB : 1 ≤ B)
    (hsmall : ∀ k, k < K₀ → c k ≤ B)
    (hbig : ∀ k, K₀ ≤ k → c k ≤ d k)
    (hD : ∀ n : ℕ, ∏ k ∈ Finset.Icc 1 n, d k ≤ C)
    (h : AxWildDescent K c) : AxLemma K (B ^ K₀ * C) :=
  axLemma_of_wildDescent_Icc K hc (prod_Icc_le_of_eventually hc hd1 hB hsmall hbig hD) h

end Exit

/-! ## §3 ★訂正の証拠 —— `axDecay p 0 = axDecay p 1` -/

section Evidence

variable {p : ℕ} [Fact p.Prime]

omit [Fact p.Prime] in
/-- ★★前波の形（任意の `Finset`）が `axDecay` に使えない理由。

`axDecay p k` の指数は `(1/(p−1))·(1/p)^{k−1}` で、★`ℕ` の引き算なので `0 − 1 = 0`。
⇒ `axDecay p 0 = axDecay p 1` である。★`{0,1}` を含む集合では同じ因子が 2 回入り、
指数の総和が `axConstant` の分（`= A/(1−r)`）を超える。 -/
theorem axDecay_zero_eq_one : axDecay p 0 = axDecay p 1 := by
  rw [axDecay, axDecay]

end Evidence

/-! ## §4 使っている公理の一覧 -/

#print axioms prod_Icc_le_of_eventually
#print axioms axLemma_of_eventually_Icc
#print axioms axDecay_zero_eq_one

end FiniteExceptionsIcc

end ABC3.Found.PGC
