import ABC3.Found.PGC.LossToExit

/-!
# [pGC] 基底段を別扱いする —— 有限個の例外は `B^{K₀}` 倍で吸収できる

## 持ち場（前波で「次の 1 点」とした点）

`(p−1)·(t k + p − 2) ≤ e k` が**浅い層で成り立たない**（`p=3, k=3` で `4 ≤ 3`）。

## ★選んだ道と、その理由

3 つ考えられた:

1. 基底段だけ別の定数で処理する ←★**これを選んだ**
2. 有限個を落として総和を取り直す
3. `+(p−2)` の余白を別の形で吸収する

★1 を選んだのは、出口 `axLemma_of_wildDescent` の仮定が
`hC : ∀ s : Finset ℕ, ∏ k ∈ s, c k ≤ C` という ★**「任意の有限積」の形**だからである。
有限個の例外は積の中で有限個の因子にしかならないので、★**定数を `B^{K₀}` 倍する**だけで済む。
★2（総和を取り直す）は `axDecay` の指数の形を変えることになり、木の
`prod_le_rpow_of_geometric` が使えなくなる。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `prod_le_of_eventually` | ★抽象核。分岐も付値も出てこない（`Finset` と `ℝ` だけ） |
| §2 | `axLemma_of_eventually` | ★★出口に差し込んだ形（`AxLemma K (B^{K₀}·C)`） |

## ★これで何が変わったか

前波の到達点 `LossToExit.axLemma_of_loss_bound` は
`hjump : ∀ k, (p−1)(t k + p − 2) ≤ e k` を **すべての `k`** に要求していた。
★本ファイルにより、★**浅い `k` は `c k ≤ B` で押さえるだけでよい**ことになった。

★`B` として何が取れるかは**測っていない**。★`AxWildDescent` の `c k` は
「1 段の損失の上界」なので、浅い層でも有限の値ではあるが、★**その具体的な値は
本日の鎖では出していない**（`ε` が有限倍しか増えないことは `axWildDescent_pow` が
`c k = p^k` で無条件に真だと言っているので、少なくとも `B = p^{K₀}` 程度は取れるはずだが、
★**確かめていない**）。

## 逸脱の記録

- §1 の `hC0 : 0 ≤ C` は**外した** —— 書いたが使わなかった（`mul_le_mul` の非負条件は
  `0 ≤ B^{K₀}` の側だけで足りる）。★「仮定は書く前に外してみる」の 4 例目。
- ★`d` の積の仮定は「任意の `Finset`」の形で受けている。木の `axLemma_of_axDecay` は
  `Icc 1 n` の形（`axLemma_of_wildDescent_Icc`）なので、★そちらと繋ぐには
  部分集合への単調性（因子が `1` 以上）が 1 本要る。★本ファイルには**無い**。
-/

namespace ABC3.Found.PGC

namespace FiniteExceptions

open Finset

/-! ## §1 抽象核 —— 有限個の例外は定数倍で吸収できる -/

section Kernel

/-- ★★抽象核（分岐も付値も出てこない）。

`k < K₀` では `c k ≤ B`、`K₀ ≤ k` では `c k ≤ d k` で、`d` の任意有限積が `C` 以下なら、
`c` の任意有限積は `B^{K₀}·C` 以下。

★これが「基底段だけ別扱い」の中身である。★定数は `B^{K₀}` 倍に悪くなるが**有限**。 -/
theorem prod_le_of_eventually {c d : ℕ → ℝ} {K₀ : ℕ} {B C : ℝ}
    (hc1 : ∀ k, 1 ≤ c k) (hB : 1 ≤ B)
    (hsmall : ∀ k, k < K₀ → c k ≤ B)
    (hbig : ∀ k, K₀ ≤ k → c k ≤ d k)
    (hD : ∀ s : Finset ℕ, ∏ k ∈ s, d k ≤ C)
    (s : Finset ℕ) : ∏ k ∈ s, c k ≤ B ^ K₀ * C := by
  classical
  have hc0 : ∀ k, (0 : ℝ) ≤ c k := fun k => le_trans zero_le_one (hc1 k)
  have hsplit : (∏ k ∈ s.filter (fun k => k < K₀), c k)
      * ∏ k ∈ s.filter (fun k => ¬ k < K₀), c k = ∏ k ∈ s, c k :=
    Finset.prod_filter_mul_prod_filter_not s _ c
  have h1 : ∏ k ∈ s.filter (fun k => k < K₀), c k ≤ B ^ K₀ := by
    have hle : ∏ k ∈ s.filter (fun k => k < K₀), c k
        ≤ ∏ _k ∈ s.filter (fun k => k < K₀), B := by
      refine Finset.prod_le_prod (fun i _ => hc0 i) ?_
      intro i hi
      exact hsmall i (Finset.mem_filter.mp hi).2
    have hcard : (s.filter (fun k => k < K₀)).card ≤ K₀ := by
      have hsub : s.filter (fun k => k < K₀) ⊆ range K₀ := by
        intro i hi
        exact Finset.mem_range.mpr (Finset.mem_filter.mp hi).2
      simpa using Finset.card_le_card hsub
    calc ∏ k ∈ s.filter (fun k => k < K₀), c k
        ≤ ∏ _k ∈ s.filter (fun k => k < K₀), B := hle
      _ = B ^ (s.filter (fun k => k < K₀)).card := by rw [Finset.prod_const]
      _ ≤ B ^ K₀ := pow_le_pow_right₀ hB hcard
  have h2 : ∏ k ∈ s.filter (fun k => ¬ k < K₀), c k ≤ C := by
    refine le_trans ?_ (hD (s.filter (fun k => ¬ k < K₀)))
    refine Finset.prod_le_prod (fun i _ => hc0 i) ?_
    intro i hi
    exact hbig i (not_lt.mp (Finset.mem_filter.mp hi).2)
  have hp1 : (0 : ℝ) ≤ ∏ k ∈ s.filter (fun k => k < K₀), c k :=
    Finset.prod_nonneg fun i _ => hc0 i
  have hp2 : (0 : ℝ) ≤ ∏ k ∈ s.filter (fun k => ¬ k < K₀), c k :=
    Finset.prod_nonneg fun i _ => hc0 i
  rw [← hsplit]
  exact mul_le_mul h1 h2 hp2 (pow_nonneg (le_trans zero_le_one hB) K₀)

end Kernel

/-! ## §2 ★出口に差し込む —— 基底段を別扱いしても `AxLemma` は出る -/

section Exit

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**基底段を別扱いした出口**。

`k < K₀` の有限個だけ `c k ≤ B` で押さえ、`K₀ ≤ k` で `c k ≤ d k`（`d` の積は `C` 以下）なら
`AxLemma K (B^{K₀}·C)`。★定数は悪くなるが**有限**で、`axLemma_of_wildDescent` が要求する
「任意の有限積が押さえられる」を満たす。 -/
theorem axLemma_of_eventually (K : PAdicLocalField p) {c d : ℕ → ℝ} {K₀ : ℕ} {B C : ℝ}
    (hc : ∀ k, 1 ≤ c k) (hB : 1 ≤ B)
    (hsmall : ∀ k, k < K₀ → c k ≤ B)
    (hbig : ∀ k, K₀ ≤ k → c k ≤ d k)
    (hD : ∀ s : Finset ℕ, ∏ k ∈ s, d k ≤ C)
    (h : AxWildDescent K c) : AxLemma K (B ^ K₀ * C) :=
  axLemma_of_wildDescent K hc (prod_le_of_eventually hc hB hsmall hbig hD) h

end Exit

/-! ## §3 使っている公理の一覧 -/

#print axioms prod_le_of_eventually
#print axioms axLemma_of_eventually

end FiniteExceptions

end ABC3.Found.PGC
