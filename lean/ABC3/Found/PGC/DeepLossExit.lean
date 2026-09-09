import ABC3.Found.PGC.BaseLayerConstant

/-!
# [pGC] 損失評価の 4 仮説を**深い段だけ**にした出口 ＋ `K₀ = 4` が**最適**であることの測定

## 持ち場（前波で「残りはちょうど 1 点」とした点）

前波の私の言葉（逐語）:

> **残りはちょうど 1 点**: `hbig : ∀ k, K₀ ≤ k → cdeep k ≤ axDecay p k` を満たす `cdeep` と
> `AxWildDescent K cdeep` を、**深さ `K₀` 以上でだけ**作ること。

## ★測って分かったこと 1 —— `K₀ = 4` は `p ≥ 3` の**すべて**で最適

`LossToExit.lean:130` の私は `p = 3` の 2 点（`k = 3` で破れ、`k = 4` で入る）しか測っていなかった。
★本波で `p` について閉じた形にした（§1、純 `ℕ` の算術）:

* `4 ≤ k` なら `2(p−1) ≤ p^{k−2}`（★`p ≥ 3` のすべてで成立）
* `k ≤ 3` なら `2(p−1) ≤ p^{k−2}` は ★**`p ≥ 3` のすべてで偽**

⇒ ★★**`K₀ = 4` は「大きめに取った安全側」ではなく、`p ≥ 3` について一様に最適**である。
★前波の `fit_fails_at_three` / `fit_holds_at_four` は `p = 3` の 1 例にすぎなかった。
★私はそこから `K₀` が `p` とともに動くかどうかを**測っていなかった**。動かない。

## ★測って分かったこと 2 —— `p = 2` は「余白 0、でも `K₀ = 3`」

★★**本ファイルの下書きで私はここを間違えた**（書く前に測って直した）。
下書きの断定は「`p = 2` では `K₀ = 0`（基底段の別扱いが要らない）」だった。★**偽**。

正しくは 2 つの別の事柄に分かれる:

1. ★`p − 2 = 0` なので `(p−1)(t + p − 2) = (p−1)t` —— ★**余白は本当に 0** で、
   要る不等式は木の `RamificationJumpBound` が与える形**そのもの**（§2b `slack_vanishes_at_two`）。
2. ★しかし円分塔の模型では `t = 2`、`e = 2^{k−2}` なので `2 ≤ 2^{k−2}`、
   すなわち ★**`k ≥ 3` が要る**（§2b、`k = 2` で `2 ≤ 1` は偽）。

⇒ ★「余白が 0」と「基底段の別扱いが要らない」は**別のこと**だった。
測定: `tools/deep-loss-k0-check.py`（厳密整数、`p ∈ {2,3,5,7,11,13,101}`, `k ∈ [0,40)`）
→ ★`p = 2` は最小 `K₀ = 3`、`p ≥ 3` はすべて最小 `K₀ = 4`。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `two_mul_pred_le_pow` / `not_two_mul_pred_le_pow` | ★**抽象核**（純 `ℕ`。素数も体も出ない） |
| §2 | `cyclotomic_jump_fits` / `cyclotomic_jump_fails` | 円分塔の模型（`t = p`, `e = p^{k−2}(p−1)`）へ代入 |
| §2b | `slack_vanishes_at_two` / `two_fails_at_two` / `two_holds_at_three` | ★`p = 2` の訂正 |
| §3 | `deep_loss_le_axDecay` | `LossToExit.axDecay_of_loss_bound` を深い段だけに絞る |
| §4 | `axLemma_of_deep_loss_bound` / `axSenTate_of_deep_loss_bound` | ★★出口（4 仮説がすべて `K₀ ≤ k` 付きになる） |
| §5 | `axSenTate_of_cyclotomic_loss` | ★★★`K₀ = 4` で `hm` / `he` / `hjump` を**消した**形 |

## ★★出口の仮説が 5 → 2 になった

`LossToExit.axSenTate_of_loss_bound` は 5 仮説（`hc`, `hform`, `hm`, `he`, `hjump`）＋ `h`。
§5 は ★**`hform`（深さ 4 以上でだけ）と `h : AxWildDescent K c` の 2 つだけ**。

* `hc : ∀ k, 1 ≤ c k` —— ★`BaseLayerConstant` で消えた（貼り合わせ後の下界が自動）
* `hm` / `he` —— 円分塔の模型で `rfl` と `positivity`
* `hjump` —— ★§2 で証明した（`K₀ = 4`）

## 逸脱の記録

- §5 は円分塔の模型（`m k = p^{k−1}`, `e k = p^{k−2}(p−1)`, `t k = p`）を**代入した**形であり、
  ★「任意の `K` について `e k` がこの値だ」と主張しているのではない。
  ★一般の `K` に対しては §4（`hm`/`he`/`hjump` を仮定に残した形）を使う。
- §1 は `Nat` の切り捨て引き算のまま書いた（`p - 1`, `p - 2`, `k - 2`）。
  ★木の `RamificationJumpBound` と `LossToExit` が同じ書き方なので合わせた。
- import は `BaseLayerConstant` 1 本。★`LossToExit` は
  `BaseLayerConstant → FiniteExceptionsIcc → FiniteExceptions → LossToExit` で推移的に入る。
-/

namespace ABC3.Found.PGC

namespace DeepLossExit

/-! ## §1 ★抽象核 —— 純 `ℕ` の算術（素数も体も分岐も出ない） -/

section Kernel

/-- ★★`4 ≤ k` かつ `3 ≤ p` なら `2(p−1) ≤ p^{k−2}`。

★これが「`+(p−2)` の余白が入る」ことの全部である。
証明は `2(p−1) ≤ 2p ≤ p·p = p² ≤ p^{k−2}`。 -/
theorem two_mul_pred_le_pow {p k : ℕ} (hp : 3 ≤ p) (hk : 4 ≤ k) :
    2 * (p - 1) ≤ p ^ (k - 2) := by
  have h1 : 2 * (p - 1) ≤ p * p := by nlinarith [Nat.sub_le p 1, Nat.sub_lt_self one_pos (by omega : 1 ≤ p)]
  have h2 : p * p ≤ p ^ (k - 2) := by
    have : p ^ 2 ≤ p ^ (k - 2) := Nat.pow_le_pow_right (by omega) (by omega)
    simpa [pow_two] using this
  exact h1.trans h2

/-- ★★逆向き —— `k ≤ 3` なら `2(p−1) ≤ p^{k−2}` は **`p ≥ 3` のすべてで偽**。

`k ≤ 3` では `k − 2 ≤ 1` なので `p^{k−2} ≤ p`。一方 `2(p−1) = 2p − 2 > p ⟺ p > 2`。
⇒ ★**`K₀ = 4` は最適**（`K₀ = 3` では `p` をどう選んでも破れる）。 -/
theorem not_two_mul_pred_le_pow {p k : ℕ} (hp : 3 ≤ p) (hk : k ≤ 3) :
    ¬ (2 * (p - 1) ≤ p ^ (k - 2)) := by
  have hle : p ^ (k - 2) ≤ p := by
    calc p ^ (k - 2) ≤ p ^ 1 := Nat.pow_le_pow_right (by omega) (by omega)
      _ = p := pow_one p
  omega

end Kernel

/-! ## §2 円分塔の模型への代入（`t = p`, `e = p^{k−2}(p−1)`） -/

section Cyclotomic

/-- ★★木の出口が要る形 `(p−1)·(t + (p−2)) ≤ e` を、円分塔（`t = p`,
`e = p^{k−2}(p−1)`）で書いたもの。★`4 ≤ k` なら `p ≥ 3` のすべてで成り立つ。 -/
theorem cyclotomic_jump_fits {p k : ℕ} (hp : 3 ≤ p) (hk : 4 ≤ k) :
    (p - 1) * (p + (p - 2)) ≤ p ^ (k - 2) * (p - 1) := by
  have hrw : p + (p - 2) = 2 * (p - 1) := by omega
  rw [hrw, mul_comm (p - 1) (2 * (p - 1))]
  exact Nat.mul_le_mul (two_mul_pred_le_pow hp hk) (le_refl _)

/-- ★★逆 —— `k ≤ 3` では `p ≥ 3` の**すべて**で破れる。⇒ `K₀ = 4` は最適。 -/
theorem cyclotomic_jump_fails {p k : ℕ} (hp : 3 ≤ p) (hk : k ≤ 3) :
    ¬ ((p - 1) * (p + (p - 2)) ≤ p ^ (k - 2) * (p - 1)) := by
  have hrw : p + (p - 2) = 2 * (p - 1) := by omega
  intro hcon
  rw [hrw, mul_comm (p - 1) (2 * (p - 1))] at hcon
  exact not_two_mul_pred_le_pow hp hk
    (Nat.le_of_mul_le_mul_right hcon (by omega : 0 < p - 1))

end Cyclotomic

/-! ## §2b ★`p = 2` の訂正 —— 「余白 0」と「基底段が要らない」は別のこと -/

section TwoCase

/-- ★`p = 2` では余白が本当に 0 —— 要る不等式が木の与える形と**同じ**になる。 -/
theorem slack_vanishes_at_two (t : ℕ) : (2 - 1) * (t + (2 - 2)) = (2 - 1) * t := by
  norm_num

/-- ★★それでも `k = 2` では破れる（`2 ≤ 1`）。
⇒ ★下書きの「`p = 2` では `K₀ = 0`」は**偽**だった。 -/
theorem two_fails_at_two : ¬ ((2 - 1) * (2 + (2 - 2)) ≤ 2 ^ (2 - 2) * (2 - 1)) := by
  norm_num

/-- ★`p = 2` は `k = 3` で入る（`2 ≤ 2`）。⇒ ★`p = 2` の最小 `K₀` は **3**。 -/
theorem two_holds_at_three : (2 - 1) * (2 + (2 - 2)) ≤ 2 ^ (3 - 2) * (2 - 1) := by
  norm_num

end TwoCase

/-! ## §3 ★★出口 —— 4 仮説がすべて「深い段だけ」になる -/

section Exit

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**損失評価の出口（深い段だけ版）**。

`LossToExit.axSenTate_of_loss_bound` との違い:

* `hc : ∀ k, 1 ≤ c k` が ★**消えた**（`BaseLayerConstant` の貼り合わせで自動）
* `hform` / `hm` / `he` / `hjump` がすべて ★**`K₀ ≤ k` 付き**になった

⇒ ★`LossToExit.lean:130` が「浅い層は別扱いが要る」と書いた点が、これで閉じる。 -/
theorem axLemma_of_deep_loss_bound (K : PAdicLocalField p) {c : ℕ → ℝ} {t m e : ℕ → ℕ}
    (K₀ : ℕ)
    (hform : ∀ k, K₀ ≤ k →
      c k ≤ (p : ℝ) ^ (((t k + (p - 2) : ℕ) : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, K₀ ≤ k → p ^ (k - 1) ≤ m k) (he : ∀ k, K₀ ≤ k → 0 < e k)
    (hjump : ∀ k, K₀ ≤ k → (p - 1) * (t k + (p - 2)) ≤ e k)
    (h : AxWildDescent K c) : AxLemma K ((p : ℝ) ^ K₀ * axConstant p) :=
  BaseLayerConstant.axLemma_of_deep_bound K K₀
    (fun k hk => (hform k hk).trans
      (LossToExit.axDecay_of_loss_bound (hm k hk) (he k hk) (hjump k hk))) h

/-- ★★そこから `AxSenTate K`。 -/
theorem axSenTate_of_deep_loss_bound (K : PAdicLocalField p) {c : ℕ → ℝ} {t m e : ℕ → ℕ}
    (K₀ : ℕ)
    (hform : ∀ k, K₀ ≤ k →
      c k ≤ (p : ℝ) ^ (((t k + (p - 2) : ℕ) : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, K₀ ≤ k → p ^ (k - 1) ≤ m k) (he : ∀ k, K₀ ≤ k → 0 < e k)
    (hjump : ∀ k, K₀ ≤ k → (p - 1) * (t k + (p - 2)) ≤ e k)
    (h : AxWildDescent K c) : AxSenTate K :=
  BaseLayerConstant.axSenTate_of_deep_bound K K₀
    (fun k hk => (hform k hk).trans
      (LossToExit.axDecay_of_loss_bound (hm k hk) (he k hk) (hjump k hk))) h

end Exit

/-! ## §4 ★★★`K₀ = 4` —— `hm` / `he` / `hjump` を**消した**形 -/

section Sharp

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**残っている仮説は 2 つだけ**。

* `hform`（★**深さ 4 以上でだけ**）—— 1 段の損失の形
* `h : AxWildDescent K c` —— 降下そのもの

`hm` は `rfl`、`he` は `p ≥ 3`、`hjump` は §2 で消えた。

★★これは円分塔の模型（`m k = p^{k−1}`, `e k = p^{k−2}(p−1)`, `t k = p`）を
**代入した**形である。一般の `K` には §3 を使う。 -/
theorem axSenTate_of_cyclotomic_loss (K : PAdicLocalField p) {c : ℕ → ℝ} (hp : 3 ≤ p)
    (hform : ∀ k, 4 ≤ k → c k ≤ (p : ℝ) ^ (((p + (p - 2) : ℕ) : ℝ) /
      (((p ^ (k - 1) : ℕ) : ℝ) * ((p ^ (k - 2) * (p - 1) : ℕ) : ℝ))))
    (h : AxWildDescent K c) : AxSenTate K :=
  axSenTate_of_deep_loss_bound K 4 (t := fun _ => p) (m := fun k => p ^ (k - 1))
    (e := fun k => p ^ (k - 2) * (p - 1)) hform (fun _ _ => le_refl _)
    (fun _ _ => mul_pos (pow_pos (by omega) _) (by omega))
    (fun k hk => cyclotomic_jump_fits hp hk) h

end Sharp

/-! ## §5 使っている公理の一覧 -/

#print axioms two_mul_pred_le_pow
#print axioms not_two_mul_pred_le_pow
#print axioms cyclotomic_jump_fits
#print axioms cyclotomic_jump_fails
#print axioms slack_vanishes_at_two
#print axioms two_fails_at_two
#print axioms two_holds_at_three
#print axioms axLemma_of_deep_loss_bound
#print axioms axSenTate_of_deep_loss_bound
#print axioms axSenTate_of_cyclotomic_loss

end DeepLossExit

end ABC3.Found.PGC
