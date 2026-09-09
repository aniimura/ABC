import ABC3.Found.PGC.BoundedJumpSum

/-!
# [pGC] ★★私の損失評価が木の出口に**そのまま入った** —— 「跳びが有界」は要らなかった

## 持ち場（前波で「残るは 1 本」とした点）

`t k ≤ T`（第 1 跳びが `k` に依らない）の供給。

## ★★開いたら、その 1 本は**要らなかった**

木の `FirstJumpRoute.axLemma_of_firstJump`（`NormalizedTraceDescent.lean:379`）の仮定は

```
(hform : ∀ k, c k ≤ (p:ℝ) ^ ((j k : ℝ) / ((m k : ℝ) * (e k : ℝ))))
(hm : ∀ k, p ^ (k - 1) ≤ m k) (he : ∀ k, 0 < e k)
(hjump : ∀ k, (p - 1) * j k ≤ e k)
```

★★**`j k` も `e k` も `k` とともに動いてよい。** 総和が閉じるのは
`hm : m k ≥ p^{k−1}` だけからである。⇒ ★私が前波で「要る」と書いた
**「跳びが `k` に依らない `T` で押さえられる」は不要**だった。
★要るのは各段ごとの `(p−1)·j k ≤ e k` —— ★これは `RamificationJumpBound` が与える形そのもの。

## ★何を差し込んだか

私の 1 段の損失は `t + p − 2`（`t` = 第 1 跳び）。★`j k := t k + (p−2)` と置くだけで
木の 3 定理にそのまま入った（§1・§2、いずれも**本体 1 行**）。

| 節 | 宣言 | 対応する木の定理 |
|---|---|---|
| §1 | `axDecay_of_loss_bound` | `JumpArith.rpow_div_le_axDecay` |
| §2 | `axLemma_of_loss_bound` | `FirstJumpRoute.axLemma_of_firstJump` |
| §2 | `axSenTate_of_loss_bound` | `FirstJumpRoute.axSenTate_of_firstJump` |

## ★★★残っているものは 1 つの不等式に落ちた

★**`(p−1)·(t k + p − 2) ≤ e k`**。

木が与えるのは `(p−1)·t k ≤ e k`（`RamificationJumpBound`）なので、
★差は `(p−1)(p−2)` の**余白**である。円分塔（`t = p`、`e = e_{E₁} = p^{k−2}(p−1)`）に入れると
`2(p−1) ≤ p^{k−2}` になり、

* `p = 3, k = 3`: `4 ≤ 3` ★**成り立たない**（§3 `fit_fails_at_three`）
* `p = 3, k = 4`: `4 ≤ 9` 成り立つ（§3 `fit_holds_at_four`）

⇒ ★**深い層では入るが、浅い層（`k` が小さいところ）は別扱いが要る。**
★これは「基底段を別に処理する」ふつうの形だが、★**本日の鎖にはその処理が無い**。
★次の 1 点はここである。

## ★配管（1 往復で当たった罠）

`NormalizedTraceDescent.lean` は ★**`NormalizedTraceDescent` という名前空間を開いていない**。
中身は `ABC3.Found.PGC.JumpArith` と `ABC3.Found.PGC.FirstJumpRoute` に在る。逐語:

```
error: Unknown identifier `NormalizedTraceDescent.JumpArith.rpow_div_le_axDecay`
```

★**ファイル名 ≠ 名前空間**。`grep -n "^namespace" <file>` で先に確かめる。

## 逸脱の記録

- §1・§2 は木の定理への**代入**であり、★新しい数学は 0。★価値は「私の `t + p − 2` が
  木の `j` の位置にちょうど嵌まる」ことを**確かめた**ことにある。
- §3 の 2 本は手計算（`norm_num`）の記録である。
-/

namespace ABC3.Found.PGC

namespace LossToExit

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 ★私の `loss ≤ t + p − 2` が木の出口の形にそのまま入る -/

section Fit

/-- ★★私の 1 段の損失 `t + p − 2`（`t` = 第 1 跳び）を、木の
`JumpArith.rpow_div_le_axDecay` に**そのまま**入れた形。

★`j := t + (p−2)` と置くだけ。★要る仮定は
`m ≥ p^{k−1}`（塔の下から `k−1` 段ぶん）と `(p−1)·(t + p − 2) ≤ e` の 2 本。 -/
theorem axDecay_of_loss_bound {k t m e : ℕ}
    (hm : p ^ (k - 1) ≤ m) (he : 0 < e)
    (hjump : (p - 1) * (t + (p - 2)) ≤ e) :
    (p : ℝ) ^ (((t + (p - 2) : ℕ) : ℝ) / ((m : ℝ) * (e : ℝ))) ≤ axDecay p k :=
  JumpArith.rpow_div_le_axDecay hm he hjump

end Fit

/-! ## §2 ★出口 —— `AxLemma` と `AxSenTate` -/

section Exit

/-- ★★★**私の損失評価を木の出口に差し込んだ形**。

`c k ≤ p^{(t k + p − 2)/(m k · e k)}` と 3 本の算術条件から `AxLemma K (axConstant p)`。
★中身は木の `FirstJumpRoute.axLemma_of_firstJump` に `j k := t k + (p−2)` を代入しただけ。 -/
theorem axLemma_of_loss_bound (K : PAdicLocalField p) {c : ℕ → ℝ} {t m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p : ℝ) ^ (((t k + (p - 2) : ℕ) : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, p ^ (k - 1) ≤ m k) (he : ∀ k, 0 < e k)
    (hjump : ∀ k, (p - 1) * (t k + (p - 2)) ≤ e k)
    (h : AxWildDescent K c) : AxLemma K (axConstant p) :=
  FirstJumpRoute.axLemma_of_firstJump K hc hform hm he hjump h

/-- ★★そこから `AxSenTate K`。 -/
theorem axSenTate_of_loss_bound (K : PAdicLocalField p) {c : ℕ → ℝ} {t m e : ℕ → ℕ}
    (hc : ∀ k, 1 ≤ c k)
    (hform : ∀ k, c k ≤ (p : ℝ) ^ (((t k + (p - 2) : ℕ) : ℝ) / ((m k : ℝ) * (e k : ℝ))))
    (hm : ∀ k, p ^ (k - 1) ≤ m k) (he : ∀ k, 0 < e k)
    (hjump : ∀ k, (p - 1) * (t k + (p - 2)) ≤ e k)
    (h : AxWildDescent K c) : AxSenTate K :=
  FirstJumpRoute.axSenTate_of_firstJump K hc hform hm he hjump h

end Exit

/-! ## §3 ★`+(p−2)` の余白は**浅い層では足りない**（手計算の記録） -/

section Slack

/-- ★木が与えるのは `(p−1)·t ≤ e`。私が要るのは `(p−1)·(t + p − 2) ≤ e` で、
差は `(p−1)(p−2)` である。円分塔（`t = p`、`e = e_{E₁} = p^{k−2}(p−1)`）で書くと
`2(p−1) ≤ p^{k−2}` になる。

★`p = 3, k = 3` では `4 ≤ 3` で**成り立たない**。 -/
theorem fit_fails_at_three : ¬ (2 * (3 - 1) ≤ 3 ^ (3 - 2)) := by norm_num

/-- ★`p = 3, k = 4` では `4 ≤ 9` で成り立つ。⇒ ★**深い層では入るが、浅い層は別扱いが要る**。 -/
theorem fit_holds_at_four : 2 * (3 - 1) ≤ 3 ^ (4 - 2) := by norm_num

end Slack

/-! ## §4 使っている公理の一覧 -/

#print axioms axDecay_of_loss_bound
#print axioms axLemma_of_loss_bound
#print axioms axSenTate_of_loss_bound
#print axioms fit_fails_at_three

end LossToExit

end ABC3.Found.PGC
