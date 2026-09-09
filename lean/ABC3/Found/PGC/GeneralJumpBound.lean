import ABC3.Found.PGC.JumpAtP

/-!
# [pGC] 一般の跳び —— ★この道は**円分塔でだけ**総和が閉じる（境界を測った）

## 持ち場（前波で私が「まだ測っていない」とした点）

跳び `t` が `p` でない層。★定数は `t + p − 2`（前波 `JumpAtP.exponent_bound_general`）。

## ★★測った境界

1 段の損失を ★**絶対付値**（`v_L` を `e_L` で割った、体に依らない目盛り）で測る。
`e_L = p·e_E` である。

| | 跳び `t` | 1 段の絶対損失 | 段を重ねると |
|---|---|---|---|
| 円分の層 | `t = p`（測定 (n1)、5 設定すべて） | `(2p−2)/(p·e_E)` ≈ `2/e_E` | ★`e_E` とともに**幾何級数的に縮む** |
| 一般の層 | `t ≤ p·e_E/(p−1)`（★古典的上界、Serre IV） | `≤ 1/(p−1) + (p−2)/(p·e_E)` | ★**第 1 項が `e_E` に依らない** |

⇒ ★★**一般の層では 1 段あたり `1/(p−1)` が縮まずに残る。**
`AxLemma K C` は `K.closure` の**すべての** `x` について 1 つの `C` を要求し、
`wildDepth x` はいくらでも大きくなるので、★**この道からは一様な `C` が出ない**
（§3 `no_uniform_constant`、中身は §1 の Archimedes）。

## ★★★これが意味すること（★誇張しないで書く）

* ★本日の鎖（`loss ≤ 2p−2` → `SharpPrizeBound.sharp_exponent_sum_le` → `axLemma_of_wildDescent`）は
  ★**円分塔については閉じる**（`t = p` で幾何級数、`SharpPrizeBound` が総和を押さえている）。
* ★**一般の `x`（＝ `AxWildDescent` の主張そのもの）には、この見積もりでは足りない。**
  ★足りないのは「1 段の損失」ではなく ★**「深い層でも跳びが小さい」ことの理由**である。
* ⇒ 出口を閉じるには、★**跳びが大きい層では別の見積もり**（古典の Ax–Sen–Tate は
  trace/different を使う）か、★**降下の際に `x′` を選んで跳びを小さく保つ**議論が要る。
  ★どちらも本日の鎖には無い。★**まだ測っていない**。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `exists_nat_mul_gt` | ★抽象核（Archimedes）。分岐も付値も出ない |
| §2 | `cyclotomic_fraction_eq` | 円分の層: `p/(p·e_E) = 1/e_E` |
| §2 | `absolute_loss_le_of_jump_bound` | 一般の層: `≤ 1/(p−1) + (p−2)/(p·e_E)` |
| §3 | `no_uniform_constant` | ★★一様な `C` は取れない |

## 逸脱の記録

- ★`t ≤ p·e_E/(p−1)`（Serre の古典的上界）は**仮説として受けた**。
  ★形式化していない。★「測って無かった」ではなく「**古典を引いた**」である。
- ★前波のファイル名は `JumpAtP.lean` である（本体の持ち場は `JumpGeneralization.lean` と
  書いていたが、そのファイルは**存在しない**。`ls` で確認した）。★台帳の訂正のため記す。
-/

namespace ABC3.Found.PGC

namespace GeneralJumpBound

/-! ## §1 抽象核 —— 正の定数を足し続ければ必ず超える -/

section Kernel

/-- ★抽象核（Archimedes）。`0 < c` なら `n·c` はどんな `C` も超える。
★分岐も付値も出てこない。★これが「一様な定数 `C` が取れない」ことの本体である。 -/
theorem exists_nat_mul_gt {c : ℝ} (hc : 0 < c) (C : ℝ) : ∃ n : ℕ, C < n * c := by
  obtain ⟨n, hn⟩ := exists_nat_gt (C / c)
  exact ⟨n, by rwa [div_lt_iff₀ hc] at hn⟩

end Kernel

/-! ## §2 1 段の損失を**絶対付値**で測る -/

section Fraction

/-- ★円分の層（跳び `t = p`）: 1 段の損失の**絶対**付値は `1/e_E` —— ★`e_E` とともに縮む。 -/
theorem cyclotomic_fraction_eq {p eE : ℝ} (hp : p ≠ 0) (he : eE ≠ 0) :
    p / (p * eE) = 1 / eE := by
  field_simp

/-- ★★一般の層: 跳びの古典的な上界 `t ≤ p·e_E/(p−1)` を入れると、
1 段の損失の**絶対**付値は `1/(p−1) + (p−2)/(p·e_E)` 以下。

★第 1 項 `1/(p−1)` は ★**`e_E` に依らない**。⇒ 段を重ねても縮まない。 -/
theorem absolute_loss_le_of_jump_bound {p eE t : ℝ} (hp : 2 ≤ p) (he : 0 < eE)
    (ht : t ≤ p * eE / (p - 1)) :
    (t + (p - 2)) / (p * eE) ≤ 1 / (p - 1) + (p - 2) / (p * eE) := by
  have hp0 : (0 : ℝ) < p := by linarith
  have hpe : (0 : ℝ) < p * eE := mul_pos hp0 he
  have hp1 : (0 : ℝ) < p - 1 := by linarith
  have hfirst : t / (p * eE) ≤ 1 / (p - 1) := by
    rw [div_le_div_iff₀ hpe hp1]
    have h1 : t * (p - 1) ≤ (p * eE / (p - 1)) * (p - 1) :=
      mul_le_mul_of_nonneg_right ht (le_of_lt hp1)
    rw [div_mul_cancel₀ _ (ne_of_gt hp1)] at h1
    linarith
  calc (t + (p - 2)) / (p * eE) = t / (p * eE) + (p - 2) / (p * eE) := by ring
    _ ≤ 1 / (p - 1) + (p - 2) / (p * eE) := by linarith

end Fraction

/-! ## §3 ★★結論 —— この道からは**一様な定数**が出ない -/

section NoUniform

/-- ★★★1 段あたり `1/(p−1)` を失う道では、★**深さに依らない一様な `C` は取れない**。

`AxLemma K C` は `K.closure` の**すべての** `x` について 1 つの `C` を要求する。
`wildDepth x` は `x` とともにいくらでも大きくなるので、
★段ごとの損失が縮まないなら合計は発散する（§1 の Archimedes）。 -/
theorem no_uniform_constant {p : ℝ} (hp : 2 ≤ p) (C : ℝ) :
    ∃ n : ℕ, C < n * (1 / (p - 1)) := by
  refine exists_nat_mul_gt ?_ C
  have : (0 : ℝ) < p - 1 := by linarith
  positivity

end NoUniform

/-! ## §4 使っている公理の一覧 -/

#print axioms exists_nat_mul_gt
#print axioms cyclotomic_fraction_eq
#print axioms absolute_loss_le_of_jump_bound
#print axioms no_uniform_constant

end GeneralJumpBound

end ABC3.Found.PGC
