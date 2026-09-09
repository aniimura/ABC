import ABC3.Found.PGC.GeneralJumpBound

/-!
# [pGC] ★★★**前波の自分の否定的結論を撤回する** —— 「第 1 跳び」は上の層の跳びではない

## 何が起きたか

前波（`GeneralJumpBound.lean`）で私はこう書いた:

> ★★**一般の層では 1 段あたり `1/(p−1)` が縮まずに残る。**
> ★**この道からは一様な `C` が出ない**

★★**これは誤りである。** 撤回する。★誤りの中身は ★**跳びをどの層のものと取るか**だった。

## ★木を読んで分かったこと（`NormalizedTraceDescent.lean` §2/§3）

`rpow_div_le_axDecay`（`:224`）の docstring はこう言っている:

> ★★分岐の言葉に戻すと `j` = 下の層 `E₁/F` の跳び、`m` = `e(L/E₁)`、`e` = `e_{E₁}`、
> `m·e = e_L` である。

★つまり木の勘定は ★**下の層の跳び `j`** を **上の `e_L = m·e`** で割る。
`m ≥ p^{k−1}` なので ★**指数は幾何級数的に縮む**。
`hjump : (p−1)·j ≤ e` は `RamificationJumpBound.lean`（★木に**在る**、561 行・`sorry` 0）の
主結果を下の層に当てた形である。

## ★★測って確かめた（`tools/zeta-tower-check.py`、`ℚ₃(ζ₂₇)`、`e_L = 18`）

```
Gal(L / Q_3(zeta_3^1)) : size 9,  min i = 3   ← ★H_K の第 1 跳び
Gal(L / Q_3(zeta_3^2)) : size 3,  min i = 9   ← 上の層 L/E₁ の跳び
```

★★**この 2 つは違う（3 と 9）。** そして ★**私の鎖の `t = v_L(σπ − π)` は前者**である ——
`hform-k3-check.py` の `HK = {a : a ≡ 1 mod p, a ≠ 1}` は `Gal(L/ℚ_p(ζ_p))` であって
`Gal(L/E₁)` ではなく、`i₁` はその**最小**（＝第 1 跳び）を取っている。
⇒ ★**私が測った `t = p` は「下の側の跳び」であり、`e_L` とともに増えない。**

## ★★訂正後の結論

| | 前波の私 | ★訂正 |
|---|---|---|
| `t` の正体 | 上の層 `L/E₁` の跳び | ★**`H_K` の第 1 跳び**（下の側） |
| `t` の上界 | `e_L/(p−1)`（`e_L` とともに増える） | ★**`k` に依らない**（円分塔では `t = p`、5 設定で測定） |
| 1 段の絶対損失 | `1/(p−1)`（縮まない） | ★`(t+p−2)/(p^{k}(p−1))`（**幾何級数**） |
| 総和 | 発散 | ★**収束**（§2、`SharpPrizeBound` と同じ形） |

★前波の `no_uniform_constant`（`GeneralJumpBound.lean` §3）は**定理としては正しい**
（`1/(p−1)` を足し続ければ発散する）が、★**その前提「1 段あたり `1/(p−1)` 失う」が誤り**
だった。★あちらの docstring は直さない（規約）。★ここに訂正として書く。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `first_jump_ne_layer_jump` / `first_jump_eq_p` | ★測定の記録（`3 ≠ 9`） |
| §2 | `sum_le_of_fixed_constant` | ★跳びが `k` に依らなければ総和は閉じる |
| §2 | `two_p_sub_two_recovers` | `c = 2p−2` で `SharpPrizeBound` の `2p/(p−1)` に一致 |

## ★残っている本当の 1 点

★**「`H_K` の第 1 跳びが `k` に依らない」ことの証明**である。円分塔では測定で `t = p`
（5 設定）だが、一般の `K` では `RamificationJumpBound` の
`(p−1)·j ≤ e_{E₁}` が要る形になっている（`e_{E₁}` は下の層の値なので `k` に依らない）。
★木の `RamificationJumpBound.lean` を**私はまだ読んでいない**。★費用は書かない。

## 逸脱の記録

- ★本ファイルは**私自身の前波を撤回する**ためのものである。★12 度目の自己訂正。
- §1 の 2 本は測定の記録（`norm_num` / `rfl`）であって数学の主張ではない。
-/

namespace ABC3.Found.PGC

namespace FirstJumpIsBottom

open Finset

/-! ## §1 ★測定の記録 —— 「第 1 跳び」は**上の層の跳びではない** -/

section Record

/-- ★★`ℚ₃(ζ₂₇)`（`e_L = 18`）で `tools/zeta-tower-check.py` が出した 2 つの数:

```
Gal(L / Q_3(zeta_3^1)) : size 9,  min i = 3   ← ★H_K の第 1 跳び
Gal(L / Q_3(zeta_3^2)) : size 3,  min i = 9   ← 上の層 L/E₁ の跳び
```

★**この 2 つは違う。** 私の鎖の `t = v_L(σπ − π)` は**前者**（`H_K` の第 1 跳び、`= p`）である。 -/
theorem first_jump_ne_layer_jump : (3 : ℕ) ≠ 9 := by norm_num

/-- ★`H_K` の第 1 跳びは `p` に等しい（`p = 3` の測定）。★上の層の跳び `9` は `e_L/2` である。 -/
theorem first_jump_eq_p : (3 : ℕ) = 3 := rfl

end Record

/-! ## §2 ★訂正後の総和 —— 跳びが `k` に依らなければ閉じる -/

section Sum

/-- ★★1 段の定数が `c/(p^k(p−1))` の形（`c` が `k` に**依らない**）なら、総和は閉じる。

★`c = 2p−2` が `SharpPrizeBound.sharp_exponent_sum_le` の場合である。 -/
theorem sum_le_of_fixed_constant {p c : ℝ} (hp : 2 ≤ p) (hc : 0 ≤ c) (n : ℕ) :
    ∑ k ∈ range n, c / (p ^ k * (p - 1)) ≤ c * p / ((p - 1) * (p - 1)) := by
  have hp1 : (1 : ℝ) < p := by linarith
  have hp0 : (0 : ℝ) < p := by linarith
  have hpm : (0 : ℝ) < p - 1 := by linarith
  have hgeom : ∑ k ∈ range n, ((1 : ℝ) / p) ^ k ≤ p / (p - 1) := by
    have h := SharpPrizeBound.partial_geom_le (r := (1 : ℝ) / p) (by positivity)
      (by rw [div_lt_one hp0]; linarith) n
    have hval : (1 : ℝ) / (1 - 1 / p) = p / (p - 1) := by field_simp
    rwa [hval] at h
  have hrw : ∑ k ∈ range n, c / (p ^ k * (p - 1))
      = (c / (p - 1)) * ∑ k ∈ range n, ((1 : ℝ) / p) ^ k := by
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [div_pow, one_pow]
    field_simp
  rw [hrw]
  calc (c / (p - 1)) * ∑ k ∈ range n, ((1 : ℝ) / p) ^ k
      ≤ (c / (p - 1)) * (p / (p - 1)) := mul_le_mul_of_nonneg_left hgeom (by positivity)
    _ = c * p / ((p - 1) * (p - 1)) := by rw [div_mul_div_comm]

/-- ★`c = 2p−2` を入れると `SharpPrizeBound` の `2p/(p−1)` に一致する。 -/
theorem two_p_sub_two_recovers {p : ℝ} (hp : 2 ≤ p) :
    (2 * p - 2) * p / ((p - 1) * (p - 1)) = 2 * p / (p - 1) := by
  have hpm : (p : ℝ) - 1 ≠ 0 := by
    intro h
    linarith [h]
  have h2 : (2 : ℝ) * p - 2 = 2 * (p - 1) := by ring
  rw [h2]
  field_simp

end Sum

/-! ## §3 使っている公理の一覧 -/

#print axioms first_jump_ne_layer_jump
#print axioms sum_le_of_fixed_constant
#print axioms two_p_sub_two_recovers

end FirstJumpIsBottom

end ABC3.Found.PGC
