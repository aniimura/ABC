import ABC3.Found.PGC.MaxMinIndex

/-!
# [pGC] 成分の構成 —— `σx − x = ρ·(形式微分) + O(‖ρ‖²)` から `B·π^{j₀} + R` まで

## 持ち場

前波で私は「残りはちょうど 1 点: 成分そのものの構成（`σx − x` の `𝒪_{E₁}` 基底展開
`Σ B_j π^j` と `B_{j₀}` の主部）」と書いた。本ファイルがそれを埋める。

## ★★★前波の自分の判定を**精密化する（訂正）**

`ComponentFormulaScope.lean` で私は「手で導いた成分の公式は**一般には偽**、`j = j₀` でだけ真」
と書いた。★**この言い方は 2 つの別の主張を混ぜていた。**

| | 主張 | 真偽 |
|---|---|---|
| (α) | **大域**の形: `σx − x = w·Σ_j ( j f_j + (j+1) f_{j+1} )·π^j + E`、`‖E‖ ≤ ‖ρ‖²` | ★**常に真**（本ファイル §3 で**証明**した。仮定は整数性と `ρ = w(1+π)` だけ） |
| (β) | **成分ごと・相対**の形: 各 `j` で `v(B_j − P_j) > v(P_j)` | ★**一般には偽**（前波の測定どおり） |

★(β) が偽になるのは、`P_j` が深いスロットでは**絶対誤差 `‖ρ‖²` の方が `P_j` より大きい**から。
`j = j₀` では `P_{j₀}` が最上位に来るので `‖ρ‖²` が下に落ち、(β) も真になる。
⇒ ★**「公式が偽」なのではなく、「相対誤差で測ると偽、絶対誤差で測ると真」**だった。
本ファイルは (α) を証明したので、★以後は (β) を使う必要がない。

## ★何が抽象核だったか

「成分の公式」の正体は ★**形式微分**だった:

```
Σ_j f_j·( (π+ρ)^j − π^j )  =  ρ · Σ_j j f_j π^{j−1}  +  O(‖ρ‖²)      （§1）
(1+π) · Σ_j j f_j π^{j−1}  =  Σ_j ( j f_j + (j+1) f_{j+1} )·π^j       （§2、可換環の恒等式）
```

⇒ `ρ = w·(1+π)`（`RhoFactorization`）を代入すると主部が出る（§3）。
★§1・§2 とも**分岐・付値・Galois・素数が 1 語も出ない**。

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `norm_sum_sub_deriv_le` | ★抽象核 1。差は `ρ×形式微分` から `‖ρ‖²` しかずれない |
| §2 | `sum_deriv_shift` | ★抽象核 2。`(1+π)×形式微分` の係数ずらし（境界項つき、可換環の恒等式） |
| §3 | `norm_expansion_main_le` | 2 つを合わせる（`f p = 0` で境界項が消える） |
| §4 | `expansion_main_at_max_min` | `j₀` を取り出して `B·π^{j₀} + R` の形に。`‖B‖ = ‖w‖·‖f_{j₀}‖` |
| §5 | `norm_moved_le` | ★抽象核 3。係数が動いてもよい形 |
| §5 | `norm_expansion_main_le_moved` / `expansion_main_at_max_min_moved` | ★塔の形（到達点） |

## ★★凍結模型は塔では**偽** —— それでも結論は変わらない（測定つき）

§1〜§4 は `σx − x = Σ_j f_j((π+ρ)^j − π^j)`、すなわち ★**`σ` が `𝒪_{E₁}` 係数 `f_j` を
動かさない**という模型である。★塔ではこれは**偽**: `σ ∈ H_K` は `a ≡ 1 (mod p)` しか
満たさないので `n ≥ 3` では `E₁ = ℚ_p(ζ_{p^{n−1}})` を動かす。
`tools/expansion-check.py` の (m5) が **360/360 件で「動いた」**と言っている
（`p=7, n=2` は `E₁ = ℚ_p` なので 84/84 件で動かない）。

⇒ ★§5 で模型を外した。動く量が `‖ρ‖²` 以下なら結論は同じ、という形にした。
その仮定 `hmv` 自体を測ったのが (m5b): ★**504/504 件で成立**。

## ★測定（`tools/expansion-check.py`、厳密整数演算）

| p | n | 標本 | v(ρ) | (m1) `ρ = w(1+π)` | (m2) §3 の主張 | (m4) 塔での結論 | (m5b) `hmv` |
|---|---|---|---|---|---|---|---|
| 3 | 3 | 60 | 3 | 成立 | 60/60（最小 v = 6 = 2v(ρ)、★sharp） | 60/60 | 180/180（最小 9） |
| 3 | 4 | 20 | 3 | 成立 | 20/20（最小 6） | 20/20 | 60/60（最小 9） |
| 5 | 3 | 12 | 5 | 成立 | 12/12（最小 10） | 12/12 | 60/60（最小 25） |
| 7 | 2 | 12 | 7 | 成立 | 12/12（最小 14） | 12/12 | 84/84（∞） |
| 2 | 4 | 60 | 2 | 成立 | 60/60（E = 0） | 60/60（最小 4） | 120/120（★最小 4 = 2v(ρ)） |

★`‖E‖ ≤ ‖ρ‖²` は**ちょうど等号で当たる**（最小 v = 2v(ρ)）。緩めれば落ちる。
★`p = 2` では `hmv` も等号で当たる（`v(σf_j − f_j) = 2v(ρ)`）。ここでも `p = 2` だけが境界。

## 逸脱の記録

- §1〜§4 は「凍結模型」の主張である。★塔で使うときは §5 の形（`norm_expansion_main_le_moved`
  / `expansion_main_at_max_min_moved`）を使うこと。§1〜§4 を直接引くと**塔では成り立たない
  仮定**を引くことになる。
- `x = Σ_{j<p} f_j π^j` という**展開の存在そのもの**は仮定である（`f` を引数で受けている）。
  材料は `TotallyRamifiedLayer.adjoin_eq_top_of_valK` / `linearIndependent_of_ne_mod`。
- ★残る数合わせ: `‖B‖ = ‖w‖·‖f_{j₀}‖` を `ResidueSeparation.loss_le_of_residue_ne` が要求する
  `‖B‖ = ‖π‖^{d+(p−1)}` に読み替えること。`v(w) = p = i₁ + 1` は上の表で測ってある
  （5 つの設定すべてで `v(w) = v(ρ) = i₁ + 1`）。
-/

namespace ABC3.Found.PGC

namespace ComponentExpansion

open Finset

/-! ## §1 抽象核 1 —— `σx − x` の主部は「形式微分 × ρ」 -/

section FirstOrder

/-- ★抽象核。`x = Σ f_j π^j` を `π ↦ π + ρ` で動かした差は、
`ρ ·（形式微分）` から `‖ρ‖²` しかずれない。★分岐・付値・Galois が 1 語も出ない。 -/
theorem norm_sum_sub_deriv_le {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (f : ℕ → M) (hf : ∀ j, ‖f j‖ ≤ 1) (n : ℕ) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1)‖ ≤ ‖ρ‖ ^ 2 := by
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
  intro j _
  match j with
  | 0 => simp
  | (m + 1) =>
      have hrw : f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1))
            - ρ * ((((m + 1 : ℕ)) : M) * f (m + 1) * π ^ (m + 1 - 1))
          = f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1) - ((m : M) + 1) * π ^ m * ρ) := by
        simp only [Nat.add_sub_cancel]
        push_cast
        ring
      rw [hrw, norm_mul]
      calc ‖f (m + 1)‖ * ‖(π + ρ) ^ (m + 1) - π ^ (m + 1) - ((m : M) + 1) * π ^ m * ρ‖
          ≤ 1 * ‖ρ‖ ^ 2 := by
            refine mul_le_mul (hf _) (BinomialFirstOrder.norm_add_pow_sub_linear_le hπ hρ m)
              (norm_nonneg _) zero_le_one
        _ = ‖ρ‖ ^ 2 := one_mul _

end FirstOrder

/-! ## §2 抽象核 2 —— `(1 + π) ×（形式微分）`の係数のずらし -/

section Shift

/-- ★抽象核。可換環の恒等式。`(1+π)` を形式微分に掛けると
係数が `j f_j + (j+1) f_{j+1}` にずれる（境界項つき）。 -/
theorem sum_deriv_shift {M : Type*} [CommRing M] (π : M) (f : ℕ → M) (n : ℕ) :
    ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = (1 + π) * (∑ j ∈ range n, (j : M) * f j * π ^ (j - 1))
        + (n : M) * f n * π ^ (n - 1) := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [Finset.sum_range_succ, ih, Finset.sum_range_succ, mul_add]
      match m with
      | 0 => simp
      | (k + 1) =>
          simp only [Nat.add_sub_cancel]
          push_cast
          ring

end Shift

/-! ## §3 2 つの核を合わせる —— 主部は `w · Σ (j f_j + (j+1) f_{j+1}) π^j` -/

section Main

theorem norm_expansion_main_le {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f : ℕ → M) (hf : ∀ j, ‖f j‖ ≤ 1) {n : ℕ} (hfn : f n = 0) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j‖
      ≤ ‖ρ‖ ^ 2 := by
  have key : w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1) := by
    rw [sum_deriv_shift, hfn, hw]
    ring
  rw [key]
  exact norm_sum_sub_deriv_le hπ hρ f hf n

end Main

/-! ## §4 `j₀` を取り出す —— `B · π^{j₀} + R` の形 -/

section AtJZero

theorem expansion_main_at_max_min {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} {p : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1)
    (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f : ℕ → M) (hf : ∀ j, ‖f j‖ ≤ 1) (hfp : f p = 0)
    {a : ℕ} (ha1 : 1 ≤ a) (ha2 : a ≤ p - 1) (hane : f a ≠ 0) :
    ∃ j₀ B R, 1 ≤ j₀ ∧ j₀ ≤ p - 1 ∧ ¬ p ∣ j₀ ∧ ‖B‖ = ‖w‖ * ‖f j₀‖ ∧
      ‖(∑ j ∈ range p, f j * ((π + ρ) ^ j - π ^ j)) - (B * π ^ j₀ + R)‖ ≤ ‖ρ‖ ^ 2 := by
  classical
  obtain ⟨j₀, h1, h2, hnd, _, _, hpair⟩ :=
    MaxMinIndex.exists_no_cancel_index hp hple f hfp ha1 ha2 hane
  refine ⟨j₀, w * ((j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)),
    w * ∑ j ∈ (range p).erase j₀, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j,
    h1, h2, hnd, ?_, ?_⟩
  · rw [norm_mul, hpair]
  · have hmem : j₀ ∈ range p := Finset.mem_range.mpr (by omega)
    have hsplit :
        w * ∑ j ∈ range p, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
          = w * ((j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)) * π ^ j₀
            + w * ∑ j ∈ (range p).erase j₀,
                ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j := by
      rw [← Finset.add_sum_erase _ _ hmem, mul_add]
      ring
    rw [← hsplit]
    exact norm_expansion_main_le hπ hρ hw f hf hfp

end AtJZero

/-! ## §5 ★凍結模型を外す —— `σ` が係数を動かしてもよい形 -/

section Moved

/-- ★抽象核。係数が `f` から `g` に動いても、動いた量が `C` 以下なら
「凍結模型」との差は `C` 以下。 -/
theorem norm_moved_le {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ : M} (hs : ‖π + ρ‖ ≤ 1) (f g : ℕ → M) {C : ℝ} (hC : 0 ≤ C)
    (h : ∀ j, ‖g j - f j‖ ≤ C) (n : ℕ) :
    ‖((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - ∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j)‖ ≤ C := by
  have hrw : ((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - ∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j)
      = ∑ j ∈ range n, (g j - f j) * (π + ρ) ^ j := by
    rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [hrw]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC ?_
  intro j _
  rw [norm_mul, norm_pow]
  calc ‖g j - f j‖ * ‖π + ρ‖ ^ j ≤ C * 1 := by
        refine mul_le_mul (h j) (pow_le_one₀ (norm_nonneg _) hs) (by positivity) hC
    _ = C := mul_one C

/-- ★塔の形。`σ` が `𝒪_{E₁}` 係数を動かしても、動きが `‖ρ‖²` 以下なら結論は同じ。 -/
theorem norm_expansion_main_le_moved {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hs : ‖π + ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f g : ℕ → M) (hf : ∀ j, ‖f j‖ ≤ 1) (hmv : ∀ j, ‖g j - f j‖ ≤ ‖ρ‖ ^ 2)
    {n : ℕ} (hfn : f n = 0) :
    ‖((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j‖
      ≤ ‖ρ‖ ^ 2 := by
  have hsplit : ((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = (((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
          - ∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        + ((∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
          - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j) := by
    ring
  rw [hsplit]
  refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
  · exact norm_moved_le hs f g (by positivity) hmv n
  · exact norm_expansion_main_le hπ hρ hw f hf hfn

/-- ★本ファイルの到達点。`σx − x` を `B·π^{j₀} + R` の形に書き、`‖B‖ = ‖w‖·‖f_{j₀}‖`
（＝先頭の対が打ち消さない）まで込みで出す。`σ` が係数を動かす場合を含む。 -/
theorem expansion_main_at_max_min_moved {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} {p : ℕ} (hp : p.Prime) (hple : ‖(p : M)‖ < 1)
    (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hs : ‖π + ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f g : ℕ → M) (hf : ∀ j, ‖f j‖ ≤ 1) (hmv : ∀ j, ‖g j - f j‖ ≤ ‖ρ‖ ^ 2)
    (hfp : f p = 0) {a : ℕ} (ha1 : 1 ≤ a) (ha2 : a ≤ p - 1) (hane : f a ≠ 0) :
    ∃ j₀ B R, 1 ≤ j₀ ∧ j₀ ≤ p - 1 ∧ ¬ p ∣ j₀ ∧ ‖B‖ = ‖w‖ * ‖f j₀‖ ∧
      ‖((∑ j ∈ range p, g j * (π + ρ) ^ j) - ∑ j ∈ range p, f j * π ^ j)
          - (B * π ^ j₀ + R)‖ ≤ ‖ρ‖ ^ 2 := by
  classical
  obtain ⟨j₀, h1, h2, hnd, _, _, hpair⟩ :=
    MaxMinIndex.exists_no_cancel_index hp hple f hfp ha1 ha2 hane
  refine ⟨j₀, w * ((j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)),
    w * ∑ j ∈ (range p).erase j₀, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j,
    h1, h2, hnd, ?_, ?_⟩
  · rw [norm_mul, hpair]
  · have hmem : j₀ ∈ range p := Finset.mem_range.mpr (by omega)
    have hsplit :
        w * ∑ j ∈ range p, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
          = w * ((j₀ : M) * f j₀ + ((j₀ + 1 : ℕ) : M) * f (j₀ + 1)) * π ^ j₀
            + w * ∑ j ∈ (range p).erase j₀,
                ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j := by
      rw [← Finset.add_sum_erase _ _ hmem, mul_add]
      ring
    rw [← hsplit]
    exact norm_expansion_main_le_moved hπ hρ hs hw f g hf hmv hfp

end Moved

/-! ## §6 使っている公理の一覧 -/

#print axioms norm_sum_sub_deriv_le
#print axioms sum_deriv_shift
#print axioms norm_expansion_main_le
#print axioms expansion_main_at_max_min
#print axioms norm_moved_le
#print axioms norm_expansion_main_le_moved
#print axioms expansion_main_at_max_min_moved

end ComponentExpansion

end ABC3.Found.PGC
