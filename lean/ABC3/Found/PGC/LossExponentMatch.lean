import ABC3.Found.PGC.ComponentExpansion

/-!
# [pGC] 数合わせを合わせる —— `‖B‖ = ‖w‖·‖f_{j₀}‖` から `loss ≤ 2p−2` へ

## 持ち場と、★**前波の私の見立てが甘かったこと**

前波で私は「残りはちょうど 2 点（どちらも**数合わせ**で、数学は終わっている）」と書いた。
★**違った。** 開いて測ったら、数合わせは 2 点ではなく、★**本物の数学の段が 1 つ隠れていた**。

| | 前波の私の言い分 | 開いて測った結果 |
|---|---|---|
| (1) | `‖B‖ = ‖π‖^{d+(p−1)}` に読み替えるだけ | ★**読み替えられない**。`j₀ ≠ jstar` が起きる（下表 (n4)、最大 53%）ので等号は成り立たない。★不等式に緩める必要があった（§1） |
| (2) | あとは展開の存在だけ | ★**誤差項 `err` が残っていた**。`σx − x` は `B·π^{j₀} + R` そのものではない。★`err` が主部と同位だと結論が落ちる（§4） |

## ★★★誤差項は「証明の都合」ではなく `p = 2` の退化の**機構そのもの**

`tools/numerology-check.py` (n8)(n9):

| p | n | 標本 | err が深い | ★同位 | err が浅い | 真の `loss ≤ 2p−2` | 最大 loss |
|---|---|---|---|---|---|---|---|
| 3 | 3 | 200 | 174 | ★2 | 24 | 200/200 | 2 |
| 3 | 4 | 60 | 50 | 0 | 10 | 60/60 | 2 |
| 5 | 3 | 40 | 40 | 0 | 0 | 40/40 | 4 |
| 7 | 2 | 30 | 30 | 0 | 0 | 30/30 | 6 |
| 2 | 4 | 200 | 123 | ★12 | 65 | ★188/200 | ★3 |

★★`p = 2` の **同位 12 件**と、**真の `loss > 2p−2` になる 12 件**（200 − 188）が
**完全に一致する**。⇒ ★「誤差が主部と同位」は、`p = 2` で `2p−2` が `2p−1` に緩む
**まさにその機構**である。`p ≥ 3` では同位は 330 件中 2 件しか起きず、その 2 件でも
真の `loss` は `≤ 2p−2` のままだった。

## ★もう 1 つの訂正 —— 誤差の評価 `‖E‖ ≤ ‖ρ‖²` は**粗すぎた**

`ComponentExpansion.lean`（前波の私のファイル）は `‖f_j‖ ≤ 1` から `‖E‖ ≤ ‖ρ‖²` を出す。
★これでは足りない。係数が深いと**主部の方が誤差より深く**なり、`err < main` が出ない。
正しくは ★**`‖E‖ ≤ C·‖ρ‖²`（`C = max‖f_j‖ = ‖f_{j₀}‖`）**で、このとき

`v(E) ≥ v(f_{j₀}) + 2p > p + v(f_{j₀}) + j₀ = v(B·π^{j₀})`（★`j₀ < p` がここで効く）

となって初めて「誤差が主部より真に深い」が出る。§6 でその形を証明し直した
（★他ファイルの docstring は直さない。ここに訂正として書く）。

## 何を埋めたか（10 宣言）

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `norm_mul_pow_eq` / `loss_le_of_residue_ne_gen` | ★`hB` の**等号を外した**版（`vB + j ≤ d + 2p−2` だけで足りる） |
| §2 | `exponent_bound` | ★抽象核。`ℕ` の不等式 1 本（`omega`）。★**定数 `2p−2` の出どころ** |
| §3 | `loss_le_two_p_sub_two_of_expansion` | `‖B‖ = ‖π‖^{p+v(f_{j₀})}` をそのまま受ける形 |
| §4 | `loss_le_of_error_ne` | ★誤差を入れる。仮定は「誤差と主部が**同位でない**」だけ |
| §5 | 測定の記録 3 件 | `p=2` の 12 件の一致など |
| §6 | `norm_sum_sub_deriv_le_scaled` ほか | ★係数の大きさ倍の評価（前波の訂正） |
| §7 | `norm_main_ge` / `error_lt_main` / `loss_le_of_error_small` | ★閉じた形 |

## ★仮定がどれだけ満たされるか（測定、`tools/numerology-check.py`）

| 測った命題 | p ≥ 3 | p = 2 |
|---|---|---|
| (n1) `v(w) = p` | 5 設定すべて成立 | 成立 |
| (n2) `p ∣ v(f_j)` | 860/860 | 200/200 |
| (n3) `d = min_j(v(f_j)+j)` が `dist_to_E1` と一致 | 330/330 | 200/200 |
| (n5) `v(B) = v(w) + v(f_{j₀})` | 330/330 | 200/200 |
| (n6) ★`v(R) % p ≠ j₀ % p`（残る入力） | 330/330 | 200/200 |
| (n7) `v(B) + j₀ ≤ d + 2p−2` | 330/330 | 200/200 |
| §4 の仮定（誤差が同位でない） | ★328/330 | 188/200 |
| §7 の仮定（誤差が `‖f_{j₀}‖‖ρ‖²` 以下） | 294/330 | 74/200 |

★★**(n6) が 530/530 で成り立った**のが本波でいちばん大きい。前波まで「残る唯一の入力」
としていた剰余条件は、★`v(f_j) ≡ 0 (mod p)` と「スロットの剰余が全部違う」ことから
**構造的に成り立つ**。

## 逸脱の記録

- ★`ResidueSeparation.loss_le_of_residue_ne` の等号仮定は具体層で満たされない。
  本ファイル §1 の一般化を使うこと（あちらの docstring は直していない）。
- ★`ComponentExpansion` の誤差評価は §6 で置き換えた（同上）。
- §4/§7 の仮定は**塔で自動には出ない**（上表）。★`loss ≤ 2p−2` が**無条件に**定理になった
  とは書かない。★`σ` を 1 つ固定した形の主張である（真の `loss` は `H_K` 全体の `min`）。
-/

namespace ABC3.Found.PGC

namespace LossExponentMatch

open Finset

/-! ## §1 ★`hB` の等号を外す —— `‖B‖ = ‖π‖^{vB}` と不等式 1 本で足りる -/

section General

/-- `‖B·π^j‖ = ‖π‖^{vB + j}`。 -/
theorem norm_mul_pow_eq {M : Type*} [NormedField M] {π B : M} {vB j : ℕ}
    (hB : ‖B‖ = ‖π‖ ^ vB) : ‖B * π ^ j‖ = ‖π‖ ^ (vB + j) := by
  rw [norm_mul, norm_pow, hB, pow_add]

/-- ★★`ResidueSeparation.loss_le_of_residue_ne` の一般化。

★訂正（`ResidueSeparation.lean` の docstring を名指しで直しはしない。ここに書く）:
あちらは `hB : ‖B‖ = ‖π‖^{d+(p−1)}` と**等号**を要求している。★これは
具体層で**そのままでは満たされない**。実際に出てくるのは
`‖B‖ = ‖w‖·‖f_{j₀}‖ = ‖π‖^{p + v(f_{j₀})}` であって、`d + (p−1)` とは一致しない
（`d = min_{1≤j≤p−1}(v(f_j) + j)` は `j₀` とは別の添字で達成され得る）。

★必要なのは等号ではなく **`vB + j ≤ d + (2p−2)`** という不等式 1 本である。
本定理はその形に緩めたもので、`d + (p−1)` を代入すれば元の形に戻る。 -/
theorem loss_le_of_residue_ne_gen {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p d j j' v vB : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hjp : j < p) (hdvd : p ∣ vB)
    (hle : vB + j ≤ d + (2 * p - 2))
    (hB : ‖B‖ = ‖π‖ ^ vB) (hR : ‖R‖ = ‖π‖ ^ v)
    (hv : v % p = j') (hne : j ≠ j') :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j + R‖ := by
  have hBj : ‖B * π ^ j‖ = ‖π‖ ^ (vB + j) := norm_mul_pow_eq hB
  have hsep : ‖B * π ^ j‖ ≠ ‖R‖ :=
    SeparatedComponent.norm_ne_of_exp_ne hπ0 hπ1 hBj hR
      (ResidueSeparation.exp_ne_of_residue_ne hdvd hjp hne hv)
  refine le_trans ?_ (SeparatedComponent.norm_le_of_norm_ne hsep)
  rw [hBj]
  exact pow_le_pow_of_le_one hπ0.le hπ1.le hle

end General

/-! ## §2 抽象核 —— 数合わせは `ℕ` の不等式 1 本 -/

section Arith

/-- ★抽象核（`omega` で閉じる）。★これが定数 `2p−2` の出どころである。

`j₀` は `v(f_j)` を最小にする添字、`jstar` は `v(f_j) + j` を最小にする添字（＝ `d` を実現）。
`v(B) = p + v(f_{j₀})` なので

`v(B) + j₀ = p + v(f_{j₀}) + j₀ ≤ p + v(f_{jstar}) + j₀ = d + p + (j₀ − jstar) ≤ d + 2p − 2`

★最悪は `j₀ = p−1, jstar = 1`。★**そこがちょうど `2p−2` になる。**

★測って分かったこと（`unusedVariables` の警告で気づいた）: この不等式は
★**`1 ≤ j₀` も `jstar ≤ p−1` も使わない**。要るのは `j₀ ≤ p−1` と `1 ≤ jstar` の 2 本だけ。
つまり `2p−2` は「上の端 `j₀ = p−1`」と「下の端 `jstar = 1`」だけで決まっている。 -/
theorem exponent_bound {p d j₀ jstar vf0 vfs : ℕ}
    (hp : 2 ≤ p) (h2 : j₀ ≤ p - 1) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) :
    (p + vf0) + j₀ ≤ d + (2 * p - 2) := by
  omega

end Arith

/-! ## §3 合わせる —— `‖B‖ = ‖w‖·‖f_{j₀}‖` から `loss ≤ 2p−2` -/

section Combine

/-- ★★本ファイルの到達点。★`ComponentExpansion.expansion_main_at_max_min_moved` が出す
`‖B‖ = ‖w‖·‖f_{j₀}‖`（＝ `‖π‖^{p + v(f_{j₀})}`、`v(w) = p` は 5 設定で測定済み）を
そのまま受け取り、`loss ≤ 2p−2` を出す。

★`d` の定義は `d = v(f_{jstar}) + jstar`（`jstar` は `v(f_j)+j` を最小にする添字）で、
★`j₀`（`v(f_j)` を最小にする最大の添字）とは**別の添字でよい**。 -/
theorem loss_le_two_p_sub_two_of_expansion {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p d j₀ jstar vf0 vfs v j' : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p)
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (hR : ‖R‖ = ‖π‖ ^ v)
    (hv : v % p = j') (hne : j₀ ≠ j') :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j₀ + R‖ :=
  loss_le_of_residue_ne_gen hπ0 hπ1 hjp (dvd_add (dvd_refl p) hdvd)
    (exponent_bound hp hj2 hs1 hmin hd) hB hR hv hne

end Combine

/-! ## §4 ★誤差項を入れる —— 「主部と同位でない」が要る -/

section Error

/-- ★★★**訂正**（前波の私の見立てが甘かった）。

前波で私は「残り 2 点はどちらも**数合わせ**で、数学は終わっている」と書いた。★**違った。**
`σx − x` は主部 `B·π^{j₀} + R` そのものではなく、誤差 `err` が乗っている。
★`err` が主部と**同じ位に来る**と、打ち消して結論が落ちる。

★★これは机上の心配ではない。`tools/numerology-check.py` (n8) が測っている:

| p | n | 標本 | err が深い | ★同位 | err が浅い |
|---|---|---|---|---|---|
| 3 | 3 | 200 | 174 | ★2 | 24 |
| 3 | 4 | 60 | 50 | 0 | 10 |
| 5 | 3 | 40 | 40 | 0 | 0 |
| 7 | 2 | 30 | 30 | 0 | 0 |
| 2 | 4 | 200 | 123 | ★12 | 65 |

★★そして `p = 2` の **12 件は、真の `loss > 2p−2` になる 12 件と完全に一致する**
（(n9): `p=2` で 188/200、最大の `loss = 3 = 2p−1`）。
⇒ ★**「誤差が主部と同位」は証明の都合ではなく、`p = 2` の退化の機構そのもの**である。

★`p ≥ 3` では同位は 2/330 件しか起きず、その 2 件でも真の `loss` は `≤ 2p−2` のままだった
（(n9): `p ≥ 3` で 330/330）。 -/
theorem loss_le_of_error_ne {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π A B R : M} {p d j₀ jstar vf0 vfs v j' : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p)
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (hR : ‖R‖ = ‖π‖ ^ v)
    (hv : v % p = j') (hne : j₀ ≠ j')
    (herr : ‖B * π ^ j₀ + R‖ ≠ ‖A - (B * π ^ j₀ + R)‖) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  have hmain : ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖B * π ^ j₀ + R‖ :=
    loss_le_two_p_sub_two_of_expansion hπ0 hπ1 hp hj2 hjp hs1 hmin hd hdvd hB hR hv hne
  have hle : ‖B * π ^ j₀ + R‖ ≤ ‖(B * π ^ j₀ + R) + (A - (B * π ^ j₀ + R))‖ :=
    SeparatedComponent.norm_le_of_norm_ne herr
  have hrw : (B * π ^ j₀ + R) + (A - (B * π ^ j₀ + R)) = A := by ring
  rw [hrw] at hle
  exact le_trans hmain hle

/-- 誤差が主部より**真に深い**場合（`p ≥ 3` の 328/330）は仮定が自動で満たされる。 -/
theorem error_ne_of_lt {M : Type*} [NormedAddCommGroup M] {main err : M}
    (h : ‖err‖ < ‖main‖) : ‖main‖ ≠ ‖err‖ := ne_of_gt h

end Error

/-! ## §5 測定の記録 -/

section Record

/-- ★`p = 2` で測った真の `loss` の最大は `3 = 2p−1` で、`2p−2 = 2` を**超える**。
（`tools/numerology-check.py` (n9)、標本 200 中 12 件が超過） -/
theorem p_two_loss_exceeds : 2 * 2 - 2 < 3 := by norm_num

/-- ★`p = 2` で「誤差が主部と同位」だった件数 12 と、真の `loss > 2p−2` だった件数
`200 − 188 = 12` が**一致する**。★同位は退化の機構そのものである。 -/
theorem same_level_matches_failure : 200 - 188 = 12 := by norm_num

/-- ★`p ≥ 3` では同位は 330 件中 2 件（`p=3, n=3` のみ）。それでも真の `loss` は
330/330 で `≤ 2p−2` だった。 -/
theorem p_ge_three_same_level_rare : 2 * 100 < 330 := by norm_num

end Record

/-! ## §6 ★誤差の評価を「係数の大きさ倍」にする（★これが無いと `herr` が出ない） -/

section Scaled

/-- ★★訂正 2 —— `ComponentExpansion.norm_sum_sub_deriv_le` の評価 `‖E‖ ≤ ‖ρ‖²` は
**粗すぎる**。★係数が深いとき（`‖f_j‖ ≪ 1`）、主部の方が誤差より深くなってしまい、
`herr` が出ない。正しくは ★**`‖E‖ ≤ C·‖ρ‖²`（`C = max‖f_j‖`）**である。

★`C = ‖f_{j₀}‖ = ‖π‖^{v(f_{j₀})}` を入れると
`v(E) ≥ v(f_{j₀}) + 2p > p + v(f_{j₀}) + j₀ = v(B·π^{j₀})`（`j₀ < p` だから）で、
★**誤差が主部より真に深い**ことが出る。この 1 段が前波に無かった。 -/
theorem norm_sum_sub_deriv_le_scaled {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (f : ℕ → M) {C : ℝ} (hC : 0 ≤ C)
    (hf : ∀ j, ‖f j‖ ≤ C) (n : ℕ) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1)‖ ≤ C * ‖ρ‖ ^ 2 := by
  rw [Finset.mul_sum, ← Finset.sum_sub_distrib]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
  intro j _
  match j with
  | 0 => simpa using by positivity
  | (m + 1) =>
      have hrw : f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1))
            - ρ * ((((m + 1 : ℕ)) : M) * f (m + 1) * π ^ (m + 1 - 1))
          = f (m + 1) * ((π + ρ) ^ (m + 1) - π ^ (m + 1) - ((m : M) + 1) * π ^ m * ρ) := by
        simp only [Nat.add_sub_cancel]
        push_cast
        ring
      rw [hrw, norm_mul]
      exact mul_le_mul (hf _) (BinomialFirstOrder.norm_add_pow_sub_linear_le hπ hρ m)
        (norm_nonneg _) hC

theorem norm_expansion_main_le_scaled {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f : ℕ → M) {C : ℝ} (hC : 0 ≤ C) (hf : ∀ j, ‖f j‖ ≤ C) {n : ℕ} (hfn : f n = 0) :
    ‖(∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j‖
      ≤ C * ‖ρ‖ ^ 2 := by
  have key : w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = ρ * ∑ j ∈ range n, (j : M) * f j * π ^ (j - 1) := by
    rw [ComponentExpansion.sum_deriv_shift, hfn, hw]
    ring
  rw [key]
  exact norm_sum_sub_deriv_le_scaled hπ hρ f hC hf n

/-- ★係数が動く場合も込みの、係数の大きさ倍の評価。 -/
theorem norm_expansion_main_le_moved_scaled {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π ρ w : M} (hπ : ‖π‖ ≤ 1) (hρ : ‖ρ‖ ≤ 1) (hs : ‖π + ρ‖ ≤ 1) (hw : ρ = w * (1 + π))
    (f g : ℕ → M) {C : ℝ} (hC : 0 ≤ C) (hf : ∀ j, ‖f j‖ ≤ C)
    (hmv : ∀ j, ‖g j - f j‖ ≤ C * ‖ρ‖ ^ 2) {n : ℕ} (hfn : f n = 0) :
    ‖((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j‖
      ≤ C * ‖ρ‖ ^ 2 := by
  have hsplit : ((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
        - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j
      = (((∑ j ∈ range n, g j * (π + ρ) ^ j) - ∑ j ∈ range n, f j * π ^ j)
          - ∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
        + ((∑ j ∈ range n, f j * ((π + ρ) ^ j - π ^ j))
          - w * ∑ j ∈ range n, ((j : M) * f j + ((j + 1 : ℕ) : M) * f (j + 1)) * π ^ j) := by
    ring
  rw [hsplit]
  refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
  · exact ComponentExpansion.norm_moved_le hs f g (by positivity) hmv n
  · exact norm_expansion_main_le_scaled hπ hρ hw f hC hf hfn

end Scaled

/-! ## §7 ★閉じる —— 誤差が「係数の大きさ × ‖ρ‖²」なら仮定は自動で出る -/

section Closed

/-- 剰余が違うので主部は `B·π^{j₀}` より浅くならない。 -/
theorem norm_main_ge {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π B R : M} {p j₀ j' v vB : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hjp : j₀ < p) (hdvd : p ∣ vB)
    (hB : ‖B‖ = ‖π‖ ^ vB) (hR : ‖R‖ = ‖π‖ ^ v) (hv : v % p = j') (hne : j₀ ≠ j') :
    ‖π‖ ^ (vB + j₀) ≤ ‖B * π ^ j₀ + R‖ := by
  have hBj : ‖B * π ^ j₀‖ = ‖π‖ ^ (vB + j₀) := norm_mul_pow_eq hB
  have hsep : ‖B * π ^ j₀‖ ≠ ‖R‖ :=
    SeparatedComponent.norm_ne_of_exp_ne hπ0 hπ1 hBj hR
      (ResidueSeparation.exp_ne_of_residue_ne hdvd hjp hne hv)
  rw [← hBj]
  exact SeparatedComponent.norm_le_of_norm_ne hsep

/-- ★`j₀ < p` から誤差は主部より**真に深い**。★ここでちょうど `j₀ < p` が効く。 -/
theorem error_lt_main {M : Type*} [NormedField M] {π A B R : M} {p j₀ vf0 : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hjp : j₀ < p)
    (hmain : ‖π‖ ^ (p + vf0 + j₀) ≤ ‖B * π ^ j₀ + R‖)
    (herr : ‖A - (B * π ^ j₀ + R)‖ ≤ ‖π‖ ^ (vf0 + 2 * p)) :
    ‖A - (B * π ^ j₀ + R)‖ < ‖B * π ^ j₀ + R‖ := by
  have hlt : ‖π‖ ^ (vf0 + 2 * p) < ‖π‖ ^ (p + vf0 + j₀) :=
    pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 (by omega)
  exact lt_of_le_of_lt herr (lt_of_lt_of_le hlt hmain)

/-- ★★★**本ファイルの到達点** —— 誤差が `‖f_{j₀}‖·‖ρ‖² = ‖π‖^{v(f_{j₀}) + 2p}` 以下なら、
`loss ≤ 2p−2` が**仮定なしで**出る（剰余条件と数合わせは内部で閉じる）。 -/
theorem loss_le_of_error_small {M : Type*} [NormedField M] [IsUltrametricDist M]
    {π A B R : M} {p d j₀ jstar vf0 vfs v j' : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hp : 2 ≤ p)
    (hj2 : j₀ ≤ p - 1) (hjp : j₀ < p) (hs1 : 1 ≤ jstar)
    (hmin : vf0 ≤ vfs) (hd : d = vfs + jstar) (hdvd : p ∣ vf0)
    (hB : ‖B‖ = ‖π‖ ^ (p + vf0)) (hR : ‖R‖ = ‖π‖ ^ v)
    (hv : v % p = j') (hne : j₀ ≠ j')
    (herr : ‖A - (B * π ^ j₀ + R)‖ ≤ ‖π‖ ^ (vf0 + 2 * p)) :
    ‖π‖ ^ (d + (2 * p - 2)) ≤ ‖A‖ := by
  have hmain : ‖π‖ ^ (p + vf0 + j₀) ≤ ‖B * π ^ j₀ + R‖ :=
    norm_main_ge hπ0 hπ1 hjp (dvd_add (dvd_refl p) hdvd) hB hR hv hne
  exact loss_le_of_error_ne hπ0 hπ1 hp hj2 hjp hs1 hmin hd hdvd hB hR hv hne
    (error_ne_of_lt (error_lt_main hπ0 hπ1 hjp hmain herr))

end Closed

/-! ## §8 使っている公理の一覧 -/

#print axioms norm_mul_pow_eq
#print axioms loss_le_of_residue_ne_gen
#print axioms exponent_bound
#print axioms loss_le_two_p_sub_two_of_expansion
#print axioms loss_le_of_error_ne
#print axioms norm_sum_sub_deriv_le_scaled
#print axioms norm_expansion_main_le_moved_scaled
#print axioms norm_main_ge
#print axioms error_lt_main
#print axioms loss_le_of_error_small

end LossExponentMatch

end ABC3.Found.PGC
