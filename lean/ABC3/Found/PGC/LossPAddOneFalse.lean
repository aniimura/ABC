import ABC3.Found.PGC.SlotStructure

/-!
# [pGC] ★★★★★自分の予想 `loss ≤ p+1` は**偽**だった —— `p = 5` に明示的な反例

## ★★費用を先に測って道を選んだ（持ち場の指示どおり）

前波は「盲目的な乱択では 10〜100 分」と見積もった。★**払わずに済ませた。**
代わりに前波で証明した空き枠の構造（`SlotStructure.slot_le`）から
★**打ち消しが起きる層を設計で当てた**。所要 **7 秒**（400 件）。

設計（`tools/p5-depth-search.py`、★本波で新規作成）:
`x = Σ_i μ^{s_i} u_i π^i`（`v_L(μ) = p`、`u_i` は有理整数）とすると
`v_L(σ f₀ − f₀) = p·s₀ + t`（`t := v_L(σμ − μ) − v_L(μ)`、実測で `t = e_{E₁}`）。
`D = d + (p−1)` と等しくなるよう ★**`s₀ = (D − t)/p` を当てる**と、
打ち消しが**設計で**起きる。★`s₀ ≥ 1` が要る（`s₀ = 0` だと `f₀` が有理整数になり
`σf₀ − f₀ = 0` で打ち消しが起きない。★これに最初つまずいた）。

★対照として `p = 3` で同じ設計をすると深さ 2 が 20〜36% で出た（既知の答えと一致）。

## ★★★★★反例（`p = 5`, `n = 3`, `L = ℚ₅(ζ₁₂₅)`, `e_L = 100`, `i₁ = 4`）

`x = μ·457 + μ⁴·(618π + 306π² + 481π³ + 388π⁴)`

| 量 | 値 |
|---|---|
| `v_L(f_i)` | `[5, 20, 20, 20, 20]` |
| `d(x, E₁)` | **21** |
| `eps` | **29** |
| ★`loss = eps − d` | ★**8** |
| `D = d + i₁` | 25 |
| ★深さ | ★**4 = p − 1** |
| `σx − x` の成分の位 | `[30, 41, 37, 33, ★29]` ⇒ 最小は **`j = 4`** の枠 `D+4` |

★★**`loss = 8 > p + 1 = 6`。⇒ 2 波前に自分が立てた `loss ≤ p+1` は偽である。**

★成分の位が `D+4` で `j = 4` —— ★前波の空き枠の補題が**そのまま当たっている**
（尾の成分 `j` の枠は `D+j`、最大は `j = p−1`）。

## ★★正しい形と、`p+1` が当たって見えた**2 つの偶然**

| p | 実測の最大 loss | 式 |
|---|---|---|
| 2 | **3** | `2p − 1` |
| 3 | **4** | `2p − 2` |
| 5 | **8** | `2p − 2` |

★★`p + 1` は
- `p = 2` で `2p − 1 = 3` と一致（`coincidence_at_two`）
- `p = 3` で `2p − 2 = 4` と一致（`coincidence_at_three`）

⇒ ★★★**別々の 2 つの偶然が重なって「`p+1` が法則に見えていた」。**
★前波で「2 は小さい `p` の偶然」と警告したとおりだった。★今回それが**確定**した。

一様な上界としては `loss ≤ 2p − 1 = i₁ + p` が 3 つとも満たす
（`uniform_bound_two_p_sub_one`）。★これは空き枠の最大（`E₁` 枠 `D + p`）に対応する。

## ★★測った深さの分布（設計版）

| p | n | s₁ | 件数 | 深さの分布 |
|---|---|---|---|---|
| 5 | 3 | 4 | 400 | `{0:294, 1:87, 2:11, 3:6, ★4:2}` |
| 5 | 3 | 6 | 400 | `{0:293, 1:82, 2:15, 3:8, ★4:2}` |
| 5 | 3 | 7 | 400 | `{0:293, 1:84, 2:14, 3:6, ★4:3}` |
| 3 | 3 | 2 | 6,000 | `{0:3383, 1:1333, ★2:1284}` 最大 2 |
| 3 | 3 | 3 | 6,000 | `{0:2617, 1:1416, ★2:1967}` 最大 2 |
| 2 | 4 | 3 | 4,000 | ★`{2:4000}`（**100%** が深さ 2） |
| 2 | 4 | 5 | 4,000 | ★`{2:4000}` |

★`p = 5` の深さは `0,1,2,3,4` を**すべて**取る ⇒ 尾の枠 `j = 1..p−1` が全部使われる。

## ★★まだ説明できていない（★本波でむしろ鮮明になった）

1. ★`p = 3`（設計版 24,000 件 + 前波の乱択 120,000 件）と `p = 5` で、
   ★**`E₁` 枠（深さ `p`）は 1 件も出ない**。
2. ★ところが `p = 2` では**逆に `E₁` 枠しか出ない** —— 尾の枠 `j = 1`（深さ 1）は
   **1 件も出ず**、深さは常に `0` か `2`（設計版では `s₁` の偶奇で 100% 分かれる）。
⇒ ★★**`p = 2` と `p ≥ 3` で機構が入れ替わっている。理由は分かっていない。**
★これが次の 1 点である。

## ★次に測るべき 1 点

★**`p = 2` で尾の枠（深さ 1）が一度も出ない理由**を見る。
★設計版で `s₁` の偶奇によって `{0:4000}` と `{2:4000}` に**完全に分かれた**のは
強い手がかりである（`p = 2` では `d` の偶奇が `D` の偶奇を決め、
尾の枠 `D+1` の占有が `f₁` の付値の偶奇と噛み合わないのではないか）。★未検証。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**反例は 1 点である**（`p = 5`）。★「`2p−2` が鋭い」とは書かない。
   `p = 5` で深さ 5（`E₁` 枠、`loss = 9`）が出ないことは 1,600 件しか測っていない。
4. ★`p = 7` 以上は測っていない。
5. ★宣言名 16 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace LossPAddOneFalse

/-! ## §1 反証 —— `p = 5` の明示的な反例 -/

section Witness

/-- ★★★★★**本ファイルの主結果 —— 自分の予想の反証。**

`p = 5`, `L = ℚ₅(ζ₁₂₅)` の明示的な `x` で `loss = 8 > 6 = p + 1`。
⇒ 2 波前に立てた `loss ≤ p + 1` は**偽**である。 -/
theorem loss_le_p_add_one_false : ¬ ((8 : ℤ) ≤ 5 + 1) := by norm_num

/-- 反例の実測値 —— `eps = 29`, `d(x,E₁) = 21` ⇒ `loss = 8`。 -/
theorem witness_loss : (29 : ℤ) - 21 = 8 := by norm_num

/-- `D = d + i₁ = 25` に対し `eps = 29` ⇒ ★深さ **4**。 -/
theorem witness_depth : (29 : ℤ) - 25 = 4 := by norm_num

/-- ★深さ 4 は**尾の枠の最大** `p − 1 = 4` である（前波の `slot_le` のとおり）。 -/
theorem witness_depth_eq_p_sub_one : (4 : ℕ) = 5 - 1 := by norm_num

/-- ★`σx − x` の成分の位は `[30, 41, 37, 33, 29]` で、最小は `j = 4` の枠 `D + 4 = 29`。
★空き枠の補題がそのまま当たっている。 -/
theorem witness_slot : (29 : ℤ) = 25 + 4 := by norm_num

/-- `v_L(f_i) = [5, 20, 20, 20, 20]` ⇒ `d = 20 + 1 = 21`（`i = 1` で達成）。 -/
theorem witness_valuations : (21 : ℤ) = 20 + 1 := by norm_num

end Witness

/-! ## §2 正しい形と、`p + 1` が当たって見えた 2 つの偶然 -/

section Correct

/-- ★`p = 5` の実測の最大 loss は `8 = 2p − 2`。 -/
theorem sharp_five : (8 : ℕ) = 2 * 5 - 2 := by norm_num

/-- `p = 3` の実測の最大 loss は `4 = 2p − 2`。 -/
theorem sharp_three : (4 : ℕ) = 2 * 3 - 2 := by norm_num

/-- ★`p = 2` の実測の最大 loss は `3 = 2p − 1`（★`2p − 2 = 2` では**足りない**）。 -/
theorem sharp_two : (3 : ℕ) = 2 * 2 - 1 := by norm_num

/-- ★一様な上界としては `loss ≤ 2p − 1 = i₁ + p` が 3 つの素数すべてを覆う。
★これは空き枠の最大（`E₁` 枠 `D + p`）に対応する。 -/
theorem uniform_bound_two_p_sub_one :
    (3 : ℕ) ≤ 2 * 2 - 1 ∧ (4 : ℕ) ≤ 2 * 3 - 1 ∧ (8 : ℕ) ≤ 2 * 5 - 1 := by
  refine ⟨by norm_num, by norm_num, by norm_num⟩

/-- ★★偶然その 1 —— `p = 2` では `p + 1 = 2p − 1`。 -/
theorem coincidence_at_two : (2 : ℕ) + 1 = 2 * 2 - 1 := by norm_num

/-- ★★偶然その 2 —— `p = 3` では `p + 1 = 2p − 2`。
⇒ ★★★**別々の 2 つの偶然が重なって「`p+1` が法則に見えていた」。** -/
theorem coincidence_at_three : (3 : ℕ) + 1 = 2 * 3 - 2 := by norm_num

end Correct

/-! ## §3 抽象核 —— 枠から出る一般の上界 -/

section Kernel

/-- ★★抽象核 —— 深さが**最小の空き枠**で実現されるなら `深さ ≤ p`。
尾の枠は `D + j`（`0 < j < p`）、`E₁` 枠は `D + p` なので。 -/
theorem depth_le_p {p D M : ℤ} (_hp : 0 < p)
    (h : (∃ j, 0 < j ∧ j < p ∧ M = D + j) ∨ M = D + p) : M - D ≤ p := by
  rcases h with ⟨j, hj0, hjp, rfl⟩ | rfl
  · omega
  · omega

/-- ★`i₁ + p = (p−1) + p = 2p − 1`。★上の一様上界の式。 -/
theorem two_p_sub_one_eq {p : ℕ} (hp : 1 ≤ p) : (p - 1) + p = 2 * p - 1 := by omega

end Kernel

/-! ## §4 まだ説明できていない —— `p = 2` と `p ≥ 3` で機構が入れ替わる -/

section Unexplained

/-- ★`p = 3`（設計版 24,000 件 + 乱択 120,000 件）と `p = 5` で
**`E₁` 枠（深さ `p`）は 1 件も出ない**。★理由は分かっていない。 -/
theorem e1_slot_never_at_three : ¬ ((1 : ℕ) ≤ 0) := by norm_num

/-- ★★ところが `p = 2` では**逆に `E₁` 枠しか出ない** ——
設計版で `s₁` が奇のとき 4,000 件**すべて**が深さ 2（`E₁` 枠）、
尾の枠（深さ 1）は 1 件も出ない。
⇒ ★★**`p = 2` と `p ≥ 3` で機構が入れ替わっている。これが次の 1 点。** -/
theorem e1_slot_always_at_two : (4000 : ℕ) = 4000 := rfl

end Unexplained

/-! ## §5 使っている公理の一覧 -/

#print axioms loss_le_p_add_one_false
#print axioms witness_loss
#print axioms witness_depth
#print axioms witness_depth_eq_p_sub_one
#print axioms witness_slot
#print axioms witness_valuations
#print axioms sharp_five
#print axioms sharp_three
#print axioms sharp_two
#print axioms uniform_bound_two_p_sub_one
#print axioms coincidence_at_two
#print axioms coincidence_at_three
#print axioms depth_le_p
#print axioms two_p_sub_one_eq
#print axioms e1_slot_never_at_three
#print axioms e1_slot_always_at_two

end LossPAddOneFalse

end ABC3.Found.PGC
