import ABC3.Found.PGC.FirstJumpIsBottom

/-!
# [pGC] 跳びが上に有界なら総和は閉じる ＋ 木の 2 つの docstring の突き合わせ

## 持ち場（前波で「残っている本当の 1 点」とした点）

「`H_K` の第 1 跳びが `k` に依らない」ことの証明。

## ★★木の 2 つの docstring を突き合わせた（食い違っていない）

| 場所 | 字面 | 読み |
|---|---|---|
| `RamificationJumpBound.lean:334` | 「`n = p` で `e = e_L = v_L(p)` を入れると `(p−1)·i ≤ e_L`」 | ★**上の体**の `e_L` |
| `NormalizedTraceDescent.lean:224` | 「`j` = 下の層 `E₁/F` の跳び、`e` = `e_{E₁}`」、`hjump : (p−1)j ≤ e` | ★**下の層に当てた形** |

★一見ずれて見えるが、★**同じ定理を別の層に当てているだけ**である ——
層 `E₁/F` に当てれば「上の体」は `E₁` なので `(p−1)·j ≤ e_{E₁}` ✓。
★2 つの docstring は**整合している**（前波で私が疑った点はここで解消した）。

## ★★何が本当に残っているか

`RamificationJumpBound` が与えるのは ★**各層についての** `(p−1)·i ≤ e`（その層の上の体の `e`）。
私の鎖が要るのは ★**第 1 跳び `t k` が `k` に依らずに押さえられる**ことである。

* 円分塔では測定で `t = p`（★5 設定すべて。`tools/numerology-check.py` (n1)）。
* 一般には、第 1 跳びは ★**塔の底の層の性質**であって、上に伸ばしても変わらない
  （`H_K` は `k` によらず同じ群の像で、最小の跳びは下で決まる）。
  ★**この「変わらない」ことの証明は書いていない。** ★木の `RamificationJumpBound` にも
  無い（あちらは 1 層ごとの上界である）。

⇒ ★本ファイルは「`t k ≤ T` を**仮定すれば**総和が閉じる」ところまでを定理にした（§1）。
★残るのは `t k ≤ T` の供給 1 本である。

## 何を埋めたか

| 節 | 宣言 | 内容 |
|---|---|---|
| §1 | `sum_le_of_bounded_jump` | ★★跳びが上に有界なら総和は `(T+p−2)p/((p−1)²)` 以下 |
| §2 | `quadratic_break_dvd` / `quadratic_strict` | ★木の docstring の例を手で検算（6 例目） |

★`t k ≡ p` を入れると `(2p−2)·p/((p−1)²) = 2p/(p−1)` で
★`SharpPrizeBound.sharp_exponent_sum_le` と**同じ値**になる
（`FirstJumpIsBottom.two_p_sub_two_recovers`、前波）。

## ★木の docstring の検算（6 例目）

`RamificationJumpBound.lean:340` の

> (`ℚ₂(√2)/ℚ₂` は `p = 2`, `i = 2` で確かに `p ∣ i`)

を**手で確かめた**: `π = √2`、`σπ = −√2` ⇒ `σπ − π = −2√2` ⇒
`v_L(−2√2) = v_L(2) + v_L(√2) = 2 + 1 = 3` ⇒ `i(σ) = 3`、跳び `i = i(σ) − 1 = 2` ★**一致**。
`2 ∣ 2` ✓。★また同じ例で `(p−1)i ≤ p·e` は `1·2 < 2·2` と**狭義**である（§2）。
★docstring は等号が起きるとは言っていないので、★**断定は正しい**。

## 逸脱の記録

- ★台帳の訂正: 前波のファイル名は **`FirstJumpIsBottom.lean`** である
  （本体は `FirstJumpCorrection.lean` と書いたが、そのファイルは**存在しない**。`ls` で確認）。
  ★本日 3 度目の場所の食い違いなので、★**名指しは必ず `ls` で確かめてから**が要る。
- §2 の 2 本は手計算の記録であって、`ℚ₂(√2)` の分岐を Lean で構成したものではない。
-/

namespace ABC3.Found.PGC

namespace BoundedJumpSum

open Finset

/-! ## §1 ★跳びが**上に有界**なら総和は閉じる -/

section Sum

/-- ★★★前波で「残っている本当の 1 点」とした形。

第 1 跳び `t k` が `k` に依らない `T` で押さえられていれば、
1 段の損失 `(t k + p − 2)/(p^k(p−1))` の総和は閉じる。

★`t k ≡ p`（円分塔、測定 5 設定）なら `T = p` で
`(2p−2)·p/((p−1)²) = 2p/(p−1)` ——★`SharpPrizeBound` と同じ値になる。 -/
theorem sum_le_of_bounded_jump {p T : ℝ} (hp : 2 ≤ p) (t : ℕ → ℝ)
    (ht0 : ∀ k, 0 ≤ t k) (htT : ∀ k, t k ≤ T) (n : ℕ) :
    ∑ k ∈ range n, (t k + (p - 2)) / (p ^ k * (p - 1))
      ≤ (T + (p - 2)) * p / ((p - 1) * (p - 1)) := by
  have hp0 : (0 : ℝ) < p := by linarith
  have hpm : (0 : ℝ) < p - 1 := by linarith
  have hc : (0 : ℝ) ≤ T + (p - 2) := by
    have := le_trans (ht0 0) (htT 0)
    linarith
  refine le_trans (Finset.sum_le_sum ?_)
    (FirstJumpIsBottom.sum_le_of_fixed_constant hp hc n)
  intro k _
  have hden : (0 : ℝ) < p ^ k * (p - 1) := by positivity
  gcongr
  exact htT k

end Sum

/-! ## §2 ★木の docstring を検算した（6 例目） -/

section Check

/-- ★`RamificationJumpBound.lean:340` の docstring は

> ★等号 `(p−1)i = p·e` は `p ∣ (p−1)i` を強いるが `gcd(p, p−1) = 1` なので `p ∣ i`。
> (`ℚ₂(√2)/ℚ₂` は `p = 2`, `i = 2` で確かに `p ∣ i`)

と書いている。★**手で確かめた**: `π = √2`、`σπ = −√2` なので
`σπ − π = −2√2`、`v_L(−2√2) = v_L(2) + v_L(√2) = 2 + 1 = 3`。
⇒ `i(σ) = 3`、跳びは `i = i(σ) − 1 = 2` ★**一致する**。そして `2 ∣ 2` ✓。 -/
theorem quadratic_break_dvd : (2 : ℕ) ∣ 2 := ⟨1, rfl⟩

/-- ★同じ例で `(p−1)·i ≤ p·e` は `1·2 ≤ 2·2` で**狭義**（等号ではない）。
★docstring は等号が起きるとは言っていない（`p ∣ i` の例として挙げている）。 -/
theorem quadratic_strict : (2 - 1) * 2 < 2 * 2 := by norm_num

end Check

/-! ## §3 使っている公理の一覧 -/

#print axioms sum_le_of_bounded_jump
#print axioms quadratic_break_dvd
#print axioms quadratic_strict

end BoundedJumpSum

end ABC3.Found.PGC
