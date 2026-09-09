import ABC3.Found.PGC.WitnessFiveRecursion

/-!
# [pGC] ★★(N) を測った —— 柱は `c` を**改善しない**。★ただし**賞金の額**は確定した

## ★(N) の問い（前波で私が「未測定」と書いた点）

> 柱が仮定ゼロになったことが `AxWildDescent` に何を与えるか。

★開いて測った。**答えは「今の形では `c` を改善しない」**である。以下、測った内容。

## ★★測定 1 —— 木の現在地

`WildDepthFieldDescent.lean:141 axWildDescent_prime` が
★**無条件に `AxWildDescent K (fun _ => p)`** を出している。同 docstring:

> ★それでも `∏_{k ∈ Icc 1 n} p = p^n` は非有界なので `AxLemma` は出ない。
> ★残っているのは**分岐による絞り込み**だけである。

⇒ ★**私の柱はまさにその「分岐による絞り込み」に居る。**目標は明確である。

## ★★★測定 2 —— しかし柱は `c` の**上界**を与えない

`AxTowerDecay.lean:473 AxWildDescent K c` は
★**`K.closure` の任意の `x`** について降下 `x′` の存在を要求する。
私の柱が言うのは

- `ρ = σπ − π = (1+π)·w` は `𝒪_{E₁}` 成分が **2 つだけ**（`R_0 = R_1 = w`, `v_L(w) = p`）
- `E₁` 上の跳びは `t = p(p−1)`

であって、★**どちらも「与えられた `x` に対する `x′` の構成」ではない。**
★上界に要る「最小の空き枠が必ず占められている」ことは、
★★**16 波かけて `SlotStructure` で枠を数え、`TailNoCancel` で原因を突き止め、
それでもなお measured のまま**である（`p=2,3` は 240,000 件、`p=5` は 2,500 件）。

⇒ ★★**柱は機構を説明したが、上界は与えない。**（★「今日はこれ以上進めない」の中身。）

## ★★★★測定 3 —— それでも**賞金の額**が確定した（本ファイルの成果）

実測の鋭い 1 段の定数は `p^{(2p−2)/e_k}`（`e_k = p^{k−1}(p−1)`）だった
（`p=3,5` で `2p−2`、`p=2` は `2p−1`）。★その**指数の総和**を計算すると:

`Σ_{k≥1} (2p−2)/(p^{k−1}(p−1)) = 2·Σ_{k≥0} p^{−k} = 2p/(p−1)`（`sharp_exponent_sum_le`）

| p | 総和 | ⇒ `AxLemma` の定数 `C` |
|---|---|---|
| 2 | **4** | `2^4 = 16` |
| 3 | **3** | `3^3 = 27` |
| 一般 | `2p/(p−1)` | `p^{2p/(p−1)}` |

⇒ ★★★**もし実測の 1 段定数が定理になれば、`AxLemma K (p^{2p/(p−1)})` が出る。**
`axLemma_of_wildDescent` は「`c` の**任意の**有限積が `C` 以下」を要求するが、
★総和が有限なので**それは満たされる**。⇒ **`AxSenTate` まで届く。**

★比較: `axConstant p = p^{p/(p−1)²}`（`p=3` で `3^{3/4}`）。
★実測の定数はそれよりずっと大きいが、★**有限であることだけが問題なので十分**である。
（`AxLemma` は定数の大きさを問わない。）

## ★★残る 1 点は**1 つに絞られた**

> ★**「1 段の損失は `p^{(2p−2)/e_k}` を超えない」を証明すること。**

★これが証明できれば、本ファイルの `sharp_exponent_sum_le` が総和を閉じ、
`axLemma_of_wildDescent` が `AxLemma` を出し、`AxSenTate` が出る。
★★**`cor_3_1` の入口の残り 1 段が、この 1 つの不等式に落ちた。**

★現状: `p=2,3` は 240,000 件、`p=5` は 2,500 件で**破れ 0**。機構も
（`f₀` が原因・空き枠は `D+j`・`p=2` は退化）**証明済み**。
★足りないのは「最小の枠が必ず占められている」ことだけである。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**本ファイルは「実測の定数が定理になれば」という条件つきの計算である。**
   ★1 段の上界は**証明していない**。★「`AxLemma` が出た」とは書かない。
4. ★`p = 2` の実測は `2p−1`（`2p−2` ではない）だが、
   `2p−1 = 3 ≤ 4 = 2p−2 + 1` で総和は `Σ (2p−1)/e_k = (2p−1)/(p−1)·... ` と少し増える。
   ★本ファイルは `2p−2` の側だけを計算した。★`p=2` は別途 `2^{(2·2−1)/(2−1)·2} ` 相当。
   ★**そこは計算していない**（`p ≥ 3` で `2p−2` が実測、という母集団を添える）。
5. ★宣言名 4 件は **実ファイル**で衝突 0 を確かめた（#348）。
6. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace SharpPrizeBound

/-! ## §1 抽象核 —— 部分等比和の上界 -/

section Geom

/-- ★核 —— 部分等比和は `1/(1−r)` を超えない（`0 ≤ r < 1`）。 -/
theorem partial_geom_le {r : ℝ} (h0 : 0 ≤ r) (h1 : r < 1) (n : ℕ) :
    ∑ k ∈ Finset.range n, r ^ k ≤ 1 / (1 - r) := by
  have hne : r ≠ 1 := ne_of_lt h1
  have hpos : (0 : ℝ) < 1 - r := by linarith
  rw [geom_sum_eq hne]
  have hrw : (r ^ n - 1) / (r - 1) = (1 - r ^ n) / (1 - r) := by
    rw [div_eq_div_iff (by linarith) (by linarith)]
    ring
  rw [hrw]
  have hrn : (0 : ℝ) ≤ r ^ n := pow_nonneg h0 n
  gcongr
  linarith

end Geom

/-! ## §2 測った鋭い定数の指数の総和 -/

section Prize

/-- ★★★**賞金の額** —— 実測の鋭い 1 段定数 `p^{(2p−2)/e_k}` の**指数の総和**は
`2p/(p−1)` を超えない（`e_k = p^{k−1}(p−1)`）。

⇒ ★もし 1 段の上界が定理になれば、`axLemma_of_wildDescent` の
「任意の有限積が `C` 以下」が満たされ **`AxLemma K (p^{2p/(p−1)})`** が出る。
★★`cor_3_1` の入口の残り 1 段が、**この 1 つの不等式**に落ちた。 -/
theorem sharp_exponent_sum_le {p : ℝ} (hp : 2 ≤ p) (n : ℕ) :
    ∑ k ∈ Finset.range n, (2 * p - 2) / (p ^ k * (p - 1)) ≤ 2 * p / (p - 1) := by
  have hp0 : (0 : ℝ) < p := by linarith
  have hp1 : (0 : ℝ) < p - 1 := by linarith
  have hrw : ∀ k : ℕ, (2 * p - 2) / (p ^ k * (p - 1)) = 2 * (1 / p) ^ k := by
    intro k
    have hpk : (p : ℝ) ^ k ≠ 0 := (pow_pos hp0 k).ne'
    rw [div_pow, one_pow]
    field_simp
  rw [Finset.sum_congr rfl (fun k _ => hrw k), ← Finset.mul_sum]
  have hgeom := partial_geom_le (r := 1 / p) (by positivity)
    (by rw [div_lt_one hp0]; linarith) n
  have hone : (1 : ℝ) - 1 / p = (p - 1) / p := by field_simp
  rw [hone] at hgeom
  have : (1 : ℝ) / ((p - 1) / p) = p / (p - 1) := by field_simp
  rw [this] at hgeom
  have h2 : 2 * p / (p - 1) = 2 * (p / (p - 1)) := by ring
  rw [h2]
  linarith

/-- `p = 3`: 総和 `3` ⇒ `C = 3^3 = 27`。★`axConstant 3 = 3^{3/4}` より大きいが、
★**有限であることだけが問題**なので十分。 -/
theorem prize_exponent_three : (2 * (3 : ℝ)) / (3 - 1) = 3 := by norm_num

/-- `p = 2`: 総和 `4` ⇒ `C = 2^4 = 16`。 -/
theorem prize_exponent_two : (2 * (2 : ℝ)) / (2 - 1) = 4 := by norm_num

end Prize

/-! ## §3 使っている公理の一覧 -/

#print axioms partial_geom_le
#print axioms sharp_exponent_sum_le
#print axioms prize_exponent_three
#print axioms prize_exponent_two

end SharpPrizeBound

end ABC3.Found.PGC
