import ABC3.Found.PGC.TopIndexSurvives

/-!
# [pGC] ★★★残り 3 つのうち **(3) 二項展開の 1 次項が定理になった**

## ★どれから当たったか（★開いてから決めた。本日 7 波連続の形）

前波で残した 3 つ:

1. `𝒪_L = 𝒪_{E₁}[π]` の相対の冪基底と成分 `B_j` の抽出
2. `v(B_j) ∈ p·ℤ`
3. ★**(2) の主項の同定（二項展開の 1 次項）**

★持ち場が挙げた `JumpStrictMono.lean:84 norm_prod_one_add_sub_one_sub_sum_le` を**開いた**。
★**積の版**（`∏(1+η) − 1 − Σ η` の残りが 2 次）で、★**冪の版ではなかった**。
私が要るのは `(a+b)^i − a^i − i·a^{i−1}·b` の側である。⇒ ★**自分で建てた。**
（★同じ「1 次の項を取り出した残りは 2 次」の思想だが、形が違う。★在庫は当たったが使えなかった。）

## ★★★主結果

> `norm_add_pow_sub_linear_le` : 超距離ノルム体で `‖a‖ ≤ 1`, `‖b‖ ≤ 1` なら
> ★**`‖(a+b)^{m+1} − a^{m+1} − (m+1)·a^m·b‖ ≤ ‖b‖²`**

★これが前波の証明の (2) と (4) の両方を支える:

- **(2) 主項の同定**: `f_i·((π+ρ)^i − π^i)` の主項は `f_i·i·π^{i−1}·ρ` で、
  ★残りは `‖ρ‖²` の位（`v ≥ 2p`）に落ちる。
- **(4) 深い項が邪魔しないこと**: `ρ²` からの寄与が `v ≥ 2p > p` であること。
  ★★**本定理がその「`ρ²` 以降」を 1 本で押さえる。**

証明は `add_pow_sub_linear_succ`（★`ring` 1 行の恒等式）

`X_{m+1} = (a+b)·X_m + (m+1)·a^m·b²`  （`X_m := (a+b)^{m+1} − a^{m+1} − (m+1)a^m b`）

に強三角不等式を当てる帰納法だけである。★`‖(m+1 : M)‖ ≤ 1` は
`IsUltrametricDist.norm_natCast_le_one`（★#297 のとおり **`M` を明示**して呼ぶ）。

## ★★残り 2 つ（★費用は書かない）

1. ★`𝒪_L = 𝒪_{E₁}[π]` の**相対の冪基底**と成分 `B_j` の抽出 —— ★**未測定**
   （`IntegerResidueBase.lean:adjoin_pi_eq_top` は**絶対版**、と持ち場が名前で挙げた。
   ★私は**開いていない**）。
2. ★`v(B_j) ∈ p·ℤ` —— ★**未測定**。

★★**この 2 つは「成分を取り出す」という 1 つのこと**である。★前波の 5 段の証明のうち、
★(0)(1)(3)(4)(5) は**部品が揃った**（(3) は本波）。★残るのは **(2) の成分抽出だけ**である。

## ★★★★次の波が見落とさないための注記（★持ち場の要請）

★★**前波 `TopIndexSurvives.lean` の docstring に 5 段の証明が入っている。**
★Lean の定理は算術の核（`not_dvd_of_pos_lt` / `delta_deeper_iff_three_le` /
`loss_bound_from_index`）だけで、★**`loss ≤ 2p−2` の Lean 定理はまだ無い。**

★★**残っているのはただ 1 つ** —— `σx − x` を `𝒪_{E₁}` 成分に分けて
`B_{j₀} = w·j₀·f_{j₀} + (深い項)` を取り出すこと。
★その `j₀` が単数であること（`not_dvd_of_pos_lt`）、
★深い項が本当に深いこと（本波 ＋ `delta_deeper_iff_three_le`）、
★上界の算術（`loss_bound_from_index`）は**すべて定理になっている**。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★在庫の測定: 持ち場が挙げた `norm_prod_one_add_sub_one_sub_sum_le` は
   ★**積の版で、冪の版ではなかった**。★「在庫は当たったが使えなかった」例である
   （★名前と形が近いので、開かずに使うと詰まる）。
4. ★宣言名 2 件は **実ファイル**で衝突 0 を確かめた（#348）。
5. ★docstring は `.py` を Write ツールで書いてから実行した（#349）。
-/

namespace ABC3.Found.PGC

namespace BinomialFirstOrder

/-! ## §1 恒等式 —— 1 次の項を取り出した残りの漸化 -/

section Identity

/-- ★核の恒等式（`ring` 1 行）——
`X_{m+1} = (a+b)·X_m + (m+1)·a^m·b²`、ここで `X_m := (a+b)^{m+1} − a^{m+1} − (m+1)a^m b`。 -/
theorem add_pow_sub_linear_succ {M : Type*} [CommRing M] (a b : M) (m : ℕ) :
    (a + b) ^ (m + 2) - a ^ (m + 2) - ((m : M) + 1 + 1) * a ^ (m + 1) * b
      = (a + b) * ((a + b) ^ (m + 1) - a ^ (m + 1) - ((m : M) + 1) * a ^ m * b)
        + ((m : M) + 1) * a ^ m * b ^ 2 := by
  ring

end Identity

/-! ## §2 ノルムの評価 —— 残りは 2 次 -/

section Norm

/-- ★★★**二項展開の 1 次項を取り出した残りは 2 次**（超距離、`‖a‖,‖b‖ ≤ 1`）——
`‖(a+b)^{m+1} − a^{m+1} − (m+1)·a^m·b‖ ≤ ‖b‖²`。

★前波（`TopIndexSurvives`）の 5 段の証明の (2)（主項の同定）と
(4)（`ρ²` 以降が深いこと）を**両方**支える。

★木の `JumpStrictMono.norm_prod_one_add_sub_one_sub_sum_le` は**積の版**で、
★冪の版はこれが初めてである。 -/
theorem norm_add_pow_sub_linear_le {M : Type*} [NormedField M] [IsUltrametricDist M]
    {a b : M} (ha : ‖a‖ ≤ 1) (hb : ‖b‖ ≤ 1) :
    ∀ m : ℕ, ‖(a + b) ^ (m + 1) - a ^ (m + 1) - ((m : M) + 1) * a ^ m * b‖ ≤ ‖b‖ ^ 2
  | 0 => by
      have h0 : (a + b) ^ (0 + 1) - a ^ (0 + 1) - (((0 : ℕ) : M) + 1) * a ^ 0 * b = 0 := by
        push_cast; ring
      rw [h0, norm_zero]
      positivity
  | (m + 1) => by
      have IH := norm_add_pow_sub_linear_le ha hb m
      have hid := add_pow_sub_linear_succ a b m
      have hcast : ((m : M) + 1 + 1) = (((m + 1 : ℕ) : M) + 1) := by push_cast; ring
      rw [hcast] at hid
      rw [hid]
      have hab : ‖a + b‖ ≤ 1 :=
        le_trans (IsUltrametricDist.norm_add_le_max a b) (max_le ha hb)
      have hm : ‖((m : M) + 1)‖ ≤ 1 := by
        have : ((m : M) + 1) = (((m + 1 : ℕ) : M)) := by push_cast; ring
        rw [this]
        exact IsUltrametricDist.norm_natCast_le_one M (m + 1)
      have hpow : ‖a‖ ^ m ≤ 1 := pow_le_one₀ (norm_nonneg a) ha
      refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
      · rw [norm_mul]
        calc ‖a + b‖ * ‖(a + b) ^ (m + 1) - a ^ (m + 1) - ((m : M) + 1) * a ^ m * b‖
            ≤ 1 * ‖b‖ ^ 2 := by
              exact mul_le_mul hab IH (norm_nonneg _) zero_le_one
          _ = ‖b‖ ^ 2 := one_mul _
      · rw [norm_mul, norm_mul, norm_pow, norm_pow]
        calc ‖((m : M) + 1)‖ * ‖a‖ ^ m * ‖b‖ ^ 2
            ≤ 1 * 1 * ‖b‖ ^ 2 := by
              gcongr
          _ = ‖b‖ ^ 2 := by ring

end Norm

/-! ## §3 使っている公理の一覧 -/

#print axioms add_pow_sub_linear_succ
#print axioms norm_add_pow_sub_linear_le

end BinomialFirstOrder

end ABC3.Found.PGC
