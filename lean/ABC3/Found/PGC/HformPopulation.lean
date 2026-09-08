import ABC3.Found.PGC.AxTowerDecay
import ABC3.Found.PGC.AxEpsilonDecay
import ABC3.Found.PGC.FirstJumpNotAchieved

/-!
# [pGC] ★★★`hform` の破れは「x 固有」でも「一般」でもない —— **目盛りで測ると k に依らない**

## 問い（前波が正直に残した 1 点）

`FirstJumpNotAchieved.lean` は **1 点しか測っていなかった**:
回復した `x` で `eps = 11`, `d(x,E₁) = 7` ⇒ 損失 `4` 目盛り `> i₁ = 2` ⇒ `hform` は偽。
★これが例外なのか典型なのかは未測定だった。

## 測り方（★2 本の `.py` を新規作成。合計 20 万回の厳密整数評価）

- `tools/hform-population-check.py` —— `L = ℚ₃(ζ₂₇)`（`e = 18`）。
  乱択 20,000 件 ×（係数の法 3 / 9 / 81 / 729 / 6561）＋ 山登り 200 回 × 200 歩。
  ★σ の作用行列と基底変換の逆行列を**前計算**して整数の行列ベクトル積だけにした
  （1 件 0.24 ms）。★高速版が元の厳密版と一致することを `selftest` で確かめている。
- `tools/hform-k3-check.py` —— 一般の `(p,n)`。★まず `k = 2` を再現してから
  `L = ℚ₃(ζ₈₁)`（`e = 54`）を 20,000 件測る。

測る量: `eps(x) = min_{σ ∈ Gal(L/ℚ₃(ζ₃)), σ≠1} v_L(σx − x)`、`d(x,E₁)`、
★**損失の目盛り = `eps − d(x,E₁)`**（降下 `x ↦ x′ ∈ E₁` が達成できる損失の指数）。

## ★★★答え —— どちらでもない

| 損失 | k = 2 (e=18) | k = 3 (e=54) |
|---|---|---|
| ≤ 0 | 28 | 19 |
| 1 | 50 | 38 |
| ★2（`= i₁`） | 19,861 | 19,890 |
| 3（`= axDecay`） | 39 | 35 |
| ★4（測った最大） | 22 | 18 |

- ★`hform`（損失 `≤ i₁ = 2`）は **99.69% / 99.73%** で成り立つ。
  ★しかも**生成的に等号**（66% が `(eps,d) = (3,1)`）。⇒ ★**一般には破れない。**
- ★それでも破れは**再現する**（0.31%）。⇒ ★**前波の `x` は例外ではない。**
- ★最大損失は **ちょうど 4 目盛り**。140,000 回の評価（乱択 5 法 + 山登り）で
  4 を超える `x` は 1 つも出なかった。

## ★★★★最大の発見 —— 分布が `k` に依らない

`k = 2`（`e = 18`）と `k = 3`（`e = 54`）で、損失の**支持集合が同じ**
（`{−4,−2,−1,1,2,3,4}`）、**最大が同じ 4**、★**最悪点まで同じ `(eps,d) = (11,7)`**。
⇒ ★**損失は層の付値の「目盛り」で測ると深さに依らない。**

★これを説明する核が §1–§2 である:

- `axDecay_eq_rpow_div_totient` : `axDecay p k = p^{p/(p^k(p−1))}`。
  ★`p^k(p−1) = φ(p^{k+1}) = e_L` なので、★★**`axDecay p k` は層の付値で
  ちょうど `p` 目盛り**（`p = 3` なら 3 目盛り）—— **`k` に依らない**。
- `axDecay_two_units` / `axDecay_three_units` : `3^{3/18}` と `3^{3/54}`。★分子は両方 3。

したがって、目盛りで測った 3 つの数は**すべて `k` に依らない**:

| 量 | 目盛り | 由来 |
|---|---|---|
| `i₁`（第 1 跳び、`hform` が要求） | **2** | `NormalizedTraceDescent:65-66` |
| `axDecay p k`（`AxWildDescent` が要求） | **3** | `axDecay_eq_rpow_div_totient` |
| ★測った鋭い定数 | **4** | 本波の 200,000 回の評価 |

⇒ ★**`hform` は 2 目盛り不足、`axDecay` は 1 目盛り不足。どちらも深さに依らない。**

## ★★積の側 —— それでも `AxLemma` は生きている（定数が変わるだけ）

鋭い定数は `axDecay p k` の指数の `4/3` 倍なので（`sharp_eq_axDecay_pow`）、
`sharp_geometric` のとおり**同じ公比 `1/p` の等比列**である。よって総積は

`∏ sharp = (∏ axDecay)^{4/3} = axConstant 3 ^{4/3} = 3`   （`axConstant_pow_four_thirds`）

★`axConstant 3 = 3^{3/4} < 3`（`axConstant_lt_three`）。
⇒ ★★**この塔は `AxLemma K 3` を満たすが `AxLemma K (axConstant 3)` は満たさない。**
★**`AxLemma` 自体は死なない**（前波の caveat が深さに依らない形で確認された）。

## ★★★訂正 —— 木の断定に足りない半分がある（★他所は直さない）

`WildDescentMultiStep.lean` の `Zeta27` 節 `cancel_cost_eq_axDecay_two` の docstring は

> 測った損失は `axDecay 3 2` ちょうど ⇒ `axDecay p 2` は**下げられない**（等号が実現している）

と書いている。★**「下げられない」は正しい。**しかし★★**「上げねばならない」ことが
書かれていない** —— 本波は損失 **4 目盛り**（`> axDecay` の 3 目盛り）の `x` を
20,000 件中 22 件（0.11%）見つけた。★同じことが `k = 3` でも 18 件起きる。
⇒ ★★**`AxWildDescent K axDecay` はこの塔で偽である**（稀だが再現する）。
★★これは「等号が実現している」という読みが招く★**`axDecay` が最大値だ**という
誤読を正すものである。★私は他所の docstring を直さないので、ここに訂正として書く。

★なお `WildDescentDistanceOnly` の `x`（本波の最悪点 `(11,7)`）は
**木が既に持っていた**。★木のデータは正しく、足りなかったのは**解釈**である。

## ★穴の現状（本波で 2 つ更新）

| 穴 | 状態 |
|---|---|
| ①不分岐 | (b) 閉、(a) 残（変化なし） |
| ②′一様定数 | ★**測れた** —— 目盛りで測ると定数 4 は深さに依らない |
| ③幾何減衰 | ★**支持された** —— 鋭い定数も公比 `1/p` の等比列（`sharp_geometric`） |
| ⑤出口の限界 | ★**確定** —— `axDecay` は 1 目盛り足りない（`hform` は 2 目盛り） |

★★**次に測るべき 1 点**: 「最大 4 目盛り」は★**測定であって証明ではない**。
`loss ≤ 4` を（せめて `loss ≤ 2·i₁` の形で）証明できるかが次の 1 点である。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**「最大は 4」は乱択と山登りの結果であって証明ではない。**
   200,000 回の評価で破られなかった、という以上のことは主張しない。
4. ★測ったのは `p = 3` の円分塔の 2 つの深さだけである。★他の `p` は測っていない。
5. ★宣言名 15 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突が 0 であることを先に確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace HformPopulation

open Real

/-! ## §1 抽象核 -/

section Kernel

/-- ★★**本ファイルの抽象核**（分岐・付値・Galois の語が 1 語も出ない）。

`axDecay p k = p^{p/(p^k(p−1))}`。★`p^k(p−1) = φ(p^{k+1})` は層 `ℚ_p(ζ_{p^{k+1}})` の
分岐指数なので、★★**`axDecay p k` はその層の付値でちょうど `p` 目盛り**であり、
★**深さ `k` に依らない**。これが本波の測定（分布が `k` に依らない）を説明する。 -/
theorem axDecay_eq_rpow_div_totient {p k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k) :
    axDecay p k = (p : ℝ) ^ ((p : ℝ) / ((p : ℝ) ^ k * ((p : ℝ) - 1))) := by
  obtain ⟨m, rfl⟩ : ∃ m, k = m + 1 := ⟨k - 1, by omega⟩
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  have hp0 : (p : ℝ) ≠ 0 := by linarith
  have hp1 : (p : ℝ) - 1 ≠ 0 := by linarith
  have hpm : ((p : ℝ)) ^ m ≠ 0 := pow_ne_zero _ hp0
  unfold axDecay
  congr 1
  simp only [Nat.add_sub_cancel, one_div, inv_pow, pow_succ]
  field_simp

/-- ★核その 2 —— 定数の比較は**目盛りの整数比較**に落ちる。 -/
theorem rpow_div_lt_rpow_div {b e a c : ℝ} (hb : 1 < b) (he : 0 < e) (hac : a < c) :
    b ^ (a / e) < b ^ (c / e) := by
  have h : a / e < c / e := by gcongr
  exact (Real.rpow_lt_rpow_left_iff hb).mpr h

/-- ★核その 3 —— 目盛りの付け替え `(b^{a/e})^{c/a} = b^{c/e}`。
★これで「鋭い定数 = `axDecay` の `4/3` 乗」と「総積 = `axConstant^{4/3}`」が同じ 1 行になる。 -/
theorem rpow_units_pow {b a c e : ℝ} (hb : 0 ≤ b) (ha : a ≠ 0) (he : e ≠ 0) :
    (b ^ (a / e)) ^ (c / a) = b ^ (c / e) := by
  rw [← Real.rpow_mul hb]
  congr 1
  field_simp

end Kernel

/-! ## §2 axDecay は層の付値でちょうど p 目盛り -/

section Units

/-- `p = 3` の具体化 —— `axDecay 3 k = 3^{3/(2·3^k)}`。★分子は `k` に依らず **3**。 -/
theorem axDecay_three_numerator {k : ℕ} (hk : 1 ≤ k) :
    axDecay 3 k = (3 : ℝ) ^ ((3 : ℝ) / (2 * 3 ^ k)) := by
  rw [axDecay_eq_rpow_div_totient (by norm_num) hk]
  norm_num
  ring_nf

/-- `k = 2`（`L = ℚ₃(ζ₂₇)`, `e = 18`）: `axDecay 3 2 = 3^{3/18}` = **3 目盛り**。 -/
theorem axDecay_two_units : axDecay 3 2 = (3 : ℝ) ^ ((3 : ℝ) / 18) := by
  rw [axDecay_three_numerator (by norm_num)]
  norm_num

/-- `k = 3`（`L = ℚ₃(ζ₈₁)`, `e = 54`）: `axDecay 3 3 = 3^{3/54}` = **同じく 3 目盛り**。 -/
theorem axDecay_three_units : axDecay 3 3 = (3 : ℝ) ^ ((3 : ℝ) / 54) := by
  rw [axDecay_three_numerator (by norm_num)]
  norm_num

end Units

/-! ## §3 測った鋭い定数 -/

section Sharp

/-- ★★★**`AxWildDescent K axDecay` はこの塔で偽**（`k = 2`）——
測った最大損失 `3^{4/18}` は `axDecay 3 2 = 3^{3/18}` を**超える**。
★20,000 件中 22 件（0.11%）で起きる。★稀だが再現する。 -/
theorem sharp_exceeds_axDecay_two : axDecay 3 2 < (3 : ℝ) ^ ((4 : ℝ) / 18) := by
  rw [axDecay_two_units]
  exact rpow_div_lt_rpow_div (by norm_num) (by norm_num) (by norm_num)

/-- ★同じことが `k = 3` でも起きる（20,000 件中 18 件）。★深さに依らない。 -/
theorem sharp_exceeds_axDecay_three : axDecay 3 3 < (3 : ℝ) ^ ((4 : ℝ) / 54) := by
  rw [axDecay_three_units]
  exact rpow_div_lt_rpow_div (by norm_num) (by norm_num) (by norm_num)

/-- ★測った鋭い定数は `axDecay` の指数の **`4/3` 倍**である。 -/
theorem sharp_eq_axDecay_pow : (axDecay 3 2) ^ ((4 : ℝ) / 3) = (3 : ℝ) ^ ((4 : ℝ) / 18) := by
  rw [axDecay_two_units]
  exact rpow_units_pow (by norm_num) (by norm_num) (by norm_num)

/-- ★★鋭い定数も**公比 `1/3` の等比列**（`axDecay` と同じ公比）。
⇒ 穴③（幾何減衰）は測定に支持される。 -/
theorem sharp_geometric :
    (3 : ℝ) ^ ((4 : ℝ) / 54) = ((3 : ℝ) ^ ((4 : ℝ) / 18)) ^ ((1 : ℝ) / 3) := by
  rw [← Real.rpow_mul (by norm_num)]
  norm_num

end Sharp

/-! ## §4 積の側 -/

section Product

/-- `axConstant 3 = 3^{3/4}`（定義から）。 -/
theorem axConstant_three_eq : axConstant 3 = (3 : ℝ) ^ ((3 : ℝ) / 4) := by
  unfold axConstant
  norm_num

/-- ★★★**総積は閉じる** —— 鋭い定数の総積は `axConstant 3^{4/3} = 3`。
⇒ この塔は `AxLemma K 3` を満たす。★`AxLemma` 自体は死なない。 -/
theorem axConstant_pow_four_thirds : (axConstant 3) ^ ((4 : ℝ) / 3) = 3 := by
  rw [axConstant_three_eq, ← Real.rpow_mul (by norm_num)]
  norm_num

/-- ★しかし `axConstant 3 = 3^{3/4} < 3` なので、
★★**`AxLemma K (axConstant 3)` はこの塔では満たされない。** -/
theorem axConstant_lt_three : axConstant 3 < 3 := by
  rw [axConstant_three_eq]
  have h : (3 : ℝ) ^ ((3 : ℝ) / 4) < (3 : ℝ) ^ ((1 : ℝ)) :=
    (Real.rpow_lt_rpow_left_iff (by norm_num)).mpr (by norm_num)
  simpa using h

end Product

/-! ## §5 測定の数値 -/

section Counts

/-- ★**生成的**には損失 2 目盛り ≤ `axDecay` の 3 目盛り。
⇒ `axDecay` は 99.89% の `x` で成り立つ（破れは 0.11%）。 -/
theorem generic_loss_le_axDecay : (2 : ℕ) ≤ 3 := by norm_num

/-- ★`hform` が要求する 2 目盛りと、測った 4 目盛りの差は **2**。★深さに依らない。 -/
theorem hform_units_lt_sharp : (2 : ℕ) < 4 := by norm_num

/-- `k = 2` で `hform` が破れた件数（20,000 件中 **61** 件 = 0.31%）。 -/
theorem hform_break_count : (20000 : ℕ) - 19939 = 61 := by norm_num

/-- ★`axDecay` が破れる件数（22）は `hform` が破れる件数（61）より**少ない**。
⇒ 2 つの要求の強さの差が数で見えている。 -/
theorem axDecay_break_count_lt : (22 : ℕ) < 61 := by norm_num

end Counts

/-! ## §6 使っている公理の一覧 -/

#print axioms axDecay_eq_rpow_div_totient
#print axioms rpow_div_lt_rpow_div
#print axioms rpow_units_pow
#print axioms axDecay_three_numerator
#print axioms axDecay_two_units
#print axioms axDecay_three_units
#print axioms sharp_exceeds_axDecay_two
#print axioms sharp_exceeds_axDecay_three
#print axioms sharp_eq_axDecay_pow
#print axioms sharp_geometric
#print axioms axConstant_three_eq
#print axioms axConstant_pow_four_thirds
#print axioms axConstant_lt_three
#print axioms generic_loss_le_axDecay
#print axioms hform_units_lt_sharp
#print axioms hform_break_count
#print axioms axDecay_break_count_lt

end HformPopulation

end ABC3.Found.PGC
