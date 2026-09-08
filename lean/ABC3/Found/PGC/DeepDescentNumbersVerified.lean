import ABC3.Found.PGC.CyclotomicJumpsVerified

/-!
# [pGC] ★★★★`DeepDescentRepair` の反例は**完全に再現できた** —— 否定的結果の土台は在る

## ★どれを選んだか、なぜか（前波の表から）

★**`DeepDescentRepair` を選んだ。** 理由:

1. ★★**掛け金が最も高い。** 本日の台帳は「`AxDeepDescent` が false であることを
   機械が計算した反例で示した」と記録しており、
   ★**反例が再現できなければその否定的結果の土台が無い**。
2. ★★**最も安い。** 同ファイル `:35` は `x` を★**明示している**
   （`x := π³ + 2·π_E + π_E⁵`）。★前波で失われていた `WildDescentDistanceOnly` の `x` と違い、
   ★**乱択探索が要らない**。しかも `d(x, ℚ₃(ζ₉))` / `d(x, ℚ₃(ζ₃))` は
   ★前波で書いた `tools/zeta27-distance-check.py` の関数が**そのまま使える**。

## ★★★★結果 —— **表の 4 値すべてが一致した**

`tools/deepdescent-repair-check.py`（★本波で新規作成。`zeta27-distance-check.py` を import）。

| `DeepDescentRepair.lean:37-42` の字面 | 本波の再導出 | 一致 |
|---|---|---|
| `[F(x):F] = 9`（⇒ `wildDepth = 2`） | 軌道の大きさ **9** | ✓ |
| `min_{σ≠1} v_L(σx − x) = 23`（8 元: `23,23,27,23,23,27,23,23`） | ★**`a=4:23, 7:23, 10:27, 13:23, 16:23, 19:27, 22:23, 25:23`**（★**並びまで一致**） | ✓ |
| `d(x, E₁) = 19` | **19** | ✓ |
| `d(x, F) = 15` | **15** | ✓ |
| `v_L(π_E) = 3` | **3** | ✓ |

⇒ 要求は `v_L(x − x') ≥ 23 − 3 = 20`（`axDecay 3 2 = 3^{3/18}`）だが
`19 < 20` かつ `15 < 20` ⇒ ★★★**`L` の中に証人は無い**（同 `:44` と一致）。

★★**したがって `deepDescent_to_E1_false`（`:276`）と `deepDescent_to_base_false`（`:282`）が
使っている指数 `19` / `15` / `23` は正しく、★否定的結果の土台は在る。**

## ★`JumpDefectTradeoff` は**スクリプト無しで検算できた**（1 箇所塞いだ）

同 `:106` は `.../scratchpad/trade/{trade,trade2,closed2}.py` を名指ししているが、
★**主結果はすべて既に Lean の定理**である（`JumpDefect.defect_eq` / `covered_iff` /
`witness_five_*` / `min_cost_fits_of_le_three`）。
★スクリプトが支えているのは**総当たりの表**（`43,188` 本で破れ 0、`p=5` で 118 件）だけで、
★**結論そのものは閉じた形で出ている**。

§2 はその閉じた形 `覆う ⟺ (p−1)(p−4) ≤ 0` を素ごとに固定した:
`p = 2, 3` は覆い、`p = 5, 7` は覆わない。⇒ ★**「`p ≤ 3` なら覆う」は式で確認できる。**

## ★掃討の現状（前波の表を更新）

| ファイル | 状態 |
|---|---|
| `WildDescentDistanceOnly.lean` | ★**塞いだ**（前々波、`x` を回復） |
| `GainedTowerDescent.lean:94` | ★跳びの数値は塞いだ（前波、`L = 8, 26`）。総当たりの表は未 |
| ★`DeepDescentRepair.lean:110` | ★★**塞いだ**（本波。4 値すべて一致） |
| ★`JumpDefectTradeoff.lean:106` | ★★**塞いだ**（本波。★結論は Lean の定理で、スクリプトは補助的な総当たりだけ） |
| `DeepDescentPairDirect.lean:37,51,147` | ★未（`ℤ[π]/(g)` 上の 95,298 件の乱択。本波の道具で届く） |
| `EquivariantProjectionDescent.lean:19` | ★未 |

⇒ ★**6 箇所中 4 箇所が塞がった。**

## ★数値が合わない箇所は**見つからなかった**（正直に）

★持ち場は「合わない箇所を見つけたら止めて報告」と書いたが、★**本波で測った範囲では
すべて一致した**。★前波（円分塔の跳び）でも一致し、前々波（`WildDescentDistanceOnly`）でも
一致した。⇒ ★**木の機械計算の記録は、測った範囲では正確である。**

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★本ファイルは**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★`p` 進の付値そのものは Lean に載せていない。★載せたのは
   再導出した数値の**算術的な帰結**（`19 < 20` など）だけである。
4. ★`DeepDescentRepair.lean:59-62` が自認する未測定点（「`x'` が `L` の外に在る可能性は
   排除していない」）は★**本波でも排除していない**。★同ファイルの結論は
   「降下先を Galois 閉包に取る限り偽」であって、それ以上ではない。
-/

namespace ABC3.Found.PGC

namespace DeepDescentVerified

/-! ## §1 `DeepDescentRepair` の反例の数値（★本波で全部再導出した） -/

section Repair

theorem eps_valuation : (23 : ℕ) = 23 := rfl

/-- ★要求される付値 `23 − 3 = 20`（`ε` の `v` が `23`、`axDecay 3 2 = 3^{3/18}`）。 -/
theorem required_valuation : (23 : ℤ) - 3 = 20 := by norm_num

/-- ★★`d(x, ℚ₃(ζ₉))` の `v` は `19 < 20` ⇒ `E₁` へは届かない
（`DeepDescentRepair.lean:276 deepDescent_to_E1_false` の数値。★本波で再導出した）。 -/
theorem to_E1_fails : (19 : ℤ) < 20 := by norm_num

/-- ★★`d(x, ℚ₃(ζ₃))` の `v` は `15 < 20` ⇒ 底へも届かない（同 `:282`）。 -/
theorem to_base_fails : (15 : ℤ) < 20 := by norm_num

theorem orbit_size_nine : (9 : ℕ) = 3 ^ 2 := by norm_num

/-- 軌道 `9 = 3²` ⇒ `wildDepth = 2`。 -/
theorem wildDepth_two : padicValNat 3 (3 ^ 2) = 2 := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  exact padicValNat.prime_pow 2

/-- 8 元の変位 `{23,23,27,23,23,27,23,23}` の最小は `23`。 -/
theorem displacement_min : min (23 : ℕ) 27 = 23 := by norm_num

end Repair

/-! ## §2 `JumpDefectTradeoff` の閉じた形（★スクリプト不要） -/

section Tradeoff

/-- ★`JumpDefect.covered_iff` の閉じた形 `(p−1)(p−4) ≤ 0` を素ごとに固定する。
★スクリプトが無くても検算できる（本波で 1 箇所塞いだ根拠）。 -/
theorem covered_two : ((2 : ℤ) - 1) * ((2 : ℤ) - 4) ≤ 0 := by norm_num

theorem covered_three : ((3 : ℤ) - 1) * ((3 : ℤ) - 4) ≤ 0 := by norm_num

/-- ★`p = 5` では覆わない（`4·1 = 4 > 0`）。 -/
theorem not_covered_five : ¬ (((5 : ℤ) - 1) * ((5 : ℤ) - 4) ≤ 0) := by norm_num

theorem not_covered_seven : ¬ (((7 : ℤ) - 1) * ((7 : ℤ) - 4) ≤ 0) := by norm_num

end Tradeoff

/-! ## §3 使っている公理の一覧 -/

#print axioms required_valuation
#print axioms to_E1_fails
#print axioms to_base_fails
#print axioms orbit_size_nine
#print axioms wildDepth_two
#print axioms displacement_min
#print axioms covered_two
#print axioms covered_three
#print axioms not_covered_five
#print axioms not_covered_seven

end DeepDescentVerified

end ABC3.Found.PGC
