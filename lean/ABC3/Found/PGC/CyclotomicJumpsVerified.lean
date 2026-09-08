import ABC3.Found.PGC.CyclotomicNumbersVerified

/-!
# [pGC] ★再現できない数値を**数えて**掃討を始めた —— 円分塔の跳びを全部検算した

## ★どれを選んだか、なぜか（持ち場が選択を任せた点）

★**(N)（同じ形の穴の掃討）を選んだ。** ★まず**数えた**:

```
grep -rn "scratchpad" lean/ABC3/Found/PGC/*.lean | wc -l     →  8 箇所
grep -rln "scratchpad" lean/ABC3/Found/PGC/*.lean            →  ★5 ファイル
  DeepDescentPairDirect / DeepDescentRepair / EquivariantProjectionDescent /
  GainedTowerDescent / JumpDefectTradeoff
ls tools/*.py  →  そこで名指しされているスクリプト（pairsearch.py / three.py /
                  sharp.py / trade.py …）は ★1 本も無い
```

⇒ ★★**前波で塞いだのと同じ形の穴が、あと 5 ファイル・8 箇所ある。**
★0〜1 箇所なら「無い」と記録して別候補に移る予定だったが、★**5 ファイルあったので掃討に入った**。

★その中で★**本波の道具（`tools/zeta27-*.py`）がそのまま効く**のは
「円分塔の跳び」の数値だったので、そこから始めた（`tools/zeta-tower-check.py` を新規作成。
`(p, n)` を引数にした一般版）。

## ★★★検算の結果 —— **3 つとも一致した**

| 木の字面 | 本波の再導出（`tools/zeta-tower-check.py`） | 一致 |
|---|---|---|
| `JumpFromValueGroup.lean:332 harith_zeta81` の `u = (2, 8, 26)` | `ℚ₃(ζ₈₁)` の wild break は **`[2, 8, 26]`** | ✓ |
| 同 `:323` の `E = 54` | `e = deg E = 54` | ✓ |
| `GainedTowerDescent.lean:411 zeta27_sharp` の `L = 8` | `ℚ₃(ζ₂₇) ⊃ ℚ₃(ζ₉)` の層の跳び **`8`** | ✓ |
| `GainedTowerDescent.lean:417 zeta81_sharp` の `L = 26` | `ℚ₃(ζ₈₁) ⊃ ℚ₃(ζ₂₇)` の層の跳び **`26`** | ✓ |
| `WildDescentDistanceOnly.lean:60-62` の `ℚ₃(ζ₂₇)` の分岐群 | 前波と同じ（`[18,9,9,3,3,3,3,3,3,1]`） | ✓ |

★★おまけ: `ℚ₂(ζ₁₆)/ℚ₂` の break は **`[1, 3, 7]`**、`e = 8`、`|G_u| = [8,8,4,4,2,2,2,2,1]`。

## ★★測った convention の食い違い（★結論の訂正ではない。名指しで記録する）

`GainedTowerDescent.lean:423 zeta16_sharp` は `sharpTwoCost 2 4 7 = 8`、
すなわち `(S, T) = (4, 7)` で `max(T, p·S) = max(7, 8) = 8` を使っている。

★本波の計算では `ℚ₂(ζ₁₆)/ℚ₂` の **break** は `[1, 3, 7]` であり、`4` は break ではない
（`4` は `Gal(L/ℚ₂(ζ₄))` の **`min i`** の方である。break は `min i − 1 = 3`）。
★`ℚ₃` の 2 行（`zeta27_sharp` の `(2,8)`、`zeta81_sharp` の `(8,26)`）は
★**どちらも break** なので、★**`ℚ₂(ζ₁₆)` の行だけ convention が違う**。

* break で揃えると `max(7, 2·3) = 7`（§2 `zeta16_break_convention`）
* 木の値は `max(7, 2·4) = 8`（§2 `zeta16_mini_convention`）
* ★2 つは異なる（§2 `zeta16_conventions_differ`）

★★**ただし木自身が同 `:421-423` で「★非巡回 `ℤ/2×ℤ/4`、本ファイルの適用外」
「★これは測定であって定理ではない」と明記している**ので、
★これは**結論の誤りではなく、適用外の行の目盛りの違い**である。★そのことを記録する。

## ★本ファイルが Lean に固定したもの

`p` 進の付値は塔を組まないと Lean に載らないので、★**算術と群の側だけ**を固定した:

* `zeta81_hasseArf_first` / `_second` —— `3 ∣ 8−2`, `9 ∣ 26−8`（★Hasse–Arf の合同）
* `zeta81_strict_mono` / `zeta81_upper` / `zeta81_first_jump_pos` ——
  `harith` の 4 条件が `u = (2,8,26)`, `e = 2` で成り立つ（★`harith_zeta81` の数値の裏づけ）
* `zeta16_hasseArf_first` / `_second` —— `2 ∣ 3−1`, `4 ∣ 7−3`（★`ℚ₂(ζ₁₆)` の break も合同を満たす）
* `zeta16_break_convention` / `_mini_convention` / `_conventions_differ` —— 上の目盛りの違い
* `zeta81_sigma_order` / `zeta81_tau_order` / `zeta16_units_not_cyclic` —— 群の側

## ★残っている同じ形の穴（★数えた。次の波のために書く）

| ファイル | 名指しされたスクリプト | 本波で塞いだか |
|---|---|---|
| `WildDescentDistanceOnly.lean` | （報告に添付、紛失） | ★前波で塞いだ |
| `GainedTowerDescent.lean:94` | `sharp/{sharp,rec,rec2,closed,probe}.py` | ★**跳びの数値だけ**塞いだ（`L = 8, 26`）。総当たりの表は未 |
| `DeepDescentPairDirect.lean:37,51,147` | `pair/pairsearch.py`, `run2.py`, `awd/cyc.py` | ★未 |
| `DeepDescentRepair.lean:110` | （再現スクリプト） | ★未 |
| `EquivariantProjectionDescent.lean:19` | `grad/three.py`, `trnorm.py`, `maxloss.py` | ★未 |
| `JumpDefectTradeoff.lean:106` | `trade/trade.py`, `trade2.py`, `closed2.py` | ★未 |

★★`DeepDescentPairDirect` / `DeepDescentRepair` は `ℤ[π]/(g)` 上の乱択探索なので、
★**本波と前波の 3 本の `.py` がほぼそのまま使える**（次の波の最安の候補である）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★本ファイルは原典の主張ではなく**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★`p` 進の付値・分岐群そのものは Lean に載せていない（塔の構成が要るため）。
4. ★`ℚ₂(ζ₁₆)` の convention の違いは★**木の結論を訂正するものではない**
   （木自身が「適用外」「測定であって定理ではない」と書いている）。
-/

namespace ABC3.Found.PGC

namespace CyclotomicJumps

/-! ## §1 `ℚ₃(ζ₈₁)` の跳び `u = (2, 8, 26)` の算術（★本波で再導出した） -/

section Zeta81

/-- ★`u = (2,8,26)` は Hasse–Arf の合同を満たす: `3 ∣ 8−2`。 -/
theorem zeta81_hasseArf_first : (3 : ℤ) ^ 1 ∣ (8 : ℤ) - 2 := by decide

/-- `9 ∣ 26−8`。 -/
theorem zeta81_hasseArf_second : (3 : ℤ) ^ 2 ∣ (26 : ℤ) - 8 := by decide

theorem zeta81_strict_mono : (2 : ℤ) < 8 ∧ (8 : ℤ) < 26 := by
  constructor <;> norm_num

/-- ★`harith` の上界 `(p−1)·u k ≤ p^{k+1}·e`（`e = 2`, `k = 2`）: `52 ≤ 54`。 -/
theorem zeta81_upper : ((3 : ℤ) - 1) * 26 ≤ (3 : ℤ) ^ (2 + 1) * 2 := by norm_num

theorem zeta81_first_jump_pos : (1 : ℤ) ≤ 2 := by norm_num

end Zeta81

/-! ## §2 `ℚ₂(ζ₁₆)` の跳び `[1, 3, 7]` と目盛りの違い -/

section Zeta16

theorem zeta16_hasseArf_first : (2 : ℤ) ^ 1 ∣ (3 : ℤ) - 1 := by decide

theorem zeta16_hasseArf_second : (2 : ℤ) ^ 2 ∣ (7 : ℤ) - 3 := by decide

/-- ★break で揃えると `max(T, p·S) = max(7, 2·3) = 7`。 -/
theorem zeta16_break_convention : max (7 : ℕ) (2 * 3) = 7 := by norm_num

/-- ★木の `zeta16_sharp` の値は `max(7, 2·4) = 8`（`4` は break ではなく `min i`）。 -/
theorem zeta16_mini_convention : max (7 : ℕ) (2 * 4) = 8 := by norm_num

/-- ★★2 つの目盛りは異なる。★ただし木自身が `GainedTowerDescent.lean:421-423` で
「非巡回・本ファイルの適用外」「測定であって定理ではない」と明記しているので、
★これは**結論の誤りではなく目盛りの違い**である。 -/
theorem zeta16_conventions_differ : max (7 : ℕ) (2 * 3) ≠ max (7 : ℕ) (2 * 4) := by norm_num

end Zeta16

/-! ## §3 群の側（`decide`） -/

section GroupSide

theorem zeta81_sigma_order : (4 : ZMod 81) ^ 27 = 1 := by decide

theorem zeta81_tau_order : (28 : ZMod 81) ^ 3 = 1 := by decide

theorem zeta16_units_not_cyclic : (7 : ZMod 16) ^ 2 = 1 ∧ (3 : ZMod 16) ^ 4 = 1 := by decide

end GroupSide

/-! ## §4 使っている公理の一覧 -/

#print axioms zeta81_hasseArf_first
#print axioms zeta81_hasseArf_second
#print axioms zeta81_strict_mono
#print axioms zeta81_upper
#print axioms zeta81_first_jump_pos
#print axioms zeta16_hasseArf_first
#print axioms zeta16_hasseArf_second
#print axioms zeta16_break_convention
#print axioms zeta16_mini_convention
#print axioms zeta16_conventions_differ
#print axioms zeta81_sigma_order
#print axioms zeta81_tau_order
#print axioms zeta16_units_not_cyclic

end CyclotomicJumps

end ABC3.Found.PGC
