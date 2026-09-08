import ABC3.Found.PGC.KrasnerCeiling

/-!
# [pGC] ★★★再現性の穴を塞いだ —— `ℚ₃(ζ₂₇)` の数値を**全部**厳密整数演算で再導出した

## ★どれを選んだか、なぜか（持ち場が選択を任せた点）

★**(R)（再現性の穴）を選んだ。** 理由:

* ★本日この鎖は `WildDescentDistanceOnly.lean:60-66` の数値
  （`e_M = 18` / `i(σ)=3` / `i(τ)=9` / 層の跳び `i=8` / `v(σx−x)=11` / `v(τx−x)=15` /
  `d(x,F)=7` / `d(x,K)=3`）を★**5 回引用**した。
* ★前波で測ったとおり、同 `:57` は「厳密な整数演算で全部当たった
  （★**スクリプトは報告に添付**）」と書いており、★`x` もスクリプトも**リポジトリに無い**。
  ⇒ ★**引用の土台が再現できない状態**であり、間違っていれば本日の測定のいくつかが崩れる。
* ★(S)（棚卸し）は各ファイルの docstring が既に `file:line` を持っており重複が大きい。
  (A)/(B)/(C) は 2〜3 波前に「1 波では入らない」と測った。

## ★★★結果 —— **8 つの数値がすべて一致した**

`tools/zeta27-ramification-check.py` と `tools/zeta27-distance-check.py`（★本波で新規作成、
★丸めゼロ・すべて Python の任意精度整数）。設定は `π = ζ₂₇ − 1`、
`E(π) = Φ₂₇(1+π)` は Eisenstein（`v₃(E₀)=1` を機械で確認）、`𝒪_M = ℤ₃[π]`、
`v_M(Σ bᵢ πⁱ) = minᵢ (18·v₃(bᵢ) + i)`。

| 木の字面（`WildDescentDistanceOnly.lean`） | 本波の再導出 | 一致 |
|---|---|---|
| `e_M = 18` | `E` は deg 18 の Eisenstein | ✓ |
| `G_0 = G`（18 個） | `{a ∈ (ℤ/27)ˣ}` | ✓ |
| `G_1 = G_2 = H`（9 個） | `{1,4,7,10,13,16,19,22,25}` ＝ `a ≡ 1 mod 3` ＝ `Gal(M/ℚ₃(ζ₃))` | ✓ |
| `G_3..G_8 = Gal(M/ℚ₃(ζ₉))` | `{1,10,19}` ＝ `a ≡ 1 mod 9` | ✓ |
| `G_9 = 1` | ✓ | ✓ |
| 位数 9 の `σ` は `i_G(σ) = 3` | `a = 4,7,13,16,22,25` すべて `3` | ✓ |
| 位数 3 の `τ` は `i_G(τ) = 9` | `a = 10,19` ともに `9` | ✓ |
| 層 `M/ℚ₃(ζ₉)` の跳び `i = 8` | `Gal(M/F)` は `G_8` に居て `G_9` に居ない | ✓ |
| `x = π` で `v(σx−x)=3`, `v(τx−x)=9` | ✓ | ✓ |
| ★実在する `x` で `v(σx−x)=11`, `v(τx−x)=15` | ★**見つかった**（下記） | ✓ |
| ★`d(x, ℚ₃(ζ₉)) = 7` | ★**7**（`𝒪_M = 𝒪_F[π]` の基底 `{μ^s πⁱ}` で厳密に分解） | ✓ |
| ★`d(x, ℚ₃(ζ₃)) = 3`（同 `:100`） | ★**3** | ✓ |
| 出現率 `37/25000`, `9/12000`, `11/12000` | `21/20000`（同じ桁） | ✓ |

★★**失われていた `x` を回復した**（π-basis、低次から）:

```
x = [71, 78, 9, 50, 30, 54, 67, 43, 40, 3, 25, 11, 40, 44, 74, 63, 38, 24]
    （すなわち x = Σ_{j<18} b_j (ζ₂₇ − 1)^j、係数は 0..80）
    v(σ₄ x − x) = 11,  v(σ₁₀ x − x) = 15,  d(x, ℚ₃(ζ₉)) = 7,  d(x, ℚ₃(ζ₃)) = 3
```

⇒ ★★★**木の数値は全部正しかった。本日の測定は崩れない。**
★そして★**今後は `tools/` の 2 本で誰でも再現できる**（`file:line` ではなくスクリプトで）。

## ★★見つけた自分のバグ（記録する）

★1 度目の実装で `μ = (1+π)^9 − 1` を「`ζ₉ − 1`」だと思い込み、基底が退化して
`AssertionError: singular` になった。★正しくは `(1+π)^9 − 1 = ζ₂₇^9 − 1 = ζ₃ − 1`
（`v_M = 9`）で、`ζ₉ − 1` は `(1+π)^3 − 1`（`v_M = 3`）である。
★**機械が `singular` で止めてくれたので気づいた**（`v_M(μ) = 9` を印字して確定した）。

## ★本ファイルが Lean に固定したもの

★`p` 進の付値そのものは塔を組まないと Lean に載らないので、
★**群の側と算術の側だけ**を `decide` / `norm_num` で固定した:

* `sigma_order_nine` / `sigma_order_nine_not_three` —— `σ₄` は位数 9
* `tau_order_three` / `tau_ne_one` —— `σ₁₀` は位数 3
* `sigma_cubed_eq_tau` —— `4³ ≡ 10 (mod 27)`（★`τ = σ³`）
* `sigma_mod_three` —— `4 ≡ 1 (mod 3)`（★`σ ∈ Gal(M/ℚ₃(ζ₃))`）
* `layer_step_fails` —— `7 < 8`（★`F` を通る 1 段は**ちょうど 1 だけ**足りない）
* `uniformizer_gain_ok` / `generic_gain_fails` —— 素元は通り、実在の `x` は通らない
* `composite_has_slack` —— `11 − 12 ≤ 3`（★合成は無傷）
* `krasner_does_not_fire` —— `8 ≤ 11`（★前波の `KrasnerCeiling` の数値条件）

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★本ファイルは原典の主張ではなく**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★`p` 進の付値・分岐群そのものは Lean に載せていない（塔の構成が要るため）。
   ★載せたのは `(ℤ/27)ˣ` の群の事実と `ℕ`/`ℤ` の不等式だけである。
4. ★乱択探索は `random.seed(20260909)` で固定した（★再現可能）。
-/

namespace ABC3.Found.PGC

namespace CyclotomicVerified

/-! ## §1 群の側（`(ℤ/27)ˣ`）—— `decide` で固定 -/

section GroupSide

/-- `σ = σ₄` は `Gal(M/ℚ₃)` で位数 9（`4⁹ ≡ 1 mod 27`）。 -/
theorem sigma_order_nine : (4 : ZMod 27) ^ 9 = 1 := by decide

theorem sigma_order_nine_not_three : (4 : ZMod 27) ^ 3 ≠ 1 := by decide

/-- `τ = σ₁₀` は位数 3（`10³ ≡ 1 mod 27`）。 -/
theorem tau_order_three : (10 : ZMod 27) ^ 3 = 1 := by decide

theorem tau_ne_one : (10 : ZMod 27) ≠ 1 := by decide

/-- ★`τ = σ³`（`4³ = 64 ≡ 10 mod 27`）。 -/
theorem sigma_cubed_eq_tau : (4 : ZMod 27) ^ 3 = 10 := by decide

/-- ★`σ₄ ∈ Gal(M/ℚ₃(ζ₃))`（`4 ≡ 1 mod 3`）。 -/
theorem sigma_mod_three : (4 : ZMod 27) - 1 = 3 := by decide

end GroupSide

/-! ## §2 算術の側 —— `norm_num` で固定 -/

section Numbers

/-- ★★`d(x, ℚ₃(ζ₉)) = 7` は要求 `8` に**ちょうど 1 だけ**足りない
（★本波で `d(x,F) = 7` を厳密整数演算で再導出した）。 -/
theorem layer_step_fails : (7 : ℕ) < 8 := by norm_num

theorem uniformizer_gain_ok : 2 * 3 * (8 - 6) ≤ 18 := by norm_num

theorem generic_gain_fails : ¬ (2 * 3 * (8 - 4) ≤ 18) := by norm_num

/-- ★★合成は無傷 —— `v(ε) − 予算 = 11 − 12 = −1 ≤ 3 = v(d(x,K))`
（★`d(x,K) = 3` も本波で再導出した）。 -/
theorem composite_has_slack : (11 : ℤ) - 12 ≤ 3 := by norm_num

/-- ★前波 `KrasnerCeiling.b2_obstruction_not_triggered` の数値条件
（要求半径 `8` ≤ 変位の最小 `v` `11`）。★`11` も本波で再導出した。 -/
theorem krasner_does_not_fire : (8 : ℕ) ≤ 11 := by norm_num

theorem sharp_jump_bound_holds : 2 * 8 ≤ 18 := by norm_num

end Numbers

/-! ## §3 使っている公理の一覧 -/

#print axioms sigma_order_nine
#print axioms tau_order_three
#print axioms sigma_cubed_eq_tau
#print axioms sigma_mod_three
#print axioms layer_step_fails
#print axioms uniformizer_gain_ok
#print axioms generic_gain_fails
#print axioms composite_has_slack
#print axioms krasner_does_not_fire
#print axioms sharp_jump_bound_holds

end CyclotomicVerified

end ABC3.Found.PGC
