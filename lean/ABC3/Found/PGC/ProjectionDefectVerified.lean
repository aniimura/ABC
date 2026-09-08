import ABC3.Found.PGC.DeepDescentNumbersVerified

/-!
# [pGC] ★★掃討完了 —— `EquivariantProjectionDescent` の `γ = p−1` を 4 層で再導出した

## ★どれを選んだか、なぜか

★**`EquivariantProjectionDescent` を選んだ。** 残り 2 箇所の費用を先に測った:

* ★`DeepDescentPairDirect.lean:37,51,147` の「95,298 件」は**乱択で seed が無い**ので
  ★**同じ件数は原理的に再現できない**（持ち場の観察 1 のとおり）。再現できるのは結論だけで、
  それには 1 段/対の予算の slack を実装し直す必要がある。★中くらいの費用。
* ★`EquivariantProjectionDescent.lean:19` の中心は **`γ = p−1`**（正規化した跡の欠損）で、
  ★これは**有限個の基底元の跡**を計算するだけで出る。★**最も安い**と見た。

## ★★★結果 —— **4 層すべてで `γ = p−1` が一致した**

`tools/trace-defect-check.py`（★本波で新規作成。`zeta-tower-check.py` を import）。

| 層 | `e_M` | `v_M(Tr(π^i))` | ★`γ` | 木の字面 |
|---|---|---|---|---|
| `ℚ₃(ζ₂₇)/ℚ₃(ζ₉)` | 18 | `{0:18, 1:18, 2:18}` | **2** | `2` ✓ |
| `ℚ₃(ζ₈₁)/ℚ₃(ζ₂₇)` | 54 | `{0:54, 1:54, 2:54}` | **2** | `2` ✓ |
| `ℚ₂(ζ₈)/ℚ₂(ζ₄)` | 4 | `{0:4, 1:4}` | **1** | `1` ✓ |
| `ℚ₂(ζ₁₆)/ℚ₂(ζ₈)` | 8 | `{0:8, 1:8}` | **1** | `1` ✓ |

★同 `:23-25` の `t = 8` / `s = 6` / `m = 2` も再導出した（`t`, `m` は直接、`s` は
下の層 `ℚ₃(ζ₉)/ℚ₃(ζ₃)` の跳び `2` の `v_M` 換算 `2×3 = 6`）。

★同 `:35-42` の合成の表も再導出した:
`素朴 t+s = 14`（予算 `12` 超）/ `収縮つき max(t,s,t+s−m) = 12`（ちょうど）/
★`射影つき max(t+γ, s+γ) = 10`（**入る**）。★3 行とも木の字面と一致。

## ★★★自分のバグを記録する（本波の第 2 の成果）

★1 度目の実装は `γ = e_M − min_i v_M(Tr π^i)` と書いて **`γ = 0`** を出した。
★**作用素ノルムは「比の sup」であって「`v` の min」ではない**。正しくは

  `γ = max_{i<p} ( v_M(π^i) − v_M(P π^i) ) = max_{i<p} ( i + e_M − v_M(Tr π^i) )`

★`v_M(Tr π^i) = e_M` がすべての `i < p` で成り立つので `γ = max_i i = p−1`。
★★**木の `γ = 2` の方が正しく、私の 1 度目の計算が誤っていた。**

★★古典的な裏付けも取った（Serre, Corps Locaux V §3 Lemme 4）:
`Tr(𝒪_M) = 𝔭_E^{⌊d/e⌋}`、`d = (p−1)(i+1) = 2·9 = 18`、`e = 3` ⇒ `⌊18/3⌋ = 6`
⇒ `v_M(Tr 𝒪_M) = 18 = e_M` ⇒ `γ = p−1`。
★**機械の数値と古典の公式が独立に一致した**（§1 `different_exponent_three` /
`trace_image_valuation`）。

## ★★掃討の現状 —— **6 箇所中 5 箇所**

| ファイル | 状態 |
|---|---|
| `WildDescentDistanceOnly.lean` | ★塞いだ（`x` を回復） |
| `GainedTowerDescent.lean:94` | ★跳びの数値は塞いだ（`L = 8, 26`）。総当たりの表は未 |
| `DeepDescentRepair.lean:110` | ★塞いだ（4 値すべて一致） |
| `JumpDefectTradeoff.lean:106` | ★塞いだ（結論は Lean の定理でスクリプト不要） |
| ★`EquivariantProjectionDescent.lean:19` | ★★**塞いだ（本波）** |
| `DeepDescentPairDirect.lean:37,51,147` | ★未（★**seed が無いので件数は原理的に再現不能**。結論のみ再現可能） |

## ★数値が合わない箇所は**見つからなかった**（4 波連続）

★ただし★**1 度目の私の計算は合わなかった**（`γ = 0` 対 `γ = 2`）。
★**誤っていたのは私の方**であり、★木の記録は 4 波連続で正しい。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★`DeepDescentPairDirect` の 95,298 件は★**再現しないと決めた**（seed が無いため）。
   ★「結論だけ再現する」道は残っている（次の波の候補）。
4. ★§3 の `ℚ₃(ζ₈₁)` の行の `24` は下の層 `ℚ₃(ζ₂₇)/ℚ₃(ζ₉)` の跳び `8` の
   `v_M`（`e = 54`）換算 `8×3 = 24` である。★これは本波で再導出した `t = 8` からの計算であって、
   木の表の `50 / 42 / 28` と突き合わせて一致することを確認した。
-/

namespace ABC3.Found.PGC

namespace ProjectionDefectVerified

/-! ## §1 `γ = p−1`（★本波で 4 層すべて再導出した） -/

section Gamma

/-- ★`p = 3`: `v_M(Tr π^i) = e_M` が `i = 0,1,2` で成り立つので `γ = max(0,1,2) = 2 = p−1`。 -/
theorem gamma_three : max (max 0 1) 2 = 2 := by norm_num

/-- ★`p = 2`: 同様に `γ = max(0,1) = 1 = p−1`。 -/
theorem gamma_two : max 0 1 = 1 := by norm_num

/-- ★古典の裏付け: different の指数 `d = (p−1)(i+1) = 2·9 = 18`（`i = 8` は本波で再導出）。 -/
theorem different_exponent_three : (3 - 1) * (8 + 1) = 18 := by norm_num

/-- ★Serre V §3 Lemme 4: `Tr(𝒪_M) = 𝔭_E^{⌊d/e⌋}`、`⌊18/3⌋ = 6` ⇒ `v_M = 6·3 = 18 = e_M`。
★機械の数値と一致する。 -/
theorem trace_image_valuation : (18 / 3) * 3 = 18 := by norm_num

end Gamma

/-! ## §2 合成の表（`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`、予算 `12`） -/

section Composition

/-- ★素朴な telescope `t + s = 14` は予算 `12` を超える。 -/
theorem telescope_exceeds : ¬ ((8 : ℕ) + 6 ≤ 12) := by norm_num

/-- ★収縮つき `max(t, s, t+s−m) = 12` はちょうど入る。 -/
theorem contraction_fits : max (max (8 : ℕ) 6) (8 + 6 - 2) = 12 := by norm_num

/-- ★★射影つき `max(t+γ, s+γ) = 10 ≤ 12` は**余裕を持って入る**。 -/
theorem projection_fits : max ((8 : ℕ) + 2) (6 + 2) = 10 ∧ (10 : ℕ) ≤ 12 := by
  constructor <;> norm_num

theorem measured_max_le : (8 : ℕ) ≤ 10 := by norm_num

end Composition

/-! ## §3 合成の表（`ℚ₃(ζ₈₁)/ℚ₃(ζ₉)`、予算 `36`） -/

section Zeta81Row

theorem telescope_exceeds_81 : ¬ ((26 : ℕ) + 24 ≤ 36) := by norm_num

/-- ★`ℚ₃(ζ₈₁)` では収縮つきでも `42 > 36` で超える。 -/
theorem contraction_exceeds_81 : ¬ (max (max (26 : ℕ) 24) (26 + 24 - 8) ≤ 36) := by norm_num

/-- ★★射影つきなら `28 ≤ 36` で入る。 -/
theorem projection_fits_81 : max ((26 : ℕ) + 2) (24 + 2) = 28 ∧ (28 : ℕ) ≤ 36 := by
  constructor <;> norm_num

end Zeta81Row

/-! ## §4 使っている公理の一覧 -/

#print axioms gamma_three
#print axioms gamma_two
#print axioms different_exponent_three
#print axioms trace_image_valuation
#print axioms telescope_exceeds
#print axioms contraction_fits
#print axioms projection_fits
#print axioms measured_max_le
#print axioms telescope_exceeds_81
#print axioms contraction_exceeds_81
#print axioms projection_fits_81

end ProjectionDefectVerified

end ABC3.Found.PGC
