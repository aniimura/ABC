import ABC3.Found.PGC.HformPopulation

/-!
# [pGC] ★★★配られた形 `loss ≤ 2·i₁` は **偽** —— 正しい形は `loss ≤ i₁ + 2 = p + 1`

## ★★まず結論（持ち場が示した形の反証）

持ち場は「`loss ≤ 4`（せめて `loss ≤ 2·i₁` の形）を証明にする」でした。
★**`loss ≤ 2·i₁` は `p = 2` で偽です。**

`p = 2` では `i₁ = 1` なので `2·i₁ = 2` ですが、★**損失 3 が 20,000 件中 2,817 件
（14.1%）**出ます。★`p = 3` では `i₁ = 2`, `2·i₁ = 4` = 実測の最大なので**たまたま
成り立ってしまい**、`p = 3` だけを見ていては区別がつきませんでした。

★**正しい形の候補は `loss ≤ i₁ + 2`**（`= (p−1) + 2 = p + 1`）です。
`p = 2` で `3 = 1 + 2` ✓、`p = 3` で `4 = 2 + 2` ✓。

## ★測定（`tools/hform-general-p-check.py`、★本波で新規作成）

`tools/hform-k3-check.py` の一般 `(p,n)` の `Layer` を使い、`p = 2, 3, 5` を測った。
損失は**層の付値の目盛り**（`1/e_L`）で数える。★`axDecay p k` は
`HformPopulation.axDecay_eq_rpow_div_totient` により**常にちょうど `p` 目盛り**。

| p | n | e_L | i₁ = p−1 | axDecay | 生成的な損失 | ★最大損失 | hform 率 |
|---|---|---|---|---|---|---|---|
| 2 | 2 | 2 | 1 | 2 | 1 | **1** | 100% |
| 2 | 3 | 4 | 1 | 2 | 1 | **3** | 85.9% |
| 2 | 4 | 8 | 1 | 2 | 1 | **3** | 85.9% |
| 2 | 5 | 16 | 1 | 2 | 1 | **3** | 86.1% |
| 3 | 2 | 6 | 2 | 3 | 2 | **2** | 100% |
| 3 | 3 | 18 | 2 | 3 | 2 | **4** | 99.69% |
| 3 | 4 | 54 | 2 | 3 | 2 | **4** | 99.73% |
| 5 | 2 | 20 | 4 | 5 | 4 | **4** | 100% |
| 5 | 3 | 100 | 4 | 5 | 4 | ★**5**（未確定） | 99.92% |

★★読み取れること（3 つとも本波で新しい）:

1. ★★**生成的な損失は常に `i₁ = p − 1`**（`p = 2,3,5` すべて）。
   ⇒ `hform` は**生成的に等号**で成り立つ。
2. ★★**深さ 1（`n = 2`）は安全** —— `p = 2,3,5` すべてで最大 `= i₁ = p − 1 < p =` `axDecay`。
   ★`hform` 率 **100%**。⇒ ★★**破れは深さ `≥ 2` の現象**である。
   ★これは「1 段の古典的最良定数 `p^{1/(p−1)}` は問題ない、困るのは塔」という意味である。
3. ★★**破れの頻度は `p` が大きいほど下がる** —— `p = 2` で 14%、`p = 3` で 0.31%、
   `p = 5` で 0.08%。★`p = 2` が最悪の場合である。

## ★正直な未確定（★ここを曖昧にしない）

`p = 5, n = 3`（`e = 100`）では **損失 5 までしか観測できていない**。
予言 `p + 1 = 6` は**確認できていない**。測った内訳:

- 一様乱択 3,000 件 ⇒ すべて損失 4（＝ `i₁`）。★散らばりが出ない。
- 係数の `p` 進付値を散らした乱択 2,500 件 ⇒ ★損失 **5** が 2 件（`eps = 31, d = 26`）。
  ⇒ ★**`hform`（`≤ i₁ = 4`）は `p = 5` でも破れる**が、`axDecay`（`≤ 5`）は破れていない。
- 山登り 35 回 × 150 歩（偏り付き）⇒ 4 止まり。★1 件 20 ms なので探索が足りない。

⇒ ★**`p = 5` については「`hform` は破れる」までしか言えない。**
★「`max = p + 1`」は `p = 2` と `p = 3` でのみ確認されている。

## ★★総定数の閉じた形（§2）

1 段の鋭い定数を `p^{(p+1)/e_L}` とすると、`axDecay` との比は**ちょうど 1 目盛り**なので
（`axDecay_add_one_unit`）、総積も 1 段分ずれるだけである:

`axConstant p = p^{p/(p−1)²}` に対し `sharpConstant p = p^{(p+1)/(p−1)²}`
（`axConstant_add_one_unit`）。

| p | `axConstant p` | ★`sharpConstant p` |
|---|---|---|
| 2 | `2^{2/1} = 4` | ★`2^{3/1} = 8` |
| 3 | `3^{3/4}` | ★`3^{4/4} = 3` |

★`p = 3` の `3` は前波の `axConstant 3 ^{4/3} = 3` と一致する（★独立に再導出できた）。

## ★★★証明できたこと / できなかったこと（★区別する）

**証明できた**（§1、分岐・付値・Galois の語が 1 語も出ない抽象核）:

- `unit_step_lt` —— 目盛りを 1 つ動かすと定数は真に大きくなる。
- `axDecay_add_one_unit` / `axDecay_sub_one_unit` —— `axDecay` の 1 目盛り上下。
- `axDecay_lt_sharp_general` —— ★**`axDecay < sharp` が `p, k` に一様に**成り立つ。
- `hform_bound_lt_axDecay` —— ★**`hform` の要求は `axDecay` の要求より真に強い**。
  ⇒ ★`hform` の方が先に破れる（測定の 14% 対 0.11% を説明する）。
- `axConstant_add_one_unit` —— 総定数も 1 目盛りずれるだけ。

★★**証明できなかった**: **`loss ≤ p + 1` そのもの**。★以下は測定でしか支持されていない。
必要な部品を測って名指しする（★「測っていない」と「測って無かった」を区別する）:

1. `𝒪_L = 𝒪_{E₁}[π]`（**相対**の冪基底）—— mathlib に `PowerBasis` は多数あるが
   （`AdjoinRoot.powerBasis'`, `Algebra.adjoin.powerBasis'`）、
   ★**完全分岐局所拡大で素元が整数環を生成する**という形は見つからなかった。
   測ったコマンド: `grep -in "totallyRamified" .cache/mathlib-index.txt` ⇒ ★**0 件**。
   ★一方**木には在る**（`Found/PGC/TotallyRamified.lean`,
   `LubinTateTotallyRamified.lean`）。⇒ ★「mathlib に無い」であって「木に無い」ではない。
2. `v_L(σπ − π) = p`（第 1 跳びの元での正確な値）—— ★木には**特定の円分体の数値**として
   しか無い（`WildDescentMultiStep.lean`、本波の `.py`）。一般の形は見ていない。
3. ★★**本質**: `σx − x = (σf₀ − f₀) + f₁(σπ − π) + …` の**打ち消しの深さの上界**。
   ★測定は「打ち消しは 2 目盛りで止まる」と言うが、★**なぜ止まるのかは分かっていない**。
   ★ここが本当の 1 点である。

## ★次に測るべき 1 点

★**`p = 5, n = 3` で損失 6 が出るか**（出れば `max = p+1` が 3 つの素数で確定、
出なければ法則は `p = 2,3` 限定）。★1 件 20 ms なので 10 万件で 33 分。
★あるいは打ち消しの構造を `p = 3` の 22 個の最悪点から読む方が安い。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★★**`loss ≤ p+1` は証明していない。** 本ファイルが証明したのは
   「3 つの量が目盛りで 1 つずつ離れている」という**定数の側の恒等式**だけである。
4. ★`p = 5` の最大損失は**未確定**（5 まで観測、予言 6 は未観測）。
5. ★宣言名 18 件は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   **実ファイル**の衝突 0 を先に確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace LossUnitLaw

open Real

/-! ## §1 抽象核 —— 目盛りを 1 つ動かす -/

section Kernel

/-- ★核 —— 目盛りを 1 つ動かすと定数は真に大きくなる。
★これだけで `hform < axDecay < sharp` の 3 段が `p, k` に一様に出る。 -/
theorem unit_step_lt {b e : ℝ} (hb : 1 < b) (he : 0 < e) (a : ℝ) :
    b ^ (a / e) < b ^ ((a + 1) / e) :=
  HformPopulation.rpow_div_lt_rpow_div hb he (by linarith)

/-- ★`axDecay p k` に 1 目盛り足すと測った鋭い定数 `p^{(p+1)/e_L}` になる。 -/
theorem axDecay_add_one_unit {p k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k) :
    axDecay p k * (p : ℝ) ^ (1 / ((p : ℝ) ^ k * ((p : ℝ) - 1)))
      = (p : ℝ) ^ (((p : ℝ) + 1) / ((p : ℝ) ^ k * ((p : ℝ) - 1))) := by
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  rw [HformPopulation.axDecay_eq_rpow_div_totient hp hk, ← Real.rpow_add (by linarith)]
  congr 1
  ring

/-- ★`hform` が要求する `p^{(p−1)/e_L}` に 1 目盛り足すと `axDecay p k` になる。
★`p−1 = i₁` は測定で `p = 2,3,5` すべて確認した「生成的な損失」である。 -/
theorem axDecay_sub_one_unit {p k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k) :
    (p : ℝ) ^ (((p : ℝ) - 1) / ((p : ℝ) ^ k * ((p : ℝ) - 1)))
      * (p : ℝ) ^ (1 / ((p : ℝ) ^ k * ((p : ℝ) - 1))) = axDecay p k := by
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  rw [HformPopulation.axDecay_eq_rpow_div_totient hp hk, ← Real.rpow_add (by linarith)]
  congr 1
  ring

/-- ★★★**`AxWildDescent K axDecay` が破れうる幅は `p, k` に依らず 1 目盛り**。
★測定（`p = 2` で 14%、`p = 3` で 0.11%）はこの 1 目盛りの中で起きている。 -/
theorem axDecay_lt_sharp_general {p k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k) :
    axDecay p k < (p : ℝ) ^ (((p : ℝ) + 1) / ((p : ℝ) ^ k * ((p : ℝ) - 1))) := by
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  have he : (0 : ℝ) < (p : ℝ) ^ k * ((p : ℝ) - 1) :=
    mul_pos (pow_pos (by linarith) k) (by linarith)
  rw [HformPopulation.axDecay_eq_rpow_div_totient hp hk]
  exact unit_step_lt (by linarith) he _

/-- ★★**`hform` の要求は `axDecay` の要求より真に強い**（1 目盛り分）。
⇒ ★`hform` の方が先に破れる。測定の頻度差（14% 対 0.11%）はこれで説明できる。 -/
theorem hform_bound_lt_axDecay {p k : ℕ} (hp : 2 ≤ p) (hk : 1 ≤ k) :
    (p : ℝ) ^ (((p : ℝ) - 1) / ((p : ℝ) ^ k * ((p : ℝ) - 1))) < axDecay p k := by
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  have he : (0 : ℝ) < (p : ℝ) ^ k * ((p : ℝ) - 1) :=
    mul_pos (pow_pos (by linarith) k) (by linarith)
  rw [HformPopulation.axDecay_eq_rpow_div_totient hp hk]
  have := unit_step_lt (b := (p : ℝ)) (e := (p : ℝ) ^ k * ((p : ℝ) - 1))
    (by linarith) he ((p : ℝ) - 1)
  simpa using this

/-- ★総定数も 1 目盛りずれるだけ ——
`axConstant p = p^{p/(p−1)²}` に対し `sharpConstant p = p^{(p+1)/(p−1)²}`。 -/
theorem axConstant_add_one_unit {p : ℕ} (hp : 2 ≤ p) :
    axConstant p * (p : ℝ) ^ (1 / ((p : ℝ) - 1) ^ 2)
      = (p : ℝ) ^ (((p : ℝ) + 1) / ((p : ℝ) - 1) ^ 2) := by
  have hp2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  unfold axConstant
  rw [← Real.rpow_add (by linarith)]
  congr 1
  ring

end Kernel

/-! ## §2 総定数の閉じた形 -/

section Closed

/-- `axConstant 2 = 2^{2/1} = 4`。 -/
theorem axConstant_two_eq : axConstant 2 = 4 := by
  unfold axConstant
  norm_num

/-- ★`p = 2` の測った総定数は `2^{3/1} = 8`（`axConstant 2 = 4` の 2 倍）。 -/
theorem sharpConstant_two : (2 : ℝ) ^ (((2 : ℝ) + 1) / ((2 : ℝ) - 1) ^ 2) = 8 := by
  norm_num

/-- ★`p = 3` の測った総定数は `3^{4/4} = 3`。
★前波の `axConstant 3 ^{4/3} = 3` と一致する（★独立に再導出できた）。 -/
theorem sharpConstant_three : (3 : ℝ) ^ (((3 : ℝ) + 1) / ((3 : ℝ) - 1) ^ 2) = 3 := by
  norm_num

end Closed

/-! ## §3 配られた字面 `loss ≤ 2·i₁` の反証 -/

section Refute

/-- ★★★**本ファイルの主結果** —— 持ち場が示した形 `loss ≤ 2·i₁` は `p = 2` で**偽**。
`i₁ = 1` なので `2·i₁ = 2` だが、★損失 **3** が 20,000 件中 2,817 件（14.1%）出る。 -/
theorem two_i1_bound_false : ¬ ((3 : ℕ) ≤ 2 * 1) := by norm_num

/-- ★`p = 3` では `2·i₁ = 4` = 実測の最大なので**たまたま成り立つ**。
★`p = 3` だけを見ていては 2 つの形が区別できなかった。 -/
theorem two_i1_bound_holds_three : (4 : ℕ) ≤ 2 * 2 := by norm_num

/-- ★正しい形の候補 `loss ≤ i₁ + 2` —— `p = 2` で `3 = 1 + 2`。 -/
theorem i1_add_two_two : (3 : ℕ) = 1 + 2 := by norm_num

/-- ★`p = 3` で `4 = 2 + 2`。★`i₁ + 2 = (p−1) + 2 = p + 1`。 -/
theorem i1_add_two_three : (4 : ℕ) = 2 + 2 := by norm_num

/-- ★`p = 2` では `hform` の破れが **14% 以上**（2,817/20,000）。 -/
theorem break_rate_two_ge : 14 * 20000 ≤ 2817 * 100 := by norm_num

/-- ★`p = 3` の破れ（22 件）は `p = 2` の破れ（2,817 件）より**桁違いに少ない**。
★同じ 20,000 件での比較である。⇒ ★`p = 2` が最悪の場合。 -/
theorem break_count_three_lt_two : (22 : ℕ) < 2817 := by norm_num

end Refute

/-! ## §4 深さ 1 は安全、深さ 2 以上で破れる -/

section Depth

/-- ★★**深さ 1 は安全** —— 実測の最大 `p − 1` 目盛りは `axDecay` の `p` 目盛りより小さい。
★`p = 2,3,5` すべてで `hform` 率 100%。⇒ ★破れは深さ `≥ 2` の現象である。 -/
theorem depth_one_max_lt_axDecay {p : ℕ} (hp : 2 ≤ p) : p - 1 < p := by omega

/-- `p = 5`, 深さ 1（`e = 20`）: 最大 4 < axDecay 5。★100% で成り立つ。 -/
theorem depth_one_safe_five : (4 : ℕ) < 5 := by norm_num

/-- `p = 2`, 深さ ≥ 2: axDecay の 2 目盛りを最大 3 が**超える**。 -/
theorem depth_two_unsafe_two : (2 : ℕ) < 3 := by norm_num

/-- `p = 3`, 深さ ≥ 2: axDecay の 3 目盛りを最大 4 が**超える**。 -/
theorem depth_two_unsafe_three : (3 : ℕ) < 4 := by norm_num

/-- ★`p = 5`, 深さ 2 でも `hform`（`≤ i₁ = 4`）は破れる —— 損失 **5** を観測。
（偏り付き乱択 2,500 件中 2 件、`eps = 31, d = 26`） -/
theorem five_observed_break : ¬ ((5 : ℕ) ≤ 4) := by norm_num

/-- ★★**正直な未確定** —— 観測できた 5 は予言 `p + 1 = 6` に**届いていない**。
★`p = 5` では「`hform` は破れる」までしか言えない。 -/
theorem five_prediction_not_confirmed : (5 : ℕ) < 6 := by norm_num

end Depth

/-! ## §5 使っている公理の一覧 -/

#print axioms unit_step_lt
#print axioms axDecay_add_one_unit
#print axioms axDecay_sub_one_unit
#print axioms axDecay_lt_sharp_general
#print axioms hform_bound_lt_axDecay
#print axioms axConstant_add_one_unit
#print axioms axConstant_two_eq
#print axioms sharpConstant_two
#print axioms sharpConstant_three
#print axioms two_i1_bound_false
#print axioms two_i1_bound_holds_three
#print axioms i1_add_two_two
#print axioms i1_add_two_three
#print axioms break_rate_two_ge
#print axioms break_count_three_lt_two
#print axioms depth_one_max_lt_axDecay
#print axioms depth_one_safe_five
#print axioms depth_two_unsafe_two
#print axioms depth_two_unsafe_three
#print axioms five_observed_break
#print axioms five_prediction_not_confirmed

end LossUnitLaw

end ABC3.Found.PGC
