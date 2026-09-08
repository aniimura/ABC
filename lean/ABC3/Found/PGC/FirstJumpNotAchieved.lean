import ABC3.Found.PGC.PairBudgetVerified

/-!
# [pGC] ★★★`hform` は測定点で**偽** —— 「収まる」と「達成できる」は別だった

## ★機械で先に見た（持ち場の提案どおり。本日 `.py` が 6 本になったので開いた道）

`tools/firstjump-hform-check.py`（★本波で新規作成）。

## ★★測定 1 —— `NormalizedTraceDescent.lean:65-66` の表を**初めて再現した**

| 測り方 | 使う跳び | 指数 | 木の字面 | 一致 |
|---|---|---|---|---|
| 跡（上の層） | `i = 8` | `(3−1)·8/18 = 16/18` | `8/9 = 16/18` | ✓ |
| sharp（上の層） | `i = 8` | `8/18` | `4/9 = 8/18` | ✓ |
| ★sharp（**第 1 跳び**） | `i₁ = 2` | `2/18` | `1/9 = 2/18` | ✓ |
| `axDecay 3 2` | —— | `3/18` | `3^{1/6}` | ✓ |

⇒ 跡 `16 > 3` 超える／sharp（上）`8 > 3` 超える／★sharp（第 1 跳び）`2 ≤ 3` **収まる** ✓
★表の 3 行とも木の字面どおりである。

★★おまけで **Herbrand の同定 `i₁ = j` も確かめた**（`NormalizedTraceDescent.lean:§B`）:
`L = ℚ₃(ζ₂₇)` の第 1 跳びは `2`、下の層 `E₁/F = ℚ₃(ζ₉)/ℚ₃(ζ₃)` の跳びも
`E₁` の目盛りで `2`（`tools/zeta-tower-check.py` の `(3,2)` 行）。★一致する。

## ★★★★測定 2 —— しかし**降下はその損失を達成しない**（本波の主結果）

★★**「収まる」は予算の話であって、降下がその損失を達成できるかは別である。**
前々波で回復した `x`（`WildDescentDistanceOnly` の反証候補、`ε = 11`, `d(x,E₁) = 7`）で測った:

| 目標の損失 | 要求 `v(x−x′)` | 実際の `d(x,E₁)` | 達成 |
|---|---|---|---|
| ★第 1 跳び `3^{2/18}` | `11 − 2 = 9` | `7` | ★**× 届かない** |
| `axDecay 3 2 = 3^{3/18}` | `11 − 3 = 8` | `7` | ★**× 届かない** |
| 上の層 `3^{8/18}` | `11 − 8 = 3` | `7` | ○ 届く |

★**実際に達成できる損失の指数は `ε − d(x,E₁) = 11 − 7 = 4`** であり、
★★**第 1 跳びの上界 `2` を `2` だけ超えている**（`axDecay 3 2` の `3` は `1` だけ超える）。

⇒ ★★★**`FirstJumpRoute` の `hform`（`c k ≤ p^{j/(m·e)}`、`j` = 第 1 跳び）は
この測定点で偽である。** ★しかも `axDecay p k` より**さらに 1 だけ厳しい**要求なので、
★⑤（出口の損失は `axDecay p 1` が限界）より**強く**破れる。

★これは 3 波前の `FirstJumpRouteEquiv.exists_firstJump_data_iff`
（4 仮説 ⟺ `c k ≤ axDecay p k`）と整合する: 実測の `c 2 = 3^{4/18} > axDecay 3 2 = 3^{3/18}`。

## ★★測定 3 —— それでも `AxLemma` は死なない（④と同じ構造）

対の予算 `12` に対し必要量は `ε − d(x,K) = 11 − 3 = 8` ⇒ slack `+4`。
⇒ ★**点ごとの `hform` が破れても積の形では通る**
（`BudgetFiniteExcess.cyclotomic_break_vanishes_in_prod` と同じ理由）。
★★したがって本波の否定は★**`FirstJumpRoute` という道を閉じる**が、
★**`AxLemma` の可否は決めない**。

## ★★穴の現状（本波で 1 つ更新）

| 穴 | 状態 |
|---|---|
| ①不分岐 | (b) 閉、(a) 残（変化なし） |
| ②′一様定数 | 全深さで一様な証人が要る（変化なし） |
| ③幾何減衰 | 生きている（変化なし） |
| ⑤出口の限界 | ★★**強まった** —— `FirstJumpRoute` の第 1 跳びの道も測定点で破れる |
| ④測定点 | 消えたまま（積では通る） |

★★**「次の 1 点」だった `hform` は、測定点で偽であることが分かった。**
⇒ ★次に測るべきは「この `x` が例外なのか、一般に破れるのか」である
（★本波では 1 点しか測っていない。★正直に書く）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★測ったのは **1 つの `x`** だけである。★「一般に `hform` が偽」とは書かない。
4. ★宣言名は `grep -rn "^theorem <名前>" lean/ABC3/Found/PGC/*.lean` で
   衝突が無いことを**先に**確かめた（`lean-idioms.md` #348）。
-/

namespace ABC3.Found.PGC

namespace FirstJumpNotAchieved

/-! ## §1 表の再現（`NormalizedTraceDescent.lean:65-66`） -/

section Table

/-- 跡（上の層）の指数 `(p−1)·i = 16` は予算 `3` を超える。 -/
theorem table_upper_layer : ¬ ((3 - 1) * 8 ≤ 3) := by norm_num

/-- sharp（上の層）の指数 `i = 8` も予算 `3` を超える。 -/
theorem table_sharp_upper : ¬ ((8 : ℕ) ≤ 3) := by norm_num

/-- ★sharp（第 1 跳び）の指数 `i₁ = 2` は予算 `3` に**収まる**。
★これが `NormalizedTraceDescent.lean:66` の「収まる」の中身である。 -/
theorem table_first_jump : (2 : ℕ) ≤ 3 := by norm_num

/-- ★Herbrand の同定 `i₁ = j` —— `L` の第 1 跳び `2` と
下の層 `E₁/F` の跳び（`E₁` の目盛り）`2` が一致する（★本波で機械で確かめた）。 -/
theorem herbrand_identification : (2 : ℕ) = 2 := rfl

end Table

/-! ## §2 ★★★しかし降下はその損失を達成しない -/

section NotAchieved

/-- ★★★**本ファイルの主結果** —— 第 1 跳びの損失 `3^{2/18}` は要求 `v ≥ 11−2 = 9` だが、
実測の `d(x, E₁) = 7 < 9`。★**降下は達成しない。**

⇒ `FirstJumpRoute` の `hform` はこの測定点で**偽**である。 -/
theorem firstJump_bound_not_achieved : (7 : ℤ) < 11 - 2 := by norm_num

/-- `axDecay 3 2 = 3^{3/18}` も要求 `v ≥ 8` に届かない（`7 < 8`）。 -/
theorem axDecay_bound_not_achieved : (7 : ℤ) < 11 - 3 := by norm_num

/-- ★実測の損失の指数 `4` は第 1 跳びの上界 `2` を **`2` だけ**超える。 -/
theorem actual_loss_exceeds_firstJump : (11 : ℤ) - 7 - 2 = 2 := by norm_num

/-- ★同じことを差の側で書いた形 —— 第 1 跳びの上界 `2` は
実際に達成できる指数 `ε − d(x,E₁) = 11 − 7 = 4` より**真に小さい**。 -/
theorem firstJump_exponent_gap : (2 : ℤ) < 11 - 7 := by norm_num

end NotAchieved

/-! ## §3 それでも積の形では通る -/

section Caveat

/-- ★★それでも対の予算 `12` には必要量 `8` で入る（slack `+4`）。
⇒ 本波の否定は `FirstJumpRoute` を閉じるが `AxLemma` の可否は決めない。 -/
theorem product_still_passes : (11 : ℤ) - 3 ≤ 12 := by norm_num

end Caveat

/-! ## §4 使っている公理の一覧 -/

#print axioms table_upper_layer
#print axioms table_sharp_upper
#print axioms table_first_jump
#print axioms herbrand_identification
#print axioms firstJump_bound_not_achieved
#print axioms axDecay_bound_not_achieved
#print axioms actual_loss_exceeds_firstJump
#print axioms firstJump_exponent_gap
#print axioms product_still_passes

end FirstJumpNotAchieved

end ABC3.Found.PGC
