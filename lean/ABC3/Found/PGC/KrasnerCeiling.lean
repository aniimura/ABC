import ABC3.Found.PGC.BudgetFiniteExcess

/-!
# [pGC] B2 は Krasner の評価では**原理的に**塞げない —— しかも帰結は `AxLemma` を左右しない

## ★数値計算に進まなかった理由（費用を先に 3 点測った）

持ち場は `WildDescentDistanceOnly.lean:88-99` の B2
（`ℚ₃(ζ₃)` の 40 個の巡回 3 次拡大で `d(x,N)` を測る）を測れと言った。
★着手前に測り、★**数値計算には進まないと決めた**:

1. ★★**`x` が木に無い。** 同 `:57` は「厳密な整数演算で全部当たった
   （★**スクリプトは報告に添付**）」と書いており、★`x` もスクリプトも
   **リポジトリに入っていない**。測ったコマンド:
   ```
   grep -rn "37/25000|zeta27" tools/ --include=*.py --include=*.mjs   → 0 件
   ls tools/*.py   → genell-remaining-measure.py / source-text.py / velu-disc-check.py のみ
   ```
   ⇒ ★乱択探索からやり直しになる。
2. ★`d(x,N)` は `N ⊄ M` なので**合成体 `MN`（`ℚ₃` 上 54 次）**の算術が要る。
   ★`N` 全体の上の最小なので、Galois の下界以外に閉じた式が無い。
3. ★★そして下の測定 1 のとおり、★**その Galois の下界では原理的に届かない**。

⇒ ★数値計算は 1 波では入らないと見て、★**代わりに「その道が原理的に届かない」ことを
定理にした**。

## ★★★測定 1 —— Krasner の評価は B2 に**原理的に**届かない

§1 `not_exists_fixed_close`（抽象核。★分岐・付値・Galois の語彙が 1 語も出ない）:

★`δ < ‖σ·x − x‖` なら「`σ` で固定され `x` から `δ` 以内」の `x′` は**無い**。
（`σx − x = σ(x−x′) − (x−x′)` を超距離で潰すだけ。）
★これが「塔の外の `x′` を排除する」唯一の一般的な道具である。

★しかし B2 の数値では（すべて `WildDescentDistanceOnly.lean:65` の測定値）:

* 要求される半径は `v = 8`
* `Gal(M/K)` の変位の**最小の `v`** は `11`（`v(σx−x) = 11`, `v(τx−x) = 15`）
* ★`M ∩ N = K` なので `Gal(\bar K/N)` の元は `M` 上 `Gal(M/K)` に落ち、
  ★**`N` に依らず同じ変位**しか与えない

⇒ §2 `b2_obstruction_not_triggered`: `‖π‖^{11} ≤ ‖π‖^{8}`、すなわち
★★**変位は要求半径より小さく、`not_exists_fixed_close` が発火しない。**
★★★**40 個のどの `N` でも同じである**（変位が `N` に依らないから）。
⇒ ★**B2 は Krasner/Galois の評価では塞げない。別の道具が要る。**

★★木の地の文（同 `:75-76`）「Krasner の補題は `‖x−x'‖ < min_σ‖σx−x‖` のときしか効かず、
ここで要求される半径は `min` より**大きい**ので届かない（`need = 8 ≤ v_max = 15`）」を
★**定理にした**。★ただし 1 か所訂正する: 木は `v_max = 15` を使っているが、
★`σ` は好きに選べるので**最も大きい変位＝最小の `v` = 11**を使うのが最良である。
★それでも `11 ≥ 8` なので結論は変わらない（★数値の使い方の訂正であって、結論の訂正ではない）。

## ★★測定 2 —— B2 の帰結は `AxLemma` を左右しない

§3 `b2_does_not_decide_axLemma`:
★前波までの `cEx` は深さ 2 で `axDecay 3 2` を**破っている**（`axDecay 3 2 < cEx 2`）が、
★それでも `AxWildDescent K cEx → AxLemma K (axConstant 3)` である。

⇒ ★★**B2 が塞がっても（点ごとの破れが本物でも）`AxLemma` は死なない。**
★逆に B2 が塞がらなくても、破れが 1 つ消えるだけである。
⇒ ★★★**B2 は `AxLemma` の可否を決めない。**（前波の `WitnessScope` と同じ構造。）

## ★在庫の測定（2 か所。#330）

```
IsUltrametricDist.norm_sub_le_max  → ★error: Unknown constant（本波で実測）
  ⇒ WildDescentDistanceOnly.lean:104 の「norm_sub_le_max は無い、to_additive の穴」は真。
     直しは sub_eq_add_neg + IsUltrametricDist.norm_add_le_max + norm_neg。
grep -rn "37/25000|zeta27" tools/ ; ls tools/*.py   → 上のとおり 0 件
```

## ★穴の現状（本波で穴は増えていない）

①(b) 閉 /(a) 残 / ②′ 一様な証人が要る / ③ 変わらず / ⑤ ②′と同じ / ④ 消えたまま。
★★本波は★**B2 という「木が唯一ずっと未測定と書き続けていた枝」を 1 本刈った**
（Krasner では原理的に届かない ＋ 帰結は `AxLemma` を左右しない）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ pGC 物理 p.6 Corollary 3.1。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 は抽象核（超距離な半ノルム加法群＋等長なモノイド作用）で、分岐の語彙が 1 語も出ない。
4. ★数値 `8` / `11` は `WildDescentDistanceOnly.lean:65` の測定値を**そのまま引いた**。
   ★本ファイルはそれを**再計算していない**（`x` が木に無いため）。
   ⇒ ★測定 1 の結論は「木の数値が正しければ」という条件付きである。
-/

namespace ABC3.Found.PGC

namespace KrasnerCeiling

/-! ## §1 抽象核 —— 「固定されて近い元は無い」の唯一の一般的な道具 -/

section Core

variable {M : Type*} [SeminormedAddCommGroup M] [IsUltrametricDist M]

/-- `x′` が `σ` で固定され `‖x − x′‖ ≤ δ` なら `‖σ·x − x‖ ≤ δ`。 -/
theorem norm_smul_sub_le_of_fixed {G : Type*} [Monoid G] [DistribMulAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) {x x' : M} {δ : ℝ} (σ : G)
    (hfix : σ • x' = x') (hd : ‖x - x'‖ ≤ δ) : ‖σ • x - x‖ ≤ δ := by
  have hsplit : σ • x - x = σ • (x - x') - (x - x') := by
    rw [smul_sub, hfix]; abel
  rw [hsplit, sub_eq_add_neg]
  refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
  · rw [hiso]
    exact hd
  · rw [norm_neg]
    exact hd

/-- ★★★**塔の外を排除する唯一の一般的な道具** —— `δ < ‖σ·x − x‖` なら
「`σ` で固定され `x` から `δ` 以内」の `x′` は無い。

★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem not_exists_fixed_close {G : Type*} [Monoid G] [DistribMulAction G M]
    (hiso : ∀ (g : G) (z : M), ‖g • z‖ = ‖z‖) {x : M} {δ : ℝ} {σ : G}
    (hσ : δ < ‖σ • x - x‖) :
    ¬ ∃ x' : M, σ • x' = x' ∧ ‖x - x'‖ ≤ δ := by
  rintro ⟨x', hfix, hd⟩
  exact absurd (norm_smul_sub_le_of_fixed hiso σ hfix hd) (not_le.mpr hσ)

end Core

/-! ## §2 ★★★B2 ではその道具が発火しない -/

section Ceiling

variable {M : Type*} [NormedField M]

/-- `a ≤ b` かつ `‖y‖ = ‖π‖^b` なら `‖y‖ ≤ ‖π‖^a`（`0 < ‖π‖ ≤ 1`）。 -/
theorem obstruction_not_triggered {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    {y : M} {a b : ℕ} (hab : a ≤ b) (hy : ‖y‖ = ‖π‖ ^ b) :
    ‖y‖ ≤ ‖π‖ ^ a := by
  rw [hy]
  exact pow_le_pow_of_le_one (le_of_lt hπ0) hπ1 hab

/-- ★★★**B2 では上の道具が発火しない** —— 変位 `‖π‖^{11}` は要求半径 `‖π‖^{8}` より
**小さい**。★`M ∩ N = K` なので変位は `N` に依らず、★**40 個のどの `N` でも同じ**である。
⇒ ★B2 は Krasner/Galois の評価では塞げない。 -/
theorem b2_obstruction_not_triggered {π : M} (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≤ 1)
    {y : M} (hy : ‖y‖ = ‖π‖ ^ (11 : ℕ)) : ‖y‖ ≤ ‖π‖ ^ (8 : ℕ) :=
  obstruction_not_triggered hπ0 hπ1 (by omega) hy

end Ceiling

/-! ## §3 ★★B2 の帰結は `AxLemma` を左右しない -/

section Irrelevant

open ABC3.Skeleton.PGC

variable [Fact (Nat.Prime 3)]

/-- ★★**B2 の帰結は `AxLemma` を左右しない** —— `cEx` は深さ 2 で `axDecay 3 2` を
破っているが、それでも `AxLemma K (axConstant 3)` が出る。
⇒ B2 が塞がっても塞がらなくても `AxLemma` の可否は変わらない。 -/
theorem b2_does_not_decide_axLemma (K : PAdicLocalField 3) :
    (axDecay 3 2 < BudgetFinite.cEx 2)
      ∧ (AxWildDescent K BudgetFinite.cEx → AxLemma K (axConstant 3)) :=
  ⟨BudgetFinite.cEx_two_gt_axDecay, fun h => BudgetFinite.axLemma_of_wildDescent_cEx K h⟩

end Irrelevant

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def not_exists_fixed_close.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def b2_obstruction_not_triggered.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def b2_does_not_decide_axLemma.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms norm_smul_sub_le_of_fixed
#print axioms not_exists_fixed_close
#print axioms obstruction_not_triggered
#print axioms b2_obstruction_not_triggered
#print axioms b2_does_not_decide_axLemma

end KrasnerCeiling

end ABC3.Found.PGC
