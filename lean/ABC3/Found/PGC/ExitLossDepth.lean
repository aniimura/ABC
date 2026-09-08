import ABC3.Found.PGC.BudgetFiniteExcess

/-!
# [pGC] 本日の出口は `cEx` を実現しない —— 1 段の損失は深さに依らず `axDecay p 1` 以上

## ★配られた 2 つの問いへの答え

**(1) 本日の出口が出す 1 段の損失は深さに依るか。★依らない。しかも `1` に落ちない。**

* `TotallyRamifiedLayer.lean:305 exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p`
  （`k = 0`、次数 `p` の層 1 枚）の損失は ★**`axDecay p 1`（定数）**。
* 一般の `k` の出口 `TotallyRamifiedLayer.lean:266
  exists_norm_sub_algebraMap_le_prod_axDecay_of_uniformizer` は `∃ y : twr g p k` を返すが、
  ★`twr g p j` は `⟨g^{p^j}⟩` の固定体（`GainedTowerModel.lean:293`）なので
  `twr g p k` は `[M : twr g p k] = p` の★**1 段下**である。
  ★その 1 段に **`∏_{j ∈ Icc 1 (k+1)} axDecay p j`（＝塔ぜんぶの予算）**を払う。
  ⇒ ★**深さが増えるほど 1 段の損失は大きくなる**（各因子 > 1）。

★どちらの読み方でも `exitLoss p d := ∏_{j ∈ Icc 1 d} axDecay p j ≥ axDecay p 1 > 1`
（§2 `axDecay_one_le_exitLoss`）で、★**`1` に落ちる段が 1 つも無い**（`exitLoss_ne_one`）。

⇒ ★★★**前波の `cEx`（`d ≥ 3` で `1`）は本日の出口では実現されない**
（§4 `exit_does_not_realize_cEx` / `degP_exit_does_not_realize_cEx`）。
`BudgetFiniteExcess.budget_of_eventually_le_one` /
`BudgetFinite.prod_Icc_le_of_eventually_one` の逃げ道（★悪い段が有限個）は**使えない**。
積は非有界である（§1 `not_forall_prod_Icc_le_of_le`）。

**(2) `c` は深さが大きいほど良くなる必要があるか。★ある。`axDecay p k` がすでにその形である。**

`axDecay p k = p^{(1/(p−1))p^{1−k}}` は `k → ∞` で **1 に収束**する
（前波 `PGroupDescentToAxWild.exists_axDecay_lt`）。
★一方 `exitLoss p d` は `d` について**増加**する。⇒ ★**向きが逆である。**

## ★★★どこで止まるか —— 「最後の跳び」では原理的に足りない

§3 `last_jump_bound_not_enough`: `d ≥ 2` のとき、`0 < e` と `(p−1)·i ≤ e` を満たす
`(i, e)`（例: `i = 1`, `e = p−1`）で ★`axDecay p d < p^{i/e}` となるものが**存在する**。
⇒ ★★**`RamificationJumpBound.lean:334` の sharp な `(p−1)i ≤ e_L` だけからは
`p^{i/e_L} ≤ axDecay p 1` しか出ず、`axDecay p d` は出ない。**

★足りないのは `NormalizedTraceDescent.lean:379 FirstJumpRoute` の
`hm : p^{k−1} ≤ m`（分母に `p^{k−1}` が入ること）であり、それは
★**最後の跳びではなく第 1 跳びで測る**ことから来る
（同 `:65-66` の表: `p = 3`・深さ 2 で「sharp(上の層) `i = 8` → `3^{4/9}` ★超える」vs
「sharp(**第 1 跳び**) `i₁ = 2` → `3^{1/9}` ★★収まる」）。
⇒ ★★★**本日の出口は「上の層の跳び」を使っており、そこが `axDecay p 1` で頭打ちになる。**

## ★①不分岐側との関係

★本波も独立である。上の否定は全分岐（`hvalK`）を仮定した出口についての主張であり、
`TotallyRamifiedLayer.lean:342 not_valK_of_norm_eq`（不分岐側）とは別である。

## ★穴の現状（本波の再測定）

| 穴 | 現状 |
|---|---|
| ①不分岐 | ★生きている |
| ②一様定数 | ★★**②′ に一般化した** —— 「定数 `> 1` で**下から**抑えられる」だけで足りる（§1） |
| ③幾何減衰 | ★生きている |
| ④測定点 | ★前波で消えた（記法の問題だった） |
| ★★⑤（本波） | ★**本日の出口の 1 段の損失は `axDecay p 1` が限界**。`cEx` は実現されない |

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ項目（pGC 物理 p.6 Corollary 3.1）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. §1 は `ℝ` と `Finset` だけ、§3 は `ℕ` と `rpow` だけで、分岐の語彙が 1 語も出ない。
4. ★`exitLoss` は「一般の `k` の出口が 1 段に払う額」を式にしたものであって、
   出口の宣言そのものではない（★宣言は `twr g p k` の元を返す形で、
   本ファイルはその**損失の値だけ**を取り出している）。
-/

namespace ABC3.Found.PGC

namespace ExitLossDepth

/-! ## §1 抽象核（②′）—— 損失が定数 `> 1` で**下から**抑えられれば積は非有界 -/

section Core

/-- ★★**②の一般化（②′）** —— 損失が `d ≥ 1` で定数 `C > 1` 以上なら
`Icc 1 n` 上の積は有界にできない。

★`PGroupDescentToAxWild.const_descent_no_go` は `c ≡ C`（一様定数）だったが、
★**下から `C` で抑えられるだけで十分**である。★これが本日の出口に効く形。
★分岐・付値・Galois・`p` 進の語彙が 1 語も出ない。 -/
theorem not_forall_prod_Icc_le_of_le {c : ℕ → ℝ} {C B : ℝ} (hC : 1 < C)
    (hc : ∀ d, 1 ≤ d → C ≤ c d) :
    ¬ (∀ n : ℕ, ∏ d ∈ Finset.Icc 1 n, c d ≤ B) := by
  intro h
  obtain ⟨n, hn⟩ := pow_unbounded_of_one_lt B hC
  have hle : ∏ d ∈ Finset.Icc 1 n, C ≤ ∏ d ∈ Finset.Icc 1 n, c d := by
    refine Finset.prod_le_prod (fun i _ => ?_) (fun i hi => hc i (Finset.mem_Icc.mp hi).1)
    linarith
  rw [Finset.prod_const, Nat.card_Icc, Nat.add_sub_cancel] at hle
  linarith [h n]

end Core

/-! ## §2 本日の出口が 1 段に払う額 -/

section Exit

variable {p : ℕ} [Fact p.Prime]

/-- ★一般の `k` の出口（`TotallyRamifiedLayer.lean:266`）が**1 段**に払う額。
★出口は `twr g p k`（`[M : twr g p k] = p`、1 段下）の元を返して
`∏_{j ∈ Icc 1 (k+1)} axDecay p j` を払う。 -/
noncomputable def exitLoss (p : ℕ) (d : ℕ) : ℝ := ∏ j ∈ Finset.Icc 1 d, axDecay p j

/-- ★★**深さに依らず `axDecay p 1` 以上** —— しかも `d` について増加する。 -/
theorem axDecay_one_le_exitLoss : ∀ d : ℕ, 1 ≤ d → axDecay p 1 ≤ exitLoss p d := by
  intro d
  induction d with
  | zero => intro h; omega
  | succ m ih =>
      intro _
      rcases Nat.eq_zero_or_pos m with hm | hm
      · subst hm
        simp [exitLoss]
      · have h1 := ih hm
        have h2 : (1:ℝ) ≤ axDecay p (m + 1) := one_le_axDecay p (m + 1)
        have h3 : (1:ℝ) ≤ axDecay p 1 := one_le_axDecay p 1
        have hstep : exitLoss p (m + 1) = exitLoss p m * axDecay p (m + 1) := by
          unfold exitLoss
          exact Finset.prod_Icc_succ_top (by omega) _
        rw [hstep]
        nlinarith

/-- ⇒ ★出口の損失を `AxWildDescent` の `c` にすると積は非有界。 -/
theorem not_forall_prod_exitLoss_le (B : ℝ) :
    ¬ (∀ n : ℕ, ∏ d ∈ Finset.Icc 1 n, exitLoss p d ≤ B) :=
  not_forall_prod_Icc_le_of_le (PGroupToAxWild.one_lt_axDecay_one (p := p))
    (fun d hd => axDecay_one_le_exitLoss (p := p) d hd)

/-- ★★`1` に落ちる段が **1 つも無い** ⇒ 前波の「悪い段が有限個」の逃げ道が使えない。 -/
theorem exitLoss_ne_one (d : ℕ) (hd : 1 ≤ d) : exitLoss p d ≠ 1 := by
  intro h
  have := axDecay_one_le_exitLoss (p := p) d hd
  rw [h] at this
  exact absurd this (not_le.mpr (PGroupToAxWild.one_lt_axDecay_one (p := p)))

end Exit

/-! ## §3 ★★★「最後の跳び」の sharp な評価では `axDecay p d` は出ない -/

section Sharp

variable {p : ℕ} [Fact p.Prime]

/-- ★★★**どこで止まるか** —— `d ≥ 2` では、sharp な `(p−1)·i ≤ e` を満たす `(i, e)` で
`axDecay p d < p^{i/e}` となるものが存在する（`i = 1`, `e = p−1`）。

⇒ ★`RamificationJumpBound.lean:334` の `(p−1)i ≤ e_L` だけからは
`p^{i/e_L} ≤ axDecay p 1` しか出ない。`axDecay p d` を出すには
`FirstJumpRoute` の `hm : p^{k−1} ≤ m`（＝**第 1 跳びで測る**こと）が要る。 -/
theorem last_jump_bound_not_enough {d : ℕ} (hd : 2 ≤ d) :
    ∃ i e : ℕ, 0 < e ∧ (p - 1) * i ≤ e ∧
      axDecay p d < (p : ℝ) ^ ((i : ℝ) / (e : ℝ)) := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  refine ⟨1, p - 1, by omega, by omega, ?_⟩
  obtain ⟨k, rfl⟩ : ∃ k, d = k + 1 := ⟨d - 1, by omega⟩
  have hk : 0 < k := by omega
  have hcast : ((1 : ℕ) : ℝ) / ((p - 1 : ℕ) : ℝ) = 1 / ((p : ℝ) - 1) := by
    rw [Nat.cast_one, Nat.cast_sub (by omega), Nat.cast_one]
  rw [hcast, ← axDecay_one p]
  exact PGroupToAxWild.axDecay_succ_lt_axDecay_one hk

end Sharp

/-! ## §4 ★★★まとめ —— `cEx` は実現されない -/

section Summary

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★**本ファイルの主結果** —— 本日の出口は前波の `cEx` を実現しない。

(a) 損失は `d ≥ 1` で常に `axDecay p 1` 以上、(b) `1` になる段が無い、
(c) 積は非有界。⇒ ★「悪い段が有限個」の逃げ道は使えない。 -/
theorem exit_does_not_realize_cEx (B : ℝ) :
    (∀ d : ℕ, 1 ≤ d → axDecay p 1 ≤ exitLoss p d)
      ∧ (∀ d : ℕ, 1 ≤ d → exitLoss p d ≠ 1)
      ∧ ¬ (∀ n : ℕ, ∏ d ∈ Finset.Icc 1 n, exitLoss p d ≤ B) :=
  ⟨fun d hd => axDecay_one_le_exitLoss (p := p) d hd,
   fun d hd => exitLoss_ne_one (p := p) d hd,
   not_forall_prod_exitLoss_le B⟩

/-- ★次数 `p` の層 1 枚の出口（`TotallyRamifiedLayer.lean:305`、損失は定数 `axDecay p 1`）
についても同じ。 -/
theorem degP_exit_does_not_realize_cEx (B : ℝ) :
    (axDecay p 1 ≠ 1) ∧ ¬ (∀ n : ℕ, ∏ _d ∈ Finset.Icc 1 n, axDecay p 1 ≤ B) :=
  ⟨ne_of_gt (PGroupToAxWild.one_lt_axDecay_one (p := p)),
   PGroupToAxWild.not_forall_prod_Icc_axDecay_one_le B⟩

end Summary

/-! ## §5 `.src` と 使っている公理の一覧 -/

section Src

def not_forall_prod_Icc_le_of_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def axDecay_one_le_exitLoss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def last_jump_bound_not_enough.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exit_does_not_realize_cEx.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms not_forall_prod_Icc_le_of_le
#print axioms axDecay_one_le_exitLoss
#print axioms not_forall_prod_exitLoss_le
#print axioms exitLoss_ne_one
#print axioms last_jump_bound_not_enough
#print axioms exit_does_not_realize_cEx
#print axioms degP_exit_does_not_realize_cEx

end ExitLossDepth

end ABC3.Found.PGC
