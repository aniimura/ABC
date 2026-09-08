import ABC3.Found.PGC.DepthLowerBound
import ABC3.Found.PGC.BudgetFiniteExcess

/-!
# [pGC] ★証人は「点ごとの目標」を殺すが「積の形」は殺さない —— ②′には一様性が要る

## ★★★持ち場が「先に測る価値がある」と書いた点を、最初に測った

問い: ★**証人が各深さに在っても、④と同じ理由（積の形なら通る）で②′が効かない
可能性があるのではないか。** ⇒ ★**そのとおりだった。**

**(1) 証人は「点ごとの目標」を殺す**（§1 `not_axWildDescent_axDecay_of_witness_at`）。

深さ `d ≠ 0` に証人（`d(x, 深さ<d) ≥ C·ε` を満たす `x`）が 1 つ在り、
`axDecay p d < C` なら ★**`AxWildDescent K (axDecay p)` は偽**である。
★`C > 1` ならそういう `d` は必ず在る（§1 `exists_depth_refuting`、
中身は前波の `PGroupDescentToAxWild.exists_axDecay_lt`）。
⇒ §3 `not_axWildDescent_axDecay_of_deep_witness`:
★★**すべての深さに一様な `C > 1` の証人が在れば、目標 `AxWildDescent K (axDecay p)` は偽。**

**(2) ★しかし「積の形」は殺さない**（§2 `witness_at_one_depth_not_enough`）。

②′（`ExitLossDepth.not_forall_prod_Icc_le_of_le`）の仮説は
★**`∀ d ≥ 1, C ≤ c d`（`C` は `d` に依らない）**である。
★有限個の深さで下界が出ても、この仮説は出ない。証拠は前波の `cEx`
（`BudgetFiniteExcess`）そのもので:

* `3^{4/18} ≤ cEx 2`（★深さ 2 で下界を満たす）
* `∀ n, ∏_{j ∈ Icc 1 n} cEx j ≤ axConstant 3`（★積は有界）
* `¬ (∀ d ≥ 1, 3^{4/18} ≤ cEx d)`（★深さ 3 以降は `1` なので偽）

⇒ ★★★**②′が効くには「すべての深さで一様な `C > 1` の証人」が要る。**
★1 つの深さ（あるいは有限個）では足りない。★これは④が消えたのと**同じ構造**である。

## ★★どこで止まるか（`file:line`）—— 証人の構成の費用を測った

証人には `∀ x' : K.closure, wildDepth K x' < d → C·ε ≤ ‖x − x'‖`、
すなわち★**塔の外にある深さ `< d` の元も含めた**下界が要る。

★これは `WildDescentDistanceOnly.lean:88-99` が
★「`M` の外の深さ ≤ 1 の `x′` を排除できていない。Krasner の補題は
`‖x−x'‖ < min_σ‖σx−x‖` のときしか効かず、ここで要求される半径は `min` より**大きい**ので
届かない」「⇒ 確定させるには `K = ℚ₃(ζ₃)` の 40 個の巡回 3 次拡大で `d(x,N)` を
実際に測る必要がある。★測っていない」と書いた★**その未測定点と同じ形**である
（★本波で突き合わせた。同ファイルは今も未測定である）。

⇒ ★★**証人の構成は「塔の中の計算」ではなく「塔の外を排除する」問題であり、
本日の塔一式（`ConcreteNormedModelK1` / `JumpFromValueGroup.lean:332 harith_zeta81` など）
では届かない。** ★これが止まる場所である。

## ★穴の現状（本波の測定）

| 穴 | 本波の測定 |
|---|---|
| ①不分岐 | (b) 閉、(a) 残（変化なし） |
| ②′ | ★★**「すべての深さで一様な証人」が要ると確定**。有限個では足りない（④と同じ構造） |
| ③幾何減衰 | 変わらず |
| ⑤出口の限界 | ②′と同じ（射程は変わらないが、効かせるには一様性が要る） |
| ④ | 消えたまま |

★★**②′/⑤ は「効く」とも「効かない」とも確定していない。**
★確定させるのに要るのは「一様な証人の存在」で、それは塔の外を排除する問題である。

## ★①不分岐側との関係

★本波は①と**独立**である（`AxWildDescent` の定義と `axDecay` の収束だけを使う）。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は上流と同じ pGC 物理 p.6 Corollary 3.1。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★証人そのものは**構成していない**（上に理由を書いた）。本ファイルが定理にしたのは
   ★「証人が在れば何が言えて、何が言えないか」だけである。
4. §2 は `p = 3` の具体列 `cEx` に依存している（前波で建てたもの）。
   一般の `p` での対応する列は作っていない。
-/

namespace ABC3.Found.PGC

namespace WitnessScope

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 証人は「点ごとの目標」を殺す -/

section Pointwise

/-- ★★★**証人 1 つで点ごとの目標が死ぬ** —— 深さ `d ≠ 0` に
`d(x, 深さ<d) ≥ C·ε` を満たす `x` が在り `axDecay p d < C` なら
`AxWildDescent K (axDecay p)` は偽。

★中身は前波の `DepthLowerBound.le_c_of_lower_bound`（測定器）1 本である。 -/
theorem not_axWildDescent_axDecay_of_witness_at (K : PAdicLocalField p) {C : ℝ} {d : ℕ}
    (hd : d ≠ 0) (hlt : axDecay p d < C)
    (hwit : ∃ (ε : ℝ) (x : K.closure), 0 < ε ∧ wildDepth K x = d ∧
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) ∧
      (∀ x' : K.closure, wildDepth K x' < d → C * ε ≤ ‖x - x'‖)) :
    ¬ AxWildDescent K (axDecay p) := by
  intro h
  obtain ⟨ε, x, hε, hdx, hx, hlow⟩ := hwit
  have hdep : wildDepth K x ≠ 0 := by rw [hdx]; exact hd
  have hle := DepthLowerBound.le_c_of_lower_bound K h hε hdep hx
    (by rw [hdx]; exact hlow)
  rw [hdx] at hle
  linarith

/-- `1 < C` なら `axDecay p d < C` となる `d ≠ 0` が在る
（`PGroupDescentToAxWild.exists_axDecay_lt`）。 -/
theorem exists_depth_refuting {C : ℝ} (hC : 1 < C) :
    ∃ d : ℕ, d ≠ 0 ∧ axDecay p d < C := by
  obtain ⟨k, hk⟩ := PGroupToAxWild.exists_axDecay_lt (p := p) hC
  exact ⟨k + 1, by omega, hk⟩

end Pointwise

/-! ## §2 ★しかし「積の形」は殺さない —— ②′には一様性が要る -/

section Product

/-- ★★★**有限個の証人では②′の仮説は出ない** —— 前波の `cEx` が証拠。

`3^{4/18} ≤ cEx 2`（深さ 2 で下界）かつ `∏ ≤ axConstant 3`（積は有界）かつ
`¬ (∀ d ≥ 1, 3^{4/18} ≤ cEx d)`（深さ 3 以降は `1`）。
⇒ ★②′には「**すべての**深さで一様な `C > 1`」が要る。★④が消えたのと同じ構造。 -/
theorem witness_at_one_depth_not_enough :
    ∃ c : ℕ → ℝ, ((3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) ≤ c 2)
      ∧ (∀ n : ℕ, ∏ j ∈ Finset.Icc 1 n, c j ≤ axConstant 3)
      ∧ ¬ (∀ d : ℕ, 1 ≤ d → (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) ≤ c d) := by
  haveI : Fact (Nat.Prime 3) := ⟨Nat.prime_three⟩
  refine ⟨BudgetFinite.cEx, ?_, BudgetFinite.cEx_prod_le, ?_⟩
  · have h : BudgetFinite.cEx 2 = (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) := by
      norm_num [BudgetFinite.cEx]
    rw [h]
  · intro hall
    have h3 := hall 3 (by omega)
    rw [BudgetFinite.cEx_eq_one_of_gt 3 (by omega)] at h3
    have hgt : (1 : ℝ) < (3 : ℝ) ^ (((4 : ℕ) : ℝ) / ((18 : ℕ) : ℝ)) := by
      refine Real.one_lt_rpow_iff_of_pos (by norm_num) |>.mpr (Or.inl ⟨by norm_num, ?_⟩)
      norm_num
    linarith

end Product

/-! ## §3 まとめ -/

section Summary

/-- ★★すべての深さに一様な `C > 1` の証人が在れば、目標は偽。 -/
theorem not_axWildDescent_axDecay_of_deep_witness (K : PAdicLocalField p) {C : ℝ} (hC : 1 < C)
    (hwit : ∀ d : ℕ, d ≠ 0 → ∃ (ε : ℝ) (x : K.closure), 0 < ε ∧ wildDepth K x = d ∧
      (∀ σ : K.absGal, ‖σ • x - x‖ ≤ ε) ∧
      (∀ x' : K.closure, wildDepth K x' < d → C * ε ≤ ‖x - x'‖)) :
    ¬ AxWildDescent K (axDecay p) := by
  obtain ⟨d, hd, hlt⟩ := exists_depth_refuting (p := p) hC
  exact not_axWildDescent_axDecay_of_witness_at K hd hlt (hwit d hd)

end Summary

/-! ## §4 `.src` と 使っている公理の一覧 -/

section Src

def not_axWildDescent_axDecay_of_witness_at.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def witness_at_one_depth_not_enough.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def not_axWildDescent_axDecay_of_deep_witness.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

end Src

#print axioms not_axWildDescent_axDecay_of_witness_at
#print axioms exists_depth_refuting
#print axioms witness_at_one_depth_not_enough
#print axioms not_axWildDescent_axDecay_of_deep_witness

end WitnessScope

end ABC3.Found.PGC
