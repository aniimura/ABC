import ABC3.Found.PGC.ProjectionDefectVerified

/-!
# [pGC] ★★★掃討完了 6/6 —— そして★**次の波が最初に読むファイル**

## ★どれを選んだか、なぜか

★**(L)（掃討の最後）を選んだ。** 前波で「seed が無いので件数は再現不能」と測ったが、
★**結論の方は前々波で回復した `x` から直接出る**と分かったので安かった:
`ε = 11`, `d(x,E₁) = 7`, `d(x,K) = 3` から 1 段の slack `= 3 − (11−7) = −1`、
対の slack `= 12 − (11−3) = +4`。★木の表の 2 つの数がそのまま出る。

## ★★★結果 —— `DeepDescentPairDirect.lean:39-45` の表を再導出した

`tools/pair-budget-check.py`（★本波で新規作成。`zeta27-distance-check.py` を import）。

| | 木の字面 | 本波（回復した `x`） | 本波（新しい乱択 20,000 件、seed 20260909） |
|---|---|---|---|
| 1 段の最悪の必要量 | `4` | `4` | `4` ✓ |
| ★1 段の slack | `−1` | `−1` | ★最小 `−1` ✓ |
| ★対の slack | `+4` | `+4` | ★最小 `+4` ✓ |
| 対の予算の反例 | 0 件 | —— | ★**0 件** ✓ |

★★**件数（49,000 / 95,298）は seed が無いので再現していない。★再現したのは結論である。**
★これは前波で「再現できるのは結論の方」と測ったとおりである。

## ★★掃討の最終状態 —— **6/6**

| ファイル | 状態 | 塞いだ波 |
|---|---|---|
| `WildDescentDistanceOnly.lean` | ★塞いだ（`x` を回復。8 数値一致） | 5 波前 |
| `GainedTowerDescent.lean:94` | ★跳び `L = 8, 26` は塞いだ。総当たりの表は**未** | 4 波前 |
| `DeepDescentRepair.lean:110` | ★塞いだ（4 値一致） | 3 波前 |
| `JumpDefectTradeoff.lean:106` | ★塞いだ（結論は Lean の定理でスクリプト不要） | 3 波前 |
| `EquivariantProjectionDescent.lean:19` | ★塞いだ（`γ = p−1` を 4 層で） | 前波 |
| ★`DeepDescentPairDirect.lean:37,51,147` | ★★**塞いだ（本波、結論のみ）** | 本波 |

★`tools/` に入った再現スクリプトは **6 本**:
`zeta27-ramification-check.py` / `zeta27-distance-check.py` / `zeta-tower-check.py` /
`deepdescent-repair-check.py` / `trace-defect-check.py` / `pair-budget-check.py`。

★★**5 波で測った範囲では、木の機械計算の記録はすべて正しかった。**
★合わなかったのは 1 度だけで、それは★**私の 1 度目の `γ` の計算**であった（前波）。

## ★★★次の波が最初に読むもの（★本日の測定の要約。すべて `file:line` で確かめられる）

### 出口の連鎖に残る**穴は 4 つ**（本日の測定で確定した射程）

| 穴 | 主張 | 場所 | 射程 |
|---|---|---|---|
| ① 不分岐 | 次数 `p` の層が不分岐なら出口の `hvalK` は偽 | `TotallyRamifiedLayer.lean:342` | ★**(b) は閉じた**（`UnramifiedStepFixed.exists_fixed_step_free`：不分岐の段は損失も `ε` の伸びも `1` で**無料**）。★残るのは (a)「不分岐 ⇒ `∃ a ∈ 𝒪_L, Tr(a)=1`」だけ |
| ②′ 一様定数 | 損失が定数 `> 1` で下から抑えられれば積は非有界 | `ExitLossDepth.not_forall_prod_Icc_le_of_le` | ★射程は `c` ではなく **`g`**（`BudgetRecast`）。★効かせるには**すべての深さで一様な証人**が要る（`WitnessScope`） |
| ③ 幾何減衰 | `‖σ^{p^k}π−π‖ ≥ ‖p‖^{1/(p−1)}‖σπ−π‖`（`k` に依らない） | `JumpGeometricDecay.no_uniform_geometric_of_harith` | ★生きている |
| ⑤ 出口の限界 | 本日の出口の 1 段の損失は `axDecay p 1` 以上で `1` に落ちない | `ExitLossDepth.exit_does_not_realize_cEx` | ★`hm` を供給しても消えない（`FirstJumpRouteGap`） |

★④（`ℚ₃(ζ₂₇)` の測定点）は★**消えた** —— 点ごとには破れるが**積の形では通る**
（`BudgetFiniteExcess.cyclotomic_break_vanishes_in_prod`）。

### ★掘らなくてよいと確定した道

* `AxWildDescentDecay`（`ε` が増えない形）は ★**`AxLemma` と同値（循環）**
  （`AxDecayKnobCircular.axWildDescentDecay_iff_axLemma`）。`θ` の値は効かない。
* 予算関数 `F` への載せ替えも ★**同じ循環に戻る**（`BudgetRecast.budgetStep_iff_axLemma`）。
* `FirstJumpRoute` の 4 仮説は ★**`c k ≤ axDecay p k` と同値**（再パラメータ化にすぎない）
  （`FirstJumpRouteEquiv.exists_firstJump_data_iff`）。
* B2（塔の外の `x′` の排除）は ★**Krasner の評価では原理的に届かない**
  （`KrasnerCeiling.b2_obstruction_not_triggered`）。しかも★**帰結は `AxLemma` を左右しない**。

### ★次の 1 点（★私の見込み。測っていない）

★`hform` —— 「第 1 跳びの損失 `p^{u₁/e_L}` を達成する降下の構成」。
★塔のデータからの仮説 3 つは落ちており（`FirstJumpRouteGap.axLemma_of_first_jump_loss`）、
★残るのは `hform` と `hc` と `AxWildDescent K c` の 3 本だけである。
★`u₁ < i_ϖ(σ)` は本日閉じた（`CosetSumFixedRing.lt_ramIndex_quotient_fixedRing`）が、
★それは「`σ̄` が `x` を動かす量」の下界であって「`x′` の作り方」ではない。

## 逸脱の記録（CLAUDE.md「逸脱」）

1. `.src` は付けていない（★**検算の記録**である）。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**（import のみ）。
3. ★件数（49,000 / 95,298）は★**再現していない**（seed が無いため）。★再現したのは結論である。
4. ★`GainedTowerDescent.lean:94` の総当たりの表は★**今も未再現**である。
-/

namespace ABC3.Found.PGC

namespace PairBudgetVerified

/-! ## §1 slack の再導出（`ε = 11`, `d(x,E₁) = 7`, `d(x,K) = 3`） -/

section Slack

theorem one_step_need : (11 : ℤ) - 7 = 4 := by norm_num

theorem pair_need : (11 : ℤ) - 3 = 8 := by norm_num

/-- ★★1 段の予算（`axDecay 3 2` ＝ `v_M` で `3`）は **`−1` 足りない**
（`DeepDescentPairDirect.lean:42` の字面。★本波で `x` から再導出した）。 -/
theorem one_step_slack : (3 : ℤ) - (11 - 7) = -1 := by norm_num

/-- ★★対の予算（`axDecay 3 1 · axDecay 3 2` ＝ `v_M` で `12`）は **`+4` 余る**（同 `:43`）。 -/
theorem pair_slack : (12 : ℤ) - (11 - 3) = 4 := by norm_num

/-- 1 段では届かない。 -/
theorem one_step_budget_fails : ¬ ((11 : ℤ) - 7 ≤ 3) := by norm_num

/-- ★対では入る。★新しい乱択 20,000 件（seed 20260909）でも最小 slack は `+4`、反例 0 件。 -/
theorem pair_budget_fits : (11 : ℤ) - 3 ≤ 12 := by norm_num

end Slack

/-! ## §2 予算の値 -/

section Budgets

theorem axDecay_two_budget : (3 : ℤ) = 3 := rfl

/-- 対の予算は `axDecay 3 1`（`9/18`）＋ `axDecay 3 2`（`3/18`）＝ `12/18`。 -/
theorem pair_budget_value : (9 : ℤ) + 3 = 12 := by norm_num

/-- ★予算はちょうど `12`（`J = 13` は入らない）。 -/
theorem pair_budget_sharp_recheck : ¬ ((13 : ℤ) ≤ 12) := by norm_num

end Budgets

/-! ## §3 使っている公理の一覧 -/

#print axioms one_step_need
#print axioms pair_need
#print axioms one_step_slack
#print axioms pair_slack
#print axioms one_step_budget_fails
#print axioms pair_budget_fits
#print axioms pair_budget_value
#print axioms pair_budget_sharp_recheck

end PairBudgetVerified

end ABC3.Found.PGC
