import Mathlib.NumberTheory.Padics.PadicNumbers
import Mathlib.Topology.Algebra.InfiniteSum.Nonarchimedean
import Mathlib.Analysis.SpecificLimits.Normed
import Mathlib.NumberTheory.Padics.ProperSpace
import Mathlib.Analysis.Normed.Group.Ultra
import Mathlib.Analysis.Normed.Group.FunctionSeries
import Mathlib.Topology.MetricSpace.Ultra.ContinuousMaps

/-!
# Track B — p進対数(`1 + 𝔪` 上)

`Found/IUTchIII/LogShell.lean` の測定で、log-shell を塞いでいる唯一の欠落が
**p進対数**だと分かった([AbsTopIII] 物理 p.3 の `log_{k̄}(𝒪^×_k)`)。
ここではそれを実際に作りにいく。

## `sorry` 無しで作れたもの

1. `one_div_le_norm_natCast` — `‖(m : ℚ_[p])‖ ≥ 1/m`(m ≠ 0)。
2. `tendsto_logTerm` / `summable_logTerm` — `‖x‖ < 1` で**総和可能**。
3. `logOneAdd` — **定義**。
4. `continuousOn_logOneAdd` — 半径 `r < 1` の閉球上で**連続**。
5. `norm_logOneAdd_eq` / **`logOneAdd_ne_zero`** — ★**非退化性**。
   `‖x‖ ≤ 1/4` なら `‖logOneAdd x‖ = ‖x‖`。ゆえに `logOneAdd` は
   `log := 0` のような退化定義ではない——**Lean で示した**。
6. `maximalIdeal_eq_closedBall` — `{x : ‖x‖ < 1} = closedBall 0 (1/p)`。
7. **`logShell`** / **`isCompact_logShell'`** — ★★**ℚ_[p] の log-shell とその
   コンパクト性**([AbsTopIII] p.5 の (L1))。`logShell_ne_singleton_zero` も。

## ★まだ作れていないもの

- ★**乗法性** `logOneAdd` が `(1+𝔪, ×) → (ℚ_p, +)` の準同型であること。
  mathlib に **`PowerSeries.log` が無い**のが直接の障害
  (`PowerSeries.exp` / `sin` / `cos` はある)。したがって形式的恒等式
  `log((1+X)(1+Y)) = log(1+X) + log(1+Y)` が使えない。
- `𝒪^×` 全体への延長(`log(ζu) := log(u)`)。これは乗法性を前提とする規約なので、
  上が済むまで意味を持たない。**したがって現状の `logShell` は
  `log(1+𝔪)` であって `log(𝒪^×)` ではない。**
- 一般の有限次拡大 K への一般化。
-/

namespace ABC3.Found.IUTchIII

open Filter

variable {p : ℕ} [Fact p.Prime]

/-- **`‖(m : ℚ_[p])‖ ≥ 1/m`**。

`‖m‖ = p^(-v_p(m))` かつ `p^(v_p(m)) ≤ m`(`Nat.ordProj_le`)による。
対数級数の項が 0 に収束することの、唯一の非自明な入力。 -/
theorem one_div_le_norm_natCast {m : ℕ} (hm : m ≠ 0) :
    (1 : ℝ) / m ≤ ‖((m : ℕ) : ℚ_[p])‖ := by
  have hne : ((m : ℕ) : ℚ_[p]) ≠ 0 := Nat.cast_ne_zero.mpr hm
  rw [Padic.norm_eq_zpow_neg_valuation hne, Padic.valuation_natCast]
  have hp1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hple : (p : ℕ) ^ (padicValNat p m) ≤ m := by
    have := Nat.ordProj_le (n := m) p hm
    rwa [Nat.factorization_def _ (Fact.out : p.Prime)] at this
  have hpleR : ((p : ℝ)) ^ (padicValNat p m) ≤ (m : ℝ) := by exact_mod_cast hple
  rw [zpow_neg, zpow_natCast]
  rw [one_div, inv_le_inv₀ (by positivity) (by positivity)]
  exact hpleR

/-- 対数級数 `log(1+x) = Σ_{n≥0} (-1)^n x^{n+1}/(n+1)` の第 n 項。 -/
noncomputable def logTerm (x : ℚ_[p]) (n : ℕ) : ℚ_[p] :=
  (-1) ^ n * x ^ (n + 1) / ((n + 1 : ℕ) : ℚ_[p])

/-- 項のノルムの上界: `‖logTerm x n‖ ≤ (n+1) * ‖x‖^(n+1)`。 -/
theorem norm_logTerm_le (x : ℚ_[p]) (n : ℕ) :
    ‖logTerm x n‖ ≤ (n + 1 : ℕ) * ‖x‖ ^ (n + 1) := by
  have hne : ((n + 1 : ℕ) : ℚ_[p]) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.succ_ne_zero n)
  have hlow : (1 : ℝ) / (n + 1 : ℕ) ≤ ‖((n + 1 : ℕ) : ℚ_[p])‖ :=
    one_div_le_norm_natCast (Nat.succ_ne_zero n)
  have hpos : (0 : ℝ) < ‖((n + 1 : ℕ) : ℚ_[p])‖ := norm_pos_iff.mpr hne
  rw [logTerm, norm_div, norm_mul, norm_pow, norm_pow, norm_neg, norm_one, one_pow, one_mul]
  rw [div_le_iff₀ hpos]
  calc ‖x‖ ^ (n + 1)
      = ‖x‖ ^ (n + 1) * ((n + 1 : ℕ) * (1 / (n + 1 : ℕ))) := by
        field_simp
    _ ≤ ‖x‖ ^ (n + 1) * ((n + 1 : ℕ) * ‖((n + 1 : ℕ) : ℚ_[p])‖) := by
        gcongr
    _ = (n + 1 : ℕ) * ‖x‖ ^ (n + 1) * ‖((n + 1 : ℕ) : ℚ_[p])‖ := by ring

/-- **`‖x‖ < 1` なら項は 0 に収束する**。 -/
theorem tendsto_logTerm {x : ℚ_[p]} (hx : ‖x‖ < 1) :
    Tendsto (logTerm x) atTop (nhds 0) := by
  have hbound : Tendsto (fun n : ℕ => (n + 1 : ℕ) * ‖x‖ ^ (n + 1)) atTop (nhds 0) := by
    have h0 : Tendsto (fun n : ℕ => (n : ℝ) * ‖x‖ ^ n) atTop (nhds 0) :=
      tendsto_self_mul_const_pow_of_lt_one (norm_nonneg x) hx
    have := (tendsto_add_atTop_iff_nat (f := fun n : ℕ => (n : ℝ) * ‖x‖ ^ n) 1).mpr h0
    simpa using this
  exact squeeze_zero_norm (fun n => norm_logTerm_le x n) hbound

/-- **総和可能**。非アルキメデス完備群では「項が 0 に行く」で足りる
(`NonarchimedeanAddGroup.summable_of_tendsto_cofinite_zero`)。 -/
theorem summable_logTerm {x : ℚ_[p]} (hx : ‖x‖ < 1) : Summable (logTerm x) :=
  NonarchimedeanAddGroup.summable_of_tendsto_cofinite_zero
    (by rw [Nat.cofinite_eq_atTop]; exact tendsto_logTerm hx)

/-- **p進対数**(`‖x‖ < 1` の範囲で `log(1+x)`)。

原文 (AbsTopIII p.3):
> — which is compact — may be thought of as defining a sort of canonical rigid
> integral structure on k. In the present paper, we shall refer to the “canonical
> rigid integral structures” obtained in this way as log-shells.

★`‖x‖ ≥ 1` では級数が総和可能とは限らないので、この `tsum` は
**`‖x‖ < 1` の外では意味を持たない**(mathlib の `tsum` は非総和可能なとき `0` を返す)。 -/
noncomputable def logOneAdd (x : ℚ_[p]) : ℚ_[p] := ∑' n : ℕ, logTerm x n

/-- `log(1 + 0) = 0`。

★非退化性(`logOneAdd` が恒等的に 0 でないこと)は本ファイル下部の
`logOneAdd_ne_zero` で証明した。 -/
@[simp] theorem logOneAdd_zero : logOneAdd (0 : ℚ_[p]) = 0 := by
  rw [logOneAdd]
  convert tsum_zero with n
  simp [logTerm]

/-- 補助: `Σ (n+1) r^(n+1)` は `0 ≤ r < 1` で総和可能。 -/
theorem summable_bound {r : ℝ} (hr0 : 0 ≤ r) (hr : r < 1) :
    Summable (fun n : ℕ => ((n : ℝ) + 1) * r ^ (n + 1)) := by
  have hnorm : ‖r‖ < 1 := by rwa [Real.norm_eq_abs, abs_of_nonneg hr0]
  have h1 : Summable (fun n : ℕ => (n : ℝ) * r ^ n) := by
    simpa using summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) 1 hnorm
  have h2 : Summable (fun n : ℕ => r ^ n) := summable_geometric_of_lt_one hr0 hr
  refine ((h1.add h2).mul_left r).congr (fun n => ?_)
  ring

/-- ★**半径 `r < 1` の閉球の上で連続**。

級数の項が `(n+1) r^{n+1}` で一様に押さえられるので `continuousOn_tsum` が使える。 -/
theorem continuousOn_logOneAdd {r : ℝ} (hr0 : 0 ≤ r) (hr : r < 1) :
    ContinuousOn (logOneAdd (p := p)) (Metric.closedBall 0 r) := by
  refine continuousOn_tsum (u := fun n : ℕ => ((n : ℝ) + 1) * r ^ (n + 1))
    (fun i => Continuous.continuousOn ?_) (summable_bound hr0 hr) (fun n x hx => ?_)
  · unfold logTerm; fun_prop
  · refine (norm_logTerm_le x n).trans ?_
    have hxr : ‖x‖ ≤ r := by simpa using hx
    have hpow : ‖x‖ ^ (n + 1) ≤ r ^ (n + 1) := by gcongr
    have : ((n : ℝ) + 1) = ((n + 1 : ℕ) : ℝ) := by push_cast; ring
    rw [this]
    gcongr

/-- ★★**本物の log-shell が1つ完成**: `log(1 + B_r)` は**コンパクト**。

`B_r`(半径 `r < 1` の閉球)はコンパクト(`ℚ_[p]` は proper)、
`logOneAdd` はその上で連続なので、像はコンパクト。

`𝔪 = {x : ‖x‖ < 1}` との同定は `maximalIdeal_eq_closedBall` で証明済み。
それを使って名前を付けたものが `logShell` / `isCompact_logShell'`。 -/
theorem isCompact_logShell {r : ℝ} (hr0 : 0 ≤ r) (hr : r < 1) :
    IsCompact (logOneAdd (p := p) '' Metric.closedBall 0 r) :=
  (isCompact_closedBall (0 : ℚ_[p]) r).image_of_continuousOn (continuousOn_logOneAdd hr0 hr)

/-! ## ★非退化性 — `logOneAdd` は恒等的に 0 ではない -/

/-- 補助: `n + 2 ≤ 2 · 4^n`。 -/
private theorem aux_le_pow (n : ℕ) : (n : ℝ) + 2 ≤ 2 * 4 ^ n := by
  induction n with
  | zero => norm_num
  | succ k ih =>
    have h4 : (0 : ℝ) ≤ 4 ^ k := by positivity
    have hpow : (4 : ℝ) ^ (k + 1) = 4 * 4 ^ k := by ring
    push_cast [hpow]
    linarith [ih, h4, Nat.cast_nonneg (α := ℝ) k]

/-- 尾部の一様評価: `‖x‖ ≤ 1/4` なら `n ≥ 1` の項はすべて `2‖x‖²` 以下。 -/
theorem norm_logTerm_succ_le {x : ℚ_[p]} (hx : ‖x‖ ≤ 1 / 4) (n : ℕ) :
    ‖logTerm x (n + 1)‖ ≤ 2 * ‖x‖ ^ 2 := by
  have h0 : (0 : ℝ) ≤ ‖x‖ := norm_nonneg x
  refine (norm_logTerm_le x (n + 1)).trans ?_
  have hstep : ‖x‖ ^ n ≤ (1 / 4 : ℝ) ^ n := by gcongr
  have hkey : ((n : ℝ) + 2) * ‖x‖ ^ n ≤ 2 := by
    calc ((n : ℝ) + 2) * ‖x‖ ^ n
        ≤ ((n : ℝ) + 2) * (1 / 4 : ℝ) ^ n := by gcongr
      _ = ((n : ℝ) + 2) / 4 ^ n := by rw [div_pow, one_pow]; ring
      _ ≤ 2 := by
          rw [div_le_iff₀ (by positivity)]
          linarith [aux_le_pow n]
  calc ((n + 1 + 1 : ℕ) : ℝ) * ‖x‖ ^ (n + 1 + 1)
      = (((n : ℝ) + 2) * ‖x‖ ^ n) * ‖x‖ ^ 2 := by push_cast; ring
    _ ≤ 2 * ‖x‖ ^ 2 := by gcongr

/-- ★**`‖x‖ ≤ 1/4` なら `‖logOneAdd x‖ = ‖x‖`**。

超距離性による: 第 0 項は `x` そのもので、残りの尾部は `2‖x‖² < ‖x‖` で押さえられる。 -/
theorem norm_logOneAdd_eq {x : ℚ_[p]} (hx : ‖x‖ ≤ 1 / 4) (hx0 : x ≠ 0) :
    ‖logOneAdd x‖ = ‖x‖ := by
  have hlt1 : ‖x‖ < 1 := lt_of_le_of_lt hx (by norm_num)
  have hs : Summable (logTerm x) := summable_logTerm hlt1
  have h0 : logTerm x 0 = x := by simp [logTerm]
  have hpos : 0 < ‖x‖ := norm_pos_iff.mpr hx0
  have htail : ‖∑' n : ℕ, logTerm x (n + 1)‖ ≤ 2 * ‖x‖ ^ 2 :=
    IsUltrametricDist.norm_tsum_le_of_forall_le_of_nonneg (by positivity)
      (fun n => norm_logTerm_succ_le hx n)
  have hsmall : 2 * ‖x‖ ^ 2 < ‖x‖ := by nlinarith
  have hne : ‖∑' n : ℕ, logTerm x (n + 1)‖ ≠ ‖x‖ := ne_of_lt (lt_of_le_of_lt htail hsmall)
  rw [logOneAdd, hs.tsum_eq_zero_add, h0]
  rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (Ne.symm hne)]
  exact max_eq_left (le_of_lt (lt_of_le_of_lt htail hsmall))

/-- ★★**非退化性**: `logOneAdd` は恒等的に 0 ではない。

`x := p²` を取ると `‖x‖ = p^{-2} ≤ 1/4` なので上の等式が使え、
`‖logOneAdd x‖ = ‖x‖ ≠ 0`。

**これが無いと `logOneAdd` は `log := 0` と区別できない**——
`Meta/Calibration.lean` の事実1がまさにこの形の失敗を扱っている。 -/
theorem logOneAdd_ne_zero : ∃ x : ℚ_[p], ‖x‖ < 1 ∧ logOneAdd x ≠ 0 := by
  have hp2 : 2 ≤ (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).two_le
  have hppos : (0 : ℝ) < (p : ℝ) := by linarith
  refine ⟨(p : ℚ_[p]) ^ 2, ?_, ?_⟩
  · rw [norm_pow, Padic.norm_p]
    rw [show ((p : ℝ))⁻¹ ^ 2 = ((p : ℝ) ^ 2)⁻¹ by rw [inv_pow]]
    rw [inv_lt_one_iff₀]
    right; nlinarith
  · have hx0 : ((p : ℚ_[p]) ^ 2) ≠ 0 :=
      pow_ne_zero 2 (Nat.cast_ne_zero.mpr (Fact.out : p.Prime).ne_zero)
    have hle : ‖(p : ℚ_[p]) ^ 2‖ ≤ 1 / 4 := by
      rw [norm_pow, Padic.norm_p, show ((p : ℝ))⁻¹ ^ 2 = ((p : ℝ) ^ 2)⁻¹ by rw [inv_pow]]
      rw [inv_le_iff_one_le_mul₀ (by positivity)]
      nlinarith
    intro hzero
    have := norm_logOneAdd_eq hle hx0
    rw [hzero, norm_zero] at this
    exact hx0 (norm_eq_zero.mp this.symm)

/-! ## ★ℚ_p の log-shell -/

/-- `ℚ_[p]` の極大イデアル `𝔪 = {x : ‖x‖ < 1}` は、付値が離散なので
半径 `1/p` の**閉**球に等しい。 -/
theorem maximalIdeal_eq_closedBall :
    {x : ℚ_[p] | ‖x‖ < 1} = Metric.closedBall 0 ((p : ℝ)⁻¹) := by
  have hp1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  ext x
  simp only [Set.mem_setOf_eq, Metric.mem_closedBall, dist_zero_right]
  constructor
  · intro hlt
    rcases eq_or_ne x 0 with rfl | hx0
    · simp
    · rw [Padic.norm_eq_zpow_neg_valuation hx0] at hlt ⊢
      have hneg : x.valuation < 0 ∨ 0 ≤ x.valuation := lt_or_ge _ _
      rcases hneg with h | h
      · -- valuation < 0 なら ‖x‖ = p^{-v} > 1、仮定に反する
        exfalso
        have : (1 : ℝ) < (p : ℝ) ^ (-x.valuation) := by
          apply one_lt_zpow₀ hp1; omega
        linarith
      · -- valuation ≥ 0。‖x‖ < 1 から valuation ≠ 0、よって valuation ≥ 1
        have hne0 : x.valuation ≠ 0 := by
          intro h0
          rw [h0] at hlt; simp at hlt
        have h1 : (1 : ℤ) ≤ x.valuation := by omega
        rw [show ((p : ℝ))⁻¹ = (p : ℝ) ^ (-1 : ℤ) by simp]
        exact zpow_le_zpow_right₀ hp1.le (by omega)
  · intro hle
    calc ‖x‖ ≤ (p : ℝ)⁻¹ := hle
      _ < 1 := by rw [inv_lt_one_iff₀]; right; exact hp1

/-- ★★**ℚ_[p] の log-shell** —— [AbsTopIII] 物理 p.3 の `log_{k̄}(𝒪^×_k)` の
`1 + 𝔪` 部分。`𝒪^×` 全体への延長は乗法性(未証明)を要するので、ここでは
`log(1 + 𝔪)` を取る(module docstring の注意を参照)。 -/
noncomputable def logShell : Set ℚ_[p] := logOneAdd '' {x : ℚ_[p] | ‖x‖ < 1}

/-- ★★**log-shell はコンパクト** —— [AbsTopIII] p.5 の (L1)。

原文 (AbsTopIII p.5):
> (L1) a log-shell is compact and hence of finite “log-volume” [cf. Corollary

**これが本課題の目的だった対象である。** `sorry` 無し。 -/
theorem isCompact_logShell' : IsCompact (logShell (p := p)) := by
  have hp1 : (1 : ℝ) < (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have h0 : (0 : ℝ) ≤ (p : ℝ)⁻¹ := by positivity
  have h1 : ((p : ℝ))⁻¹ < 1 := by rw [inv_lt_one_iff₀]; right; exact hp1
  rw [logShell, maximalIdeal_eq_closedBall]
  exact isCompact_logShell h0 h1

/-- log-shell は `{0}` ではない —— 非退化性の帰結。 -/
theorem logShell_ne_singleton_zero : (logShell (p := p)) ≠ {0} := by
  obtain ⟨x, hx, hne⟩ := logOneAdd_ne_zero (p := p)
  intro h
  have : logOneAdd x ∈ logShell (p := p) := ⟨x, hx, rfl⟩
  rw [h, Set.mem_singleton_iff] at this
  exact hne this

end ABC3.Found.IUTchIII
