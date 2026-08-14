import ABC3.Found.IUTchIII.PadicLog
import ABC3.Found.IUTchIII.PowerSeriesLog
import Mathlib.RingTheory.PowerSeries.Substitution

namespace ABC3.Found.IUTchIII

open PowerSeries

variable {p : ℕ} [Fact p.Prime]

/-! ### 部品1: 形式的対数の係数のノルム -/

/-- `‖coeff d (log ℚ_p)‖ ≤ d`。分母は `d` ただ1つなので、`‖1/d‖_p ≤ d` で足りる。 -/
theorem norm_coeff_log_le (d : ℕ) : ‖coeff d (log ℚ_[p])‖ ≤ d := by
  rcases Nat.eq_zero_or_pos d with rfl | hd
  · simp
  · rw [coeff_log, if_neg hd.ne']
    have hcast : (algebraMap ℚ ℚ_[p]) ((-1 : ℚ) ^ (d + 1) / d)
        = (-1 : ℚ_[p]) ^ (d + 1) / (d : ℚ_[p]) := by
      rw [map_div₀, map_pow, map_neg, map_one, map_natCast]
    rw [hcast, norm_div, norm_pow, norm_neg, norm_one, one_pow]
    have h1 : (1 : ℝ) / d ≤ ‖((d : ℕ) : ℚ_[p])‖ := one_div_le_norm_natCast hd.ne'
    have hdpos : (0 : ℝ) < d := by exact_mod_cast hd
    have hnpos : 0 < ‖((d : ℕ) : ℚ_[p])‖ := lt_of_lt_of_le (by positivity) h1
    rw [div_le_iff₀ hnpos]
    calc (1 : ℝ) = d * (1 / d) := by field_simp
      _ ≤ d * ‖((d : ℕ) : ℚ_[p])‖ := by
          exact mul_le_mul_of_nonneg_left h1 hdpos.le

/-! ### 部品2: 多項式の冪の係数のノルム -/

/-- 係数がすべてノルム `≤ r` なら、`P^d` の係数はノルム `≤ r^d`。超距離だから成り立つ。 -/
theorem norm_coeff_pow_le {P : Polynomial ℚ_[p]} {r : ℝ} (hr0 : 0 ≤ r)
    (hr : ∀ n, ‖P.coeff n‖ ≤ r) (d : ℕ) : ∀ n, ‖(P ^ d).coeff n‖ ≤ r ^ d := by
  induction d with
  | zero =>
    intro n
    simp only [pow_zero, Polynomial.coeff_one]
    split_ifs <;> simp
  | succ k ih =>
    intro n
    rw [pow_succ, Polynomial.coeff_mul]
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    rintro ⟨i, j⟩ -
    rw [norm_mul, pow_succ]
    exact mul_le_mul (ih i) (hr j) (norm_nonneg _) (by positivity)

/-! ### 部品3: 係数の総和 = `eval 1` -/

/-- 多項式の係数の(無限和としての)総和は `eval 1` に等しい。有限和なので解析は要らない。 -/
theorem tsum_coeff_eq_eval_one (Q : Polynomial ℚ_[p]) : ∑' n, Q.coeff n = Q.eval 1 := by
  rw [tsum_eq_sum (s := Finset.range (Q.natDegree + 1)) ?_, Polynomial.eval_eq_sum_range]
  · simp only [one_pow, mul_one]
  · intro n hn
    have h : Q.natDegree + 1 ≤ n := by simpa using hn
    exact Q.coeff_eq_zero_of_natDegree_lt h

/-! ### 部品4: `logOneAdd` を形式的対数の係数で書き直す -/

/-- `logOneAdd w = ∑' d, coeff d (log ℚ_p) * w ^ d`。単なる添字のずらし。 -/
theorem logOneAdd_eq_tsum_coeff_log {w : ℚ_[p]} (hw : ‖w‖ < 1) :
    logOneAdd w = ∑' d : ℕ, coeff d (log ℚ_[p]) * w ^ d := by
  have key : ∀ d : ℕ, coeff (d + 1) (log ℚ_[p]) * w ^ (d + 1) = logTerm w d := by
    intro d
    rw [coeff_log, if_neg (Nat.succ_ne_zero d), logTerm, map_div₀, map_pow, map_neg,
      map_one, map_natCast]
    ring
  have hsummable : Summable fun d : ℕ => coeff d (log ℚ_[p]) * w ^ d := by
    rw [← summable_nat_add_iff 1]
    simpa only [key] using summable_logTerm hw
  rw [hsummable.tsum_eq_zero_add]
  have h0 : coeff 0 (log ℚ_[p]) * w ^ 0 = 0 := by simp
  rw [h0, zero_add, logOneAdd]
  exact (tsum_congr key).symm

/-! ### 部品5: `P^d` の係数の台 -/

/-- 定数項が 0 なら `P^d` は次数 `d` 未満の係数を持たない。 -/
theorem coeff_pow_eq_zero_of_lt {P : Polynomial ℚ_[p]} (hP0 : P.coeff 0 = 0) :
    ∀ (d n : ℕ), n < d → (P ^ d).coeff n = 0 := by
  intro d
  induction d with
  | zero => intro n hn; exact absurd hn (Nat.not_lt_zero n)
  | succ k ih =>
    intro n hn
    rw [pow_succ, Polynomial.coeff_mul]
    refine Finset.sum_eq_zero ?_
    rintro ⟨i, j⟩ hij
    rw [Finset.mem_antidiagonal] at hij
    rcases Nat.eq_zero_or_pos j with rfl | hj
    · rw [hP0, mul_zero]
    · rw [ih i (by omega), zero_mul]

/-! ### 部品6: ★橋 —— 形式的対数の係数の総和 = p進対数 -/

/-- ★**評価の橋**。定数項 0・係数のノルムが `r < 1` 以下の多項式 `P` に対し、
形式的対数 `logOf (1 + P)` の**係数の総和**は `logOneAdd (P.eval 1)` に等しく、
その総和は収束する。

これは mathlib の `eval₂` / `aeval` を**使っていない**。あれらは
`IsLinearTopology S S` を要求し、体である ℚ_p は満たさない(2026-08-15 実測)。
ここでは形式的な側(`ℚ_p⟦X⟧`、X進位相なので線形位相)で係数の等式を取り、
ℚ_p 側では2重級数を手で並べ替えている。 -/
theorem summable_and_tsum_coeff_logOf (P : Polynomial ℚ_[p]) (hP0 : P.coeff 0 = 0)
    {r : ℝ} (hr0 : 0 ≤ r) (hr1 : r < 1) (hr : ∀ n, ‖P.coeff n‖ ≤ r) :
    (Summable fun n => PowerSeries.coeff n (logOf (1 + (P : ℚ_[p]⟦X⟧)))) ∧
      ∑' n, PowerSeries.coeff n (logOf (1 + (P : ℚ_[p]⟦X⟧))) = logOneAdd (P.eval 1) := by
  classical
  set G : ℕ → ℕ → ℚ_[p] := fun d n => coeff d (log ℚ_[p]) * (P ^ d).coeff n with hG
  -- 台: `n < d` でも `d * deg P < n` でも消える
  have hGzero_d : ∀ d n, n < d → G d n = 0 := fun d n h => by
    simp [hG, coeff_pow_eq_zero_of_lt hP0 d n h]
  have hGzero_n : ∀ d n, d * P.natDegree < n → G d n = 0 := fun d n h => by
    have h2 : (P ^ d).coeff n = 0 :=
      Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_le_of_lt Polynomial.natDegree_pow_le h)
    simp [hG, h2]
  -- ノルム評価: `‖G d n‖ ≤ d * r ^ d`
  have hGbound : ∀ d n, ‖G d n‖ ≤ d * r ^ d := by
    intro d n
    rw [hG, norm_mul]
    exact mul_le_mul (norm_coeff_log_le d) (norm_coeff_pow_le hr0 hr d n)
      (norm_nonneg _) (Nat.cast_nonneg d)
  -- 各方向の総和可能性(有限台)
  have h₁ : ∀ d, Summable (G d) := fun d =>
    (hasSum_sum_of_ne_finset_zero (s := Finset.range (d * P.natDegree + 1))
      (fun n hn => hGzero_n d n (by simpa using hn))).summable
  have h₂ : ∀ n, Summable fun d => G d n := fun n =>
    (hasSum_sum_of_ne_finset_zero (s := Finset.range (n + 1))
      (fun d hd => hGzero_d d n (by simpa using hd))).summable
  -- ★2重族の総和可能性: 非アルキメデスなので「余有限で 0 に行く」だけでよい
  have hF : Summable (Function.uncurry G) := by
    refine NonarchimedeanAddGroup.summable_of_tendsto_cofinite_zero ?_
    rw [NormedAddGroup.tendsto_nhds_zero]
    intro ε hε
    obtain ⟨D, hD⟩ := (Metric.tendsto_atTop.mp
      (tendsto_self_mul_const_pow_of_lt_one hr0 hr1) ε hε)
    rw [Filter.eventually_cofinite]
    refine Set.Finite.subset
      ((Finset.range D ×ˢ Finset.range (D * P.natDegree + 1)).finite_toSet) ?_
    rintro ⟨d, n⟩ hq
    simp only [Set.mem_setOf_eq, not_lt] at hq
    have hne : G d n ≠ 0 := fun h => by simp [Function.uncurry, h] at hq; linarith
    have hdlt : d < D := by
      by_contra hcon
      have := hD d (not_lt.mp hcon)
      rw [Real.dist_eq, sub_zero] at this
      have hle : ‖G d n‖ ≤ d * r ^ d := hGbound d n
      have : (d : ℝ) * r ^ d < ε := lt_of_le_of_lt (le_abs_self _) this
      exact absurd hq (by simp only [Function.uncurry]; linarith)
    have hnle : n ≤ d * P.natDegree := by
      by_contra hcon
      exact hne (hGzero_n d n (not_le.mp hcon))
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_range]
    exact ⟨hdlt, by nlinarith [Nat.mul_le_mul_right P.natDegree (Nat.le_of_lt_succ
      (Nat.lt_succ_of_lt hdlt))]⟩
  -- 係数の等式
  have hsub : HasSubst ((P : ℚ_[p]⟦X⟧)) :=
    HasSubst.of_constantCoeff_zero'
      (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, Polynomial.coeff_coe]; exact hP0)
  have hcoeff : ∀ n, PowerSeries.coeff n (logOf (1 + (P : ℚ_[p]⟦X⟧))) = ∑' d, G d n := by
    intro n
    have hlog : logOf (1 + (P : ℚ_[p]⟦X⟧)) = (log ℚ_[p]).subst (P : ℚ_[p]⟦X⟧) := by
      rw [logOf_eq, add_sub_cancel_left]
    rw [hlog, coeff_subst' hsub]
    have hterm : ∀ d, coeff d (log ℚ_[p]) • PowerSeries.coeff n ((P : ℚ_[p]⟦X⟧) ^ d) = G d n := by
      intro d
      rw [smul_eq_mul, hG, ← Polynomial.coe_pow, Polynomial.coeff_coe]
    simp only [hterm]
    rw [finsum_eq_sum_of_support_subset (s := Finset.range (n + 1)),
      tsum_eq_sum (s := Finset.range (n + 1)) (fun d hd => hGzero_d d n (by simpa using hd))]
    intro d hd
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_range, Set.mem_Iio]
    by_contra hcon
    exact hd (hGzero_d d n (by omega))
  -- `‖P.eval 1‖ < 1`
  have hEval : ‖P.eval 1‖ < 1 := by
    rw [Polynomial.eval_eq_sum_range]
    refine lt_of_le_of_lt (IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hr0 ?_) hr1
    intro i _
    simpa using hr i
  refine ⟨?_, ?_⟩
  · have hswap : Summable (Function.uncurry fun n d => G d n) := hF.prod_symm
    have := (hswap.hasSum.prod_fiberwise (fun n => (h₂ n).hasSum)).summable
    simpa only [← hcoeff] using this
  · calc ∑' n, PowerSeries.coeff n (logOf (1 + (P : ℚ_[p]⟦X⟧)))
        = ∑' n, ∑' d, G d n := tsum_congr hcoeff
      _ = ∑' d, ∑' n, G d n := hF.tsum_comm' h₁ h₂
      _ = ∑' d, coeff d (log ℚ_[p]) * (P.eval 1) ^ d := by
            refine tsum_congr fun d => ?_
            rw [hG]
            rw [tsum_mul_left, tsum_coeff_eq_eval_one, Polynomial.eval_pow]
      _ = logOneAdd (P.eval 1) := (logOneAdd_eq_tsum_coeff_log hEval).symm

/-! ### 部品7: ★p進対数の乗法性 -/

/-- `C a * X^k` の係数のノルムは `‖a‖` 以下。 -/
private theorem norm_coeff_CX_le (a : ℚ_[p]) (k n : ℕ) :
    ‖(Polynomial.C a * Polynomial.X ^ k).coeff n‖ ≤ ‖a‖ := by
  rw [Polynomial.coeff_C_mul, norm_mul, Polynomial.coeff_X_pow]
  split_ifs <;> simp

/-- ★★**p進対数の乗法性**。`‖x‖, ‖y‖ < 1` のとき
`log((1+x)(1+y)) = log(1+x) + log(1+y)`。

`(1+x)(1+y) = 1 + (x + y + xy)` なので、左辺は `logOneAdd (x + y + x*y)`。

証明の骨組み: 帳簿用の変数 `X` を1つ入れた `ℚ_[p]⟦X⟧` の中で
`(1 + C x·X)(1 + C y·X) = 1 + C(x+y)·X + C(xy)·X²` とし、
`PowerSeriesLog.logOf_mul`(形式的加法公式)を当ててから、
`summable_and_tsum_coeff_logOf`(評価の橋)で係数の総和を取る。
`x`, `y` は**係数**であって代入する点ではないので、
mathlib の評価 API(`IsLinearTopology` を要求し ℚ_p は満たさない)は一度も使わない。 -/
theorem logOneAdd_mul {x y : ℚ_[p]} (hx : ‖x‖ < 1) (hy : ‖y‖ < 1) :
    logOneAdd (x + y + x * y) = logOneAdd x + logOneAdd y := by
  classical
  set r : ℝ := max ‖x‖ ‖y‖ with hrdef
  have hr0 : 0 ≤ r := le_trans (norm_nonneg x) (le_max_left _ _)
  have hr1 : r < 1 := max_lt hx hy
  set Px : Polynomial ℚ_[p] := Polynomial.C x * Polynomial.X ^ 1 with hPx
  set Py : Polynomial ℚ_[p] := Polynomial.C y * Polynomial.X ^ 1 with hPy
  set Pz : Polynomial ℚ_[p] :=
    Polynomial.C (x + y) * Polynomial.X ^ 1 + Polynomial.C (x * y) * Polynomial.X ^ 2 with hPz
  -- 定数項
  have hPx0 : Px.coeff 0 = 0 := by simp [hPx]
  have hPy0 : Py.coeff 0 = 0 := by simp [hPy]
  have hPz0 : Pz.coeff 0 = 0 := by simp [hPz]
  -- 係数のノルム
  have hbx : ∀ n, ‖Px.coeff n‖ ≤ r := fun n =>
    le_trans (by rw [hPx]; exact norm_coeff_CX_le x 1 n) (le_max_left _ _)
  have hby : ∀ n, ‖Py.coeff n‖ ≤ r := fun n =>
    le_trans (by rw [hPy]; exact norm_coeff_CX_le y 1 n) (le_max_right _ _)
  have hbz : ∀ n, ‖Pz.coeff n‖ ≤ r := by
    intro n
    rw [hPz, Polynomial.coeff_add]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
    · exact le_trans (norm_coeff_CX_le (x + y) 1 n)
        (le_trans (IsUltrametricDist.norm_add_le_max x y) (le_of_eq hrdef.symm))
    · refine le_trans (norm_coeff_CX_le (x * y) 2 n) ?_
      rw [norm_mul]
      calc ‖x‖ * ‖y‖ ≤ ‖x‖ * 1 := by
            exact mul_le_mul_of_nonneg_left hy.le (norm_nonneg x)
        _ = ‖x‖ := mul_one _
        _ ≤ r := le_max_left _ _
  -- 多項式の等式
  have hPoly : Px + Py + Px * Py = Pz := by
    simp only [hPx, hPy, hPz, Polynomial.C_add, Polynomial.C_mul]
    ring
  have hPS : (1 + (Px : ℚ_[p]⟦X⟧)) * (1 + (Py : ℚ_[p]⟦X⟧)) = 1 + (Pz : ℚ_[p]⟦X⟧) := by
    rw [← hPoly]
    push_cast
    ring
  -- 定数項が 1
  have hcx : PowerSeries.constantCoeff (1 + (Px : ℚ_[p]⟦X⟧)) = 1 := by
    rw [map_add, map_one, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
      Polynomial.coeff_coe, hPx0, add_zero]
  have hcy : PowerSeries.constantCoeff (1 + (Py : ℚ_[p]⟦X⟧)) = 1 := by
    rw [map_add, map_one, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
      Polynomial.coeff_coe, hPy0, add_zero]
  -- 形式的加法公式
  have hform : logOf (1 + (Pz : ℚ_[p]⟦X⟧))
      = logOf (1 + (Px : ℚ_[p]⟦X⟧)) + logOf (1 + (Py : ℚ_[p]⟦X⟧)) := by
    rw [← hPS, logOf_mul hcx hcy]
  -- 橋を渡す
  obtain ⟨hsx, htx⟩ := summable_and_tsum_coeff_logOf Px hPx0 hr0 hr1 hbx
  obtain ⟨hsy, hty⟩ := summable_and_tsum_coeff_logOf Py hPy0 hr0 hr1 hby
  obtain ⟨_, htz⟩ := summable_and_tsum_coeff_logOf Pz hPz0 hr0 hr1 hbz
  have hze : Pz.eval 1 = x + y + x * y := by simp [hPz]
  have hxe : Px.eval 1 = x := by simp [hPx]
  have hye : Py.eval 1 = y := by simp [hPy]
  rw [← hze, ← htz, ← hxe, ← htx, ← hye, ← hty, hform]
  simp only [map_add]
  exact hsx.tsum_add hsy

/-! ### 部品8: 乗法性の帰結 —— log-shell は加法で閉じる -/

/-- ★`logShell` は**加法で閉じる**。

これは乗法性の直接の帰結であり、「`log` が乗法群 `1 + 𝔪` から加法群 `ℚ_p` への
準同型である」ことの集合レベルの言明。[AbsTopIII] p.3 が log-shell を
「canonical rigid **integral structure** on k」と呼ぶときの「structure」の一部である。 -/
theorem logShell_add_mem {a b : ℚ_[p]}
    (ha : a ∈ logShell (p := p)) (hb : b ∈ logShell (p := p)) :
    a + b ∈ logShell (p := p) := by
  obtain ⟨u, hu, rfl⟩ := ha
  obtain ⟨v, hv, rfl⟩ := hb
  refine ⟨u + v + u * v, ?_, logOneAdd_mul hu hv⟩
  simp only [Set.mem_setOf_eq] at hu hv ⊢
  refine lt_of_le_of_lt (IsUltrametricDist.norm_add_le_max _ _) (max_lt ?_ ?_)
  · exact lt_of_le_of_lt (IsUltrametricDist.norm_add_le_max u v) (max_lt hu hv)
  · rw [norm_mul]
    calc ‖u‖ * ‖v‖ ≤ ‖u‖ * 1 := mul_le_mul_of_nonneg_left hv.le (norm_nonneg u)
      _ = ‖u‖ := mul_one _
      _ < 1 := hu

end ABC3.Found.IUTchIII
