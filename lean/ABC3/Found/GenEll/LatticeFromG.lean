import ABC3.Found.GenEll.PexcODE

/-!
# GenEll 第 353 ブロック —— **★★★★★★★★`(g₂,g₃)` は束を決める**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★到達点——節点が閉じた

> **`g₂(L) = g₂(L')` かつ `g₃(L) = g₃(L')` ならば `Λ = Λ'`**(`lattice_eq_of_g₂_g₃_eq`)

★これで `Skeleton/GenEll/LatticeFromInvariants.lean` の `sorry` が消え、
第 350 の `archInv_congr` と合わせて **`archInv` は `(g₂,g₃)` で決まる**——
すなわち**曲線の関数**になる。

## ★★★★★★★見積もりより短くなった——**漸化式を立てずに済んだ**

スケルトンの段取りは「係数の漸化式 `(i−4)(i+3)a_i = 6Σa_ja_k − (g₂/2)[i=2]` を立てて
帰納法で `G` を決める」だった(見積もり 4-8 ブロック)。
★★★実際には**差を取ると漸化式が要らなくなる**:

`D := ℘[L−0] − ℘[L'−0]` と置くと、2 つの恒等式(第 352)の差から

    z²·D'' = (12 + 6z²(f+h))·D

★`D` の `0` での最小非零係数の番号を `m` とすると、両辺に `iteratedDeriv m` を当てて
`0` で評価する。**Leibniz の和が両側とも 1 項に潰れる**:

| 辺 | 潰れる理由 | 残る項 |
|---|---|---|
| 左 | `iteratedDeriv i (z²) 0` は `i = 2` 以外 0 | `C(m,2)·2·a_m` |
| 右 | `i ≥ 1` では `m−i < m` なので `a_{m−i} = 0` | `A(0)·a_m = 12·a_m` |

★★したがって `2·C(m,2) = 12`、すなわち **`C(m,2) = 6`、`m = 4`**。
★★★★しかし `(g₂,g₃)` が等しければ `a_0,…,a_4` は一致する
(`a_i = (i+1)!·G(i+2)`、`G` の奇数番は 0、`G₄ = g₂/60`、`G₆ = g₃/140`)ので **`m ≥ 5`**。
★★★★★矛盾。ゆえに `D ≡ 0`(`0` の近傍)。

★**`i = 4` で係数 `(i−4)(i+3)` が消える**という漸化式の特異性が、
ここでは「`m = 4` が唯一の可能性」という形で現れている——同じ事実の別の顔である。

## ★★★★★★束の一致まで

1. `℘ = ℘[L−0] + 1/z²`(第 351)なので `℘_L =ᶠ[𝓝 0] ℘_{L'}`。
2. `0` の近くの**束に属さない点** `z₀` を取り、**一致の定理**で
   `(Λ ∪ Λ')ᶜ`(可算集合の補集合なので連結)全体へ延ばす。
3. `x ∈ Λ` かつ `x ∉ Λ'` なら、`℘_{L'}` は `x` で解析的だから有界、
   しかし `℘_L` は `x` の近くで**非有界**(`eventually_norm_weierstrassP_gt`)。矛盾。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `iteratedDeriv_pexc` | ★★`a_n = (n+1)!·G(n+2)` |
| `iteratedDeriv_pexc_le_four` | ★★★`n ≤ 4` の係数は `(g₂,g₃)` で決まる |
| `iteratedDeriv_sq` | ★`z²` の高階微分 |
| `diff_ode` | ★★★★★**差の方程式 `z²D'' = (12+6z²(f+h))D`** |
| `key_coeff` | ★★★★★★**`C(m,2)·2·a_m = 12·a_m`** |
| `pexc_eventuallyEq` | ★★★★★★★`℘[L−0] =ᶠ ℘[L'−0]` |
| `eventually_norm_weierstrassP_gt` | ★★★★★`℘` は束の点の近くで非有界 |
| `eqOn_weierstrassP` | ★★★★★★`(Λ∪Λ')ᶜ` 上で `℘` が一致 |
| `lattice_eq_of_g₂_g₃_eq` | ★★★★★★★★**`(g₂,g₃)` は束を決める** |
-/

namespace ABC3.Found.GenEll

open Complex PeriodPair Filter Topology Nat

/-! ## ★★`℘[L−0]` の `0` での高階微分 -/

/-- ★★**`a_n = (n+1)!·G(n+2)`**(mathlib の冪級数展開)。 -/
theorem iteratedDeriv_pexc (L : PeriodPair) (n : ℕ) :
    iteratedDeriv n ℘[L - (0:ℂ)] 0
      = if n = 0 then 0 else ((n + 1)! : ℂ) * L.G (n + 2) := by
  rw [L.iteratedDeriv_weierstrassPExcept_self 0]
  split
  · simp [L.weierstrassPExcept_zero]
  · rw [L.sumInvPow_zero]

theorem G_four (L : PeriodPair) : L.G 4 = L.g₂ / 60 := by
  rw [PeriodPair.g₂]; ring

theorem G_six (L : PeriodPair) : L.G 6 = L.g₃ / 140 := by
  rw [PeriodPair.g₃]; ring

/-- ★★★`n ≤ 4` の係数は `(g₂,g₃)` で決まる——奇数番は消え、`G₄`・`G₆` が効く。 -/
theorem iteratedDeriv_pexc_le_four {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃)
    {n : ℕ} (hn : n ≤ 4) :
    iteratedDeriv n ℘[L - (0:ℂ)] 0 = iteratedDeriv n ℘[L' - (0:ℂ)] 0 := by
  rw [iteratedDeriv_pexc, iteratedDeriv_pexc]
  interval_cases n
  · simp
  · simp [L.G_eq_zero_of_odd 3 (by decide), L'.G_eq_zero_of_odd 3 (by decide)]
  · simp [G_four, h₂]
  · simp [L.G_eq_zero_of_odd 5 (by decide), L'.G_eq_zero_of_odd 5 (by decide)]
  · simp [G_six, h₃]

/-! ## ★`z²` の高階微分と反復微分の合成 -/

theorem iteratedDeriv_comp_iteratedDeriv (f : ℂ → ℂ) (m n : ℕ) :
    iteratedDeriv (n + m) f = iteratedDeriv m (iteratedDeriv n f) := by
  induction m with
  | zero => simp
  | succ k ih => rw [← add_assoc, iteratedDeriv_succ, ih, iteratedDeriv_succ]

theorem deriv_sq_fun : deriv (fun w : ℂ => w ^ 2) = fun w => 2 * w := by
  funext w; simp

theorem iteratedDeriv_three_sq : iteratedDeriv 3 (fun w : ℂ => w ^ 2) = 0 := by
  rw [show (3 : ℕ) = 1 + 1 + 1 from rfl, iteratedDeriv_succ, iteratedDeriv_succ,
    iteratedDeriv_one, deriv_sq_fun]
  have h : deriv (fun w : ℂ => 2 * w) = fun _ => 2 := by funext w; simp
  rw [h]
  funext w
  simp

/-- ★**`z²` の高階微分は `i = 2` のときだけ `2`**。 -/
theorem iteratedDeriv_sq (i : ℕ) :
    iteratedDeriv i (fun w : ℂ => w ^ 2) 0 = if i = 2 then 2 else 0 := by
  match i with
  | 0 => simp
  | 1 => simp
  | 2 =>
      rw [show (2 : ℕ) = 1 + 1 from rfl, iteratedDeriv_succ, iteratedDeriv_one, deriv_sq_fun]
      have h : deriv (fun w : ℂ => 2 * w) = fun _ => 2 := by funext w; simp
      rw [h]
      simp
  | (n + 3) =>
      rw [show n + 3 = 3 + n from by omega, iteratedDeriv_comp_iteratedDeriv,
        iteratedDeriv_three_sq]
      rw [if_neg (by omega)]
      simp

/-! ## ★★★★★差の方程式 -/

theorem eventually_analyticAt_pexc (L : PeriodPair) :
    ∀ᶠ w in 𝓝 (0:ℂ), AnalyticAt ℂ ℘[L - (0:ℂ)] w :=
  (isOpen_analyticAt ℂ _).eventually_mem (L.analyticAt_weierstrassPExcept 0)

theorem eventually_deriv_sub (L L' : PeriodPair) :
    (deriv (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]))
      =ᶠ[𝓝 (0:ℂ)] (fun w => deriv ℘[L - (0:ℂ)] w - deriv ℘[L' - (0:ℂ)] w) := by
  filter_upwards [eventually_analyticAt_pexc L, eventually_analyticAt_pexc L'] with w h1 h2
  exact deriv_sub h1.differentiableAt h2.differentiableAt

theorem eventually_deriv2_sub (L L' : PeriodPair) :
    ∀ᶠ w in 𝓝 (0:ℂ), deriv (deriv (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)])) w
      = deriv (deriv ℘[L - (0:ℂ)]) w - deriv (deriv ℘[L' - (0:ℂ)]) w := by
  filter_upwards [(eventually_deriv_sub L L').eventuallyEq_nhds,
    eventually_analyticAt_pexc L, eventually_analyticAt_pexc L'] with w hw h1 h2
  rw [hw.deriv_eq]
  exact deriv_sub h1.deriv.differentiableAt h2.deriv.differentiableAt

/-- ★★★★★**差 `D = f − h` の満たす方程式** `z²D'' = (12 + 6z²(f+h))·D`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★第 352 の 2 つの恒等式の差を取るだけである。 -/
theorem diff_ode {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) :
    ∀ᶠ w in 𝓝 (0:ℂ),
      w ^ 2 * iteratedDeriv 2 (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) w
        = (12 + 6 * w ^ 2 * (℘[L - (0:ℂ)] w + ℘[L' - (0:ℂ)] w))
          * (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) w := by
  filter_upwards [pexc_ode L, pexc_ode L', eventually_deriv2_sub L L'] with w hw1 hw2 hw3
  have h2 : iteratedDeriv 2 (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) w
      = deriv (deriv ℘[L - (0:ℂ)]) w - deriv (deriv ℘[L' - (0:ℂ)]) w := by
    rw [show (2:ℕ) = 1 + 1 from rfl, iteratedDeriv_succ, iteratedDeriv_one, hw3]
  rw [h2]
  show w ^ 2 * (deriv (deriv ℘[L - (0:ℂ)]) w - deriv (deriv ℘[L' - (0:ℂ)]) w) = _
  show _ = (12 + 6 * w ^ 2 * (℘[L - (0:ℂ)] w + ℘[L' - (0:ℂ)] w))
    * (℘[L - (0:ℂ)] w - ℘[L' - (0:ℂ)] w)
  rw [h₂] at hw1
  linear_combination hw1 - hw2

/-! ## ★★★★★★最小次数の係数を突き合わせる -/

/-- ★★★★★★**`C(m,2)·2·a_m = 12·a_m`**——Leibniz の和が両側とも 1 項に潰れる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★左辺は `iteratedDeriv i (z²) 0` が `i = 2` 以外で消えるから、
右辺は `i ≥ 1` で `m − i < m` だから潰れる。 -/
theorem key_coeff {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (m : ℕ) (hm2 : 2 ≤ m)
    (hmin : ∀ i < m, iteratedDeriv i (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) 0 = 0) :
    ((m.choose 2 : ℂ)) * 2 * iteratedDeriv m (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) 0
      = 12 * iteratedDeriv m (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) 0 := by
  set D : ℂ → ℂ := ℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)] with hD
  set A : ℂ → ℂ := fun w => 12 + 6 * w ^ 2 * (℘[L - (0:ℂ)] w + ℘[L' - (0:ℂ)] w) with hA
  have hDana : AnalyticAt ℂ D 0 :=
    (L.analyticAt_weierstrassPExcept 0).sub (L'.analyticAt_weierstrassPExcept 0)
  have hD2ana : AnalyticAt ℂ (iteratedDeriv 2 D) 0 := by
    rw [show (2:ℕ) = 1 + 1 from rfl, iteratedDeriv_succ, iteratedDeriv_one]
    exact hDana.deriv.deriv
  have hAana : AnalyticAt ℂ A 0 := by
    have h1 := L.analyticAt_weierstrassPExcept 0
    have h2 := L'.analyticAt_weierstrassPExcept 0
    exact analyticAt_const.add ((analyticAt_const.mul ((analyticAt_id).pow 2)).mul (h1.add h2))
  have hsq : ContDiffAt ℂ (m : ℕ∞) (fun w : ℂ => w ^ 2) 0 :=
    (contDiff_id.pow 2).contDiffAt.of_le (by exact_mod_cast le_top)
  have heq := Filter.EventuallyEq.iteratedDeriv_eq m (diff_ode (L := L) (L' := L') h₂)
  rw [iteratedDeriv_fun_mul hsq (hD2ana.contDiffAt),
    iteratedDeriv_fun_mul (hAana.contDiffAt) (hDana.contDiffAt)] at heq
  have hstep : iteratedDeriv (m - 2) (iteratedDeriv 2 D) 0 = iteratedDeriv m D 0 := by
    rw [← iteratedDeriv_comp_iteratedDeriv D (m - 2) 2, show 2 + (m - 2) = m from by omega]
  have hLHS : (∑ i ∈ Finset.range (m + 1), ((m.choose i : ℂ)) *
      iteratedDeriv i (fun w : ℂ => w ^ 2) 0 * iteratedDeriv (m - i) (iteratedDeriv 2 D) 0)
      = ((m.choose 2 : ℂ)) * 2 * iteratedDeriv m D 0 := by
    rw [Finset.sum_eq_single 2]
    · rw [iteratedDeriv_sq 2, if_pos rfl, hstep]
    · intro b _ hb
      rw [iteratedDeriv_sq b, if_neg hb]
      ring
    · intro hmem
      exact absurd (Finset.mem_range.2 (by omega)) hmem
  have hRHS : (∑ i ∈ Finset.range (m + 1), ((m.choose i : ℂ)) *
      iteratedDeriv i A 0 * iteratedDeriv (m - i) D 0) = 12 * iteratedDeriv m D 0 := by
    rw [Finset.sum_eq_single 0]
    · simp only [Nat.choose_zero_right, Nat.cast_one, one_mul, Nat.sub_zero, iteratedDeriv_zero]
      rw [hA]
      norm_num
    · intro b _ hb
      rw [hmin (m - b) (by omega)]
      ring
    · intro hmem
      exact absurd (Finset.mem_range.2 (by omega)) hmem
  rw [hLHS, hRHS] at heq
  exact heq

/-! ## ★★★★★★★`℘[L−0]` の一致 -/

theorem iteratedDeriv_diff_le_four {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃)
    {n : ℕ} (hn : n ≤ 4) : iteratedDeriv n (℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)]) 0 = 0 := by
  have hc1 : ContDiffAt ℂ (n : ℕ∞) ℘[L - (0:ℂ)] 0 := (L.analyticAt_weierstrassPExcept 0).contDiffAt
  have hc2 : ContDiffAt ℂ (n : ℕ∞) ℘[L' - (0:ℂ)] 0 :=
    (L'.analyticAt_weierstrassPExcept 0).contDiffAt
  have h := iteratedDeriv_fun_sub (n := n) (x := (0:ℂ)) hc1 hc2
  show iteratedDeriv n (fun w => ℘[L - (0:ℂ)] w - ℘[L' - (0:ℂ)] w) 0 = 0
  rw [h, iteratedDeriv_pexc_le_four h₂ h₃ hn, sub_self]

/-- ★★★★★★★**`(g₂,g₃)` が等しければ `℘[L−0]` は `0` の近傍で一致する**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★最小非零係数の番号 `m` は `C(m,2) = 6` を強いるので `m = 4` だが、
`(g₂,g₃)` が等しければ `m ≥ 5` である。 -/
theorem pexc_eventuallyEq {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) :
    ℘[L - (0:ℂ)] =ᶠ[𝓝 (0:ℂ)] ℘[L' - (0:ℂ)] := by
  classical
  set D : ℂ → ℂ := ℘[L - (0:ℂ)] - ℘[L' - (0:ℂ)] with hD
  have hDana : AnalyticAt ℂ D 0 :=
    (L.analyticAt_weierstrassPExcept 0).sub (L'.analyticAt_weierstrassPExcept 0)
  have hall : ∀ n : ℕ, iteratedDeriv n D 0 = 0 := by
    by_contra! hcon
    obtain ⟨n₀, hn₀⟩ := hcon
    have hex : ∃ n, iteratedDeriv n D 0 ≠ 0 := ⟨n₀, hn₀⟩
    set m := Nat.find hex with hm_def
    have hm : iteratedDeriv m D 0 ≠ 0 := Nat.find_spec hex
    have hmin : ∀ i < m, iteratedDeriv i D 0 = 0 := fun i hi => not_not.1 (Nat.find_min hex hi)
    have hm5 : 5 ≤ m := by
      by_contra! hlt
      exact hm (iteratedDeriv_diff_le_four h₂ h₃ (by omega))
    have hkey := key_coeff h₂ m (by omega) hmin
    have h6 : ((m.choose 2 : ℂ)) * 2 = 12 := mul_right_cancel₀ hm hkey
    have h6' : (m.choose 2 : ℕ) = 6 := by
      have hc : ((m.choose 2 : ℂ)) = 6 := by linear_combination h6 / 2
      exact_mod_cast hc
    have h10 : (5:ℕ).choose 2 ≤ m.choose 2 := Nat.choose_le_choose 2 hm5
    have h5 : (5:ℕ).choose 2 = 10 := by decide
    omega
  have htop : analyticOrderAt D 0 = ⊤ :=
    ENat.eq_top_iff_forall_ge.mpr fun n =>
      (natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hDana).mpr fun i _ => hall i
  have hD0 := analyticOrderAt_eq_top.1 htop
  filter_upwards [hD0] with z hz
  have hzz : ℘[L - (0:ℂ)] z - ℘[L' - (0:ℂ)] z = 0 := hz
  linear_combination hzz

/-! ## ★★★★★束の点の近くでの非有界性 -/

theorem tendsto_inv_norm_sq_at (x : ℂ) :
    Tendsto (fun z : ℂ => (‖z - x‖ ^ 2)⁻¹) (𝓝[≠] x) atTop := by
  refine Filter.Tendsto.inv_tendsto_nhdsGT_zero ?_
  refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
  · have h1 : Tendsto (fun z : ℂ => z - x) (𝓝 x) (𝓝 (x - x)) :=
      (continuous_id.sub continuous_const).tendsto x
    rw [sub_self] at h1
    have h : Tendsto (fun z : ℂ => ‖z - x‖ ^ 2) (𝓝 x) (𝓝 0) := by
      simpa using (h1.norm).pow 2
    exact h.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin] with z hz
    exact pow_pos (norm_pos_iff.2 (sub_ne_zero.2 hz)) 2

/-- ★★★★★**`℘` は束の点の近くで非有界**——`℘ = (解析的) + 1/(z−x)²`。 -/
theorem eventually_norm_weierstrassP_gt (L : PeriodPair) {x : ℂ} (hx : x ∈ L.lattice) (C : ℝ) :
    ∀ᶠ z in 𝓝[≠] x, C < ‖℘[L] z‖ := by
  set g : ℂ → ℂ := fun z => ℘[L - x] z - 1 / x ^ 2 with hg
  have hgc : ContinuousAt g x :=
    ((L.analyticAt_weierstrassPExcept x).continuousAt).sub continuousAt_const
  set M := ‖g x‖ with hM
  have e1 : ∀ᶠ z in 𝓝[≠] x, ‖g z‖ < M + 1 :=
    (hgc.norm.tendsto.eventually_lt_const (by linarith : M < M + 1)).filter_mono nhdsWithin_le_nhds
  have e2 := (tendsto_inv_norm_sq_at x).eventually_gt_atTop (C + M + 1)
  filter_upwards [e1, e2, self_mem_nhdsWithin] with z hz1 hz2 hz3
  have hzx : z - x ≠ 0 := sub_ne_zero.2 hz3
  have hsplit : ℘[L] z = g z + 1 / (z - x) ^ 2 := by
    have h := L.weierstrassPExcept_add ⟨x, hx⟩ z
    simp only [hg]
    linear_combination -h
  have hnorm : ‖(1 : ℂ) / (z - x) ^ 2‖ = (‖z - x‖ ^ 2)⁻¹ := by
    rw [norm_div, norm_one, norm_pow, one_div]
  have hge : (‖z - x‖ ^ 2)⁻¹ ≤ ‖℘[L] z‖ + ‖g z‖ := by
    rw [← hnorm, hsplit]
    calc ‖(1:ℂ)/(z-x)^2‖ = ‖(g z + 1/(z-x)^2) - g z‖ := by ring_nf
      _ ≤ ‖g z + 1/(z-x)^2‖ + ‖g z‖ := norm_sub_le _ _
  linarith

/-! ## ★★★★★★★★束の一致 -/

theorem weierstrassP_eventuallyEq {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) :
    ℘[L] =ᶠ[𝓝 (0:ℂ)] ℘[L'] := by
  filter_upwards [pexc_eventuallyEq h₂ h₃] with z hz
  rw [weierstrassP_eq_except_add L z, weierstrassP_eq_except_add L' z, hz]

theorem countable_lattice (L : PeriodPair) : (L.lattice : Set ℂ).Countable :=
  Set.countable_coe_iff.1 (countable_of_Lindelof_of_discrete (X := L.lattice))

/-- ★★★★★★`(g₂,g₃)` が等しければ `℘` は共通の定義域で一致する(一致の定理)。 -/
theorem eqOn_weierstrassP {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) :
    Set.EqOn ℘[L] ℘[L'] (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) := by
  obtain ⟨V, hVsub, hVopen, hV0⟩ := mem_nhds_iff.1 (weierstrassP_eventuallyEq h₂ h₃)
  have hpick : ∀ᶠ z in 𝓝[≠] (0:ℂ),
      z ∈ V ∧ z ∉ (L.lattice : Set ℂ) ∧ z ∉ (L'.lattice : Set ℂ) := by
    filter_upwards [nhdsWithin_le_nhds (hVopen.mem_nhds hV0),
      nhdsWithin_le_nhds (L.compl_lattice_sdiff_singleton_mem_nhds 0),
      nhdsWithin_le_nhds (L'.compl_lattice_sdiff_singleton_mem_nhds 0),
      self_mem_nhdsWithin] with z h1 h2 h3 h4
    exact ⟨h1, fun hm => h2 ⟨hm, h4⟩, fun hm => h3 ⟨hm, h4⟩⟩
  obtain ⟨z₀, hz₀V, hz₀L, hz₀L'⟩ := hpick.exists
  have hz₀U : z₀ ∈ (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) := by
    intro h
    rcases h with h | h
    · exact hz₀L h
    · exact hz₀L' h
  have heq : ℘[L] =ᶠ[𝓝 z₀] ℘[L'] := by
    filter_upwards [hVopen.mem_nhds hz₀V] with z hz using hVsub hz
  have hUconn : IsPreconnected (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) :=
    (Set.Countable.isConnected_compl_of_one_lt_rank (by simp)
      ((countable_lattice L).union (countable_lattice L'))).2
  have hana1 : AnalyticOnNhd ℂ ℘[L] (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) :=
    fun z hz => L.analyticOnNhd_weierstrassP z (fun h => hz (Or.inl h))
  have hana2 : AnalyticOnNhd ℂ ℘[L'] (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) :=
    fun z hz => L'.analyticOnNhd_weierstrassP z (fun h => hz (Or.inr h))
  exact hana1.eqOn_of_preconnected_of_eventuallyEq hana2 hUconn hz₀U heq

theorem lattice_subset_of_eqOn {L L' : PeriodPair}
    (heq : Set.EqOn ℘[L] ℘[L'] (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ)) :
    (L.lattice : Set ℂ) ⊆ (L'.lattice : Set ℂ) := by
  intro x hx
  by_contra hx'
  have hana : AnalyticAt ℂ ℘[L'] x := L'.analyticOnNhd_weierstrassP x hx'
  set K := ‖℘[L'] x‖ + 1 with hK
  have hbound : ∀ᶠ z in 𝓝 x, ‖℘[L'] z‖ < K :=
    hana.continuousAt.norm.tendsto.eventually_lt_const (by linarith : ‖℘[L'] x‖ < K)
  have hinU : ∀ᶠ z in 𝓝[≠] x, z ∈ (((L.lattice : Set ℂ) ∪ (L'.lattice : Set ℂ))ᶜ) := by
    filter_upwards [nhdsWithin_le_nhds (L.compl_lattice_sdiff_singleton_mem_nhds x),
      nhdsWithin_le_nhds (L'.isClosed_lattice.isOpen_compl.mem_nhds hx'),
      self_mem_nhdsWithin] with z h1 h2 h3
    rintro (hm | hm)
    · exact h1 ⟨hm, h3⟩
    · exact h2 hm
  have hbig := eventually_norm_weierstrassP_gt L hx K
  obtain ⟨z, hz1, hz2, hz3⟩ :=
    (hbig.and ((hbound.filter_mono nhdsWithin_le_nhds).and hinU)).exists
  rw [heq hz3] at hz1
  linarith

/-- ★★★★★★★★**`(g₂,g₃)` は束を決める**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★これで `archInv` は曲線の関数になる(第 350 `archInv_congr` と合わせて)。 -/
theorem lattice_eq_of_g₂_g₃_eq {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) :
    L.lattice = L'.lattice := by
  apply SetLike.coe_injective
  exact Set.Subset.antisymm (lattice_subset_of_eqOn (eqOn_weierstrassP h₂ h₃))
    (lattice_subset_of_eqOn (eqOn_weierstrassP h₂.symm h₃.symm))

/-- ★★★★★★★★**アルキメデス不変量は `(g₂,g₃)` で決まる**。 -/
theorem archInv_eq_of_g₂_g₃_eq {L L' : PeriodPair} (h₂ : L.g₂ = L'.g₂) (h₃ : L.g₃ = L'.g₃) :
    archInv L = archInv L' :=
  archInv_congr (lattice_eq_of_g₂_g₃_eq h₂ h₃)

/-! ## ★出典の紐付け(`.src`) -/

def lattice_eq_of_g₂_g₃_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def key_coeff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def pexc_eventuallyEq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
