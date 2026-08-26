import ABC3.Found.GenEll.LatticeInvariance

/-!
# GenEll 第 351 ブロック —— **★★★★★★`℘'' = 6℘² − g₂/2`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★節点「`(g₂,g₃)` は束を決める」の第 1 歩

`Skeleton/GenEll/LatticeFromInvariants.lean` の段取りは
**係数の漸化式**を `z²f'' = 6z²f² + 12f − (g₂/2)z²` から取り出すことだった。
★その出発点が本ブロックの **`℘'' = 6℘² − g₂/2`**(mathlib に無い)である。

## ★★★★★段取り——**微分して `℘'` で割る**

mathlib は `derivWeierstrassP_sq : ℘'² = 4℘³ − g₂℘ − g₃` を持つ。
★両辺を微分すると `2℘'·℘'' = 12℘²·℘' − g₂·℘'`、すなわち

    ℘' · (2℘'' − 12℘² + g₂) = 0   (`Λᶜ` 上)

★★`℘'` は `Λᶜ` 上で恒等的に 0 ではないので、**一致の定理**で第 2 因子が消える。

### ★★★★★★`℘'` が恒等的に 0 でないこと——**0 の近くで `℘` は非有界**

★もし `℘' ≡ 0` なら `℘` は `Λᶜ`(連結開集合)上で定数になる。
★★しかし `℘[L] z = ℘[L−0] z + 1/z²`(`weierstrassPExcept_add` を `l₀ = 0` で)であり、
`℘[L−0]` は `0` で解析的だから、`z → 0` で `‖℘[L] z‖ → ∞` である。
★★★定数は有界なので矛盾する。

★`Λᶜ` の連結性は mathlib の `Set.Countable.isConnected_compl_of_one_lt_rank`
(可算集合の補集合は連結)で出る——束は離散なので可算である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isPreconnected_compl_lattice` | ★★`Λᶜ` は連結 |
| `hasDerivAt_weierstrassP`・`hasDerivAt_derivWeierstrassP` | ★微分の形 |
| `weierstrassP_eq_except_add` | ★★**`℘ = ℘[L−0] + 1/z²`** |
| `tendsto_inv_norm_sq` | ★`‖z‖⁻² → ∞` |
| `deriv_relation` | ★★★★微分方程式を 1 回微分した形 |
| `exists_derivWeierstrassP_ne_zero` | ★★★★★**`℘'` は恒等的に 0 でない** |
| `deriv_derivWeierstrassP` | ★★★★★★**`℘'' = 6℘² − g₂/2`** |
-/

namespace ABC3.Found.GenEll

open Complex PeriodPair Filter Topology

/-! ## ★★`Λᶜ` の連結性 -/

/-- ★★**`Λᶜ` は連結**——束は離散なので可算であり、可算集合の補集合は連結。 -/
theorem isPreconnected_compl_lattice (L : PeriodPair) : IsPreconnected ((L.lattice : Set ℂ)ᶜ) :=
  (Set.Countable.isConnected_compl_of_one_lt_rank (by simp)
    (countable_of_Lindelof_of_discrete (X := L.lattice))).2

/-! ## ★微分の形 -/

theorem hasDerivAt_weierstrassP (L : PeriodPair) {z : ℂ} (hz : z ∉ L.lattice) :
    HasDerivAt ℘[L] (℘'[L] z) z := by
  have h := (L.analyticOnNhd_weierstrassP z hz).differentiableAt.hasDerivAt
  rwa [L.deriv_weierstrassP] at h

theorem hasDerivAt_derivWeierstrassP (L : PeriodPair) {z : ℂ} (hz : z ∉ L.lattice) :
    HasDerivAt ℘'[L] (deriv ℘'[L] z) z :=
  (L.analyticOnNhd_derivWeierstrassP z hz).differentiableAt.hasDerivAt

/-! ## ★★`℘` の主要部を切り出す -/

/-- ★★**`℘[L] z = ℘[L−0] z + 1/z²`**——`0` の項だけを取り出した形。 -/
theorem weierstrassP_eq_except_add (L : PeriodPair) (z : ℂ) :
    ℘[L] z = ℘[L - (0:ℂ)] z + 1 / z ^ 2 := by
  have h := L.weierstrassPExcept_add ⟨0, zero_mem _⟩ z
  simpa using h.symm

/-- ★`‖z‖⁻² → ∞`(`z → 0`、`z ≠ 0`)。 -/
theorem tendsto_inv_norm_sq : Tendsto (fun z : ℂ => (‖z‖ ^ 2)⁻¹) (𝓝[≠] (0:ℂ)) atTop := by
  refine Filter.Tendsto.inv_tendsto_nhdsGT_zero ?_
  refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
  · have h : Tendsto (fun z : ℂ => ‖z‖ ^ 2) (𝓝 (0:ℂ)) (𝓝 0) := by
      have := (continuous_norm (E := ℂ)).tendsto (0 : ℂ)
      simpa using this.pow 2
    exact h.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin] with z hz
    exact pow_pos (norm_pos_iff.2 hz) 2

/-! ## ★★★★微分方程式を 1 回微分する -/

/-- ★★★★**`2℘'℘'' = 12℘²℘' − g₂℘'`**——`℘'² = 4℘³ − g₂℘ − g₃` の微分。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★両辺は `Λᶜ`(開集合)で一致するので導関数も一致する。 -/
theorem deriv_relation (L : PeriodPair) {z : ℂ} (hz : z ∉ L.lattice) :
    2 * ℘'[L] z * deriv ℘'[L] z
      = 12 * ℘[L] z ^ 2 * ℘'[L] z - L.g₂ * ℘'[L] z := by
  have hopen := L.isClosed_lattice.isOpen_compl
  have hmem : (L.lattice : Set ℂ)ᶜ ∈ nhds z := hopen.mem_nhds hz
  have heq : (fun w => ℘'[L] w ^ 2)
      =ᶠ[nhds z] (fun w => 4 * ℘[L] w ^ 3 - L.g₂ * ℘[L] w - L.g₃) := by
    filter_upwards [hmem] with w hw using L.derivWeierstrassP_sq w hw
  have hL : HasDerivAt (fun w => ℘'[L] w ^ 2)
      ((2 : ℕ) * ℘'[L] z ^ (2 - 1) * deriv ℘'[L] z) z :=
    (hasDerivAt_derivWeierstrassP L hz).pow 2
  have hR : HasDerivAt (fun w => 4 * ℘[L] w ^ 3 - L.g₂ * ℘[L] w - L.g₃)
      (4 * ((3 : ℕ) * ℘[L] z ^ (3 - 1) * ℘'[L] z) - L.g₂ * ℘'[L] z - 0) z :=
    ((((hasDerivAt_weierstrassP L hz).pow 3).const_mul 4).sub
      ((hasDerivAt_weierstrassP L hz).const_mul L.g₂)).sub (hasDerivAt_const z L.g₃)
  have hu := (hL.congr_of_eventuallyEq heq.symm).unique hR
  push_cast at hu
  linear_combination hu

/-! ## ★★★★★`℘'` は恒等的に 0 でない -/

/-- ★★★★★**`℘'` は `Λᶜ` 上で恒等的に 0 ではない**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★もし恒等的に 0 なら `℘` は連結な `Λᶜ` 上で定数になるが、
`℘ = ℘[L−0] + 1/z²` は `0` の近くで非有界である。 -/
theorem exists_derivWeierstrassP_ne_zero (L : PeriodPair) :
    ∃ z, z ∉ L.lattice ∧ ℘'[L] z ≠ 0 := by
  by_contra! hcon
  have hopen := L.isClosed_lattice.isOpen_compl
  have hdiff : DifferentiableOn ℂ ℘[L] ((L.lattice : Set ℂ)ᶜ) := L.differentiableOn_weierstrassP
  have hfd : Set.EqOn (fderiv ℂ ℘[L]) 0 ((L.lattice : Set ℂ)ᶜ) := by
    intro z hz
    have h := ((L.analyticOnNhd_weierstrassP z hz).differentiableAt).hasDerivAt
    rw [L.deriv_weierstrassP, hcon z hz] at h
    have h1 : HasFDerivAt ℘[L] (0 : ℂ →L[ℂ] ℂ) z := by simpa using h.hasFDerivAt
    exact h1.fderiv
  obtain ⟨z₀, hz₀⟩ : ∃ z₀ : ℂ, z₀ ∉ L.lattice := ⟨L.ω₁ / 2, L.ω₁_div_two_notMem_lattice⟩
  set c := ℘[L] z₀ with hc
  have hconst : ∀ z, z ∉ L.lattice → ℘[L] z = c :=
    fun z hz => hopen.is_const_of_fderiv_eq_zero (isPreconnected_compl_lattice L) hdiff hfd hz hz₀
  set M := ‖℘[L - (0:ℂ)] 0‖ with hM
  have hf : ContinuousAt ℘[L - (0:ℂ)] 0 := (L.analyticAt_weierstrassPExcept 0).continuousAt
  have e1 : ∀ᶠ z in 𝓝[≠] (0:ℂ), ‖℘[L - (0:ℂ)] z‖ < M + 1 :=
    (hf.norm.tendsto.eventually_lt_const (by linarith : M < M + 1)).filter_mono nhdsWithin_le_nhds
  have e2 := tendsto_inv_norm_sq.eventually_gt_atTop (‖c‖ + M + 1)
  have e3 : ∀ᶠ z in 𝓝[≠] (0:ℂ), z ∉ L.lattice := by
    filter_upwards [nhdsWithin_le_nhds (L.compl_lattice_sdiff_singleton_mem_nhds 0),
      self_mem_nhdsWithin] with z hz hz0 hmem
    exact hz ⟨hmem, hz0⟩
  obtain ⟨z, hz1, hz2, hz3⟩ := (e1.and (e2.and e3)).exists
  have hzc := hconst z hz3
  have hsplit := weierstrassP_eq_except_add L z
  have hnorm : ‖(1 : ℂ) / z ^ 2‖ = (‖z‖ ^ 2)⁻¹ := by
    rw [norm_div, norm_one, norm_pow, one_div]
  have hle : (‖z‖ ^ 2)⁻¹ ≤ ‖℘[L] z‖ + ‖℘[L - (0:ℂ)] z‖ := by
    rw [← hnorm, hsplit]
    calc ‖(1:ℂ)/z^2‖ = ‖(℘[L - (0:ℂ)] z + 1/z^2) - ℘[L - (0:ℂ)] z‖ := by ring_nf
      _ ≤ ‖℘[L - (0:ℂ)] z + 1/z^2‖ + ‖℘[L - (0:ℂ)] z‖ := norm_sub_le _ _
  rw [hzc] at hle
  linarith

/-! ## ★★★★★★`℘'' = 6℘² − g₂/2` -/

/-- ★★★★★★**`℘'' = 6℘² − g₂/2`**(`Λᶜ` 上)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`℘'·(2℘'' − 12℘² + g₂) = 0` で、`℘'` が消えない点の近傍から一致の定理で延ばす。 -/
theorem deriv_derivWeierstrassP (L : PeriodPair) {z : ℂ} (hz : z ∉ L.lattice) :
    deriv ℘'[L] z = 6 * ℘[L] z ^ 2 - L.g₂ / 2 := by
  set S : ℂ → ℂ := fun w => 2 * deriv ℘'[L] w - 12 * ℘[L] w ^ 2 + L.g₂ with hS
  have hopen : IsOpen ((L.lattice : Set ℂ)ᶜ) := L.isClosed_lattice.isOpen_compl
  have hd : AnalyticOnNhd ℂ (deriv ℘'[L]) ((L.lattice : Set ℂ)ᶜ) :=
    AnalyticOnNhd.deriv (L.analyticOnNhd_derivWeierstrassP)
  have hSana : AnalyticOnNhd ℂ S ((L.lattice : Set ℂ)ᶜ) := by
    intro w hw
    have h1 : AnalyticAt ℂ (deriv ℘'[L]) w := hd w hw
    have h2 : AnalyticAt ℂ ℘[L] w := L.analyticOnNhd_weierstrassP w hw
    exact AnalyticAt.add (AnalyticAt.sub (analyticAt_const.mul h1)
      (analyticAt_const.mul (h2.pow 2))) analyticAt_const
  have hprod : ∀ w ∈ ((L.lattice : Set ℂ)ᶜ), ℘'[L] w * S w = 0 := by
    intro w hw
    have h := deriv_relation L (z := w) hw
    simp only [hS]
    linear_combination h
  obtain ⟨z₀, hz₀, hz₀ne⟩ := exists_derivWeierstrassP_ne_zero L
  have hnbhd : ∀ᶠ w in nhds z₀, ℘'[L] w ≠ 0 :=
    ContinuousAt.eventually_ne ((L.analyticOnNhd_derivWeierstrassP z₀ hz₀).continuousAt) hz₀ne
  have hSzero : S =ᶠ[nhds z₀] 0 := by
    filter_upwards [hnbhd, hopen.mem_nhds hz₀] with w hw hwU
    exact (mul_eq_zero.1 (hprod w hwU)).resolve_left hw
  have hEq := hSana.eqOn_zero_of_preconnected_of_eventuallyEq_zero
    (isPreconnected_compl_lattice L) hz₀ hSzero
  have h0 : S z = 0 := hEq hz
  simp only [hS] at h0
  linear_combination h0 / 2

/-! ## ★出典の紐付け(`.src`) -/

def deriv_derivWeierstrassP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def exists_derivWeierstrassP_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def deriv_relation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
