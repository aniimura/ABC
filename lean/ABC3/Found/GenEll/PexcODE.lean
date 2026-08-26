import ABC3.Found.GenEll.WeierstrassODE

/-!
# GenEll 第 352 ブロック —— **★★★★★★`f = ℘[L−0]` の解析的な恒等式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★節点「`(g₂,g₃)` は束を決める」の第 2 歩

第 351 で `℘'' = 6℘² − g₂/2`(`Λᶜ` 上)を出した。
★本ブロックはこれを **`0` で解析的な `f := ℘[L−0]`** の言葉に書き直す:

> **`z²·f'' = 6z²·f² + 12·f − (g₂/2)·z²`**(`0` の近傍で、`pexc_ode`)

★★これが**係数の漸化式**を取り出す土台である
(`f` の `0` での冪級数係数は mathlib の `hasFPowerSeriesAt_weierstrassPExcept` で
`a_i = (i+1)·G(i+2)` と分かっている)。

## ★★★★★段取り

`℘ = f + z⁻²`(第 351 `weierstrassP_eq_except_add`)を 2 回微分する:

    ℘'  = f' − 2/z³        (`derivWeierstrassP_eq_pexc`)
    ℘'' = f'' + 6/z⁴       (`deriv_derivWeierstrassP_eq_pexc`)

★これを `℘'' = 6℘² − g₂/2` に入れると `6/z⁴` が両辺で打ち消し、
`f'' = 6f² + 12f/z² − g₂/2`、すなわち `z²f'' = 6z²f² + 12f − (g₂/2)z²`(`Λᶜ` 上)。

★★★**`0` を含む近傍へ延ばす段**は、値の突き合わせだけで済んだ:
両辺は `w = 0` でともに `12·f(0) = 0`(`weierstrassPExcept_zero`)なので、
`𝓝[≠] 0 ⊔ pure 0 = 𝓝 0`(`nhdsNE_sup_pure`)で貼り合わせればよい。
★連続性や一致の定理を持ち出す必要はなかった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hasDerivAt_inv_sq`・`hasDerivAt_inv_cube` | ★`1/z²`・`−2/z³` の微分 |
| `derivWeierstrassP_eq_pexc` | ★★★**`℘' = f' − 2/z³`** |
| `deriv_derivWeierstrassP_eq_pexc` | ★★★**`℘'' = f'' + 6/z⁴`** |
| `pexc_ode_of_notMem` | ★★★★★`Λᶜ` 上の恒等式 |
| `pexc_ode` | ★★★★★★**`0` の近傍での恒等式** |
-/

namespace ABC3.Found.GenEll

open Complex PeriodPair Filter Topology

/-! ## ★負冪の微分 -/

theorem hasDerivAt_inv_sq {w : ℂ} (hw : w ≠ 0) :
    HasDerivAt (fun t : ℂ => 1 / t ^ 2) (-2 / w ^ 3) w := by
  have h := hasDerivAt_zpow (-2 : ℤ) w (Or.inl hw)
  have he : (fun t : ℂ => t ^ (-2 : ℤ)) = fun t : ℂ => 1 / t ^ 2 := by
    funext t
    rw [zpow_neg, show ((2:ℤ)) = ((2:ℕ):ℤ) from rfl, zpow_natCast]
    ring
  have hval : ((-2 : ℤ) : ℂ) * w ^ ((-2 : ℤ) - 1) = -2 / w ^ 3 := by
    rw [show ((-2 : ℤ) - 1 : ℤ) = -((3 : ℕ) : ℤ) by norm_num, zpow_neg, zpow_natCast]
    push_cast
    ring
  rw [he, hval] at h
  exact h

theorem hasDerivAt_inv_cube {w : ℂ} (hw : w ≠ 0) :
    HasDerivAt (fun t : ℂ => -2 / t ^ 3) (6 / w ^ 4) w := by
  have h := hasDerivAt_zpow (-3 : ℤ) w (Or.inl hw)
  have he : (fun t : ℂ => ((-2 : ℂ)) * t ^ (-3 : ℤ)) = fun t : ℂ => -2 / t ^ 3 := by
    funext t
    rw [zpow_neg, show ((3:ℤ)) = ((3:ℕ):ℤ) from rfl, zpow_natCast]
    ring
  have h2 := h.const_mul (-2 : ℂ)
  rw [he] at h2
  have hval : (-2 : ℂ) * (((-3 : ℤ) : ℂ) * w ^ ((-3 : ℤ) - 1)) = 6 / w ^ 4 := by
    rw [show ((-3 : ℤ) - 1 : ℤ) = -((4 : ℕ) : ℤ) by norm_num, zpow_neg, zpow_natCast]
    push_cast
    ring
  rw [hval] at h2
  exact h2

/-! ## ★★★`℘` の微分を `f` で書く -/

theorem hasDerivAt_pexc (L : PeriodPair) {w : ℂ} (hw : w ∉ (L.lattice : Set ℂ) \ {0}) :
    HasDerivAt ℘[L - (0:ℂ)] (deriv ℘[L - (0:ℂ)] w) w :=
  ((L.analyticOnNhd_weierstrassPExcept 0) w hw).differentiableAt.hasDerivAt

theorem hasDerivAt_deriv_pexc (L : PeriodPair) {w : ℂ} (hw : w ∉ (L.lattice : Set ℂ) \ {0}) :
    HasDerivAt (deriv ℘[L - (0:ℂ)]) (deriv (deriv ℘[L - (0:ℂ)]) w) w :=
  ((AnalyticOnNhd.deriv (L.analyticOnNhd_weierstrassPExcept 0)) w hw).differentiableAt.hasDerivAt

/-- ★★★**`℘' = f' − 2/z³`**(`f := ℘[L−0]`)。 -/
theorem derivWeierstrassP_eq_pexc (L : PeriodPair) {w : ℂ} (hw : w ∉ L.lattice) :
    ℘'[L] w = deriv ℘[L - (0:ℂ)] w - 2 / w ^ 3 := by
  have hw0 : w ≠ 0 := fun h => hw (h ▸ zero_mem _)
  have hsd : w ∉ (L.lattice : Set ℂ) \ {0} := fun h => hw h.1
  have hsum : HasDerivAt (fun t : ℂ => ℘[L - (0:ℂ)] t + 1 / t ^ 2)
      (deriv ℘[L - (0:ℂ)] w + -2 / w ^ 3) w :=
    (hasDerivAt_pexc L hsd).add (hasDerivAt_inv_sq hw0)
  have hfun : (fun t : ℂ => ℘[L - (0:ℂ)] t + 1 / t ^ 2) = ℘[L] := by
    funext t; exact (weierstrassP_eq_except_add L t).symm
  rw [hfun] at hsum
  have h := (hasDerivAt_weierstrassP L hw).unique hsum
  rw [h]
  ring

/-- ★★★**`℘'' = f'' + 6/z⁴`**。 -/
theorem deriv_derivWeierstrassP_eq_pexc (L : PeriodPair) {w : ℂ} (hw : w ∉ L.lattice) :
    deriv ℘'[L] w = deriv (deriv ℘[L - (0:ℂ)]) w + 6 / w ^ 4 := by
  have hw0 : w ≠ 0 := fun h => hw (h ▸ zero_mem _)
  have hsd : w ∉ (L.lattice : Set ℂ) \ {0} := fun h => hw h.1
  have hopen := L.isClosed_lattice.isOpen_compl
  have heq : ℘'[L] =ᶠ[nhds w] (fun t => deriv ℘[L - (0:ℂ)] t + (-2 / t ^ 3)) := by
    filter_upwards [hopen.mem_nhds hw] with t ht
    rw [derivWeierstrassP_eq_pexc L ht]; ring
  rw [heq.deriv_eq]
  exact ((hasDerivAt_deriv_pexc L hsd).add (hasDerivAt_inv_cube hw0)).deriv

/-! ## ★★★★★★`f` の恒等式 -/

/-- ★★★★★**`z²f'' = 6z²f² + 12f − (g₂/2)z²`**(`Λᶜ` 上)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`℘'' = 6℘² − g₂/2` に `℘ = f + z⁻²` を入れると `6/z⁴` が打ち消す。 -/
theorem pexc_ode_of_notMem (L : PeriodPair) {w : ℂ} (hw : w ∉ L.lattice) :
    w ^ 2 * deriv (deriv ℘[L - (0:ℂ)]) w
      = 6 * w ^ 2 * ℘[L - (0:ℂ)] w ^ 2 + 12 * ℘[L - (0:ℂ)] w - L.g₂ / 2 * w ^ 2 := by
  have hw0 : w ≠ 0 := fun h => hw (h ▸ zero_mem _)
  have h1 := deriv_derivWeierstrassP L hw
  rw [deriv_derivWeierstrassP_eq_pexc L hw, weierstrassP_eq_except_add L w] at h1
  field_simp at h1
  have h2w : (2 : ℂ) * w ^ 2 ≠ 0 := by simp [hw0]
  refine mul_left_cancel₀ h2w ?_
  linear_combination h1

/-- ★★★★★★**`0` の近傍での恒等式**——係数の漸化式の土台。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`w = 0` では両辺とも `12·f(0) = 0` なので、`𝓝[≠]0 ⊔ pure 0 = 𝓝 0` で貼り合わせる。 -/
theorem pexc_ode (L : PeriodPair) :
    (fun w : ℂ => w ^ 2 * deriv (deriv ℘[L - (0:ℂ)]) w)
      =ᶠ[nhds 0] (fun w : ℂ => 6 * w ^ 2 * ℘[L - (0:ℂ)] w ^ 2 + 12 * ℘[L - (0:ℂ)] w
        - L.g₂ / 2 * w ^ 2) := by
  have hpunct : ∀ᶠ w in 𝓝[≠] (0:ℂ),
      w ^ 2 * deriv (deriv ℘[L - (0:ℂ)]) w
        = 6 * w ^ 2 * ℘[L - (0:ℂ)] w ^ 2 + 12 * ℘[L - (0:ℂ)] w - L.g₂ / 2 * w ^ 2 := by
    filter_upwards [nhdsWithin_le_nhds (L.compl_lattice_sdiff_singleton_mem_nhds 0),
      self_mem_nhdsWithin] with w hw hw0
    exact pexc_ode_of_notMem L (fun hmem => hw ⟨hmem, hw0⟩)
  have hzero : (0:ℂ) ^ 2 * deriv (deriv ℘[L - (0:ℂ)]) 0
      = 6 * (0:ℂ) ^ 2 * ℘[L - (0:ℂ)] 0 ^ 2 + 12 * ℘[L - (0:ℂ)] 0 - L.g₂ / 2 * (0:ℂ) ^ 2 := by
    simp [L.weierstrassPExcept_zero]
  show ∀ᶠ w in nhds (0:ℂ), _
  rw [← nhdsNE_sup_pure (0:ℂ), Filter.eventually_sup]
  exact ⟨hpunct, Filter.eventually_pure.2 hzero⟩

/-! ## ★出典の紐付け(`.src`) -/

def pexc_ode.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def derivWeierstrassP_eq_pexc.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def pexc_ode_of_notMem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
