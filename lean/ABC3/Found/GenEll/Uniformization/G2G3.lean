/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Uniformization.Sublattice

/-!
# 一様化 —— `g₂`・`g₃` の変換・代数側との突き合わせ

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`g₂` の比較 -/

/-- ★★★★★★★★`℘[Λ − 0]` の原点での 2 階微分は `g₂/10`。

★mathlib の `iteratedDeriv_weierstrassPExcept_self`（`= 3!·sumInvPow 0 4 = 6·G₄`）
と `g₂ = 60 G₄` から。 -/
theorem iteratedDeriv_two_weierstrassPExcept (P : PeriodPair) :
    iteratedDeriv 2 (P.weierstrassPExcept 0) 0 = P.g₂ / 10 := by
  rw [P.iteratedDeriv_weierstrassPExcept_self 0 (n := 2)]
  have hG : P.sumInvPow 0 4 = P.G 4 := by rw [PeriodPair.sumInvPow_zero]
  simp only [if_neg (by decide : ¬(2 : ℕ) = 0), hG, PeriodPair.g₂]
  norm_num
  ring

/-- ★★★★★★★★`℘[Λ − 0]` の原点での 4 階微分は `6g₃/7`。

★`= 5!·sumInvPow 0 6 = 120·G₆` と `g₃ = 140 G₆` から。 -/
theorem iteratedDeriv_four_weierstrassPExcept (P : PeriodPair) :
    iteratedDeriv 4 (P.weierstrassPExcept 0) 0 = 6 * P.g₃ / 7 := by
  rw [P.iteratedDeriv_weierstrassPExcept_self 0 (n := 4)]
  have hG : P.sumInvPow 0 6 = P.G 6 := by rw [PeriodPair.sumInvPow_zero]
  simp only [if_neg (by decide : ¬(4 : ℕ) = 0), hG, PeriodPair.g₃]
  norm_num
  ring

/-- ★★★★★★★★`℘″ = 6℘² − g₂/2`——`iteratedDeriv` の言葉で。 -/
theorem iteratedDeriv_two_weierstrassP (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice) :
    iteratedDeriv 2 P.weierstrassP w = 6 * P.weierstrassP w ^ 2 - P.g₂ / 2 := by
  rw [iteratedDeriv_succ, iteratedDeriv_one, PeriodPair.deriv_weierstrassP]
  exact deriv_derivWeierstrassP P hw

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`g₂` の同種写像公式**

    `g₂(Λ′) = g₂(Λ) + 10·Σ_{w ∈ T∖{0}} (6·℘_Λ(w)² − g₂(Λ)/2)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★第 671 の等式

    `℘[Λ′−0](z) − ℘[Λ−0](z) = Σ_{w ∈ T∖{0}} (℘_Λ(z+w) − ℘_Λ(w))`

の両辺を原点で 2 回微分する。左辺は `g₂′/10 − g₂/10`（mathlib の Taylor 係数）、
右辺は `Σ ℘_Λ″(w) = Σ (6℘_Λ(w)² − g₂/2)`。

★★★★☆**これが Vélu の商 `E/H` の `a₄` の解析版である**——代数側
（`Found/GenEll/Velu.lean` の `veluQuotient`：`a₄ ↦ a₄ − 5v`）と突き合わせれば
`latticeCurve P′ = veluQuotient (latticeCurve P) H`、すなわち `α = 1` が出る。 -/
theorem g₂_isogeny (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) :
    P'.g₂ = P.g₂ + 10 * ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2) := by
  have hfun : (fun z => P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z)
      = fun z => ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) :=
    funext (weierstrassPExcept_sub_eq_sum P P' T h0T hvelu)
  have hL : iteratedDeriv 2
      (fun z => P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z) 0
      = P'.g₂ / 10 - P.g₂ / 10 := by
    rw [iteratedDeriv_fun_sub (P'.analyticAt_weierstrassPExcept 0).contDiffAt
      (P.analyticAt_weierstrassPExcept 0).contDiffAt,
      iteratedDeriv_two_weierstrassPExcept, iteratedDeriv_two_weierstrassPExcept]
  have hterm : ∀ w ∈ T.erase 0,
      iteratedDeriv 2 (fun z : ℂ => P.weierstrassP (z + w) - P.weierstrassP w) 0
        = 6 * P.weierstrassP w ^ 2 - P.g₂ / 2 := by
    intro w hw
    have hwn : w ∉ P.lattice := rep_notMem_lattice P P' T h0T hT hrep hw
    have hshift : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 :=
      shifted_analyticAt P 0 w (by rw [zero_add]; exact hwn)
    rw [iteratedDeriv_fun_sub hshift.contDiffAt contDiffAt_const]
    have h1 : iteratedDeriv 2 (fun z : ℂ => P.weierstrassP (z + w)) 0
        = iteratedDeriv 2 P.weierstrassP w := by
      rw [iteratedDeriv_comp_add_const]
      simp
    rw [h1, iteratedDeriv_two_weierstrassP P hwn,
      show iteratedDeriv 2 (fun _ : ℂ => P.weierstrassP w) 0 = 0 by
        rw [iteratedDeriv_succ, iteratedDeriv_one]
        simp, sub_zero]
  have hR : iteratedDeriv 2
      (fun z => ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w)) 0
      = ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2) := by
    rw [iteratedDeriv_fun_sum (fun w hw => by
      have hwn : w ∉ P.lattice := rep_notMem_lattice P P' T h0T hT hrep hw
      have hshift : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 :=
        shifted_analyticAt P 0 w (by rw [zero_add]; exact hwn)
      exact hshift.contDiffAt.sub contDiffAt_const)]
    exact Finset.sum_congr rfl hterm
  rw [hfun, hR] at hL
  linear_combination (-10 : ℂ) * hL

def g₂_isogeny.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(g₂ の同種写像公式——Vélu の商の a₄ の解析版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def iteratedDeriv_two_weierstrassPExcept.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘[Λ−0] の原点での 2 階微分は g₂/10。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`g₃` の比較 -/

/-- ★★★★★★`℘‴ = 12·℘·℘′`。 -/
theorem iteratedDeriv_three_weierstrassP (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice) :
    iteratedDeriv 3 P.weierstrassP w
      = 12 * P.weierstrassP w * P.derivWeierstrassP w := by
  have hopen : IsOpen ((P.lattice : Set ℂ)ᶜ) := P.isClosed_lattice.isOpen_compl
  have hiter2 : iteratedDeriv 2 P.weierstrassP = deriv P.derivWeierstrassP := by
    rw [iteratedDeriv_succ, iteratedDeriv_one, PeriodPair.deriv_weierstrassP]
  have heq : deriv P.derivWeierstrassP
      =ᶠ[nhds w] fun z => 6 * P.weierstrassP z ^ 2 - P.g₂ / 2 := by
    filter_upwards [hopen.mem_nhds hw] with z hz
    exact deriv_derivWeierstrassP P hz
  rw [iteratedDeriv_succ, hiter2, heq.deriv_eq]
  have h1 := hasDerivAt_weierstrassP P hw
  have hb := ((h1.pow 2).const_mul (6 : ℂ)).sub_const (P.g₂ / 2)
  have h2 : HasDerivAt (fun z : ℂ => 6 * P.weierstrassP z ^ 2 - P.g₂ / 2) _ w :=
    hb.congr_of_eventuallyEq (by filter_upwards with z; simp only [Pi.pow_apply])
  rw [h2.deriv]
  push_cast
  ring

/-- ★★★★★★★★`℘⁗ = 120℘³ − 18g₂℘ − 12g₃`。

★`℘⁗ = 12(℘′² + ℘·℘″) = 12(4℘³ − g₂℘ − g₃ + 6℘³ − g₂℘/2)`。 -/
theorem iteratedDeriv_four_weierstrassP (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice) :
    iteratedDeriv 4 P.weierstrassP w
      = 120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃ := by
  have hopen : IsOpen ((P.lattice : Set ℂ)ᶜ) := P.isClosed_lattice.isOpen_compl
  have heq : iteratedDeriv 3 P.weierstrassP
      =ᶠ[nhds w] fun z => 12 * P.weierstrassP z * P.derivWeierstrassP z := by
    filter_upwards [hopen.mem_nhds hw] with z hz
    exact iteratedDeriv_three_weierstrassP P hz
  rw [iteratedDeriv_succ, heq.deriv_eq]
  have h1 := hasDerivAt_weierstrassP P hw
  have h2 : HasDerivAt P.derivWeierstrassP (6 * P.weierstrassP w ^ 2 - P.g₂ / 2) w := by
    have h := hasDerivAt_derivWeierstrassP P hw
    rwa [deriv_derivWeierstrassP P hw] at h
  have h3 := ((h1.mul h2).const_mul (12 : ℂ))
  have h4 : HasDerivAt (fun z => 12 * P.weierstrassP z * P.derivWeierstrassP z)
      (12 * (P.derivWeierstrassP w * P.derivWeierstrassP w
        + P.weierstrassP w * (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))) w := by
    refine h3.congr_of_eventuallyEq ?_
    filter_upwards with z
    simp only [Pi.mul_apply]
    ring
  rw [h4.deriv]
  have hsq := P.derivWeierstrassP_sq w hw
  linear_combination 12 * hsq

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`g₃` の同種写像公式**

    `g₃(Λ′) = g₃(Λ) + (7/6)·Σ_{w ∈ T∖{0}} (120℘_Λ(w)³ − 18g₂℘_Λ(w) − 12g₃)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★第 671 の等式の両辺を原点で 4 回微分する。左辺は `6g₃′/7 − 6g₃/7`
（mathlib の Taylor 係数 `5!·sumInvPow 0 6 = 120 G₆`、`g₃ = 140 G₆`）、
右辺は `Σ ℘_Λ⁗(w)`。

★★★★☆**第 672 の `g₂` と合わせて、`latticeCurve P′` の係数が
`latticeCurve P` と代表系だけで決まった**——これを代数側の `veluQuotient`
（`Found/GenEll/Velu.lean`）と突き合わせれば `α = 1` が出る。 -/
theorem g₃_isogeny (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) :
    P'.g₃ = P.g₃ + (7 / 6) * ∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃) := by
  have hfun : (fun z => P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z)
      = fun z => ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) :=
    funext (weierstrassPExcept_sub_eq_sum P P' T h0T hvelu)
  have hconst : ∀ c : ℂ, iteratedDeriv 4 (fun _ : ℂ => c) 0 = 0 := by
    intro c
    rw [iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_succ, iteratedDeriv_one]
    simp
  have hL : iteratedDeriv 4
      (fun z => P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z) 0
      = 6 * P'.g₃ / 7 - 6 * P.g₃ / 7 := by
    rw [iteratedDeriv_fun_sub (P'.analyticAt_weierstrassPExcept 0).contDiffAt
      (P.analyticAt_weierstrassPExcept 0).contDiffAt,
      iteratedDeriv_four_weierstrassPExcept, iteratedDeriv_four_weierstrassPExcept]
  have hterm : ∀ w ∈ T.erase 0,
      iteratedDeriv 4 (fun z : ℂ => P.weierstrassP (z + w) - P.weierstrassP w) 0
        = 120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃ := by
    intro w hw
    have hwn : w ∉ P.lattice := rep_notMem_lattice P P' T h0T hT hrep hw
    have hshift : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 :=
      shifted_analyticAt P 0 w (by rw [zero_add]; exact hwn)
    rw [iteratedDeriv_fun_sub hshift.contDiffAt contDiffAt_const]
    have h1 : iteratedDeriv 4 (fun z : ℂ => P.weierstrassP (z + w)) 0
        = iteratedDeriv 4 P.weierstrassP w := by
      rw [iteratedDeriv_comp_add_const]
      simp
    rw [h1, iteratedDeriv_four_weierstrassP P hwn, hconst, sub_zero]
  have hR : iteratedDeriv 4
      (fun z => ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w)) 0
      = ∑ w ∈ T.erase 0,
        (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃) := by
    rw [iteratedDeriv_fun_sum (fun w hw => by
      have hwn : w ∉ P.lattice := rep_notMem_lattice P P' T h0T hT hrep hw
      have hshift : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 :=
        shifted_analyticAt P 0 w (by rw [zero_add]; exact hwn)
      exact hshift.contDiffAt.sub contDiffAt_const)]
    exact Finset.sum_congr rfl hterm
  rw [hfun, hR] at hL
  linear_combination (-7 / 6 : ℂ) * hL

def iteratedDeriv_four_weierstrassP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘⁗ = 120℘³ − 18g₂℘ − 12g₃。★無条件)",
    sectionId := "genell-lemma-3-5" }

def g₃_isogeny.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(g₃ の同種写像公式——Vélu の商の a₆ の解析版。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（解析側・完全形）——位数 `l` の点から `E/H` の係数まで**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、次がすべて同時に取れる:

* `z₀`（`Φ(z₀) = Q`）と周期対 `P′`（格子は `Λ′ = Λ + ℤz₀`）
* 整数 `A, B, C, D` で `ω₁ = Aω₁′ + Bω₂′`・`ω₂ = Cω₁′ + Dω₂′`・`|AD − BC| = l`
* 代表系 `T`（`|T| = l`、`0 ∈ T`）と Vélu の公式
* **`E/H` の係数**:

      g₂(Λ′) = g₂(Λ) + 10·Σ_{w∈T∖0} (6℘(w)² − g₂/2)
      g₃(Λ′) = g₃(Λ) + (7/6)·Σ_{w∈T∖0} (120℘(w)³ − 18g₂℘(w) − 12g₃)

★★★★★★☆**`latticeCurve P′` の係数が `latticeCurve P` と代表系だけで決まった。**
☆残るのは、これを代数側の Vélu の商（`Found/GenEll/Velu.lean` の `veluQuotient`）
と突き合わせて `α = 1`（`u′ = u`）を出すことである。 -/
theorem exists_isogeny_data_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ) (T : Finset ℂ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (0 : ℂ) ∈ T ∧ T.card = l ∧
      (∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) ∧
      (∑ w ∈ T.erase 0, P.derivWeierstrassP w = 0) ∧
      P'.g₂ = P.g₂ + 10 * ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2) ∧
      P'.g₃ = P.g₃ + (7 / 6) * ∑ w ∈ T.erase 0,
        (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃) := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  exact ⟨z₀, P', A, B, C, D, T, hz₀, h1, h2, hdet, hP', h0T, hcard, hvelu,
    sum_derivWeierstrassP_rep_eq_zero P P' T h0T hT hrep,
    g₂_isogeny P P' T h0T hT hrep hvelu,
    g₃_isogeny P P' T h0T hT hrep hvelu⟩

def exists_isogeny_data_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側・完全形——位数 l の点から E/H の係数まで。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★代数側との突き合わせ -/

/-- ★★★★★★★★★★**Vélu の `v_Q` は `℘″(w)`**。

`latticeCurve` では `a₁ = a₂ = a₃ = 0`・`a₄ = −g₂/4` なので

    v_Q = 2·g^x_Q = 2(3℘(w)² − g₂/4) = 6℘(w)² − g₂/2 = ℘″(w) -/
theorem veluV_latticePoint (P : PeriodPair) (w : ℂ) :
    veluV (latticeCurve P) (latticePointX P w) (latticePointY P w)
      = 6 * P.weierstrassP w ^ 2 - P.g₂ / 2 := by
  simp only [veluV, veluGx, veluGy, latticeCurve, latticePointX, latticePointY]
  ring

/-- ★★★★★★★★★★**Vélu の `w_Q`**——`u_Q = ℘′(w)²` と微分方程式から

    w_Q = ℘′(w)² + ℘″(w)·℘(w) = 10℘(w)³ − (3/2)g₂℘(w) − g₃ -/
theorem veluW_latticePoint (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice) :
    veluW (latticeCurve P) (latticePointX P w) (latticePointY P w)
      = 10 * P.weierstrassP w ^ 3 - (3 / 2) * P.g₂ * P.weierstrassP w - P.g₃ := by
  have hsq := P.derivWeierstrassP_sq w hw
  simp only [veluW, veluU, veluV, veluGx, veluGy, latticeCurve, latticePointX,
    latticePointY]
  linear_combination hsq

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**解析側の `Λ′` と代数側の Vélu の商が一致する**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`S` が `(H∖{O})/±` の代表系であること（＝仮説 `hvS`・`hwS`：`T∖{0}` にわたる和が
`S` にわたる和の 2 倍であること）を認めれば

    latticeCurve Λ′ = veluQuotient (latticeCurve Λ) S

★★★数値の照合:

* `a₄`: `g₂′ = g₂ + 10·Σ_{T∖0}(6℘²−g₂/2) = g₂ + 20v` ⟺ `−g₂′/4 = −g₂/4 − 5v` ✓
* `a₆`: `g₃′ = g₃ + (7/6)·Σ_{T∖0}(120℘³−18g₂℘−12g₃) = g₃ + 28w`
  ⟺ `−g₃′/4 = −g₃/4 − 7w` ✓（`b₂ = 0`）

★★★★★★☆**これで `latticeCurve P′` が `E/H` の Weierstrass モデルそのもので
あることが確定した**——変数変換は要らない、すなわち `α = 1`。 -/
theorem latticeCurve_eq_veluQuotient (P P' : PeriodPair) (T : Finset ℂ)
    (S : Finset (ℂ × ℂ))
    (hg₂ : P'.g₂ = P.g₂ + 10 * ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))
    (hg₃ : P'.g₃ = P.g₃ + (7 / 6) * ∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃))
    (hvS : (2 : ℂ) * veluVSum (latticeCurve P) S
      = ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))
    (hwS : (2 : ℂ) * veluWSum (latticeCurve P) S
      = ∑ w ∈ T.erase 0,
        (10 * P.weierstrassP w ^ 3 - (3 / 2) * P.g₂ * P.weierstrassP w - P.g₃)) :
    latticeCurve P' = veluQuotient (latticeCurve P) S := by
  have hsum : (∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃))
      = 12 * ∑ w ∈ T.erase 0,
        (10 * P.weierstrassP w ^ 3 - (3 / 2) * P.g₂ * P.weierstrassP w - P.g₃) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun w _ => by ring)
  have hb₂ : (latticeCurve P).b₂ = 0 := by
    simp [latticeCurve, WeierstrassCurve.b₂]
  have ha₄ : -P'.g₂ / 4 = (latticeCurve P).a₄ - 5 * veluVSum (latticeCurve P) S := by
    show -P'.g₂ / 4 = -P.g₂ / 4 - 5 * veluVSum (latticeCurve P) S
    rw [hg₂]
    linear_combination (5 / 2 : ℂ) * hvS
  have ha₆ : -P'.g₃ / 4
      = (latticeCurve P).a₆ - (latticeCurve P).b₂ * veluVSum (latticeCurve P) S
        - 7 * veluWSum (latticeCurve P) S := by
    rw [hb₂]
    show -P'.g₃ / 4
      = -P.g₃ / 4 - 0 * veluVSum (latticeCurve P) S - 7 * veluWSum (latticeCurve P) S
    rw [hg₃, hsum]
    linear_combination (7 / 2 : ℂ) * hwS
  have hP'eq : latticeCurve P' = ⟨0, 0, 0, -P'.g₂ / 4, -P'.g₃ / 4⟩ := rfl
  simp only [veluQuotient, veluCurve]
  rw [hP'eq, WeierstrassCurve.mk.injEq]
  exact ⟨rfl, rfl, rfl, ha₄, ha₆⟩

def veluV_latticePoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の v_Q は ℘″(w)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluW_latticePoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の w_Q = 10℘³ − (3/2)g₂℘ − g₃。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeCurve_eq_veluQuotient.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側の Λ′ と代数側の Vélu の商が一致する——α = 1)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**解析側の `Λ′` は代数側の Vélu の商そのもの——`±` 代表系を使わない形**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★Vélu の `v = Σ_S v_Q`（`S` は `(H∖{O})/±` の代表系、`v_Q = 2g^x_Q`）は、
`℘` が偶なので `H∖{O}` 全体にわたる `g^x_Q` の和に等しい:

    v = Σ_S 2g^x_Q = Σ_{H∖{O}} g^x_Q

同様に `w = Σ_S (u_Q + v_Q x_Q) = Σ_{H∖{O}} (u_Q/2 + g^x_Q·x_Q)`。
★★★★☆**これで代表系を選ぶ必要がなくなり、仮説なしで**

    latticeCurve Λ′ = veluCurve (latticeCurve Λ) v w

**が書ける。すなわち `α = 1`。** -/
theorem latticeCurve_eq_veluCurve (P P' : PeriodPair) (T : Finset ℂ)
    (hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice)
    (hg₂ : P'.g₂ = P.g₂ + 10 * ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))
    (hg₃ : P'.g₃ = P.g₃ + (7 / 6) * ∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃)) :
    latticeCurve P' = veluCurve (latticeCurve P)
      (∑ w ∈ T.erase 0,
        veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w))
      (∑ w ∈ T.erase 0,
        (veluU (latticeCurve P) (latticePointX P w) (latticePointY P w) / 2
          + veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w)
            * latticePointX P w)) := by
  have hVsum : (∑ w ∈ T.erase 0,
      veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w))
      = ∑ w ∈ T.erase 0, (3 * P.weierstrassP w ^ 2 - P.g₂ / 4) :=
    Finset.sum_congr rfl (fun w _ => by
      simp only [veluV2, veluGx, latticeCurve, latticePointX, latticePointY]
      ring)
  have hWsum : (∑ w ∈ T.erase 0,
      (veluU (latticeCurve P) (latticePointX P w) (latticePointY P w) / 2
        + veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w)
          * latticePointX P w))
      = ∑ w ∈ T.erase 0,
        (5 * P.weierstrassP w ^ 3 - (3 / 4) * P.g₂ * P.weierstrassP w - P.g₃ / 2) :=
    Finset.sum_congr rfl (fun w hw => by
      have hsq := P.derivWeierstrassP_sq w (hnot w hw)
      simp only [veluU, veluV2, veluGx, veluGy, latticeCurve, latticePointX,
        latticePointY]
      linear_combination hsq / 2)
  have hsum₂ : (∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))
      = 2 * ∑ w ∈ T.erase 0, (3 * P.weierstrassP w ^ 2 - P.g₂ / 4) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun w _ => by ring)
  have hsum₃ : (∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃))
      = 24 * ∑ w ∈ T.erase 0,
        (5 * P.weierstrassP w ^ 3 - (3 / 4) * P.g₂ * P.weierstrassP w - P.g₃ / 2) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl (fun w _ => by ring)
  have hb₂ : (latticeCurve P).b₂ = 0 := by
    simp [latticeCurve, WeierstrassCurve.b₂]
  have ha₄ : -P'.g₂ / 4 = (latticeCurve P).a₄ - 5 * ∑ w ∈ T.erase 0,
      veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w) := by
    rw [hVsum]
    show -P'.g₂ / 4
      = -P.g₂ / 4 - 5 * ∑ w ∈ T.erase 0, (3 * P.weierstrassP w ^ 2 - P.g₂ / 4)
    rw [hg₂, hsum₂]
    ring
  have ha₆ : -P'.g₃ / 4
      = (latticeCurve P).a₆ - (latticeCurve P).b₂ * (∑ w ∈ T.erase 0,
          veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w))
        - 7 * ∑ w ∈ T.erase 0,
          (veluU (latticeCurve P) (latticePointX P w) (latticePointY P w) / 2
            + veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w)
              * latticePointX P w) := by
    rw [hb₂, hWsum]
    show -P'.g₃ / 4
      = -P.g₃ / 4 - 0 * _ - 7 * ∑ w ∈ T.erase 0,
        (5 * P.weierstrassP w ^ 3 - (3 / 4) * P.g₂ * P.weierstrassP w - P.g₃ / 2)
    rw [hg₃, hsum₃]
    ring
  have hP'eq : latticeCurve P' = ⟨0, 0, 0, -P'.g₂ / 4, -P'.g₃ / 4⟩ := rfl
  simp only [veluCurve]
  rw [hP'eq, WeierstrassCurve.mk.injEq]
  exact ⟨rfl, rfl, rfl, ha₄, ha₆⟩

def latticeCurve_eq_veluCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側の Λ′ は代数側の Vélu の商そのもの——± 代表系なし。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（解析側・最終形）——位数 `l` の点から `E/H` の Weierstrass モデルまで**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ) = latticeCurve Λ` の位数ちょうど `l` の点 `Q` に対し、次が**仮説なしで**取れる:

* `z₀`（`Φ(z₀) = Q`）と周期対 `P′`（格子は `Λ′ = Λ + ℤz₀`）
* 整数 `A, B, C, D` で `ω₁ = Aω₁′ + Bω₂′`・`ω₂ = Cω₁′ + Dω₂′`・`|AD − BC| = l`
* 代表系 `T`（`|T| = l`）と、**`latticeCurve Λ′` が Vélu の商そのものであること**:

      latticeCurve Λ′ = veluCurve (latticeCurve Λ) v w
      v = Σ_{w ∈ T∖{0}} g^x_{Φ(w)},   w = Σ_{w ∈ T∖{0}} (u_{Φ(w)}/2 + g^x_{Φ(w)}·x_{Φ(w)})

★★★★★★★☆**変数変換は要らない。すなわち `α = 1`。**
☆`htFalt_isogeny_le_of_analytic_minimal`（第 617）の残る仮説 `α`・`hu` は
これで `α = 1`・`u′ = u` として満たせる。 -/
theorem exists_velu_model_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ) (T : Finset ℂ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (0 : ℂ) ∈ T ∧ T.card = l ∧
      latticeCurve P' = veluCurve (latticeCurve P)
        (∑ w ∈ T.erase 0,
          veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w))
        (∑ w ∈ T.erase 0,
          (veluU (latticeCurve P) (latticePointX P w) (latticePointY P w) / 2
            + veluV2 (latticeCurve P) (latticePointX P w) (latticePointY P w)
              * latticePointX P w)) := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  refine ⟨z₀, P', A, B, C, D, T, hz₀, h1, h2, hdet, hP', h0T, hcard, ?_⟩
  exact latticeCurve_eq_veluCurve P P' T
    (fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw)
    (g₂_isogeny P P' T h0T hT hrep hvelu)
    (g₃_isogeny P P' T h0T hT hrep hvelu)

def exists_velu_model_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側・最終形——位数 l の点から E/H の Weierstrass モデルまで。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★**点が一致すれば `Λ` を法として一致**——第 624 の言い換え。 -/
theorem sub_mem_lattice_of_point_eq (P : PeriodPair) {w v : ℂ} (hw : w ∉ P.lattice)
    (hv : v ∉ P.lattice) (hx : latticePointX P w = latticePointX P v)
    (hy : latticePointY P w = latticePointY P v) : w - v ∈ P.lattice := by
  refine mem_lattice_of_shift_eq P (w - v) hv ?_ ?_ ?_
  · rw [show v + (w - v) = w by ring]; exact hw
  · rw [show v + (w - v) = w by ring]; exact hx
  · rw [show v + (w - v) = w by ring]
    simp only [latticePointY] at hy
    linear_combination 2 * hy

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`latticeCurve Λ′` は `H∖{O}` 全体で書いた Vélu の商そのもの**

    latticeCurve Λ′ = veluQuotientFull (latticeCurve Λ) S,
    S = { (℘(w), ℘′(w)/2) : w ∈ T∖{0} }

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`w ↦ (℘(w), ℘′(w)/2)` は `T∖{0}` 上で単射（第 624 の単射性 ＋
代表系の元が `Λ` を法として相異なること）なので、`T∖{0}` にわたる和は
`S` にわたる和に等しい。

★★★★☆**これで解析側の結論が代数側の語彙（`veluQuotientFull`）で書けた**——
`Found/GenEll/Velu.lean` の `veluQuotientFull_map`（第 679）と合わせれば
`L` 上の `E/H` と各 `σ` での一意化が結びつく。 -/
theorem latticeCurve_eq_veluQuotientFull (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hg₂ : P'.g₂ = P.g₂ + 10 * ∑ w ∈ T.erase 0, (6 * P.weierstrassP w ^ 2 - P.g₂ / 2))
    (hg₃ : P'.g₃ = P.g₃ + (7 / 6) * ∑ w ∈ T.erase 0,
      (120 * P.weierstrassP w ^ 3 - 18 * P.g₂ * P.weierstrassP w - 12 * P.g₃)) :
    latticeCurve P' = veluQuotientFull (latticeCurve P)
      ((T.erase 0).image (fun w => (latticePointX P w, latticePointY P w))) := by
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  rw [veluQuotientFull, veluVFull, veluWFull,
    Finset.sum_image hinj, Finset.sum_image hinj]
  exact latticeCurve_eq_veluCurve P P' T hnot hg₂ hg₃

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（解析側・代数語彙）——位数 `l` の点から `E/H` まで**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ) = latticeCurve Λ` の位数ちょうど `l` の点 `Q` に対し、**仮説なしで**

* `z₀`（`Φ(z₀) = Q`）と周期対 `P′`（格子は `Λ′ = Λ + ℤz₀`）
* 整数 `A, B, C, D` で `ω₁ = Aω₁′ + Bω₂′`・`ω₂ = Cω₁′ + Dω₂′`・`|AD − BC| = l`
* 点の集合 `S`（`|S| = l − 1`）で

      latticeCurve Λ′ = veluQuotientFull (latticeCurve Λ) S

★★★★★★★☆**変数変換は要らない（`α = 1`）。これが
`htFalt_isogeny_le_of_velu`（第 678）に渡す形である。** -/
theorem exists_veluQuotientFull_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ) (S : Finset (ℂ × ℂ)),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      S.card = l - 1 ∧
      latticeCurve P' = veluQuotientFull (latticeCurve P) S := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep, hTdef⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  have hvelu : ∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z
  have hnot : ∀ w ∈ T.erase 0, w ∉ P.lattice :=
    fun w hw => rep_notMem_lattice P P' T h0T hT hrep hw
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0,
      (latticePointX P w, latticePointY P w)
        = (latticePointX P v, latticePointY P v) → w = v := by
    intro w hw v hv he
    refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw)
      (Finset.mem_of_mem_erase hv) ?_
    exact sub_mem_lattice_of_point_eq P (hnot w hw) (hnot v hv)
      (congrArg Prod.fst he) (congrArg Prod.snd he)
  refine ⟨z₀, P', A, B, C, D,
    (T.erase 0).image (fun w => (latticePointX P w, latticePointY P w)),
    hz₀, h1, h2, hdet, hP', ?_, ?_⟩
  · rw [Finset.card_image_of_injOn (fun w hw v hv he =>
      hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he),
      Finset.card_erase_of_mem h0T, hcard]
  · exact latticeCurve_eq_veluQuotientFull P P' T h0T hT hrep
      (g₂_isogeny P P' T h0T hT hrep hvelu) (g₃_isogeny P P' T h0T hT hrep hvelu)

def latticeCurve_eq_veluQuotientFull.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(latticeCurve Λ′ は H∖{O} 全体で書いた Vélu の商そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_veluQuotientFull_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側・代数語彙——位数 l の点から E/H まで。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
