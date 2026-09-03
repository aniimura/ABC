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
import ABC3.Found.GenEll.Uniformization.FilledPole

/-!
# 一様化 —— `y` 座標の加法公式・mathlib の `Point` への橋・倍加公式

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★`y` 座標の加法公式へ -/

/-- ★★★★★★★★★★★★★★**`y` 座標の加法公式の代数の核**——2 つの ODE だけから。

    `℘″(z)·D − N·℘′(z) = −N²/2 + (4℘(z) + 2℘(w))·D²`

（`N = ℘′(z) − ℘′(w)`、`D = ℘(z) − ℘(w)`、`℘″ = 6℘² − g₂/2`）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`℘′² = 4℘³ − g₂℘ − g₃` を `z` と `w` の両方に当てると**両辺とも**

    `2℘z³ − 6℘z²℘w + (g₂/2)℘z + (g₂/2)℘w + g₃ + ℘′z·℘′w`

になる。☆これが `R′ = −R²/2 + 4℘(z) + 2℘(w)` の中身であり、
加法定理を微分して得られる `℘′(z+w) = R·R′/2 − ℘′(z)` を
群法則の形 `℘′(z+w) = −R·(℘(z+w) − ℘(z)) − ℘′(z)` に直す鍵である。 -/
theorem y_addition_algebraic (P : PeriodPair) (z w : ℂ)
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) :
    (6 * P.weierstrassP z ^ 2 - P.g₂ / 2) * (P.weierstrassP z - P.weierstrassP w)
        - (P.derivWeierstrassP z - P.derivWeierstrassP w) * P.derivWeierstrassP z
      = -(P.derivWeierstrassP z - P.derivWeierstrassP w) ^ 2 / 2
        + (4 * P.weierstrassP z + 2 * P.weierstrassP w)
          * (P.weierstrassP z - P.weierstrassP w) ^ 2 := by
  have h1 := P.derivWeierstrassP_sq z hz
  have h2 := P.derivWeierstrassP_sq w hw
  linear_combination (-1/2 : ℂ) * h1 + (1/2 : ℂ) * h2

def y_addition_algebraic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(y 座標の加法公式の代数の核——2 つの ODE だけから。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`y` 座標の加法公式**

    `℘′(z+w) = −R·(℘(z+w) − ℘(z)) − ℘′(z)`,  `R = (℘′(z) − ℘′(w))/(℘(z) − ℘(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これで一様化 `(℘, ℘′/2) : ℂ/Λ → E(ℂ)` は群準同型である**
——`latticeCurve` の群法則（`a₁ = a₂ = a₃ = 0`）は
`x₃ = λ² − x₁ − x₂`・`y₃ = −λ(x₃ − x₁) − y₁`（`λ = R/2`）だから。

★機構: 加法定理（第 641）を `z` で微分すると `℘′(z+w) = R·R′/2 − ℘′(z)`。
第 644 の代数の核から `R′ = −R²/2 + 4℘(z) + 2℘(w) = −2(℘(z+w) − ℘(z))`。 -/
theorem derivWeierstrassP_addition (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    P.derivWeierstrassP (z + w)
      = -((P.derivWeierstrassP z - P.derivWeierstrassP w)
            / (P.weierstrassP z - P.weierstrassP w))
          * (P.weierstrassP (z + w) - P.weierstrassP z)
        - P.derivWeierstrassP z := by
  set D : ℂ := P.weierstrassP z - P.weierstrassP w with hD
  set N : ℂ := P.derivWeierstrassP z - P.derivWeierstrassP w with hN
  -- 左辺の微分
  have hL : HasDerivAt (fun s : ℂ => P.weierstrassP (s + w))
      (P.derivWeierstrassP (z + w)) z :=
    HasDerivAt.comp_add_const z w (hasDerivAt_weierstrassP P hzw)
  -- 右辺の微分
  have hnum : HasDerivAt (fun s : ℂ => P.derivWeierstrassP s - P.derivWeierstrassP w)
      (6 * P.weierstrassP z ^ 2 - P.g₂ / 2) z := by
    have h := (hasDerivAt_derivWeierstrassP P hz).sub_const (P.derivWeierstrassP w)
    rwa [deriv_derivWeierstrassP P hz] at h
  have hden : HasDerivAt (fun s : ℂ => P.weierstrassP s - P.weierstrassP w)
      (P.derivWeierstrassP z) z :=
    (hasDerivAt_weierstrassP P hz).sub_const _
  have hR : HasDerivAt (fun s : ℂ => (P.derivWeierstrassP s - P.derivWeierstrassP w)
        / (P.weierstrassP s - P.weierstrassP w))
      (((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) * D - N * P.derivWeierstrassP z) / D ^ 2) z :=
    hnum.div hden hne
  have hRHS : HasDerivAt (fun s : ℂ => ((P.derivWeierstrassP s - P.derivWeierstrassP w)
        / (P.weierstrassP s - P.weierstrassP w)) ^ 2 / 4
      - P.weierstrassP s - P.weierstrassP w)
      (2 * (N / D) * (((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) * D
          - N * P.derivWeierstrassP z) / D ^ 2) / 4 - P.derivWeierstrassP z) z := by
    have h1 : HasDerivAt (fun s : ℂ => ((P.derivWeierstrassP s - P.derivWeierstrassP w)
          / (P.weierstrassP s - P.weierstrassP w)) ^ 2 / 4)
        (2 * (N / D) * (((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) * D
          - N * P.derivWeierstrassP z) / D ^ 2) / 4) z := by
      have h2 := (hR.pow 2).div_const 4
      simpa using h2
    simpa using (h1.sub (hasDerivAt_weierstrassP P hz)).sub_const (P.weierstrassP w)
  -- 加法定理は開集合で成り立つ
  have hEq : (fun s : ℂ => P.weierstrassP (s + w))
      =ᶠ[nhds z] fun s : ℂ => ((P.derivWeierstrassP s - P.derivWeierstrassP w)
        / (P.weierstrassP s - P.weierstrassP w)) ^ 2 / 4
      - P.weierstrassP s - P.weierstrassP w := by
    filter_upwards [(isOpen_goodSet P w).mem_nhds ⟨hz, hzw, hne⟩] with s hs
    exact weierstrassP_addition P w hw h2w hs.1 hs.2.1 hs.2.2
  have hderiv : P.derivWeierstrassP (z + w)
      = 2 * (N / D) * (((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) * D
          - N * P.derivWeierstrassP z) / D ^ 2) / 4 - P.derivWeierstrassP z := by
    rw [← hL.deriv, hEq.deriv_eq, hRHS.deriv]
  -- 代数の核（第 644）を入れる
  have halg := y_addition_algebraic P z w hz hw
  rw [← hD, ← hN] at halg
  have hadd := weierstrassP_addition P w hw h2w hz hzw hne
  rw [← hD, ← hN] at hadd
  rw [hderiv, hadd, halg]
  field_simp
  ring

def derivWeierstrassP_addition.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(y 座標の加法公式——一様化は群準同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★mathlib の `Point` へのパッケージ -/

/-- ★★★★★★`(℘(z), ℘′(z)/2)` は非特異——`Δ ≠ 0` だから。 -/
theorem nonsingular_latticePoint (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ)
    (hz : z ∉ P.lattice) :
    (latticeCurve P).toAffine.Nonsingular (latticePointX P z) (latticePointY P z) := by
  refine ((latticeCurve P).toAffine.equation_iff_nonsingular_of_Δ_ne_zero ?_).1
    (latticeCurve_equation P z hz)
  rw [latticeCurve_Δ]
  exact hΔ

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**一様化は群準同型**
——mathlib の `Affine.Point` の加法と一致する。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`latticeCurve P` は `a₁ = a₂ = a₃ = 0` なので

* `slope = (y₁ − y₂)/(x₁ − x₂) = R/2`
* `addX = ℓ² − x₁ − x₂ = R²/4 − ℘z − ℘w = ℘(z+w)`（第 641 の加法定理）
* `addY = −ℓ(addX − x₁) − y₁ = −(R/2)(℘(z+w) − ℘z) − ℘′z/2 = ℘′(z+w)/2`
  （第 645 の `y` 座標の公式）

★★★★☆**これで `Φ : ℂ → E(ℂ)` は群準同型であり、
第 603-604・624 の全単射性と合わせて群同型である。** -/
theorem latticePoint_add (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (w z : ℂ)
    (hw : w ∉ P.lattice) (h2w : 2 * w ∉ P.lattice) (hz : z ∉ P.lattice)
    (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    (WeierstrassCurve.Affine.Point.some (latticePointX P z) (latticePointY P z)
        (nonsingular_latticePoint P hΔ z hz)
      + WeierstrassCurve.Affine.Point.some (latticePointX P w) (latticePointY P w)
        (nonsingular_latticePoint P hΔ w hw) : (latticeCurve P).toAffine.Point)
      = WeierstrassCurve.Affine.Point.some (latticePointX P (z + w))
          (latticePointY P (z + w)) (nonsingular_latticePoint P hΔ (z + w) hzw) := by
  have hxne : latticePointX P z ≠ latticePointX P w := by
    intro hc
    exact hne (by simp only [latticePointX] at hc; rw [hc]; ring)
  have hxy : ¬(latticePointX P z = latticePointX P w
      ∧ latticePointY P z
        = (latticeCurve P).toAffine.negY (latticePointX P w) (latticePointY P w)) :=
    fun h => hxne h.1
  have hslope : (latticeCurve P).toAffine.slope (latticePointX P z) (latticePointX P w)
      (latticePointY P z) (latticePointY P w)
      = (P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w) / 2 := by
    rw [WeierstrassCurve.Affine.slope, if_neg hxne]
    simp only [latticePointX, latticePointY]
    field_simp
  have hadd := weierstrassP_addition P w hw h2w hz hzw hne
  have hadd' := derivWeierstrassP_addition P w hw h2w hz hzw hne
  have hX : (latticeCurve P).toAffine.addX (latticePointX P z) (latticePointX P w)
      ((latticeCurve P).toAffine.slope (latticePointX P z) (latticePointX P w)
        (latticePointY P z) (latticePointY P w))
      = latticePointX P (z + w) := by
    rw [hslope, WeierstrassCurve.Affine.addX]
    simp only [latticeCurve, latticePointX]
    rw [hadd]
    ring
  have hY : (latticeCurve P).toAffine.addY (latticePointX P z) (latticePointX P w)
      (latticePointY P z)
      ((latticeCurve P).toAffine.slope (latticePointX P z) (latticePointX P w)
        (latticePointY P z) (latticePointY P w))
      = latticePointY P (z + w) := by
    rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY, hX, hslope]
    simp only [WeierstrassCurve.Affine.negY, latticeCurve, latticePointX, latticePointY]
    linear_combination (-1/2 : ℂ) * hadd'
  rw [WeierstrassCurve.Affine.Point.add_some hxy]
  congr 1

def nonsingular_latticePoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5((℘(z), ℘′(z)/2) は非特異——Δ ≠ 0 から。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticePoint_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化は群準同型——mathlib の Point の加法と一致する。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★倍加公式 -/

open Filter Topology in
/-- ★★★★★★`℘(y) ≠ ℘(z)` は `z` の除いた近傍で成り立つ（`2z ∉ Λ`）。

★第 629 の「`℘ − ℘(z)` は `z` で 1 位の零点」から、零点が孤立するので従う。 -/
theorem eventually_weierstrassP_ne_self (P : PeriodPair) {z : ℂ} (hz : z ∉ P.lattice)
    (h2z : 2 * z ∉ P.lattice) :
    ∀ᶠ y in 𝓝[≠] z, P.weierstrassP y - P.weierstrassP z ≠ 0 := by
  have hana : AnalyticAt ℂ (fun y : ℂ => P.weierstrassP y - P.weierstrassP z) z :=
    (P.analyticOnNhd_weierstrassP z hz).sub analyticAt_const
  have ho : analyticOrderAt (fun y : ℂ => P.weierstrassP y - P.weierstrassP z) z ≠ ⊤ := by
    rw [analyticOrderAt_weierstrassP_sub_self P z hz h2z]
    decide
  rcases hana.eventually_eq_zero_or_eventually_ne_zero with h1 | h2
  · exact absurd (analyticOrderAt_eq_top.2 h1) ho
  · exact h2

open Filter Topology in
/-- ★★★★★★★★★★`w → z` で `R(z,w) → ℘″(z)/℘′(z)`。 -/
theorem tendsto_slopeRatio (P : PeriodPair) {z : ℂ} (hz : z ∉ P.lattice)
    (h2z : 2 * z ∉ P.lattice) :
    Tendsto (fun w : ℂ => (P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w)) (𝓝[≠] z)
      (nhds ((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) / P.derivWeierstrassP z)) := by
  have hpne : P.derivWeierstrassP z ≠ 0 := fun hc =>
    h2z ((derivWeierstrassP_eq_zero_iff P z hz).1 hc)
  have hnum : Tendsto (slope P.derivWeierstrassP z) (𝓝[≠] z)
      (nhds (6 * P.weierstrassP z ^ 2 - P.g₂ / 2)) := by
    have h := hasDerivAt_derivWeierstrassP P hz
    rw [deriv_derivWeierstrassP P hz] at h
    exact hasDerivAt_iff_tendsto_slope.1 h
  have hden : Tendsto (slope P.weierstrassP z) (𝓝[≠] z) (nhds (P.derivWeierstrassP z)) :=
    hasDerivAt_iff_tendsto_slope.1 (hasDerivAt_weierstrassP P hz)
  refine (hnum.div hden hpne).congr' ?_
  filter_upwards [self_mem_nhdsWithin, eventually_weierstrassP_ne_self P hz h2z]
    with w hw hne2
  have hwz : w - z ≠ 0 := sub_ne_zero.2 (by simpa using hw)
  have hne3 : P.weierstrassP z - P.weierstrassP w ≠ 0 := fun hc =>
    hne2 (by linear_combination -hc)
  simp only [Pi.div_apply, slope_def_field]
  rw [div_div_div_eq, div_eq_div_iff (mul_ne_zero hwz hne2) hne3]
  ring

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★**倍加公式**

    `℘(z+z) = (℘″(z)/℘′(z))²/4 − 2℘(z)`,  `℘″ = 6℘² − g₂/2`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★加法定理（第 641）で `w → z` の極限を取る。
☆`R(z,w) → ℘″(z)/℘′(z)`・`℘(z+w) → ℘(z+z)`・`℘(w) → ℘(z)`。

★★★これで mathlib の `slope`（`x₁ = x₂` の枝、`(3℘z² − g₂/4)/℘′z = ℘″/(2℘′)`）と
一致し、`Φ : ℂ → E(ℂ)` の群準同型が**倍加の場合**でも書ける。 -/
theorem weierstrassP_duplication (P : PeriodPair) {z : ℂ} (hz : z ∉ P.lattice)
    (h2z : 2 * z ∉ P.lattice) :
    P.weierstrassP (z + z)
      = ((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) / P.derivWeierstrassP z) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP z := by
  have hzz : z + z ∉ P.lattice := by
    intro hc; exact h2z (by simpa [two_mul] using hc)
  have hL : Tendsto (fun w : ℂ => P.weierstrassP (z + w)) (𝓝[≠] z)
      (nhds (P.weierstrassP (z + z))) := by
    have hf : Tendsto (fun w : ℂ => z + w) (𝓝[≠] z) (nhds (z + z)) := by
      have ht : Tendsto (fun w : ℂ => z + w) (nhds z) (nhds (z + z)) :=
        (continuous_const.add continuous_id).tendsto _
      exact ht.mono_left nhdsWithin_le_nhds
    exact ((P.analyticOnNhd_weierstrassP (z + z) hzz).continuousAt).tendsto.comp hf
  have hp : Tendsto P.weierstrassP (𝓝[≠] z) (nhds (P.weierstrassP z)) :=
    ((P.analyticOnNhd_weierstrassP z hz).continuousAt).continuousWithinAt.tendsto
  have hR : Tendsto (fun w : ℂ => ((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
      - P.weierstrassP z - P.weierstrassP w) (𝓝[≠] z)
      (nhds (((6 * P.weierstrassP z ^ 2 - P.g₂ / 2) / P.derivWeierstrassP z) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP z)) :=
    ((((tendsto_slopeRatio P hz h2z).pow 2).div_const 4).sub_const
      (P.weierstrassP z)).sub hp
  have hEq : (fun w : ℂ => P.weierstrassP (z + w)) =ᶠ[𝓝[≠] z]
      fun w : ℂ => ((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
      - P.weierstrassP z - P.weierstrassP w := by
    have hL1 : ∀ᶠ w in 𝓝[≠] z, w ∉ P.lattice :=
      mem_nhdsWithin_of_mem_nhds ((P.isClosed_lattice.isOpen_compl).mem_nhds hz)
    have hL2 : ∀ᶠ w in 𝓝[≠] z, z + w ∉ P.lattice := by
      have hopen : IsOpen {w : ℂ | z + w ∉ P.lattice} := by
        have he : {w : ℂ | z + w ∉ P.lattice}
            = (fun w : ℂ => z + w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
        rw [he]
        exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
      exact mem_nhdsWithin_of_mem_nhds (hopen.mem_nhds hzz)
    have hL3 : ∀ᶠ w in 𝓝[≠] z, P.weierstrassP z - P.weierstrassP w ≠ 0 := by
      filter_upwards [eventually_weierstrassP_ne_self P hz h2z] with y hy
      intro hc
      exact hy (by linear_combination -hc)
    have hL4 : ∀ᶠ w in 𝓝[≠] z, 2 * w ∉ P.lattice := by
      have hopen : IsOpen {w : ℂ | 2 * w ∉ P.lattice} := by
        have he : {w : ℂ | 2 * w ∉ P.lattice}
            = (fun w : ℂ => 2 * w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
        rw [he]
        exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
      exact mem_nhdsWithin_of_mem_nhds (hopen.mem_nhds h2z)
    filter_upwards [hL1, hL2, hL3, hL4] with w hw1 hw2 hw3 hw4
    exact weierstrassP_addition P w hw1 hw4 hz hw2 hw3
  exact tendsto_nhds_unique (hL.congr' hEq) hR

def weierstrassP_duplication.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(倍加公式——加法定理の w → z 極限。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
