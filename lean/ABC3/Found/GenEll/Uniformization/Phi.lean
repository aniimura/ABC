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
import ABC3.Found.GenEll.Uniformization.AdditionFormula

/-!
# 一様化 —— 一様化写像 `Φ : ℂ → E(ℂ)`

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★一様化写像 `Φ : ℂ → E(ℂ)` -/

/-- ★★★★★`Point.some` の合同——`Nonsingular` は `Prop` だから。 -/
theorem point_some_congr {P : PeriodPair} {x₁ y₁ x₂ y₂ : ℂ}
    {h₁ : (latticeCurve P).toAffine.Nonsingular x₁ y₁}
    {h₂ : (latticeCurve P).toAffine.Nonsingular x₂ y₂}
    (hx : x₁ = x₂) (hy : y₁ = y₂) :
    (WeierstrassCurve.Affine.Point.some x₁ y₁ h₁ : (latticeCurve P).toAffine.Point)
      = WeierstrassCurve.Affine.Point.some x₂ y₂ h₂ := by
  subst hx; subst hy; rfl

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**一様化写像**
`Φ(z) = (℘(z), ℘′(z)/2)`（格子点では `O`）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★☆これが `ℂ/Λ ≅ E(ℂ)` の実体である——全射（第 603-604）・単射（第 624）・
群準同型（第 648・650）。 -/
noncomputable def uniformMap (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    (latticeCurve P).toAffine.Point :=
  if h : z ∈ P.lattice then 0
  else WeierstrassCurve.Affine.Point.some (latticePointX P z) (latticePointY P z)
    (nonsingular_latticePoint P hΔ z h)

open scoped Classical in
@[simp] theorem uniformMap_of_mem (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z : ℂ}
    (h : z ∈ P.lattice) : uniformMap P hΔ z = 0 := by
  simp [uniformMap, h]

open scoped Classical in
theorem uniformMap_of_notMem (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z : ℂ}
    (h : z ∉ P.lattice) :
    uniformMap P hΔ z = WeierstrassCurve.Affine.Point.some (latticePointX P z)
      (latticePointY P z) (nonsingular_latticePoint P hΔ z h) := by
  simp [uniformMap, h]

/-- ★★★★★★★★**`Φ` は `Λ`-周期的**。 -/
theorem uniformMap_periodic (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) {l : ℂ}
    (hl : l ∈ P.lattice) : uniformMap P hΔ (z + l) = uniformMap P hΔ z := by
  by_cases hz : z ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ (P.lattice.add_mem hz hl), uniformMap_of_mem P hΔ hz]
  · have hzl : z + l ∉ P.lattice := fun hc => hz (by simpa using P.lattice.sub_mem hc hl)
    rw [uniformMap_of_notMem P hΔ hzl, uniformMap_of_notMem P hΔ hz]
    refine point_some_congr ?_ ?_
    · simp only [latticePointX]
      exact P.weierstrassP_add_coe z ⟨l, hl⟩
    · simp only [latticePointY]
      rw [P.derivWeierstrassP_add_coe z ⟨l, hl⟩]

def uniformMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化写像 Φ(z) = (℘(z), ℘′(z)/2))",
    sectionId := "genell-lemma-3-5" }

def uniformMap_periodic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は Λ-周期的。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★**比 `R` は `z` と `w` の入れ替えで不変**——分子・分母がともに符号を変えるから。 -/
theorem slopeRatio_symm (P : PeriodPair) (z w : ℂ) :
    (P.derivWeierstrassP w - P.derivWeierstrassP z)
        / (P.weierstrassP w - P.weierstrassP z)
      = (P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w) := by
  rw [show P.derivWeierstrassP w - P.derivWeierstrassP z
      = -(P.derivWeierstrassP z - P.derivWeierstrassP w) by ring,
    show P.weierstrassP w - P.weierstrassP z
      = -(P.weierstrassP z - P.weierstrassP w) by ring, neg_div_neg_eq]

/-- ★★★★★★★★★★★★★★★★**加法定理の対称版**——`2z ∉ Λ` から。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 641 は `2w ∉ Λ` を仮定するが、`℘(z+w) = ℘(w+z)` と `R` の対称性
（`slopeRatio_symm`）で `2z ∉ Λ` からも同じ結論が出る。
☆これで「`z` か `w` のどちらかが 2-捩れでない」場合が尽きる（第 653 の記録）。 -/
theorem weierstrassP_addition' (P : PeriodPair) (z : ℂ) (hz : z ∉ P.lattice)
    (h2z : 2 * z ∉ P.lattice) {w : ℂ}
    (hw : w ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    P.weierstrassP (z + w)
      = ((P.derivWeierstrassP z - P.derivWeierstrassP w)
          / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP w := by
  have hne' : P.weierstrassP w - P.weierstrassP z ≠ 0 := fun hc =>
    hne (by linear_combination -hc)
  have hwz : w + z ∉ P.lattice := by rw [add_comm]; exact hzw
  have h := weierstrassP_addition P z hz h2z hw hwz hne'
  rw [add_comm w z] at h
  rw [h, slopeRatio_symm]
  ring

def weierstrassP_addition'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(加法定理の対称版——2z ∉ Λ から。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★**`4x³ − g₂x − g₃` の 3 つの相異なる根の和は `0`**（Vieta）。

★2 つの方程式を引くと `(e₁−e₂)(4(e₁²+e₁e₂+e₂²) − g₂) = 0`、
さらに 2 つの関係を引くと `(e₂−e₃)(e₁+e₂+e₃) = 0` になる。 -/
theorem sum_roots_eq_zero (g₂ g₃ e₁ e₂ e₃ : ℂ)
    (h1 : 4 * e₁ ^ 3 - g₂ * e₁ - g₃ = 0) (h2 : 4 * e₂ ^ 3 - g₂ * e₂ - g₃ = 0)
    (h3 : 4 * e₃ ^ 3 - g₂ * e₃ - g₃ = 0)
    (h12 : e₁ ≠ e₂) (h13 : e₁ ≠ e₃) (h23 : e₂ ≠ e₃) :
    e₁ + e₂ + e₃ = 0 := by
  have d12 : e₁ - e₂ ≠ 0 := sub_ne_zero.2 h12
  have d13 : e₁ - e₃ ≠ 0 := sub_ne_zero.2 h13
  have d23 : e₂ - e₃ ≠ 0 := sub_ne_zero.2 h23
  have k1 : 4 * (e₁ ^ 2 + e₁ * e₂ + e₂ ^ 2) - g₂ = 0 := by
    refine mul_left_cancel₀ d12 ?_
    rw [mul_zero]
    linear_combination h1 - h2
  have k2 : 4 * (e₁ ^ 2 + e₁ * e₃ + e₃ ^ 2) - g₂ = 0 := by
    refine mul_left_cancel₀ d13 ?_
    rw [mul_zero]
    linear_combination h1 - h3
  refine mul_left_cancel₀ d23 ?_
  rw [mul_zero]
  linear_combination (k1 - k2) / 4

/-- ★★★★★★★★★★★★★★★★★★**両方 2-捩れの場合の加法定理**

    `2z, 2w ∈ Λ`・`℘z ≠ ℘w` なら `℘(z+w) = −℘z − ℘w`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★このとき `℘′z = ℘′w = 0`（第 605）なので群法則の `slope` は `0`、
すなわち `x₃ = −x₁ − x₂` である。
☆`℘z`・`℘w`・`℘(z+w)` は `4x³ − g₂x − gₙ` の**相異なる 3 根**なので
（第 605 の `cubic_eq_zero_of_two_mem` と第 624 の単射性）、Vieta で和が `0`。

★★これで第 653 に記録した 3 つの場合がすべて済んだ。 -/
theorem weierstrassP_addition_two_torsion (P : PeriodPair) {z w : ℂ}
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (h2z : 2 * z ∈ P.lattice) (h2w : 2 * w ∈ P.lattice)
    (hne : P.weierstrassP z ≠ P.weierstrassP w) :
    P.weierstrassP (z + w) = -P.weierstrassP z - P.weierstrassP w := by
  have h2zw : 2 * (z + w) ∈ P.lattice := by
    have : 2 * (z + w) = 2 * z + 2 * w := by ring
    rw [this]
    exact P.lattice.add_mem h2z h2w
  have e1 := cubic_eq_zero_of_two_mem P z hz h2z
  have e2 := cubic_eq_zero_of_two_mem P w hw h2w
  have e3 := cubic_eq_zero_of_two_mem P (z + w) hzw h2zw
  simp only [latticePointX] at e1 e2 e3
  -- 3 つの値は相異なる
  have hd13 : P.weierstrassP z ≠ P.weierstrassP (z + w) := by
    intro hc
    refine hw ?_
    have hmem := mem_lattice_of_shift_eq P w hz (by simpa using hzw) (by simpa using hc.symm)
      (by
        rw [show z + w = z + w from rfl,
          derivWeierstrassP_eq_zero_of_two_mem P (z + w) h2zw,
          derivWeierstrassP_eq_zero_of_two_mem P z h2z])
    exact hmem
  have hd23 : P.weierstrassP w ≠ P.weierstrassP (z + w) := by
    intro hc
    refine hz ?_
    have hwzc : w + z ∉ P.lattice := by rw [add_comm]; exact hzw
    have hmem := mem_lattice_of_shift_eq P z hw (by simpa using hwzc)
      (by rw [show w + z = z + w by ring]; exact hc.symm)
      (by
        rw [show w + z = z + w by ring,
          derivWeierstrassP_eq_zero_of_two_mem P (z + w) h2zw,
          derivWeierstrassP_eq_zero_of_two_mem P w h2w])
    exact hmem
  have := sum_roots_eq_zero P.g₂ P.g₃ (P.weierstrassP z) (P.weierstrassP w)
    (P.weierstrassP (z + w)) (by linear_combination e1) (by linear_combination e2)
    (by linear_combination e3) hne hd13 hd23
  linear_combination this

def sum_roots_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(4x³ − g₂x − g₃ の 3 根の和は 0——Vieta。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_addition_two_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(両方 2-捩れの場合の加法定理。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**加法定理（一般形）**
——2-捩れの仮定なし。

    `℘(z+w) = (1/4)·((℘′(z) − ℘′(w))/(℘(z) − ℘(w)))² − ℘(z) − ℘(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★3 つの場合に分ける（第 653 の記録）:

* `2w ∉ Λ` —— 第 641
* `2w ∈ Λ` かつ `2z ∉ Λ` —— 第 654（`R` の対称性）
* 両方 2-捩れ —— 第 655（`R = 0` になり `℘(z+w) = −℘z − ℘w`）

★★★☆**これで `℘` の加法定理が無条件（`z, w, z+w ∉ Λ` と `℘z ≠ ℘w` のみ）になった**。 -/
theorem weierstrassP_addition_general (P : PeriodPair) {z w : ℂ}
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    P.weierstrassP (z + w)
      = ((P.derivWeierstrassP z - P.derivWeierstrassP w)
          / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP w := by
  by_cases h2w : 2 * w ∈ P.lattice
  · by_cases h2z : 2 * z ∈ P.lattice
    · -- 両方 2-捩れ
      have hpz : P.derivWeierstrassP z = 0 :=
        derivWeierstrassP_eq_zero_of_two_mem P z h2z
      have hpw : P.derivWeierstrassP w = 0 :=
        derivWeierstrassP_eq_zero_of_two_mem P w h2w
      have hxne : P.weierstrassP z ≠ P.weierstrassP w := fun hc =>
        hne (by linear_combination hc)
      rw [weierstrassP_addition_two_torsion P hz hw hzw h2z h2w hxne, hpz, hpw]
      simp
    · -- 2z ∉ Λ
      exact weierstrassP_addition' P z hz h2z hw hzw hne
  · exact weierstrassP_addition P w hw h2w hz hzw hne

def weierstrassP_addition_general.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の加法定理・一般形——2-捩れの仮定なし。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**`y` 座標の加法公式（一般形）**
——2-捩れの仮定なし。

    `℘′(z+w) = −R·(℘(z+w) − ℘(z)) − ℘′(z)`,  `R = (℘′(z) − ℘′(w))/(℘(z) − ℘(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★3 つの場合:

* `2w ∉ Λ` —— 第 645
* `2w ∈ Λ` かつ `2z ∉ Λ` —— 役を入れ替えた第 645 と `R(℘z − ℘w) = ℘′z`（`℘′w = 0`）
* 両方 2-捩れ —— `℘′(z+w) = ℘′z = ℘′w = 0` で両辺 `0`

★★★☆**これで一様化の群準同型が 2-捩れも込めて書ける**。 -/
theorem derivWeierstrassP_addition_general (P : PeriodPair) {z w : ℂ}
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    P.derivWeierstrassP (z + w)
      = -((P.derivWeierstrassP z - P.derivWeierstrassP w)
            / (P.weierstrassP z - P.weierstrassP w))
          * (P.weierstrassP (z + w) - P.weierstrassP z)
        - P.derivWeierstrassP z := by
  by_cases h2w : 2 * w ∈ P.lattice
  · by_cases h2z : 2 * z ∈ P.lattice
    · have h2zw : 2 * (z + w) ∈ P.lattice := by
        have he : 2 * (z + w) = 2 * z + 2 * w := by ring
        rw [he]
        exact P.lattice.add_mem h2z h2w
      rw [derivWeierstrassP_eq_zero_of_two_mem P (z + w) h2zw,
        derivWeierstrassP_eq_zero_of_two_mem P z h2z,
        derivWeierstrassP_eq_zero_of_two_mem P w h2w]
      simp
    · have hne' : P.weierstrassP w - P.weierstrassP z ≠ 0 := fun hc =>
        hne (by linear_combination -hc)
      have hwz : w + z ∉ P.lattice := by rw [add_comm]; exact hzw
      have h := derivWeierstrassP_addition P z hz h2z hw hwz hne'
      rw [add_comm w z] at h
      rw [h, slopeRatio_symm, derivWeierstrassP_eq_zero_of_two_mem P w h2w]
      field_simp
      ring
  · exact derivWeierstrassP_addition P w hw h2w hz hzw hne

def derivWeierstrassP_addition_general.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(y 座標の加法公式・一般形——2-捩れの仮定なし。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★**`Point` の加法との一致（一般形）**
——2-捩れの仮定なし。 -/
theorem latticePoint_add_general (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z w : ℂ}
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice) (hzw : z + w ∉ P.lattice)
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
  have hadd := weierstrassP_addition_general P hz hw hzw hne
  have hadd' := derivWeierstrassP_addition_general P hz hw hzw hne
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

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`Φ` は群準同型**

    `Φ(z + w) = Φ(z) + Φ(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★5 つの場合で尽きる（第 651 の記録）:

1. `z ∈ Λ` —— `Φz = 0`、周期性
2. `w ∈ Λ` —— 対称
3. `z + w ∈ Λ` —— `w ≡ −z` なので `x₁ = x₂`・`y₁ = negY x₂ y₂`、`Point.add_of_Y_eq`
4. `℘z ≠ ℘w` —— 第 658 の `latticePoint_add_general`
5. `℘z = ℘w`（かつ `z+w ∉ Λ`）—— 単射性（第 624）で `z ≡ w`、すなわち**倍加**

★★★★★☆**これで一様化 `ℂ/Λ → E(ℂ)` は全単射（第 603-604・624）かつ群準同型
——群同型である。** -/
theorem uniformMap_add_of_ne (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z w : ℂ}
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0 ∨ z ∈ P.lattice ∨ w ∈ P.lattice
      ∨ z + w ∈ P.lattice) :
    uniformMap P hΔ (z + w) = uniformMap P hΔ z + uniformMap P hΔ w := by
  by_cases hz : z ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hz, zero_add, add_comm z w]
    exact uniformMap_periodic P hΔ w hz
  by_cases hw : w ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hw, add_zero]
    exact uniformMap_periodic P hΔ z hw
  by_cases hzw : z + w ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hzw, uniformMap_of_notMem P hΔ hz,
      uniformMap_of_notMem P hΔ hw]
    have hpw : P.weierstrassP w = P.weierstrassP z := by
      have h1 : P.weierstrassP ((-z) + (z + w)) = P.weierstrassP (-z) :=
        P.weierstrassP_add_coe (-z) ⟨z + w, hzw⟩
      rw [show (-z) + (z + w) = w by ring, P.weierstrassP_neg] at h1
      exact h1
    have hpdw : P.derivWeierstrassP w = -P.derivWeierstrassP z := by
      have h1 : P.derivWeierstrassP ((-z) + (z + w)) = P.derivWeierstrassP (-z) :=
        P.derivWeierstrassP_add_coe (-z) ⟨z + w, hzw⟩
      rw [show (-z) + (z + w) = w by ring, P.derivWeierstrassP_neg] at h1
      exact h1
    refine (WeierstrassCurve.Affine.Point.add_of_Y_eq ?_ ?_).symm
    · simp only [latticePointX]
      exact hpw.symm
    · simp only [latticePointY, WeierstrassCurve.Affine.negY, latticeCurve, latticePointX]
      rw [hpdw]
      ring
  · rcases hne with hne | hne | hne | hne
    · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_notMem P hΔ hw,
        uniformMap_of_notMem P hΔ hzw]
      exact (latticePoint_add_general P hΔ hz hw hzw hne).symm
    · exact absurd hne hz
    · exact absurd hne hw
    · exact absurd hne hzw

def latticePoint_add_general.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Point の加法との一致・一般形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def uniformMap_add_of_ne.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は群準同型——倍加以外の場合)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★**倍加の `y` 座標の公式**

    `℘′(w+w) = −(℘″(w)/℘′(w))·(℘(w+w) − ℘(w)) − ℘′(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 657（`y` 座標の一般形）で `z → w` の極限を取る。
☆`R(z,w) = R(w,z) → ℘″(w)/℘′(w)`（第 649・第 654）、
`℘′(z+w) → ℘′(w+w)`・`℘(z+w) → ℘(w+w)`・`℘z → ℘w`・`℘′z → ℘′w`。

★★★これで mathlib の `addY`（`x₁ = x₂` の枝、`slope = ℘″/(2℘′)`）と一致し、
`Φ` の群準同型が**倍加の場合**でも書ける。 -/
theorem derivWeierstrassP_duplication (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    P.derivWeierstrassP (w + w)
      = -((6 * P.weierstrassP w ^ 2 - P.g₂ / 2) / P.derivWeierstrassP w)
          * (P.weierstrassP (w + w) - P.weierstrassP w)
        - P.derivWeierstrassP w := by
  have hww : w + w ∉ P.lattice := by
    intro hc; exact h2w (by simpa [two_mul] using hc)
  have hshift : Tendsto (fun z : ℂ => z + w) (𝓝[≠] w) (nhds (w + w)) := by
    have ht : Tendsto (fun z : ℂ => z + w) (nhds w) (nhds (w + w)) :=
      (continuous_id.add continuous_const).tendsto _
    exact ht.mono_left nhdsWithin_le_nhds
  have hLd : Tendsto (fun z : ℂ => P.derivWeierstrassP (z + w)) (𝓝[≠] w)
      (nhds (P.derivWeierstrassP (w + w))) :=
    ((P.analyticOnNhd_derivWeierstrassP (w + w) hww).continuousAt).tendsto.comp hshift
  have hLp : Tendsto (fun z : ℂ => P.weierstrassP (z + w)) (𝓝[≠] w)
      (nhds (P.weierstrassP (w + w))) :=
    ((P.analyticOnNhd_weierstrassP (w + w) hww).continuousAt).tendsto.comp hshift
  have hp : Tendsto P.weierstrassP (𝓝[≠] w) (nhds (P.weierstrassP w)) :=
    ((P.analyticOnNhd_weierstrassP w hw).continuousAt).continuousWithinAt.tendsto
  have hpd : Tendsto P.derivWeierstrassP (𝓝[≠] w) (nhds (P.derivWeierstrassP w)) :=
    ((P.analyticOnNhd_derivWeierstrassP w hw).continuousAt).continuousWithinAt.tendsto
  have hRatio : Tendsto (fun z : ℂ => (P.derivWeierstrassP z - P.derivWeierstrassP w)
      / (P.weierstrassP z - P.weierstrassP w)) (𝓝[≠] w)
      (nhds ((6 * P.weierstrassP w ^ 2 - P.g₂ / 2) / P.derivWeierstrassP w)) := by
    refine (tendsto_slopeRatio P hw h2w).congr ?_
    intro z
    exact (slopeRatio_symm P w z).symm
  have hR : Tendsto (fun z : ℂ => -((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w))
        * (P.weierstrassP (z + w) - P.weierstrassP z) - P.derivWeierstrassP z) (𝓝[≠] w)
      (nhds (-((6 * P.weierstrassP w ^ 2 - P.g₂ / 2) / P.derivWeierstrassP w)
        * (P.weierstrassP (w + w) - P.weierstrassP w) - P.derivWeierstrassP w)) :=
    ((hRatio.neg).mul (hLp.sub hp)).sub hpd
  have hEq : (fun z : ℂ => P.derivWeierstrassP (z + w)) =ᶠ[𝓝[≠] w]
      fun z : ℂ => -((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w))
        * (P.weierstrassP (z + w) - P.weierstrassP z) - P.derivWeierstrassP z := by
    have hL1 : ∀ᶠ z in 𝓝[≠] w, z ∉ P.lattice :=
      mem_nhdsWithin_of_mem_nhds ((P.isClosed_lattice.isOpen_compl).mem_nhds hw)
    have hL2 : ∀ᶠ z in 𝓝[≠] w, z + w ∉ P.lattice := by
      have hopen : IsOpen {z : ℂ | z + w ∉ P.lattice} := by
        have he : {z : ℂ | z + w ∉ P.lattice}
            = (fun z : ℂ => z + w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
        rw [he]
        exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
      exact mem_nhdsWithin_of_mem_nhds (hopen.mem_nhds hww)
    filter_upwards [hL1, hL2, eventually_weierstrassP_ne_self P hw h2w] with z hz1 hz2 hz3
    exact derivWeierstrassP_addition_general P hz1 hw hz2 hz3
  exact tendsto_nhds_unique (hLd.congr' hEq) hR

def derivWeierstrassP_duplication.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(倍加の y 座標の公式——第 657 の z → w 極限。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★**`Point` の倍加との一致**。

★mathlib の `slope` は `x₁ = x₂`・`y₁ ≠ negY` のとき
`(3x₁² + 2a₂x₁ + a₄ − a₁y₁)/(y₁ − negY x₁ y₁)`、`latticeCurve` では
`(3℘w² − g₂/4)/℘′w = ℘″(w)/(2℘′w)`。
☆`addX` は第 650、`addY` は第 659 と一致する。 -/
theorem latticePoint_double (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {w : ℂ}
    (hw : w ∉ P.lattice) (h2w : 2 * w ∉ P.lattice)
    (hww : w + w ∉ P.lattice) :
    (WeierstrassCurve.Affine.Point.some (latticePointX P w) (latticePointY P w)
        (nonsingular_latticePoint P hΔ w hw)
      + WeierstrassCurve.Affine.Point.some (latticePointX P w) (latticePointY P w)
        (nonsingular_latticePoint P hΔ w hw) : (latticeCurve P).toAffine.Point)
      = WeierstrassCurve.Affine.Point.some (latticePointX P (w + w))
          (latticePointY P (w + w)) (nonsingular_latticePoint P hΔ (w + w) hww) := by
  have hpne : P.derivWeierstrassP w ≠ 0 := fun hc =>
    h2w ((derivWeierstrassP_eq_zero_iff P w hw).1 hc)
  have hnegY : (latticeCurve P).toAffine.negY (latticePointX P w) (latticePointY P w)
      = -latticePointY P w := by
    simp [WeierstrassCurve.Affine.negY, latticeCurve]
  have hyne : latticePointY P w
      ≠ (latticeCurve P).toAffine.negY (latticePointX P w) (latticePointY P w) := by
    rw [hnegY]
    intro hc
    refine hpne ?_
    simp only [latticePointY] at hc
    linear_combination hc
  have hxy : ¬(latticePointX P w = latticePointX P w
      ∧ latticePointY P w
        = (latticeCurve P).toAffine.negY (latticePointX P w) (latticePointY P w)) :=
    fun h => hyne h.2
  have hslope : (latticeCurve P).toAffine.slope (latticePointX P w) (latticePointX P w)
      (latticePointY P w) (latticePointY P w)
      = (3 * P.weierstrassP w ^ 2 - P.g₂ / 4) / P.derivWeierstrassP w := by
    rw [WeierstrassCurve.Affine.slope, if_pos rfl, if_neg hyne]
    simp only [latticeCurve, latticePointX, latticePointY, WeierstrassCurve.Affine.negY]
    congr 1
    · ring
    · ring
  have hadd := weierstrassP_duplication P hw h2w
  have hadd' := derivWeierstrassP_duplication P hw h2w
  have hX : (latticeCurve P).toAffine.addX (latticePointX P w) (latticePointX P w)
      ((latticeCurve P).toAffine.slope (latticePointX P w) (latticePointX P w)
        (latticePointY P w) (latticePointY P w))
      = latticePointX P (w + w) := by
    rw [hslope, WeierstrassCurve.Affine.addX]
    simp only [latticeCurve, latticePointX]
    rw [hadd]
    field_simp
    ring
  have hY : (latticeCurve P).toAffine.addY (latticePointX P w) (latticePointX P w)
      (latticePointY P w)
      ((latticeCurve P).toAffine.slope (latticePointX P w) (latticePointX P w)
        (latticePointY P w) (latticePointY P w))
      = latticePointY P (w + w) := by
    rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY, hX, hslope]
    simp only [WeierstrassCurve.Affine.negY, latticeCurve, latticePointX, latticePointY]
    rw [hadd']
    field_simp
    ring
  rw [WeierstrassCurve.Affine.Point.add_some hxy]
  congr 1

/-- ★★★★★★★★★★★★**`℘z = ℘w` なら `z ≡ ±w`**——第 624 の言い換え。 -/
theorem sub_or_add_mem_of_weierstrassP_eq (P : PeriodPair) {z w : ℂ}
    (hz : z ∉ P.lattice) (hw : w ∉ P.lattice)
    (hpz : P.weierstrassP z = P.weierstrassP w) :
    z - w ∈ P.lattice ∨ z + w ∈ P.lattice := by
  have hnw : -w ∉ P.lattice := fun hm => hw (by simpa using neg_mem hm)
  have hsq : P.derivWeierstrassP z ^ 2 = P.derivWeierstrassP w ^ 2 := by
    rw [P.derivWeierstrassP_sq z hz, P.derivWeierstrassP_sq w hw, hpz]
  rcases sq_eq_sq_iff_eq_or_eq_neg.1 hsq with hcase | hcase
  · left
    refine mem_lattice_of_shift_eq P (z - w) hw ?_ ?_ ?_
    · rw [show w + (z - w) = z by ring]; exact hz
    · rw [show w + (z - w) = z by ring]; exact hpz
    · rw [show w + (z - w) = z by ring]; exact hcase
  · right
    refine mem_lattice_of_shift_eq P (z + w) hnw ?_ ?_ ?_
    · rw [show -w + (z + w) = z by ring]; exact hz
    · rw [show -w + (z + w) = z by ring, P.weierstrassP_neg]; exact hpz
    · rw [show -w + (z + w) = z by ring, P.derivWeierstrassP_neg, hcase]

def latticePoint_double.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Point の倍加との一致。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sub_or_add_mem_of_weierstrassP_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘z = ℘w なら z ≡ ±w——単射性の言い換え。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`Φ` は群準同型**

    `Φ(z + w) = Φ(z) + Φ(w)`（すべての `z, w`）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★5 つの場合で尽きる:

1. `z ∈ Λ`・2. `w ∈ Λ` —— 周期性（第 652）
3. `z + w ∈ Λ` —— `w ≡ −z` なので `Point.add_of_Y_eq`
4. `℘z ≠ ℘w` —— 第 658 の `latticePoint_add_general`
5. `℘z = ℘w`（かつ `z+w ∉ Λ`）—— 単射性（第 660 の
   `sub_or_add_mem_of_weierstrassP_eq`）で `z ≡ w`、
   すなわち**倍加**（第 660 の `latticePoint_double`）

★★★★★☆**これで一様化 `ℂ/Λ → E(ℂ)` は全単射（第 603-604・624）かつ群準同型
——群同型である。** ☆どの部品も mathlib に無い。 -/
theorem uniformMap_add (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z w : ℂ) :
    uniformMap P hΔ (z + w) = uniformMap P hΔ z + uniformMap P hΔ w := by
  by_cases hz : z ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hz, zero_add, add_comm z w]
    exact uniformMap_periodic P hΔ w hz
  by_cases hw : w ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hw, add_zero]
    exact uniformMap_periodic P hΔ z hw
  by_cases hzw : z + w ∈ P.lattice
  · rw [uniformMap_of_mem P hΔ hzw, uniformMap_of_notMem P hΔ hz,
      uniformMap_of_notMem P hΔ hw]
    have hpw : P.weierstrassP w = P.weierstrassP z := by
      have h1 : P.weierstrassP ((-z) + (z + w)) = P.weierstrassP (-z) :=
        P.weierstrassP_add_coe (-z) ⟨z + w, hzw⟩
      rw [show (-z) + (z + w) = w by ring, P.weierstrassP_neg] at h1
      exact h1
    have hpdw : P.derivWeierstrassP w = -P.derivWeierstrassP z := by
      have h1 : P.derivWeierstrassP ((-z) + (z + w)) = P.derivWeierstrassP (-z) :=
        P.derivWeierstrassP_add_coe (-z) ⟨z + w, hzw⟩
      rw [show (-z) + (z + w) = w by ring, P.derivWeierstrassP_neg] at h1
      exact h1
    refine (WeierstrassCurve.Affine.Point.add_of_Y_eq ?_ ?_).symm
    · simp only [latticePointX]
      exact hpw.symm
    · simp only [latticePointY, WeierstrassCurve.Affine.negY, latticeCurve, latticePointX]
      rw [hpdw]
      ring
  by_cases hne : P.weierstrassP z - P.weierstrassP w = 0
  · have hpz : P.weierstrassP z = P.weierstrassP w := by linear_combination hne
    rcases sub_or_add_mem_of_weierstrassP_eq P hz hw hpz with hsub | hadd
    · have hzeq : uniformMap P hΔ z = uniformMap P hΔ w := by
        have h1 := uniformMap_periodic P hΔ w hsub
        rw [show w + (z - w) = z by ring] at h1
        exact h1
      have h2w : 2 * w ∉ P.lattice := by
        intro hc
        refine hzw ?_
        have he : z + w = (w + w) + (z - w) := by ring
        rw [he]
        exact P.lattice.add_mem (by simpa [two_mul] using hc) hsub
      have hww : w + w ∉ P.lattice := fun hc => h2w (by simpa [two_mul] using hc)
      have hzw' : uniformMap P hΔ (z + w) = uniformMap P hΔ (w + w) := by
        have h1 := uniformMap_periodic P hΔ (w + w) hsub
        rw [show w + w + (z - w) = z + w by ring] at h1
        exact h1
      rw [hzw', hzeq, uniformMap_of_notMem P hΔ hw, uniformMap_of_notMem P hΔ hww]
      exact (latticePoint_double P hΔ hw h2w hww).symm
    · exact absurd hadd hzw
  · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_notMem P hΔ hw,
      uniformMap_of_notMem P hΔ hzw]
    exact (latticePoint_add_general P hΔ hz hw hzw hne).symm

def uniformMap_add.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は群準同型——一様化は群同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
