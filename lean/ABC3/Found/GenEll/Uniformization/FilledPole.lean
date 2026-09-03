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
import ABC3.Found.GenEll.Uniformization.AdditionODE

/-!
# 一様化 —— 極を埋めた `F_w`

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★極を埋めた `F_w` -/

open Filter Topology in
/-- ★★★★★★★★**極を埋めた `F_w`**——`Ext(z) ≔ limUnder (𝓝[≠] z) F_w`。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★mathlib の `℘` は格子点で junk value を取るので `addDefect` はそこで連続でない。
`elliptic_liouville`（第 598）は `Differentiable` を要求するので、
★★**除去可能特異点を埋めた関数**を作る必要がある（道具は第 630）。 -/
noncomputable def addDefectExt (P : PeriodPair) (w : ℂ) : ℂ → ℂ :=
  fun z => limUnder (𝓝[≠] z) (addDefect P w)

open Filter Topology in
/-- ★★★★★★★★★★**良い点では `Ext` は解析的**（場合 (a)）。 -/
theorem analyticAt_addDefectExt_of_good (P : PeriodPair) (w : ℂ) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    AnalyticAt ℂ (addDefectExt P w) z := by
  refine analyticAt_limUnder_of_analyticAt _ z ?_
  filter_upwards [(isOpen_goodSet P w).mem_nhds ⟨hz, hzw, hne⟩] with y hy
  exact analyticAt_addDefect P w hy.1 hy.2.1 hy.2.2

open Filter Topology in
/-- ★★★★★★★★**良い点では `Ext` は `F_w` そのもの**。 -/
theorem addDefectExt_eq_of_good (P : PeriodPair) (w : ℂ) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    addDefectExt P w z = addDefect P w z :=
  limUnder_eq_self_of_continuousAt _ z (analyticAt_addDefect P w hz hzw hne).continuousAt

def addDefectExt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極を埋めた F_w——Liouville に当てるための整関数)",
    sectionId := "genell-lemma-3-5" }

def analyticAt_addDefectExt_of_good.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い点では Ext は解析的——組み立ての場合 (a)。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology Bornology in
/-- ★★★★★★原点の近く（除いた）では `℘ z ≠ ℘ w`——`℘ z → ∞` だから。 -/
theorem eventually_weierstrassP_ne (P : PeriodPair) (c : ℂ) :
    ∀ᶠ z in 𝓝[≠] (0:ℂ), P.weierstrassP z ≠ c := by
  have hord : meromorphicOrderAt P.weierstrassP 0 < 0 := by
    rw [P.order_weierstrassP 0 P.lattice.zero_mem]; decide
  have h1 : Tendsto P.weierstrassP (𝓝[≠] (0:ℂ)) (cobounded ℂ) :=
    tendsto_cobounded_of_meromorphicOrderAt_neg hord
  exact h1.eventually (eventually_ne_cobounded c)

open Filter Topology in
/-- ★★★★★★★★原点の近く（除いた）はすべて「良い点」（`w ∉ Λ`）。 -/
theorem eventually_good_near_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    ∀ᶠ z in 𝓝[≠] (0:ℂ), z ∉ P.lattice ∧ z + w ∉ P.lattice
      ∧ P.weierstrassP z - P.weierstrassP w ≠ 0 := by
  have h1 : ∀ᶠ z in 𝓝[≠] (0:ℂ), z ∈ ((P.lattice : Set ℂ) \ {0})ᶜ :=
    mem_nhdsWithin_of_mem_nhds (P.isOpen_compl_lattice_sdiff.mem_nhds (by simp))
  have h2 : ∀ᶠ z in 𝓝[≠] (0:ℂ), z + w ∉ P.lattice := by
    have hopen : IsOpen {z : ℂ | z + w ∉ P.lattice} := by
      have he : {z : ℂ | z + w ∉ P.lattice}
          = (fun z : ℂ => z + w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
      rw [he]
      exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
    exact mem_nhdsWithin_of_mem_nhds (hopen.mem_nhds (by simpa using hw))
  filter_upwards [h1, h2, eventually_weierstrassP_ne P (P.weierstrassP w),
    self_mem_nhdsWithin] with z hz1 hz2 hz3 hz4
  refine ⟨fun hc => hz1 ⟨hc, by simpa using hz4⟩, hz2, ?_⟩
  intro hc
  exact hz3 (by linear_combination hc)

open Filter Topology in
/-- ★★★★★★★★★★★★★★**`Ext` は原点で解析的**（場合 (b)）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 610 の `addDefect_eq_near`（`z ≠ 0` で `F_w = addDefectNear`）と
`analyticAt_addDefectNear` に、第 630 の道具を当てる。 -/
theorem analyticAt_addDefectExt_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    AnalyticAt ℂ (addDefectExt P w) 0 := by
  refine analyticAt_limUnder_of_eventuallyEq _ (addDefectNear P w) 0
    (analyticAt_addDefectNear P w hw) ?_ ?_
  · filter_upwards [eventually_good_near_zero P w hw] with z hz
    exact (analyticAt_addDefect P w hz.1 hz.2.1 hz.2.2).continuousAt
  · filter_upwards [eventually_good_near_zero P w hw, self_mem_nhdsWithin] with z hz hz0
    refine addDefect_eq_near P w z ?_ hz.2.2
    simpa using hz0

def analyticAt_addDefectExt_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は原点で解析的——組み立ての場合 (b)。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★`q(t) ≔ 2(℘(t−w) − ℘(w)) − t(℘′(t−w) − ℘′(w))` は `t = 0` で解析的。 -/
theorem analyticAt_addQ (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    AnalyticAt ℂ (addQ P w) 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := by
    intro hc; exact hw (by simpa using neg_mem hc)
  have hf : AnalyticAt ℂ (fun t : ℂ => t - w) 0 := analyticAt_id.sub analyticAt_const
  have hp : AnalyticAt ℂ (fun t : ℂ => P.weierstrassP (t - w)) 0 :=
    AnalyticAt.comp (f := fun t : ℂ => t - w) (x := 0)
      (P.analyticOnNhd_weierstrassP _ hnw) hf
  have hpd : AnalyticAt ℂ (fun t : ℂ => P.derivWeierstrassP (t - w)) 0 :=
    AnalyticAt.comp (f := fun t : ℂ => t - w) (x := 0)
      (P.analyticOnNhd_derivWeierstrassP _ hnw) hf
  show AnalyticAt ℂ (fun t : ℂ => 2 * (P.weierstrassP (t - w) - P.weierstrassP w)
      - t * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)) 0
  exact (analyticAt_const.mul (hp.sub analyticAt_const)).sub
    (analyticAt_id.mul (hpd.sub analyticAt_const))

/-- ★★★★★★★★**`q = t³·g`**（`w ∉ Λ` だけから）——第 626 の仮説を外した形。 -/
theorem exists_addQ_factor' (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    ∃ g : ℂ → ℂ, AnalyticAt ℂ g 0 ∧
      ∀ᶠ t in nhds (0:ℂ), addQ P w t = t ^ 3 * g t :=
  exists_addQ_factor P w hw (analyticAt_addQ P w hw)

def analyticAt_addQ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q は t = 0 で解析的。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★**`Ext` は `−w` で解析的**（場合 (d)）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 627 の `addDefect_eq_nearNeg`（局所形）と第 635・636 の因数分解
（`u = t·û` で `û(0) ≠ 0`、`q = t³·g`）に、第 630 の道具を当てる。 -/
theorem analyticAt_addDefectExt_negW (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    AnalyticAt ℂ (addDefectExt P w) (-w) := by
  obtain ⟨û, hûana, hû0, hûeq⟩ := exists_shiftU_factor P w hw h2w
  obtain ⟨g, hgana, hgeq⟩ := exists_addQ_factor' P w hw
  have hnw : -w ∉ P.lattice := fun hm => hw (by simpa using neg_mem hm)
  have hshift : Tendsto (fun z : ℂ => z + w) (nhds (-w)) (nhds (0:ℂ)) := by
    have ht : Tendsto (fun z : ℂ => z + w) (nhds (-w)) (nhds ((-w) + w)) :=
      (continuous_id.add continuous_const).tendsto _
    simpa using ht
  have hcompA : AnalyticAt ℂ (fun z : ℂ => z + w) (-w) := analyticAt_id.add analyticAt_const
  have hz0 : ((fun z : ℂ => z + w) (-w)) = 0 := by ring
  -- 局所解析関数 A
  have hgA : AnalyticAt ℂ (fun z : ℂ => g (z + w)) (-w) :=
    AnalyticAt.comp (f := fun z : ℂ => z + w) (x := -w) (by simpa using hgana) hcompA
  have hûA : AnalyticAt ℂ (fun z : ℂ => û (z + w)) (-w) :=
    AnalyticAt.comp (f := fun z : ℂ => z + w) (x := -w) (by simpa using hûana) hcompA
  have heA : AnalyticAt ℂ (fun z : ℂ => P.weierstrassPExcept 0 (z + w)) (-w) := by
    refine AnalyticAt.comp (g := P.weierstrassPExcept 0) (f := fun z : ℂ => z + w)
      (x := -w) ?_ hcompA
    have he0 : AnalyticAt ℂ (P.weierstrassPExcept 0) 0 :=
      ((P.differentiableOn_weierstrassPExcept 0).analyticOnNhd
        P.isOpen_compl_lattice_sdiff) 0 (by simp)
    simpa using he0
  have hpA : AnalyticAt ℂ P.weierstrassP (-w) := P.analyticOnNhd_weierstrassP _ hnw
  have hpdA : AnalyticAt ℂ P.derivWeierstrassP (-w) := P.analyticOnNhd_derivWeierstrassP _ hnw
  have hdenne : (fun z : ℂ => 4 * û (z + w) ^ 2) (-w) ≠ 0 := by
    simp only [neg_add_cancel]
    exact mul_ne_zero (by norm_num) (pow_ne_zero _ hû0)
  have hAana : AnalyticAt ℂ (fun z : ℂ =>
      g (z + w) * (2 * û (z + w) + (P.derivWeierstrassP z - P.derivWeierstrassP w))
        / (4 * û (z + w) ^ 2)
      + P.weierstrassPExcept 0 (z + w) + P.weierstrassP z + P.weierstrassP w) (-w) :=
    (((hgA.mul ((analyticAt_const.mul hûA).add (hpdA.sub analyticAt_const))).div
      (analyticAt_const.mul (hûA.pow 2)) hdenne).add heA).add hpA |>.add analyticAt_const
  -- 良い点であることと局所形
  have hgood : ∀ᶠ z in 𝓝[≠] (-w), z ∉ P.lattice ∧ z + w ∉ P.lattice
      ∧ P.weierstrassP z - P.weierstrassP w ≠ 0 := by
    have hL : ∀ᶠ z in 𝓝[≠] (-w), z ∉ P.lattice :=
      mem_nhdsWithin_of_mem_nhds ((P.isClosed_lattice.isOpen_compl).mem_nhds hnw)
    have hM : ∀ᶠ z in 𝓝[≠] (-w), z + w ∈ ((P.lattice : Set ℂ) \ {0})ᶜ := by
      refine mem_nhdsWithin_of_mem_nhds (hshift.eventually ?_)
      exact P.isOpen_compl_lattice_sdiff.mem_nhds (by simp)
    have hU : ∀ᶠ z in 𝓝[≠] (-w),
        P.weierstrassP (z + w - w) - P.weierstrassP w = (z + w) * û (z + w) :=
      mem_nhdsWithin_of_mem_nhds (hshift.eventually hûeq)
    have hUne : ∀ᶠ z in 𝓝[≠] (-w), û (z + w) ≠ 0 :=
      mem_nhdsWithin_of_mem_nhds (hshift.eventually (hûana.continuousAt.eventually_ne hû0))
    filter_upwards [hL, hM, hU, hUne, self_mem_nhdsWithin] with z hz1 hz2 hz3 hz4 hz5
    have hzne : z + w ≠ 0 := by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hz5
      intro hc; exact hz5 (by linear_combination hc)
    refine ⟨hz1, fun hc => hz2 ⟨hc, by simpa using hzne⟩, ?_⟩
    have hzz : z + w - w = z := by ring
    rw [hzz] at hz3
    rw [hz3]
    exact mul_ne_zero hzne hz4
  refine analyticAt_limUnder_of_eventuallyEq _ _ (-w) hAana ?_ ?_
  · filter_upwards [hgood] with z hz
    exact (analyticAt_addDefect P w hz.1 hz.2.1 hz.2.2).continuousAt
  · have hU : ∀ᶠ z in 𝓝[≠] (-w),
        P.weierstrassP (z + w - w) - P.weierstrassP w = (z + w) * û (z + w) :=
      mem_nhdsWithin_of_mem_nhds (hshift.eventually hûeq)
    have hG : ∀ᶠ z in 𝓝[≠] (-w), addQ P w (z + w) = (z + w) ^ 3 * g (z + w) :=
      mem_nhdsWithin_of_mem_nhds (hshift.eventually hgeq)
    have hUne : ∀ᶠ z in 𝓝[≠] (-w), û (z + w) ≠ 0 :=
      mem_nhdsWithin_of_mem_nhds (hshift.eventually (hûana.continuousAt.eventually_ne hû0))
    filter_upwards [hU, hG, hUne, self_mem_nhdsWithin] with z hz1 hz2 hz3 hz5
    have hzne : z + w ≠ 0 := by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hz5
      intro hc; exact hz5 (by linear_combination hc)
    have hkey := addDefect_eq_nearNeg P w û g (P.weierstrassPExcept 0) (z + w) hzne hz3
      hz1 hz2 (by rw [← weierstrassP_sub_invSq]; ring)
    have hzz : z + w - w = z := by ring
    rw [hzz] at hkey
    exact hkey

def analyticAt_addDefectExt_negW.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は −w で解析的——組み立ての場合 (d)。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★`℘′ − ℘′(w)` は `w` で `1` 位以上の零点。 -/
theorem one_le_analyticOrderAt_derivSub (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    ((1 : ℕ) : ℕ∞)
      ≤ analyticOrderAt (fun z : ℂ => P.derivWeierstrassP z - P.derivWeierstrassP w) w := by
  refine (natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero
    ((P.analyticOnNhd_derivWeierstrassP w hw).sub analyticAt_const)).2 ?_
  intro i hi
  interval_cases i
  simp

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★**`Ext` は `w` で解析的**（場合 (c)）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★分母 `℘ − ℘(w)` は `w` でちょうど 1 位の零点（第 629）、
分子 `℘′ − ℘′(w)` も `w` で消えるので、比 `n/d` は解析的。 -/
theorem analyticAt_addDefectExt_atW (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    AnalyticAt ℂ (addDefectExt P w) w := by
  have hsubana : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP z - P.weierstrassP w) w :=
    (P.analyticOnNhd_weierstrassP w hw).sub analyticAt_const
  have hdsubana : AnalyticAt ℂ
      (fun z : ℂ => P.derivWeierstrassP z - P.derivWeierstrassP w) w :=
    (P.analyticOnNhd_derivWeierstrassP w hw).sub analyticAt_const
  obtain ⟨d, hd, hd0, hdeq⟩ := (AnalyticAt.analyticOrderAt_eq_natCast hsubana (n := 1)).1
      (by rw [analyticOrderAt_weierstrassP_sub_self P w hw h2w]; norm_num)
  obtain ⟨n, hn, hneq⟩ := (natCast_le_analyticOrderAt hdsubana (n := 1)).1
      (one_le_analyticOrderAt_derivSub P w hw)
  have hww : w + w ∉ P.lattice := by
    intro hc; exact h2w (by simpa [two_mul] using hc)
  have hshiftA : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) w := by
    have hf : AnalyticAt ℂ (fun z : ℂ => z + w) w := analyticAt_id.add analyticAt_const
    exact AnalyticAt.comp (g := P.weierstrassP) (f := fun z : ℂ => z + w) (x := w)
      (P.analyticOnNhd_weierstrassP _ hww) hf
  have hAana : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w) + P.weierstrassP z
      + P.weierstrassP w - (n z / d z) ^ 2 / 4) w :=
    ((hshiftA.add (P.analyticOnNhd_weierstrassP w hw)).add analyticAt_const).sub
      (((hn.div hd hd0).pow 2).div analyticAt_const (by norm_num))
  have hdne : ∀ᶠ z in 𝓝[≠] w, d z ≠ 0 :=
    mem_nhdsWithin_of_mem_nhds (hd.continuousAt.eventually_ne hd0)
  have hLat : ∀ᶠ z in 𝓝[≠] w, z ∉ P.lattice :=
    mem_nhdsWithin_of_mem_nhds ((P.isClosed_lattice.isOpen_compl).mem_nhds hw)
  have hLat2 : ∀ᶠ z in 𝓝[≠] w, z + w ∉ P.lattice := by
    have hopen : IsOpen {z : ℂ | z + w ∉ P.lattice} := by
      have he : {z : ℂ | z + w ∉ P.lattice}
          = (fun z : ℂ => z + w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
      rw [he]
      exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
    exact mem_nhdsWithin_of_mem_nhds (hopen.mem_nhds hww)
  have hdE : ∀ᶠ z in 𝓝[≠] w, P.weierstrassP z - P.weierstrassP w = (z - w) * d z := by
    filter_upwards [mem_nhdsWithin_of_mem_nhds hdeq] with z hz
    simpa using hz
  have hnE : ∀ᶠ z in 𝓝[≠] w,
      P.derivWeierstrassP z - P.derivWeierstrassP w = (z - w) * n z := by
    filter_upwards [mem_nhdsWithin_of_mem_nhds hneq] with z hz
    simpa using hz
  have hgood : ∀ᶠ z in 𝓝[≠] w, z ∉ P.lattice ∧ z + w ∉ P.lattice
      ∧ P.weierstrassP z - P.weierstrassP w ≠ 0 := by
    filter_upwards [hLat, hLat2, hdE, hdne, self_mem_nhdsWithin] with z h1 h2 h3 h4 h5
    refine ⟨h1, h2, ?_⟩
    rw [h3]
    exact mul_ne_zero (sub_ne_zero.2 (by simpa using h5)) h4
  refine analyticAt_limUnder_of_eventuallyEq _ _ w hAana ?_ ?_
  · filter_upwards [hgood] with z hz
    exact (analyticAt_addDefect P w hz.1 hz.2.1 hz.2.2).continuousAt
  · filter_upwards [hdE, hnE, hdne, self_mem_nhdsWithin] with z h3 h4 h5 h6
    have hzw : z - w ≠ 0 := sub_ne_zero.2 (by simpa using h6)
    show addDefect P w z = _
    simp only [addDefect]
    rw [h3, h4]
    have : (z - w) * n z / ((z - w) * d z) = n z / d z := by
      field_simp
    rw [this]

def analyticAt_addDefectExt_atW.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は w で解析的——組み立ての場合 (c)。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★`𝓝[≠]` は平行移動で移る。 -/
theorem map_add_nhdsNE (z l : ℂ) :
    Filter.map (fun y : ℂ => y + l) (𝓝[≠] z) = 𝓝[≠] (z + l) := by
  have h1 : Filter.map (fun y : ℂ => y + l) (nhds z) = nhds (z + l) :=
    (Homeomorph.addRight l).map_nhds_eq z
  rw [nhdsWithin, nhdsWithin, Filter.map_inf (add_left_injective l), h1]
  congr 1
  rw [Filter.map_principal]
  congr 1
  ext y
  simp only [Set.mem_image, Set.mem_compl_iff, Set.mem_singleton_iff]
  constructor
  · rintro ⟨x, hx, rfl⟩
    intro hc
    exact hx (by linear_combination hc)
  · intro hy
    exact ⟨y - l, fun hc => hy (by linear_combination hc), by ring⟩

/-- ★★★★★★★★**`F_w` は `Λ`-周期的**。 -/
theorem addDefect_periodic (P : PeriodPair) (w : ℂ) (z : ℂ) (l : ℂ) (hl : l ∈ P.lattice) :
    addDefect P w (z + l) = addDefect P w z := by
  simp only [addDefect]
  rw [show z + l + w = (z + w) + l by ring, P.weierstrassP_add_coe (z + w) ⟨l, hl⟩,
    P.weierstrassP_add_coe z ⟨l, hl⟩, P.derivWeierstrassP_add_coe z ⟨l, hl⟩]

open Filter Topology in
/-- ★★★★★★★★★★**`Ext` も `Λ`-周期的**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`limUnder` は `map` で書けるので、`𝓝[≠]` の平行移動と `F_w` の周期性から従う。
☆これで第 634（原点）の解析性が**格子の全点**へ延びる。 -/
theorem addDefectExt_periodic (P : PeriodPair) (w : ℂ) (z : ℂ) (l : ℂ)
    (hl : l ∈ P.lattice) :
    addDefectExt P w (z + l) = addDefectExt P w z := by
  simp only [addDefectExt, Filter.limUnder]
  congr 1
  rw [← map_add_nhdsNE z l, Filter.map_map]
  refine Filter.map_congr ?_
  filter_upwards with y
  simp only [Function.comp_apply]
  exact addDefect_periodic P w y l hl

open Filter Topology in
/-- ★★★★★★★★★★★★**`Ext` は格子の全点で解析的**——周期性で原点から延ばす。 -/
theorem analyticAt_addDefectExt_lattice (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    {p : ℂ} (hp : p ∈ P.lattice) :
    AnalyticAt ℂ (addDefectExt P w) p := by
  have h0 := analyticAt_addDefectExt_zero P w hw
  have hshift : AnalyticAt ℂ (fun z : ℂ => addDefectExt P w (z - p + p)) p := by
    have hf : AnalyticAt ℂ (fun z : ℂ => z - p) p := analyticAt_id.sub analyticAt_const
    have hg : AnalyticAt ℂ (addDefectExt P w) ((fun z : ℂ => z - p) p) := by
      simpa using h0
    have hcomp := AnalyticAt.comp (g := addDefectExt P w) (f := fun z : ℂ => z - p)
      (x := p) hg hf
    refine hcomp.congr ?_
    filter_upwards with z
    simp only [Function.comp_apply]
    exact (addDefectExt_periodic P w (z - p) p hp).symm
  refine hshift.congr ?_
  filter_upwards with z
  ring_nf

def addDefectExt_periodic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は Λ-周期的。★無条件)",
    sectionId := "genell-lemma-3-5" }

def analyticAt_addDefectExt_lattice.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は格子の全点で解析的——周期性で原点から延ばす。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★**解析性は `Λ` の平行移動で移る**。 -/
theorem analyticAt_addDefectExt_of_shift (P : PeriodPair) (w : ℂ) {p q : ℂ}
    (hpq : q - p ∈ P.lattice) (h : AnalyticAt ℂ (addDefectExt P w) p) :
    AnalyticAt ℂ (addDefectExt P w) q := by
  have hf : AnalyticAt ℂ (fun z : ℂ => z - (q - p)) q := analyticAt_id.sub analyticAt_const
  have hg : AnalyticAt ℂ (addDefectExt P w) ((fun z : ℂ => z - (q - p)) q) := by
    simpa using h
  have hcomp := AnalyticAt.comp (g := addDefectExt P w)
    (f := fun z : ℂ => z - (q - p)) (x := q) hg hf
  refine hcomp.congr ?_
  filter_upwards with z
  simp only [Function.comp_apply]
  have := addDefectExt_periodic P w (z - (q - p)) (q - p) hpq
  rw [show z - (q - p) + (q - p) = z by ring] at this
  exact this.symm

open Filter Topology in
/-- ★★★★★★★★★★**`Ext 0 = 0`**——第 610 の `addDefectNear_zero` から。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★☆**これが「加法定理が原点で成り立つ」ことであり、Liouville の定数を決める**。 -/
theorem addDefectExt_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    addDefectExt P w 0 = 0 := by
  have heq : addDefect P w =ᶠ[𝓝[≠] (0:ℂ)] addDefectNear P w := by
    filter_upwards [eventually_good_near_zero P w hw, self_mem_nhdsWithin] with z hz hz0
    exact addDefect_eq_near P w z (by simpa using hz0) hz.2.2
  have h1 : Tendsto (addDefectNear P w) (𝓝[≠] (0:ℂ)) (nhds (addDefectNear P w 0)) :=
    (analyticAt_addDefectNear P w hw).continuousAt.continuousWithinAt.tendsto
  have h2 : addDefectExt P w 0 = addDefectNear P w 0 :=
    (h1.congr' heq.symm).limUnder_eq
  rw [h2, addDefectNear_zero]

def addDefectExt_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext 0 = 0——Liouville の定数を決める。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`Ext` は整関数**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★任意の点が次の 4 通りのいずれかである:
`Λ`（第 639）・`z ≡ −w`（第 637）・`z ≡ w`（第 638）・良い点（第 633）。
☆`℘(z) = ℘(w)` なら `℘′(z)² = ℘′(w)²` なので `℘′(z) = ±℘′(w)`、
どちらの場合も第 624 の単射性で `z ≡ ±w` になる。 -/
theorem differentiable_addDefectExt (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    Differentiable ℂ (addDefectExt P w) := by
  intro z
  refine AnalyticAt.differentiableAt ?_
  by_cases hzl : z ∈ P.lattice
  · exact analyticAt_addDefectExt_lattice P w hw hzl
  by_cases hzw : z + w ∈ P.lattice
  · refine analyticAt_addDefectExt_of_shift P w (p := -w) ?_
      (analyticAt_addDefectExt_negW P w hw h2w)
    simpa using hzw
  by_cases hne : P.weierstrassP z - P.weierstrassP w = 0
  · have hpz : P.weierstrassP z = P.weierstrassP w := by linear_combination hne
    have hsq : P.derivWeierstrassP z ^ 2 = P.derivWeierstrassP w ^ 2 := by
      rw [P.derivWeierstrassP_sq z hzl, P.derivWeierstrassP_sq w hw, hpz]
    have hnw : -w ∉ P.lattice := fun hm => hw (by simpa using neg_mem hm)
    rcases sq_eq_sq_iff_eq_or_eq_neg.1 hsq with hcase | hcase
    · have hmem : z - w ∈ P.lattice := by
        refine mem_lattice_of_shift_eq P (z - w) hw ?_ ?_ ?_
        · rw [show w + (z - w) = z by ring]; exact hzl
        · rw [show w + (z - w) = z by ring]; exact hpz
        · rw [show w + (z - w) = z by ring]; exact hcase
      exact analyticAt_addDefectExt_of_shift P w (p := w) hmem
        (analyticAt_addDefectExt_atW P w hw h2w)
    · exfalso
      refine hzw ?_
      have hmem := mem_lattice_of_shift_eq P (z + w) (z₀ := -w) hnw
        (by rw [show -w + (z + w) = z by ring]; exact hzl)
        (by rw [show -w + (z + w) = z by ring, P.weierstrassP_neg]; exact hpz)
        (by rw [show -w + (z + w) = z by ring, P.derivWeierstrassP_neg, hcase])
      exact hmem
  · exact analyticAt_addDefectExt_of_good P w hzl hzw hne

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`F_w ≡ 0`**——Liouville で閉じる。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Ext` は整で `Λ`-周期的（第 639）なので `elliptic_liouville`（第 598）で定数、
その値は `Ext 0 = 0`（第 640）。★★良い点では `Ext = F_w`（第 633）である。 -/
theorem addDefect_eq_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    addDefect P w z = 0 := by
  have hper : ∀ (y : ℂ), ∀ l ∈ P.lattice, addDefectExt P w (y + l) = addDefectExt P w y :=
    fun y l hl => addDefectExt_periodic P w y l hl
  have hconst := elliptic_liouville P (addDefectExt P w)
    (differentiable_addDefectExt P w hw h2w) hper z 0
  rw [addDefectExt_zero P w hw] at hconst
  rw [← addDefectExt_eq_of_good P w hz hzw hne]
  exact hconst

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`℘` の加法定理**

    `℘(z+w) = (1/4)·((℘′(z) − ℘′(w))/(℘(z) − ℘(w)))² − ℘(z) − ℘(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★★☆**mathlib に無い**（`Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` は
`℘` の理論を 1080 行ぶん持つが加法定理は無い、2026-08-29 に測定）。

★★機構は Liouville（第 598）で、極は 4 か所とも塞がっている:
`Λ`（第 610・639）・`z ≡ w`（第 638）・`z ≡ −w`（第 637）・
その他は単射性（第 624-625）で存在しない。
☆★単射性は**零点勘定を使わず**、線型 2 階 ODE `h″ = 6(℘(·+a)+℘)h` の
一意性（第 622-624）から出ている。 -/
theorem weierstrassP_addition (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    P.weierstrassP (z + w)
      = ((P.derivWeierstrassP z - P.derivWeierstrassP w)
          / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4
        - P.weierstrassP z - P.weierstrassP w := by
  have h := addDefect_eq_zero P w hw h2w hz hzw hne
  simp only [addDefect] at h
  linear_combination h

def differentiable_addDefectExt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Ext は整関数——4 通りの場合分けが尽きる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_addition.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の加法定理。★無条件——mathlib に無い)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
