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
import ABC3.Found.GenEll.Uniformization.AdditionEntry

/-!
# 一様化 —— ODE の一意性で零点勘定を回避する・因数分解・極を埋める道具

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★ODE の一意性で零点勘定を回避する -/

/-- ★★★★★平行移動した `℘` の微分。 -/
theorem hasDerivAt_weierstrassP_shift (P : PeriodPair) (a : ℂ) {z : ℂ}
    (hz : z + a ∉ P.lattice) :
    HasDerivAt (fun s : ℂ => P.weierstrassP (s + a)) (P.derivWeierstrassP (z + a)) z :=
  HasDerivAt.comp_add_const z a (hasDerivAt_weierstrassP P hz)

/-- ★★★★★平行移動した `℘′` の微分。 -/
theorem hasDerivAt_derivWeierstrassP_shift (P : PeriodPair) (a : ℂ) {z : ℂ}
    (hz : z + a ∉ P.lattice) :
    HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s + a))
      (deriv P.derivWeierstrassP (z + a)) z :=
  HasDerivAt.comp_add_const z a (hasDerivAt_derivWeierstrassP P hz)

/-- ★★★★`{s | s ∉ Λ ∧ s + a ∉ Λ}` は開集合。 -/
theorem isOpen_shiftDomain (P : PeriodPair) (a : ℂ) :
    IsOpen {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} := by
  have h1 : IsOpen {s : ℂ | s ∉ P.lattice} := P.isClosed_lattice.isOpen_compl
  have h2 : IsOpen {s : ℂ | s + a ∉ P.lattice} := by
    have he : {s : ℂ | s + a ∉ P.lattice}
        = (fun s : ℂ => s + a) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
    rw [he]
    exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
  exact h1.inter h2

/-- ★★★★★★**差 `h(z) ≔ ℘(z+a) − ℘(z)` の 1 階導関数**。 -/
theorem hasDerivAt_shiftDiff (P : PeriodPair) (a : ℂ) {z : ℂ}
    (hz : z ∉ P.lattice) (hza : z + a ∉ P.lattice) :
    HasDerivAt (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s)
      (P.derivWeierstrassP (z + a) - P.derivWeierstrassP z) z :=
  (hasDerivAt_weierstrassP_shift P a hza).sub (hasDerivAt_weierstrassP P hz)

/-- ★★★★★★★★★★★★★★★★★★★★★★**差は線型の 2 階 ODE を満たす**

    `h″(z) = 6·(℘(z+a) + ℘(z))·h(z)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Found/GenEll/WeierstrassODE.lean` の `deriv_derivWeierstrassP`
（`℘″ = 6℘² − g₂/2`）を `℘(·+a)` と `℘` の両方に当てて引き算するだけである
——`g₂/2` が消えて `6(u² − v²) = 6(u+v)(u−v)` になる。

★★★☆**これが零点勘定を回避する鍵である**。`h(z₂) = h′(z₂) = 0` なら、
解析的位数の算術で `h` は `z₂` の近傍で恒等的に `0` になる:
位数を `n`（有限、`≥ 2`）とすると `h″` の位数は `n − 2`、
一方 `6(u+v)·h` の位数は `≥ n`——矛盾。 -/
theorem deriv_deriv_shiftDiff (P : PeriodPair) (a : ℂ) {z : ℂ}
    (hz : z ∉ P.lattice) (hza : z + a ∉ P.lattice) :
    deriv (deriv (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s)) z
      = 6 * (P.weierstrassP (z + a) + P.weierstrassP z)
          * (P.weierstrassP (z + a) - P.weierstrassP z) := by
  have hnhds : {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} ∈ nhds z :=
    (isOpen_shiftDomain P a).mem_nhds ⟨hz, hza⟩
  have heq : deriv (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s)
      =ᶠ[nhds z] fun s : ℂ => P.derivWeierstrassP (s + a) - P.derivWeierstrassP s := by
    filter_upwards [hnhds] with s hs
    exact (hasDerivAt_shiftDiff P a hs.1 hs.2).deriv
  rw [heq.deriv_eq]
  have h2 : HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s + a) - P.derivWeierstrassP s)
      (deriv P.derivWeierstrassP (z + a) - deriv P.derivWeierstrassP z) z :=
    (hasDerivAt_derivWeierstrassP_shift P a hza).sub (hasDerivAt_derivWeierstrassP P hz)
  rw [h2.deriv, deriv_derivWeierstrassP P hza, deriv_derivWeierstrassP P hz]
  ring

def deriv_deriv_shiftDiff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘(z+a) − ℘(z) は線型の 2 階 ODE h″ = 6(u+v)h を満たす。★無条件)",
    sectionId := "genell-lemma-3-5" }

def deriv_deriv_shiftDiff.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★到達点(2026-08-29、第 622): これで零点勘定を回避できる見込みが立った。" ++
       "☆h(z₂) = h′(z₂) = 0 なら解析的位数の算術で h は z₂ の近傍で恒等的に 0:" ++
       "位数を n(有限、≥ 2)とすると h″ の位数は n − 2(analyticOrderAt_deriv_of_pos を 2 回)、" ++
       "一方 6(u+v)·h の位数は ≥ n(analyticOrderAt_smul)——矛盾。" ++
       "★したがって位数は ⊤ で h は近傍で 0。" ++
       "★★一致の定理で ℂ ∖ (Λ ∪ (−a+Λ)) 全体へ延ばし、" ++
       "mem_lattice_of_weierstrassP_periodic(第 620)で a ∈ Λ、すなわち単射性") 13 ]

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★**線型 2 階 ODE の一意性（解析的位数版）**。

`h(z₀) = h′(z₀) = 0` かつ `h″ = c·h`（`c` は `z₀` で解析的）なら `h` は `z₀` の近傍で `0`。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★機構は**解析的位数の算術だけ**——留数定理も偏角の原理も要らない:

* `h(z₀) = h′(z₀) = 0` から位数は `≥ 2`
* 位数が有限 `= m + 2` なら `h″` の位数は `m`（`analyticOrderAt_deriv_of_pos` を 2 回）
* 一方 `c·h` の位数は `order(c) + (m+2) ≥ m + 2`（`analyticOrderAt_smul`）
* `m ≥ m + 2` は偽——★したがって位数は `⊤`、すなわち `h` は近傍で `0` -/
theorem eventually_eq_zero_of_linear_ode (h c : ℂ → ℂ) (z₀ : ℂ)
    (hana : AnalyticAt ℂ h z₀) (hc : AnalyticAt ℂ c z₀)
    (h0 : h z₀ = 0) (h1 : deriv h z₀ = 0)
    (hode : deriv (deriv h) =ᶠ[nhds z₀] c • h) :
    ∀ᶠ s in nhds z₀, h s = 0 := by
  rw [← analyticOrderAt_eq_top]
  by_contra hfin
  obtain ⟨n, hn0⟩ := Option.ne_none_iff_exists'.1 hfin
  have hn : analyticOrderAt h z₀ = (n : ℕ∞) := hn0
  have h2le : ((2 : ℕ) : ℕ∞) ≤ analyticOrderAt h z₀ := by
    rw [natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hana]
    intro i hi
    interval_cases i
    · simpa using h0
    · rw [iteratedDeriv_one]; exact h1
  rw [hn] at h2le
  have h2n : 2 ≤ n := by exact_mod_cast h2le
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 2 := ⟨n - 2, by omega⟩
  have hd1 : analyticOrderAt (deriv h) z₀ = ((m + 1 : ℕ) : ℕ∞) := by
    refine analyticOrderAt_deriv_of_pos hana ?_
    rw [hn]; push_cast; ring
  have hd2 : analyticOrderAt (deriv (deriv h)) z₀ = ((m : ℕ) : ℕ∞) := by
    refine analyticOrderAt_deriv_of_pos hana.deriv ?_
    rw [hd1]; push_cast; ring
  rw [analyticOrderAt_congr hode, analyticOrderAt_smul hc hana, hn] at hd2
  have hfinal : (((m + 2 : ℕ)) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by
    calc (((m + 2 : ℕ)) : ℕ∞) ≤ analyticOrderAt c z₀ + (((m + 2 : ℕ)) : ℕ∞) := le_add_self
    _ = ((m : ℕ) : ℕ∞) := hd2
  have hle : (m + 2 : ℕ) ≤ m := by exact_mod_cast hfinal
  omega

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**`℘(z₀+a) = ℘(z₀)` かつ
`℘′(z₀+a) = ℘′(z₀)` なら `℘(·+a) = ℘` は `z₀` の近傍で一致する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**これが一様化の単射性の核である**——第 622 の線型 2 階 ODE
`h″ = 6(℘(·+a) + ℘)·h` に上の `eventually_eq_zero_of_linear_ode` を当てる。

☆残るのは (1) 一致の定理で `ℂ ∖ (Λ ∪ (−a+Λ))` 全体へ延ばすこと、
(2) 第 620 の `mem_lattice_of_weierstrassP_periodic` で `a ∈ Λ` を出すこと。 -/
theorem weierstrassP_shift_eventually_eq (P : PeriodPair) (a : ℂ) {z₀ : ℂ}
    (hz : z₀ ∉ P.lattice) (hza : z₀ + a ∉ P.lattice)
    (h0 : P.weierstrassP (z₀ + a) = P.weierstrassP z₀)
    (h1 : P.derivWeierstrassP (z₀ + a) = P.derivWeierstrassP z₀) :
    ∀ᶠ s in nhds z₀, P.weierstrassP (s + a) - P.weierstrassP s = 0 := by
  set h : ℂ → ℂ := fun s => P.weierstrassP (s + a) - P.weierstrassP s with hh
  set c : ℂ → ℂ := fun s => 6 * (P.weierstrassP (s + a) + P.weierstrassP s) with hcdef
  have hshiftAna : AnalyticAt ℂ (fun s : ℂ => P.weierstrassP (s + a)) z₀ := by
    have hf : AnalyticAt ℂ (fun s : ℂ => s + a) z₀ := analyticAt_id.add analyticAt_const
    have hg : AnalyticAt ℂ P.weierstrassP ((fun s : ℂ => s + a) z₀) :=
      P.analyticOnNhd_weierstrassP (z₀ + a) hza
    exact AnalyticAt.comp (f := fun s : ℂ => s + a) (x := z₀) hg hf
  have hpAna : AnalyticAt ℂ P.weierstrassP z₀ := P.analyticOnNhd_weierstrassP z₀ hz
  have hana : AnalyticAt ℂ h z₀ := hshiftAna.sub hpAna
  have hcAna : AnalyticAt ℂ c z₀ := analyticAt_const.mul (hshiftAna.add hpAna)
  have hh0 : h z₀ = 0 := by simp only [hh]; rw [h0]; ring
  have hh1 : deriv h z₀ = 0 := by
    rw [(hasDerivAt_shiftDiff P a hz hza).deriv, h1]
    ring
  have hnhds : {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} ∈ nhds z₀ :=
    (isOpen_shiftDomain P a).mem_nhds ⟨hz, hza⟩
  have hode : deriv (deriv h) =ᶠ[nhds z₀] c • h := by
    filter_upwards [hnhds] with s hs
    show deriv (deriv (fun t : ℂ => P.weierstrassP (t + a) - P.weierstrassP t)) s
        = 6 * (P.weierstrassP (s + a) + P.weierstrassP s)
          * (P.weierstrassP (s + a) - P.weierstrassP s)
    exact deriv_deriv_shiftDiff P a hs.1 hs.2
  exact eventually_eq_zero_of_linear_ode h c z₀ hana hcAna hh0 hh1 hode

def eventually_eq_zero_of_linear_ode.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(線型 2 階 ODE の一意性——解析的位数の算術だけ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_shift_eventually_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化の単射性の核——℘ と ℘′ が一致すれば近傍で一致。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★`{s | s ∉ Λ ∧ s + a ∉ Λ}` は連結——2 つの可算集合の補集合だから。 -/
theorem isPreconnected_shiftDomain (P : PeriodPair) (a : ℂ) :
    IsPreconnected {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} := by
  have hcount : ((P.lattice : Set ℂ) ∪ ((fun w : ℂ => w - a) '' (P.lattice : Set ℂ))).Countable := by
    refine Set.Countable.union ?_ (Set.Countable.image ?_ _) <;>
      exact countable_of_Lindelof_of_discrete (X := P.lattice)
  have hset : {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice}
      = ((P.lattice : Set ℂ) ∪ ((fun w : ℂ => w - a) '' (P.lattice : Set ℂ)))ᶜ := by
    ext s
    simp only [Set.mem_setOf_eq, Set.mem_compl_iff, Set.mem_union, Set.mem_image,
      not_or, not_exists, not_and]
    constructor
    · rintro ⟨h1, h2⟩
      refine ⟨h1, ?_⟩
      rintro w hw rfl
      exact h2 (by simpa using hw)
    · rintro ⟨h1, h2⟩
      refine ⟨h1, fun hc => ?_⟩
      exact h2 (s + a) hc (by ring)
  rw [hset]
  exact (Set.Countable.isConnected_compl_of_one_lt_rank (by simp) hcount).2

open Filter Topology Bornology Metric in
/-- ★★★★★★★★★★**`℘` の周期群はちょうど `Λ`（局所版）**——`−a` の近くで
`℘(z+a) = ℘(z)` が成り立てば `a ∈ Λ`。

★第 620 の `mem_lattice_of_weierstrassP_periodic` を「`−a` の除いた近傍で」に弱めた形。
☆一致の定理から出てくるのはこちらの形である。 -/
theorem mem_lattice_of_eventually_shift_eq (P : PeriodPair) (a : ℂ)
    (h : ∀ᶠ z in 𝓝[≠] (-a), P.weierstrassP (z + a) = P.weierstrassP z) :
    a ∈ P.lattice := by
  by_contra hc
  have hna : -a ∉ P.lattice := fun hm => hc (by simpa using neg_mem hm)
  have hcont : ContinuousAt P.weierstrassP (-a) :=
    (P.analyticOnNhd_weierstrassP (-a) hna).continuousAt
  have hord : meromorphicOrderAt P.weierstrassP 0 < 0 := by
    rw [P.order_weierstrassP 0 P.lattice.zero_mem]; decide
  have h1 : Tendsto P.weierstrassP (𝓝[≠] (0:ℂ)) (cobounded ℂ) :=
    tendsto_cobounded_of_meromorphicOrderAt_neg hord
  have hshift : Tendsto (fun z : ℂ => z + a) (𝓝[≠] (-a)) (𝓝[≠] (0:ℂ)) := by
    rw [tendsto_nhdsWithin_iff]
    refine ⟨?_, ?_⟩
    · have ht : Tendsto (fun z : ℂ => z + a) (𝓝 (-a)) (𝓝 ((-a) + a)) :=
        (continuous_id.add continuous_const).tendsto _
      simpa using ht.mono_left nhdsWithin_le_nhds
    · filter_upwards [self_mem_nhdsWithin] with z hz
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hz ⊢
      intro hcc
      exact hz (by linear_combination hcc)
  have h2 : Tendsto (fun z : ℂ => P.weierstrassP (z + a)) (𝓝[≠] (-a)) (cobounded ℂ) :=
    h1.comp hshift
  have h3 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (cobounded ℂ) := h2.congr' h
  have h4 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (𝓝 (P.weierstrassP (-a))) :=
    hcont.continuousWithinAt
  exact (h4.not_tendsto (disjoint_nhds_cobounded _)) h3

def isPreconnected_shiftDomain.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5({s | s ∉ Λ ∧ s + a ∉ Λ} は連結。★無条件)",
    sectionId := "genell-lemma-3-5" }

def mem_lattice_of_eventually_shift_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の周期群はちょうど Λ——局所版。★無条件)",
    sectionId := "genell-lemma-3-5" }

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**一様化の単射性**

    `℘(z₀+a) = ℘(z₀)` かつ `℘′(z₀+a) = ℘′(z₀)` なら `a ∈ Λ`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**零点勘定（偏角の原理）を使わずに出た**。機構は 3 段:

1. `h ≔ ℘(·+a) − ℘` は線型 2 階 ODE `h″ = 6(℘(·+a) + ℘)·h` を満たす（第 622）
2. `h(z₀) = h′(z₀) = 0` なので解析的位数の算術で `h` は `z₀` の近傍で `0`（第 623）
3. 一致の定理で `{s | s ∉ Λ ∧ s+a ∉ Λ}`（連結）全体へ延ばし、
   `−a` の近くで `℘(z+a) = ℘(z)`。★`℘` は `−a` で解析的なのに
   `℘(z+a)` は `z → −a` で発散する——`a ∉ Λ` なら矛盾

★★★★これで **`(℘, ℘′/2) : ℂ/Λ → E(ℂ)` は単射**であり、
第 603-604 の全射性と合わせて**一様化は全単射**である。 -/
theorem mem_lattice_of_shift_eq (P : PeriodPair) (a : ℂ) {z₀ : ℂ}
    (hz : z₀ ∉ P.lattice) (hza : z₀ + a ∉ P.lattice)
    (h0 : P.weierstrassP (z₀ + a) = P.weierstrassP z₀)
    (h1 : P.derivWeierstrassP (z₀ + a) = P.derivWeierstrassP z₀) :
    a ∈ P.lattice := by
  by_contra hc
  have hna : -a ∉ P.lattice := fun hm => hc (by simpa using neg_mem hm)
  have hana : AnalyticOnNhd ℂ (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s)
      {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} := by
    intro s hs
    have hf : AnalyticAt ℂ (fun t : ℂ => t + a) s := analyticAt_id.add analyticAt_const
    have hg : AnalyticAt ℂ P.weierstrassP ((fun t : ℂ => t + a) s) :=
      P.analyticOnNhd_weierstrassP (s + a) hs.2
    exact (AnalyticAt.comp (f := fun t : ℂ => t + a) (x := s) hg hf).sub
      (P.analyticOnNhd_weierstrassP s hs.1)
  have hloc : (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s) =ᶠ[nhds z₀] 0 :=
    weierstrassP_shift_eventually_eq P a hz hza h0 h1
  have hEq : Set.EqOn (fun s : ℂ => P.weierstrassP (s + a) - P.weierstrassP s) 0
      {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} :=
    hana.eqOn_zero_of_preconnected_of_eventuallyEq_zero
      (isPreconnected_shiftDomain P a) ⟨hz, hza⟩ hloc
  have hVopen : IsOpen ({s : ℂ | s ∉ P.lattice}
      ∩ {s : ℂ | s + a ∉ (P.lattice : Set ℂ) \ {0}}) := by
    refine IsOpen.inter (P.isClosed_lattice.isOpen_compl) ?_
    have he : {s : ℂ | s + a ∉ (P.lattice : Set ℂ) \ {0}}
        = (fun s : ℂ => s + a) ⁻¹' (((P.lattice : Set ℂ) \ {0})ᶜ) := rfl
    rw [he]
    exact P.isOpen_compl_lattice_sdiff.preimage (by fun_prop)
  have hVmem : (-a) ∈ ({s : ℂ | s ∉ P.lattice}
      ∩ {s : ℂ | s + a ∉ (P.lattice : Set ℂ) \ {0}}) := by
    refine ⟨hna, ?_⟩
    simp
  have hev : ∀ᶠ z in 𝓝[≠] (-a), P.weierstrassP (z + a) = P.weierstrassP z := by
    filter_upwards [mem_nhdsWithin_of_mem_nhds (hVopen.mem_nhds hVmem), self_mem_nhdsWithin]
      with z hzV hzne
    have hzne' : z + a ≠ 0 := by
      simp only [Set.mem_compl_iff, Set.mem_singleton_iff] at hzne
      intro hcc
      exact hzne (by linear_combination hcc)
    have hzU : z ∈ {s : ℂ | s ∉ P.lattice ∧ s + a ∉ P.lattice} :=
      ⟨hzV.1, fun hcc => hzV.2 ⟨hcc, by simpa using hzne'⟩⟩
    have hh := hEq hzU
    simp only [Pi.zero_apply] at hh
    linear_combination hh
  exact hc (mem_lattice_of_eventually_shift_eq P a hev)

def mem_lattice_of_shift_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化の単射性——零点勘定を使わずに出た。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★★★**`℘′(w) = 0` なら `w` は 2-捩れ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★第 624 の単射性の系である: `℘′(w) = 0` なら `℘(−w) = ℘(w)`（偶）かつ
`℘′(−w) = −℘′(w) = 0 = ℘′(w)`（奇）なので、`a = −2w` に単射性を当てて `−2w ∈ Λ`。

☆古典的には「`℘′` はちょうど 3 つの零点（非自明な 2-捩れ点）を持つ」という
**零点勘定**から出る事実であるが、★**ODE の一意性から零点勘定なしで出た**。 -/
theorem two_mem_lattice_of_derivWeierstrassP_eq_zero (P : PeriodPair) (w : ℂ)
    (hw : w ∉ P.lattice) (h : P.derivWeierstrassP w = 0) : 2 * w ∈ P.lattice := by
  have hnw : w + (-2 * w) ∉ P.lattice := by
    have he : w + (-2 * w) = -w := by ring
    rw [he]
    exact fun hm => hw (by simpa using neg_mem hm)
  have h0 : P.weierstrassP (w + (-2 * w)) = P.weierstrassP w := by
    have he : w + (-2 * w) = -w := by ring
    rw [he, P.weierstrassP_neg]
  have h1 : P.derivWeierstrassP (w + (-2 * w)) = P.derivWeierstrassP w := by
    have he : w + (-2 * w) = -w := by ring
    rw [he, P.derivWeierstrassP_neg, h]
    ring
  have hm := mem_lattice_of_shift_eq P (-2 * w) hw hnw h0 h1
  have he : (2 : ℂ) * w = -(-2 * w) := by ring
  rw [he]
  exact neg_mem hm

/-- ★★★★★★★★★★★★★★★★★★★★★★**2-捩れの完全な特徴づけ**

    `℘′(w) = 0  ⟺  2w ∈ Λ`（`w ∉ Λ`）

★`⟸` は第 605（`℘′` が奇であることから 3 行）、`⟹` は第 624 の単射性の系。
☆★**mathlib に無い**（`Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` は
`℘′` の零点を同定していない、2026-08-29 に測定）。 -/
theorem derivWeierstrassP_eq_zero_iff (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    P.derivWeierstrassP w = 0 ↔ 2 * w ∈ P.lattice :=
  ⟨two_mem_lattice_of_derivWeierstrassP_eq_zero P w hw,
   derivWeierstrassP_eq_zero_of_two_mem P w⟩

def two_mem_lattice_of_derivWeierstrassP_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘′(w) = 0 なら w は 2-捩れ——単射性の系。★無条件)",
    sectionId := "genell-lemma-3-5" }

def derivWeierstrassP_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2-捩れの完全な特徴づけ ℘′(w) = 0 ⟺ 2w ∈ Λ。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★`z ≡ −w` の側——因数分解 -/

/-- ★★★★★`u(t) ≔ ℘(t−w) − ℘(w)` は `t = 0` で解析的。 -/
theorem analyticAt_shiftU (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    AnalyticAt ℂ (fun t : ℂ => P.weierstrassP (t - w) - P.weierstrassP w) 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := by
    intro hc; exact hw (by simpa using neg_mem hc)
  have hf : AnalyticAt ℂ (fun t : ℂ => t - w) 0 := analyticAt_id.sub analyticAt_const
  have hg : AnalyticAt ℂ P.weierstrassP ((fun t : ℂ => t - w) 0) :=
    P.analyticOnNhd_weierstrassP _ hnw
  exact (AnalyticAt.comp (f := fun t : ℂ => t - w) (x := 0) hg hf).sub analyticAt_const

/-- ★★★★★★★★★★★★★★**`u(t) ≔ ℘(t−w) − ℘(w)` は `t = 0` でちょうど 1 位の零点**
（`2w ∉ Λ` のとき）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`u(0) = ℘(−w) − ℘(w) = 0`（偶）、`u′(0) = ℘′(−w) = −℘′(w) ≠ 0`
（第 625 の `derivWeierstrassP_eq_zero_iff` で `2w ∉ Λ` から `℘′(w) ≠ 0`）。

☆これが `z ≡ −w` での極の打ち消しに要る `û ≔ u/t` の可逆性を与える。 -/
theorem analyticOrderAt_shiftU (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    analyticOrderAt (fun t : ℂ => P.weierstrassP (t - w) - P.weierstrassP w) 0 = 1 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := by
    intro hc; exact hw (by simpa using neg_mem hc)
  have hana := analyticAt_shiftU P w hw
  refine hana.analyticOrderAt_eq_one_of_zero_deriv_ne_zero ?_ ?_
  · simp only [zero_sub, P.weierstrassP_neg]
    ring
  · have hder : HasDerivAt (fun t : ℂ => P.weierstrassP (t - w) - P.weierstrassP w)
        (P.derivWeierstrassP ((0:ℂ) - w)) 0 :=
      (HasDerivAt.comp_sub_const 0 w (hasDerivAt_weierstrassP P hnw)).sub_const _
    rw [hder.deriv]
    have hne : P.derivWeierstrassP w ≠ 0 := by
      intro hcc
      exact h2w ((derivWeierstrassP_eq_zero_iff P w hw).1 hcc)
    simp only [zero_sub, P.derivWeierstrassP_neg]
    exact neg_ne_zero.2 hne

/-- ★★★★★★★★**`u = t·û`（`û` は解析的、`û(0) ≠ 0`）**。 -/
theorem exists_shiftU_factor (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (h2w : 2 * w ∉ P.lattice) :
    ∃ û : ℂ → ℂ, AnalyticAt ℂ û 0 ∧ û 0 ≠ 0 ∧
      ∀ᶠ t in nhds (0:ℂ), P.weierstrassP (t - w) - P.weierstrassP w = t * û t := by
  have hord : analyticOrderAt (fun t : ℂ => P.weierstrassP (t - w) - P.weierstrassP w) 0
      = ((1 : ℕ) : ℕ∞) := by
    rw [analyticOrderAt_shiftU P w hw h2w]
    norm_num
  obtain ⟨g, hg, hg0, hgeq⟩ :=
    (AnalyticAt.analyticOrderAt_eq_natCast (analyticAt_shiftU P w hw)).1 hord
  refine ⟨g, hg, hg0, ?_⟩
  filter_upwards [hgeq] with t ht
  simpa using ht

/-- ★★★★★★★★★★**`q = t³·g`（`g` は解析的）**——第 614 の位数 `≥ 3` から。

☆これと `u = t·û` を合わせると

    `F_w(t−w) = g·(2û+v)/(4û²) + e(t) + ℘(t−w) + ℘(w)`

となり、★`û(0) = −℘′(w) ≠ 0` なので**右辺は `t = 0` で解析的**
——`z ≡ −w` の極が消える。 -/
theorem exists_addQ_factor (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (hana : AnalyticAt ℂ (addQ P w) 0) :
    ∃ g : ℂ → ℂ, AnalyticAt ℂ g 0 ∧
      ∀ᶠ t in nhds (0:ℂ), addQ P w t = t ^ 3 • g t := by
  obtain ⟨g, hg, hgeq⟩ :=
    (natCast_le_analyticOrderAt hana).1 (three_le_analyticOrderAt_addQ P w hw hana)
  exact ⟨g, hg, by simpa using hgeq⟩

def analyticOrderAt_shiftU.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘(t−w) − ℘(w) は t = 0 でちょうど 1 位の零点。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_addQ_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q = t³·g——z ≡ −w の極が消える理由)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★★★★★**`z ≡ −w` の極が消える（等式）**

    `F_w(t−w) = g·(2û + v)/(4û²) + e(t) + ℘(t−w) + ℘(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★導出は 4 行:

    `F_w(t−w) = 1/t² + e + ℘(t−w) + ℘(w) − v²/(4t²û²)`   （`u = tû`）
    `         = [4û² − v²]/(4t²û²) + e + ℘(t−w) + ℘(w)`
    `4û² − v² = (2û−v)(2û+v)`,  `2û − v = q/t = t²g`      （`q = t³g`）
    `         = g(2û+v)/(4û²) + e + ℘(t−w) + ℘(w)`

☆右辺は `û(0) ≠ 0`（第 626）なので `t = 0` で解析的——**極が消える**。 -/
theorem addDefect_eq_nearNeg (P : PeriodPair) (w : ℂ) (û g e : ℂ → ℂ) (t : ℂ)
    (ht : t ≠ 0) (hû : û t ≠ 0)
    (hu : P.weierstrassP (t - w) - P.weierstrassP w = t * û t)
    (hq : addQ P w t = t ^ 3 * g t)
    (he : P.weierstrassP t = 1 / t ^ 2 + e t) :
    addDefect P w (t - w)
      = g t * (2 * û t + (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w))
          / (4 * û t ^ 2)
        + e t + P.weierstrassP (t - w) + P.weierstrassP w := by
  have ht2 : t ^ 2 ≠ 0 := pow_ne_zero _ ht
  have hune : P.weierstrassP (t - w) - P.weierstrassP w ≠ 0 := by
    rw [hu]; exact mul_ne_zero ht hû
  -- `2û − v = t²g`
  have hkey : 2 * û t - (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)
      = t ^ 2 * g t := by
    have h1 : t * (2 * û t - (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w))
        = t * (t ^ 2 * g t) := by
      have h2 : addQ P w t
          = t * (2 * û t - (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)) := by
        simp only [addQ]
        linear_combination 2 * hu
      rw [← h2, hq]
      ring
    exact mul_left_cancel₀ ht h1
  have hV : P.derivWeierstrassP (t - w) - P.derivWeierstrassP w
      = 2 * û t - t ^ 2 * g t := by linear_combination -hkey
  simp only [addDefect]
  have hzw : t - w + w = t := by ring
  rw [hzw, he, hu, hV]
  field_simp
  ring

def addDefect_eq_nearNeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z ≡ −w の極が消える——u = tû と q = t³g で t² が約される。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★**`℘ − ℘(w)` は `w` でちょうど 1 位の零点**（`2w ∉ Λ`）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`℘′(w) ≠ 0`（第 625 の `derivWeierstrassP_eq_zero_iff` で `2w ∉ Λ` から）だから。

☆★これが `F_w` の `z ≡ w` での極が立たない理由である——
分母 `℘(z) − ℘(w)` が 1 位、分子 `℘′(z) − ℘′(w)` も `w` で消えるので比は有界。 -/
theorem analyticOrderAt_weierstrassP_sub_self (P : PeriodPair) (w : ℂ)
    (hw : w ∉ P.lattice) (h2w : 2 * w ∉ P.lattice) :
    analyticOrderAt (fun z : ℂ => P.weierstrassP z - P.weierstrassP w) w = 1 := by
  have hana : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP z - P.weierstrassP w) w :=
    (P.analyticOnNhd_weierstrassP w hw).sub analyticAt_const
  refine hana.analyticOrderAt_eq_one_of_zero_deriv_ne_zero (by simp) ?_
  have hd : deriv (fun z : ℂ => P.weierstrassP z - P.weierstrassP w) w
      = P.derivWeierstrassP w := ((hasDerivAt_weierstrassP P hw).sub_const _).deriv
  rw [hd]
  exact fun hcc => h2w ((derivWeierstrassP_eq_zero_iff P w hw).1 hcc)

def analyticOrderAt_weierstrassP_sub_self.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ − ℘(w) は w で 1 位の零点——2w ∉ Λ のとき。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★極を埋めた関数を作る道具 -/

open Filter Topology in
/-- ★★★★★★連続な点では `limUnder` は値そのもの。 -/
theorem limUnder_eq_self_of_continuousAt (f : ℂ → ℂ) (z : ℂ) (h : ContinuousAt f z) :
    limUnder (𝓝[≠] z) f = f z :=
  h.continuousWithinAt.tendsto.limUnder_eq

open Filter Topology in
/-- ★★★★★★★★★★★★★★★★★★★★**除去可能特異点を埋める一般の道具**。

`f` が `z` の除いた近傍で解析関数 `A` と一致し、その近傍の各点で連続なら、

    `Ext ≔ fun y => limUnder (𝓝[≠] y) f`

は `z` で解析的（そして `Ext = A` が `z` の近傍で成り立つ）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**これが `℘` の加法定理を Liouville に繋ぐのに要る道具である**——
mathlib の `℘` は格子点で junk value を取るので `addDefect` はそこで連続でなく、
`elliptic_liouville`（`Differentiable` を要求）に直接は当てられない。
☆`MeromorphicAt.analyticAt` は `ContinuousAt` を要求するので使えない。 -/
theorem analyticAt_limUnder_of_eventuallyEq (f A : ℂ → ℂ) (z : ℂ)
    (hA : AnalyticAt ℂ A z)
    (hpunct : ∀ᶠ y in 𝓝[≠] z, ContinuousAt f y)
    (heq : f =ᶠ[𝓝[≠] z] A) :
    AnalyticAt ℂ (fun y => limUnder (𝓝[≠] y) f) z := by
  refine hA.congr ?_
  have hz : limUnder (𝓝[≠] z) f = A z := by
    have h1 : Tendsto A (𝓝[≠] z) (𝓝 (A z)) := hA.continuousAt.continuousWithinAt.tendsto
    exact (h1.congr' heq.symm).limUnder_eq
  have hpt : ∀ᶠ y in 𝓝[≠] z, A y = limUnder (𝓝[≠] y) f := by
    filter_upwards [hpunct, heq] with y hy hyeq
    exact (hy.continuousWithinAt.tendsto.limUnder_eq.trans hyeq).symm
  rw [Filter.EventuallyEq, ← nhdsNE_sup_pure z, Filter.eventually_sup]
  exact ⟨hpt, by simpa using hz.symm⟩

open Filter Topology in
/-- ★★★★★★★★★★解析的な点では `Ext` は `f` と一致し、そこで解析的。 -/
theorem analyticAt_limUnder_of_analyticAt (f : ℂ → ℂ) (z : ℂ)
    (hana : ∀ᶠ y in 𝓝 z, AnalyticAt ℂ f y) :
    AnalyticAt ℂ (fun y => limUnder (𝓝[≠] y) f) z := by
  have hz : AnalyticAt ℂ f z := hana.self_of_nhds
  refine analyticAt_limUnder_of_eventuallyEq f f z hz ?_ (by rfl)
  filter_upwards [mem_nhdsWithin_of_mem_nhds hana] with y hy
  exact hy.continuousAt

def analyticAt_limUnder_of_eventuallyEq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(除去可能特異点を埋める一般の道具——limUnder で整関数を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★**`addDefect` は「良い点」で解析的**

    `z ∉ Λ`・`z + w ∉ Λ`・`℘(z) ≠ ℘(w)` なら `F_w` は `z` で解析的

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが `℘` の加法定理の組み立ての場合 (a) である。
☆残る 3 か所（`Λ`・`z ≡ w`・`z ≡ −w`）は第 610・629・627 で局所形が取れており、
それ以外の点は第 624-625 の単射性で存在しない。 -/
theorem analyticAt_addDefect (P : PeriodPair) (w : ℂ) {z : ℂ}
    (hz : z ∉ P.lattice) (hzw : z + w ∉ P.lattice)
    (hne : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    AnalyticAt ℂ (addDefect P w) z := by
  have hshift : AnalyticAt ℂ (fun s : ℂ => P.weierstrassP (s + w)) z := by
    have hf : AnalyticAt ℂ (fun s : ℂ => s + w) z := analyticAt_id.add analyticAt_const
    have hg : AnalyticAt ℂ P.weierstrassP ((fun s : ℂ => s + w) z) :=
      P.analyticOnNhd_weierstrassP (z + w) hzw
    exact AnalyticAt.comp (f := fun s : ℂ => s + w) (x := z) hg hf
  have hp : AnalyticAt ℂ P.weierstrassP z := P.analyticOnNhd_weierstrassP z hz
  have hpd : AnalyticAt ℂ P.derivWeierstrassP z := P.analyticOnNhd_derivWeierstrassP z hz
  have hratio : AnalyticAt ℂ
      (fun s : ℂ => (P.derivWeierstrassP s - P.derivWeierstrassP w)
        / (P.weierstrassP s - P.weierstrassP w)) z :=
    (hpd.sub analyticAt_const).div (hp.sub analyticAt_const) hne
  show AnalyticAt ℂ (fun s : ℂ => P.weierstrassP (s + w) + P.weierstrassP s
      + P.weierstrassP w
      - ((P.derivWeierstrassP s - P.derivWeierstrassP w)
          / (P.weierstrassP s - P.weierstrassP w)) ^ 2 / 4) z
  exact ((hshift.add hp).add analyticAt_const).sub
    ((hratio.pow 2).div analyticAt_const (by norm_num))

open Filter Topology in
/-- ★★★★★★「良い点」の集合は開——`℘` は格子の外で連続だから。 -/
theorem isOpen_goodSet (P : PeriodPair) (w : ℂ) :
    IsOpen {z : ℂ | z ∉ P.lattice ∧ z + w ∉ P.lattice
      ∧ P.weierstrassP z - P.weierstrassP w ≠ 0} := by
  rw [isOpen_iff_mem_nhds]
  intro z hz
  have h1 : {y : ℂ | y ∉ P.lattice} ∈ nhds z :=
    (P.isClosed_lattice.isOpen_compl).mem_nhds hz.1
  have h2 : {y : ℂ | y + w ∉ P.lattice} ∈ nhds z := by
    have hopen : IsOpen {y : ℂ | y + w ∉ P.lattice} := by
      have he : {y : ℂ | y + w ∉ P.lattice}
          = (fun y : ℂ => y + w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
      rw [he]
      exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
    exact hopen.mem_nhds hz.2.1
  have hcont : ContinuousAt (fun y : ℂ => P.weierstrassP y - P.weierstrassP w) z :=
    ((P.analyticOnNhd_weierstrassP z hz.1).sub analyticAt_const).continuousAt
  have h3 : ∀ᶠ y in nhds z, P.weierstrassP y - P.weierstrassP w ≠ 0 :=
    hcont.eventually_ne hz.2.2
  filter_upwards [h1, h2, h3] with y hy1 hy2 hy3
  exact ⟨hy1, hy2, hy3⟩

def analyticAt_addDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(F_w は良い点で解析的——加法定理の組み立ての場合 (a)。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
