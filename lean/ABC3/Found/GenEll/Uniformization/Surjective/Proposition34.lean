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
import ABC3.Found.GenEll.Uniformization.VeluAnalytic
import ABC3.Found.GenEll.Uniformization.Surjective.Lemma35

/-!
# Surjective —— `[GenEll] Proposition 3.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`℘` は全射 -/

/-- ★★★★★`ω₁/2` は格子に入らない（`ω₁`・`ω₂` の ℝ-一次独立から）。 -/
theorem half_omega1_notMem (P : PeriodPair) : (P.ω₁ / 2) ∉ P.lattice := by
  intro h
  rw [PeriodPair.lattice, Submodule.mem_span_pair] at h
  obtain ⟨m, n, hmn⟩ := h
  have hind := LinearIndependent.pair_iff.1 P.indep ((m : ℝ) - 1/2) (n : ℝ) ?_
  · have h1 := hind.1
    have h2 : (2 * m : ℤ) = 1 := by exact_mod_cast (by linarith : (2 * (m:ℝ)) = 1)
    omega
  · have hc : (m : ℂ) • P.ω₁ + (n : ℂ) • P.ω₂ = P.ω₁ / 2 := by
      simpa [zsmul_eq_mul] using hmn
    have h2 : ((m : ℝ) - 1/2) • P.ω₁ + ((n : ℝ)) • P.ω₂
        = ((m : ℂ) * P.ω₁ + (n : ℂ) * P.ω₂) - P.ω₁ / 2 := by
      push_cast [Complex.real_smul]
      ring
    rw [h2, ← hc]
    ring

open Classical Filter Topology Bornology in
/-- ★★★★★★格子の外では `(℘ − x₀)⁻¹` は微分可能。 -/
theorem wp_inv_differentiableAt_of_notMem (P : PeriodPair) (x₀ : ℂ)
    (hcon : ∀ z ∉ P.lattice, P.weierstrassP z ≠ x₀) (p : ℂ) (hp : p ∉ P.lattice) :
    DifferentiableAt ℂ (fun z => if z ∈ P.lattice then (0:ℂ)
      else (P.weierstrassP z - x₀)⁻¹) p := by
  have hopen : IsOpen ((P.lattice : Set ℂ)ᶜ) := P.isClosed_lattice.isOpen_compl
  have hA : AnalyticAt ℂ (fun z => (P.weierstrassP z - x₀)⁻¹) p :=
    ((P.analyticOnNhd_weierstrassP p hp).sub analyticAt_const).inv (sub_ne_zero.2 (hcon p hp))
  refine hA.differentiableAt.congr_of_eventuallyEq ?_
  filter_upwards [hopen.mem_nhds hp] with z hz
  simp only [Set.mem_compl_iff, SetLike.mem_coe] at hz
  simp [hz]

open Classical Filter Topology Bornology in
/-- ★★★★★★★★★★格子点でも `(℘ − x₀)⁻¹` は微分可能——★**除去可能特異点**。

★`℘ → ∞` なので `(℘ − x₀)⁻¹ → 0`。値を `0` に決めれば連続になり、
Riemann の除去可能特異点定理（mathlib の
`Complex.differentiableOn_compl_singleton_and_continuousAt_iff`）で微分可能になる。 -/
theorem wp_inv_differentiableAt_of_mem (P : PeriodPair) (x₀ : ℂ)
    (hcon : ∀ z ∉ P.lattice, P.weierstrassP z ≠ x₀) (p : ℂ) (hp : p ∈ P.lattice) :
    DifferentiableAt ℂ (fun z => if z ∈ P.lattice then (0:ℂ)
      else (P.weierstrassP z - x₀)⁻¹) p := by
  set g : ℂ → ℂ := fun z => if z ∈ P.lattice then (0:ℂ) else (P.weierstrassP z - x₀)⁻¹ with hg
  set s : Set ℂ := ((P.lattice : Set ℂ) \ {p})ᶜ with hs
  have hsnhds : s ∈ 𝓝 p := P.isOpen_compl_lattice_sdiff.mem_nhds (by simp)
  have hoff : ∀ z ∈ s \ {p}, z ∉ P.lattice := by
    rintro z ⟨hz1, hz2⟩
    intro hc
    exact hz1 ⟨hc, by simpa using hz2⟩
  have hcont : ContinuousAt g p := by
    rw [← continuousWithinAt_compl_self]
    have hord : meromorphicOrderAt P.weierstrassP p < 0 := by
      rw [P.order_weierstrassP p hp]; decide
    have h1 : Tendsto P.weierstrassP (𝓝[≠] p) (cobounded ℂ) :=
      tendsto_cobounded_of_meromorphicOrderAt_neg hord
    have hsub : Tendsto (fun w : ℂ => w - x₀) (cobounded ℂ) (cobounded ℂ) := by
      simpa using (tendsto_sub_cobounded_right (α := ℂ) x₀)
    have h3 : Tendsto (fun z => (P.weierstrassP z - x₀)⁻¹) (𝓝[≠] p) (𝓝 0) :=
      tendsto_inv₀_cobounded.comp (hsub.comp h1)
    have hgp : g p = 0 := by simp [hg, hp]
    rw [ContinuousWithinAt, hgp]
    refine h3.congr' ?_
    filter_upwards [self_mem_nhdsWithin, mem_nhdsWithin_of_mem_nhds hsnhds] with z hz1 hz2
    have hz : z ∉ P.lattice := hoff z ⟨hz2, by simpa using hz1⟩
    simp [hg, hz]
  have hdon : DifferentiableOn ℂ g s := by
    rw [← Complex.differentiableOn_compl_singleton_and_continuousAt_iff hsnhds]
    exact ⟨fun z hz =>
      (wp_inv_differentiableAt_of_notMem P x₀ hcon z (hoff z hz)).differentiableWithinAt, hcont⟩
  exact hdon.differentiableAt hsnhds

open Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`℘` は全射**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★機構は第 598 の Liouville そのもの: もし `℘` が `x₀` を取らないなら

    `g(z) ≔ 1/(℘(z) − x₀)`（格子点では `0`）

は**整で二重周期的**なので定数。`g(0) = 0` だから `g ≡ 0`、
すなわち `℘(ω₁/2) = x₀`——仮定に反する。

★★★★☆**これが一様化 `ℂ/Λ → E(ℂ)` の全射性の `x` 座標の段である**
（`§9-1039`（第 597）で「mathlib に無い」と測ったもの）。
☆残るのは `y` 座標（`℘′` の符号の選択）と単射性である。 -/
theorem weierstrassP_surjective (P : PeriodPair) (x₀ : ℂ) :
    ∃ z, z ∉ P.lattice ∧ P.weierstrassP z = x₀ := by
  by_contra hcon0
  push_neg at hcon0
  set g : ℂ → ℂ := fun z => if z ∈ P.lattice then (0:ℂ) else (P.weierstrassP z - x₀)⁻¹ with hg
  have hper : ∀ z : ℂ, ∀ l ∈ P.lattice, g (z + l) = g z := by
    intro z l hl
    by_cases hz : z ∈ P.lattice
    · have hzl : z + l ∈ P.lattice := P.lattice.add_mem hz hl
      simp [hg, hz, hzl]
    · have hzl : z + l ∉ P.lattice := fun hc => hz (by simpa using P.lattice.sub_mem hc hl)
      simp only [hg, if_neg hz, if_neg hzl, P.weierstrassP_add_coe z ⟨l, hl⟩]
  have hdiff : Differentiable ℂ g := by
    intro p
    by_cases hp : p ∈ P.lattice
    · exact wp_inv_differentiableAt_of_mem P x₀ hcon0 p hp
    · exact wp_inv_differentiableAt_of_notMem P x₀ hcon0 p hp
  have hhalf : P.ω₁ / 2 ∉ P.lattice := half_omega1_notMem P
  have hconst := elliptic_liouville P g hdiff hper (P.ω₁ / 2) 0
  have h0 : g 0 = 0 := by simp [hg, P.lattice.zero_mem]
  have hval : (P.weierstrassP (P.ω₁ / 2) - x₀)⁻¹ = 0 := by
    have hgz : g (P.ω₁ / 2) = 0 := by rw [hconst, h0]
    simpa [hg, hhalf] using hgz
  exact hcon0 (P.ω₁ / 2) hhalf (by
    have hz := inv_eq_zero.1 hval
    linear_combination hz)

def weierstrassP_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(℘ は全射——一様化の全射性の x 座標の段。★無条件)",
    sectionId := "genell-prop-3-4" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**一様化は全射**。

    `latticeCurve P` の上の任意の点 `(x₀, y₀)` は `(℘(z), ℘′(z)/2)` の形である

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★機構は 2 段:

1. `℘` は全射（第 603）なので `℘(z₀) = x₀` となる `z₀ ∉ Λ` がある
2. `℘′(z₀)² = 4x₀³ − g₂x₀ − g₃ = (2y₀)²` なので `℘′(z₀) = ±2y₀`。
   ★符号が合わなければ `−z₀` を取る（`℘` は偶・`℘′` は奇）

★★★★☆**これが `§9-1039`（第 597）で「mathlib に無い」と測った
一様化 `ℂ/Λ → E(ℂ)` の全射性である**——第 603 と合わせて塞がった。

☆残るのは単射性（`℘(z) = ℘(w)` かつ `℘′(z) = ℘′(w)` なら `z ≡ w mod Λ`）。 -/
theorem latticePoint_surjective (P : PeriodPair) (x₀ y₀ : ℂ)
    (h : (latticeCurve P).toAffine.Equation x₀ y₀) :
    ∃ z, z ∉ P.lattice ∧ latticePointX P z = x₀ ∧ latticePointY P z = y₀ := by
  obtain ⟨z₀, hz₀, hx⟩ := weierstrassP_surjective P x₀
  have hsq := P.derivWeierstrassP_sq z₀ hz₀
  rw [WeierstrassCurve.Affine.equation_iff] at h
  simp only [latticeCurve] at h
  have hy : (P.derivWeierstrassP z₀) ^ 2 = (2 * y₀) ^ 2 := by
    rw [hsq, hx]
    linear_combination -4 * h
  rcases sq_eq_sq_iff_eq_or_eq_neg.1 hy with hcase | hcase
  · exact ⟨z₀, hz₀, hx, by simp only [latticePointY, hcase]; ring⟩
  · refine ⟨-z₀, fun hc => hz₀ (by simpa using neg_mem hc), ?_, ?_⟩
    · simp only [latticePointX, P.weierstrassP_neg z₀]; exact hx
    · simp only [latticePointY, P.derivWeierstrassP_neg z₀, hcase]; ring

def latticePoint_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(一様化は全射——ℂ/Λ → E(ℂ) の全射性。★無条件)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
