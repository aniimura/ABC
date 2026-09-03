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

/-!
# 一様化 —— 平行移動が誘導する置換・`℘` は全射・2-捩れ点

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★平行移動が誘導する置換 -/

open Classical in
/-- ★★★★★平行移動 `w ↦ w + w₀` が代表系 `T` に誘導する写像。

★`-(w + w₀)` の代表を選ぶ（`Λ` を法として `w + w₀` と合同な `T` の元）。 -/
noncomputable def shiftRep (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (w : ℂ) : ℂ :=
  if h : ∃ v ∈ T, -(w + w₀) + v ∈ P.lattice then h.choose else 0

/-- ★★★★★★**代表系の一意性**——`hrep` の言い換え。 -/
theorem rep_unique (P P' : PeriodPair) (T : Finset ℂ)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (p : ℂ) (hp : p ∈ P'.lattice) (v v' : ℂ) (hv : v ∈ T) (hv' : v' ∈ T)
    (h1 : p + v ∈ P.lattice) (h2 : p + v' ∈ P.lattice) : v = v' := by
  obtain ⟨u, hu, -, huniq⟩ := hrep p hp
  have e1 : v = u := by by_contra hc; exact huniq v hv hc h1
  have e2 : v' = u := by by_contra hc; exact huniq v' hv' hc h2
  rw [e1, e2]

/-- ★★★★★★★★★★★★★★★★**平行移動は代表系を置換する**——`veluAnalyticX_shift` の
仮定 `σ` が**代表系の定義から作れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★単射性は代表系の一意性から、全射性は `T` が有限だから出る。 -/
theorem exists_veluShiftPerm (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (w₀ : ℂ) (hw₀ : w₀ ∈ P'.lattice) :
    ∃ σ : ℂ → ℂ, (∀ w ∈ T, σ w ∈ T) ∧ (∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
      ∧ (∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v) ∧ (∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) := by
  classical
  have hex : ∀ w ∈ T, ∃ v ∈ T, -(w + w₀) + v ∈ P.lattice := by
    intro w hw
    have hmem : -(w + w₀) ∈ P'.lattice := neg_mem (P'.lattice.add_mem (hT w hw) hw₀)
    obtain ⟨v, hv, hv2, -⟩ := hrep _ hmem
    exact ⟨v, hv, hv2⟩
  have hmemT : ∀ w ∈ T, shiftRep P T w₀ w ∈ T := by
    intro w hw
    rw [shiftRep, dif_pos (hex w hw)]
    exact (hex w hw).choose_spec.1
  have hinj : ∀ w ∈ T, ∀ w' ∈ T, shiftRep P T w₀ w = shiftRep P T w₀ w' → w = w' := by
    intro w hw w' hw' h
    have s1 := (hex w hw).choose_spec.2
    have s2 := (hex w' hw').choose_spec.2
    rw [shiftRep, dif_pos (hex w hw)] at h
    rw [shiftRep, dif_pos (hex w' hw')] at h
    rw [h] at s1
    have hdiff : -w + w' ∈ P.lattice := by
      have hs := P.lattice.sub_mem s2 s1
      have he : (-(w' + w₀) + (hex w' hw').choose) - (-(w + w₀) + (hex w' hw').choose)
          = -w' + w := by ring
      rw [he] at hs
      have hn := neg_mem hs
      simpa using hn
    have hzero : -w + w ∈ P.lattice := by simpa using P.lattice.zero_mem
    exact rep_unique P P' T hrep (-w) (neg_mem (hT w hw)) w w' hw hw' hzero hdiff
  refine ⟨shiftRep P T w₀, hmemT, hinj, ?_, ?_⟩
  · intro v hv
    obtain ⟨a, ha, hav⟩ := Finset.surj_on_of_inj_on_of_card_le (s := T) (t := T)
      (fun a _ => shiftRep P T w₀ a) (fun a ha => hmemT a ha)
      (fun a₁ a₂ ha₁ ha₂ h => hinj a₁ ha₁ a₂ ha₂ h) le_rfl v hv
    exact ⟨a, ha, hav.symm⟩
  · intro w hw
    have s1 := (hex w hw).choose_spec.2
    rw [shiftRep, dif_pos (hex w hw)]
    have he : w + w₀ - (hex w hw).choose = -(-(w + w₀) + (hex w hw).choose) := by ring
    rw [he]
    exact neg_mem s1

/-- ★★★★★★★★★★★★★★★★★★★★**`veluAnalyticX` は `Λ′`-周期的**——★**無条件**
（代表系であることだけから）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 599 では `σ` を仮定として受けていたが、第 602 でそれが**代表系の定義から作れた**。 -/
theorem veluAnalyticX_periodic (P P' : PeriodPair) (T : Finset ℂ) (c : ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (z : ℂ) (l : ℂ) (hl : l ∈ P'.lattice) :
    veluAnalyticX P T c (z + l) = veluAnalyticX P T c z := by
  obtain ⟨σ, hmem, hinj, hsurj, hshift⟩ := exists_veluShiftPerm P P' T hT hrep l hl
  exact veluAnalyticX_shift P T c l σ hmem hinj hsurj hshift z

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側
——代表系であることだけから**。

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**仮定は `T` が `Λ′/Λ` の代表系であること、それだけである**:

| 仮定 | 内容 |
|---|---|
| `hle` | `Λ ⊆ Λ′` |
| `h0T` | `0 ∈ T` |
| `hT` | `T ⊆ Λ′` |
| `hrep` | `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` |

★★★★★周期性・極の打ち消し・正規化はすべて塞がった
（第 598 Liouville、第 599 周期性、第 600 極、第 601 正規化、第 602 置換の構成）。

☆残るのは、この `℘` の等式を `Found/GenEll/Velu.lean` の**代数側 `veluXGen`**
（第 591）へ翻訳することだけである——`℘` の加法定理が要る（mathlib に無い）。 -/
theorem weierstrassP_eq_velu_of_rep (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
  weierstrassP_eq_velu P P' hle T h0T hT hrep
    (fun y l hl => veluAnalyticX_periodic P P' T _ hT hrep y l hl) z

def exists_veluShiftPerm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(平行移動は代表系を置換する——σ は代表系の定義から作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu_of_rep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——仮定は T が Λ′/Λ の代表系であることだけ)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu_of_rep.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆この ℘ の等式を Found/GenEll/Velu.lean の代数側 veluXGen(第 591)へ" ++
       "翻訳すること。★℘ の加法定理が要る(mathlib に無い、2026-08-29 に測定)") 8,
    .implicitStep
      ("★★★★到達点(2026-08-29、第 602): Vélu の公式の解析側が" ++
       "「T が Λ′/Λ の代表系である」だけから従う形になった。" ++
       "第 598 Liouville・第 599 周期性・第 600 極・第 601 正規化・第 602 置換の構成") 9 ]

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

/-! ## ★★★★★★★★★★★★★★★★2-捩れ点 -/

/-- ★★★★★★★★★★★★**2-捩れ点では `℘′` が消える**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★証明は 3 行: `℘′` は奇であり、`2z ∈ Λ` なら `−z ≡ z (mod Λ)` なので
`℘′(z) = ℘′(−z) = −℘′(z)`。

★★`Found/GenEll/Velu.lean` の `veluV2`（2-捩れの場合分け）と
`velu2_omega` の仮定 `2y₀ + a₁x₀ + a₃ = 0` は、`latticeCurve` では
`a₁ = a₃ = 0` なので `2·(℘′(z)/2) = ℘′(z) = 0`——**本定理そのもの**である。 -/
theorem derivWeierstrassP_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ)
    (h2 : 2 * z ∈ P.lattice) : P.derivWeierstrassP z = 0 := by
  have hneg : P.derivWeierstrassP (-z) = -P.derivWeierstrassP z := P.derivWeierstrassP_neg z
  have hper : P.derivWeierstrassP (z + (-(2 * z))) = P.derivWeierstrassP z :=
    P.derivWeierstrassP_add_coe z ⟨-(2 * z), neg_mem h2⟩
  have hz : z + (-(2 * z)) = -z := by ring
  rw [hz, hneg] at hper
  linear_combination -hper / 2

/-- ★★★★★★★★★★**2-捩れ点の `y` 座標は `0`**——`Velu.lean` の 2-捩れの場合分けの中身。 -/
theorem latticePointY_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ)
    (h2 : 2 * z ∈ P.lattice) : latticePointY P z = 0 := by
  simp [latticePointY, derivWeierstrassP_eq_zero_of_two_mem P z h2]

/-- ★★★★★★★★★★★★**2-捩れ点の `x` 座標は `4x³ − g₂x − g₃` の根**。

★`℘′² = 4℘³ − g₂℘ − g₃` の左辺が `0` になるから。
☆`latticeCurve P = ⟨0,0,0,−g₂/4,−g₃/4⟩` の 2-捩れ点はちょうどここである。 -/
theorem cubic_eq_zero_of_two_mem (P : PeriodPair) (z : ℂ) (hz : z ∉ P.lattice)
    (h2 : 2 * z ∈ P.lattice) :
    4 * (latticePointX P z) ^ 3 - P.g₂ * (latticePointX P z) - P.g₃ = 0 := by
  have hsq := P.derivWeierstrassP_sq z hz
  rw [derivWeierstrassP_eq_zero_of_two_mem P z h2] at hsq
  simp only [latticePointX]
  linear_combination -hsq

def derivWeierstrassP_eq_zero_of_two_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2-捩れ点では ℘′ が消える——℘′ は奇だから。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
