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

/-!
# 一様化 —— 一様化の座標・周期性と対称性・`ω`-正規化・同種写像の定義的性質

☆`Found/GenEll/Uniformization.lean`(292 KB / 325 宣言)を
**ファイル内の見出し**で割った 1 枚である(2026-09-03、第 1456)。
★論文のセクションでは割れない——このファイルは [GenEll] §3 の
`Lemma 3.5` と `Proposition 3.4` の 2 項目しか持たず、割っても 146 KB のままだからである。
☆読む順序は `Basic → VeluAnalytic → Surjective → AdditionEntry → AdditionODE
→ FilledPole → AdditionFormula → Phi → GroupIso → Sublattice → G2G3 → Assemble`。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★一様化の座標 -/

/-- ★★★★★一様化の `x` 座標 `℘(z)`。 -/
noncomputable def latticePointX (P : PeriodPair) (z : ℂ) : ℂ := P.weierstrassP z

/-- ★★★★★一様化の `y` 座標 `℘′(z)/2`。

★`latticeCurve P` は `y² = x³ − (g₂/4)x − (g₃/4)` なので、
`℘′² = 4℘³ − g₂℘ − g₃` を `4` で割った形に合わせるために `/2` が要る。 -/
noncomputable def latticePointY (P : PeriodPair) (z : ℂ) : ℂ := P.derivWeierstrassP z / 2

/-- ★★★★★★★★★★★★★★★★**`(℘(z), ℘′(z)/2)` は `latticeCurve P` の上にある**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★mathlib の `derivWeierstrassP_sq`（`℘′² = 4℘³ − g₂℘ − g₃`）を `4` で割るだけである。 -/
theorem latticeCurve_equation (P : PeriodPair) (z : ℂ) (hz : z ∉ P.lattice) :
    (latticeCurve P).toAffine.Equation (latticePointX P z) (latticePointY P z) := by
  have h := P.derivWeierstrassP_sq z hz
  rw [WeierstrassCurve.Affine.equation_iff]
  simp only [latticeCurve, latticePointX, latticePointY]
  linear_combination h / 4

/-! ## ★★★★★周期性と対称性 -/

/-- ★★★★★`x` 座標は周期的。 -/
theorem latticePointX_add_coe (P : PeriodPair) (z : ℂ) (l : P.lattice) :
    latticePointX P (z + (l : ℂ)) = latticePointX P z :=
  P.weierstrassP_add_coe z l

/-- ★★★★★`y` 座標は周期的。 -/
theorem latticePointY_add_coe (P : PeriodPair) (z : ℂ) (l : P.lattice) :
    latticePointY P (z + (l : ℂ)) = latticePointY P z := by
  simp only [latticePointY, P.derivWeierstrassP_add_coe z l]

/-- ★★★★`x` 座標は偶——`℘(−z) = ℘(z)`。 -/
theorem latticePointX_neg (P : PeriodPair) (z : ℂ) :
    latticePointX P (-z) = latticePointX P z :=
  P.weierstrassP_neg z

/-- ★★★★`y` 座標は奇——`℘′(−z) = −℘′(z)`。

★`latticeCurve P` は `a₁ = a₃ = 0` なので `negY x y = −y`。
すなわち `−z` は `z` の点の `−` である。 -/
theorem latticePointY_neg (P : PeriodPair) (z : ℂ) :
    latticePointY P (-z) = -latticePointY P z := by
  simp only [latticePointY, P.derivWeierstrassP_neg z]
  ring

/-! ## ★★★★★★★★★★★★★★★★★★★★`ω`-正規化は自動である -/

/-- ★★★★★★★★★★★★★★★★★★★★**一様化は自動で `ω`-正規化されている**。

    `d(℘(z))/dz = 2·(℘′(z)/2) = 2y`

★`latticeCurve P` の不変微分は `dx/(2y + a₁x + a₃) = dx/(2y)`（`a₁ = a₃ = 0`）なので、
上の等式は **`dx/(2y) = dz`** と同じことである。

★★★これが `Found/GaloisRep/VeluNormalized.lean` の `archDefect_isogeny`（第 596）で
`α_σ` が現れる場所であり、Vélu の正規化（`Found/GenEll/Velu.lean` の
`velu_omega_gen`、第 591）が `α_σ = 1` を与える理由である。 -/
theorem deriv_latticePointX (P : PeriodPair) :
    deriv (latticePointX P) = fun z => 2 * latticePointY P z := by
  have hfun : latticePointX P = P.weierstrassP := rfl
  rw [hfun, P.deriv_weierstrassP]
  funext z
  simp only [latticePointY]
  ring

/-! ## ★★★★★★★★★★★★同種写像の定義的性質 -/

/-- ★★★★★★★★★★★★**`Λ ⊆ Λ′` なら `℘_{Λ′}` は `Λ`-周期的**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★これが「`z ↦ z` が `ℂ/Λ → ℂ/Λ′` を誘導する」ことの中身であり、
同種写像 `E_Λ → E_{Λ′}` が存在する理由である。

☆残るのは「`℘_{Λ′}` は `℘_Λ` の有理式で書ける」こと——それが Vélu の公式の解析側で
あり、`Found/GenEll/Velu.lean` が代数側で与えている式そのものである。 -/
theorem weierstrassP_periodic_of_le (P P' : PeriodPair) (h : P.lattice ≤ P'.lattice)
    (z : ℂ) (l : P.lattice) :
    P'.weierstrassP (z + (l : ℂ)) = P'.weierstrassP z :=
  P'.weierstrassP_add_coe z ⟨(l : ℂ), h l.2⟩

/-- ★★★★★★★★★★★★**`Λ ⊆ Λ′` なら `℘′_{Λ′}` も `Λ`-周期的**。 -/
theorem derivWeierstrassP_periodic_of_le (P P' : PeriodPair) (h : P.lattice ≤ P'.lattice)
    (z : ℂ) (l : P.lattice) :
    P'.derivWeierstrassP (z + (l : ℂ)) = P'.derivWeierstrassP z :=
  P'.derivWeierstrassP_add_coe z ⟨(l : ℂ), h l.2⟩

/-- ★★★★★★★★★★★★★★**`Λ ⊆ Λ′` なら `z ↦ (℘_{Λ′}(z), ℘′_{Λ′}(z)/2)` は
`Λ`-周期的で `latticeCurve P′` の上にある**——★同種写像そのもの。 -/
theorem latticePoint_isogeny (P P' : PeriodPair) (h : P.lattice ≤ P'.lattice)
    (z : ℂ) (hz : z ∉ P'.lattice) (l : P.lattice) :
    (latticeCurve P').toAffine.Equation (latticePointX P' z) (latticePointY P' z)
      ∧ latticePointX P' (z + (l : ℂ)) = latticePointX P' z
      ∧ latticePointY P' (z + (l : ℂ)) = latticePointY P' z :=
  ⟨latticeCurve_equation P' z hz,
   weierstrassP_periodic_of_le P P' h z l,
   by simp only [latticePointY, derivWeierstrassP_periodic_of_le P P' h z l]⟩

/-! ## ★出典の紐付け（`.src`）——★★条つき（一様化の全射性は含まない） -/

def latticeCurve_equation.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4((℘(z), ℘′(z)/2) は latticeCurve の上にある。★無条件)",
    sectionId := "genell-prop-3-4" }

def deriv_latticePointX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(一様化は自動で ω-正規化されている——dx/(2y) = dz。★無条件)",
    sectionId := "genell-prop-3-4" }

def latticePoint_isogeny.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Λ ⊆ Λ′ なら z ↦ z が同種写像を誘導する——解析側。★無条件)",
    sectionId := "genell-prop-3-4" }

def latticePoint_isogeny.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆一様化写像 ℂ/Λ → E(ℂ) が全単射であること(全射性)は mathlib に無い" ++
       "(2026-08-29 に測定: Analysis/SpecialFunctions/Elliptic/Weierstrass.lean は" ++
       "℘・℘′・g₂・g₃・周期性・偶奇・℘′² = 4℘³ − g₂℘ − g₃ まで)") 9,
    .implicitStep
      ("☆℘_{Λ′} を ℘_Λ の有理式で書くこと(Vélu の公式の解析側)。" ++
       "★代数側は Found/GenEll/Velu.lean が持っている(第 586-593)。" ++
       "★★繋ぐには「極を持たない楕円関数は定数」(Liouville)が要る") 9,
    .implicitStep
      ("★★★到達点(2026-08-29、第 597): mathlib が ℘ を丸ごと持っていることが分かり、" ++
       "解析側の第一の煉瓦が置けた。★一様化が自動で ω-正規化されている" ++
       "(deriv_latticePointX)ので、第 596 の α_σ が 1 になる仕組みが Lean 上に載った") 8 ]

end ABC3.Found.GenEll
