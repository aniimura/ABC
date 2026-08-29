/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★一様化 `z ↦ (℘(z), ℘′(z)/2)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★これは何か

`Found/GaloisRep/VeluNormalized.lean` の `htFalt_isogeny_le_of_analytic`（第 596）に
残った入力は「**次数 `l` の解析的同種写像のデータ**」であった:
`α_σ`・`u′_σ = α_σ·u_σ`・`α_σ·Λ_σ ⊆ Λ′_σ` の指数が `l`。

★本ファイルはその**解析側の第一の煉瓦**である。

## ★★★★★★mathlib の在庫（2026-08-29 に測定）

★★**mathlib は Weierstrass `℘` の理論を丸ごと持っている**
（`Mathlib/Analysis/SpecialFunctions/Elliptic/Weierstrass.lean`、1080 行）:

| | |
|---|---|
| `PeriodPair` | ★本プロジェクトが使っている構造そのもの |
| `weierstrassP` `℘[L]`・`derivWeierstrassP` `℘′[L]` | 定義・局所一様収束・微分可能性 |
| `weierstrassP_add_coe`・`derivWeierstrassP_add_coe` | ★**周期性** |
| `weierstrassP_neg`・`derivWeierstrassP_neg` | ★偶・奇 |
| `deriv_weierstrassP` | `deriv ℘ = ℘′` |
| `g₂`・`g₃`・`G n` | Eisenstein 級数 |
| `derivWeierstrassP_sq` | ★★★**`℘′(z)² = 4℘(z)³ − g₂℘(z) − g₃`** |

☆**無いもの**: 一様化写像が全単射であること（`ℂ/Λ ≅ E(ℂ)`）、
群法則との両立、`℘_{Λ′}` を `℘_Λ` の有理式で書くこと（＝Vélu の公式の解析側）。

## ★★★本ファイルが取るもの

| 定理 | 内容 |
|---|---|
| `latticeCurve_equation` | ★★★★`(℘(z), ℘′(z)/2)` は `latticeCurve P` の上にある |
| `latticePointX/Y_add_coe` | 周期性 |
| `latticePointX/Y_neg` | `−z ↦ −(点)` |
| `deriv_latticePointX` | ★★★★★**一様化は自動で `ω`-正規化されている** |
| `weierstrassP_periodic_of_le` | ★★★★★★`Λ ⊆ Λ′` なら `℘_{Λ′}` は `Λ`-周期的（同種写像の定義的性質） |

★★`latticeCurve P = ⟨0, 0, 0, −g₂/4, −g₃/4⟩` なので不変微分は `dx/(2y)`。
`x = ℘`、`y = ℘′/2` を入れると `dx = ℘′dz`、`2y = ℘′` で **`dx/(2y) = dz`**
——☆これが `§9-1038`（第 596）の `α_σ` が `1` になる仕組みである。
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

/-! ## ★★★★★★★★★★★★★★★★★★★★★★楕円関数の Liouville -/

/-- ★★★★★`{ω₁, ω₂}` を ℝ-基底として読んだもの（`indep` と `dim_ℝ ℂ = 2` から）。 -/
noncomputable def realBasis (P : PeriodPair) : Module.Basis (Fin 2) ℝ ℂ :=
  basisOfLinearIndependentOfCardEqFinrank P.indep (by simp)

theorem coe_realBasis (P : PeriodPair) : ⇑(realBasis P) = ![P.ω₁, P.ω₂] :=
  coe_basisOfLinearIndependentOfCardEqFinrank _ _

/-- ★★★★座標で書き直す。 -/
theorem realBasis_repr_sum (P : PeriodPair) (z : ℂ) :
    ((realBasis P).repr z 0 : ℝ) • P.ω₁ + ((realBasis P).repr z 1 : ℝ) • P.ω₂ = z := by
  have h := (realBasis P).sum_repr z
  rw [Fin.sum_univ_two, coe_realBasis] at h
  simpa using h

/-- ★★★★★★**どの `z` も格子を引けば閉平行四辺形に入る**——`Int.fract` を取るだけ。 -/
theorem exists_mem_box (P : PeriodPair) (z : ℂ) :
    ∃ l ∈ P.lattice, ∃ c : Fin 2 → ℝ, c ∈ Set.Icc (0 : Fin 2 → ℝ) 1 ∧
      z - l = c 0 • P.ω₁ + c 1 • P.ω₂ := by
  set a := (realBasis P).repr z 0 with ha
  set b := (realBasis P).repr z 1 with hb
  refine ⟨(⌊a⌋ : ℤ) • P.ω₁ + (⌊b⌋ : ℤ) • P.ω₂, ?_, ![Int.fract a, Int.fract b], ?_, ?_⟩
  · exact Submodule.add_mem _
      (Submodule.smul_mem _ _ (Submodule.subset_span (by simp)))
      (Submodule.smul_mem _ _ (Submodule.subset_span (by simp)))
  · rw [Set.mem_Icc]
    refine ⟨fun i => ?_, fun i => ?_⟩ <;> fin_cases i <;>
      simp [Int.fract_nonneg, (Int.fract_lt_one _).le]
  · have h := realBasis_repr_sum P z
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Int.fract]
    rw [← h]
    module

/-- ★★★★閉平行四辺形をパラメータ付ける写像。 -/
noncomputable def boxMap (P : PeriodPair) : (Fin 2 → ℝ) → ℂ :=
  fun c => c 0 • P.ω₁ + c 1 • P.ω₂

theorem continuous_boxMap (P : PeriodPair) : Continuous (boxMap P) := by
  unfold boxMap; fun_prop

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**楕円関数の Liouville**
——**整で二重周期的な関数は定数**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★機構は 3 行:

1. どの `z` も格子を引けば閉平行四辺形 `K` に入る（`exists_mem_box`——`Int.fract`）
2. `K` はコンパクト（`Set.Icc 0 1` の連続像）なので `f` は `K` 上で有界、
   周期性から `Set.range f ⊆ f '' K` で**全平面で有界**
3. mathlib の Liouville（`Differentiable.apply_eq_apply_of_bounded`）

★★★☆**これが mathlib に無く、Vélu の公式の解析側に要る道具である**
（2026-08-29 に測定: `Analysis/SpecialFunctions/Elliptic/Weierstrass.lean` は
`℘` の理論を持つが、楕円関数の Liouville は無い）。

★★`℘_{Λ′}(z) = ℘_Λ(z) + Σ_w [℘_Λ(z+w) − ℘_Λ(w)]`（Vélu の `X` の解析側）は、
両辺の差が整で `Λ′`-周期的であること＋原点で `0` になることから従う。
☆本定理はその「整で `Λ′`-周期的なら定数」の段である。 -/
theorem elliptic_liouville (P : PeriodPair) (f : ℂ → ℂ) (hf : Differentiable ℂ f)
    (hper : ∀ (z : ℂ), ∀ l ∈ P.lattice, f (z + l) = f z) (z w : ℂ) : f z = f w := by
  refine hf.apply_eq_apply_of_bounded ?_ z w
  have hcpt : IsCompact (boxMap P '' Set.Icc 0 1) :=
    isCompact_Icc.image (continuous_boxMap P)
  have hsub : Set.range f ⊆ f '' (boxMap P '' Set.Icc 0 1) := by
    rintro _ ⟨x, rfl⟩
    obtain ⟨l, hl, c, hc, hx⟩ := exists_mem_box P x
    refine ⟨boxMap P c, ⟨c, hc, rfl⟩, ?_⟩
    rw [boxMap, ← hx]
    have h2 := hper (x - l) l hl
    rw [sub_add_cancel] at h2
    exact h2.symm
  exact ((hcpt.image hf.continuous).isBounded).subset hsub

/-- ★★★★★★★★★★★★★★★★★★**整で二重周期的なら `f` は `f 0` に等しい**（使いやすい形）。 -/
theorem elliptic_liouville_eq_zero (P : PeriodPair) (f : ℂ → ℂ) (hf : Differentiable ℂ f)
    (hper : ∀ (z : ℂ), ∀ l ∈ P.lattice, f (z + l) = f z) (h0 : f 0 = 0) (z : ℂ) :
    f z = 0 := by
  rw [elliptic_liouville P f hf hper z 0, h0]

def elliptic_liouville.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(楕円関数の Liouville——整で二重周期的なら定数。★無条件)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★Vélu の公式の解析側 -/

/-- ★★★★★★**Vélu の `X` の解析側**

    `X(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

★`T` は `Λ′/Λ` の代表系（`0` を含む）、`c = Σ_{w ∈ T∖{0}} ℘_Λ(w)` のつもりである。
原文の形 `℘(z) + Σ_{w≠0}[℘(z+w) − ℘(w)]` を、和をひとまとめにして書いたものである。 -/
noncomputable def veluAnalyticX (P : PeriodPair) (T : Finset ℂ) (c : ℂ) (z : ℂ) : ℂ :=
  (∑ w ∈ T, P.weierstrassP (z + w)) - c

/-- ★★★★★★★★★★★★**代表系は平行移動で置換される**——`℘` 側。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`σ` は「`w ↦ w + w₀` が `Λ′/Λ` に誘導する置換」を代表系の上に持ち上げたものである
（`hshift`: `w + w₀ − σ(w) ∈ Λ`）。 -/
theorem veluAnalyticSum_shift (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    ∑ w ∈ T, P.weierstrassP (z + w + w₀) = ∑ w ∈ T, P.weierstrassP (z + w) := by
  refine Finset.sum_bij (fun w _ => σ w) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) hsurj ?_
  intro w hw
  have hz : z + w + w₀ = (z + σ w) + (w + w₀ - σ w) := by ring
  rw [hz]
  exact P.weierstrassP_add_coe _ ⟨_, hshift w hw⟩

/-- ★★★★★★★★★★★★**代表系は平行移動で置換される**——`℘′` 側。 -/
theorem derivVeluAnalyticSum_shift (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    ∑ w ∈ T, P.derivWeierstrassP (z + w + w₀) = ∑ w ∈ T, P.derivWeierstrassP (z + w) := by
  refine Finset.sum_bij (fun w _ => σ w) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) hsurj ?_
  intro w hw
  have hz : z + w + w₀ = (z + σ w) + (w + w₀ - σ w) := by ring
  rw [hz]
  exact P.derivWeierstrassP_add_coe _ ⟨_, hshift w hw⟩

/-- ★★★★★★★★★★★★★★★★**Vélu の `X` は `Λ′`-周期的**。

★これが「`X` が `E/H` の座標である」ことの中身である。 -/
theorem veluAnalyticX_shift (P : PeriodPair) (T : Finset ℂ) (c w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    veluAnalyticX P T c (z + w₀) = veluAnalyticX P T c z := by
  simp only [veluAnalyticX]
  congr 1
  rw [← veluAnalyticSum_shift P T w₀ σ hmem hinj hsurj hshift z]
  exact Finset.sum_congr rfl fun w _ => by ring_nf

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側——Liouville で閉じる形**。

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★仮定は 3 つ:

* `hper`  : 右辺が `Λ′`-周期的（★`veluAnalyticX_shift` で取れる——**代表系の置換だけ**）
* `hdiff` : ☆**差が整である**（極が打ち消し合うこと——極の解析が残る）
* `h0`    : ☆差が原点で `0`

★★★`hper` は本ファイルで無条件に取れており、残るのは `hdiff`（極の打ち消し）と `h0`。
☆すなわち **Vélu の公式の解析側は「極の解析」1 点に絞られた**。

★★★★これが `Found/GenEll/Velu.lean`（第 586-593）の代数側と対になる。
両者を繋げば `Found/GaloisRep/VeluNormalized.lean` の
`htFalt_isogeny_le_of_analytic`（第 596）の入力が揃う。 -/
theorem weierstrassP_eq_of_liouville (P P' : PeriodPair) (T : Finset ℂ) (c : ℂ)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice, veluAnalyticX P T c (z + l) = veluAnalyticX P T c z)
    (hdiff : Differentiable ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z))
    (h0 : P'.weierstrassP 0 - veluAnalyticX P T c 0 = 0) (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T c z := by
  have hkey : ∀ (y : ℂ), ∀ l ∈ P'.lattice,
      (fun z => P'.weierstrassP z - veluAnalyticX P T c z) (y + l)
        = (fun z => P'.weierstrassP z - veluAnalyticX P T c z) y := by
    intro y l hl
    simp only
    rw [P'.weierstrassP_add_coe y ⟨l, hl⟩, hper y l hl]
  have h := elliptic_liouville_eq_zero P' _ hdiff hkey h0 z
  simpa [sub_eq_zero] using h

def veluAnalyticX_shift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の X は Λ′-周期的——代表系は平行移動で置換される。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_of_liouville.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——Liouville で閉じる形。残るのは極の解析)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_of_liouville.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hdiff: ℘_{Λ′}(z) − Σ_{w∈T} ℘_Λ(z+w) + c が整であること。" ++
       "★両辺とも Λ′ の各点で 2 位の極をもち、主要部が一致するので差は整になる。" ++
       "★★mathlib は order_weierstrassP(各格子点で 2 位の極)を持っているので" ++
       "道具はある(2026-08-29 に測定)") 8,
    .implicitStep
      ("☆仮定 h0: 差が原点で 0 であること。★c の取り方(c = Σ_{w≠0} ℘_Λ(w))で決まる") 6,
    .implicitStep
      ("★★★到達点(2026-08-29、第 599): Vélu の公式の解析側が" ++
       "「代表系の置換」(無条件で取れた)と「極の解析」(残り)に分離した。" ++
       "★代数側は Found/GenEll/Velu.lean が持っている(第 586-593)") 8 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★極が打ち消し合うこと -/

/-- ★★★★★★★★★★**`℘` から主要部を引くと格子点で解析的**。

    `℘_Λ(z) − 1/(z−p)²` は `p ∈ Λ` で解析的

★mathlib の `weierstrassPExcept`（`l₀` 項を抜いた `℘`）がそのまま使える:
`℘[Λ] z = ℘[Λ − p] z + (1/(z−p)² − 1/p²)`。 -/
theorem weierstrassP_sub_pole_analyticAt (P : PeriodPair) (p : ℂ) (hp : p ∈ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP z - 1 / (z - p) ^ 2) p := by
  have hEq : (fun z => P.weierstrassP z - 1 / (z - p) ^ 2)
      = fun z => P.weierstrassPExcept p z - 1 / p ^ 2 := by
    funext z; rw [← P.weierstrassPExcept_add ⟨p, hp⟩]; ring
  rw [hEq]
  exact (((P.differentiableOn_weierstrassPExcept p).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) p (by simp)).sub analyticAt_const

/-- ★★★★★★★★平行移動した `℘` からも同じ主要部を引ける——`(z+w) − (p+w) = z − p` だから。 -/
theorem shifted_sub_pole_analyticAt (P : PeriodPair) (p w : ℂ) (h : p + w ∈ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP (z + w) - 1 / (z - p) ^ 2) p := by
  have hf : AnalyticAt ℂ (fun z : ℂ => z + w) p := analyticAt_id.add analyticAt_const
  have hg : AnalyticAt ℂ (fun z => P.weierstrassP z - 1 / (z - (p + w)) ^ 2)
      ((fun z : ℂ => z + w) p) := weierstrassP_sub_pole_analyticAt P (p + w) h
  refine (AnalyticAt.comp (f := fun z : ℂ => z + w) (x := p) hg hf).congr ?_
  filter_upwards with z
  simp only [Function.comp_apply]
  ring_nf

/-- ★★★★★格子の外では平行移動した `℘` は解析的。 -/
theorem shifted_analyticAt (P : PeriodPair) (p w : ℂ) (h : p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP (z + w)) p := by
  have hf : AnalyticAt ℂ (fun z : ℂ => z + w) p := analyticAt_id.add analyticAt_const
  have hg : AnalyticAt ℂ P.weierstrassP ((fun z : ℂ => z + w) p) :=
    P.analyticOnNhd_weierstrassP (p + w) h
  exact AnalyticAt.comp (f := fun z : ℂ => z + w) (x := p) hg hf

/-- ★★★★★★★★★★★★★★★★★★**極が打ち消し合う**——`p ∈ Λ′` のとき。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Λ′` の点 `p` では、`℘_{Λ′}` と（代表系のうちちょうど 1 つの）`℘_Λ(·+w₀)` が
**同じ主要部 `1/(z−p)²`** を持つので、差は解析的になる。 -/
theorem veluDiff_analyticAt_of_mem (P P' : PeriodPair) (T : Finset ℂ) (c p : ℂ)
    (hp : p ∈ P'.lattice) (w₀ : ℂ) (hw₀ : w₀ ∈ T) (hpw₀ : p + w₀ ∈ P.lattice)
    (hother : ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) p := by
  have hEq : (fun z => P'.weierstrassP z - veluAnalyticX P T c z)
      = fun z => ((P'.weierstrassP z - 1 / (z - p) ^ 2)
            - (P.weierstrassP (z + w₀) - 1 / (z - p) ^ 2))
          - (∑ w ∈ T.erase w₀, P.weierstrassP (z + w)) + c := by
    funext z
    simp only [veluAnalyticX]
    rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (z + w)) hw₀]
    ring
  rw [hEq]
  refine (((weierstrassP_sub_pole_analyticAt P' p hp).sub
    (shifted_sub_pole_analyticAt P p w₀ hpw₀)).sub ?_).add analyticAt_const
  refine Finset.analyticAt_fun_sum _ fun w hw => ?_
  exact shifted_analyticAt P p w
    (hother w (Finset.mem_of_mem_erase hw) (Finset.ne_of_mem_erase hw))

/-- ★★★★★★★★★★★★格子の外では両方とも解析的。 -/
theorem veluDiff_analyticAt_of_notMem (P P' : PeriodPair) (T : Finset ℂ) (c p : ℂ)
    (hp : p ∉ P'.lattice) (hall : ∀ w ∈ T, p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) p := by
  have hEq : (fun z => P'.weierstrassP z - veluAnalyticX P T c z)
      = fun z => (P'.weierstrassP z - ∑ w ∈ T, P.weierstrassP (z + w)) + c := by
    funext z; simp only [veluAnalyticX]; ring
  rw [hEq]
  refine (((P'.analyticOnNhd_weierstrassP p hp).sub ?_)).add analyticAt_const
  exact Finset.analyticAt_fun_sum _ fun w hw => shifted_analyticAt P p w (hall w hw)

/-- ★★★★★★`p ∉ Λ′` なら `p + w ∉ Λ`（`w ∈ T ⊆ Λ′`・`Λ ⊆ Λ′` だから）。 -/
theorem notMem_of_notMem (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice) (p : ℂ) (hp : p ∉ P'.lattice) (w : ℂ) (hw : w ∈ T) :
    p + w ∉ P.lattice := by
  intro hc
  refine hp ?_
  have h1 : (p + w) - w ∈ P'.lattice := P'.lattice.sub_mem (hle hc) (hT w hw)
  simpa using h1

/-- ★★★★★★★★★★★★★★★★★★★★★★**差は整である**——★極の打ち消しが済んだ形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★仮定は **`T` が `Λ′/Λ` の代表系である**ことだけ（`hT`・`hrep`）:

* `hT`  : `T ⊆ Λ′`
* `hrep`: `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` を満たす

☆すなわち `weierstrassP_eq_of_liouville`（第 599）の仮定 `hdiff` は
**代表系の性質だけから従う**。★残るのは `h0`（差が原点で `0`）＝`c` の取り方である。 -/
theorem veluDiff_differentiable (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (c : ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    Differentiable ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) := by
  intro p
  by_cases hp : p ∈ P'.lattice
  · obtain ⟨w₀, hw₀, hpw₀, hother⟩ := hrep p hp
    exact (veluDiff_analyticAt_of_mem P P' T c p hp w₀ hw₀ hpw₀ hother).differentiableAt
  · exact (veluDiff_analyticAt_of_notMem P P' T c p hp
      (fun w hw => notMem_of_notMem P P' hle T hT p hp w hw)).differentiableAt

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側——代表系だけから**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

★★★仮定は **`T` が `Λ′/Λ` の代表系であること**（`hT`・`hrep`）、
**平行移動が代表系を置換すること**（`hshift`）、
そして `h0`（`c` の取り方で決まる正規化）だけになった。

☆☆**極の解析は済んだ**（`veluDiff_differentiable`）——第 599 で残っていた `hdiff` である。 -/
theorem weierstrassP_eq_veluAnalyticX (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (c : ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice, veluAnalyticX P T c (z + l) = veluAnalyticX P T c z)
    (h0 : P'.weierstrassP 0 - veluAnalyticX P T c 0 = 0) (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T c z :=
  weierstrassP_eq_of_liouville P P' T c hper
    (veluDiff_differentiable P P' hle T c hT hrep) h0 z

def veluDiff_differentiable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極が打ち消し合う——代表系の性質だけから従う。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_veluAnalyticX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——代表系と正規化だけから)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★正規化定数と最終形 -/

/-- ★★★★★**`℘(0) = 0`**——mathlib の定義（`∑' l, (1/(z−l)² − 1/l²)`）では
`z = 0` の各項が `1/l² − 1/l² = 0` になる（`l = 0` の項も junk value で `0`）。

★これは「極での値」ではなく**除去可能特異点を埋めた関数の値**として整合している
（`weierstrassPExcept_add` が `z = l₀` でも成り立つのはそのためである）。 -/
theorem weierstrassP_zero (P : PeriodPair) : P.weierstrassP 0 = 0 := by
  simp [PeriodPair.weierstrassP]

/-- ★★★★★★**Vélu の正規化定数** `c = Σ_{w ∈ T∖{0}} ℘_Λ(w)`。

★原文の形 `℘(z) + Σ_{w≠0}[℘(z+w) − ℘(w)]` の第 2 項の定数部分である。 -/
noncomputable def veluAnalyticC (P : PeriodPair) (T : Finset ℂ) : ℂ :=
  ∑ w ∈ T.erase 0, P.weierstrassP w

/-- ★★★★★★★★**正規化定数を入れると原点で `0`**。 -/
theorem veluAnalyticX_zero (P : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) :
    veluAnalyticX P T (veluAnalyticC P T) 0 = 0 := by
  simp only [veluAnalyticX, veluAnalyticC]
  rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (0 + w)) h0T]
  simp

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側（最終形）**

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★**仮定は `T` が `Λ′/Λ` の代表系であることだけ**である:

| 仮定 | 内容 |
|---|---|
| `hle` | `Λ ⊆ Λ′` |
| `h0T` | `0 ∈ T`（自明な剰余類の代表） |
| `hT` | `T ⊆ Λ′` |
| `hrep` | `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` |
| `hper` | 平行移動が代表系を置換する（★`veluAnalyticX_shift` で取れる） |

★★★★**極の解析も正規化も済んだ**——第 598（Liouville）・第 599（周期性）・
第 600（極の打ち消し）・第 601（正規化定数）で塞いだ。

☆残るのは `hper` を代表系の定義から出すこと（`veluAnalyticX_shift` に `σ` を与えること）と、
この `℘` の等式を `Found/GenEll/Velu.lean` の**代数側の `veluXGen`** に翻訳することである
——後者には `℘` の加法定理（mathlib に無い）が要る。 -/
theorem weierstrassP_eq_velu (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice,
      veluAnalyticX P T (veluAnalyticC P T) (z + l)
        = veluAnalyticX P T (veluAnalyticC P T) z)
    (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
  weierstrassP_eq_veluAnalyticX P P' hle T _ hT hrep hper
    (by rw [weierstrassP_zero, veluAnalyticX_zero P T h0T]; ring) z

def weierstrassP_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘(0) = 0——除去可能特異点を埋めた値。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側・最終形——仮定は代表系であることだけ)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hper は veluAnalyticX_shift(第 599)で取れるが、" ++
       "代表系の定義から σ(平行移動が誘導する置換)を作る段が残る") 5,
    .implicitStep
      ("☆この ℘ の等式を Found/GenEll/Velu.lean の代数側 veluXGen(第 591)に" ++
       "翻訳すること。★℘ の加法定理が要る(mathlib に無い、2026-08-29 に測定)") 8,
    .implicitStep
      ("★★★到達点(2026-08-29、第 601): Vélu の公式の解析側が" ++
       "「T が Λ′/Λ の代表系である」だけから従う形になった。" ++
       "★極の解析(第 600)・Liouville(第 598)・周期性(第 599)・正規化(第 601)で塞いだ") 9 ]

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

/-! ## ★★★★★★★★★★★★★★★★★★Laurent の入口——加法定理へ -/

/-- ★★★★★**`℘(z) − 1/z²` は mathlib の `℘[Λ − 0]` そのもの**（原点で解析的）。 -/
theorem weierstrassP_sub_invSq (P : PeriodPair) (z : ℂ) :
    P.weierstrassP z - 1 / z ^ 2 = P.weierstrassPExcept 0 z := by
  have h := P.weierstrassPExcept_add ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  rw [← h]
  simp

/-- ★★★★★**`℘′(z) + 2/z³` は `℘′[Λ − 0]`**。 -/
theorem derivWeierstrassP_add_invCube (P : PeriodPair) (z : ℂ) :
    P.derivWeierstrassP z + 2 / z ^ 3 = P.derivWeierstrassPExcept 0 z := by
  have h := P.derivWeierstrassPExcept_sub ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  rw [← h]
  ring

/-- ★★★★★★**`z²·℘(z)` の解析接続**——原点で `1`。

★これが `℘(z) = z⁻² + O(z²)` の Lean 上の姿である
（`weierstrassPExcept` は原点で解析的で値 `0`）。 -/
noncomputable def laurentB (P : PeriodPair) (z : ℂ) : ℂ :=
  1 + z ^ 2 * P.weierstrassPExcept 0 z

/-- ★★★★★★**`z³·℘′(z)` の解析接続**——原点で `−2`。 -/
noncomputable def laurentA (P : PeriodPair) (z : ℂ) : ℂ :=
  -2 + z ^ 3 * P.derivWeierstrassPExcept 0 z

@[simp] theorem laurentB_zero (P : PeriodPair) : laurentB P 0 = 1 := by simp [laurentB]

@[simp] theorem laurentA_zero (P : PeriodPair) : laurentA P 0 = -2 := by simp [laurentA]

/-- ★★★★★★`z ≠ 0` では `laurentB P z = z²·℘(z)`。 -/
theorem laurentB_eq (P : PeriodPair) (z : ℂ) (hz : z ≠ 0) :
    laurentB P z = z ^ 2 * P.weierstrassP z := by
  have h := P.weierstrassPExcept_add ⟨0, P.lattice.zero_mem⟩ z
  simp only [sub_zero] at h
  simp only [laurentB, ← h]
  have h0 : (1 : ℂ) / 0 ^ 2 = 0 := by norm_num
  rw [h0]
  field_simp
  ring

/-- ★★★★★★`z ≠ 0` では `laurentA P z = z³·℘′(z)`。 -/
theorem laurentA_eq (P : PeriodPair) (z : ℂ) (hz : z ≠ 0) :
    laurentA P z = z ^ 3 * P.derivWeierstrassP z := by
  simp only [laurentA, ← derivWeierstrassP_add_invCube]
  field_simp
  ring

/-- ★★★★★★★★`laurentA`・`laurentB` は原点で解析的。 -/
theorem analyticAt_laurentB (P : PeriodPair) : AnalyticAt ℂ (laurentB P) 0 := by
  refine analyticAt_const.add (((analyticAt_id).pow 2).mul ?_)
  exact ((P.differentiableOn_weierstrassPExcept 0).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) 0 (by simp)

theorem analyticAt_laurentA (P : PeriodPair) : AnalyticAt ℂ (laurentA P) 0 := by
  refine analyticAt_const.add (((analyticAt_id).pow 3).mul ?_)
  exact ((P.differentiableOn_derivWeierstrassPExcept 0).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) 0 (by simp)

def laurentB.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z²℘(z) の解析接続——加法定理の Laurent の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

def laurentA.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z³℘′(z) の解析接続——加法定理の Laurent の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★`M(z)`——`laurent_cancel` の右辺の因子。 -/
noncomputable def laurentM (P : PeriodPair) (x y z : ℂ) : ℂ :=
  4 * P.weierstrassPExcept 0 z + 8 * (P.weierstrassPExcept 0 z - x)
  + 4 * z ^ 2 * (P.weierstrassPExcept 0 z - x) ^ 2
  + 8 * z ^ 2 * P.weierstrassPExcept 0 z * (P.weierstrassPExcept 0 z - x)
  + 4 * z ^ 4 * P.weierstrassPExcept 0 z * (P.weierstrassPExcept 0 z - x) ^ 2
  + 4 * z * (P.derivWeierstrassPExcept 0 z - y)
  - z ^ 4 * (P.derivWeierstrassPExcept 0 z - y) ^ 2

/-- ★★★★★★★★★★★★★★★★★★★★★★**加法定理の極の打ち消し（原点）**——★純粋な恒等式。

    `4·B·(B − z²·x) − (A − z³·y)² = z² · M(z)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/AdditionTheorem.lean`（第 606）の証明の核である。
`x = ℘(w)`・`y = ℘′(w)` と置くと、`z ≠ 0` では `B = z²℘(z)`・`A = z³℘′(z)` なので

    `℘(z) + ℘(w) + ℘(z+w) − (1/4)·((℘′(z)−℘′(w))/(℘(z)−℘(w)))²`
      `= ℘(z+w) + ℘(w) + [4·B·(B − z²℘(w))² − (A − z³℘′(w))²] / (4·z²·(B − z²℘(w))²)`

であり、★**分子が `z²` でくくれる**（本定理）ので `z²` が約されて原点で有界になる。

☆`B(0) = 1`・`A(0) = −2` なので `4·1·1 − 4 = 0`——**それが打ち消しの正体**である。

☆残るのは (1) この式から「`F` は原点の除いた近傍で解析関数と一致する」を出すこと、
(2) `z ≡ −w` での極の打ち消し（そちらは `℘` の 2 階の Taylor が要る）。 -/

theorem laurent_cancel (P : PeriodPair) (z x y : ℂ) :
    4 * laurentB P z * (laurentB P z - z ^ 2 * x) ^ 2
        - (laurentA P z - z ^ 3 * y) ^ 2
      = z ^ 2 * laurentM P x y z := by
  simp only [laurentA, laurentB, laurentM]
  ring

/-- ★★★★★`M(0) = −8·x`——これが原点での値を `0` にする。 -/
@[simp] theorem laurentM_zero (P : PeriodPair) (x y : ℂ) : laurentM P x y 0 = -8 * x := by
  simp [laurentM]

/-- ★★★★★★★★★★★★**打ち消しの正体**——`z = 0` での値。

`4·B(0)·B(0)² − A(0)² = 4·1·1 − (−2)² = 0`。 -/
theorem laurent_cancel_zero (P : PeriodPair) (x y : ℂ) :
    4 * laurentB P 0 * (laurentB P 0 - (0:ℂ) ^ 2 * x) ^ 2
      - (laurentA P 0 - (0:ℂ) ^ 3 * y) ^ 2 = 0 := by
  simp
  norm_num

def laurent_cancel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(加法定理の極の打ち消し——分子が z² でくくれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★加法定理の欠損関数 -/

/-- ★★★★★★★★**加法定理の欠損** `F_w(z)`。

    `F_w(z) ≔ ℘(z+w) + ℘(z) + ℘(w) − (1/4)·((℘′(z)−℘′(w))/(℘(z)−℘(w)))²`

★`Skeleton/GenEll/AdditionTheorem.lean`（第 606）の `weierstrassP_add` は
**`F_w ≡ 0`** と同値である。 -/
noncomputable def addDefect (P : PeriodPair) (w z : ℂ) : ℂ :=
  P.weierstrassP (z + w) + P.weierstrassP z + P.weierstrassP w
    - ((P.derivWeierstrassP z - P.derivWeierstrassP w)
        / (P.weierstrassP z - P.weierstrassP w)) ^ 2 / 4

/-- ★★★★★★★★**原点の近くでの解析的な姿**——`z²` が約された形。 -/
noncomputable def addDefectNear (P : PeriodPair) (w z : ℂ) : ℂ :=
  P.weierstrassP (z + w) + P.weierstrassP w
    + laurentM P (P.weierstrassP w) (P.derivWeierstrassP w) z
      / (4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2)

/-- ★★★★★★★★★★★★★★★★★★**`z ≠ 0` では両者は一致する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`℘(z) = B/z²`・`℘′(z) = A/z³` を入れると分母の `z²` が `laurent_cancel` の
`z²` と約されるからである。☆**これが「原点の極が消える」ことの中身**である。 -/
theorem addDefect_eq_near (P : PeriodPair) (w z : ℂ) (hz : z ≠ 0)
    (hpz : P.weierstrassP z - P.weierstrassP w ≠ 0) :
    addDefect P w z = addDefectNear P w z := by
  have hB := laurentB_eq P z hz
  have hA := laurentA_eq P z hz
  have hcan := laurent_cancel P z (P.weierstrassP w) (P.derivWeierstrassP w)
  rw [hB, hA] at hcan
  have hz2 : z ^ 2 ≠ 0 := pow_ne_zero _ hz
  have hkey : z ^ 4 * (4 * P.weierstrassP z
        * (P.weierstrassP z - P.weierstrassP w) ^ 2
      - (P.derivWeierstrassP z - P.derivWeierstrassP w) ^ 2)
      = laurentM P (P.weierstrassP w) (P.derivWeierstrassP w) z := by
    refine mul_left_cancel₀ hz2 ?_
    linear_combination hcan
  simp only [addDefect, addDefectNear, hB]
  field_simp
  linear_combination hkey

/-- ★★★★★★★★★★★★★★★★**原点での値は `0`**——`M(0) = −8℘(w)`・`B(0) = 1` から

    `℘(w) + ℘(w) + (−8℘(w))/4 = 2℘(w) − 2℘(w) = 0`。

★★☆**これが「加法定理が原点で成り立つ」ことである**——`F_w` の解析接続は原点で消える。 -/
@[simp] theorem addDefectNear_zero (P : PeriodPair) (w : ℂ) :
    addDefectNear P w 0 = 0 := by
  simp only [addDefectNear, laurentM_zero, laurentB_zero, zero_add]
  norm_num
  ring

/-- ★★★★★★★★★★★★**`addDefectNear` は原点で解析的**（`w ∉ Λ` のとき）。 -/
theorem analyticAt_addDefectNear (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    AnalyticAt ℂ (addDefectNear P w) 0 := by
  have hpw : AnalyticAt ℂ (fun z : ℂ => P.weierstrassP (z + w)) 0 := by
    have hf : AnalyticAt ℂ (fun z : ℂ => z + w) 0 := analyticAt_id.add analyticAt_const
    have hg : AnalyticAt ℂ P.weierstrassP ((fun z : ℂ => z + w) 0) := by
      simpa using P.analyticOnNhd_weierstrassP w hw
    exact AnalyticAt.comp (f := fun z : ℂ => z + w) (x := 0) hg hf
  have he : AnalyticAt ℂ (P.weierstrassPExcept 0) 0 :=
    ((P.differentiableOn_weierstrassPExcept 0).analyticOnNhd
      P.isOpen_compl_lattice_sdiff) 0 (by simp)
  have hf' : AnalyticAt ℂ (P.derivWeierstrassPExcept 0) 0 :=
    ((P.differentiableOn_derivWeierstrassPExcept 0).analyticOnNhd
      P.isOpen_compl_lattice_sdiff) 0 (by simp)
  have hM : AnalyticAt ℂ (laurentM P (P.weierstrassP w) (P.derivWeierstrassP w)) 0 := by
    unfold laurentM
    fun_prop (disch := assumption)
  have hD : AnalyticAt ℂ (fun z : ℂ => 4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2) 0 := by
    unfold laurentB
    fun_prop (disch := assumption)
  have hDne : (fun z : ℂ => 4 * (laurentB P z - z ^ 2 * P.weierstrassP w) ^ 2) 0 ≠ 0 := by
    simp
  exact (hpw.add analyticAt_const).add (hM.div hD hDne)

def addDefect.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(加法定理の欠損 F_w——F_w ≡ 0 が加法定理そのもの)",
    sectionId := "genell-lemma-3-5" }

def addDefect_eq_near.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(原点の極が消える——z² が約される。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★`z ≡ −w` の側——`q` の零点の位数 -/

/-- ★★★★★★★★**`z ≡ −w` の側の鍵となる関数**

    `q(t) ≔ 2·(℘(t−w) − ℘(w)) − t·(℘′(t−w) − ℘′(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`t ≔ z + w` と置くと `F_w(t−w)` の極は `q` が `t = 0` で 3 位の零点を持つことで消える
（`u ≔ ℘(t−w) − ℘(w)`・`v ≔ ℘′(t−w) − ℘′(w)`・`û ≔ u/t` として `2û − v = q/t`、
`4û² − v² = (2û−v)(2û+v)`）。 -/
noncomputable def addQ (P : PeriodPair) (w t : ℂ) : ℂ :=
  2 * (P.weierstrassP (t - w) - P.weierstrassP w)
    - t * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)

/-- ★★★★★★**`q(0) = 0`**——`℘` は偶だから `℘(−w) = ℘(w)`。 -/
@[simp] theorem addQ_zero (P : PeriodPair) (w : ℂ) : addQ P w 0 = 0 := by
  simp only [addQ, zero_sub, zero_mul, sub_zero, P.weierstrassP_neg]
  ring

/-- ★★★★★★★★**`q` の 1 階導関数** `q′(t) = ℘′(t−w) + ℘′(w) − t·℘″(t−w)`。 -/
theorem hasDerivAt_addQ (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ P w)
      (P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
        - t * deriv P.derivWeierstrassP (t - w)) t := by
  have h1 : HasDerivAt (fun s : ℂ => P.weierstrassP (s - w))
      (P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_weierstrassP P ht)
  have h2 : HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s - w))
      (deriv P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_derivWeierstrassP P ht)
  have h3 : HasDerivAt (fun s : ℂ => s * (P.derivWeierstrassP (s - w) - P.derivWeierstrassP w))
      (1 * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)
        + t * deriv P.derivWeierstrassP (t - w)) t :=
    (hasDerivAt_id t).mul (h2.sub_const _)
  have h4 := ((h1.sub_const (P.weierstrassP w)).const_mul (2:ℂ)).sub h3
  have hval : (2:ℂ) * P.derivWeierstrassP (t - w)
      - (1 * (P.derivWeierstrassP (t - w) - P.derivWeierstrassP w)
        + t * deriv P.derivWeierstrassP (t - w))
      = P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
        - t * deriv P.derivWeierstrassP (t - w) := by ring
  rw [← hval]
  exact h4

/-- ★★★★★★★★★★**`q′(0) = 0`**——`℘′` が奇であることがちょうど効く。

    `q′(0) = ℘′(−w) + ℘′(w) − 0 = −℘′(w) + ℘′(w) = 0` -/
theorem deriv_addQ_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    deriv (addQ P w) 0 = 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := by
    intro hc
    exact hw (by simpa using neg_mem hc)
  have h := hasDerivAt_addQ P w 0 hnw
  rw [h.deriv]
  simp only [zero_sub, zero_mul, sub_zero, P.derivWeierstrassP_neg]
  ring

def addQ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(z ≡ −w の側の鍵——q は t = 0 で 3 位の零点)",
    sectionId := "genell-lemma-3-5" }

def deriv_addQ_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q′(0) = 0——℘′ が奇であることが効く。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★**`q` の 1 階導関数（閉じた形）**——`Found/GenEll/WeierstrassODE.lean` の
`deriv_derivWeierstrassP`（`deriv ℘′ = 6℘² − g₂/2`）を入れた形。 -/
noncomputable def addQ' (P : PeriodPair) (w t : ℂ) : ℂ :=
  P.derivWeierstrassP (t - w) + P.derivWeierstrassP w
    - t * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)

theorem hasDerivAt_addQ_closed (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ P w) (addQ' P w t) t := by
  have h := hasDerivAt_addQ P w t ht
  rw [deriv_derivWeierstrassP P ht] at h
  exact h

/-- ★★★★★★★★★★★★**`q″(t) = −12·t·℘(t−w)·℘′(t−w)`**——★**`t` の因子が残る**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`deriv ℘′ = 6℘² − g₂/2` の項がちょうど打ち消し合うので `t` の因子だけが残る。
☆これが第 611 で「自動的に消える」と書いた段である。 -/
theorem hasDerivAt_addQ' (P : PeriodPair) (w t : ℂ) (ht : t - w ∉ P.lattice) :
    HasDerivAt (addQ' P w)
      (-(12 * t * P.weierstrassP (t - w) * P.derivWeierstrassP (t - w))) t := by
  have hp : HasDerivAt (fun s : ℂ => P.weierstrassP (s - w))
      (P.derivWeierstrassP (t - w)) t :=
    HasDerivAt.comp_sub_const t w (hasDerivAt_weierstrassP P ht)
  have h2 : HasDerivAt (fun s : ℂ => P.derivWeierstrassP (s - w))
      (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2) t := by
    have h := HasDerivAt.comp_sub_const t w (hasDerivAt_derivWeierstrassP P ht)
    rwa [deriv_derivWeierstrassP P ht] at h
  have hsq : HasDerivAt (fun s : ℂ => 6 * P.weierstrassP (s - w) ^ 2 - P.g₂ / 2)
      (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w))) t :=
    ((hp.pow 2).const_mul 6).sub_const _
  have hmul : HasDerivAt (fun s : ℂ => s * (6 * P.weierstrassP (s - w) ^ 2 - P.g₂ / 2))
      (1 * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
        + t * (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w)))) t :=
    (hasDerivAt_id t).mul hsq
  have h4 := (h2.add_const (P.derivWeierstrassP w)).sub hmul
  have hval : (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
      - (1 * (6 * P.weierstrassP (t - w) ^ 2 - P.g₂ / 2)
        + t * (6 * (2 * P.weierstrassP (t - w) ^ 1 * P.derivWeierstrassP (t - w))))
      = -(12 * t * P.weierstrassP (t - w) * P.derivWeierstrassP (t - w)) := by
    ring
  rw [← hval]
  exact h4

/-- ★★★★★★★★★★★★**`q″(0) = 0`**。 -/
theorem iteratedDeriv_two_addQ_zero (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice) :
    iteratedDeriv 2 (addQ P w) 0 = 0 := by
  have hnw : (0 : ℂ) - w ∉ P.lattice := fun hc => hw (by simpa using neg_mem hc)
  have hopen : IsOpen {t : ℂ | t - w ∉ P.lattice} := by
    have : {t : ℂ | t - w ∉ P.lattice} = (fun t : ℂ => t - w) ⁻¹' ((P.lattice : Set ℂ)ᶜ) := rfl
    rw [this]
    exact (P.isClosed_lattice.isOpen_compl).preimage (by fun_prop)
  have hnhds : {t : ℂ | t - w ∉ P.lattice} ∈ nhds (0 : ℂ) := hopen.mem_nhds hnw
  have heq : deriv (addQ P w) =ᶠ[nhds (0:ℂ)] addQ' P w := by
    filter_upwards [hnhds] with t ht
    exact (hasDerivAt_addQ_closed P w t ht).deriv
  rw [iteratedDeriv_succ, iteratedDeriv_one, heq.deriv_eq,
    (hasDerivAt_addQ' P w 0 hnw).deriv]
  ring

/-- ★★★★★★★★★★★★★★★★★★**`q` は `t = 0` で 3 位の零点**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`q(0) = q′(0) = q″(0) = 0`（第 612・第 614）から、mathlib の
`natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero` で位数が `≥ 3` と分かる。

☆これで `2û − v = q/t` が `2` 位の零点になり、
`4û² − v² = (2û−v)(2û+v)` が `2` 位で消えて **`z ≡ −w` の極が打ち消される**。 -/
theorem three_le_analyticOrderAt_addQ (P : PeriodPair) (w : ℂ) (hw : w ∉ P.lattice)
    (hana : AnalyticAt ℂ (addQ P w) 0) :
    (3 : ℕ) ≤ analyticOrderAt (addQ P w) 0 := by
  rw [natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hana]
  intro i hi
  interval_cases i
  · simpa using addQ_zero P w
  · rw [iteratedDeriv_one]; exact deriv_addQ_zero P w hw
  · exact iteratedDeriv_two_addQ_zero P w hw

def hasDerivAt_addQ'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q″(t) = −12t·℘(t−w)·℘′(t−w)——t の因子が残る。★無条件)",
    sectionId := "genell-lemma-3-5" }

def three_le_analyticOrderAt_addQ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(q は t = 0 で 3 位の零点——z ≡ −w の極が打ち消される)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★零点勘定への第一の煉瓦 -/

/-- ★★★★★★★★★★★★★★★★★★★★**楕円関数の周期平行四辺形の境界積分は消える**。

    `∮_{∂D} f dz = 0`（`D` は `a` を頂点とする周期平行四辺形）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★機構は**周期性だけ**——Cauchy の定理は要らない:

* `f(a + ω₂ + tω₁) = f(a + tω₁)`（`ω₂ ∈ Λ`）なので上辺と下辺が打ち消し合う
* `f(a + ω₁ + tω₂) = f(a + tω₂)`（`ω₁ ∈ Λ`）なので左辺と右辺が打ち消し合う

★★★★☆**これが楕円関数の零点勘定の第一の煉瓦である**。
残る半分は**留数定理**（`∮ f = 2πi·Σ res`）であり、それと合わせて

    `Σ_{D} res(f) = 0`

が出る。`f = g′/g` に当てると **`#零点 = #極`**（偏角の原理）になり、
`g = ℘ − c` なら「`℘` は各値をちょうど 2 回取る」が従う
——それが `Skeleton/GenEll/AdditionTheorem.lean`（第 615）で同定した唯一の入口である。

☆mathlib は軸平行な長方形での Cauchy（`integral_boundary_rect_eq_zero_of_...`）を持つが、
一般の格子の平行四辺形には当たらない（2026-08-29 に測定）。 -/
theorem elliptic_boundary_integral_zero (P : PeriodPair) (f : ℂ → ℂ)
    (hper : ∀ (z : ℂ), ∀ l ∈ P.lattice, f (z + l) = f z) (a : ℂ) :
    (∫ t in (0:ℝ)..1, f (a + (t : ℂ) * P.ω₁) * P.ω₁)
      + (∫ t in (0:ℝ)..1, f (a + P.ω₁ + (t : ℂ) * P.ω₂) * P.ω₂)
      - (∫ t in (0:ℝ)..1, f (a + P.ω₂ + (t : ℂ) * P.ω₁) * P.ω₁)
      - (∫ t in (0:ℝ)..1, f (a + (t : ℂ) * P.ω₂) * P.ω₂) = 0 := by
  have hω₁ : P.ω₁ ∈ P.lattice := Submodule.subset_span (by simp)
  have hω₂ : P.ω₂ ∈ P.lattice := Submodule.subset_span (by simp)
  have h1 : ∀ t : ℝ, f (a + P.ω₂ + (t : ℂ) * P.ω₁) = f (a + (t : ℂ) * P.ω₁) := by
    intro t
    have hz : a + P.ω₂ + (t : ℂ) * P.ω₁ = (a + (t : ℂ) * P.ω₁) + P.ω₂ := by ring
    rw [hz, hper _ _ hω₂]
  have h2 : ∀ t : ℝ, f (a + P.ω₁ + (t : ℂ) * P.ω₂) = f (a + (t : ℂ) * P.ω₂) := by
    intro t
    have hz : a + P.ω₁ + (t : ℂ) * P.ω₂ = (a + (t : ℂ) * P.ω₂) + P.ω₁ := by ring
    rw [hz, hper _ _ hω₁]
  simp only [h1, h2]
  ring

def elliptic_boundary_integral_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(楕円関数の境界積分は消える——周期性だけから。★無条件)",
    sectionId := "genell-lemma-3-5" }

def elliptic_boundary_integral_zero.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残る半分は留数定理(∮ f = 2πi·Σ res)である。" ++
       "★合わせて Σ_D res(f) = 0 が出て、f = g′/g に当てると #零点 = #極(偏角の原理)、" ++
       "g = ℘ − c なら「℘ は各値をちょうど 2 回取る」が従う" ++
       "——Skeleton/GenEll/AdditionTheorem.lean(第 615)で同定した唯一の入口である") 21,
    .implicitStep
      ("☆mathlib の在庫(2026-08-29): 軸平行な長方形での Cauchy" ++
       "(Complex.integral_boundary_rect_eq_zero_of_differentiable_on_off_countable)・" ++
       "Jensen の公式・Nevanlinna 理論(ValueDistribution/)・MeromorphicOn.divisor はある。" ++
       "☆偏角の原理と楕円関数の零点和 = 極和は無い。" ++
       "★一般の格子の平行四辺形には長方形版の Cauchy は当たらない" ++
       "(変数変換が ℝ-線型で正則性を壊す)") 13 ]

/-- ★★★★★★★★**`2f + a ∈ Λ` なら `℘(f+a) = ℘(f)`**——`℘` が偶であることから。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★`G(z) ≔ ℘(z+a) − ℘(z)` は `z ↦ −a−z` について奇であり、その**不動点**
（`2f ≡ −a (mod Λ)` を満たす `f`——`Λ/2Λ` の分だけ **4 個**ある）で消える。
☆これが「`℘(z+a) = ℘(z)` の解は `a ∉ Λ` でも 4 個ある」という事実であり、
★★★零点勘定（`#零点 = #極`）と組み合わせると
**`℘` が各値を 2 回しか取らない**ことの矛盾を作る材料になる
（`G` の極は `Λ` と `−a+Λ` に 2 位ずつ＝計 4）。 -/
theorem weierstrassP_shift_eq_of_two_add_mem (P : PeriodPair) (f a : ℂ)
    (h : 2 * f + a ∈ P.lattice) :
    P.weierstrassP (f + a) = P.weierstrassP f := by
  have h1 : P.weierstrassP (f + a) = P.weierstrassP (-(f + a)) := (P.weierstrassP_neg _).symm
  have h2 : -(f + a) = f + (-(2 * f + a)) := by ring
  rw [h1, h2, P.weierstrassP_add_coe f ⟨-(2 * f + a), neg_mem h⟩]

open Filter Topology Bornology Metric in
/-- ★★★★★★★★★★★★**`℘` の周期群はちょうど `Λ`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★機構: `℘(z+a) = ℘(z)` が恒等的に成り立つなら、`℘` は `−a` でも極を持つ
（`℘(z+a)` が `z = −a` で極だから）。★★しかし `−a ∉ Λ` なら `℘` は `−a` で解析的
——連続なので有界な極限を持ち、`cobounded` へ発散することと両立しない。

☆これで `G(z) ≔ ℘(z+a) − ℘(z)` が `a ∉ Λ` のとき恒等的に `0` でないことが分かる
——零点勘定の議論で「`G ≢ 0`」を言うのに要る。 -/
theorem mem_lattice_of_weierstrassP_periodic (P : PeriodPair) (a : ℂ)
    (h : ∀ z, P.weierstrassP (z + a) = P.weierstrassP z) : a ∈ P.lattice := by
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
  have h3 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (cobounded ℂ) := by
    have h2 := h1.comp hshift
    have hfun : (P.weierstrassP ∘ fun z : ℂ => z + a) = P.weierstrassP := by
      funext z; exact h z
    rwa [hfun] at h2
  have h4 : Tendsto P.weierstrassP (𝓝[≠] (-a)) (𝓝 (P.weierstrassP (-a))) :=
    hcont.continuousWithinAt
  exact (h4.not_tendsto (disjoint_nhds_cobounded _)) h3

def weierstrassP_shift_eq_of_two_add_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2f + a ∈ Λ なら ℘(f+a) = ℘(f)——℘ が偶だから。★無条件)",
    sectionId := "genell-lemma-3-5" }

def mem_lattice_of_weierstrassP_periodic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘ の周期群はちょうど Λ。★無条件)",
    sectionId := "genell-lemma-3-5" }

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

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★部分群の原像 -/

open scoped Classical in
/-- ★★★★★★`Φ(0) = 0`。 -/
@[simp] theorem uniformMap_zero (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    uniformMap P hΔ 0 = 0 :=
  uniformMap_of_mem P hΔ P.lattice.zero_mem

open scoped Classical in
/-- ★★★★★★★★`Φ(−z) = −Φ(z)`。 -/
theorem uniformMap_neg (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformMap P hΔ (-z) = -uniformMap P hΔ z := by
  have h := uniformMap_add P hΔ z (-z)
  rw [add_neg_cancel, uniformMap_zero] at h
  exact (neg_eq_of_add_eq_zero_right h.symm).symm

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★**部分群 `H ⊆ E(ℂ)` の原像は ℂ の部分群**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★`Φ` が群準同型（第 661）だから。☆これが `Lemma 3.5` の `Λ′` である。 -/
noncomputable def preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) : AddSubgroup ℂ where
  carrier := {z : ℂ | uniformMap P hΔ z ∈ H}
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq, uniformMap_add] at *
    exact H.add_mem ha hb
  zero_mem' := by
    simp only [Set.mem_setOf_eq, uniformMap_zero]
    exact H.zero_mem
  neg_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq, uniformMap_neg] at *
    exact H.neg_mem ha

open scoped Classical in
@[simp] theorem mem_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) (z : ℂ) :
    z ∈ preimageSubgroup P hΔ H ↔ uniformMap P hΔ z ∈ H := Iff.rfl

open scoped Classical in
/-- ★★★★★★★★★★**`Λ ⊆ Λ′`**——格子は `Φ` で `0` に落ちるから。 -/
theorem lattice_le_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) {z : ℂ} (hz : z ∈ P.lattice) :
    z ∈ preimageSubgroup P hΔ H := by
  simp only [mem_preimageSubgroup, uniformMap_of_mem P hΔ hz]
  exact H.zero_mem

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`Φ` は `Λ` を法として単射**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 624 の単射性を `Φ` の言葉に直したもの。☆`Λ′/Λ → H` の単射性に要る。 -/
theorem sub_mem_lattice_of_uniformMap_eq (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {z w : ℂ} (h : uniformMap P hΔ z = uniformMap P hΔ w) : z - w ∈ P.lattice := by
  by_cases hz : z ∈ P.lattice
  · by_cases hw : w ∈ P.lattice
    · exact P.lattice.sub_mem hz hw
    · rw [uniformMap_of_mem P hΔ hz, uniformMap_of_notMem P hΔ hw] at h
      exact absurd h.symm (by simp)
  · by_cases hw : w ∈ P.lattice
    · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_mem P hΔ hw] at h
      exact absurd h (by simp)
    · rw [uniformMap_of_notMem P hΔ hz, uniformMap_of_notMem P hΔ hw] at h
      have hx : latticePointX P z = latticePointX P w := by
        injection h with hx hy
      have hy : latticePointY P z = latticePointY P w := by
        injection h with hx hy
      have hpx : P.weierstrassP z = P.weierstrassP w := hx
      have hpy : P.derivWeierstrassP z = P.derivWeierstrassP w := by
        simp only [latticePointY] at hy
        linear_combination 2 * hy
      refine mem_lattice_of_shift_eq P (z - w) hw ?_ ?_ ?_
      · rw [show w + (z - w) = z by ring]; exact hz
      · rw [show w + (z - w) = z by ring]; exact hpx
      · rw [show w + (z - w) = z by ring]; exact hpy

def preimageSubgroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(部分群 H ⊆ E(ℂ) の原像は ℂ の部分群——これが Λ′ である)",
    sectionId := "genell-lemma-3-5" }

def sub_mem_lattice_of_uniformMap_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は Λ を法として単射。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★群同型 `ℂ/Λ ≅ E(ℂ)` -/

open scoped Classical in
/-- ★★★★★★★★★★`Φ(z) = 0 ⟺ z ∈ Λ`——`Φ` の核はちょうど `Λ`。 -/
theorem uniformMap_eq_zero_iff (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformMap P hΔ z = 0 ↔ z ∈ P.lattice := by
  refine ⟨fun h => ?_, fun h => uniformMap_of_mem P hΔ h⟩
  by_contra hz
  rw [uniformMap_of_notMem P hΔ hz] at h
  exact absurd h (by simp)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**一様化写像を加法群準同型として束ねたもの**。

★中身は第 661 の `uniformMap_add`。 -/
noncomputable def uniformHom (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    ℂ →+ (latticeCurve P).toAffine.Point :=
  AddMonoidHom.mk' (uniformMap P hΔ) (uniformMap_add P hΔ)

open scoped Classical in
@[simp] theorem uniformHom_apply (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformHom P hΔ z = uniformMap P hΔ z := rfl

open scoped Classical in
/-- ★★★★★★★★★★★★`Φ` の核は `Λ`。 -/
theorem ker_uniformHom (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    (uniformHom P hΔ).ker = P.lattice.toAddSubgroup := by
  ext z
  simp only [AddMonoidHom.mem_ker, uniformHom_apply, uniformMap_eq_zero_iff,
    Submodule.mem_toAddSubgroup]

open scoped Classical in
/-- ★★★★★★★★★★★★★★**`Φ` は全射**——第 604 の `latticePoint_surjective` を
`Point` の言葉に直したもの。 -/
theorem uniformMap_surjective (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    Function.Surjective (uniformMap P hΔ) := by
  intro Q
  cases Q with
  | zero => exact ⟨0, uniformMap_zero P hΔ⟩
  | some x y h =>
      obtain ⟨z, hz, hx, hy⟩ := latticePoint_surjective P x y h.left
      refine ⟨z, ?_⟩
      rw [uniformMap_of_notMem P hΔ hz]
      exact point_some_congr hx hy

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**一様化定理——`ℂ/Λ ≅ E(ℂ)` は加法群の同型**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★三つの部品がすべて揃った:

* **全射** —— 第 603（`weierstrassP_surjective`）＋ 第 604
* **単射** —— 第 624（`mem_lattice_of_shift_eq`）＋ 第 662
* **準同型** —— 第 661（`uniformMap_add`）

★★★★★☆**どの部品も mathlib に無い**（`§9-1039` で測った通り）。
☆これで `Lemma 3.5` の「`l`-捻れの部分群 ↔ `Λ` を含む指数 `l` の格子」が
純粋に群論の言葉で書けるようになる。 -/
noncomputable def uniformEquiv (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) :
    (ℂ ⧸ P.lattice.toAddSubgroup) ≃+ (latticeCurve P).toAffine.Point :=
  AddEquiv.ofBijective
    (QuotientAddGroup.lift P.lattice.toAddSubgroup (uniformHom P hΔ)
      (fun z hz => uniformMap_of_mem P hΔ hz))
    ⟨by
      intro a b hab
      obtain ⟨z, rfl⟩ := QuotientAddGroup.mk_surjective a
      obtain ⟨w, rfl⟩ := QuotientAddGroup.mk_surjective b
      simp only [QuotientAddGroup.lift_mk, uniformHom_apply] at hab
      rw [QuotientAddGroup.eq, Submodule.mem_toAddSubgroup,
        show -z + w = -(z - w) by ring]
      exact P.lattice.neg_mem (sub_mem_lattice_of_uniformMap_eq P hΔ hab),
     by
      intro Q
      obtain ⟨z, hz⟩ := uniformMap_surjective P hΔ Q
      exact ⟨(z : ℂ ⧸ P.lattice.toAddSubgroup), by simpa using hz⟩⟩

open scoped Classical in
@[simp] theorem uniformEquiv_mk (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ) :
    uniformEquiv P hΔ (z : ℂ ⧸ P.lattice.toAddSubgroup) = uniformMap P hΔ z := rfl

def uniformEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(一様化定理——ℂ/Λ ≅ E(ℂ) は加法群の同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def uniformMap_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ は全射。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★階数 1 の部分群（`Lemma 3.5`） -/

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**位数 `l` の点 `Q` を生む `z₀`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Φ` は全射（第 663）だから `Φ(z₀) = Q` なる `z₀` が取れ、`Q ≠ 0` だから
`z₀ ∉ Λ`、`l·Q = 0` だから `Φ(l z₀) = 0`、すなわち `l z₀ ∈ Λ`。
☆この `z₀` が「`Λ` に `1/l` の分母で足す元」である。 -/
theorem exists_generator_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} (hQ : Q ≠ 0) {l : ℕ} (hl : l • Q = 0) :
    ∃ z₀ : ℂ, z₀ ∉ P.lattice ∧ uniformMap P hΔ z₀ = Q ∧ (l : ℂ) * z₀ ∈ P.lattice := by
  obtain ⟨z₀, hz₀⟩ := uniformMap_surjective P hΔ Q
  refine ⟨z₀, fun hc => hQ ?_, hz₀, ?_⟩
  · rw [← hz₀, uniformMap_of_mem P hΔ hc]
  · have hzero : uniformMap P hΔ ((l : ℂ) * z₀) = 0 := by
      have h1 : ((l : ℂ) * z₀) = l • z₀ := by simp [nsmul_eq_mul]
      rw [h1, ← uniformHom_apply, map_nsmul, uniformHom_apply, hz₀, hl]
    exact (uniformMap_eq_zero_iff P hΔ _).1 hzero

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**巡回部分群の原像は `Λ′ = Λ + ℤ z₀`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`⊇` は準同型性（第 661）から、`⊆` は単射性（第 662）から:
`Φ(z) = k·Φ(z₀) = Φ(k z₀)` なら `z − k z₀ ∈ Λ`。

★★★★☆**これが `Lemma 3.5` の「階数 1」の内容である**——`E(ℂ)` の巡回部分群
`⟨Q⟩` は `Λ` に元を 1 つ添加した格子に対応する。 -/
theorem preimageSubgroup_zmultiples (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z₀ : ℂ) :
    preimageSubgroup P hΔ (AddSubgroup.zmultiples (uniformMap P hΔ z₀))
      = P.lattice.toAddSubgroup ⊔ AddSubgroup.zmultiples z₀ := by
  ext z
  simp only [mem_preimageSubgroup, AddSubgroup.mem_zmultiples_iff, AddSubgroup.mem_sup,
    Submodule.mem_toAddSubgroup]
  constructor
  · rintro ⟨k, hk⟩
    have hzk : uniformMap P hΔ (k • z₀) = uniformMap P hΔ z := by
      rw [← uniformHom_apply, map_zsmul, uniformHom_apply, hk]
    exact ⟨z - k • z₀, sub_mem_lattice_of_uniformMap_eq P hΔ hzk.symm, k • z₀,
      ⟨k, rfl⟩, by abel⟩
  · rintro ⟨y, hy, w, ⟨k, rfl⟩, rfl⟩
    refine ⟨k, ?_⟩
    have h2 : uniformMap P hΔ (y + k • z₀) = k • uniformMap P hΔ z₀ := by
      rw [← uniformHom_apply, map_add, map_zsmul, uniformHom_apply, uniformHom_apply,
        uniformMap_of_mem P hΔ hy, zero_add]
    exact h2.symm

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**`Λ ⊆ Λ′` かつ `l·Λ′ ⊆ Λ`**——`Λ′` が `Λ` の「`1/l` 倍の中」にある。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆これで `Λ′` は `Λ` と `(1/l)Λ` に挟まれ、有限指数の格子であることが確定する
（`Λ ⊆ Λ′ ⊆ (1/l)Λ`、`[(1/l)Λ : Λ] = l²`）。 -/
theorem smul_preimageSubgroup_le (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {H : AddSubgroup (latticeCurve P).toAffine.Point} {l : ℕ}
    (hH : ∀ Q ∈ H, l • Q = 0) {z : ℂ} (hz : z ∈ preimageSubgroup P hΔ H) :
    (l : ℂ) * z ∈ P.lattice := by
  have hzero : uniformMap P hΔ ((l : ℂ) * z) = 0 := by
    have h1 : ((l : ℂ) * z) = l • z := by simp [nsmul_eq_mul]
    rw [h1, ← uniformHom_apply, map_nsmul, uniformHom_apply]
    exact hH _ ((mem_preimageSubgroup P hΔ H z).1 hz)
  exact (uniformMap_eq_zero_iff P hΔ _).1 hzero

def exists_generator_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数 l の点 Q を生む z₀ の存在。★無条件)",
    sectionId := "genell-lemma-3-5" }

def preimageSubgroup_zmultiples.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(巡回部分群の原像は Λ′ = Λ + ℤz₀——階数 1 の内容。★無条件)",
    sectionId := "genell-lemma-3-5" }

def smul_preimageSubgroup_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l·Λ′ ⊆ Λ——Λ′ は Λ と (1/l)Λ に挟まれる。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★指数 `[Λ′ : Λ] = |H|` -/

open scoped Classical in
/-- ★★★★★★★★★★`Φ` を `Λ′ → H` に制限したもの。 -/
noncomputable def preimageToH (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    preimageSubgroup P hΔ H →+ H where
  toFun z := ⟨uniformMap P hΔ z.1, (mem_preimageSubgroup P hΔ H z.1).1 z.2⟩
  map_zero' := Subtype.ext (by simpa using uniformMap_zero P hΔ)
  map_add' := fun x y => Subtype.ext (by simpa using uniformMap_add P hΔ x.1 y.1)

open scoped Classical in
@[simp] theorem preimageToH_coe (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) (z : preimageSubgroup P hΔ H) :
    (preimageToH P hΔ H z : (latticeCurve P).toAffine.Point) = uniformMap P hΔ z.1 := rfl

open scoped Classical in
/-- ★★★★★★★★★★★★`Λ′ → H` は全射——`Φ` が全射（第 663）だから。 -/
theorem preimageToH_surjective (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    Function.Surjective (preimageToH P hΔ H) := by
  rintro ⟨Q, hQ⟩
  obtain ⟨z, hz⟩ := uniformMap_surjective P hΔ Q
  refine ⟨⟨z, ?_⟩, Subtype.ext ?_⟩
  · rw [mem_preimageSubgroup, hz]; exact hQ
  · exact hz

open scoped Classical in
/-- ★★★★★★★★★★★★`Λ′ → H` の核は `Λ`（`Λ′` の中で見たもの）。 -/
theorem ker_preimageToH (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    (preimageToH P hΔ H).ker
      = P.lattice.toAddSubgroup.addSubgroupOf (preimageSubgroup P hΔ H) := by
  ext z
  simp only [AddMonoidHom.mem_ker, AddSubgroup.mem_addSubgroupOf, Submodule.mem_toAddSubgroup]
  constructor
  · intro h
    exact (uniformMap_eq_zero_iff P hΔ _).1 (congrArg Subtype.val h)
  · intro h
    exact Subtype.ext (by simpa using uniformMap_of_mem P hΔ h)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Λ′/Λ ≅ H`**——原像の商はもとの部分群と同型。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★全射（第 663）＋核 `= Λ`（本ブロック）＋第一同型定理。 -/
noncomputable def preimageQuotEquiv (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    (preimageSubgroup P hΔ H
        ⧸ P.lattice.toAddSubgroup.addSubgroupOf (preimageSubgroup P hΔ H)) ≃+ H :=
  (QuotientAddGroup.quotientAddEquivOfEq (ker_preimageToH P hΔ H).symm).trans
    (QuotientAddGroup.quotientKerEquivOfSurjective (preimageToH P hΔ H)
      (preimageToH_surjective P hΔ H))

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`[Λ′ : Λ] = |H|`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★☆**これで `Lemma 3.5` の「位数 `l` の部分群 ↔ 指数 `l` の格子」の
指数の部分が塞がった**——`H` が位数 `l` なら `Λ ⊆ Λ′` は指数 `l`。
☆残るのは `Λ′` の基底を取り出して行列式 `= l` を書くこと
（`Λ ⊆ Λ′ ⊆ (1/l)Λ` は第 664 で押さえてある）。 -/
theorem relIndex_preimageSubgroup (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (H : AddSubgroup (latticeCurve P).toAffine.Point) :
    P.lattice.toAddSubgroup.relIndex (preimageSubgroup P hΔ H) = Nat.card H := by
  rw [AddSubgroup.relIndex, AddSubgroup.index]
  exact Nat.card_congr (preimageQuotEquiv P hΔ H).toEquiv

def preimageQuotEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ′/Λ ≅ H——原像の商はもとの部分群と同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def relIndex_preimageSubgroup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5([Λ′ : Λ] = |H|——指数は部分群の位数。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`Λ′` の基底と行列式 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**指数 `l` の格子 `Λ′ = Λ + ℤz₀` の基底と行列式**

`l·z₀ = a·ω₁ + b·ω₂`（`gcd(a, b, l) = 1`）のとき、`Λ′` の基底 `ω₁′, ω₂′` と
整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★構成は初等的である。`h = gcd(a,b)`、`a = h a₁`、`b = h b₁` として
`a₁p + b₁q = 1` を取り

    η₁ ≔ a₁ω₁ + b₁ω₂,   η₂ ≔ −qω₁ + pω₂

とすると `(η₁, η₂)` は `Λ` の基底（行列式 `a₁p + b₁q = 1`）で `l z₀ = h η₁`。
`gcd(h, l) = 1`（`gcd(a,b,l) = 1` から）なので `xh + yl = 1` が取れ、

    ω₁′ ≔ η₁/l,   ω₂′ ≔ η₂

とすれば `z₀ = h·ω₁′`・`ω₁′ = x·z₀ + y·η₁` となって
`Λ′ = ℤω₁′ + ℤω₂′`。行列式は

    (pl)·a₁ − (−b₁)·(ql) = l·(pa₁ + b₁q) = l

★★★★☆**これで `Lemma 3.5` の格子側——「位数 `l` の巡回部分群 ↔
指数 `l` の格子」——が完全に閉じた。**
☆`htFalt_isogeny_le_of_analytic_minimal`（第 617）が要求する
`h₁`・`h₂`・`hdet` はこの `A, B, C, D` そのものである。 -/
theorem exists_lattice_basis_of_cyclic (P : PeriodPair) (z₀ : ℂ) (l : ℕ) (hl : 0 < l)
    (a b : ℤ) (hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂)
    (hgcd : Nat.gcd (Int.gcd a b) l = 1) :
    ∃ (ω₁' ω₂' : ℂ) (A B C D : ℤ),
      P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂' ∧
      P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂' ∧
      (A * D - B * C).natAbs = l ∧
      Submodule.span ℤ ({ω₁', ω₂'} : Set ℂ)
        = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  have hlC : (l : ℂ) ≠ 0 := Nat.cast_ne_zero.2 hl.ne'
  by_cases hab : Int.gcd a b = 0
  · -- `a = b = 0` の退化した場合。`l = 1`・`z₀ = 0`。
    have ha : a = 0 := (Int.gcd_eq_zero_iff.1 hab).1
    have hb : b = 0 := (Int.gcd_eq_zero_iff.1 hab).2
    have hl1 : l = 1 := by simpa [hab] using hgcd
    have hz0 : z₀ = 0 := by
      have := hz
      rw [ha, hb] at this
      simp only [Int.cast_zero, zero_mul, add_zero] at this
      exact (mul_eq_zero.1 this).resolve_left hlC
    refine ⟨P.ω₁, P.ω₂, 1, 0, 0, 1, by push_cast; ring, by push_cast; ring, by
      simp [hl1], ?_⟩
    rw [hz0, Submodule.span_zero_singleton, sup_bot_eq]
    rfl
  · -- 本体。`h = gcd(a,b) ≠ 0`。
    set h : ℤ := (Int.gcd a b : ℤ) with hh
    have hhpos : 0 < h := by
      simpa [hh] using Nat.pos_of_ne_zero hab
    have hhne : h ≠ 0 := hhpos.ne'
    set a₁ : ℤ := a / h with ha₁
    set b₁ : ℤ := b / h with hb₁
    have hae : a = h * a₁ := by
      rw [ha₁, Int.mul_ediv_cancel' (Int.gcd_dvd_left a b)]
    have hbe : b = h * b₁ := by
      rw [hb₁, Int.mul_ediv_cancel' (Int.gcd_dvd_right a b)]
    have hg1 : Int.gcd a₁ b₁ = 1 := by
      rw [ha₁, hb₁, hh]
      exact Int.gcd_div_gcd_div_gcd (Nat.pos_of_ne_zero hab)
    -- `a₁ p + b₁ q = 1`
    set p : ℤ := Int.gcdA a₁ b₁ with hp
    set q : ℤ := Int.gcdB a₁ b₁ with hq
    have hbez1 : a₁ * p + b₁ * q = 1 := by
      have := Int.gcd_eq_gcd_ab a₁ b₁
      rw [hg1] at this
      simpa [hp, hq] using this.symm
    -- `gcd(h, l) = 1` から `x h + y l = 1`
    have hghl : Int.gcd h (l : ℤ) = 1 := by
      simpa [hh] using hgcd
    set x : ℤ := Int.gcdA h (l : ℤ) with hx
    set y : ℤ := Int.gcdB h (l : ℤ) with hy
    have hbez2 : h * x + (l : ℤ) * y = 1 := by
      have := Int.gcd_eq_gcd_ab h (l : ℤ)
      rw [hghl] at this
      simpa [hx, hy] using this.symm
    -- 新しい基底
    set η₁ : ℂ := (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ with hη₁
    set η₂ : ℂ := (-q : ℤ) * P.ω₁ + (p : ℂ) * P.ω₂ with hη₂
    set ω₁' : ℂ := η₁ / (l : ℂ) with hω₁'
    have hlω : (l : ℂ) * ω₁' = η₁ := by
      rw [hω₁']; field_simp
    clear_value ω₁'
    have hlz : (l : ℂ) * z₀ = (h : ℂ) * η₁ := by
      rw [hz, hη₁, hae, hbe]
      push_cast
      ring
    have hz0 : z₀ = (h : ℂ) * ω₁' := by
      have : (l : ℂ) * z₀ = (l : ℂ) * ((h : ℂ) * ω₁') := by
        rw [hlz, ← hlω]; ring
      exact mul_left_cancel₀ hlC this
    have hbez1C : (a₁ : ℂ) * (p : ℂ) + (b₁ : ℂ) * (q : ℂ) = 1 := by
      exact_mod_cast congrArg (fun n : ℤ => (n : ℂ)) hbez1
    have hbez2C : (h : ℂ) * (x : ℂ) + (l : ℂ) * (y : ℂ) = 1 := by
      exact_mod_cast congrArg (fun n : ℤ => (n : ℂ)) hbez2
    have hω1 : P.ω₁ = ((p * l : ℤ) : ℂ) * ω₁' + ((-b₁ : ℤ) : ℂ) * η₂ := by
      rw [hη₂]
      push_cast
      have hlω' : (l : ℂ) * ω₁' = (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ := by
        rw [hlω, hη₁]
      linear_combination (-(p : ℂ)) * hlω' - P.ω₁ * hbez1C
    have hω2 : P.ω₂ = ((q * l : ℤ) : ℂ) * ω₁' + ((a₁ : ℤ) : ℂ) * η₂ := by
      rw [hη₂]
      push_cast
      have hlω' : (l : ℂ) * ω₁' = (a₁ : ℂ) * P.ω₁ + (b₁ : ℂ) * P.ω₂ := by
        rw [hlω, hη₁]
      linear_combination (-(q : ℂ)) * hlω' - P.ω₂ * hbez1C
    refine ⟨ω₁', η₂, p * l, -b₁, q * l, a₁, hω1, hω2, ?_, ?_⟩
    · have : (p * l) * a₁ - (-b₁) * (q * l) = (l : ℤ) * (a₁ * p + b₁ * q) := by ring
      rw [this, hbez1, mul_one, Int.natAbs_natCast]
    · -- `span {ω₁′, η₂} = Λ ⊔ span {z₀}`
      have hη₁mem : η₁ ∈ P.lattice := by
        rw [PeriodPair.mem_lattice]
        exact ⟨a₁, b₁, hη₁.symm⟩
      have hη₂mem : η₂ ∈ P.lattice := by
        rw [PeriodPair.mem_lattice]
        exact ⟨-q, p, hη₂.symm⟩
      refine le_antisymm ?_ ?_
      · -- `ω₁′ = x·z₀ + y·η₁`
        have hval : ω₁' = (x : ℂ) * z₀ + (y : ℂ) * η₁ := by
          rw [hz0, ← hlω]
          linear_combination (-ω₁') * hbez2C
        have hm1 : ω₁' ∈ P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
          rw [hval]
          refine Submodule.add_mem _ ?_ ?_
          · exact Submodule.mem_sup_right (by
              rw [show (x : ℂ) * z₀ = x • z₀ by simp [zsmul_eq_mul]]
              exact Submodule.smul_mem _ _ (Submodule.mem_span_singleton_self z₀))
          · exact Submodule.mem_sup_left (by
              rw [show (y : ℂ) * η₁ = y • η₁ by simp [zsmul_eq_mul]]
              exact Submodule.smul_mem _ _ hη₁mem)
        have hm2 : η₂ ∈ P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) :=
          Submodule.mem_sup_left hη₂mem
        refine Submodule.span_le.2 ?_
        intro w hw
        simp only [Set.mem_insert_iff, Set.mem_singleton_iff] at hw
        rcases hw with h1 | h1
        · rw [SetLike.mem_coe, h1]; exact hm1
        · rw [SetLike.mem_coe, h1]; exact hm2
      · refine sup_le ?_ ?_
        · intro w hw
          obtain ⟨m, n, rfl⟩ := PeriodPair.mem_lattice.1 hw
          have h1 : P.ω₁ ∈ Submodule.span ℤ ({ω₁', η₂} : Set ℂ) := by
            rw [Submodule.mem_span_pair]
            exact ⟨p * l, -b₁, by rw [hω1]; push_cast; simp [zsmul_eq_mul]⟩
          have h2 : P.ω₂ ∈ Submodule.span ℤ ({ω₁', η₂} : Set ℂ) := by
            rw [Submodule.mem_span_pair]
            exact ⟨q * l, a₁, by rw [hω2]; push_cast; simp [zsmul_eq_mul]⟩
          refine Submodule.add_mem _ ?_ ?_
          · rw [show (m : ℂ) * P.ω₁ = m • P.ω₁ by simp [zsmul_eq_mul]]
            exact Submodule.smul_mem _ _ h1
          · rw [show (n : ℂ) * P.ω₂ = n • P.ω₂ by simp [zsmul_eq_mul]]
            exact Submodule.smul_mem _ _ h2
        · rw [Submodule.span_le, Set.singleton_subset_iff, SetLike.mem_coe, hz0]
          rw [show (h : ℂ) * ω₁' = h • ω₁' by simp [zsmul_eq_mul]]
          exact Submodule.smul_mem _ _ (Submodule.subset_span (by simp))

def exists_lattice_basis_of_cyclic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(指数 l の格子 Λ′ = Λ + ℤz₀ の基底と行列式 = l。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★**位数がちょうど `l` なら `gcd(a, b, l) = 1`**。

`l z₀ = a ω₁ + b ω₂` のとき、`g ≔ gcd(a, b, l)` とすると
`(l/g)·z₀ = (a/g)ω₁ + (b/g)ω₂ ∈ Λ`、すなわち `(l/g)·Q = 0`。
`Q` の位数が `l` なら `l ∣ l/g` なので `g = 1`。 -/
theorem gcd_eq_one_of_addOrderOf (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) {a b : ℤ}
    (hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂) :
    Nat.gcd (Int.gcd a b) l = 1 := by
  set g : ℕ := Nat.gcd (Int.gcd a b) l with hg
  have hgpos : 0 < g := Nat.gcd_pos_of_pos_right _ hl
  have hgC : (g : ℂ) ≠ 0 := Nat.cast_ne_zero.2 hgpos.ne'
  obtain ⟨n, hn⟩ := (Nat.gcd_dvd_right (Int.gcd a b) l : g ∣ l)
  have hnpos : 0 < n := by
    rcases Nat.eq_zero_or_pos n with hc | hc
    · rw [hc, Nat.mul_zero] at hn; omega
    · exact hc
  have hga : (g : ℤ) ∣ a :=
    dvd_trans (Int.natCast_dvd_natCast.2 (Nat.gcd_dvd_left _ _)) (Int.gcd_dvd_left a b)
  have hgb : (g : ℤ) ∣ b :=
    dvd_trans (Int.natCast_dvd_natCast.2 (Nat.gcd_dvd_left _ _)) (Int.gcd_dvd_right a b)
  have hA : ((a / (g : ℤ) : ℤ) : ℂ) * (g : ℂ) = (a : ℂ) := by
    exact_mod_cast congrArg (fun t : ℤ => (t : ℂ)) (Int.ediv_mul_cancel hga)
  have hB : ((b / (g : ℤ) : ℤ) : ℂ) * (g : ℂ) = (b : ℂ) := by
    exact_mod_cast congrArg (fun t : ℤ => (t : ℂ)) (Int.ediv_mul_cancel hgb)
  have hnz : (n : ℂ) * z₀
      = ((a / (g : ℤ) : ℤ) : ℂ) * P.ω₁ + ((b / (g : ℤ) : ℤ) : ℂ) * P.ω₂ := by
    refine mul_left_cancel₀ hgC ?_
    calc (g : ℂ) * ((n : ℂ) * z₀) = (l : ℂ) * z₀ := by
          rw [hn]; push_cast; ring
      _ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂ := hz
      _ = (g : ℂ) * (((a / (g : ℤ) : ℤ) : ℂ) * P.ω₁
            + ((b / (g : ℤ) : ℤ) : ℂ) * P.ω₂) := by
          rw [← hA, ← hB]; ring
  have hnzmem : (n : ℂ) * z₀ ∈ P.lattice :=
    PeriodPair.mem_lattice.2 ⟨a / (g : ℤ), b / (g : ℤ), hnz.symm⟩
  have hnQ : n • Q = 0 := by
    have hzero : uniformMap P hΔ ((n : ℂ) * z₀) = 0 := uniformMap_of_mem P hΔ hnzmem
    rw [show ((n : ℂ) * z₀) = n • z₀ by simp [nsmul_eq_mul], ← uniformHom_apply,
      map_nsmul, uniformHom_apply, hz₀] at hzero
    exact hzero
  have hdvd : l ∣ n := hQ ▸ addOrderOf_dvd_of_nsmul_eq_zero hnQ
  have hle : l ≤ n := Nat.le_of_dvd hnpos hdvd
  have hne : n ≤ l := by
    rw [hn]
    exact Nat.le_mul_of_pos_left n hgpos
  have hnl : n = l := le_antisymm hne hle
  have hn' : l = g * l := by
    conv_lhs => rw [hn]
    rw [hnl]
  refine Nat.eq_of_mul_eq_mul_right hl ?_
  rw [← hn', Nat.one_mul]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（格子側）——位数 `l` の点から指数 `l` の格子を作る**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`（`Φ(z₀) = Q`）と
`Λ′ = Λ + ℤz₀` の基底 `ω₁′, ω₂′` と整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

★★★★★☆**これが `htFalt_isogeny_le_of_analytic_minimal`（第 617）の
`h₁`・`h₂`・`hdet` そのものである。**

★部品:

* 全射（第 663）で `z₀` を取る
* `l·Q = 0` と核 `= Λ`（第 663）で `l z₀ ∈ Λ`、すなわち `l z₀ = aω₁ + bω₂`
* 位数がちょうど `l` だから `gcd(a, b, l) = 1`（本ブロック）
* Hermite 標準形（第 666）で基底と行列式 -/
theorem exists_isogeny_lattice_basis (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ ω₁' ω₂' : ℂ) (A B C D : ℤ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂' ∧
      P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂' ∧
      (A * D - B * C).natAbs = l ∧
      Submodule.span ℤ ({ω₁', ω₂'} : Set ℂ)
        = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  obtain ⟨z₀, hz₀⟩ := uniformMap_surjective P hΔ Q
  have hlQ : l • Q = 0 := hQ ▸ addOrderOf_nsmul_eq_zero Q
  have hlz : ((l : ℂ) * z₀) ∈ P.lattice := by
    refine (uniformMap_eq_zero_iff P hΔ _).1 ?_
    rw [show ((l : ℂ) * z₀) = l • z₀ by simp [nsmul_eq_mul], ← uniformHom_apply,
      map_nsmul, uniformHom_apply, hz₀, hlQ]
  obtain ⟨a, b, hab⟩ := PeriodPair.mem_lattice.1 hlz
  have hz : (l : ℂ) * z₀ = (a : ℂ) * P.ω₁ + (b : ℂ) * P.ω₂ := hab.symm
  obtain ⟨ω₁', ω₂', A, B, C, D, h1, h2, hdet, hspan⟩ :=
    exists_lattice_basis_of_cyclic P z₀ l hl a b hz
      (gcd_eq_one_of_addOrderOf P hΔ hl hQ hz₀ hz)
  exact ⟨z₀, ω₁', ω₂', A, B, C, D, hz₀, h1, h2, hdet, hspan⟩

def gcd_eq_one_of_addOrderOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(位数がちょうど l なら gcd(a, b, l) = 1。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isogeny_lattice_basis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子側——位数 l の点から指数 l の格子を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★★★★★★★★★★★★★**基底変換は ℝ-線型独立性を保つ**。

`ω₁ = Aω₁′ + Bω₂′`・`ω₂ = Cω₁′ + Dω₂′`・`AD − BC ≠ 0` なら
`ω₁′, ω₂′` も ℝ 上独立。

★逆行列は `Δω₁′ = Dω₁ − Bω₂`・`Δω₂′ = −Cω₁ + Aω₂`（`Δ = AD − BC`）。
`rω₁′ + tω₂′ = 0` に `Δ` を掛けて `ω₁, ω₂` の独立性に帰着させる。 -/
theorem linearIndependent_of_basis_change (P : PeriodPair) {ω₁' ω₂' : ℂ} {A B C D : ℤ}
    (h1 : P.ω₁ = (A : ℂ) * ω₁' + (B : ℂ) * ω₂')
    (h2 : P.ω₂ = (C : ℂ) * ω₁' + (D : ℂ) * ω₂')
    (hdet : A * D - B * C ≠ 0) :
    LinearIndependent ℝ ![ω₁', ω₂'] := by
  have hP := LinearIndependent.pair_iff.1 P.indep
  have hΔR : ((A : ℝ) * D - B * C) ≠ 0 := by
    have : ((A * D - B * C : ℤ) : ℝ) ≠ 0 := Int.cast_ne_zero.2 hdet
    push_cast at this
    exact this
  refine LinearIndependent.pair_iff.2 ?_
  intro r t hrt
  have hrt' : (r : ℂ) * ω₁' + (t : ℂ) * ω₂' = 0 := by
    simpa [Complex.real_smul] using hrt
  have hkey : (r * D - t * C : ℝ) • P.ω₁ + (-(r * B) + t * A : ℝ) • P.ω₂ = 0 := by
    simp only [Complex.real_smul]
    push_cast
    rw [h1, h2]
    linear_combination ((A : ℂ) * (D : ℂ) - (B : ℂ) * (C : ℂ)) * hrt'
  obtain ⟨e1, e2⟩ := hP _ _ hkey
  constructor
  · have hr : r * ((A : ℝ) * D - B * C) = 0 := by
      linear_combination (A : ℝ) * e1 + (C : ℝ) * e2
    exact (mul_eq_zero.1 hr).resolve_right hΔR
  · have ht : t * ((A : ℝ) * D - B * C) = 0 := by
      linear_combination (B : ℝ) * e1 + (D : ℝ) * e2
    exact (mul_eq_zero.1 ht).resolve_right hΔR

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（格子側・完成形）——位数 `l` の点から `PeriodPair` `P′` を作る**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`（`Φ(z₀) = Q`）と
**周期対** `P′`（格子は `Λ′ = Λ + ℤz₀`）と整数 `A, B, C, D` が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

★★★★★☆**これで `htFalt_isogeny_le_of_analytic_minimal`（第 617）が要求する
`P′`・`h₁`・`h₂`・`hdet` がすべて揃った。**
☆残るのは `α`（`u′ = α·u`）——代数的な同種写像 `E → E/H` と
この解析的な `Λ ⊆ Λ′` を突き合わせることである。 -/
theorem exists_isogeny_periodPair (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) := by
  obtain ⟨z₀, ω₁', ω₂', A, B, C, D, hz₀, h1, h2, hdet, hspan⟩ :=
    exists_isogeny_lattice_basis P hΔ hl hQ
  have hdet0 : A * D - B * C ≠ 0 := by
    intro hc
    rw [hc] at hdet
    simp at hdet
    omega
  exact ⟨z₀, ⟨ω₁', ω₂', linearIndependent_of_basis_change P h1 h2 hdet0⟩,
    A, B, C, D, hz₀, h1, h2, hdet, hspan⟩

def linearIndependent_of_basis_change.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(基底変換は ℝ-線型独立性を保つ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_isogeny_periodPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子側・完成形——位数 l の点から PeriodPair P′ を作る。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★巡回の場合の代表系 -/

/-- ★`|k| < l` で `l ∣ k` なら `k = 0`。 -/
theorem eq_zero_of_dvd_of_abs_lt {l k : ℤ} (hl : 0 < l) (h : l ∣ k)
    (h1 : -l < k) (h2 : k < l) : k = 0 := by
  obtain ⟨c, rfl⟩ := h
  rcases lt_trichotomy c 0 with hc | hc | hc
  · nlinarith
  · simp [hc]
  · nlinarith

open scoped Classical in
/-- ★★★★★★★★★★★★★★`k·z₀ ∈ Λ ⟺ l ∣ k`——`Q` の位数がちょうど `l` のとき。 -/
theorem intCast_mul_mem_lattice_iff (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) (k : ℤ) :
    (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k := by
  rw [← uniformMap_eq_zero_iff P hΔ,
    show ((k : ℂ) * z₀) = k • z₀ by simp [zsmul_eq_mul],
    ← uniformHom_apply, map_zsmul, uniformHom_apply, hz₀, ← hQ,
    addOrderOf_dvd_iff_zsmul_eq_zero]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**巡回の場合の代表系 `T = {0, z₀, 2z₀, …, (l−1)z₀}`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`Λ′ = Λ + ℤz₀` で `z₀` の `Λ` を法とした位数がちょうど `l` のとき、
`T` は `Λ′/Λ` の代表系である（第 601 `weierstrassP_eq_velu_of_rep` の仮説）。

★★★`p = y + n z₀ ∈ Λ′` に対し `w₀ ≔ ((−n) mod l)·z₀` を取れば
`p + w₀ ∈ Λ`。一意性は `|k − m| < l` かつ `l ∣ k − m` から。 -/
theorem exists_velu_rep (P P' : PeriodPair) (z₀ : ℂ) (l : ℕ) (hl : 0 < l)
    (hP' : P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ))
    (hord : ∀ k : ℤ, (k : ℂ) * z₀ ∈ P.lattice ↔ (l : ℤ) ∣ k) :
    ∃ T : Finset ℂ, (0 : ℂ) ∈ T ∧ T.card = l ∧ (∀ w ∈ T, w ∈ P'.lattice) ∧
      (∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
        ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) := by
  have hlZ : (0 : ℤ) < (l : ℤ) := by exact_mod_cast hl
  have hinj : ∀ k ∈ Finset.range l, ∀ k' ∈ Finset.range l,
      (k : ℂ) * z₀ = (k' : ℂ) * z₀ → k = k' := by
    intro k hk k' hk' he
    simp only [Finset.mem_range] at hk hk'
    have h0 : (((k : ℤ) - (k' : ℤ) : ℤ) : ℂ) * z₀ ∈ P.lattice := by
      have hz : (((k : ℤ) - (k' : ℤ) : ℤ) : ℂ) * z₀ = 0 := by
        push_cast
        linear_combination he
      rw [hz]
      exact P.lattice.zero_mem
    have hd := (hord _).1 h0
    have hzero := eq_zero_of_dvd_of_abs_lt hlZ hd (by omega) (by omega)
    omega
  refine ⟨(Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀), ?_, ?_, ?_, ?_⟩
  · exact Finset.mem_image.2 ⟨0, Finset.mem_range.2 hl, by simp⟩
  · have hio : Set.InjOn (fun k : ℕ => (k : ℂ) * z₀) ↑(Finset.range l) := by
      intro k hk k' hk' he
      exact hinj k (by simpa using hk) k' (by simpa using hk') he
    rw [Finset.card_image_of_injOn hio, Finset.card_range]
  · intro w hw
    obtain ⟨k, -, rfl⟩ := Finset.mem_image.1 hw
    rw [hP']
    refine Submodule.mem_sup_right ?_
    rw [show ((k : ℂ) * z₀) = (k : ℤ) • z₀ by simp [zsmul_eq_mul]]
    exact Submodule.smul_mem _ _ (Submodule.mem_span_singleton_self z₀)
  · intro p hp
    rw [hP', Submodule.mem_sup] at hp
    obtain ⟨y, hy, w, hw, rfl⟩ := hp
    obtain ⟨n, rfl⟩ := Submodule.mem_span_singleton.1 hw
    set r : ℤ := (-n) % (l : ℤ) with hrdef
    have hr0 : 0 ≤ r := Int.emod_nonneg _ hlZ.ne'
    have hrl : r < (l : ℤ) := Int.emod_lt_of_pos _ hlZ
    set m : ℕ := r.toNat with hmdef
    have hmr : (m : ℤ) = r := Int.toNat_of_nonneg hr0
    have hml : m < l := by omega
    have hdvd : (l : ℤ) ∣ (n + (m : ℤ)) := by
      rw [hmr]
      refine ⟨-((-n) / (l : ℤ)), ?_⟩
      have hq := Int.emod_add_ediv (-n) (l : ℤ)
      linear_combination hq
    refine ⟨(m : ℂ) * z₀, Finset.mem_image.2 ⟨m, Finset.mem_range.2 hml, rfl⟩, ?_, ?_⟩
    · have he : y + (n • z₀) + (m : ℂ) * z₀ = y + ((n + (m : ℤ) : ℤ) : ℂ) * z₀ := by
        push_cast [zsmul_eq_mul]
        ring
      rw [he]
      exact P.lattice.add_mem hy ((hord _).2 hdvd)
    · intro v hv hvne hmem
      obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hv
      rw [Finset.mem_range] at hk
      have he : y + (n • z₀) + (k : ℂ) * z₀ = y + ((n + (k : ℤ) : ℤ) : ℂ) * z₀ := by
        push_cast [zsmul_eq_mul]
        ring
      rw [he] at hmem
      have hk2 : ((n + (k : ℤ) : ℤ) : ℂ) * z₀ ∈ P.lattice := by
        have hd2 := P.lattice.sub_mem hmem hy
        simpa using hd2
      have hsub : (l : ℤ) ∣ ((k : ℤ) - (m : ℤ)) := by
        have h2 := (hord (n + (k : ℤ))).1 hk2
        have h3 : (l : ℤ) ∣ ((n + (k : ℤ)) - (n + (m : ℤ))) := dvd_sub h2 hdvd
        have h4 : (n + (k : ℤ)) - (n + (m : ℤ)) = (k : ℤ) - (m : ℤ) := by ring
        rwa [h4] at h3
      have hzero : ((k : ℤ) - (m : ℤ)) = 0 :=
        eq_zero_of_dvd_of_abs_lt hlZ hsub (by omega) (by omega)
      exact hvne (by rw [show k = m by omega])

def exists_velu_rep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(巡回の場合の代表系 T = {0, z₀, …, (l−1)z₀}。★無条件)",
    sectionId := "genell-lemma-3-5" }

def intCast_mul_mem_lattice_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(k·z₀ ∈ Λ ⟺ l ∣ k——位数がちょうど l のとき。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`Lemma 3.5`（解析側・完成形）——位数 `l` の点から Vélu の公式まで**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E(ℂ)` の位数ちょうど `l` の点 `Q` に対し、`z₀`・周期対 `P′`・整数 `A, B, C, D`・
代表系 `T`（`|T| = l`）が取れて

    ω₁ = A·ω₁′ + B·ω₂′,  ω₂ = C·ω₁′ + D·ω₂′,  |AD − BC| = l

    ℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)

★★★★★★☆**これで「位数 `l` の点 → 指数 `l` の格子 → その `℘` は
`Λ` の `℘` の Vélu 和で書ける」が一本につながった。**

☆残るのは、この解析的な等式を代数的な Vélu の商
（`Found/GenEll/Velu.lean` の `veluQuotient`）と突き合わせて
`α`（`u′ = α·u`）を作ることである。 -/
theorem exists_velu_formula_of_torsion (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    ∃ (z₀ : ℂ) (P' : PeriodPair) (A B C D : ℤ) (T : Finset ℂ),
      uniformMap P hΔ z₀ = Q ∧
      P.ω₁ = (A : ℂ) * P'.ω₁ + (B : ℂ) * P'.ω₂ ∧
      P.ω₂ = (C : ℂ) * P'.ω₁ + (D : ℂ) * P'.ω₂ ∧
      (A * D - B * C).natAbs = l ∧
      P'.lattice = P.lattice ⊔ Submodule.span ℤ ({z₀} : Set ℂ) ∧
      (0 : ℂ) ∈ T ∧ T.card = l ∧
      (∀ z : ℂ, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) := by
  obtain ⟨z₀, P', A, B, C, D, hz₀, h1, h2, hdet, hP'⟩ :=
    exists_isogeny_periodPair P hΔ hl hQ
  obtain ⟨T, h0T, hcard, hT, hrep⟩ :=
    exists_velu_rep P P' z₀ l hl hP' (intCast_mul_mem_lattice_iff P hΔ hQ hz₀)
  have hle : P.lattice ≤ P'.lattice := by rw [hP']; exact le_sup_left
  exact ⟨z₀, P', A, B, C, D, T, hz₀, h1, h2, hdet, hP', h0T, hcard,
    fun z => weierstrassP_eq_velu_of_rep P P' hle T h0T hT hrep z⟩

def exists_velu_formula_of_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(解析側・完成形——位数 l の点から Vélu の公式まで。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★代表系の対称性 -/

/-- ★★★★★★★★★★★★**代表系の元は `Λ` を法として相異なる**。 -/
theorem rep_sub_mem_lattice_imp_eq (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {w v : ℂ} (hw : w ∈ T) (hv : v ∈ T) (hd : w - v ∈ P.lattice) : w = v := by
  obtain ⟨w₀, hw₀T, hw₀Λ, hw₀u⟩ := hrep (-v) (neg_mem (hT v hv))
  have hv0 : v = w₀ := by
    by_contra hc
    exact hw₀u v hv hc (by simpa using P.lattice.zero_mem)
  have hw0 : w = w₀ := by
    by_contra hc
    exact hw₀u w hw hc (by rw [show -v + w = w - v by ring]; exact hd)
  rw [hw0, hv0]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系の上での `℘′` の和は消える**

    `Σ_{w ∈ T∖{0}} ℘′_Λ(w) = 0`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★`w ↦ ν(w)`（`w + ν w ∈ Λ` を満たす唯一の代表元、つまり `ν w ≡ −w`）は
`T∖{0}` の対合であり、`℘′(ν w) = ℘′(−w) = −℘′(w)`。
したがって `S = Σ ℘′(ν w) = −S`、すなわち `S = 0`。

★★☆これが Vélu の `ω`-正規化（第 586-593 の代数側 `velu_omega_gen`）の解析版である。 -/
theorem sum_derivWeierstrassP_rep_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    ∑ w ∈ T.erase 0, P.derivWeierstrassP w = 0 := by
  classical
  have huniq : ∀ {w v : ℂ}, w ∈ T → v ∈ T → w - v ∈ P.lattice → w = v :=
    fun hw hv hd => rep_sub_mem_lattice_imp_eq P P' T hT hrep hw hv hd
  have hex : ∀ w ∈ T, ∃ v, v ∈ T ∧ w + v ∈ P.lattice := by
    intro w hw
    obtain ⟨v, hv, hv2, -⟩ := hrep w (hT w hw)
    exact ⟨v, hv, hv2⟩
  choose! ν hνT hνΛ using hex
  have hνe : ∀ w ∈ T.erase 0, ν w ∈ T.erase 0 := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have hw0 : w ≠ 0 := Finset.ne_of_mem_erase hw
    refine Finset.mem_erase.2 ⟨?_, hνT w hw'⟩
    intro hc
    refine hw0 (huniq hw' h0T ?_)
    have hz := hνΛ w hw'
    rw [hc, add_zero] at hz
    simpa using hz
  have hinvol : ∀ w ∈ T.erase 0, ν (ν w) = w := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have h1 := hνΛ w hw'
    have h2 := hνΛ (ν w) (hνT w hw')
    refine huniq (hνT (ν w) (hνT w hw')) hw' ?_
    have hd := P.lattice.sub_mem h2 h1
    rw [show ν w + ν (ν w) - (w + ν w) = ν (ν w) - w by ring] at hd
    exact hd
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0, ν w = ν v → w = v := by
    intro w hw v hv he
    rw [← hinvol w hw, ← hinvol v hv, he]
  have hodd : ∀ w ∈ T, P.derivWeierstrassP (ν w) = -P.derivWeierstrassP w := by
    intro w hw
    have hl : w + ν w ∈ P.lattice := hνΛ w hw
    have he : ν w = -w + (w + ν w) := by ring
    rw [he, P.derivWeierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.derivWeierstrassP_neg]
  have hinjOn : Set.InjOn ν ↑(T.erase 0) := fun w hw v hv he =>
    hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he
  have himg : (T.erase 0).image ν = T.erase 0 :=
    Finset.eq_of_subset_of_card_le
      (fun v hv => by
        obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hv
        exact hνe w hw)
      (le_of_eq (Finset.card_image_of_injOn hinjOn).symm)
  have h1 : ∑ v ∈ T.erase 0, P.derivWeierstrassP v
      = ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w) := by
    conv_lhs => rw [← himg]
    exact Finset.sum_image (fun w hw v hv he => hinj w hw v hv he)
  have h2 : ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w)
      = -∑ w ∈ T.erase 0, P.derivWeierstrassP w := by
    have hc : ∑ w ∈ T.erase 0, P.derivWeierstrassP (ν w)
        = ∑ w ∈ T.erase 0, (-P.derivWeierstrassP w) :=
      Finset.sum_congr rfl (fun w hw => hodd w (Finset.mem_of_mem_erase hw))
    rw [hc, Finset.sum_neg_distrib]
  have hSS : ∑ w ∈ T.erase 0, P.derivWeierstrassP w
      = -∑ w ∈ T.erase 0, P.derivWeierstrassP w := h1.trans h2
  have h3 : (2 : ℂ) * ∑ w ∈ T.erase 0, P.derivWeierstrassP w = 0 := by
    linear_combination hSS
  simpa using h3

def rep_sub_mem_lattice_imp_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の元は Λ を法として相異なる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_derivWeierstrassP_rep_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の上での ℘′ の和は消える——Vélu の ω-正規化の解析版。★無条件)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★Laurent 比較の入口 -/

/-- ★★★★★★★★★★**代表系の `0` 以外の元は格子の外**——`w ≡ 0` なら `w = 0` だから。 -/
theorem rep_notMem_lattice (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {w : ℂ} (hw : w ∈ T.erase 0) : w ∉ P.lattice := by
  intro hc
  refine Finset.ne_of_mem_erase hw ?_
  refine rep_sub_mem_lattice_imp_eq P P' T hT hrep (Finset.mem_of_mem_erase hw) h0T ?_
  simpa using hc

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`℘_{Λ′} − ℘_Λ` は `Λ` の `℘` の平行移動の和**

    `℘_{Λ′}(z) − ℘_Λ(z) = Σ_{w ∈ T∖{0}} (℘_Λ(z + w) − ℘_Λ(w))`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 669 の Vélu の公式から `w = 0` の項を分離しただけ。
☆左辺は原点の `z⁻²` の極が打ち消し合っており、右辺は原点で解析的
（`T∖{0}` の元は格子の外）。**これが Laurent 係数の比較の入口である**——
`z²` の係数から `g₂′ − g₂`、`z⁴` の係数から `g₃′ − g₃` が出る。 -/
theorem weierstrassP_sub_eq_sum (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) (z : ℂ) :
    P'.weierstrassP z - P.weierstrassP z
      = ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) := by
  rw [hvelu z]
  simp only [veluAnalyticX, veluAnalyticC]
  rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (z + w)) h0T]
  simp only [add_zero]
  rw [Finset.sum_sub_distrib]
  ring

/-- ★★★★★★★★★★★★**平行移動の和は原点で解析的**。 -/
theorem analyticAt_veluShiftSum (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0 : ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => ∑ w ∈ T.erase 0,
      (P.weierstrassP (z + w) - P.weierstrassP w)) 0 := by
  have h : AnalyticAt ℂ (∑ w ∈ T.erase 0,
      fun z : ℂ => P.weierstrassP (z + w) - P.weierstrassP w) 0 := by
    refine Finset.analyticAt_sum _ ?_
    intro w hw
    refine AnalyticAt.sub ?_ analyticAt_const
    refine shifted_analyticAt P 0 w ?_
    rw [zero_add]
    exact rep_notMem_lattice P P' T h0T hT hrep hw
  refine h.congr ?_
  filter_upwards with z
  simp [Finset.sum_apply]

/-- ★★★★★★★★★★★★★★★★
**`℘[Λ′ − 0] − ℘[Λ − 0]` は平行移動の和に等しい**——`z⁻²` が消える形。

☆左辺は mathlib の `iteratedDeriv_weierstrassPExcept_self` で Taylor 係数が
`sumInvPow`（＝`G n`、したがって `g₂`・`g₃`）で書けている。 -/
theorem weierstrassPExcept_sub_eq_sum (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T)
    (hvelu : ∀ z, P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z) (z : ℂ) :
    P'.weierstrassPExcept 0 z - P.weierstrassPExcept 0 z
      = ∑ w ∈ T.erase 0, (P.weierstrassP (z + w) - P.weierstrassP w) := by
  rw [← weierstrassP_sub_invSq P' z, ← weierstrassP_sub_invSq P z,
    ← weierstrassP_sub_eq_sum P P' T h0T hvelu z]
  ring

def weierstrassP_sub_eq_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘_{Λ′} − ℘_Λ は ℘_Λ の平行移動の和——Laurent 比較の入口。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassPExcept_sub_eq_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘[Λ′−0] − ℘[Λ−0] は平行移動の和——z⁻² が消える形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def rep_notMem_lattice.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の 0 以外の元は格子の外。★無条件)",
    sectionId := "genell-lemma-3-5" }

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
  obtain ⟨T, h0T, hcard, hT, hrep⟩ :=
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
  obtain ⟨T, h0T, hcard, hT, hrep⟩ :=
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
  obtain ⟨T, h0T, hcard, hT, hrep⟩ :=
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

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★代表系の像は `⟨Q⟩` -/

open scoped Classical in
/-- ★★★★★★★★`Φ(k·z) = k • Φ(z)`（`k : ℕ`）。 -/
theorem uniformMap_natCast_mul (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z : ℂ)
    (k : ℕ) : uniformMap P hΔ ((k : ℂ) * z) = k • uniformMap P hΔ z := by
  rw [show ((k : ℂ) * z) = k • z by simp [nsmul_eq_mul], ← uniformHom_apply,
    map_nsmul, uniformHom_apply]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系 `T = {0, z₀, …, (l−1)z₀}` の像は `⟨Q⟩ = {0, Q, …, (l−1)Q}`**

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★☆**これが Galois 降下の鍵である**——解析側で選んだ代表系 `T` の像は
`Q` が生成する巡回部分群そのものなので、`veluQuotientFull` に渡す点集合 `S` は
**`Q` だけで決まる**（`z₀` の選び方に依らない）。
したがって `Q` が `L`-有理なら `S` も `L`-有理である。 -/
theorem image_uniformMap_veluRep (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) (z₀ : ℂ)
    (l : ℕ) :
    ((Finset.range l).image (fun k : ℕ => (k : ℂ) * z₀)).image (uniformMap P hΔ)
      = (Finset.range l).image (fun k : ℕ => k • uniformMap P hΔ z₀) := by
  rw [Finset.image_image]
  refine Finset.image_congr ?_
  intro k _
  exact uniformMap_natCast_mul P hΔ z₀ k

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★**`T∖{0}` の像は `⟨Q⟩∖{O}`**。

★`0 < k < l` で `k • Q ≠ 0`（`Q` の位数がちょうど `l` だから）。 -/
theorem uniformMap_ne_zero_of_mem_erase (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    {Q : (latticeCurve P).toAffine.Point} {l : ℕ} (hQ : addOrderOf Q = l)
    {z₀ : ℂ} (hz₀ : uniformMap P hΔ z₀ = Q) {k : ℕ} (hk0 : k ≠ 0) (hkl : k < l) :
    uniformMap P hΔ ((k : ℂ) * z₀) ≠ 0 := by
  rw [uniformMap_natCast_mul, hz₀]
  intro hc
  have hdvd : addOrderOf Q ∣ k := addOrderOf_dvd_of_nsmul_eq_zero hc
  rw [hQ] at hdvd
  have := Nat.le_of_dvd (Nat.pos_of_ne_zero hk0) hdvd
  omega

def image_uniformMap_veluRep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系 T の像は ⟨Q⟩——Galois 降下の鍵。★無条件)",
    sectionId := "genell-lemma-3-5" }

def uniformMap_ne_zero_of_mem_erase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(0 < k < l なら k·z₀ の像は O でない。★無条件)",
    sectionId := "genell-lemma-3-5" }

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**代表系の上で `f(℘)·℘′` の和は消える**

    `Σ_{w ∈ T∖{0}} f(℘_Λ(w))·℘′_Λ(w) = 0`   （任意の `f : ℂ → ℂ`）

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★第 670 の一般化。`ν w ≡ −w` は `T∖{0}` の対合で、
`℘` は偶（`℘(ν w) = ℘(w)`）・`℘′` は奇（`℘′(ν w) = −℘′(w)`）だから
`S = Σ f(℘(ν w))·℘′(ν w) = −S`。

★★☆`f = 1` なら第 670（`Σ ℘′ = 0`）、`f = id` なら `Σ ℘·℘′ = 0`。
後者は Vélu の `Σ B·x = 0`（第 688 の仮説）に要る。 -/
theorem sum_mul_derivWeierstrassP_rep_eq_zero (P P' : PeriodPair) (T : Finset ℂ)
    (h0T : (0 : ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) (f : ℂ → ℂ) :
    ∑ w ∈ T.erase 0, f (P.weierstrassP w) * P.derivWeierstrassP w = 0 := by
  classical
  have huniq : ∀ {w v : ℂ}, w ∈ T → v ∈ T → w - v ∈ P.lattice → w = v :=
    fun hw hv hd => rep_sub_mem_lattice_imp_eq P P' T hT hrep hw hv hd
  have hex : ∀ w ∈ T, ∃ v, v ∈ T ∧ w + v ∈ P.lattice := by
    intro w hw
    obtain ⟨v, hv, hv2, -⟩ := hrep w (hT w hw)
    exact ⟨v, hv, hv2⟩
  choose! ν hνT hνΛ using hex
  have hνe : ∀ w ∈ T.erase 0, ν w ∈ T.erase 0 := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have hw0 : w ≠ 0 := Finset.ne_of_mem_erase hw
    refine Finset.mem_erase.2 ⟨?_, hνT w hw'⟩
    intro hc
    refine hw0 (huniq hw' h0T ?_)
    have hz := hνΛ w hw'
    rw [hc, add_zero] at hz
    simpa using hz
  have hinvol : ∀ w ∈ T.erase 0, ν (ν w) = w := by
    intro w hw
    have hw' : w ∈ T := Finset.mem_of_mem_erase hw
    have h1 := hνΛ w hw'
    have h2 := hνΛ (ν w) (hνT w hw')
    refine huniq (hνT (ν w) (hνT w hw')) hw' ?_
    have hd := P.lattice.sub_mem h2 h1
    rw [show ν w + ν (ν w) - (w + ν w) = ν (ν w) - w by ring] at hd
    exact hd
  have hinj : ∀ w ∈ T.erase 0, ∀ v ∈ T.erase 0, ν w = ν v → w = v := by
    intro w hw v hv he
    rw [← hinvol w hw, ← hinvol v hv, he]
  have hodd : ∀ w ∈ T, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
      = -(f (P.weierstrassP w) * P.derivWeierstrassP w) := by
    intro w hw
    have hl : w + ν w ∈ P.lattice := hνΛ w hw
    have he : ν w = -w + (w + ν w) := by ring
    have hp : P.weierstrassP (ν w) = P.weierstrassP w := by
      rw [he, P.weierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.weierstrassP_neg]
    have hpd : P.derivWeierstrassP (ν w) = -P.derivWeierstrassP w := by
      rw [he, P.derivWeierstrassP_add_coe (-w) ⟨w + ν w, hl⟩, P.derivWeierstrassP_neg]
    rw [hp, hpd]
    ring
  have hinjOn : Set.InjOn ν ↑(T.erase 0) := fun w hw v hv he =>
    hinj w (Finset.mem_coe.1 hw) v (Finset.mem_coe.1 hv) he
  have himg : (T.erase 0).image ν = T.erase 0 :=
    Finset.eq_of_subset_of_card_le
      (fun v hv => by
        obtain ⟨w, hw, rfl⟩ := Finset.mem_image.1 hv
        exact hνe w hw)
      (le_of_eq (Finset.card_image_of_injOn hinjOn).symm)
  have h1 : ∑ v ∈ T.erase 0, f (P.weierstrassP v) * P.derivWeierstrassP v
      = ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w) := by
    conv_lhs => rw [← himg]
    exact Finset.sum_image (fun w hw v hv he => hinj w hw v hv he)
  have h2 : ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
      = -∑ w ∈ T.erase 0, f (P.weierstrassP w) * P.derivWeierstrassP w := by
    have hc : ∑ w ∈ T.erase 0, f (P.weierstrassP (ν w)) * P.derivWeierstrassP (ν w)
        = ∑ w ∈ T.erase 0, (-(f (P.weierstrassP w) * P.derivWeierstrassP w)) :=
      Finset.sum_congr rfl (fun w hw => hodd w (Finset.mem_of_mem_erase hw))
    rw [hc, Finset.sum_neg_distrib]
  have h3 : (2 : ℂ) * ∑ w ∈ T.erase 0,
      f (P.weierstrassP w) * P.derivWeierstrassP w = 0 := by
    linear_combination h1.trans h2
  simpa using h3

def sum_mul_derivWeierstrassP_rep_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の上で f(℘)·℘′ の和は消える——第 670 の一般化。★無条件)",
    sectionId := "genell-lemma-3-5" }

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
