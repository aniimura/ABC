/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
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
