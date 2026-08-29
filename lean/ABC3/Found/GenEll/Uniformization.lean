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
