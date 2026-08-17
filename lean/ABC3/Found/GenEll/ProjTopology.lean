import ABC3.Meta.Claim
import Mathlib.LinearAlgebra.Projectivization.Basic
import Mathlib.Analysis.Normed.Module.FiniteDimension
import Mathlib.Analysis.Complex.Basic
import Mathlib.LinearAlgebra.Complex.FiniteDimensional

/-!
# [GenEll] Definition 1.1 の足場 —— `ℙ(V)` の位相とコンパクト性(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

## ★★なぜこれが要るのか —— 2026-08-17 の測り直し

`Definition 1.1` は「`X^arc`(コンパクト正規複素解析空間)上のエルミート計量」を要求する。
★本日の測定で、**§1 が実際に使うのは複素解析空間ではない**ことが分かった:

| 項目 | 要るもの |
|---|---|
| `ht` の定義そのもの | 可逆層 + **1 点ごとの**ノルム。★位相すら要らない |
| `Proposition 1.4, (ii)(iii)` / `Proposition 1.6` のアルキメデス側 | ★**コンパクト空間上の連続関数は有界** |

★★**正則構造・連接層・GAGA はどこにも現れない。**
現れるのは**位相・コンパクト性・連続性**だけである。

★`X` は ℤ-固有なので `X ⊆ ℙⁿ` と埋め込め、`X(ℂ)` は `ℙⁿ(ℂ)` の閉集合になる。
**コンパクト空間の閉部分集合はコンパクト**(Hausdorff 性は要らない)。
★したがって出発点は **`ℙⁿ(ℂ)` のコンパクト性**である。

## ★mathlib 実測(2026-08-17)

`Projectivization` に**位相の instance は無い**——`Action.lean` / `Basic.lean` の
instance は代数と群作用だけである。`ℙⁿ(ℂ)` のコンパクト性も無い。

★ただし **`Projectivization K V` は `Quotient` として定義されている**ので、
商位相はそのまま載る。

## ★機構

単位球面 `S(0,1) ⊆ V` は有限次元ノルム空間ではコンパクトであり、
**`ℙ(V)` はその連続像である**(`v ↦ ‖v‖⁻¹ • v` で正規化できるから)。
★連続像なのでコンパクト。★★**Hausdorff 性は 1 度も使っていない**——要らない。
-/

namespace ABC3.Found.GenEll

open scoped LinearAlgebra.Projectivization

/-! ## ★商位相 -/

section Topology

variable (K V : Type*) [DivisionRing K] [AddCommGroup V] [Module K V]
  [TopologicalSpace V]

/-- ★**`ℙ(V)` の商位相**。

★mathlib に無い(2026-08-17 実測)。`Projectivization K V` は
`Quotient (projectivizationSetoid K V)` なので、そのまま商位相が載る。 -/
instance instTopologicalSpaceProjectivization : TopologicalSpace (ℙ K V) :=
  inferInstanceAs (TopologicalSpace (Quotient (projectivizationSetoid K V)))

variable {K V}

/-- ★商写像 `{v // v ≠ 0} → ℙ(V)` は連続。 -/
theorem continuous_projectivization_mk' :
    Continuous (Projectivization.mk' K : {v : V // v ≠ 0} → ℙ K V) :=
  continuous_quot_mk

end Topology

/-! ## ★★コンパクト性 -/

section Compact

variable {V : Type*} [NormedAddCommGroup V] [NormedSpace ℂ V] [FiniteDimensional ℂ V]

/-- ★単位球面の点は非零。 -/
theorem sphere_ne_zero (v : Metric.sphere (0 : V) 1) : (v : V) ≠ 0 := by
  intro hc
  have h : ‖(v : V)‖ = 1 := mem_sphere_zero_iff_norm.1 v.2
  rw [hc, norm_zero] at h
  exact zero_ne_one h

/-- ★単位球面から `ℙ(V)` への写像。 -/
noncomputable def sphereToProj (v : Metric.sphere (0 : V) 1) : ℙ ℂ V :=
  Projectivization.mk ℂ (v : V) (sphere_ne_zero v)

theorem continuous_sphereToProj : Continuous (sphereToProj (V := V)) := by
  have hc : Continuous (fun v : Metric.sphere (0 : V) 1 =>
      (⟨(v : V), sphere_ne_zero v⟩ : {w : V // w ≠ 0})) :=
    continuous_subtype_val.subtype_mk _
  exact continuous_projectivization_mk'.comp hc

/-- ★★**単位球面は `ℙ(V)` を覆う** —— `v ↦ ‖v‖⁻¹ • v` で正規化できるから。 -/
theorem sphereToProj_surjective : Function.Surjective (sphereToProj (V := V)) := by
  intro p
  have hv : p.rep ≠ 0 := p.rep_nonzero
  have hn : (0 : ℝ) < ‖p.rep‖ := norm_pos_iff.2 hv
  set c : ℂ := ((‖p.rep‖ : ℝ) : ℂ)⁻¹ with hcdef
  have hcne : c ≠ 0 := by
    rw [hcdef]
    simp only [ne_eq, inv_eq_zero, Complex.ofReal_eq_zero]
    exact hn.ne'
  have hnorm : ‖c • p.rep‖ = 1 := by
    rw [norm_smul, hcdef, norm_inv, Complex.norm_real, Real.norm_eq_abs,
      abs_of_pos hn]
    field_simp
  refine ⟨⟨c • p.rep, mem_sphere_zero_iff_norm.2 hnorm⟩, ?_⟩
  have hkey : Projectivization.mk ℂ (c • p.rep)
        (sphere_ne_zero ⟨c • p.rep, mem_sphere_zero_iff_norm.2 hnorm⟩)
      = Projectivization.mk ℂ p.rep p.rep_nonzero :=
    (Projectivization.mk_eq_mk_iff ℂ _ _ _ _).2 ⟨Units.mk0 c hcne, rfl⟩
  rw [sphereToProj, hkey, p.mk_rep]

/-- ★★★**`ℙ(V)` はコンパクト**(`V` は `ℂ` 上有限次元)。

★★**Hausdorff 性は使っていない。** 単位球面(コンパクト)の連続像だからである。
`Proposition 1.6` のアルキメデス側が要求するのは
「コンパクト空間上の連続関数は有界」であり、それにはこれで足りる。 -/
theorem compactSpace_projectivization : CompactSpace (ℙ ℂ V) := by
  haveI : NormedSpace ℝ V := NormedSpace.complexToReal
  haveI : FiniteDimensional ℝ V := FiniteDimensional.trans ℝ ℂ V
  haveI : CompactSpace (Metric.sphere (0 : V) 1) :=
    isCompact_iff_compactSpace.1 (isCompact_sphere (0 : V) 1)
  rw [← isCompact_univ_iff, ← Set.range_eq_univ.2 sphereToProj_surjective]
  exact isCompact_range continuous_sphereToProj

end Compact

/-! ## ★出典の紐付け(`.src`) -/

def compactSpace_projectivization.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(X^arc のコンパクト性の足場のみ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
