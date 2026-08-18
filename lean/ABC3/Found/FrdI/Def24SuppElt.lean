/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24Supp
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Prop32Equiv

/-!
# [FrdI] 台の理論を **`M` の元**へ降ろす —— `M ≅ M^pf`(perfect ＋ divisorial)

原文 (FrdI p.11):
> which is injective if M is torsion-free, integral, and saturated, hence, in particular, if

`Def24Supp.lean` は台の理論を **`M^pf` の元**について作った。
★★ところが `Proposition 4.1` が扱うのは **`Φ(A)` の元**である。

## ★★橋 —— `Pf.of : M → M^pf` は全単射になる

| 向き | 要る条件 | 宣言 |
|---|---|---|
| 単射 | ★`divisorial`(`n •` が単射になる) | `Pf.of_injective_of_divisorial` |
| 全射 | ★`perfect`(`mk a n = of b`、`n • b = a`) | `Pf.of_surjective_of_perfect` |

★★**単射の中身**: `n • a = n • b` から `MLe (n•a) (n•b)` を両向きに取り、
`mle_of_nsmul_mle`(saturated ＋ integral)で `MLe a b` を両向きに落とし、
`mle_antisymm`(integral ＋ sharp)で潰す。
★`divisorial` の 4 条のうち **3 条(integral・saturated・sharp)がここで効く**。

★`Proposition 4.1` は `𝒞` が **perfect かつ isotropic 型**、`Φ` が **perf-factorial**
であることを仮定するので、`Φ(A)` は perfect(`Proposition 1.10, (iii)`)かつ
divisorial(`PreFrobenioid.divisorial`)であり、両条件が揃う。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe w

variable {M : Type w} [AddCommMonoid M]

/-! ## ★1. `M ≅ M^pf` -/

/-- ★★**divisorial なら `n •` は単射**。

★`divisorial` の 3 条(saturated・integral・sharp)がここで効く。 -/
theorem nsmul_inj_of_divisorial (hdiv : IsDivisorial M) {n : ℕ} (hn : 0 < n) {a b : M}
    (h : n • a = n • b) : a = b :=
  mle_antisymm hdiv.1.1 hdiv.2
    (mle_of_nsmul_mle hdiv.1.2.1 hdiv.1.1 hn ⟨0, by rw [h, add_zero]⟩)
    (mle_of_nsmul_mle hdiv.1.2.1 hdiv.1.1 hn ⟨0, by rw [← h, add_zero]⟩)

/-- ★**divisorial なら `M → M^pf` は単射**。 -/
theorem Pf.of_injective_of_divisorial (hdiv : IsDivisorial M) :
    Function.Injective (Pf.of : M →+ Pf M) := by
  intro a b h
  obtain ⟨j, hj⟩ := Pf.of_eq_of_iff.mp h
  exact nsmul_inj_of_divisorial hdiv j.2 hj

/-- ★★**perfect なら `M → M^pf` は全射** —— `mk a n` は `n • b = a` なる `b` の像。 -/
theorem Pf.of_surjective_of_perfect (hperf : IsPerfectMonoid M) :
    Function.Surjective (Pf.of : M →+ Pf M) := by
  intro x
  induction x using Pf.inductionOn with
  | _ a n =>
    obtain ⟨b, hb⟩ := (hperf n).2 a
    refine ⟨b, ?_⟩
    rw [Pf.of_apply, ← hb]
    exact (Pf.sound 1 (by simp)).symm

/-! ## ★2. `M` の元の台 -/

variable {ι : Prime M → Pf M → ℝ≥0}

theorem suppElt_eq (a : M) : SuppElt ι a = Supp (factorMap ι (Pf.of a)) := rfl

/-- ★元の台は加法的。 -/
theorem suppElt_add (H : IsPerfFactorialWith M ι) (a b : M) :
    SuppElt ι (a + b) = SuppElt ι a ∪ SuppElt ι b := by
  rw [suppElt_eq, suppElt_eq, suppElt_eq, map_add, supp_factorMap_add H]

/-- ★`MLe` は元の台の包含を与える。 -/
theorem suppElt_mono_of_mle (H : IsPerfFactorialWith M ι) {a b : M} (h : MLe a b) :
    SuppElt ι a ⊆ SuppElt ι b := by
  obtain ⟨c, rfl⟩ := h
  rw [suppElt_add H]
  exact Set.subset_union_left

/-- ★台が空なら `0`(`M` が divisorial なら)。 -/
theorem eq_zero_of_suppElt_empty (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {a : M} (h : SuppElt ι a = ∅) : a = 0 := by
  refine Pf.of_injective_of_divisorial hdiv ?_
  rw [map_zero]
  exact eq_zero_of_supp_empty H h

/-- ★★★**素点の集合 `S` で `M` の元を割る**(`M` が perfect かつ divisorial なら)。

★★これが原文の「primary factorizations」の、`Φ(A)` の元についての形である。 -/
theorem exists_split_suppElt (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) (a : M) (S : Set (Prime M)) :
    ∃ y z : M, a = y + z ∧ SuppElt ι y = S ∩ SuppElt ι a ∧
      SuppElt ι z = Sᶜ ∩ SuppElt ι a := by
  obtain ⟨y₀, z₀, hsum, hy, hz⟩ := exists_split_supp H (Pf.of a) S
  obtain ⟨y, rfl⟩ := Pf.of_surjective_of_perfect hperf y₀
  obtain ⟨z, rfl⟩ := Pf.of_surjective_of_perfect hperf z₀
  refine ⟨y, z, Pf.of_injective_of_divisorial hdiv (by rw [map_add]; exact hsum), ?_, ?_⟩
  · show Supp (factorMap ι (Pf.of y)) = S ∩ Supp (factorMap ι (Pf.of a))
    rw [hy, supp_indicator]
  · show Supp (factorMap ι (Pf.of z)) = Sᶜ ∩ Supp (factorMap ι (Pf.of a))
    rw [hz, supp_indicator]

/-- ★★**台が交わらない 2 元に共通の非零下界は無い**(`M` の元の版)。

★★`Proposition 4.1, (ii)(iii)` の「disjoint supports」の条の骨である。 -/
theorem eq_zero_of_mle_of_suppElt_disjoint (H : IsPerfFactorialWith M ι)
    (hdiv : IsDivisorial M) {a b d : M} (hda : MLe d a) (hdb : MLe d b)
    (hdisj : SuppElt ι a ∩ SuppElt ι b = ∅) : d = 0 := by
  refine eq_zero_of_suppElt_empty H hdiv (Set.eq_empty_iff_forall_notMem.mpr fun p hp => ?_)
  exact absurd (Set.eq_empty_iff_forall_notMem.mp hdisj p)
    (by simp only [not_not]
        exact ⟨suppElt_mono_of_mle H hda hp, suppElt_mono_of_mle H hdb hp⟩)

/-! ## ★出典の紐付け(条つき) -/

def Pf.of_injective_of_divisorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11,
    item := "§0 Monoids — perfection M^pf の単射性",
    sectionId := "frdi-s0-perfect" }

end ABC3.Found.FrdI
