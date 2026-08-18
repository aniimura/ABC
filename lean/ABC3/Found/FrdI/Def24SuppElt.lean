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

/-! ## ★3. ★★★**primary ⟺ 台が 1 点**

★★これが原文の「the structure of the `Φ(A)_p`」の中身であり、
`Proposition 4.1, (ii)–(v)` の単系層はすべてここに帰着する。 -/

theorem suppElt_zero (H : IsPerfFactorialWith M ι) : SuppElt ι (0 : M) = ∅ := by
  rw [suppElt_eq, map_zero, factorMap_zero H]
  ext p
  simp [Supp]

theorem suppElt_nsmul (H : IsPerfFactorialWith M ι) (a : M) {n : ℕ} (hn : 0 < n) :
    SuppElt ι (n • a) = SuppElt ι a := by
  induction n with
  | zero => omega
  | succ m ih =>
    rcases Nat.eq_zero_or_pos m with hm | hm
    · subst hm; simp
    · rw [succ_nsmul, suppElt_add H, ih hm, Set.union_self]

theorem suppElt_mono_of_mprec (H : IsPerfFactorialWith M ι) {a b : M} (h : MPrec a b) :
    SuppElt ι a ⊆ SuppElt ι b := by
  obtain ⟨n, hn, c, hc⟩ := h
  have h1 : SuppElt ι a ⊆ SuppElt ι (n • b) := suppElt_mono_of_mle H ⟨c, hc⟩
  rwa [suppElt_nsmul H b hn] at h1

theorem suppElt_ne_empty (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {a : M} (ha : a ≠ 0) : SuppElt ι a ≠ ∅ :=
  fun h => ha (eq_zero_of_suppElt_empty H hdiv h)

/-- ★★★**primary な元の台は 1 点**。

★2 点あれば片方だけを残す分解 `a = y + z` が取れ、`y ⪯ a ⪯ y` から
`Supp a ⊆ Supp y = {P}` になって矛盾する。 -/
theorem suppElt_singleton_of_primary (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {a : M} (hp : IsPrimaryElt a) :
    ∃ P : Prime M, SuppElt ι a = {P} := by
  obtain ⟨P, hP⟩ : ∃ P, P ∈ SuppElt ι a :=
    Set.nonempty_iff_ne_empty.mpr (suppElt_ne_empty H hdiv hp.1)
  refine ⟨P, Set.eq_singleton_iff_unique_mem.mpr ⟨hP, fun P' hP' => ?_⟩⟩
  by_contra hne
  obtain ⟨y, z, hsum, hy, hz⟩ := exists_split_suppElt H hperf hdiv a {P}
  have hyP : P ∈ SuppElt ι y := by rw [hy]; exact ⟨rfl, hP⟩
  have hy0 : y ≠ 0 := by
    intro h0
    rw [h0, suppElt_zero H] at hyP
    exact hyP
  have hay : MPrec a y := hp.2 y hy0 ⟨1, one_pos, z, by rw [one_smul]; exact hsum.symm⟩
  have hsub : SuppElt ι a ⊆ SuppElt ι y := suppElt_mono_of_mprec H hay
  rw [hy] at hsub
  exact hne (hsub hP').1

/-- ★`n • a` が primary なら `a` も primary(`a ⪯ n•a ⪯ a` だから)。 -/
theorem isPrimaryElt_of_nsmul {a : M} {n : ℕ} (hn : 0 < n) (ha : a ≠ 0)
    (h : IsPrimaryElt (n • a)) : IsPrimaryElt a := by
  have h1 : MPrec a (n • a) := ⟨1, one_pos, by rw [one_smul]; exact mle_nsmul_self' hn⟩
  exact ⟨ha, fun c hc hca => mprec_trans h1 (h.2 c hc (mprec_trans hca h1))⟩

/-- ★★★**台が 1 点なら primary**(`suppElt_singleton_of_primary` の逆)。

★★手: `factorMem` で `p` 成分を実現する `x ∈ pCarrierPf M P` を取り、
`factorMap_pCarrier_self` / `factorMap_pCarrier_other` で
**`factorMap x = factorMap (of a)`** を出し、`factorInj` で `x = of a`。
すると `n • a` が(類 `P` の)primary な元になり、`isPrimaryElt_of_nsmul` で降りる。 -/
theorem isPrimaryElt_of_suppElt_singleton (H : IsPerfFactorialWith M ι)
    (hdiv : IsDivisorial M) {a : M} (ha : a ≠ 0) {P : Prime M}
    (hsupp : SuppElt ι a = {P}) : IsPrimaryElt a := by
  have hP : factorMap ι (Pf.of a) P ≠ 0 := by
    have h : P ∈ SuppElt ι a := by rw [hsupp]; rfl
    exact h
  obtain ⟨x, hxmem, hxval⟩ := H.factorMem (Pf.of a) P
  have heq : factorMap ι x = factorMap ι (Pf.of a) := by
    funext Q
    by_cases hQ : Q = P
    · subst hQ
      rw [factorMap_pCarrier_self H hxmem, hxval]
    · rw [factorMap_pCarrier_other H hdiv.2 hQ hxmem]
      by_contra hc
      refine hQ ?_
      have h : Q ∈ SuppElt ι a := fun h0 => hc h0.symm
      rw [hsupp] at h
      exact h
  have hxa : x = Pf.of a := H.factorInj heq
  have hx0 : x ≠ 0 := by
    intro h0
    rw [h0, factorMap_zero H] at heq
    exact hP (by rw [← congrFun heq P]; rfl)
  rcases hxmem with ⟨n, b, hb, hnx⟩ | hx0'
  · obtain ⟨hbp, -⟩ := hb
    rw [hxa, ← map_nsmul] at hnx
    have hnab : ((n : ℕ+) : ℕ) • a = b := Pf.of_injective_of_divisorial hdiv hnx
    exact isPrimaryElt_of_nsmul n.2 ha (hnab ▸ hbp)
  · exact absurd hx0' hx0

/-- ★★★★**primary ⟺ 台が 1 点**。 -/
theorem isPrimaryElt_iff_suppElt_singleton (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {a : M} (ha : a ≠ 0) :
    IsPrimaryElt a ↔ ∃ P : Prime M, SuppElt ι a = {P} :=
  ⟨suppElt_singleton_of_primary H hperf hdiv,
    fun ⟨_, hP⟩ => isPrimaryElt_of_suppElt_singleton H hdiv ha hP⟩

/-! ## ★出典の紐付け(条つき) -/

def Pf.of_injective_of_divisorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11,
    item := "§0 Monoids — perfection M^pf の単射性",
    sectionId := "frdi-s0-perfect" }

end ABC3.Found.FrdI
