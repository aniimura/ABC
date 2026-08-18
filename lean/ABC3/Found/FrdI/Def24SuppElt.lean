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

/-! ## ★4. 同じ 1 点の台を持つ 2 元 —— `Proposition 4.1, (iii)` の「disjoint supports」 -/

/-- ★★**台が 1 点なら `M^pf` の像は `p` の担い手に入る**。 -/
theorem of_mem_pCarrierPf_of_suppElt_singleton (H : IsPerfFactorialWith M ι)
    (hdiv : IsDivisorial M) {a : M} {P : Prime M} (hsupp : SuppElt ι a = {P}) :
    Pf.of a ∈ pCarrierPf M P := by
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
  rw [← H.factorInj heq]
  exact hxmem

/-- ★★**同じ 1 点の台を持つ 2 元は `⪯`-同値**。

★どちらも `p` の担い手に入るので、`n • z` と `m • w` が類 `p` の primary 元になり、
`toPrime_eq_iff` で繋がる。 -/
theorem mprec_of_suppElt_eq_singleton (H : IsPerfFactorialWith M ι) (hdiv : IsDivisorial M)
    {z w : M} {P : Prime M} (hzs : SuppElt ι z = {P}) (hws : SuppElt ι w = {P}) :
    MPrec z w := by
  rcases of_mem_pCarrierPf_of_suppElt_singleton H hdiv hzs with ⟨n, bz, hbz, hnz⟩ | hz0
  · rcases of_mem_pCarrierPf_of_suppElt_singleton H hdiv hws with ⟨m, bw, hbw, hmw⟩ | hw0
    · obtain ⟨hbzp, hbze⟩ := hbz
      obtain ⟨hbwp, hbwe⟩ := hbw
      have hcls : MPrec bz bw := (toPrime_eq_iff hbzp hbwp).mp (by rw [hbze, hbwe])
      have h1 : MPrec (Pf.of z) (Pf.of bz) := by
        refine ⟨1, one_pos, ?_⟩
        rw [one_smul, ← hnz]
        exact mle_nsmul_self' n.2
      have h2 : MPrec (Pf.of bz) (Pf.of bw) := mprec_pf_of_iff.mpr hcls
      have h3 : MPrec (Pf.of bw) (Pf.of w) :=
        ⟨((m : ℕ+) : ℕ), m.2, ⟨0, by rw [add_zero]; exact hmw.symm⟩⟩
      exact mprec_pf_of_iff.mp (mprec_trans h1 (mprec_trans h2 h3))
    · exfalso
      have hemp : SuppElt ι w = ∅ := by
        show Supp (factorMap ι (Pf.of w)) = ∅
        rw [hw0, factorMap_zero H]
        ext p; simp [Supp]
      rw [hws] at hemp
      exact absurd hemp.symm (by simp)
  · exfalso
    have hemp : SuppElt ι z = ∅ := by
      show Supp (factorMap ι (Pf.of z)) = ∅
      rw [hz0, factorMap_zero H]
      ext p; simp [Supp]
    rw [hzs] at hemp
    exact absurd hemp.symm (by simp)

/-- ★★**`z ⪯ w` なら共通の非零下界がある**(perfect で `n` 割り)。 -/
theorem exists_common_lower_of_mprec (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M)
    {z w : M} (hz : z ≠ 0) (h : MPrec z w) :
    ∃ d : M, d ≠ 0 ∧ MLe d z ∧ MLe d w := by
  obtain ⟨n, hn, c, hc⟩ := h
  obtain ⟨d, hd⟩ := (hperf ⟨n, hn⟩).2 z
  have hnd : n • d = z := hd
  refine ⟨d, ?_, ?_, ?_⟩
  · intro h0
    exact hz (by rw [← hnd, h0, smul_zero])
  · exact ⟨(n - 1) • d, by rw [← hnd, ← succ_nsmul']; congr 1; omega⟩
  · obtain ⟨c', hc'⟩ := (hperf ⟨n, hn⟩).2 c
    have hnc' : n • c' = c := hc'
    refine ⟨c', nsmul_inj_of_divisorial hdiv hn ?_⟩
    rw [smul_add, hnd, hnc']
    exact hc

/-- ★★★**台が交わらない ⟺ 共通の非零下界が無い**。

★★これが `Proposition 4.1, (iii)` の「disjoint supports」の単系層である。 -/
theorem suppElt_disjoint_iff (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) (a b : M) :
    SuppElt ι a ∩ SuppElt ι b = ∅ ↔ ∀ d : M, MLe d a → MLe d b → d = 0 := by
  constructor
  · intro hdisj d hda hdb
    exact eq_zero_of_mle_of_suppElt_disjoint H hdiv hda hdb hdisj
  · intro hcond
    refine Set.eq_empty_iff_forall_notMem.mpr (fun P hP => ?_)
    obtain ⟨y, z, hsum, hy, hz⟩ := exists_split_suppElt H hperf hdiv a {P}
    obtain ⟨y', z', hsum', hy', hz'⟩ := exists_split_suppElt H hperf hdiv b {P}
    have hzs : SuppElt ι y = {P} := by
      rw [hy]
      ext Q
      simp only [Set.mem_inter_iff, Set.mem_singleton_iff]
      constructor
      · exact fun h => h.1
      · rintro rfl
        exact ⟨rfl, hP.1⟩
    have hws : SuppElt ι y' = {P} := by
      rw [hy']
      ext Q
      simp only [Set.mem_inter_iff, Set.mem_singleton_iff]
      constructor
      · exact fun h => h.1
      · rintro rfl
        exact ⟨rfl, hP.2⟩
    have hy0 : y ≠ 0 := by
      intro h0
      rw [h0, suppElt_zero H] at hzs
      exact absurd hzs.symm (by simp)
    obtain ⟨d, hd0, hdy, hdy'⟩ := exists_common_lower_of_mprec hperf hdiv hy0
      (mprec_of_suppElt_eq_singleton H hdiv hzs hws)
    exact hd0 (hcond d (mle_trans hdy ⟨z, hsum.symm⟩) (mle_trans hdy' ⟨z', hsum'.symm⟩))


/-! ## ★5. ★★★★`Proposition 4.1, (iv)(v)` の単系層

原文 (FrdI p.77):
> For every primary element x ∈/ p, x ≤ x if and only if x ≤ xδ + x .

★原文の条件は
「**`p` に属さない任意の primary 元 `y` について `y ≼ e ⟺ y ≼ d + e`**」
であり、これは **`d` の台が 1 点 `{p}` である**ことと同値である。
★★これが (iv)(v) の中身であり、圏の側では
(iv) は**後置**の圏同値、(v) は**前置**の圏同値でこの単系の条件に翻訳される。 -/

/-- ★★素点の集合 `S` での分解を、**因子の値まで込めて**取る。 -/
theorem exists_split_factorMap (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) (a : M) (S : Set (Prime M)) :
    ∃ y z : M, a = y + z ∧
      factorMap ι (Pf.of y) = Set.indicator S (factorMap ι (Pf.of a)) ∧
      factorMap ι (Pf.of z) = Set.indicator Sᶜ (factorMap ι (Pf.of a)) := by
  obtain ⟨y₀, z₀, hsum, hy, hz⟩ := exists_split_supp H (Pf.of a) S
  obtain ⟨y, rfl⟩ := Pf.of_surjective_of_perfect hperf y₀
  obtain ⟨z, rfl⟩ := Pf.of_surjective_of_perfect hperf z₀
  exact ⟨y, z, Pf.of_injective_of_divisorial hdiv (by rw [map_add]; exact hsum), hy, hz⟩

/-- ★★**台が `S` に収まる元は、`x` の下にあるなら `x` の `S` 成分の下にある**。

★★これが「素点ごとに比べる」ことの実体である —— 差 `c` を `S` で割り、
`S` の上では `y` がそのまま残ることを使う。 -/
theorem mle_of_restrict (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) {x y xS : M} {S : Set (Prime M)}
    (hy : SuppElt ι y ⊆ S) (h : MLe y x)
    (hxS : factorMap ι (Pf.of xS) = Set.indicator S (factorMap ι (Pf.of x))) :
    MLe y xS := by
  classical
  obtain ⟨c, hc⟩ := h
  obtain ⟨cS, cSc, hcsum, hcS, -⟩ := exists_split_factorMap H hperf hdiv c S
  refine ⟨cS, Pf.of_injective_of_divisorial hdiv (H.factorInj ?_)⟩
  rw [map_add, H.factorAdd, hcS, hxS, ← hc, map_add, H.factorAdd]
  have hself : Set.indicator S (factorMap ι (Pf.of y)) = factorMap ι (Pf.of y) :=
    Set.indicator_eq_self.mpr hy
  have hadd := Set.indicator_add (M := ℝ≥0) S (factorMap ι (Pf.of y)) (factorMap ι (Pf.of c))
  rw [show (fun r => factorMap ι (Pf.of y) r + factorMap ι (Pf.of c) r)
      = factorMap ι (Pf.of y) + factorMap ι (Pf.of c) from rfl] at hadd
  rw [hadd, hself]
  rfl

/-- ★★**`d` の台と交わらない台を持つ `y` は、`d + e` の下にあれば `e` の下にある**。 -/
theorem mle_of_mle_add_of_suppElt_disjoint (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {d e y : M} {q : Prime M}
    (hy : SuppElt ι y ⊆ {q}) (hq : q ∉ SuppElt ι d) (h : MLe y (d + e)) : MLe y e := by
  classical
  obtain ⟨u, v, hsum, hu, -⟩ :=
    exists_split_factorMap H hperf hdiv (d + e) ({q} : Set (Prime M))
  obtain ⟨u', v', hsum', hu', -⟩ := exists_split_factorMap H hperf hdiv e ({q} : Set (Prime M))
  have hdq : factorMap ι (Pf.of d) q = 0 := by
    by_contra hc
    exact hq (show q ∈ Supp (factorMap ι (Pf.of d)) from hc)
  have heq : factorMap ι (Pf.of u) = factorMap ι (Pf.of u') := by
    rw [hu, hu']
    funext r
    by_cases hr : r ∈ ({q} : Set (Prime M))
    · have hrq : r = q := hr
      subst hrq
      rw [Set.indicator_of_mem hr, Set.indicator_of_mem hr, map_add, H.factorAdd,
        Pi.add_apply, hdq, zero_add]
    · rw [Set.indicator_of_notMem hr, Set.indicator_of_notMem hr]
  have huu' : u = u' := Pf.of_injective_of_divisorial hdiv (H.factorInj heq)
  exact mle_trans (mle_of_restrict H hperf hdiv hy h hu) (huu' ▸ ⟨v', hsum'.symm⟩)

/-- ★★★**原文の条件が成り立てば `d` の台は `{p}` に収まる**。

★`d + e` の `q` 成分だけを残した `u` が primary になり、
`u ≼ e` から `q` 成分で `d` の分が消えることが出る。 -/
theorem suppElt_subset_of_cond (H : IsPerfFactorialWith M ι) (hperf : IsPerfectMonoid M)
    (hdiv : IsDivisorial M) {d e : M} {p : Prime M}
    (hcond : ∀ y : M, IsPrimaryElt y → SuppElt ι y ≠ {p} → (MLe y e ↔ MLe y (d + e))) :
    SuppElt ι d ⊆ {p} := by
  classical
  intro q hq
  by_contra hqp
  obtain ⟨u, v, hsum, hu, -⟩ :=
    exists_split_factorMap H hperf hdiv (d + e) ({q} : Set (Prime M))
  have hqde : q ∈ SuppElt ι (d + e) := by
    rw [suppElt_add H]; exact Or.inl hq
  have hsuppu : SuppElt ι u = ({q} : Set (Prime M)) := by
    show Supp (factorMap ι (Pf.of u)) = _
    rw [hu, supp_indicator]
    exact Set.inter_eq_left.mpr (fun r hr => by rw [show r = q from hr]; exact hqde)
  have hu0 : u ≠ 0 := by
    intro h0
    rw [h0, suppElt_zero H] at hsuppu
    exact absurd hsuppu.symm (by simp)
  have hup : IsPrimaryElt u := isPrimaryElt_of_suppElt_singleton H hdiv hu0 hsuppu
  have hne : SuppElt ι u ≠ {p} := by
    rw [hsuppu]
    intro hcontra
    exact hqp (by rw [← Set.singleton_eq_singleton_iff.mp hcontra]; rfl)
  obtain ⟨c, hc⟩ := (hcond u hup hne).mpr ⟨v, hsum.symm⟩
  have h1 : factorMap ι (Pf.of u) q + factorMap ι (Pf.of c) q = factorMap ι (Pf.of e) q := by
    rw [← hc, map_add, H.factorAdd, Pi.add_apply]
  have h2 : factorMap ι (Pf.of u) q
      = factorMap ι (Pf.of d) q + factorMap ι (Pf.of e) q := by
    rw [hu, Set.indicator_of_mem (Set.mem_singleton q), map_add, H.factorAdd, Pi.add_apply]
  have h3 : factorMap ι (Pf.of d) q + factorMap ι (Pf.of c) q = 0 := by
    refine add_right_cancel (b := factorMap ι (Pf.of e) q) ?_
    rw [zero_add, add_right_comm, ← h2, h1]
  exact hq (add_eq_zero.mp h3).1

/-- ★★★★**[FrdI] Proposition 4.1, (iv)(v) の単系層**。

原文 (FrdI p.77):
> The necessity and sufficiency of this condition then follow immediately by com-

★**`d` が primary ⟺ ある素点 `p` があって、`p` に属さない primary 元 `y` については
`y ≼ e` と `y ≼ d + e` が同値**。★(iv) と (v) は、圏の側で
後置／前置のどちらの圏同値を使うかだけが違い、単系の内容はこの 1 本である。 -/
theorem isPrimaryElt_iff_exists_prime_cond (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {d e : M} (hd : d ≠ 0) :
    IsPrimaryElt d ↔ ∃ p : Prime M, ∀ y : M, IsPrimaryElt y → SuppElt ι y ≠ {p} →
      (MLe y e ↔ MLe y (d + e)) := by
  classical
  constructor
  · intro hp
    obtain ⟨P, hP⟩ := suppElt_singleton_of_primary H hperf hdiv hp
    refine ⟨P, fun y hy hne => ⟨fun h => mle_trans h ⟨d, add_comm _ _⟩, fun h => ?_⟩⟩
    obtain ⟨q, hq⟩ := suppElt_singleton_of_primary H hperf hdiv hy
    have hqP : q ≠ P := by rintro rfl; exact hne hq
    refine mle_of_mle_add_of_suppElt_disjoint H hperf hdiv (q := q) (by rw [hq]) ?_ h
    rw [hP]
    exact hqP
  · rintro ⟨p, hcond⟩
    have hsub : SuppElt ι d ⊆ {p} := suppElt_subset_of_cond H hperf hdiv hcond
    obtain ⟨q, hq⟩ : ∃ q, q ∈ SuppElt ι d :=
      Set.nonempty_iff_ne_empty.mpr (suppElt_ne_empty H hdiv hd)
    refine isPrimaryElt_of_suppElt_singleton H hdiv hd (P := p) ?_
    exact Set.eq_singleton_iff_unique_mem.mpr ⟨hsub hq ▸ hq, fun r hr => hsub hr⟩

/-- ★★★**primary な元の台は、その元の属する素点 1 点**。

★★これにより、原文の条件「`x_ϵ′ ∉ p`」が
`SuppElt ι x_ϵ′ ≠ {p}` と**同じこと**であることが分かる。 -/
theorem suppElt_eq_singleton_toPrime (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {a : M} (ha : IsPrimaryElt a) :
    SuppElt ι a = {toPrime M a ha} := by
  obtain ⟨P, hP⟩ := suppElt_singleton_of_primary H hperf hdiv ha
  obtain ⟨n, b, hb, hnb⟩ | h0 := of_mem_pCarrierPf_of_suppElt_singleton H hdiv hP
  · obtain ⟨hbp, hbe⟩ := hb
    rw [← map_nsmul] at hnb
    have hnab : ((n : ℕ+) : ℕ) • a = b := Pf.of_injective_of_divisorial hdiv hnb
    have hab : MPrec a b := by
      refine ⟨1, one_pos, ?_⟩
      rw [one_smul, ← hnab]
      exact mle_nsmul_self' n.2
    rw [hP, ← hbe, (toPrime_eq_iff ha hbp).mpr hab]
  · exact absurd (by rw [suppElt_eq, h0, factorMap_zero H]; ext r; simp [Supp] :
      SuppElt ι a = (∅ : Set (Prime M))) (by rw [hP]; simp)


/-! ## ★6. 和の下界 —— `Proposition 4.1, (iii)` の cartesian 性

原文 (FrdI p.77):
> follows from the fact that "for xU ∈ Φ(F ), x + xι ≤ xU if and only if x ≤ xU ,

★原文が cartesian 性の理由として挙げるのがこの 1 行である。 -/

/-- ★★★**台が交わらない 2 元では、和が下にあることと両方が下にあることが同値**。

★`c` を `Supp a` とその補集合で割り、それぞれに `a`, `b` が収まることを使う。 -/
theorem mle_add_iff_of_suppElt_disjoint (H : IsPerfFactorialWith M ι)
    (hperf : IsPerfectMonoid M) (hdiv : IsDivisorial M) {a b c : M}
    (hdisj : SuppElt ι a ∩ SuppElt ι b = ∅) :
    MLe (a + b) c ↔ MLe a c ∧ MLe b c := by
  classical
  refine ⟨fun h => ⟨mle_trans ⟨b, rfl⟩ h, mle_trans ⟨a, add_comm b a⟩ h⟩, fun ⟨ha, hb⟩ => ?_⟩
  obtain ⟨cS, cSc, hsum, hcS, hcSc⟩ :=
    exists_split_factorMap H hperf hdiv c (SuppElt ι a)
  have h1 : MLe a cS := mle_of_restrict H hperf hdiv (subset_refl _) ha hcS
  have h2 : MLe b cSc := by
    refine mle_of_restrict H hperf hdiv (fun r hr => ?_) hb hcSc
    intro hra
    exact Set.eq_empty_iff_forall_notMem.mp hdisj r ⟨hra, hr⟩
  obtain ⟨d1, hd1⟩ := h1
  obtain ⟨d2, hd2⟩ := h2
  exact ⟨d1 + d2, by rw [hsum, ← hd1, ← hd2]; abel⟩


/-! ## ★出典の紐付け(条つき) -/

def Pf.of_injective_of_divisorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 11,
    item := "§0 Monoids — perfection M^pf の単射性",
    sectionId := "frdi-s0-perfect" }

end ABC3.Found.FrdI
