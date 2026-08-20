/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.CartierPf
import Mathlib.RingTheory.Int.Basic

/-!
# `Φ(L)` の各素点での局所化 `M_p` は monoprime

★★★[FrdI] `Definition 2.4, (i)` の条件 **(b)**「各 `p ∈ Prime(M)` で `M_p` は
monoprime」を、`Φ := Γ ∩ ℤ≥0[S]` について閉じる
(`Example 6.1` の「one verifies immediately that `Φ(L)` is perf-factorial」の一部)。

## ★中身

1. `Mp_eq_suppAtSub` —— **`M_p` はちょうど「台が `{s}` に入る `Φ` の元」**
   (`s` は `p` に対応する素因子)。`primeCarrier` が「台がちょうど `{s}`」
   (`CartierPrime.lean` の primary の判定)なので、その生成する部分単系がこれになる。
2. `sGrpAt Γ s := {c : ℤ | c·s ∈ Γ}` は `ℤ` の**部分群**で、`Q`-Cartier だから
   `n_s > 0` を含む。したがって**正の生成元 `d`** を持つ
   (`ℤ` の部分群は単項 —— `Submodule.IsPrincipal`)。
3. よって `M_p ≃+ ℕ`(`k ↦ k·d·s`)、すなわち **`ℤ`-monoprime**。

★原文の `Φ(L)^pf_p = ℚ≥0` に当たるのは、これを完全化した姿である
(`CartierPf.lean` の `pfEquivNonneg` が `Φ^pf ≃ ℚ≥0[D_L]` を与えている)。
-/

namespace ABC3.Found.FrdI

open Finsupp

/-! ## ★1. `ℤ` の部分群は単項 -/

/-- ★`ℤ` の部分群は単項 —— `c ∈ H ↔ a ∣ c`。 -/
theorem exists_gen_addSubgroupInt (H : AddSubgroup ℤ) : ∃ a : ℤ, ∀ c : ℤ, c ∈ H ↔ a ∣ c := by
  set S := AddSubgroup.toIntSubmodule H with hS
  set g := Submodule.IsPrincipal.generator S with hg
  have hspan : S = Submodule.span ℤ {g} :=
    (Submodule.IsPrincipal.span_singleton_generator S).symm
  refine ⟨g, fun c => ?_⟩
  constructor
  · intro hc
    have hcS : c ∈ S := hc
    rw [hspan, Submodule.mem_span_singleton] at hcS
    obtain ⟨b, hb⟩ := hcS
    exact ⟨b, by rw [← hb, smul_eq_mul, mul_comm]⟩
  · rintro ⟨b, hb⟩
    have hcS : c ∈ S := by
      rw [hspan, Submodule.mem_span_singleton]
      exact ⟨b, by rw [smul_eq_mul, mul_comm]; exact hb.symm⟩
    exact hcS

/-- ★正の元を持つ `ℤ` の部分群は、**正の**生成元を持つ。 -/
theorem exists_pos_gen_addSubgroupInt (H : AddSubgroup ℤ) {n : ℤ} (hn : 0 < n) (hnH : n ∈ H) :
    ∃ d : ℤ, 0 < d ∧ ∀ c : ℤ, c ∈ H ↔ d ∣ c := by
  obtain ⟨a, ha⟩ := exists_gen_addSubgroupInt H
  have hane : a ≠ 0 := by
    intro h0
    have h1 := (ha n).mp hnH
    rw [h0] at h1
    have hz : n = 0 := zero_dvd_iff.mp h1
    omega
  refine ⟨|a|, abs_pos.mpr hane, fun c => ?_⟩
  rw [abs_dvd]
  exact ha c

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℤ)}

/-! ## ★2. `M_p` の同定 -/

/-- ★`s` に台を持つ `Φ` の元がなす部分単系。 -/
def suppAtSub (Γ : AddSubgroup (S →₀ ℤ)) (s : S) : AddSubmonoid (effSub Γ) where
  carrier := {a : effSub Γ | ((a : effSub Γ) : S →₀ ℤ).support ⊆ {s}}
  add_mem' {a b} ha hb := by
    classical
    have hsub : (((a + b : effSub Γ)) : S →₀ ℤ).support
        ⊆ ((a : effSub Γ) : S →₀ ℤ).support ∪ ((b : effSub Γ) : S →₀ ℤ).support := by
      show (((a : effSub Γ) : S →₀ ℤ) + ((b : effSub Γ) : S →₀ ℤ)).support ⊆ _
      exact Finsupp.support_add
    exact hsub.trans (Finset.union_subset ha hb)
  zero_mem' := by
    show ((0 : effSub Γ) : S →₀ ℤ).support ⊆ {s}
    simp

theorem mem_suppAtSub {s : S} {a : effSub Γ} :
    a ∈ suppAtSub Γ s ↔ ((a : effSub Γ) : S →₀ ℤ).support ⊆ {s} := Iff.rfl

/-- ★`primeCarrier` は「台がちょうど `{s}`」。 -/
theorem mem_primeCarrier_iff (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) (a : effSub Γ) :
    a ∈ primeCarrier (effSub Γ) p
      ↔ ((a : effSub Γ) : S →₀ ℤ).support = {effSubPrimeEquiv hQ p} := by
  constructor
  · rintro ⟨ha, hp⟩
    obtain ⟨t, ht⟩ := support_singleton_of_isPrimaryElt hQ ha
    have hpt : effSubPrimeEquiv hQ p = t := by
      rw [← hp]
      exact primePt_toPrime hQ ha ht
    rw [hpt, ht]
  · intro h
    have hprim : IsPrimaryElt a := isPrimaryElt_of_support_singleton h
    refine ⟨hprim, ?_⟩
    have h1 : effSubPrimeEquiv hQ (toPrime _ a hprim) = effSubPrimeEquiv hQ p :=
      primePt_toPrime hQ hprim h
    exact (effSubPrimeEquiv hQ).injective h1

/-- ★★**`M_p` はちょうど「台が `{s}` に入る元」**。 -/
theorem Mp_eq_suppAtSub (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) :
    Mp (effSub Γ) p = suppAtSub Γ (effSubPrimeEquiv hQ p) := by
  refine le_antisymm (AddSubmonoid.closure_le.mpr ?_) ?_
  · intro a ha
    exact mem_suppAtSub.mpr ((mem_primeCarrier_iff hQ p a).mp ha).le
  · intro a ha
    classical
    by_cases h0 : a = 0
    · rw [h0]; exact AddSubmonoid.zero_mem _
    · have hane : ((a : effSub Γ) : S →₀ ℤ) ≠ 0 := fun hz => h0 ((effSub_eq_zero_iff a).mpr hz)
      obtain ⟨t, ht⟩ := Finsupp.support_nonempty_iff.mpr hane
      have hsub := mem_suppAtSub.mp ha
      have hts : t = effSubPrimeEquiv hQ p := Finset.mem_singleton.mp (hsub ht)
      have heq : ((a : effSub Γ) : S →₀ ℤ).support = {effSubPrimeEquiv hQ p} :=
        Finset.Subset.antisymm hsub (by rw [← hts]; exact Finset.singleton_subset_iff.mpr ht)
      exact mem_Mp_of_mem_primeCarrier ((mem_primeCarrier_iff hQ p a).mpr heq)

/-! ## ★3. `M_p ≃+ ℕ` -/

/-- ★`{c : ℤ | c·s ∈ Γ}` —— `ℤ` の部分群。 -/
noncomputable def sGrpAt (Γ : AddSubgroup (S →₀ ℤ)) (s : S) : AddSubgroup ℤ :=
  Γ.comap (Finsupp.singleAddHom s)

theorem mem_sGrpAt {s : S} {c : ℤ} : c ∈ sGrpAt Γ s ↔ single s c ∈ Γ := Iff.rfl

/-- ★`Q`-Cartier なら `sGrpAt` は**正の生成元**を持つ。 -/
theorem exists_pos_gen_sGrpAt (hQ : IsQCartierSubgroup Γ) (s : S) :
    ∃ d : ℤ, 0 < d ∧ ∀ c : ℤ, single s c ∈ Γ ↔ d ∣ c := by
  obtain ⟨n, hn, hnmem⟩ := hQ s
  have hnZ : (0 : ℤ) < (n : ℤ) := by exact_mod_cast hn
  obtain ⟨d, hd, hdiff⟩ := exists_pos_gen_addSubgroupInt (sGrpAt Γ s) hnZ (mem_sGrpAt.mpr hnmem)
  exact ⟨d, hd, fun c => (mem_sGrpAt (Γ := Γ) (s := s) (c := c)).symm.trans (hdiff c)⟩

/-- ★台が `{s}` に入る元は `single s (a s)` である。 -/
theorem eq_single_of_support_subset {x : S →₀ ℤ} {s : S} (h : x.support ⊆ {s}) :
    x = single s (x s) := by
  classical
  ext t
  rcases eq_or_ne s t with rfl | hst
  · simp
  · have hns : t ∉ x.support := fun ht => hst (Finset.mem_singleton.mp (h ht)).symm
    rw [Finsupp.notMem_support_iff.mp hns, Finsupp.single_apply, if_neg hst]

/-- ★★**`{a ∈ Φ | supp a ⊆ {s}} ≃+ ℕ`** —— `Q`-Cartier のもとで。

★`{c | c·s ∈ Γ}` は `ℤ` の部分群で `n_s > 0` を含むから、正の生成元 `d` を持つ。 -/
theorem nonempty_suppAtSubEquivNat (hQ : IsQCartierSubgroup Γ) (s : S) :
    Nonempty (suppAtSub Γ s ≃+ ℕ) := by
  classical
  obtain ⟨d, hd, hdiv⟩ := exists_pos_gen_sGrpAt hQ s
  have hmem : ∀ k : ℕ, (single s ((k : ℤ) * d) : S →₀ ℤ) ∈ effSub Γ := by
    intro k
    refine mem_effSub.mpr ⟨(hdiv _).mpr ⟨(k : ℤ), by ring⟩, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simp only [Finsupp.single_eq_same]
      positivity
    · simp [hst]
  have hsupp : ∀ k : ℕ, (⟨_, hmem k⟩ : effSub Γ) ∈ suppAtSub Γ s :=
    fun k => mem_suppAtSub.mpr Finsupp.support_single_subset
  refine ⟨(AddEquiv.ofBijective
    (⟨⟨fun k : ℕ => (⟨⟨_, hmem k⟩, hsupp k⟩ : suppAtSub Γ s), ?_⟩, ?_⟩ : ℕ →+ suppAtSub Γ s)
    ?_).symm⟩
  · refine Subtype.ext (Subtype.ext ?_)
    show (single s (((0 : ℕ) : ℤ) * d) : S →₀ ℤ) = 0
    simp
  · intro k l
    refine Subtype.ext (Subtype.ext ?_)
    show (single s (((k + l : ℕ) : ℤ) * d) : S →₀ ℤ)
      = single s ((k : ℤ) * d) + single s ((l : ℤ) * d)
    rw [← Finsupp.single_add]
    congr 1
    push_cast
    ring
  · constructor
    · intro k l hkl
      have h1 : (single s ((k : ℤ) * d) : S →₀ ℤ) = single s ((l : ℤ) * d) :=
        congrArg (fun t : suppAtSub Γ s => ((t : effSub Γ) : S →₀ ℤ)) hkl
      have h2 := congrArg (fun t : S →₀ ℤ => t s) h1
      simp only [Finsupp.single_eq_same] at h2
      have h3 : (k : ℤ) = (l : ℤ) := mul_right_cancel₀ (ne_of_gt hd) h2
      exact_mod_cast h3
    · rintro ⟨⟨a, hae⟩, hsu⟩
      have hsub : a.support ⊆ {s} := mem_suppAtSub.mp hsu
      have haeq : a = single s (a s) := eq_single_of_support_subset hsub
      have hann : 0 ≤ a s := (mem_effSub.mp hae).2 s
      have hdvd : d ∣ a s := (hdiv (a s)).mp (haeq ▸ (mem_effSub.mp hae).1)
      obtain ⟨q, hq⟩ := hdvd
      have hq0 : 0 ≤ q := by
        by_contra hqn
        push_neg at hqn
        have hlt : a s < 0 := by
          rw [hq]
          exact mul_neg_of_pos_of_neg hd hqn
        omega
      refine ⟨q.toNat, Subtype.ext (Subtype.ext ?_)⟩
      show (single s ((q.toNat : ℤ) * d) : S →₀ ℤ) = a
      rw [Int.toNat_of_nonneg hq0, haeq, hq, mul_comm]

/-- ★★★**`M_p` は monoprime**(しかも `ℤ`-monoprime、すなわち `≃+ ℕ`)——
[FrdI] `Definition 2.4, (i)` の条件 **(b)**。 -/
theorem isMonoprime_Mp_effSub (hQ : IsQCartierSubgroup Γ) (p : Prime (effSub Γ)) :
    IsMonoprime (Mp (effSub Γ) p) := by
  refine ⟨MonoidType.int, ?_⟩
  rw [Mp_eq_suppAtSub hQ p]
  exact nonempty_suppAtSubEquivNat hQ _

/-! ### ★出典の紐付け -/

/-- ★locator —— `Definition 2.4, (i)` の条件 (b) を `Example 6.1` の `Φ(L)` で。 -/
def isMonoprime_Mp_effSub.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L) の M_p は monoprime (Definition 2.4, (i) の (b))",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
