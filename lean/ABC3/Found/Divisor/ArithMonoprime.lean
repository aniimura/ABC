/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPhiMonoid
import ABC3.Found.Divisor.CartierMonoprime

/-!
# 実係数の `M_p` は monoprime —— 非アルキメデスで `ℕ`、アルキメデスで `ℝ≥0`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.112。

原文 (FrdI p.112):
> plex archimedean valuations with their complex conjugates]; OF for the ring of

## ★★幾何版との違いはここだけ

`Cartier*.lean`(`Example 6.1`)では `M_p ≃+ ℕ` が**つねに**成り立つ ——
`ℤ` の部分群はすべて単項だからである。

★算術(`Example 6.3`)では素点が 2 種類ある:

* **非アルキメデス** `v`: `Γ_v = ℤ·log(N v)` は離散 —— `M_p ≃+ ℕ`(`ℤ`-monoprime)
* **アルキメデス** `v`: `Γ_v = ℝ` —— `M_p ≃+ ℝ≥0`(**`ℝ`-monoprime**)

★★これが原文の `ord(O_v^▷) ≅ ℤ≥0` / `≅ ℝ≥0` の分岐そのものである。
`ℝ` の部分群は「離散(巡回)」か「稠密」の二択だが、**稠密なものは monoprime でない**
(可除でないから `ℝ≥0` と同型にならない)。★そこで
「各素点で離散か全体か」を仮定 `IsLocallyMonoprimeR` として明示に置く ——
`Example 6.3` の `Γ = arithDivGroup L` はこれを満たす。
-/

namespace ABC3.Found.Divisor

open Finsupp ABC3.Found.FrdI

open scoped NNReal

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℝ)}

/-! ## ★1. 素点ごとの局所群 -/

/-- ★`{c : ℝ | c·s ∈ Γ}` —— `ℝ` の部分群。 -/
noncomputable def sGrpAtR (Γ : AddSubgroup (S →₀ ℝ)) (s : S) : AddSubgroup ℝ :=
  Γ.comap (Finsupp.singleAddHom s)

theorem mem_sGrpAtR {s : S} {c : ℝ} : c ∈ sGrpAtR Γ s ↔ single s c ∈ Γ := Iff.rfl

/-- ★★**各素点で局所群が「離散」か「全体」**。

★`Example 6.3` の `Γ = arithDivGroup L` では
非アルキメデス `v` で `ℤ·log(N v)`、アルキメデス `v` で `ℝ` になる。 -/
def IsLocallyMonoprimeR (Γ : AddSubgroup (S →₀ ℝ)) : Prop :=
  ∀ s : S, (∃ d : ℝ, 0 < d ∧ ∀ c : ℝ, single s c ∈ Γ ↔ ∃ n : ℤ, c = (n : ℝ) * d)
    ∨ (∀ c : ℝ, single s c ∈ Γ)

/-- ★局所群が離散か全体なら、各素点で正の元がある。 -/
theorem isGenSubgroupR_of_isLocallyMonoprimeR (hL : IsLocallyMonoprimeR Γ) :
    IsGenSubgroupR Γ := by
  intro s
  rcases hL s with ⟨d, hd, hspec⟩ | hfull
  · exact ⟨d, hd, (hspec d).mpr ⟨1, by ring⟩⟩
  · exact ⟨1, one_pos, hfull 1⟩

/-! ## ★2. `M_p` の同定 -/

/-- ★`s` に台を持つ `Φ` の元がなす部分単系。 -/
def suppAtSubR (Γ : AddSubgroup (S →₀ ℝ)) (s : S) : AddSubmonoid (effR Γ) where
  carrier := {a : effR Γ | ((a : effR Γ) : S →₀ ℝ).support ⊆ {s}}
  add_mem' {a b} ha hb := by
    classical
    have hsub : (((a + b : effR Γ)) : S →₀ ℝ).support
        ⊆ ((a : effR Γ) : S →₀ ℝ).support ∪ ((b : effR Γ) : S →₀ ℝ).support := by
      show (((a : effR Γ) : S →₀ ℝ) + ((b : effR Γ) : S →₀ ℝ)).support ⊆ _
      exact Finsupp.support_add
    exact hsub.trans (Finset.union_subset ha hb)
  zero_mem' := by
    show ((0 : effR Γ) : S →₀ ℝ).support ⊆ {s}
    simp

theorem mem_suppAtSubR {s : S} {a : effR Γ} :
    a ∈ suppAtSubR Γ s ↔ ((a : effR Γ) : S →₀ ℝ).support ⊆ {s} := Iff.rfl

/-- ★`primeCarrier` は「台がちょうど `{s}`」。 -/
theorem mem_primeCarrier_iff_R (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) (a : effR Γ) :
    a ∈ primeCarrier (effR Γ) p
      ↔ ((a : effR Γ) : S →₀ ℝ).support = {effRPrimeEquiv hG p} := by
  constructor
  · rintro ⟨ha, hp⟩
    obtain ⟨t, ht⟩ := support_singleton_of_isPrimaryElt_R hG ha
    have hpt : effRPrimeEquiv hG p = t := by
      rw [← hp]
      exact primePtR_toPrime hG ha ht
    rw [hpt, ht]
  · intro h
    have hprim : IsPrimaryElt a := isPrimaryElt_of_support_singleton_R h
    refine ⟨hprim, ?_⟩
    have h1 : effRPrimeEquiv hG (toPrime _ a hprim) = effRPrimeEquiv hG p :=
      primePtR_toPrime hG hprim h
    exact (effRPrimeEquiv hG).injective h1

/-- ★★**`M_p` はちょうど「台が `{s}` に入る元」**。 -/
theorem Mp_eq_suppAtSubR (hG : IsGenSubgroupR Γ) (p : Prime (effR Γ)) :
    Mp (effR Γ) p = suppAtSubR Γ (effRPrimeEquiv hG p) := by
  refine le_antisymm (AddSubmonoid.closure_le.mpr ?_) ?_
  · intro a ha
    exact mem_suppAtSubR.mpr ((mem_primeCarrier_iff_R hG p a).mp ha).le
  · intro a ha
    classical
    by_cases h0 : a = 0
    · rw [h0]; exact AddSubmonoid.zero_mem _
    · have hane : ((a : effR Γ) : S →₀ ℝ) ≠ 0 := fun hz => h0 ((effR_eq_zero_iff a).mpr hz)
      obtain ⟨t, ht⟩ := Finsupp.support_nonempty_iff.mpr hane
      have hsub := mem_suppAtSubR.mp ha
      have hts : t = effRPrimeEquiv hG p := Finset.mem_singleton.mp (hsub ht)
      have heq : ((a : effR Γ) : S →₀ ℝ).support = {effRPrimeEquiv hG p} :=
        Finset.Subset.antisymm hsub (by rw [← hts]; exact Finset.singleton_subset_iff.mpr ht)
      exact mem_Mp_of_mem_primeCarrier ((mem_primeCarrier_iff_R hG p a).mpr heq)

/-- ★台が `{s}` に入る元は `single s (a s)` である。 -/
theorem eq_single_of_support_subset_R {x : S →₀ ℝ} {s : S} (h : x.support ⊆ {s}) :
    x = single s (x s) := by
  classical
  ext t
  rcases eq_or_ne s t with rfl | hst
  · simp
  · have hns : t ∉ x.support := fun ht => hst (Finset.mem_singleton.mp (h ht)).symm
    rw [Finsupp.notMem_support_iff.mp hns, Finsupp.single_apply, if_neg hst]

/-! ## ★3. 離散な素点 —— `M_p ≃+ ℕ` -/

/-- ★★**局所群が `ℤ·d` なら `{a ∈ Φ | supp a ⊆ {s}} ≃+ ℕ`**(非アルキメデスの場合)。 -/
theorem nonempty_suppAtSubR_equiv_nat {s : S} {d : ℝ} (hd : 0 < d)
    (hspec : ∀ c : ℝ, single s c ∈ Γ ↔ ∃ n : ℤ, c = (n : ℝ) * d) :
    Nonempty (suppAtSubR Γ s ≃+ ℕ) := by
  classical
  have hmem : ∀ k : ℕ, (single s ((k : ℝ) * d) : S →₀ ℝ) ∈ effR Γ := by
    intro k
    refine mem_effR.mpr ⟨(hspec _).mpr ⟨(k : ℤ), by push_cast; ring⟩, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simp only [Finsupp.single_eq_same]
      positivity
    · simp [hst]
  have hsupp : ∀ k : ℕ, (⟨_, hmem k⟩ : effR Γ) ∈ suppAtSubR Γ s :=
    fun k => mem_suppAtSubR.mpr Finsupp.support_single_subset
  refine ⟨(AddEquiv.ofBijective
    (⟨⟨fun k : ℕ => (⟨⟨_, hmem k⟩, hsupp k⟩ : suppAtSubR Γ s), ?_⟩, ?_⟩ : ℕ →+ suppAtSubR Γ s)
    ?_).symm⟩
  · refine Subtype.ext (Subtype.ext ?_)
    show (single s (((0 : ℕ) : ℝ) * d) : S →₀ ℝ) = 0
    simp
  · intro k l
    refine Subtype.ext (Subtype.ext ?_)
    show (single s (((k + l : ℕ) : ℝ) * d) : S →₀ ℝ)
      = single s ((k : ℝ) * d) + single s ((l : ℝ) * d)
    rw [← Finsupp.single_add]
    congr 1
    push_cast
    ring
  · constructor
    · intro k l hkl
      have h1 : (single s ((k : ℝ) * d) : S →₀ ℝ) = single s ((l : ℝ) * d) :=
        congrArg (fun t : suppAtSubR Γ s => ((t : effR Γ) : S →₀ ℝ)) hkl
      have h2 := congrArg (fun t : S →₀ ℝ => t s) h1
      simp only [Finsupp.single_eq_same] at h2
      have h3 : (k : ℝ) = (l : ℝ) := mul_right_cancel₀ (ne_of_gt hd) h2
      exact_mod_cast h3
    · rintro ⟨⟨a, hae⟩, hsu⟩
      have hsub : a.support ⊆ {s} := mem_suppAtSubR.mp hsu
      have haeq : a = single s (a s) := eq_single_of_support_subset_R hsub
      have hann : 0 ≤ a s := (mem_effR.mp hae).2 s
      obtain ⟨n, hn⟩ := (hspec (a s)).mp (haeq ▸ (mem_effR.mp hae).1)
      have hn0 : 0 ≤ n := by
        by_contra hneg
        push_neg at hneg
        have : (n : ℝ) < 0 := by exact_mod_cast hneg
        nlinarith
      refine ⟨n.toNat, Subtype.ext (Subtype.ext ?_)⟩
      show (single s ((n.toNat : ℝ) * d) : S →₀ ℝ) = a
      have hcast : ((n.toNat : ℤ) : ℝ) = (n : ℝ) := by
        rw [Int.toNat_of_nonneg hn0]
      rw [haeq, hn]
      congr 1
      push_cast at hcast ⊢
      rw [hcast]

/-! ## ★4. 稠密(全体)な素点 —— `M_p ≃+ ℝ≥0` -/

/-- ★★★**局所群が `ℝ` 全体なら `{a ∈ Φ | supp a ⊆ {s}} ≃+ ℝ≥0`**(アルキメデスの場合)。

★★これが `Example 6.1`(幾何)には現れない、`Example 6.3` に固有の分岐である。 -/
theorem nonempty_suppAtSubR_equiv_nnreal {s : S} (hfull : ∀ c : ℝ, single s c ∈ Γ) :
    Nonempty (suppAtSubR Γ s ≃+ ℝ≥0) := by
  classical
  have hmem : ∀ r : ℝ≥0, (single s (r : ℝ) : S →₀ ℝ) ∈ effR Γ := by
    intro r
    refine mem_effR.mpr ⟨hfull _, fun t => ?_⟩
    rcases eq_or_ne s t with rfl | hst
    · simpa using r.2
    · simp [hst]
  have hsupp : ∀ r : ℝ≥0, (⟨_, hmem r⟩ : effR Γ) ∈ suppAtSubR Γ s :=
    fun r => mem_suppAtSubR.mpr Finsupp.support_single_subset
  refine ⟨(AddEquiv.ofBijective
    (⟨⟨fun r : ℝ≥0 => (⟨⟨_, hmem r⟩, hsupp r⟩ : suppAtSubR Γ s), ?_⟩, ?_⟩ :
      ℝ≥0 →+ suppAtSubR Γ s) ?_).symm⟩
  · refine Subtype.ext (Subtype.ext ?_)
    show (single s (((0 : ℝ≥0) : ℝ)) : S →₀ ℝ) = 0
    simp
  · intro r t
    refine Subtype.ext (Subtype.ext ?_)
    show (single s (((r + t : ℝ≥0) : ℝ)) : S →₀ ℝ)
      = single s ((r : ℝ)) + single s ((t : ℝ))
    rw [← Finsupp.single_add, NNReal.coe_add]
  · constructor
    · intro r t hrt
      have h1 : (single s ((r : ℝ)) : S →₀ ℝ) = single s ((t : ℝ)) :=
        congrArg (fun z : suppAtSubR Γ s => ((z : effR Γ) : S →₀ ℝ)) hrt
      have h2 := congrArg (fun z : S →₀ ℝ => z s) h1
      simp only [Finsupp.single_eq_same] at h2
      exact NNReal.coe_injective h2
    · rintro ⟨⟨a, hae⟩, hsu⟩
      have hsub : a.support ⊆ {s} := mem_suppAtSubR.mp hsu
      have haeq : a = single s (a s) := eq_single_of_support_subset_R hsub
      have hann : 0 ≤ a s := (mem_effR.mp hae).2 s
      refine ⟨⟨a s, hann⟩, Subtype.ext (Subtype.ext ?_)⟩
      show (single s ((⟨a s, hann⟩ : ℝ≥0) : ℝ) : S →₀ ℝ) = a
      exact haeq.symm

/-! ## ★5. まとめ -/

/-- ★★★★**`M_p` は monoprime**(実係数版)——
非アルキメデスで `ℤ`-monoprime、アルキメデスで **`ℝ`-monoprime**。 -/
theorem isMonoprime_Mp_effR (hL : IsLocallyMonoprimeR Γ) (p : Prime (effR Γ)) :
    IsMonoprime (Mp (effR Γ) p) := by
  have hG := isGenSubgroupR_of_isLocallyMonoprimeR hL
  rw [Mp_eq_suppAtSubR hG p]
  rcases hL (effRPrimeEquiv hG p) with ⟨d, hd, hspec⟩ | hfull
  · exact ⟨MonoidType.int, nonempty_suppAtSubR_equiv_nat hd hspec⟩
  · exact ⟨MonoidType.real, nonempty_suppAtSubR_equiv_nnreal hfull⟩

/-- ★★**アルキメデス素点では `M_p` は `ℝ`-monoprime**。 -/
theorem isLambdaMonoprime_real_Mp_effR (hL : IsLocallyMonoprimeR Γ)
    (p : Prime (effR Γ))
    (hfull : ∀ c : ℝ, single (effRPrimeEquiv (isGenSubgroupR_of_isLocallyMonoprimeR hL) p) c ∈ Γ) :
    IsLambdaMonoprime (Mp (effR Γ) p) MonoidType.real := by
  rw [Mp_eq_suppAtSubR (isGenSubgroupR_of_isLocallyMonoprimeR hL) p]
  exact nonempty_suppAtSubR_equiv_nnreal hfull

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Example 6.3` の `ord(O_v^▷) ≅ ℤ≥0 / ℝ≥0`(`M_p` の monoprime 性)。 -/
def isMonoprime_Mp_effR.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 112,
    item := "Example 6.3 — Φ(L) の M_p は monoprime(非アルキメデスで ℕ、アルキメデスで ℝ≥0)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
