/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.IdealCountSplit
import ABC3.Found.NumberField.PrimeDensity
import Mathlib.NumberTheory.NumberField.Discriminant.Different
import Mathlib.RingTheory.RamificationInertia.Ramification

/-!
# 完全分解する素数の密度は `1/[L:ℚ]`(鎖 `cheb` の `cheb-split-density`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★★何をするか —— **類体論を使わない Chebotarev(完全分解の場合)**

`L/ℚ` が Galois のとき、`L` で完全分解する有理素数の Dirichlet 密度は `1/[L:ℚ]`。

★★これは Artin の相互法則も Hecke L 関数も**使わない**。使うのは

* `PrimeDensity.lean`: `Σ_p a_L(p)p^{-s} / log(1/(s−1)) → 1`(解析側)
* `IdealCountSplit.lean`: `a_L(p)` を分解の様子で読む(代数側)
* mathlib: `NumberField.not_dvd_discr_iff_forall_liesOver`(分岐は `discr` を割る素数だけ)

の 3 つだけである。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `ramificationIdx_eq_one_of_not_dvd_discr` | `p ∤ discr L` なら `p` の上の素はすべて不分岐 |
| `finite_ramSet` | ★分岐しうる素数は有限個(`discr L ≠ 0` の約数) |
| `SplQ` | `L` で完全分解する有理素数の集合 |
| `tsum_compl_le` | ★★完全分解しない素数の寄与は**一様有界**(`d·#ram`) |
| `tendsto_splSum_div_log` | `Σ_{p ∈ Spl} a_L(p)p^{-s} / log(1/(s−1)) → 1` |
| `tendsto_splQ_div_log` | ★★★★★★**密度は `1/[L:ℚ]`** |

★★要点は「不分岐なら `a_L(p)` は `0` か `[L:ℚ]` しかない」——
分岐する素数は有限個なので、その寄与は `s → 1+` で消える。
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal Filter Topology
open scoped _root_.NumberField

/-! ## ★1. 分岐する素数は有限個 -/

section Ram

variable {L : Type*} [Field L] [NumberField L]

theorem ne_bot_of_liesOver_span {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L))
    (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) : P ≠ ⊥ := by
  intro h
  refine span_intCast_ne_bot hp ?_
  rw [hLO.over, h]
  exact Ideal.comap_bot_of_injective _ (FaithfulSMul.algebraMap_injective ℤ (𝓞 L))

/-- ★★**`p ∤ discr L` なら `p` の上の素はすべて不分岐**。

★mathlib の `not_dvd_discr_iff_forall_liesOver` と
`Ideal.ramificationIdx'_eq_one` を繋ぐだけ。 -/
theorem ramificationIdx_eq_one_of_not_dvd_discr {p : ℕ} (hp : p.Prime)
    (hdvd : ¬ ((p : ℤ) ∣ discr L)) (P : Ideal (𝓞 L)) (hP : P.IsPrime)
    (hLO : P.LiesOver (Ideal.span {(p : ℤ)})) :
    ramificationIdx (Ideal.span {(p : ℤ)}) P = 1 := by
  haveI : (Ideal.span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  haveI := hP
  haveI := hLO
  haveI : P.IsMaximal := hP.isMaximal (ne_bot_of_liesOver_span hp P hLO)
  haveI : Algebra.IsUnramifiedAt ℤ P :=
    (not_dvd_discr_iff_forall_liesOver L (𝓞 L) (Nat.prime_iff_prime_int.mp hp)).mp hdvd P
      ‹P.IsMaximal› hLO
  rw [Ideal.ramificationIdx_eq_ramificationIdx' _ P (span_intCast_ne_bot hp)]
  exact Ideal.ramificationIdx'_eq_one P ℤ

end Ram

variable (L : Type*) [Field L] [NumberField L]

/-- ★分岐しうる素数の集合(`discr L` を割る素数)。 -/
def ramSet : Set ℕ := {p | p.Prime ∧ ((p : ℤ) ∣ discr L)}

/-- ★★**分岐しうる素数は有限個** —— `discr L ≠ 0` の約数だから。 -/
theorem finite_ramSet : (ramSet L).Finite := by
  refine Set.Finite.subset (Set.finite_Icc 0 (discr L).natAbs) ?_
  rintro p ⟨hp, hdvd⟩
  refine ⟨Nat.zero_le _, ?_⟩
  have h1 : p ∣ (discr L).natAbs := Int.natAbs_dvd_natAbs.mpr (by simpa using hdvd)
  exact Nat.le_of_dvd (Int.natAbs_pos.mpr (discr_ne_zero L)) h1

/-! ## ★2. 完全分解する素数の集合 -/

/-- ★**`L` で完全分解する有理素数の集合**。 -/
def SplQ : Set Nat.Primes :=
  {p | (primesOver (Ideal.span {((p : ℕ) : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L}

/-- ★例外項(分岐する素数の寄与)の一様定数。 -/
noncomputable def splConst : ℝ := (Module.finrank ℚ L : ℝ) * (finite_ramSet L).toFinset.card

variable {L}

theorem splConst_nonneg : 0 ≤ splConst L := by
  rw [splConst]; positivity

/-- ★★**不分岐で完全分解しなければ `a_L(p) = 0`**。 -/
theorem zetaSummandR_eq_zero_of_notMem [IsGalois ℚ L] {s : ℝ} {p : Nat.Primes}
    (hnot : p ∉ SplQ L) (hram : (p : ℕ) ∉ ramSet L) :
    zetaSummandR L s (p : ℕ) = 0 := by
  have hp : (p : ℕ).Prime := p.2
  have hdvd : ¬ (((p : ℕ) : ℤ) ∣ discr L) := fun h => hram ⟨hp, h⟩
  have h0 : idealCount L (p : ℕ) = 0 :=
    idealCount_eq_zero_of_unramified_of_not_split hp
      (fun P hP hLO => ramificationIdx_eq_one_of_not_dvd_discr hp hdvd P hP hLO) hnot
  rw [zetaSummandR, h0]
  simp

theorem zetaSummandR_le_finrank {s : ℝ} (hs : 1 ≤ s) {p : Nat.Primes} :
    zetaSummandR L s (p : ℕ) ≤ (Module.finrank ℚ L : ℝ) := by
  have hp : (p : ℕ).Prime := p.2
  have h1 : zetaSummandR L s (p : ℕ) ≤ (idealCount L (p : ℕ) : ℝ) / ((p : ℕ) : ℝ) :=
    zetaSummandR_le_div L hs (le_of_lt hp.one_lt)
  have h2 : (idealCount L (p : ℕ) : ℝ) ≤ (Module.finrank ℚ L : ℝ) := by
    exact_mod_cast idealCount_le_finrank hp
  have h3 : (1 : ℝ) ≤ ((p : ℕ) : ℝ) := by exact_mod_cast hp.one_lt.le
  calc zetaSummandR L s (p : ℕ) ≤ (idealCount L (p : ℕ) : ℝ) / ((p : ℕ) : ℝ) := h1
    _ ≤ (idealCount L (p : ℕ) : ℝ) / 1 :=
        div_le_div_of_nonneg_left (Nat.cast_nonneg _) zero_lt_one h3
    _ = (idealCount L (p : ℕ) : ℝ) := div_one _
    _ ≤ (Module.finrank ℚ L : ℝ) := h2

/-- ★★★**完全分解しない素数の寄与は一様有界** —— 非ゼロなのは分岐する有限個だけ。 -/
theorem tsum_compl_le [IsGalois ℚ L] {s : ℝ} (hs : 1 < s) :
    ∑' p : ((SplQ L)ᶜ : Set Nat.Primes), zetaSummandR L s (p : ℕ) ≤ splConst L := by
  classical
  have hsummable : Summable (fun p : ((SplQ L)ᶜ : Set Nat.Primes) =>
      zetaSummandR L s ((p : Nat.Primes) : ℕ)) :=
    (summable_zetaSummandR_primes L hs).subtype _
  refine hsummable.tsum_le_of_sum_le (fun u => ?_)
  set T : Finset ((SplQ L)ᶜ : Set Nat.Primes) :=
    u.filter (fun i => ((i : Nat.Primes) : ℕ) ∈ ramSet L) with hT
  have hsplit : ∑ i ∈ u, zetaSummandR L s ((i : Nat.Primes) : ℕ)
      = ∑ i ∈ T, zetaSummandR L s ((i : Nat.Primes) : ℕ) := by
    rw [hT, Finset.sum_filter]
    refine Finset.sum_congr rfl (fun i _ => ?_)
    by_cases hi : ((i : Nat.Primes) : ℕ) ∈ ramSet L
    · rw [if_pos hi]
    · rw [if_neg hi]
      exact zetaSummandR_eq_zero_of_notMem i.2 hi
  have hcard : T.card ≤ (finite_ramSet L).toFinset.card := by
    refine Finset.card_le_card_of_injOn (fun i => ((i : Nat.Primes) : ℕ)) ?_ ?_
    · intro i hi
      simp only [hT, Finset.coe_filter, Set.mem_setOf_eq] at hi
      simpa using hi.2
    · intro i _ j _ hij
      exact Subtype.ext (Subtype.ext hij)
  calc ∑ i ∈ u, zetaSummandR L s ((i : Nat.Primes) : ℕ)
      = ∑ i ∈ T, zetaSummandR L s ((i : Nat.Primes) : ℕ) := hsplit
    _ ≤ ∑ _i ∈ T, (Module.finrank ℚ L : ℝ) :=
        Finset.sum_le_sum (fun i _ => zetaSummandR_le_finrank (le_of_lt hs))
    _ = (T.card : ℝ) * (Module.finrank ℚ L : ℝ) := by rw [Finset.sum_const, nsmul_eq_mul]
    _ ≤ ((finite_ramSet L).toFinset.card : ℝ) * (Module.finrank ℚ L : ℝ) :=
        mul_le_mul_of_nonneg_right (by exact_mod_cast hcard) (Nat.cast_nonneg _)
    _ = splConst L := by rw [splConst]; ring

theorem abs_tsum_sub_splSum_le [IsGalois ℚ L] {s : ℝ} (hs : 1 < s) :
    |(∑' p : Nat.Primes, zetaSummandR L s (p : ℕ))
      - (∑' p : SplQ L, zetaSummandR L s ((p : Nat.Primes) : ℕ))| ≤ splConst L := by
  have hsum := summable_zetaSummandR_primes L hs
  have hsplit := hsum.tsum_subtype_add_tsum_subtype_compl (SplQ L)
  have heq : (∑' p : Nat.Primes, zetaSummandR L s (p : ℕ))
      - (∑' p : SplQ L, zetaSummandR L s ((p : Nat.Primes) : ℕ))
      = ∑' p : ((SplQ L)ᶜ : Set Nat.Primes), zetaSummandR L s ((p : Nat.Primes) : ℕ) := by
    rw [← hsplit]; ring
  rw [heq, abs_of_nonneg (tsum_nonneg (fun p => zetaSummandR_nonneg L _))]
  exact tsum_compl_le hs

/-! ## ★3. 密度 -/

/-- ★★★★★**完全分解する素数の和も `log(1/(s−1))` と同値**。 -/
theorem tendsto_splSum_div_log [IsGalois ℚ L] :
    Tendsto (fun s : ℝ ↦
        (∑' p : SplQ L, zetaSummandR L s ((p : Nat.Primes) : ℕ)) / Real.log (1 / (s - 1)))
      (𝓝[>] 1) (𝓝 1) := by
  set D : ℝ → ℝ := fun s => Real.log (1 / (s - 1)) with hD
  set A : ℝ → ℝ := fun s => ∑' p : Nat.Primes, zetaSummandR L s (p : ℕ) with hA
  set B : ℝ → ℝ := fun s => ∑' p : SplQ L, zetaSummandR L s ((p : Nat.Primes) : ℕ) with hB
  have hDtop : Tendsto D (𝓝[>] (1:ℝ)) atTop := tendsto_logInv_atTop
  have h1 : Tendsto (fun s => A s / D s) (𝓝[>] (1:ℝ)) (𝓝 1) :=
    tendsto_tsum_primes_div_log L
  have hconst : Tendsto (fun _ : ℝ => splConst L) (𝓝[>] (1:ℝ)) (𝓝 (splConst L)) :=
    tendsto_const_nhds
  have h2 : Tendsto (fun s => (A s - B s) / D s) (𝓝[>] (1:ℝ)) (𝓝 0) := by
    refine squeeze_zero_norm' ?_ (hconst.div_atTop hDtop)
    filter_upwards [self_mem_nhdsWithin, hDtop.eventually_gt_atTop 0] with s hs hDs
    have hs1 : (1 : ℝ) < s := hs
    rw [Real.norm_eq_abs, abs_div, abs_of_pos hDs]
    exact div_le_div_of_nonneg_right (abs_tsum_sub_splSum_le hs1) (le_of_lt hDs)
  have hcomb : Tendsto (fun s => A s / D s - (A s - B s) / D s) (𝓝[>] (1:ℝ)) (𝓝 (1 - 0)) :=
    h1.sub h2
  rw [sub_zero] at hcomb
  refine hcomb.congr' ?_
  filter_upwards [hDtop.eventually_gt_atTop 0] with s hDs
  have hne : D s ≠ 0 := ne_of_gt hDs
  show A s / D s - (A s - B s) / D s = B s / D s
  field_simp
  ring

theorem zetaSummandR_of_mem_SplQ [IsGalois ℚ L] {s : ℝ} {p : Nat.Primes} (h : p ∈ SplQ L) :
    zetaSummandR L s (p : ℕ) = (Module.finrank ℚ L : ℝ) * ((p : ℕ) : ℝ) ^ (-s) := by
  rw [zetaSummandR, idealCount_eq_finrank_of_splitsCompletely p.2 h]

theorem summable_primes_rpow {s : ℝ} (hs : 1 < s) :
    Summable (fun p : Nat.Primes ↦ ((p : ℕ) : ℝ) ^ (-s)) := by
  have h : Summable (fun n : ℕ ↦ ((n : ℝ) ^ s)⁻¹) := Real.summable_nat_rpow_inv.mpr hs
  refine (h.comp_injective (fun _ _ hij => Subtype.ext hij)).congr fun p ↦ ?_
  show ((((p : Nat.Primes) : ℕ) : ℝ) ^ s)⁻¹ = _
  rw [← Real.rpow_neg (Nat.cast_nonneg _)]

theorem splSum_eq [IsGalois ℚ L] {s : ℝ} :
    (∑' p : SplQ L, zetaSummandR L s ((p : Nat.Primes) : ℕ))
      = (Module.finrank ℚ L : ℝ) * ∑' p : SplQ L, (((p : Nat.Primes) : ℕ) : ℝ) ^ (-s) := by
  rw [← tsum_mul_left]
  exact tsum_congr fun p => zetaSummandR_of_mem_SplQ p.2

/-- ★★★★★★**[cheb-split-density] 完全分解する素数の Dirichlet 密度は `1/[L:ℚ]`**
—— **類体論を使わない**。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★★使うのは `Σ_p a_L(p)p^{-s} ~ log(1/(s−1))`(Dedekind ζ の単純極)と、
「不分岐なら `a_L(p)` は `0` か `[L:ℚ]`」「分岐する素数は有限個」だけである。
Artin の相互法則も Hecke L 関数も要らない。 -/
theorem tendsto_splQ_div_log [IsGalois ℚ L] :
    Tendsto (fun s : ℝ ↦
        (∑' p : SplQ L, (((p : Nat.Primes) : ℕ) : ℝ) ^ (-s)) / Real.log (1 / (s - 1)))
      (𝓝[>] 1) (𝓝 (1 / (Module.finrank ℚ L : ℝ))) := by
  have hd : (Module.finrank ℚ L : ℝ) ≠ 0 := by
    have h0 : 0 < Module.finrank ℚ L := Module.finrank_pos
    positivity
  have h := (tendsto_splSum_div_log (L := L)).const_mul ((Module.finrank ℚ L : ℝ)⁻¹)
  rw [mul_one, ← one_div] at h
  refine h.congr' ?_
  filter_upwards with s
  rw [splSum_eq]
  field_simp

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— 類体論を使わない Chebotarev(完全分解の場合)。 -/
def tendsto_splQ_div_log.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 完全分解する素数の密度は 1/[L:ℚ]",
    sectionId := "frdi-thm-6-4" }

def tendsto_splQ_div_log.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "tendsto_tsum_primes_div_log(解析側)"
      (.inProject "ABC3" "ABC3.Found.NF.tendsto_tsum_primes_div_log") 116,
    .citation "[ABC3]" "idealCount_eq_finrank_of_splitsCompletely(代数側)"
      (.inProject "ABC3" "ABC3.Found.NF.idealCount_eq_finrank_of_splitsCompletely") 116,
    .citation "[mathlib]" "not_dvd_discr_iff_forall_liesOver(分岐は discr の約数だけ)"
      (.inMathlib "NumberField.not_dvd_discr_iff_forall_liesOver") 116 ]

end ABC3.Found.NF
