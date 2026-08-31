/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicSparse
import ABC3.Found.GaloisRep.AdicFinsetSum
import ABC3.Found.GaloisRep.DualTate

/-!
# Galois (G6) 第 858 ブロック —— **★★★★★★★★★★`∑_ζ T_{D²f}(ζ) = l⁴s₃(q^l) − s₃(q)`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★これは何か

`T_{D²f}(u) = ∑_n q^n ∑_{d∣n} d³ u^d`（第 856、証明済み）を `μ_l∖{1}` 上で足すと

    `∑_{i≠0} T_{D²f}(ζ^i) = ∑_n q^n ∑_{d∣n} d³ (l[l∣d] − 1)`
                          `= ∑_n q^n (l⁴σ₃(n/l)[l∣n] − σ₃(n))`
                          `= l⁴·s₃(q^l) − s₃(q)`

★★最後の段が `adicSum_sparse`（第 857）である——**`q` の級数を `q^l` の級数に
付け替える**のがここであり、これが `σ₁(q) → σ₁(q^l)` を生む機構の `σ₃` 版である。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★★★★指標和を通す -/

open Finset in
/-- ★★★★★★**`∑_{i≠0} T_{D²f}(ζ^i)` の係数形**。 -/
theorem sum_mu_d2xtail [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ i) q hq
      = adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
            (d : R) ^ 3 * ((if l ∣ d then (l : R) else 0) - 1))
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  rw [Finset.sum_congr rfl (fun i _ => tateD2Xtail_eq_divisorSum (ζ ^ i) q hq)]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => q ^ n * ∑ d ∈ n.divisors, (d : R) ^ 3 * (ζ ^ i) ^ d)
      (fun i n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n))]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [← Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun d _ => ?_)
  rw [← Finset.mul_sum]
  congr 1
  rw [← sum_mu_pow_erase_zero hl hζ d]
  exact Finset.sum_congr rfl (fun i _ => by rw [← pow_mul])

/-- ★★★★**係数を `σ₃` の言葉に直す**。 -/
theorem d2x_coeff_sigma {l : ℕ} (hl : l.Prime) (n : ℕ) :
    ∑ d ∈ n.divisors, (d : R) ^ 3 * ((if l ∣ d then (l : R) else 0) - 1)
      = (l : R) ^ 4 * (if l ∣ n then ((∑ e ∈ (n / l).divisors, e ^ 3 : ℕ) : R) else 0)
        - ((∑ d ∈ n.divisors, d ^ 3 : ℕ) : R) := by
  classical
  have hsplit : ∑ d ∈ n.divisors, (d : R) ^ 3 * ((if l ∣ d then (l : R) else 0) - 1)
      = (l : R) * (∑ d ∈ n.divisors.filter (fun d => l ∣ d), (d : R) ^ 3)
        - ∑ d ∈ n.divisors, (d : R) ^ 3 := by
    rw [Finset.mul_sum, Finset.sum_filter, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun d _ => ?_)
    by_cases hd : l ∣ d
    · rw [if_pos hd, if_pos hd]
      ring
    · rw [if_neg hd, if_neg hd]
      ring
  rw [hsplit]
  congr 1
  · by_cases hn : n = 0
    · subst hn
      simp
    · have h := sum_divisors_dvd_pow 3 l n hl hn
      have hcast : (∑ d ∈ n.divisors.filter (fun d => l ∣ d), (d : R) ^ 3)
          = ((∑ d ∈ n.divisors.filter (fun d => l ∣ d), d ^ 3 : ℕ) : R) := by
        push_cast
        rfl
      rw [hcast, h]
      by_cases hln : l ∣ n
      · rw [if_pos hln, if_pos hln]
        push_cast
        ring
      · rw [if_neg hln, if_neg hln]
        simp
  · push_cast
    rfl

/-! ## ★★★★★★★★`q^l` の級数への付け替え -/

open Finset in
/-- ★★★★★★★★★★**`∑_{i≠0} T_{D²f}(ζ^i) = l⁴·s₃(q^l) − s₃(q)`**。 -/
theorem sum_mu_d2xtail_sigma [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) :
    ∑ i ∈ (range l).erase 0, tateD2Xtail (ζ ^ i) q hq
      = (l : R) ^ 4 * evalAdic (sigmaSeries 3) (q ^ l) hql
        - evalAdic (sigmaSeries 3) q hq := by
  classical
  set A : ℕ → R := fun n => q ^ n *
    ((l : R) ^ 4 * (if l ∣ n then ((∑ e ∈ (n / l).divisors, e ^ 3 : ℕ) : R) else 0)) with hA
  set B : ℕ → R := fun n => q ^ n * ((∑ d ∈ n.divisors, d ^ 3 : ℕ) : R) with hB
  have hAmem : ∀ n, A n ∈ I ^ n := fun n =>
    Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  have hBmem : ∀ n, B n ∈ I ^ n := fun n =>
    Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)
  rw [sum_mu_d2xtail hl hζ q hq]
  have hstep : adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
        (d : R) ^ 3 * ((if l ∣ d then (l : R) else 0) - 1))
        (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n))
      = adicSum (fun n => A n - B n) (fun n => Submodule.sub_mem _ (hAmem n) (hBmem n)) := by
    refine adicSum_congr _ _ (fun n => ?_)
    rw [d2x_coeff_sigma hl n, hA, hB]
    ring
  rw [hstep, adicSum_sub A B hAmem hBmem]
  congr 1
  · -- 疎な級数の付け替え
    have hb : ∀ m : ℕ, (l : R) ^ 4 * ((q ^ l) ^ m * ((∑ e ∈ m.divisors, e ^ 3 : ℕ) : R))
        ∈ I ^ m := by
      intro m
      refine Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ ?_)
      rw [← pow_mul]
      exact Ideal.pow_le_pow_right (Nat.le_mul_of_pos_left m hl.pos)
        (Ideal.pow_mem_pow hq (l * m))
    have hzero : ∀ n, ¬ l ∣ n → A n = 0 := by
      intro n hn
      rw [hA]
      simp [if_neg hn]
    have hshift : ∀ m : ℕ, A (l * m)
        = (l : R) ^ 4 * ((q ^ l) ^ m * ((∑ e ∈ m.divisors, e ^ 3 : ℕ) : R)) := by
      intro m
      rw [hA]
      have hdvd : l ∣ l * m := Dvd.intro m rfl
      have hdiv : l * m / l = m := Nat.mul_div_cancel_left _ hl.pos
      simp only [if_pos hdvd, hdiv]
      rw [← pow_mul]
      ring
    rw [adicSum_sparse l hl.pos A _ hAmem hb hzero hshift,
      adicSum_smul ((l : R) ^ 4) _ (fun m => Ideal.mul_mem_right _ _
        (Ideal.pow_le_pow_right (Nat.le_mul_of_pos_left m hl.pos)
          (by rw [← pow_mul]; exact Ideal.pow_mem_pow hq (l * m))))]
    congr 1
    rw [evalAdic_eq_adicSum]
    refine adicSum_congr _ _ (fun n => ?_)
    rw [coeff_sigmaSeries]
    by_cases hn : n = 0
    · subst hn
      simp
    · rw [if_neg hn]
      rw [ArithmeticFunction.sigma_apply]
      push_cast
      ring
  · rw [evalAdic_eq_adicSum]
    refine adicSum_congr _ _ (fun n => ?_)
    rw [hB, coeff_sigmaSeries]
    by_cases hn : n = 0
    · subst hn
      simp
    · rw [if_neg hn, ArithmeticFunction.sigma_apply]
      push_cast
      ring

/-! ## ★出典の紐付け(`.src`) -/

def sum_mu_d2xtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} T_{D²f}(ζ^i) の係数形。★l は素数)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_d2xtail_sigma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{i≠0} T_{D²f}(ζ^i) = l⁴s₃(q^l) − s₃(q)。★l は素数)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
