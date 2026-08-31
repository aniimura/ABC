/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicSeries
import ABC3.Found.GaloisRep.TateSigma
import ABC3.Found.GaloisRep.MuCharSum
import ABC3.Meta.Claim

/-!
# `adicSum` は有限和と可換（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`tateModel_of_quot_mu`（葉 1）では `μ_l` の点にわたる **Vélu の和**を取る。
その各項は `adicSum`（`I` 進和）で書かれているので、

    `∑_{i ∈ s} adicSum (a i) = adicSum (fun n => ∑_{i ∈ s} a i n)`

がないと `Found/GaloisRep/MuCharSum.lean` の指標和が使えない。★本ファイルはそれを与える。
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★`0` の `I` 進和は `0`。 -/
theorem adicSum_zero [IsAdicComplete I R] :
    adicSum (I := I) (fun _ : ℕ => (0 : R)) (fun n => Submodule.zero_mem _) = 0 := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  simp [partialSum]

/-- ★★★★★★★★★★★★★★★★**`adicSum` は有限和と可換**。

★これで `μ_l` の点にわたる和を `q` の係数ごとの和に直せる。 -/
theorem adicSum_finsetSum [IsAdicComplete I R] {ι : Type*} [DecidableEq ι]
    (s : Finset ι) (a : ι → ℕ → R) (ha : ∀ i, ∀ n, a i n ∈ I ^ n) :
    adicSum (I := I) (fun n => ∑ i ∈ s, a i n)
        (fun n => Submodule.sum_mem _ (fun i _ => ha i n))
      = ∑ i ∈ s, adicSum (a i) (ha i) := by
  refine adicSum_unique _ _ _ (fun n => ?_)
  have hps : partialSum (fun n => ∑ i ∈ s, a i n) n
      = ∑ i ∈ s, partialSum (a i) n := by
    simp only [partialSum]
    exact Finset.sum_comm
  rw [hps, SModEq.sub_mem, ← Finset.sum_sub_distrib]
  refine Submodule.sum_mem _ (fun i _ => ?_)
  exact SModEq.sub_mem.1 (adicSum_spec (a i) (ha i) n)

/-! ## ★★★★★★★★★★★★★★★★`μ_l` にわたる Tate の尾の和 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`∑_{ζ ∈ μ_l∖{1}} tateXtail(ζ)` は `ζ` を含まない**。

    `∑_ζ tateXtail(ζ, q) = adicSum (n ↦ q^n · ∑_{d ∣ n} d·([l ∣ d]·l − 1))`

★★★これが `tateModel_of_quot_mu`（葉 1）の手順 2 を
**実際の Tate 級数に適用した形**である。 -/
theorem sum_mu_tateXtail [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXtail (ζ ^ i) q hq
      = adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
            (d : R) * ((if l ∣ d then (l : R) else 0) - 1))
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  rw [Finset.sum_congr rfl (fun i _ => tateXtail_eq_divisorSum (ζ ^ i) q hq)]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => q ^ n * ∑ d ∈ n.divisors, (d : R) * (ζ ^ i) ^ d)
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

open Finset in
/-- ★★★★★★★★★★★★★★★★★★**`Y` 側も同じ**。 -/
theorem sum_mu_tateYtail [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateYtail (ζ ^ i) q hq
      = adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
            ((d.choose 2 : ℕ) : R) * ((if l ∣ d then (l : R) else 0) - 1))
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  rw [Finset.sum_congr rfl (fun i _ => tateYtail_eq_divisorSum (ζ ^ i) q hq)]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => q ^ n * ∑ d ∈ n.divisors, ((d.choose 2 : ℕ) : R) * (ζ ^ i) ^ d)
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

/-! ## ★出典の紐付け(`.src`) -/

def sum_mu_tateXtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l にわたる Tate の X の尾の和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateYtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l にわたる Tate の Y の尾の和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def adicSum_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(0 の I 進和は 0。★無条件)",
    sectionId := "genell-lemma-3-2" }

def adicSum_finsetSum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(adicSum は有限和と可換。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
