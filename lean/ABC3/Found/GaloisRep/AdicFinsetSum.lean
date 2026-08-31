/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicSeries
import ABC3.Found.GaloisRep.TateSigma
import ABC3.Found.GaloisRep.TatePair
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

/-! ## ★★★★★★★★★★★★★★★★`w = qζ⁻¹` 側の和 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★
**`∑_ζ tateXterm(qζ⁻¹)` も `ζ` を含まない**。

`tateXpair ζ (qζ⁻¹) q` の `w` 側の項である。 -/
theorem sum_mu_tateXterm_w [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXterm (q * ζ ^ (i * (l - 1)))
      = adicSum (I := I)
          (fun n => (n : R) * q ^ n * ((if l ∣ n then (l : R) else 0) - 1))
          (fun n => Ideal.mul_mem_right _ _
            (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n))) := by
  classical
  have hmem : ∀ i : ℕ, q * ζ ^ (i * (l - 1)) ∈ I := fun i => Ideal.mul_mem_right _ _ hq
  rw [Finset.sum_congr rfl (fun i _ => tateXterm_eq_adicSum (hmem i))]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => (n : R) * (q * ζ ^ (i * (l - 1))) ^ n)
      (fun i n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow (hmem i) n))]
  refine adicSum_congr _ _ (fun n => ?_)
  have hexp : ∀ i : ℕ, (n : R) * (q * ζ ^ (i * (l - 1))) ^ n
      = ((n : R) * q ^ n) * (ζ ^ (i * (l - 1))) ^ n := by
    intro i; rw [mul_pow]; ring
  rw [Finset.sum_congr rfl (fun i _ => hexp i), ← Finset.mul_sum,
    sum_mu_neg_pow hl hζ n, mul_assoc]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**`Y` 側の `w` の項**。 -/
theorem sum_mu_tateYterm_w [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateYterm (q * ζ ^ (i * (l - 1)))
      = adicSum (I := I)
          (fun n => ((n.choose 2 : ℕ) : R) * q ^ n * ((if l ∣ n then (l : R) else 0) - 1))
          (fun n => Ideal.mul_mem_right _ _
            (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n))) := by
  classical
  have hmem : ∀ i : ℕ, q * ζ ^ (i * (l - 1)) ∈ I := fun i => Ideal.mul_mem_right _ _ hq
  rw [Finset.sum_congr rfl (fun i _ => tateYterm_eq_adicSum (hmem i))]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => ((n.choose 2 : ℕ) : R) * (q * ζ ^ (i * (l - 1))) ^ n)
      (fun i n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow (hmem i) n))]
  refine adicSum_congr _ _ (fun n => ?_)
  have hexp : ∀ i : ℕ, ((n.choose 2 : ℕ) : R) * (q * ζ ^ (i * (l - 1))) ^ n
      = (((n.choose 2 : ℕ) : R) * q ^ n) * (ζ ^ (i * (l - 1))) ^ n := by
    intro i; rw [mul_pow]; ring
  rw [Finset.sum_congr rfl (fun i _ => hexp i), ← Finset.mul_sum,
    sum_mu_neg_pow hl hζ n, mul_assoc]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**`w` 側の尾**。 -/
theorem sum_mu_tateXtail_w [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXtail (q * ζ ^ (i * (l - 1))) q hq
      = adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
            (d : R) * (q ^ d * ((if l ∣ d then (l : R) else 0) - 1)))
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  classical
  rw [Finset.sum_congr rfl
    (fun i _ => tateXtail_eq_divisorSum (q * ζ ^ (i * (l - 1))) q hq)]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => q ^ n * ∑ d ∈ n.divisors,
        (d : R) * (q * ζ ^ (i * (l - 1))) ^ d)
      (fun i n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n))]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [← Finset.mul_sum]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun d _ => ?_)
  have hexp : ∀ i : ℕ, (d : R) * (q * ζ ^ (i * (l - 1))) ^ d
      = ((d : R) * q ^ d) * (ζ ^ (i * (l - 1))) ^ d := by
    intro i; rw [mul_pow]; ring
  rw [Finset.sum_congr rfl (fun i _ => hexp i), ← Finset.mul_sum,
    sum_mu_neg_pow hl hζ d, mul_assoc]

/-! ## ★★★★★★★★★★★★★★★★★★★★`tateXpair` を `μ_l` にわたって足す -/

open Finset in
/-- ★`ζ^i` と `qζ^{−i}` の積は `q`——`tateXpair` の対の条件。 -/
theorem mu_pair_mul {l : ℕ} (hl : 0 < l) {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (i : ℕ) :
    ζ ^ i * (q * ζ ^ (i * (l - 1))) = q := by
  rw [← mul_assoc, mul_comm (ζ ^ i) q, mul_assoc, ← pow_add]
  have hkey : i + i * (l - 1) = l * i := by
    cases l with
    | zero => omega
    | succ n => simp only [Nat.add_sub_cancel]; ring
  rw [hkey, pow_mul, hζ.pow_eq_one, one_pow, mul_one]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★
**`∑_{ζ∈μ_l∖{1}} X(ζ, q)` を 4 つの部分に分ける**。

★`tateXpair a w q = (term a + tail a) + (term w + tail w) − 2 s₁(q)` なので、
有限和を分配するだけである。★★右辺の 4 つはすべて
本ファイル・`MuCharSum.lean` で **`ζ`-free に計算済み**である。 -/
theorem sum_mu_tateXpair [IsAdicComplete I R] {l : ℕ} {ζ : R} (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * ζ ^ (i * (l - 1))) q hq
      = ((∑ i ∈ (range l).erase 0, tateXterm (ζ ^ i))
          + (∑ i ∈ (range l).erase 0, tateXtail (ζ ^ i) q hq))
        + ((∑ i ∈ (range l).erase 0, tateXterm (q * ζ ^ (i * (l - 1))))
          + (∑ i ∈ (range l).erase 0, tateXtail (q * ζ ^ (i * (l - 1))) q hq))
        - (((range l).erase 0).card : R) * (2 * evalAdic (sigmaSeries 1) q hq) := by
  classical
  simp only [tateXpair]
  rw [Finset.sum_sub_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
    Finset.sum_add_distrib, Finset.sum_const, nsmul_eq_mul]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★**`Y` 側も同じ分配**。 -/
theorem sum_mu_tateYpair [IsAdicComplete I R] {l : ℕ} {ζ : R} (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateYpair (ζ ^ i) (q * ζ ^ (i * (l - 1))) q hq
      = ((∑ i ∈ (range l).erase 0, tateYterm (ζ ^ i))
          + (∑ i ∈ (range l).erase 0, tateYtail (ζ ^ i) q hq))
        - ((∑ i ∈ (range l).erase 0, tateXterm (q * ζ ^ (i * (l - 1))))
          + (∑ i ∈ (range l).erase 0, tateXtail (q * ζ ^ (i * (l - 1))) q hq))
        - ((∑ i ∈ (range l).erase 0, tateYterm (q * ζ ^ (i * (l - 1))))
          + (∑ i ∈ (range l).erase 0, tateYtail (q * ζ ^ (i * (l - 1))) q hq))
        + (((range l).erase 0).card : R) * evalAdic (sigmaSeries 1) q hq := by
  classical
  simp only [tateYpair]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.sum_sub_distrib,
    Finset.sum_add_distrib, Finset.sum_add_distrib, Finset.sum_add_distrib,
    Finset.sum_const, nsmul_eq_mul]

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**`∑_{ζ∈μ_l∖{1}} X(ζ, q)` は定数項を除いて `ζ` を含まない**。

★★★これが `tateModel_of_quot_mu`（葉 1）の手順 1・2 の到達点である。
右辺の 3 つの `adicSum` はどれも係数に `[l ∣ d]·l − 1` を持ち、
`l ∣ d` の項だけが残る——これが `q` 展開を `q^l` 展開に付け替える。

☆定数項 `∑_ζ tateXterm(ζ)` は `MuCharSum.lean` の
`twelve_mul_sum_mu_ringInverse` で `12·(それ) = −(l²−1)` と定まっている。 -/
theorem sum_mu_tateXpair_eq [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * ζ ^ (i * (l - 1))) q hq
      = (∑ i ∈ (range l).erase 0, tateXterm (ζ ^ i))
        + adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
              (d : R) * ((if l ∣ d then (l : R) else 0) - 1))
            (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n))
        + adicSum (I := I)
            (fun n => (n : R) * q ^ n * ((if l ∣ n then (l : R) else 0) - 1))
            (fun n => Ideal.mul_mem_right _ _
              (Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)))
        + adicSum (I := I) (fun n => q ^ n * ∑ d ∈ n.divisors,
              (d : R) * (q ^ d * ((if l ∣ d then (l : R) else 0) - 1)))
            (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n))
        - ((l - 1 : ℕ) : R) * (2 * evalAdic (sigmaSeries 1) q hq) := by
  classical
  have hcard : ((range l).erase 0).card = l - 1 := by
    rw [Finset.card_erase_of_mem (Finset.mem_range.2 hl.pos), Finset.card_range]
  rw [sum_mu_tateXpair q hq, sum_mu_tateXtail hl hζ q hq,
    sum_mu_tateXterm_w hl hζ q hq, sum_mu_tateXtail_w hl hζ q hq, hcard]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def sum_mu_tateXtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l にわたる Tate の X の尾の和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateYtail.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l にわたる Tate の Y の尾の和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateXterm_w.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(w = qζ⁻¹ 側の X の項の μ_l 和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateYterm_w.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(w 側の Y の項の μ_l 和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateXtail_w.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(w 側の尾の μ_l 和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def mu_pair_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ^i と qζ^{−i} の積は q。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateXpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_ζ X(ζ,q) を 4 つの部分に分ける。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateYpair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_ζ Y(ζ,q) の分配。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_tateXpair_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_ζ X(ζ,q) は定数項を除いて ζ-free。★無条件)",
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
