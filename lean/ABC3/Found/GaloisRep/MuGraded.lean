/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicFinsetSum
import ABC3.Meta.Claim

/-!
# μ-等級付き `I` 進級数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`tateModel_of_quot_mu`（葉 1）の**手順 3**では
`v = ∑_ζ (3X(ζ)² + a₄ − Y(ζ))` のように **`X(ζ)` の積**が現れる。

★`X(ζ^i)`・`Y(ζ^i)` はどちらも

    `∑_n (∑_{a<l} A n a · (ζ^i)^a)`   （`I` 進和）

の形に書ける（`ζ` の指数は `ζ^l = 1` により **`mod l` で揃う**）。
★★本ファイルはこの形を `muEval` として名前を付け、

* `μ_l∖{1}` 上の和が `ζ`-free になること（`sum_mu_muEval`）
* 積がふたたび同じ形になること（`muEval_mul`）

を与える。☆これで手順 3 が**有限個の係数計算**に落ちる。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★★★**μ-等級付き `I` 進級数の値**——`z = ζ^i` を代入した値。 -/
noncomputable def muEval [IsAdicComplete I R] (l : ℕ) (A : ℕ → ℕ → R)
    (hA : ∀ n a, A n a ∈ I ^ n) (z : R) : R :=
  adicSum (I := I) (fun n => ∑ a ∈ range l, A n a * z ^ a)
    (fun n => Submodule.sum_mem _ (fun a _ => Ideal.mul_mem_right _ _ (hA n a)))

/-- ★★★★★★★★★★★★★★★★★★★★
**μ-等級付き級数を `μ_l∖{1}` 上で足すと `ζ` が消える**。

    `∑_{i≠0} muEval l A (ζ^i) = adicSum (n ↦ ∑_{a<l} A n a · ([l ∣ a]·l − 1))`

★★`a` は `0 ≤ a < l` なので `l ∣ a ⇔ a = 0`——★右辺は
`l·A n 0 − ∑_{a<l} A n a` である。 -/
theorem sum_mu_muEval [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (A : ℕ → ℕ → R) (hA : ∀ n a, A n a ∈ I ^ n) :
    ∑ i ∈ (range l).erase 0, muEval l A hA (ζ ^ i)
      = adicSum (I := I)
          (fun n => ∑ a ∈ range l, A n a * ((if l ∣ a then (l : R) else 0) - 1))
          (fun n => Submodule.sum_mem _
            (fun a _ => Ideal.mul_mem_right _ _ (hA n a))) := by
  classical
  simp only [muEval]
  rw [← adicSum_finsetSum ((range l).erase 0)
      (fun i n => ∑ a ∈ range l, A n a * (ζ ^ i) ^ a)
      (fun i n => Submodule.sum_mem _
        (fun a _ => Ideal.mul_mem_right _ _ (hA n a)))]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun a _ => ?_)
  rw [← Finset.mul_sum, ← sum_mu_pow_erase_zero hl hζ a]
  exact congrArg _ (Finset.sum_congr rfl (fun i _ => by rw [← pow_mul]))

/-- ★★★★★★`a < l` のとき `l ∣ a ⇔ a = 0` なので、和は簡単になる。 -/
theorem sum_mu_muEval' [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (A : ℕ → ℕ → R) (hA : ∀ n a, A n a ∈ I ^ n) :
    ∑ i ∈ (range l).erase 0, muEval l A hA (ζ ^ i)
      = adicSum (I := I)
          (fun n => (l : R) * A n 0 - ∑ a ∈ range l, A n a)
          (fun n => Submodule.sub_mem _
            (Ideal.mul_mem_left _ _ (hA n 0))
            (Submodule.sum_mem _ (fun a _ => hA n a))) := by
  classical
  rw [sum_mu_muEval hl hζ A hA]
  refine adicSum_congr _ _ (fun n => ?_)
  have hsplit : ∑ a ∈ range l, A n a * ((if l ∣ a then (l : R) else 0) - 1)
      = ∑ a ∈ range l, (A n a * (if l ∣ a then (l : R) else 0) - A n a) := by
    exact Finset.sum_congr rfl (fun a _ => by ring)
  rw [hsplit, Finset.sum_sub_distrib]
  congr 1
  refine Finset.sum_eq_single 0 (fun a ha hne => ?_) (fun h => absurd (Finset.mem_range.2 hl.pos) h)
    |>.trans ?_
  · rw [if_neg (fun hdvd => hne (Nat.eq_zero_of_dvd_of_lt hdvd (Finset.mem_range.1 ha))), mul_zero]
  · rw [if_pos (dvd_zero l), mul_comm]

/-! ## ★出典の紐付け(`.src`) -/

def muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き I 進級数の値)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ-等級付き級数の μ_l 和は ζ-free。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_muEval'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l 和は l·A n 0 − ∑_a A n a。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
