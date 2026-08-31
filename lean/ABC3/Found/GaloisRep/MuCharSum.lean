/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.RootsOfUnity.PrimitiveRoots
import ABC3.Meta.Claim

/-!
# `μ_l` 上の指標和（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`Skeleton/GenEll/TateIsogeny.lean` の `tateModel_of_quot_mu`（葉 1）の
**手順 2**——「`ζ` について足すと `ζ` の指数が `l` の倍数の項だけ残る」——である。

    `∑_{ζ ∈ μ_l} ζ^k = l`  （`l ∣ k`）
    `∑_{ζ ∈ μ_l} ζ^k = 0`  （`l ∤ k`）

★これが `σ_k(q) → σ_k(q^l)` を生む機構であり、Tate 曲線の `q` 展開恒等式
`a₄(q^l) = a₄(q) − 5v`・`a₆(q^l) = a₆(q) − v − 7w` の核である。

☆材料は mathlib の `IsPrimitiveRoot.geom_sum_eq_zero` だけである。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R : Type*} [CommRing R] [IsDomain R]

/-- ★★★★★★★★★★★★**`μ_l` 上の指標和**——`∑_{i<l} ζ^{ik} = l` か `0`。

`ζ` が原始 `l` 乗根で `l` が素数なら、`ζ^k` は `l ∣ k` のとき `1`、
そうでないとき再び原始 `l` 乗根である。 -/
theorem sum_mu_pow {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    ∑ i ∈ range l, ζ ^ (i * k) = if l ∣ k then (l : R) else 0 := by
  have hkey : ∑ i ∈ range l, ζ ^ (i * k) = ∑ i ∈ range l, (ζ ^ k) ^ i := by
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [← pow_mul, mul_comm k i]
  rw [hkey]
  by_cases hdvd : l ∣ k
  · obtain ⟨m, rfl⟩ := hdvd
    have hone : ζ ^ (l * m) = 1 := by
      rw [pow_mul, hζ.pow_eq_one, one_pow]
    rw [if_pos (dvd_mul_right l m), hone]
    simp
  · rw [if_neg hdvd]
    have hcop : Nat.Coprime k l := ((Nat.Prime.coprime_iff_not_dvd hl).2 hdvd).symm
    have hprim : IsPrimitiveRoot (ζ ^ k) l := hζ.pow_of_coprime k hcop
    exact hprim.geom_sum_eq_zero hl.one_lt

/-- ★★★★★★★★**非自明な `l` 乗根にわたる和**——`∑_{i≠0} ζ^{ik}`。

`veluVFull`・`veluWFull` は `H∖{O}` 全体にわたる和なので、こちらの形で使う。 -/
theorem sum_mu_pow_erase_zero {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    ∑ i ∈ (range l).erase 0, ζ ^ (i * k) = (if l ∣ k then (l : R) else 0) - 1 := by
  have hmem : (0 : ℕ) ∈ range l := mem_range.2 hl.pos
  have hsplit : ∑ i ∈ range l, ζ ^ (i * k)
      = ζ ^ (0 * k) + ∑ i ∈ (range l).erase 0, ζ ^ (i * k) :=
    (Finset.add_sum_erase _ _ hmem).symm
  rw [sum_mu_pow hl hζ k] at hsplit
  simp only [zero_mul, pow_zero] at hsplit
  rw [hsplit]
  ring

/-- ★★★★★★`ζ^i`（`i` が `l` を割らない）は再び原始 `l` 乗根。 -/
theorem isPrimitiveRoot_pow_of_not_dvd {l : ℕ} (hl : l.Prime) {ζ : R}
    (hζ : IsPrimitiveRoot ζ l) {i : ℕ} (hi : ¬ l ∣ i) : IsPrimitiveRoot (ζ ^ i) l :=
  hζ.pow_of_coprime i (((Nat.Prime.coprime_iff_not_dvd hl).2 hi).symm)

/-! ## ★出典の紐付け(`.src`) -/

def sum_mu_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(μ_l 上の指標和——∑ ζ^{ik} は l か 0。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_pow_erase_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(非自明な l 乗根にわたる指標和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def isPrimitiveRoot_pow_of_not_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ^i は再び原始 l 乗根。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
