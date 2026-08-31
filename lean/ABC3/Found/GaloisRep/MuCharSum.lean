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

/-! ## ★★★★★★★★`1/(1−ζ)` を `ζ` の多項式で書く -/

/-- ★★★★★★★★★★★★**`(1−η)·∑_{k<l} kη^k = −l`**。

`η^l = 1` と `∑_{k<l} η^k = 0` だけから出る。

☆機構: `T = ∑_{k<l} kη^k` と置くと
`∑_{k<l}(k+1)η^{k+1} = ηT + η·∑η^k = ηT` であり、
同じ和を左へずらすと `T + lη^l = T + l` になる。 -/
theorem one_sub_mul_sum_nsmul {l : ℕ} {η : R} (hpow : η ^ l = 1)
    (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (1 - η) * ∑ k ∈ range l, (k : R) * η ^ k = -(l : R) := by
  set T : R := ∑ k ∈ range l, (k : R) * η ^ k with hT
  have hshift : ∑ k ∈ range l, ((k : R) + 1) * η ^ (k + 1) = η * T := by
    have hpt : ∀ k : ℕ, ((k : R) + 1) * η ^ (k + 1)
        = η * ((k : R) * η ^ k) + η * η ^ k := by
      intro k; rw [pow_succ]; ring
    rw [Finset.sum_congr rfl (fun k _ => hpt k), Finset.sum_add_distrib,
      ← Finset.mul_sum, ← Finset.mul_sum, hsum, mul_zero, add_zero, hT]
  have hsplit : ∑ k ∈ range l, ((k : R) + 1) * η ^ (k + 1) = T + (l : R) := by
    have h1 := Finset.sum_range_succ' (fun j => (j : R) * η ^ j) l
    have h2 := Finset.sum_range_succ (fun j => (j : R) * η ^ j) l
    rw [h2, hpow, mul_one] at h1
    push_cast at h1
    simpa using h1.symm
  have hkey : η * T = T + (l : R) := by rw [← hshift, hsplit]
  rw [sub_mul, one_mul, hkey]
  ring

/-- ★★★★★★★★★★**`1/(1−η) = −(1/l)·∑_{k<l} kη^k`**。 -/
theorem inv_one_sub_eq {F : Type*} [Field F] {l : ℕ} (hl : (l : F) ≠ 0)
    {η : F} (hpow : η ^ l = 1) (hsum : ∑ k ∈ range l, η ^ k = 0) :
    (1 - η)⁻¹ = -(l : F)⁻¹ * ∑ k ∈ range l, (k : F) * η ^ k := by
  have hcore : (1 - η) * ∑ k ∈ range l, (k : F) * η ^ k = -(l : F) := by
    set T : F := ∑ k ∈ range l, (k : F) * η ^ k with hT
    have hshift : ∑ k ∈ range l, ((k : F) + 1) * η ^ (k + 1) = η * T := by
      have hpt : ∀ k : ℕ, ((k : F) + 1) * η ^ (k + 1)
          = η * ((k : F) * η ^ k) + η * η ^ k := by
        intro k; rw [pow_succ]; ring
      rw [Finset.sum_congr rfl (fun k _ => hpt k), Finset.sum_add_distrib,
        ← Finset.mul_sum, ← Finset.mul_sum, hsum, mul_zero, add_zero, hT]
    have hsplit : ∑ k ∈ range l, ((k : F) + 1) * η ^ (k + 1) = T + (l : F) := by
      have h1 := Finset.sum_range_succ' (fun j => (j : F) * η ^ j) l
      have h2 := Finset.sum_range_succ (fun j => (j : F) * η ^ j) l
      rw [h2, hpow, mul_one] at h1
      push_cast at h1
      simpa using h1.symm
    have hkey : η * T = T + (l : F) := by rw [← hshift, hsplit]
    rw [sub_mul, one_mul, hkey]
    ring
  refine inv_eq_of_mul_eq_one_right ?_
  rw [show (1 - η) * (-(l : F)⁻¹ * ∑ k ∈ range l, (k : F) * η ^ k)
      = -(l : F)⁻¹ * ((1 - η) * ∑ k ∈ range l, (k : F) * η ^ k) from by ring, hcore]
  field_simp

/-! ## ★★★★★★★★★★★★定数項 `∑_{ζ≠1} ζ/(1−ζ)² = −(l²−1)/12` -/

variable {F : Type*} [Field F]

/-- ★非自明な `l` 乗根の周りの事実をまとめて取る。 -/
theorem zeta_pow_facts {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l)
    {i : ℕ} (hi : i ∈ (range l).erase 0) :
    (ζ ^ i) ^ l = 1 ∧ ∑ k ∈ range l, (ζ ^ i) ^ k = 0 := by
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  have hnd : ¬ l ∣ i := fun h => hi0 (Nat.eq_zero_of_dvd_of_lt h hil)
  refine ⟨?_, ?_⟩
  · rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
  · exact (isPrimitiveRoot_pow_of_not_dvd (R := F) hl hζ hnd).geom_sum_eq_zero hl.one_lt

/-- ★★★★★★★★★★★★★★**各項を二重和に直す**。 -/
theorem term_eq_double_sum {l : ℕ} (hl : l.Prime) (hlF : (l : F) ≠ 0) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) {i : ℕ} (hi : i ∈ (range l).erase 0) :
    ζ ^ i * ((1 - ζ ^ i)⁻¹) ^ 2
      = ((l : F)⁻¹) ^ 2 * ∑ k ∈ range l, ∑ m ∈ range l,
          (k : F) * (m : F) * ζ ^ (i * (k + m + 1)) := by
  obtain ⟨hpow, hsum⟩ := zeta_pow_facts hl hζ hi
  have hRHS : ∑ k ∈ range l, ∑ m ∈ range l, (k : F) * (m : F) * ζ ^ (i * (k + m + 1))
      = ζ ^ i * ((∑ k ∈ range l, (k : F) * (ζ ^ i) ^ k)
          * ∑ m ∈ range l, (m : F) * (ζ ^ i) ^ m) := by
    rw [Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    have hexp : ζ ^ (i * (k + m + 1)) = ζ ^ i * ((ζ ^ i) ^ k * (ζ ^ i) ^ m) := by
      rw [← pow_mul, ← pow_mul, ← pow_add, ← pow_add]
      congr 1
      ring
    rw [hexp]
    ring
  rw [inv_one_sub_eq hlF hpow hsum, hRHS]
  ring

/-! ## ★★★★★★初等な和の公式 -/

theorem sum_range_cast [CharZero F] (n : ℕ) :
    ∑ k ∈ range n, (k : F) = (n : F) * ((n : F) - 1) / 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      push_cast
      field_simp
      ring

theorem sum_range_cast_sq [CharZero F] (n : ℕ) :
    ∑ k ∈ range n, (k : F) ^ 2 = (n : F) * ((n : F) - 1) * (2 * (n : F) - 1) / 6 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih]
      push_cast
      field_simp
      ring

/-! ## ★★★★★★★★★★★★★★★★定数項の値 -/

/-- ★★★★★★**内側の和**——`l ∣ k+m+1` なる `m` は `l−1−k` だけ。 -/
theorem sum_indicator_inner {l : ℕ} (hl : 0 < l) {k : ℕ} (hk : k ∈ range l) :
    ∑ m ∈ range l, (m : F) * (if l ∣ k + m + 1 then (1 : F) else 0)
      = ((l : F) - 1 - (k : F)) := by
  have hkl : k < l := Finset.mem_range.1 hk
  have hmem : l - 1 - k ∈ range l := Finset.mem_range.2 (by omega)
  have hunique : ∀ m ∈ range l, m ≠ l - 1 - k →
      (m : F) * (if l ∣ k + m + 1 then (1 : F) else 0) = 0 := by
    intro m hm hne
    have hml : m < l := Finset.mem_range.1 hm
    have hnd : ¬ (l ∣ k + m + 1) := by
      intro hdvd
      have heq : k + m + 1 = l :=
        Nat.eq_of_dvd_of_lt_two_mul (by omega) hdvd (by omega)
      omega
    rw [if_neg hnd, mul_zero]
  rw [Finset.sum_eq_single (l - 1 - k) hunique (fun h => absurd hmem h)]
  have hdvd : l ∣ k + (l - 1 - k) + 1 := ⟨1, by omega⟩
  rw [if_pos hdvd, mul_one]
  have h2 : k ≤ l - 1 := by omega
  have h1 : (1 : ℕ) ≤ l := hl
  rw [Nat.cast_sub h2, Nat.cast_sub h1]
  push_cast
  ring

/-- ★★★★★★★★★★★★★★★★★★★★
**`∑_{ζ ∈ μ_l, ζ≠1} ζ/(1−ζ)² = −(l²−1)/12`**。

★★Tate 曲線の Vélu の和の**定数項**である。

☆証明は `1/(1−ζ) = −(1/l)∑ kζ^k`（`inv_one_sub_eq`）で多項式に直し、
`μ_l` 上の指標和（`sum_mu_pow_erase_zero`）を使うだけ。
**微分も対称式も要らない。** -/
theorem sum_mu_frac [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, ζ ^ i * ((1 - ζ ^ i)⁻¹) ^ 2
      = -((l : F) ^ 2 - 1) / 12 := by
  have hlF : (l : F) ≠ 0 := Nat.cast_ne_zero.2 hl.ne_zero
  rw [Finset.sum_congr rfl (fun i hi => term_eq_double_sum hl hlF hζ hi), ← Finset.mul_sum]
  have hswap : ∑ i ∈ (range l).erase 0, ∑ k ∈ range l, ∑ m ∈ range l,
        (k : F) * (m : F) * ζ ^ (i * (k + m + 1))
      = ∑ k ∈ range l, ∑ m ∈ range l,
        (k : F) * (m : F) * ((if l ∣ k + m + 1 then (l : F) else 0) - 1) := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun k _ => ?_)
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl (fun m _ => ?_)
    rw [← Finset.mul_sum, sum_mu_pow_erase_zero hl hζ (k + m + 1)]
  rw [hswap]
  have hval : ∀ k ∈ range l, ∑ m ∈ range l,
        (k : F) * (m : F) * ((if l ∣ k + m + 1 then (l : F) else 0) - 1)
      = (l : F) * ((k : F) * ((l : F) - 1 - (k : F)))
        - (k : F) * ((l : F) * ((l : F) - 1) / 2) := by
    intro k hk
    have h1 : ∀ m : ℕ, (k : F) * (m : F) * ((if l ∣ k + m + 1 then (l : F) else 0) - 1)
        = (l : F) * ((k : F) * ((m : F) * (if l ∣ k + m + 1 then (1 : F) else 0)))
          - (k : F) * (m : F) := by
      intro m
      by_cases h : l ∣ k + m + 1 <;> simp [h] <;> ring
    rw [Finset.sum_congr rfl (fun m _ => h1 m), Finset.sum_sub_distrib,
      ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum,
      sum_indicator_inner hl.pos hk, sum_range_cast]
  rw [Finset.sum_congr rfl hval]
  have hexp : ∑ k ∈ range l, ((l : F) * ((k : F) * ((l : F) - 1 - (k : F)))
        - (k : F) * ((l : F) * ((l : F) - 1) / 2))
      = ((l : F) * ((l : F) - 1) - (l : F) * ((l : F) - 1) / 2) * (∑ k ∈ range l, (k : F))
        - (l : F) * ∑ k ∈ range l, (k : F) ^ 2 := by
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun k _ => by ring)
  rw [hexp, sum_range_cast, sum_range_cast_sq]
  field_simp
  ring

/-! ## ★★★★★★★★★★★★多項式を `μ_l∖{1}` 上で足す -/

/-- ★★★★★★★★★★★★★★★★**一般形**——多項式 `∑_j c_j X^j` を
`μ_l∖{1}` 上で足すと、係数に `l·[l ∣ j] − 1` を掛けた和になる。

★これが `tateModel_of_quot_mu` の手順 2 を**任意の係数列に対して**書いた形である。 -/
theorem sum_mu_poly {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l)
    (c : ℕ → F) (N : ℕ) :
    ∑ i ∈ (range l).erase 0, ∑ j ∈ range N, c j * ζ ^ (i * j)
      = ∑ j ∈ range N, c j * ((if l ∣ j then (l : F) else 0) - 1) := by
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [← Finset.mul_sum, sum_mu_pow_erase_zero hl hζ j]

/-- ★`k ∈ range l` のとき `l ∣ k ⇔ k = 0`。 -/
theorem dvd_iff_eq_zero_of_mem_range {l k : ℕ} (hl : 0 < l) (hk : k ∈ range l) :
    l ∣ k ↔ k = 0 := by
  constructor
  · intro h
    exact Nat.eq_zero_of_dvd_of_lt h (Finset.mem_range.1 hk)
  · rintro rfl
    exact dvd_zero l

/-- ★★★★★★★★★★★★**`∑_{ζ≠1} 1/(1−ζ) = (l−1)/2`**。 -/
theorem sum_mu_inv_one_sub [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, (1 - ζ ^ i)⁻¹ = ((l : F) - 1) / 2 := by
  have hlF : (l : F) ≠ 0 := Nat.cast_ne_zero.2 hl.ne_zero
  have hterm : ∀ i ∈ (range l).erase 0,
      (1 - ζ ^ i)⁻¹ = -(l : F)⁻¹ * ∑ k ∈ range l, (k : F) * ζ ^ (i * k) := by
    intro i hi
    obtain ⟨hpow, hsum⟩ := zeta_pow_facts hl hζ hi
    rw [inv_one_sub_eq hlF hpow hsum]
    congr 1
    exact Finset.sum_congr rfl (fun k _ => by rw [← pow_mul])
  rw [Finset.sum_congr rfl hterm, ← Finset.mul_sum,
    sum_mu_poly hl hζ (fun k => (k : F)) l]
  have hval : ∑ k ∈ range l, (k : F) * ((if l ∣ k then (l : F) else 0) - 1)
      = -((l : F) * ((l : F) - 1) / 2) := by
    have h1 : ∀ k ∈ range l, (k : F) * ((if l ∣ k then (l : F) else 0) - 1) = -(k : F) := by
      intro k hk
      by_cases hk0 : k = 0
      · subst hk0; simp
      · rw [if_neg (fun h => hk0 ((dvd_iff_eq_zero_of_mem_range hl.pos hk).1 h))]
        ring
    rw [Finset.sum_congr rfl h1, Finset.sum_neg_distrib, sum_range_cast]
  rw [hval]
  field_simp

/-- ★★★★★★★★★★★★**`∑_{ζ≠1} 1/(1−ζ)² = (l−1)(5−l)/12`**。

`ζ/(1−ζ)² = 1/(1−ζ)² − 1/(1−ζ)` なので `sum_mu_frac` と
`sum_mu_inv_one_sub` から出る。 -/
theorem sum_mu_inv_one_sub_sq [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, ((1 - ζ ^ i)⁻¹) ^ 2 = ((l : F) - 1) * (5 - (l : F)) / 12 := by
  have hsplit : ∀ i ∈ (range l).erase 0,
      ((1 - ζ ^ i)⁻¹) ^ 2 = ζ ^ i * ((1 - ζ ^ i)⁻¹) ^ 2 + (1 - ζ ^ i)⁻¹ := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hnd : ¬ l ∣ i := fun h => hi0 (Nat.eq_zero_of_dvd_of_lt h hil)
    have hne : (1 : F) - ζ ^ i ≠ 0 := by
      refine sub_ne_zero.2 (Ne.symm ?_)
      intro h
      exact (isPrimitiveRoot_pow_of_not_dvd (R := F) hl hζ hnd).ne_one hl.one_lt h
    field_simp
    ring
  rw [Finset.sum_congr rfl hsplit, Finset.sum_add_distrib, sum_mu_frac hl hζ,
    sum_mu_inv_one_sub hl hζ]
  ring

/-! ## ★★★★★★★★★★★★★★逆元側の指標和と Tate 係数 -/

/-- ★`ζ^i` の逆元は `ζ^(l−i)`。 -/
theorem inv_zeta_pow {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) {i : ℕ} (hi : i ≤ l) :
    (ζ ^ i)⁻¹ = ζ ^ (l - i) := by
  have hne : ζ ^ i ≠ 0 := by
    intro h
    have : (1 : F) = 0 := by
      rw [← hζ.pow_eq_one, ← Nat.sub_add_cancel hi, pow_add, h, mul_zero]
    exact one_ne_zero this
  refine inv_eq_of_mul_eq_one_right ?_
  rw [← pow_add, Nat.add_sub_cancel' hi, hζ.pow_eq_one]

/-- ★★★★★★★★★★★★**逆元側の指標和**——`∑_{ζ≠1} ζ^{−d}` も同じ値。 -/
theorem sum_mu_inv_pow {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (d : ℕ) :
    ∑ i ∈ (range l).erase 0, ((ζ ^ i)⁻¹) ^ d = (if l ∣ d then (l : F) else 0) - 1 := by
  have hrw : ∀ i ∈ (range l).erase 0, ((ζ ^ i)⁻¹) ^ d = ζ ^ ((l - i) * d) := by
    intro i hi
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    rw [inv_zeta_pow hζ hil.le, ← pow_mul]
  rw [Finset.sum_congr rfl hrw]
  rw [← sum_mu_pow_erase_zero (R := F) hl hζ d]
  refine Finset.sum_nbij' (fun i => l - i) (fun i => l - i) ?_ ?_ ?_ ?_ ?_
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    exact Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    exact Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    omega
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    omega
  · intro a _
    rfl

/-- ★★★★★★★★★★★★★★★★★★★★
**Tate の `X` の `q^N` 係数を `μ_l∖{1}` 上で足す**。

古典的な `q` 展開 `X(u) = u/(1−u)² + ∑_{N≥1} c_N(u) q^N`、
`c_N(u) = ∑_{d ∣ N} d(u^d + u^{−d} − 2)` の係数を、
`μ_l∖{1}` 上で足すと

    `∑_ζ c_N(ζ) = 2l·∑_{d ∣ N} d([l ∣ d] − 1)`

になる。★★右辺は `ζ` を含まない——これが
`σ₁(q) → σ₁(q^l)` を生む機構である。 -/
theorem sum_mu_coeff {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l)
    (D : Finset ℕ) :
    ∑ i ∈ (range l).erase 0, ∑ d ∈ D, (d : F) * ((ζ ^ i) ^ d + ((ζ ^ i)⁻¹) ^ d - 2)
      = 2 * (l : F) * ∑ d ∈ D, (d : F) * ((if l ∣ d then (1 : F) else 0) - 1) := by
  rw [Finset.sum_comm]
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun d _ => ?_)
  have hexp : ∀ i, (d : F) * ((ζ ^ i) ^ d + ((ζ ^ i)⁻¹) ^ d - 2)
      = (d : F) * (ζ ^ (i * d)) + (d : F) * (((ζ ^ i)⁻¹) ^ d) - (d : F) * 2 := by
    intro i; rw [← pow_mul]; ring
  rw [Finset.sum_congr rfl (fun i _ => hexp i), Finset.sum_sub_distrib,
    Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum, ← Finset.sum_mul,
    sum_mu_pow_erase_zero hl hζ d, sum_mu_inv_pow hl hζ d]
  have hcard : ((range l).erase 0).card = l - 1 := by
    rw [Finset.card_erase_of_mem (Finset.mem_range.2 hl.pos), Finset.card_range]
  rw [Finset.sum_const, hcard, nsmul_eq_mul]
  have hlc : ((l - 1 : ℕ) : F) = (l : F) - 1 := by
    have h1 : (1 : ℕ) ≤ l := hl.pos
    rw [Nat.cast_sub h1, Nat.cast_one]
  rw [hlc]
  by_cases h : l ∣ d <;> simp [h] <;> ring

/-! ## ★出典の紐付け(`.src`) -/

def one_sub_mul_sum_nsmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)((1−η)·∑ kη^k = −l。★無条件)",
    sectionId := "genell-lemma-3-2" }

def inv_one_sub_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(1/(1−η) を η の多項式で書く。★無条件)",
    sectionId := "genell-lemma-3-2" }

def zeta_pow_facts.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(非自明な l 乗根の事実。★無条件)",
    sectionId := "genell-lemma-3-2" }

def term_eq_double_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ/(1−ζ)² を二重和に直す。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_frac.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{ζ≠1} ζ/(1−ζ)² = −(l²−1)/12。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_range_cast.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{k<n} k = n(n−1)/2。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_range_cast_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{k<n} k² = n(n−1)(2n−1)/6。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_indicator_inner.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(l ∣ k+m+1 なる m は l−1−k だけ。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_poly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(多項式を μ_l∖{1} 上で足す一般形。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_inv_one_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{ζ≠1} 1/(1−ζ) = (l−1)/2。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_inv_one_sub_sq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{ζ≠1} 1/(1−ζ)² = (l−1)(5−l)/12。★無条件)",
    sectionId := "genell-lemma-3-2" }

def dvd_iff_eq_zero_of_mem_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(range l の中で l ∣ k ⇔ k = 0。★無条件)",
    sectionId := "genell-lemma-3-2" }

def inv_zeta_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(ζ^i の逆元は ζ^(l−i)。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_inv_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(逆元側の指標和。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_coeff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate の X の q^N 係数を μ_l∖{1} 上で足す。★無条件)",
    sectionId := "genell-lemma-3-2" }

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
