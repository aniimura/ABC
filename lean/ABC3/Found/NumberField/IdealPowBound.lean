/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.DedekindEuler
import Mathlib.RingTheory.UniqueFactorizationDomain.Finsupp

/-!
# 素冪でのイデアル計数の多項式上界(鎖 `cheb` の `cheb-log-zeta` の前段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何のためか

`log ζ_L(s) = Σ_{f(𝔭)=1} N𝔭^{-s} + O(1)`(`cheb-log-zeta`)の **`O(1)` を
`s → 1+` で一様に押さえる**には、局所因子

  `L_p(s) = Σ_e a(p^e) p^{-es}`

の `e ≥ 2` の部分が `Σ_p` について一様に有界でなければならない。
★これには `a(p^e)` の **多項式増大**が要る。

★★★**実測(2026-08-25)**: 部分和の漸近から出る `a(n) ≤ C·n` では**弱すぎる**
(`Σ_p Σ_{e≥2} C p^e p^{-es}` は `s → 1+` で発散する)。
本ファイルで `a(p^e) ≤ (e·d + 1)^d`(`d = [K:ℚ]`)を証明する。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `factorization_le_of_dvd` | 約数の `factorization` は元の `factorization` 以下 |
| `card_dvd_le` | ★★**約数の個数の上界**(`factorization` が単射だから) |
| `le_absNorm_of_mem_normalizedFactors` | `(p)` の素因子のノルムは `p` 以上 |
| `card_normalizedFactors_span_le` | ★`(p)` の素因子の個数(重複込み)は `d` 以下 |
| `idealCount_prime_pow_le` | ★★★★**`a(p^e) ≤ (e·d + 1)^d`** |

★★要は「`(p)` のノルムは `p^d`、各素因子のノルムは `p` 以上、
ゆえに素因子は重複込みで高々 `d` 個」である。
-/

namespace ABC3.Found.NF

open UniqueFactorizationMonoid NumberField Ideal

/-! ## ★1. 約数の個数の上界(一般の一意分解モノイド) -/

open scoped Classical in
/-- ★約数の `factorization` は元の `factorization` 以下。 -/
theorem factorization_le_of_dvd {α : Type*} [CommMonoidWithZero α]
    [UniqueFactorizationMonoid α] [NormalizationMonoid α] {I J : α}
    (hJ : J ≠ 0) (h : I ∣ J) (q : α) : factorization I q ≤ factorization J q := by
  obtain ⟨c, rfl⟩ := h
  have hI : I ≠ 0 := fun h => hJ (by simp [h])
  have hc : c ≠ 0 := fun h => hJ (by simp [h])
  rw [factorization_mul hI hc]
  simp

open scoped Classical in
/-- ★★**約数の個数の上界** —— `factorization` が単射だから。

★`P` を満たすものがすべて `J` の約数なら、`factorization` の値の組
`(support J) → Fin (N+1)` への写像が単射になる。 -/
theorem card_dvd_le {α : Type*} [CommMonoidWithZero α] [UniqueFactorizationMonoid α]
    [NormalizationMonoid α] [Subsingleton αˣ] {J : α} (hJ : J ≠ 0) (N : ℕ)
    (hN : ∀ q, factorization J q ≤ N) (P : α → Prop) (hP : ∀ I, P I → I ∣ J)
    (hP0 : ∀ I, P I → I ≠ 0) :
    Nat.card {I : α // P I} ≤ (N + 1) ^ (factorization J).support.card := by
  classical
  have hinj : Function.Injective
      (fun I : {I : α // P I} => fun q : ((factorization J).support : Finset α) =>
        (⟨factorization I.1 q, Nat.lt_succ_of_le (le_trans
          (factorization_le_of_dvd hJ (hP I.1 I.2) q) (hN q))⟩ : Fin (N + 1))) := by
    intro I I' h
    refine Subtype.ext (associated_iff_eq.mp (associated_of_factorization_eq _ _
      (hP0 I.1 I.2) (hP0 I'.1 I'.2) (Finsupp.ext fun q => ?_)))
    by_cases hq : q ∈ (factorization J).support
    · exact congrArg Fin.val (congrFun h ⟨q, hq⟩)
    · have h0 : factorization J q = 0 := Finsupp.notMem_support_iff.mp hq
      have h1 := factorization_le_of_dvd hJ (hP I.1 I.2) q
      have h2 := factorization_le_of_dvd hJ (hP I'.1 I'.2) q
      omega
  calc Nat.card {I : α // P I}
      ≤ Nat.card (((factorization J).support : Finset α) → Fin (N + 1)) :=
        Nat.card_le_card_of_injective _ hinj
    _ = (N + 1) ^ (factorization J).support.card := by
        rw [Nat.card_eq_fintype_card, Fintype.card_fun, Fintype.card_fin, Fintype.card_coe]

open scoped Classical in
theorem support_nsmul_eq_nat {α : Type*} {e : ℕ} (he : e ≠ 0) (f : α →₀ ℕ) :
    (e • f).support = f.support := by
  ext q
  simp [Finsupp.mem_support_iff, he]

/-! ## ★2. `(p)` の素因子 -/

variable (K : Type*) [Field K] [NumberField K]

theorem natCast_ne_zero_ringOfIntegers {p : ℕ} (hp : p ≠ 0) : (p : 𝓞 K) ≠ 0 :=
  Nat.cast_ne_zero.mpr hp

theorem span_natCast_ne_zero {p : ℕ} (hp : p ≠ 0) : Ideal.span {(p : 𝓞 K)} ≠ 0 := by
  rw [Ideal.zero_eq_bot, Ne, Ideal.span_singleton_eq_bot]
  exact natCast_ne_zero_ringOfIntegers K hp

/-- ★`(p)` の素因子のノルムは `p` 以上 —— ノルムは `p^d` を割る `p` 冪で、`1` ではない。 -/
theorem le_absNorm_of_mem_normalizedFactors {p : ℕ} (hp : p.Prime)
    {q : Ideal (𝓞 K)} (hq : q ∈ normalizedFactors (Ideal.span {(p : 𝓞 K)})) :
    p ≤ absNorm q := by
  obtain ⟨c, hc⟩ := dvd_of_mem_normalizedFactors hq
  have hdvd : absNorm q ∣ p ^ Module.finrank ℤ (𝓞 K) := by
    rw [← absNorm_span_natCast (R := 𝓞 K) p, hc, map_mul]
    exact Dvd.intro _ rfl
  obtain ⟨j, hj, hqj⟩ := (Nat.dvd_prime_pow hp).mp hdvd
  have hne : absNorm q ≠ 1 := by
    intro h
    exact (irreducible_of_normalized_factor q hq).not_isUnit
      (Ideal.isUnit_iff.mpr (absNorm_eq_one_iff.mp h))
  rcases Nat.eq_zero_or_pos j with rfl | hjpos
  · rw [pow_zero] at hqj; exact absurd hqj hne
  · rw [hqj]
    calc p = p ^ 1 := (pow_one p).symm
      _ ≤ p ^ j := Nat.pow_le_pow_right hp.pos hjpos

/-- ★★`(p)` の素因子の個数(重複込み)は `d` 以下 —— `p^{個数} ≤ N((p)) = p^d`。 -/
theorem card_normalizedFactors_span_le {p : ℕ} (hp : p.Prime) :
    Multiset.card (normalizedFactors (Ideal.span {(p : 𝓞 K)}))
      ≤ Module.finrank ℤ (𝓞 K) := by
  set J : Ideal (𝓞 K) := Ideal.span {(p : 𝓞 K)} with hJdef
  have hJ0 : J ≠ 0 := span_natCast_ne_zero K hp.ne_zero
  have hprod : (normalizedFactors J).prod = J :=
    associated_iff_eq.mp (prod_normalizedFactors hJ0)
  have h1 : p ^ Multiset.card (normalizedFactors J)
      ≤ ((normalizedFactors J).map absNorm).prod := by
    refine le_trans (le_of_eq (by rw [Multiset.card_map])) (Multiset.pow_card_le_prod ?_)
    intro x hx
    obtain ⟨q, hq, rfl⟩ := Multiset.mem_map.mp hx
    exact le_absNorm_of_mem_normalizedFactors K hp hq
  have h2 : ((normalizedFactors J).map absNorm).prod = p ^ Module.finrank ℤ (𝓞 K) := by
    rw [← map_multiset_prod, hprod, hJdef, absNorm_span_natCast]
  rw [h2] at h1
  exact (Nat.pow_le_pow_iff_right hp.one_lt).mp h1

/-! ## ★3. 素冪でのイデアル計数の多項式上界 -/

open scoped Classical in
/-- ★★★★**`a(p^e) ≤ (e·d + 1)^d`**(`d = [K:ℚ]`)。

★★これが `cheb-log-zeta` の `e ≥ 2` の項を `s → 1+` で一様に押さえる鍵である。 -/
theorem idealCount_prime_pow_le {p : ℕ} (hp : p.Prime) (e : ℕ) :
    idealCount K (p ^ e)
      ≤ (e * Module.finrank ℤ (𝓞 K) + 1) ^ Module.finrank ℤ (𝓞 K) := by
  set d := Module.finrank ℤ (𝓞 K) with hd
  rcases Nat.eq_zero_or_pos e with rfl | he
  · simp [idealCount_one]
  set Jp : Ideal (𝓞 K) := Ideal.span {(p : 𝓞 K)} with hJp
  set J : Ideal (𝓞 K) := Jp ^ e with hJ
  have hJp0 : Jp ≠ 0 := span_natCast_ne_zero K hp.ne_zero
  have hJ0 : J ≠ 0 := pow_ne_zero _ hJp0
  have hfac : _root_.factorization J = e • _root_.factorization Jp := factorization_pow
  have hcard := card_normalizedFactors_span_le K hp
  have hN : ∀ q, _root_.factorization J q ≤ e * d := by
    intro q
    rw [hfac]
    show e * _root_.factorization Jp q ≤ e * d
    refine Nat.mul_le_mul_left _ ?_
    rw [factorization_eq_count]
    exact le_trans (Multiset.count_le_card _ _) hcard
  have hsupp : (_root_.factorization J).support.card ≤ d := by
    rw [hfac, support_nsmul_eq_nat he.ne' _, support_factorization]
    exact le_trans (Multiset.toFinset_card_le _) hcard
  have hdvd : ∀ I : Ideal (𝓞 K), absNorm I = p ^ e → I ∣ J := by
    intro I hI
    rw [Ideal.dvd_iff_le, hJ, hJp, Ideal.span_singleton_pow]
    refine Ideal.span_le.mpr ?_
    rintro x hx
    rw [Set.mem_singleton_iff] at hx
    subst hx
    have h := Ideal.absNorm_mem I
    rw [hI] at h
    simpa using h
  have hne : ∀ I : Ideal (𝓞 K), absNorm I = p ^ e → I ≠ 0 := by
    intro I hI h
    rw [h] at hI
    simp only [Ideal.zero_eq_bot, absNorm_bot] at hI
    exact absurd hI.symm (pow_ne_zero e hp.ne_zero)
  calc idealCount K (p ^ e)
      ≤ (e * d + 1) ^ (_root_.factorization J).support.card :=
        card_dvd_le hJ0 (e * d) hN (fun I => absNorm I = p ^ e) hdvd hne
    _ ≤ (e * d + 1) ^ d := Nat.pow_le_pow_right (Nat.succ_pos _) hsupp

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `cheb-log-zeta` の一様評価の鍵。 -/
def idealCount_prime_pow_le.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — 素冪でのイデアル計数の多項式上界",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
