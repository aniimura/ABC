/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.DedekindEuler
import ABC3.Found.NumberField.SplitInfinite
import Mathlib.NumberTheory.RamificationInertia.Basic
import Mathlib.NumberTheory.RamificationInertia.Galois

/-!
# イデアル計数 `a_L(p)` を分解の様子で読む(鎖 `cheb` の `cheb-split-density` の代数側)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何をするか

`PrimeDensity.lean` で解析側 `Σ_p a_L(p)p^{-s} / log(1/(s−1)) → 1` が済んだ。
残るのは `a_L(p) = #{I : N I = p}` を**分解の様子**で読むことである:

| 状況 | `a_L(p)` |
|---|---|
| `p` が完全分解 | `[L:ℚ]` |
| `p` が不分岐で完全分解しない(`L/ℚ` Galois) | `0` |
| つねに | `≤ [L:ℚ]` |

## ★★★要点 —— ノルムが `p` のイデアルは「`p` の上の剰余次数 1 の素イデアル」

* `N I = p`(素数)なら `I` は**既約**(ノルムが素数だから分解できない)、
  ゆえに素イデアル。
* `N I = p ∈ I` なので `I` は `(p)` の上にある。
* `N P = p^{f(P)}`(mathlib の `Ideal.absNorm_eq_pow_inertiaDeg`)なので
  `N P = p ⟺ f(P) = 1`。

★`L/ℚ` が Galois なら `f` は `p` の上のすべての素で等しいので、
`a_L(p)` は「全部」か「ゼロ」しかない。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `isPrime_of_absNorm_prime` | ノルムが素数なら素イデアル |
| `liesOver_of_absNorm_prime` | ノルムが `p` なら `(p)` の上にある |
| `absNorm_eq_iff_inertiaDeg_eq_one` | `N P = p ⟺ f(P) = 1` |
| `idealCount_eq_ncard` | ★★`a_L(p) = #{P ∣ p : f(P) = 1}` |
| `ncard_primesOver_eq_finrank_iff_int` | 完全分解 ⟺ すべての `P` で `e·f = 1`(底が `ℤ`) |
| `idealCount_eq_finrank_of_splitsCompletely` | ★★★完全分解なら `a_L(p) = [L:ℚ]` |
| `idealCount_eq_zero_of_unramified_of_not_split` | ★★★不分岐で非完全分解なら `a_L(p) = 0` |
| `idealCount_le_finrank` | つねに `a_L(p) ≤ [L:ℚ]` |
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal UniqueFactorizationMonoid
open scoped _root_.NumberField

variable {L : Type*} [Field L] [NumberField L]

/-! ## ★1. `(p)` は `ℤ` の極大イデアル -/

theorem span_intCast_isMaximal {p : ℕ} (hp : p.Prime) :
    (Ideal.span {(p : ℤ)}).IsMaximal :=
  ((Ideal.span_singleton_prime (by exact_mod_cast hp.ne_zero)).mpr
    (Nat.prime_iff_prime_int.mp hp)).isMaximal (by simpa using hp.ne_zero)

theorem span_intCast_ne_bot {p : ℕ} (hp : p.Prime) : Ideal.span {(p : ℤ)} ≠ ⊥ := by
  rw [Ne, Ideal.span_singleton_eq_bot]
  exact_mod_cast hp.ne_zero

/-! ## ★2. ノルムが素数のイデアル -/

/-- ★★**ノルムが素数のイデアルは素イデアル** —— ノルムが素数だから分解できない。 -/
theorem isPrime_of_absNorm_prime {p : ℕ} (hp : p.Prime) {I : Ideal (𝓞 L)}
    (h : absNorm I = p) : I.IsPrime := by
  have hirr : Irreducible I := by
    constructor
    · intro hu
      rw [Ideal.isUnit_iff] at hu
      rw [hu] at h
      simp at h
      exact hp.ne_one h.symm
    · rintro a b rfl
      rw [map_mul] at h
      rcases (Nat.Prime.eq_one_or_self_of_dvd hp (absNorm a) ⟨absNorm b, h.symm⟩) with h1 | h1
      · left; exact Ideal.isUnit_iff.mpr (absNorm_eq_one_iff.mp h1)
      · right
        rw [h1] at h
        have hb : absNorm b = 1 := by
          have hppos : 0 < p := hp.pos
          nlinarith [Nat.one_le_iff_ne_zero.mpr (show absNorm b ≠ 0 by
            intro hz; rw [hz, mul_zero] at h; omega)]
        exact Ideal.isUnit_iff.mpr (absNorm_eq_one_iff.mp hb)
  rw [← Ideal.prime_iff_isPrime]
  · exact UniqueFactorizationMonoid.irreducible_iff_prime.mp hirr
  · intro hz
    rw [hz] at h
    simp at h
    exact hp.ne_zero h.symm

/-- ★★**ノルムが `p` のイデアルは `(p)` の上にある** —— `N I = p ∈ I` だから。 -/
theorem liesOver_of_absNorm_prime {p : ℕ} (hp : p.Prime) {I : Ideal (𝓞 L)}
    (h : absNorm I = p) : I.LiesOver (Ideal.span {(p : ℤ)}) := by
  haveI := isPrime_of_absNorm_prime hp h
  constructor
  have hmem : ((p : ℤ)) ∈ Ideal.under ℤ I := by
    show algebraMap ℤ (𝓞 L) (p : ℤ) ∈ I
    have hI := Ideal.absNorm_mem I
    rw [h] at hI
    have hcast : algebraMap ℤ (𝓞 L) (p : ℤ) = ((p : ℕ) : 𝓞 L) := by push_cast; ring
    rw [hcast]
    exact hI
  have hle : Ideal.span {(p : ℤ)} ≤ Ideal.under ℤ I := Ideal.span_le.mpr (by
    rintro x hx; rw [Set.mem_singleton_iff] at hx; exact hx ▸ hmem)
  have hne : Ideal.under ℤ I ≠ ⊤ := by
    intro htop
    have h1 : (1 : ℤ) ∈ Ideal.under ℤ I := htop ▸ Submodule.mem_top
    have hone : (1 : 𝓞 L) ∈ I := by simpa using h1
    exact (isPrime_of_absNorm_prime hp h).ne_top (Ideal.eq_top_of_isUnit_mem I hone isUnit_one)
  exact (span_intCast_isMaximal hp).eq_of_le hne hle

/-! ## ★3. 基本等式(底が `ℤ` の形) -/

theorem one_le_ef_int {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L))
    (hP : P ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L)) :
    1 ≤ ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [← Finset.mem_coe,
    IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L)] at hP
  haveI : P.IsPrime := hP.1
  haveI : P.LiesOver (span {(p : ℤ)}) := hP.2
  have h2 : 0 < inertiaDeg (span {(p : ℤ)}) P := Ideal.inertiaDeg_pos' _ P
  have h1 : ramificationIdx (span {(p : ℤ)}) P ≠ 0 :=
    Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P (span_intCast_ne_bot hp)
  exact Nat.one_le_iff_ne_zero.mpr (Nat.mul_ne_zero h1 h2.ne')

theorem sum_ef_int {p : ℕ} (hp : p.Prime) :
    ∑ P ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L),
      ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P
      = Module.finrank ℚ L := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  exact Ideal.sum_ramification_inertia (𝓞 L) ℚ L (span_intCast_ne_bot hp)

theorem finite_primesOver_int {p : ℕ} (hp : p.Prime) :
    (primesOver (span {(p : ℤ)}) (𝓞 L)).Finite := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [← IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L)]
  exact (IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L)).finite_toSet

theorem ncard_primesOver_le_finrank_int {p : ℕ} (hp : p.Prime) :
    (primesOver (span {(p : ℤ)}) (𝓞 L)).ncard ≤ Module.finrank ℚ L := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [← IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L),
    Set.ncard_coe_finset, ← sum_ef_int (L := L) hp]
  have h := Finset.card_nsmul_le_sum
    (IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L))
    (fun P => ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P) 1
    (fun P hP => one_le_ef_int hp P hP)
  simpa using h

/-- ★★**完全分解 ⟺ すべての `P` で `e·f = 1`**(底が `ℤ` の形)。 -/
theorem ncard_primesOver_eq_finrank_iff_int {p : ℕ} (hp : p.Prime) :
    (primesOver (span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L ↔
      ∀ P ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L),
        ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P = 1 := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  have hsum := sum_ef_int (L := L) hp
  rw [← IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L),
    Set.ncard_coe_finset]
  constructor
  · intro hcard P hP
    by_contra hne
    have hge : 2 ≤ ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P := by
      have := one_le_ef_int hp P hP
      omega
    have hlt : (IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L)).card
        < ∑ Q ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L),
            ramificationIdx (span {(p : ℤ)}) Q * inertiaDeg (span {(p : ℤ)}) Q := by
      have hle : ∀ Q ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L),
          1 ≤ ramificationIdx (span {(p : ℤ)}) Q * inertiaDeg (span {(p : ℤ)}) Q :=
        fun Q hQ => one_le_ef_int hp Q hQ
      calc (IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L)).card
          = ∑ _Q ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L), 1 := by simp
        _ < _ := Finset.sum_lt_sum hle ⟨P, hP, by omega⟩
    omega
  · intro hall
    rw [← hsum, Finset.sum_congr rfl hall]
    simp

/-! ## ★4. `a_L(p) = #{P ∣ p : f(P) = 1}` -/

theorem absNorm_eq_pow_inertiaDeg' {p : ℕ} (hp : p.Prime) (P : Ideal (𝓞 L))
    [P.LiesOver (span {(p : ℤ)})] :
    absNorm P = p ^ inertiaDeg (span {(p : ℤ)}) P := by
  have h := Ideal.absNorm_eq_pow_inertiaDeg (R := 𝓞 L) P (Nat.prime_iff_prime_int.mp hp)
  rwa [Int.natAbs_natCast] at h

/-- ★★**`N P = p ⟺ f(P) = 1`**。 -/
theorem absNorm_eq_iff_inertiaDeg_eq_one {p : ℕ} (hp : p.Prime) {P : Ideal (𝓞 L)}
    [P.LiesOver (span {(p : ℤ)})] :
    absNorm P = p ↔ inertiaDeg (span {(p : ℤ)}) P = 1 := by
  rw [absNorm_eq_pow_inertiaDeg' hp P]
  constructor
  · intro h
    refine Nat.pow_right_injective hp.two_le (a₁ := inertiaDeg (span {(p : ℤ)}) P) (a₂ := 1) ?_
    show p ^ inertiaDeg (span {(p : ℤ)}) P = p ^ 1
    rw [pow_one]; exact h
  · intro h; rw [h, pow_one]

/-- ★★★**ノルムが `p` のイデアル = `p` の上の剰余次数 1 の素イデアル**。 -/
theorem setOf_absNorm_eq_prime {p : ℕ} (hp : p.Prime) :
    {I : Ideal (𝓞 L) | absNorm I = p}
      = {P | P ∈ primesOver (span {(p : ℤ)}) (𝓞 L) ∧ inertiaDeg (span {(p : ℤ)}) P = 1} := by
  ext I
  simp only [Set.mem_setOf_eq, Ideal.primesOver, Set.mem_setOf_eq]
  constructor
  · intro h
    haveI := isPrime_of_absNorm_prime hp h
    haveI := liesOver_of_absNorm_prime hp h
    exact ⟨⟨‹_›, ‹_›⟩, (absNorm_eq_iff_inertiaDeg_eq_one hp).mp h⟩
  · rintro ⟨⟨hP, hLO⟩, hf⟩
    haveI := hP
    haveI := hLO
    exact (absNorm_eq_iff_inertiaDeg_eq_one hp).mpr hf

/-- ★★★★**`a_L(p) = #{P ∣ p : f(P) = 1}`**。 -/
theorem idealCount_eq_ncard {p : ℕ} (hp : p.Prime) :
    idealCount L p
      = {P | P ∈ primesOver (span {(p : ℤ)}) (𝓞 L)
          ∧ inertiaDeg (span {(p : ℤ)}) P = 1}.ncard := by
  rw [idealCount, ← setOf_absNorm_eq_prime (L := L) hp, Set.ncard]
  rfl

/-! ## ★5. 分解の様子で読む -/

theorem idealCount_le_finrank {p : ℕ} (hp : p.Prime) :
    idealCount L p ≤ Module.finrank ℚ L := by
  rw [idealCount_eq_ncard hp]
  exact le_trans (Set.ncard_le_ncard (fun P hP => hP.1) (finite_primesOver_int hp))
    (ncard_primesOver_le_finrank_int hp)

/-- ★★★**`p` が完全分解すれば `a_L(p) = [L:ℚ]`**。 -/
theorem idealCount_eq_finrank_of_splitsCompletely [IsGalois ℚ L] {p : ℕ} (hp : p.Prime)
    (h : (primesOver (span {(p : ℤ)}) (𝓞 L)).ncard = Module.finrank ℚ L) :
    idealCount L p = Module.finrank ℚ L := by
  haveI : (span {(p : ℤ)}).IsMaximal := span_intCast_isMaximal hp
  rw [idealCount_eq_ncard hp]
  have hall := (ncard_primesOver_eq_finrank_iff_int (L := L) hp).mp h
  have hset : {P : Ideal (𝓞 L) | P ∈ primesOver (span {(p : ℤ)}) (𝓞 L)
      ∧ inertiaDeg (span {(p : ℤ)}) P = 1} = primesOver (span {(p : ℤ)}) (𝓞 L) := by
    ext P
    refine ⟨fun hh => hh.1, fun hP => ⟨hP, ?_⟩⟩
    have hPmem : P ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L) := by
      rw [← Finset.mem_coe, IsDedekindDomain.coe_primesOverFinset (span_intCast_ne_bot hp) (𝓞 L)]
      exact hP
    have hef := hall P hPmem
    exact Nat.eq_one_of_mul_eq_one_right (by rwa [mul_comm] at hef)
  rw [hset, h]

/-- ★★★**`p` が不分岐で完全分解しなければ `a_L(p) = 0`**(`L/ℚ` が Galois)。

★もし剰余次数 1 の素が 1 つでもあれば、不分岐と合わせて
`ncard_primesOver_eq_finrank_of_int` が完全分解を与えて矛盾する。 -/
theorem idealCount_eq_zero_of_unramified_of_not_split [IsGalois ℚ L] {p : ℕ} (hp : p.Prime)
    (hunram : ∀ P : Ideal (𝓞 L), P.IsPrime → P.LiesOver (span {(p : ℤ)}) →
      ramificationIdx (span {(p : ℤ)}) P = 1)
    (hnot : (primesOver (span {(p : ℤ)}) (𝓞 L)).ncard ≠ Module.finrank ℚ L) :
    idealCount L p = 0 := by
  haveI : Fact p.Prime := ⟨hp⟩
  by_contra hne
  rw [idealCount_eq_ncard hp] at hne
  obtain ⟨P, hP, hf⟩ := Set.nonempty_of_ncard_ne_zero hne
  haveI : P.IsPrime := hP.1
  haveI : P.LiesOver (span {(p : ℤ)}) := hP.2
  exact hnot (ncard_primesOver_eq_finrank_of_int P (hunram P hP.1 hP.2) hf)

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.4, (iv)` の密度論法の代数側。 -/
def idealCount_eq_ncard.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — a_L(p) = #{P ∣ p : f(P) = 1}",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
