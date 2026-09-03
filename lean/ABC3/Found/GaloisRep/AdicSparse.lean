/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuCharSum
import ABC3.Found.GaloisRep.AdicSeries

/-!
# Galois (G6) 第 857 ブロック —— **★★★★★★★★疎な級数と `σ_k` の `l` 倍**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か

残る葉 `sum_mu_d2xpair` の右辺には `evalAdic (sigmaSeries 3) (q^l)` が出る。
`∑_{i≠0} T_{D²f}(ζ^i)` は `q` の級数として

    `adicSum (n ↦ q^n · (l⁴σ₃(n/l)[l∣n] − σ₃(n)))`

であり、その前半は `l∤n` で消える**疎な級数**である。★これを `q^l` の級数に
付け替えるのが `adicSum_sparse` であり、`∑_{d∣N, l∣d} d^k = l^k σ_k(N/l)` が
`sum_divisors_dvd_pow` である。

| 定理 | 内容 |
|---|---|
| `partialSum_sparse` | `l∤n` で消える級数の部分和は `q^l` 側の部分和に等しい |
| `adicSum_sparse` | ★★★**疎な級数の和は付け替えられる** |
| `sum_divisors_dvd_pow` | ★★`∑_{d∣N, l∣d} d^k = l^k·σ_k(N/l)`（`l∣N` のとき） |
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★疎な級数 -/

theorem partialSum_sparse (l N : ℕ) (hl : 0 < l) (a b : ℕ → R)
    (h0 : ∀ n, ¬ l ∣ n → a n = 0) (hab : ∀ m, a (l * m) = b m) :
    partialSum a (l * N) = partialSum b N := by
  classical
  simp only [partialSum]
  rw [← Finset.sum_filter_add_sum_filter_not (Finset.range (l * N)) (fun n => l ∣ n)]
  have hz : ∑ n ∈ (Finset.range (l * N)).filter (fun n => ¬ l ∣ n), a n = 0 :=
    Finset.sum_eq_zero (fun n hn => h0 n (Finset.mem_filter.1 hn).2)
  rw [hz, add_zero]
  refine Finset.sum_nbij' (fun n => n / l) (fun m => l * m) ?_ ?_ ?_ ?_ ?_
  · intro n hn
    rw [Finset.mem_filter] at hn
    have hnl : n < l * N := Finset.mem_range.1 hn.1
    refine Finset.mem_range.2 ?_
    exact Nat.div_lt_of_lt_mul (by linarith [hnl])
  · intro m hm
    have hmN : m < N := Finset.mem_range.1 hm
    refine Finset.mem_filter.2 ⟨Finset.mem_range.2 ?_, Dvd.intro m rfl⟩
    exact (Nat.mul_lt_mul_left hl).2 hmN
  · intro n hn
    rw [Finset.mem_filter] at hn
    exact Nat.mul_div_cancel' hn.2
  · intro m _
    exact Nat.mul_div_cancel_left _ hl
  · intro n hn
    rw [Finset.mem_filter] at hn
    have := hab (n / l)
    rw [Nat.mul_div_cancel' hn.2] at this
    exact this

/-- ★★★★★★**疎な級数の和は `q^l` 側に付け替えられる**。 -/
theorem adicSum_sparse [IsAdicComplete I R] (l : ℕ) (hl : 0 < l) (a b : ℕ → R)
    (ha : ∀ n, a n ∈ I ^ n) (hb : ∀ n, b n ∈ I ^ n)
    (h0 : ∀ n, ¬ l ∣ n → a n = 0) (hab : ∀ m, a (l * m) = b m) :
    adicSum a ha = adicSum b hb := by
  refine (IsHausdorff.eq_iff_smodEq (I := I)).2 (fun N => ?_)
  rw [SModEq.sub_mem]
  have h1 := adicSum_spec a ha (l * N)
  have h2 := adicSum_spec b hb N
  rw [SModEq.sub_mem] at h1 h2
  have m1 : partialSum a (l * N) - adicSum a ha ∈ I ^ N := by
    refine Ideal.pow_le_pow_right (Nat.le_mul_of_pos_left N hl) ?_
    simpa using h1
  have m2 : partialSum b N - adicSum b hb ∈ I ^ N := by simpa using h2
  have hps : partialSum a (l * N) = partialSum b N := partialSum_sparse l N hl a b h0 hab
  have hsplit : adicSum a ha - adicSum b hb
      = (partialSum b N - adicSum b hb) - (partialSum a (l * N) - adicSum a ha) := by
    rw [hps]
    ring
  rw [hsplit]
  simpa using Submodule.sub_mem _ m2 m1

/-! ## ★★★★★`σ_k` の `l` 倍 -/

/-- ★★★**`∑_{d∣N, l∣d} d^k = l^k·σ_k(N/l)`**（`l∣N` のとき、そうでなければ `0`）。 -/
theorem sum_divisors_dvd_pow (k l N : ℕ) (hl : l.Prime) (hN : N ≠ 0) :
    ∑ d ∈ N.divisors.filter (fun d => l ∣ d), d ^ k
      = if l ∣ N then l ^ k * ∑ e ∈ (N / l).divisors, e ^ k else 0 := by
  classical
  by_cases hlN : l ∣ N
  · rw [if_pos hlN]
    obtain ⟨M, hM⟩ := hlN
    have hM0 : M ≠ 0 := by rintro rfl; exact hN (by simp [hM])
    have hNl : N / l = M := by rw [hM, Nat.mul_div_cancel_left _ hl.pos]
    rw [hNl, Finset.mul_sum]
    refine Finset.sum_nbij' (fun d => d / l) (fun e => l * e) ?_ ?_ ?_ ?_ ?_
    · intro d hd
      rw [Finset.mem_filter] at hd
      obtain ⟨hdN, hld⟩ := hd
      obtain ⟨c, hc⟩ := hld
      have hdvd : d ∣ N := (Nat.mem_divisors.1 hdN).1
      refine Nat.mem_divisors.2 ⟨?_, hM0⟩
      refine ⟨N / d, ?_⟩
      have hdl : d / l = c := by rw [hc, Nat.mul_div_cancel_left _ hl.pos]
      rw [hdl]
      have hmm : l * (c * (N / d)) = l * M := by
        rw [← hM, ← mul_assoc, ← hc, Nat.mul_div_cancel' hdvd]
      exact Nat.eq_of_mul_eq_mul_left hl.pos hmm.symm
    · intro e he
      have heM : e ∣ M := (Nat.mem_divisors.1 he).1
      refine Finset.mem_filter.2 ⟨Nat.mem_divisors.2 ⟨?_, hN⟩, Dvd.intro e rfl⟩
      rw [hM]
      exact mul_dvd_mul_left l heM
    · intro d hd
      rw [Finset.mem_filter] at hd
      exact Nat.mul_div_cancel' hd.2
    · intro e _
      exact Nat.mul_div_cancel_left _ hl.pos
    · intro d hd
      rw [Finset.mem_filter] at hd
      rw [← mul_pow, Nat.mul_div_cancel' hd.2]
  · rw [if_neg hlN]
    have hempty : N.divisors.filter (fun d => l ∣ d) = ∅ := by
      refine Finset.filter_eq_empty_iff.2 (fun {d} hd hld => ?_)
      exact hlN (hld.trans (Nat.mem_divisors.1 hd).1)
    rw [hempty, Finset.sum_empty]

/-! ## ★出典の紐付け(`.src`) -/

def adicSum_sparse.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(疎な級数の和は q^l 側に付け替えられる。★無条件)",
    sectionId := "genell-lemma-3-2" }

def sum_divisors_dvd_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_{d∣N, l∣d} d^k = l^k·σ_k(N/l)。★l は素数)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
