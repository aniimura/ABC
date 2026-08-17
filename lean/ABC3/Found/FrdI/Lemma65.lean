import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Data.Nat.Factorization.Basic
import ABC3.Meta.Claim

/-!
# [FrdI] Lemma 6.5, (i) —— 素数の対数は ℚ 上一次独立

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> Lemma 6.5.

原文 (FrdI p.116):
> linearly independent over Q.

## ★★この補題の 2 主張は**難度がまったく違う**(測定、2026-08-17)

| 条 | 内容 | 状況 |
|---|---|---|
| (i) | `log p`(`p` 素数)は ℚ 上一次独立 | ★**本ファイルで証明**。原文の言うとおり「`ℤ` が UFD であること」だけで出る |
| (ii) | 相異なる 6 素数についての比の関係式が不可能 | ★★**Lang の定理**が要る。mathlib に無い |

原文 (FrdI p.117):
> Assertion (i) is a formal consequence of the fact that Z is a unique factor-

★★**(ii) は超越論**である —— 原文は「a consequence of a theorem of Lang」と書き、
Baker の本の p.119 を指す。★これは **six exponentials theorem** であり、
mathlib には**影も無い**(`sixExponentials` も `Baker` も 0 件、2026-08-17 に検索)。
★したがって (ii) は `Gap/` 側の案件で、**本ファイルには置かない**。

## ★(i) の証明の骨

1. 分母を払って**係数を整数にする**(`N := ∏ (g p).den`)
2. 正の部分と負の部分に分け、**`log` の外に出して積にする**
   —— `∑ a_p log p = log (∏ p^{a_p})`
3. `log` は正の実数の上で単射なので `∏ p^{a_p} = ∏ p^{b_p}`(自然数として)
4. ★★**素因数分解の一意性**(`Nat.factorization`)で `a_p = b_p`

★**唯一の技術的な段は 1**(有理係数を整数係数に直す)で、
数学的な中身は 4 の一意分解だけである。

★★**添字は最初から `Nat.Primes` に取る** —— `Finset ℕ` に落として
「元が相異なる素数である」を持ち回ると、像・単射性の配管が増える(2026-08-17 に実測)。
-/

namespace ABC3.Found.FrdI

open Finset

/-! ## ★段 4 —— 素因数分解の一意性 -/

/-- ★**相異なる素数の冪の積は指数を決める**。

★`Nat.factorization` を素数 `q` で評価すると、添字が相異なるので
`Finsupp.single` の対角成分だけが残る。 -/
theorem prod_prime_pow_inj {t : Finset Nat.Primes} (a b : Nat.Primes → ℕ)
    (h : ∏ p ∈ t, (p : ℕ) ^ a p = ∏ p ∈ t, (p : ℕ) ^ b p) : ∀ q ∈ t, a q = b q := by
  have hfac : ∀ c : Nat.Primes → ℕ, (∏ p ∈ t, (p : ℕ) ^ c p).factorization
      = ∑ p ∈ t, Finsupp.single (p : ℕ) (c p) := by
    intro c
    rw [Nat.factorization_prod (fun p _ => pow_ne_zero _ p.2.ne_zero)]
    exact Finset.sum_congr rfl (fun p _ => Nat.Prime.factorization_pow p.2)
  have hsingle : ∀ (c : Nat.Primes → ℕ) (q : Nat.Primes), q ∈ t →
      (∑ p ∈ t, Finsupp.single (p : ℕ) (c p)) (q : ℕ) = c q := by
    intro c q hq
    rw [Finset.sum_apply']
    rw [Finset.sum_eq_single q]
    · simp
    · intro p _ hpq
      refine Finsupp.single_eq_of_ne ?_
      exact fun hc => hpq (Subtype.ext hc.symm)
    · intro hqt
      exact absurd hq hqt
  intro q hq
  have hq' := congrArg (fun m : ℕ => m.factorization (q : ℕ)) h
  simp only [hfac] at hq'
  rw [hsingle a q hq, hsingle b q hq] at hq'
  exact hq'

/-! ## ★段 2・3 —— 整数係数の場合 -/

/-- ★★**整数係数版** —— 相異なる素数の `log` の整数一次結合が `0` なら係数はすべて `0`。

★正の部分と負の部分に分けて `log` の外に出し、`log` の単射性(正の実数上)と
素因数分解の一意性で潰す。 -/
theorem int_sum_log_eq_zero {t : Finset Nat.Primes} (n : Nat.Primes → ℤ)
    (h : ∑ p ∈ t, (n p : ℝ) * Real.log (p : ℕ) = 0) : ∀ q ∈ t, n q = 0 := by
  set a : Nat.Primes → ℕ := fun p => (n p).toNat with ha
  set b : Nat.Primes → ℕ := fun p => (-(n p)).toNat with hb
  have hab : ∀ p, (a p : ℤ) - (b p : ℤ) = n p := fun p => Int.toNat_sub_toNat_neg (n p)
  -- ★段 2: 正の部分と負の部分に分ける
  have hsplit : ∑ p ∈ t, (a p : ℝ) * Real.log (p : ℕ)
      = ∑ p ∈ t, (b p : ℝ) * Real.log (p : ℕ) := by
    have hz : ∑ p ∈ t, ((a p : ℝ) * Real.log (p : ℕ) - (b p : ℝ) * Real.log (p : ℕ)) = 0 := by
      refine Eq.trans (Finset.sum_congr rfl (fun p _ => ?_)) h
      have hc : ((a p : ℝ) - (b p : ℝ)) = (n p : ℝ) := by exact_mod_cast hab p
      rw [← sub_mul, hc]
    rw [Finset.sum_sub_distrib] at hz
    linarith [hz]
  -- ★段 2': `log` の外へ出す
  have hlog : ∀ c : Nat.Primes → ℕ,
      Real.log ((∏ p ∈ t, (p : ℕ) ^ c p : ℕ) : ℝ)
        = ∑ p ∈ t, (c p : ℝ) * Real.log (p : ℕ) := by
    intro c
    push_cast
    rw [Real.log_prod (fun p _ => pow_ne_zero _ (Nat.cast_ne_zero.mpr p.2.ne_zero))]
    exact Finset.sum_congr rfl (fun p _ => Real.log_pow _ _)
  have hpos : ∀ c : Nat.Primes → ℕ, (0 : ℝ) < ((∏ p ∈ t, (p : ℕ) ^ c p : ℕ) : ℝ) :=
    fun c => Nat.cast_pos.mpr (Finset.prod_pos (fun p _ => pow_pos p.2.pos _))
  -- ★段 3: `log` は正の実数上で単射
  have hprod : ((∏ p ∈ t, (p : ℕ) ^ a p : ℕ) : ℝ)
      = ((∏ p ∈ t, (p : ℕ) ^ b p : ℕ) : ℝ) := by
    have hl : Real.log ((∏ p ∈ t, (p : ℕ) ^ a p : ℕ) : ℝ)
        = Real.log ((∏ p ∈ t, (p : ℕ) ^ b p : ℕ) : ℝ) := by
      rw [hlog a, hlog b]; exact hsplit
    calc ((∏ p ∈ t, (p : ℕ) ^ a p : ℕ) : ℝ)
        = Real.exp (Real.log ((∏ p ∈ t, (p : ℕ) ^ a p : ℕ) : ℝ)) := (Real.exp_log (hpos a)).symm
      _ = Real.exp (Real.log ((∏ p ∈ t, (p : ℕ) ^ b p : ℕ) : ℝ)) := by rw [hl]
      _ = ((∏ p ∈ t, (p : ℕ) ^ b p : ℕ) : ℝ) := Real.exp_log (hpos b)
  have hprodN : (∏ p ∈ t, (p : ℕ) ^ a p) = (∏ p ∈ t, (p : ℕ) ^ b p) := Nat.cast_injective hprod
  -- ★段 4: 一意分解
  intro q hq
  have hqab := prod_prime_pow_inj a b hprodN q hq
  have hd := hab q
  rw [hqab] at hd
  simpa using hd.symm

/-! ## ★段 1 —— 有理係数を整数係数に直す -/

/-- ★**分母を払える** —— `N` が `q` の分母で割り切れるなら `N * q` は整数。 -/
theorem exists_int_of_den_dvd {q : ℚ} {N : ℕ} (hdvd : q.den ∣ N) :
    ∃ z : ℤ, (z : ℚ) = (N : ℚ) * q := by
  obtain ⟨k, hk⟩ := hdvd
  have hd : ((q.den : ℚ)) ≠ 0 := Nat.cast_ne_zero.mpr q.den_ne_zero
  have h : (q.num : ℚ) = q * (q.den : ℚ) := by
    rw [← div_eq_iff hd]
    exact Rat.num_div_den q
  refine ⟨(k : ℤ) * q.num, ?_⟩
  rw [hk]
  push_cast
  rw [h]
  ring

/-- ★★★**[FrdI] Lemma 6.5, (i)** —— 素数の対数は ℚ 上一次独立。

原文 (FrdI p.116):
> linearly independent over Q.

★分母の積 `N` を掛けて整数係数に直し、整数版(`int_sum_log_eq_zero`)へ渡す。 -/
theorem log_primes_linearIndependent :
    LinearIndependent ℚ (fun p : Nat.Primes => Real.log (p : ℕ)) := by
  rw [linearIndependent_iff']
  intro t g hsum i hi
  set N : ℕ := ∏ p ∈ t, (g p).den with hN
  have hNpos : 0 < N := Finset.prod_pos (fun p _ => (g p).pos)
  set n : Nat.Primes → ℤ := fun p => ((N : ℚ) * g p).num with hn
  have hcast : ∀ p ∈ t, ((n p : ℤ) : ℚ) = (N : ℚ) * g p := by
    intro p hp
    obtain ⟨z, hz⟩ := exists_int_of_den_dvd
      (Finset.dvd_prod_of_mem (fun q : Nat.Primes => (g q).den) hp)
    have hnum : ((N : ℚ) * g p).num = z := by rw [← hz]; exact Rat.num_intCast z
    show (((N : ℚ) * g p).num : ℚ) = (N : ℚ) * g p
    rw [hnum, hz]
  have hg : ∑ p ∈ t, (g p : ℝ) * Real.log (p : ℕ) = 0 := by
    refine Eq.trans (Finset.sum_congr rfl (fun p _ => ?_)) hsum
    exact (Rat.smul_def (g p) (Real.log (p : ℕ))).symm
  have hsum' : ∑ p ∈ t, ((n p : ℤ) : ℝ) * Real.log (p : ℕ) = 0 := by
    have hstep : ∑ p ∈ t, ((n p : ℤ) : ℝ) * Real.log (p : ℕ)
        = (N : ℝ) * ∑ p ∈ t, (g p : ℝ) * Real.log (p : ℕ) := by
      rw [Finset.mul_sum]
      refine Finset.sum_congr rfl (fun p hp => ?_)
      have h1 : (((n p : ℤ)) : ℝ) = (N : ℝ) * ((g p : ℚ) : ℝ) := by
        have := congrArg (fun x : ℚ => (x : ℝ)) (hcast p hp)
        push_cast at this ⊢
        linarith [this]
      rw [h1]
      ring
    rw [hstep, hg, mul_zero]
  have hzero := int_sum_log_eq_zero n hsum' i hi
  have h0 : (N : ℚ) * g i = 0 := by rw [← hcast i hi, hzero]; simp
  rcases mul_eq_zero.mp h0 with h | h
  · exact absurd h (by exact_mod_cast hNpos.ne')
  · exact h

/-! ## ★(ii) は本ファイルに置かない

原文 (FrdI p.117):
> Assertion (i) is a formal consequence of the fact that Z is a unique factor-

★★**(ii) は Lang の定理**を要し、mathlib に無い。
★`Gap/FrdI/Section6.lean` に**分類 ②(missingMath)**として記録する。 -/

/-- ★★**locator** —— `Lemma 6.5, (i)` のみ。★★**条つき**である。 -/
def log_primes_linearIndependent.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116, item := "Lemma 6.5, (i)",
    sectionId := "frdi-lemma-6-5" }

end ABC3.Found.FrdI
