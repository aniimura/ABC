import ABC3.Found.GaloisRep.AdicFubini

/-!
# Galois (G6) 第 110 ブロック —— **★★★★★片側の尾の q 展開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★第 108・109 を合わせる

第 108 ブロック: `f(t) = ∑_{d≥0} d·tᵈ`。
第 109 ブロック: 添字の付け替えと二重和の入れ替え。

★これを合わせると片側の尾が **`q` の冪ごとに**まとまる:

    ∑_{m≥1} f(qᵐu) = ∑_{m≥1} ∑_{d≥1} d·q^{md}·u^d
                   = ∑_{n≥1} qⁿ·(∑_{e|n} (n/e)·u^{n/e})
                   = ∑_{n≥1} qⁿ·(∑_{d|n} d·u^d)

★★最後の等式は約数の対応 `e ↔ n/e` である(本ブロックでは
`∑_{m<n} [(m+1)|n] (n/(m+1))·u^{n/(m+1)}` の形まで出す)。

## ★★★★手順

| 段 | 内容 |
|---|---|
| 1 | `tateXcoef u q m n` —— `f(q^{m+1}u)` の `qⁿ` への寄与 |
| 2 | `tateXterm_pow_eq` —— 第 108 + 添字の付け替え |
| 3 | `tateXtail_eq_adicSum` —— 第 109 の Fubini |
| 4 | `sum_tateXcoef` —— 係数から `qⁿ` を括り出す |

★`Y` の側も同じ形で通る(`C(d,2)` に置き換わるだけ)。
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★`X` の側 -/

/-- ★`f(q^{m+1}u)` の `qⁿ` 係数への寄与。 -/
noncomputable def tateXcoef (u q : R) (m n : ℕ) : R :=
  if (m + 1) ∣ n then ((n / (m + 1) : ℕ) : R) * (q ^ (m + 1) * u) ^ (n / (m + 1)) else 0

theorem tateXpow_mem {u q : R} (hq : q ∈ I) (m d : ℕ) :
    (d : R) * (q ^ (m + 1) * u) ^ d ∈ I ^ ((m + 1) * d) := by
  have h1 : (q ^ (m + 1) * u) ^ d = q ^ ((m + 1) * d) * u ^ d := by
    rw [mul_pow, ← pow_mul]
  rw [h1]
  exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq _))

theorem tateXcoef_mem {u q : R} (hq : q ∈ I) (m n : ℕ) : tateXcoef u q m n ∈ I ^ n :=
  reindex_mem (m + 1) (tateXpow_mem hq m) n

theorem tateXcoef_vanish (u q : R) (m n : ℕ) (hn : n ≤ m) : tateXcoef u q m n = 0 := by
  rw [tateXcoef]
  by_cases h : (m + 1) ∣ n
  · rw [if_pos h]
    have hn0 : n = 0 := by
      rcases Nat.eq_zero_or_pos n with h0 | h0
      · exact h0
      · exact absurd (Nat.le_of_dvd h0 h) (by omega)
    subst hn0
    simp
  · rw [if_neg h]

theorem pow_succ_mul_mem {u q : R} (hq : q ∈ I) (m : ℕ) : q ^ (m + 1) * u ∈ I := by
  refine Ideal.mul_mem_right _ _ ?_
  rw [pow_succ']
  exact Ideal.mul_mem_right _ _ hq

/-- ★★★**第 108 + 添字の付け替え**——`f(q^{m+1}u)` を `q` の冪で並べ直す。 -/
theorem tateXterm_pow_eq [IsAdicComplete I R] {u q : R} (hq : q ∈ I) (m : ℕ) :
    tateXterm (q ^ (m + 1) * u) = adicSum (tateXcoef u q m) (tateXcoef_mem hq m) := by
  rw [tateXterm_eq_adicSum (pow_succ_mul_mem hq m)]
  exact adicSum_reindex_mul (m + 1) (by omega) _ (tateXpow_mem hq m)

/-- ★★★★★**片側の尾の q 展開**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXtail_eq_adicSum [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail u q hq
      = adicSum (fun n => ∑ m ∈ Finset.range n, tateXcoef u q m n)
          (diag_mem _ (fun m n => tateXcoef_mem hq m n)) := by
  rw [tateXtail,
    adicSum_congr (tateXtail_aux hq)
      (adicSum_mem_pow_of_vanish _ (fun m n => tateXcoef_mem hq m n)
        (fun m n h => tateXcoef_vanish u q m n h))
      (fun m => tateXterm_pow_eq hq m)]
  exact adicSum_fubini (tateXcoef u q) (fun m n => tateXcoef_mem hq m n)
    (fun m n h => tateXcoef_vanish u q m n h)

/-- ★★係数から `qⁿ` を括り出す。 -/
theorem tateXcoef_eq (u q : R) (m n : ℕ) (h : (m + 1) ∣ n) :
    tateXcoef u q m n = ((n / (m + 1) : ℕ) : R) * q ^ n * u ^ (n / (m + 1)) := by
  rw [tateXcoef, if_pos h, mul_pow, ← pow_mul, Nat.mul_div_cancel' h]
  ring

/-- ★★★**`qⁿ` の係数は約数の和である**。 -/
theorem sum_tateXcoef (u q : R) (n : ℕ) :
    ∑ m ∈ Finset.range n, tateXcoef u q m n
      = q ^ n * ∑ m ∈ Finset.range n,
          (if (m + 1) ∣ n then ((n / (m + 1) : ℕ) : R) * u ^ (n / (m + 1)) else 0) := by
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  by_cases h : (m + 1) ∣ n
  · rw [tateXcoef_eq u q m n h, if_pos h]
    ring
  · rw [tateXcoef, if_neg h, if_neg h]
    ring

/-! ## ★`Y` の側 -/

/-- ★`g(q^{m+1}u)` の `qⁿ` 係数への寄与。 -/
noncomputable def tateYcoef (u q : R) (m n : ℕ) : R :=
  if (m + 1) ∣ n then (((n / (m + 1)).choose 2 : ℕ) : R) * (q ^ (m + 1) * u) ^ (n / (m + 1))
  else 0

theorem tateYpow_mem {u q : R} (hq : q ∈ I) (m d : ℕ) :
    ((d.choose 2 : ℕ) : R) * (q ^ (m + 1) * u) ^ d ∈ I ^ ((m + 1) * d) := by
  have h1 : (q ^ (m + 1) * u) ^ d = q ^ ((m + 1) * d) * u ^ d := by
    rw [mul_pow, ← pow_mul]
  rw [h1]
  exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq _))

theorem tateYcoef_mem {u q : R} (hq : q ∈ I) (m n : ℕ) : tateYcoef u q m n ∈ I ^ n :=
  reindex_mem (m + 1) (tateYpow_mem hq m) n

theorem tateYcoef_vanish (u q : R) (m n : ℕ) (hn : n ≤ m) : tateYcoef u q m n = 0 := by
  rw [tateYcoef]
  by_cases h : (m + 1) ∣ n
  · rw [if_pos h]
    have hn0 : n = 0 := by
      rcases Nat.eq_zero_or_pos n with h0 | h0
      · exact h0
      · exact absurd (Nat.le_of_dvd h0 h) (by omega)
    subst hn0
    simp
  · rw [if_neg h]

theorem tateYterm_pow_eq [IsAdicComplete I R] {u q : R} (hq : q ∈ I) (m : ℕ) :
    tateYterm (q ^ (m + 1) * u) = adicSum (tateYcoef u q m) (tateYcoef_mem hq m) := by
  rw [tateYterm_eq_adicSum (pow_succ_mul_mem hq m)]
  exact adicSum_reindex_mul (m + 1) (by omega) _ (tateYpow_mem hq m)

/-- ★★★★★**`Y` 側の片側の尾の q 展開**。 -/
theorem tateYtail_eq_adicSum [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail u q hq
      = adicSum (fun n => ∑ m ∈ Finset.range n, tateYcoef u q m n)
          (diag_mem _ (fun m n => tateYcoef_mem hq m n)) := by
  rw [tateYtail,
    adicSum_congr (tateYtail_aux hq)
      (adicSum_mem_pow_of_vanish _ (fun m n => tateYcoef_mem hq m n)
        (fun m n h => tateYcoef_vanish u q m n h))
      (fun m => tateYterm_pow_eq hq m)]
  exact adicSum_fubini (tateYcoef u q) (fun m n => tateYcoef_mem hq m n)
    (fun m n h => tateYcoef_vanish u q m n h)

theorem tateYcoef_eq (u q : R) (m n : ℕ) (h : (m + 1) ∣ n) :
    tateYcoef u q m n = (((n / (m + 1)).choose 2 : ℕ) : R) * q ^ n * u ^ (n / (m + 1)) := by
  rw [tateYcoef, if_pos h, mul_pow, ← pow_mul, Nat.mul_div_cancel' h]
  ring

theorem sum_tateYcoef (u q : R) (n : ℕ) :
    ∑ m ∈ Finset.range n, tateYcoef u q m n
      = q ^ n * ∑ m ∈ Finset.range n,
          (if (m + 1) ∣ n then (((n / (m + 1)).choose 2 : ℕ) : R) * u ^ (n / (m + 1)) else 0) := by
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl (fun m _ => ?_)
  by_cases h : (m + 1) ∣ n
  · rw [tateYcoef_eq u q m n h, if_pos h]
    ring
  · rw [tateYcoef, if_neg h, if_neg h]
    ring

/-! ## ★出典の紐付け(`.src`) -/

def tateXtail_eq_adicSum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 級数の片側の尾の q 展開——葉 (b))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
