import ABC3.Found.GaloisRep.TateQExp

/-!
# Galois (G6) 第 111 ブロック —— **★★★★★約数和の形と `s₁(q)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★古典的な形にする

第 110 ブロックで片側の尾の `qⁿ` 係数が
`∑_{m<n} [(m+1)|n]·(n/(m+1))·u^{n/(m+1)}` と出た。
★これは `range n` 上の約数フィルタだが、**`n.divisors` 上の和**に書き直せる:

    ∑_{m<n} [(m+1)|n]·f(m+1) = ∑_{e|n} f(e)          (`sum_range_ite_dvd`)

★★さらに約数の対応 `e ↔ n/e`(mathlib の `Nat.sum_div_divisors`)で

    ∑_{m≥1} f(qᵐu) = ∑_{n≥1} qⁿ·(∑_{d|n} d·u^d)
    ∑_{m≥1} g(qᵐu) = ∑_{n≥1} qⁿ·(∑_{d|n} C(d,2)·u^d)

——★★★**古典的な q 展開そのもの**である。

## ★★★★★★そして `u = 1` で `s₁(q)` が出る

`u = 1` を入れると `∑_{d|n} d = σ₁(n)` なので

    ∑_{m≥1} f(qᵐ) = ∑_{n≥1} σ₁(n) qⁿ = s₁(q)        (`tateXtail_one`)

★`s₁(q)` は `X` の定義に現れる正規化定数(`− 2 s₁(q)`)であった。
★★**級数の側と `ℤ⟦q⟧` の側が同じものを指していることが確かめられた**
——第 94 ブロックの `sigmaSeries` と第 105–110 の `adicSum` が噛み合った。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `sum_range_ite_dvd` | ★★★`range n` の約数フィルタ = `n.divisors` 上の和 |
| `sum_range_div_eq_sigma` | ★★★`∑_{d|n} (n/d) = σ₁(n)` |
| `sum_tateXcoef_divisors` | ★★★★`qⁿ` の係数 = `∑_{d|n} d·u^d` |
| `tateXtail_eq_divisorSum` | ★★★★★**古典的な q 展開** |
| `tateXtail_one` | ★★★★★★**`u = 1` の尾は `s₁(q)`** |
-/

namespace ABC3.Found.GaloisRep

/-! ## ★★★約数フィルタ -/

/-- ★★★**`range n` 上の約数フィルタは `n.divisors` 上の和である**。 -/
theorem sum_range_ite_dvd {M : Type} [AddCommMonoid M] (n : ℕ) (f : ℕ → M) :
    ∑ m ∈ Finset.range n, (if (m + 1) ∣ n then f (m + 1) else 0)
      = ∑ e ∈ n.divisors, f e := by
  have h1 : ∑ e ∈ Finset.Ico 1 (n + 1), (if e ∣ n then f e else 0)
      = ∑ m ∈ Finset.range n, (if (m + 1) ∣ n then f (m + 1) else 0) := by
    rw [Finset.sum_Ico_eq_sum_range]
    simp only [Nat.add_sub_cancel]
    exact Finset.sum_congr rfl (fun m _ => by rw [add_comm])
  rw [← h1, ← Finset.sum_filter]
  rfl

/-- ★★★**`∑_{d|n} (n/d) = σ₁(n)`**。 -/
theorem sum_range_div_eq_sigma (n : ℕ) :
    ∑ m ∈ Finset.range n, (if (m + 1) ∣ n then n / (m + 1) else 0)
      = ArithmeticFunction.sigma 1 n := by
  rw [sum_range_ite_dvd n (fun e => n / e), Nat.sum_div_divisors n (fun x => x),
    ArithmeticFunction.sigma_apply]
  exact Finset.sum_congr rfl (fun d _ => (pow_one d).symm)

/-! ## ★★★★約数和の形 -/

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★**`qⁿ` の係数は約数和 `∑_{d|n} d·u^d` である**。 -/
theorem sum_tateXcoef_divisors (u q : R) (n : ℕ) :
    ∑ m ∈ Finset.range n, tateXcoef u q m n
      = q ^ n * ∑ d ∈ n.divisors, (d : R) * u ^ d := by
  rw [sum_tateXcoef, sum_range_ite_dvd n (fun e => ((n / e : ℕ) : R) * u ^ (n / e))]
  congr 1
  exact Nat.sum_div_divisors n (fun d => (d : R) * u ^ d)

/-- ★★★★**`Y` 側**: `qⁿ` の係数は `∑_{d|n} C(d,2)·u^d`。 -/
theorem sum_tateYcoef_divisors (u q : R) (n : ℕ) :
    ∑ m ∈ Finset.range n, tateYcoef u q m n
      = q ^ n * ∑ d ∈ n.divisors, ((d.choose 2 : ℕ) : R) * u ^ d := by
  rw [sum_tateYcoef, sum_range_ite_dvd n (fun e => (((n / e).choose 2 : ℕ) : R) * u ^ (n / e))]
  congr 1
  exact Nat.sum_div_divisors n (fun d => ((d.choose 2 : ℕ) : R) * u ^ d)

/-- ★★★★★**片側の尾の古典的な q 展開**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tateXtail_eq_divisorSum [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, (d : R) * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  rw [tateXtail_eq_adicSum u q hq]
  exact adicSum_congr _ _ (fun n => sum_tateXcoef_divisors u q n)

theorem tateYtail_eq_divisorSum [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail u q hq
      = adicSum (fun n => q ^ n * ∑ d ∈ n.divisors, ((d.choose 2 : ℕ) : R) * u ^ d)
          (fun n => Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)) := by
  rw [tateYtail_eq_adicSum u q hq]
  exact adicSum_congr _ _ (fun n => sum_tateYcoef_divisors u q n)

/-! ## ★★★★★★`u = 1` の尾は `s₁(q)` -/

/-- ★★★★★★**`u = 1` の尾は `s₁(q)` そのものである**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★`s₁(q)` は `X` の定義に現れる正規化定数(`− 2 s₁(q)`)であった。
★★級数の側(第 105–110 の `adicSum`)と `ℤ⟦q⟧` の側(第 94 の `sigmaSeries`)が
同じものを指していることが確かめられた。 -/
theorem tateXtail_one [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    tateXtail 1 q hq = evalAdic (sigmaSeries 1) q hq := by
  rw [tateXtail_eq_adicSum 1 q hq, evalAdic_eq_adicSum]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [sum_tateXcoef]
  have hsum : ∑ m ∈ Finset.range n,
      (if (m + 1) ∣ n then ((n / (m + 1) : ℕ) : R) * (1 : R) ^ (n / (m + 1)) else 0)
      = ((ArithmeticFunction.sigma 1 n : ℕ) : R) := by
    rw [← sum_range_div_eq_sigma n]
    push_cast
    refine Finset.sum_congr rfl (fun m _ => ?_)
    by_cases h : (m + 1) ∣ n
    · rw [if_pos h, if_pos h]
      simp
    · rw [if_neg h, if_neg h]
  have hc : ((PowerSeries.coeff n (sigmaSeries 1) : ℤ) : R)
      = ((ArithmeticFunction.sigma 1 n : ℕ) : R) := by
    rw [coeff_sigmaSeries]
    by_cases h : n = 0
    · subst h
      simp
    · rw [if_neg h]
      push_cast
      ring
  rw [hsum, hc]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def tateXtail_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 級数の約数和展開と s₁(q) との一致)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
