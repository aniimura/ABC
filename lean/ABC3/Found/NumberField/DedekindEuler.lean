/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.DedekindZeta
import Mathlib.NumberTheory.EulerProduct.Basic
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas

/-!
# Dedekind ζ の Euler 積(鎖 `cheb` の `cheb-euler`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何のためか

`Theorem 6.4, (iv)` が使う Tchebotarev 密度定理を、**類体論を使わずに**
「完全分解する素点の密度は `1/[L:K]`」の形で出すための一段目である
(`ResearchPaper/frdi-decomposition.json` の鎖 `cheb`)。

## ★★★mathlib の在庫(2026-08-25 実測)

* `NumberField.dedekindZeta K s = LSeries (fun n ↦ #{I : N I = n}) s` —— **在る**
* `NumberField.tendsto_sub_one_mul_dedekindZeta_nhdsGT`(類数公式・`s = 1` の留数)—— **在る**
* `EulerProduct.eulerProduct_hasProd`(互いに素で乗法的な `f` の Euler 積)—— **在る**
* ★**イデアル計数関数 `a_K(n) = #{I : N I = n}` の乗法性** —— **無い**。本ファイルで証明する。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `sup_eq_top_of_coprime_absNorm` | ノルムが互いに素な 2 つのイデアルは comaximal |
| `mul_sup_span_absNorm` | ★**復元** `I = I·J + (N I)`(Bézout 1 行) |
| `mul_split_eq` | ★**分解** `A = (A + (m))·(A + (n))` |
| `absNorm_sup_span_eq` | ★分解した各因子のノルムはちょうど `m`・`n` |
| `coprimeIdealEquiv` | ★★★★`{N I = m} × {N J = n} ≃ {N A = m·n}` |
| `idealCount_mul` | ★★★★★**計数関数は乗法的** |
| `sum_card_absNorm_isBigO` | 部分和は `O(n)`(mathlib の類数公式の漸近から) |
| `summable_zetaSummand` | `Re s > 1` で絶対収束 |
| `dedekindZeta_eulerProduct_hasProd` | ★★★★★★**Euler 積** |

★★要は **Bézout 1 本**である —— `u·m + v·n = 1` を `x ∈ I` に掛けると
`x = (u·m)·x + v·(x·n)` で、前者は `(m)`、後者は `I·J` に入る。
これで `I = I·J + (m)` となり、積からの復元ができる。
-/

namespace ABC3.Found.NF

open NumberField Ideal Filter Finset Topology Asymptotics Complex
open scoped nonZeroDivisors

/-! ## ★1. ノルムが互いに素なイデアルの積と分解 -/

section Dedekind

variable {R : Type*} [CommRing R] [IsDedekindDomain R] [Module.Free ℤ R] [Module.Finite ℤ R]

omit [Module.Finite ℤ R] in
/-- ★ノルムが互いに素なら 2 つのイデアルは comaximal。 -/
theorem sup_eq_top_of_coprime_absNorm {I J : Ideal R}
    (h : Nat.Coprime (absNorm I) (absNorm J)) : I ⊔ J = ⊤ :=
  absNorm_eq_one_iff.mp (Nat.eq_one_of_dvd_coprimes h
    (absNorm_dvd_absNorm_of_le le_sup_left) (absNorm_dvd_absNorm_of_le le_sup_right))

omit [IsDedekindDomain R] [Module.Free ℤ R] [Module.Finite ℤ R] in
/-- ★ℕ の互いに素から `R` の元としての `IsCoprime` を得る。 -/
theorem isCoprime_cast_of_coprime {m n : ℕ} (h : Nat.Coprime m n) :
    IsCoprime (m : R) (n : R) := by
  have hz : IsCoprime (m : ℤ) (n : ℤ) := Int.isCoprime_iff_gcd_eq_one.mpr (by simpa using h)
  simpa using hz.map (Int.castRingHom R)

/-- ★★**互いに素な部分の復元** —— `I = I·J + (N I)`。

★中身は Bézout 1 本 —— `u·(N I) + v·(N J) = 1` を `x ∈ I` に掛けると
`x = (u·(N I))·x + v·(x·(N J))` で、前者は `((N I))`、後者は `I·J` に入る。 -/
theorem mul_sup_span_absNorm {I J : Ideal R}
    (h : Nat.Coprime (absNorm I) (absNorm J)) :
    I * J ⊔ Ideal.span {(absNorm I : R)} = I := by
  refine le_antisymm (sup_le Ideal.mul_le_right (Ideal.span_le.mpr ?_)) ?_
  · rintro x hx
    rw [Set.mem_singleton_iff] at hx
    exact hx ▸ Ideal.absNorm_mem I
  · intro x hx
    obtain ⟨u, v, huv⟩ := isCoprime_cast_of_coprime (R := R) h
    have hx1 : x = u * (absNorm I : R) * x + v * (x * (absNorm J : R)) :=
      calc x = (u * (absNorm I : R) + v * (absNorm J : R)) * x := by rw [huv, one_mul]
        _ = u * (absNorm I : R) * x + v * (x * (absNorm J : R)) := by ring
    rw [hx1]
    refine Ideal.add_mem _ (Ideal.mem_sup_right ?_) (Ideal.mem_sup_left ?_)
    · exact Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ (Ideal.subset_span rfl))
    · exact Ideal.mul_mem_left _ _ (Ideal.mul_mem_mul hx (Ideal.absNorm_mem J))

/-- ★整数を生成元とする単項イデアルのノルムは `m ^ d`。 -/
theorem absNorm_span_natCast (m : ℕ) :
    absNorm (Ideal.span {(m : R)}) = m ^ Module.finrank ℤ R := by
  rw [Ideal.absNorm_span_singleton]
  have h : ((m : ℤ) : R) = (m : R) := by push_cast; ring
  rw [← h, show ((m : ℤ) : R) = algebraMap ℤ R (m : ℤ) from rfl,
    Algebra.norm_algebraMap (R := ℤ) (S := R) (m : ℤ)]
  simp

/-- ★★**分解** —— `N A = m·n`(互いに素)なら `A = (A + (m))·(A + (n))`。 -/
theorem mul_split_eq {A : Ideal R} {m n : ℕ} (hmn : Nat.Coprime m n)
    (hA : absNorm A = m * n) :
    (A ⊔ Ideal.span {(m : R)}) * (A ⊔ Ideal.span {(n : R)}) = A := by
  have htop : (A ⊔ Ideal.span {(m : R)}) ⊔ (A ⊔ Ideal.span {(n : R)}) = ⊤ := by
    refine top_le_iff.mp ?_
    rw [← Ideal.sup_eq_top_iff_isCoprime (m : R) (n : R) |>.mpr
      (isCoprime_cast_of_coprime hmn)]
    exact sup_le (le_sup_of_le_left le_sup_right) (le_sup_of_le_right le_sup_right)
  refine le_antisymm ?_ ?_
  · refine Ideal.mul_le.mpr fun r hr s hs => ?_
    obtain ⟨a, ha, z, hz, rfl⟩ := Submodule.mem_sup.mp hr
    obtain ⟨b, hb, w, hw, rfl⟩ := Submodule.mem_sup.mp hs
    obtain ⟨c, rfl⟩ := Ideal.mem_span_singleton'.mp hz
    obtain ⟨d, rfl⟩ := Ideal.mem_span_singleton'.mp hw
    have hmnA : ((m : R) * (n : R)) ∈ A := by
      have h := Ideal.absNorm_mem A
      rw [hA] at h
      simpa using h
    have hexp : (a + c * (m : R)) * (b + d * (n : R))
        = a * (b + d * (n : R)) + (c * b) * (m : R) + (c * d) * ((m : R) * (n : R)) := by ring
    rw [hexp]
    exact Ideal.add_mem _ (Ideal.add_mem _ (Ideal.mul_mem_right _ _ ha)
      (Ideal.mul_mem_right _ _ (Ideal.mul_mem_left _ _ hb))) (Ideal.mul_mem_left _ _ hmnA)
  · rw [Ideal.mul_eq_inf_of_isCoprime (Ideal.isCoprime_iff_sup_eq.mpr htop)]
    exact le_inf le_sup_left le_sup_left

/-- ★分解した因子のノルムは `m` を割る。 -/
theorem absNorm_sup_span_dvd {A : Ideal R} {m n : ℕ} (hmn : Nat.Coprime m n)
    (hA : absNorm A = m * n) :
    absNorm (A ⊔ Ideal.span {(m : R)}) ∣ m := by
  have h1 : absNorm (A ⊔ Ideal.span {(m : R)}) ∣ m ^ Module.finrank ℤ R := by
    have h := absNorm_dvd_absNorm_of_le (I := A ⊔ Ideal.span {(m : R)})
      (J := Ideal.span {(m : R)}) le_sup_right
    rwa [absNorm_span_natCast] at h
  have h2 : absNorm (A ⊔ Ideal.span {(m : R)}) ∣ m * n := by
    have h := absNorm_dvd_absNorm_of_le (I := A ⊔ Ideal.span {(m : R)}) (J := A) le_sup_left
    rwa [hA] at h
  exact (Nat.Coprime.coprime_dvd_left h1 (Nat.Coprime.pow_left _ hmn)).dvd_of_dvd_mul_right h2

/-- ★★★**ノルムはちょうど `m`**。 -/
theorem absNorm_sup_span_eq {A : Ideal R} {m n : ℕ} (hmn : Nat.Coprime m n)
    (hA : absNorm A = m * n) (hn : n ≠ 0) :
    absNorm (A ⊔ Ideal.span {(m : R)}) = m := by
  have hIm := absNorm_sup_span_dvd hmn hA
  have hJn := absNorm_sup_span_dvd hmn.symm (by rw [hA, Nat.mul_comm])
  have hprod :
      absNorm (A ⊔ Ideal.span {(m : R)}) * absNorm (A ⊔ Ideal.span {(n : R)}) = m * n := by
    rw [← map_mul, mul_split_eq hmn hA, hA]
  refine Nat.dvd_antisymm hIm ?_
  have hd : m * n ∣ absNorm (A ⊔ Ideal.span {(m : R)}) * n := by
    rw [← hprod]
    exact Nat.mul_dvd_mul_left _ hJn
  exact (Nat.mul_dvd_mul_iff_right (Nat.pos_of_ne_zero hn)).mp hd

/-- ★★★★**ノルム `m` × ノルム `n` ≃ ノルム `m·n`**(`m`, `n` は互いに素)。 -/
def coprimeIdealEquiv {m n : ℕ} (hmn : Nat.Coprime m n) (hm : m ≠ 0) (hn : n ≠ 0) :
    {I : Ideal R // absNorm I = m} × {J : Ideal R // absNorm J = n}
      ≃ {A : Ideal R // absNorm A = m * n} where
  toFun p := ⟨p.1.1 * p.2.1, by rw [map_mul, p.1.2, p.2.2]⟩
  invFun A := (⟨A.1 ⊔ Ideal.span {(m : R)}, absNorm_sup_span_eq hmn A.2 hn⟩,
    ⟨A.1 ⊔ Ideal.span {(n : R)},
      absNorm_sup_span_eq hmn.symm (by rw [A.2, Nat.mul_comm]) hm⟩)
  left_inv := by
    rintro ⟨⟨I, hI⟩, ⟨J, hJ⟩⟩
    have hco : Nat.Coprime (absNorm I) (absNorm J) := by rw [hI, hJ]; exact hmn
    refine Prod.ext (Subtype.ext ?_) (Subtype.ext ?_)
    · show I * J ⊔ Ideal.span {(m : R)} = I
      rw [← hI]
      exact mul_sup_span_absNorm hco
    · show I * J ⊔ Ideal.span {(n : R)} = J
      rw [← hJ, mul_comm I J]
      exact mul_sup_span_absNorm hco.symm
  right_inv := by
    rintro ⟨A, hA⟩
    exact Subtype.ext (mul_split_eq hmn hA)

/-- ★★★★★**イデアル計数関数は乗法的**。 -/
theorem card_absNorm_eq_mul {m n : ℕ} (hmn : Nat.Coprime m n) (hm : m ≠ 0) (hn : n ≠ 0) :
    Nat.card {A : Ideal R // absNorm A = m * n}
      = Nat.card {I : Ideal R // absNorm I = m} * Nat.card {J : Ideal R // absNorm J = n} := by
  rw [← Nat.card_prod]
  exact (Nat.card_congr (coprimeIdealEquiv hmn hm hn)).symm

end Dedekind

/-! ## ★2. イデアル計数関数 -/

/-- ★★**イデアル計数関数** `a_K(n) = #{I : N I = n}`。 -/
noncomputable def idealCount (K : Type*) [Field K] [NumberField K] (n : ℕ) : ℕ :=
  Nat.card {I : Ideal (𝓞 K) // Ideal.absNorm I = n}

variable (K : Type*) [Field K] [NumberField K]

@[simp] theorem idealCount_one : idealCount K 1 = 1 := by
  rw [idealCount, Nat.card_eq_one_iff_unique]
  refine ⟨⟨fun a b => Subtype.ext ?_⟩, ⟨⟨⊤, by simp⟩⟩⟩
  rw [absNorm_eq_one_iff.mp a.2, absNorm_eq_one_iff.mp b.2]

@[simp] theorem idealCount_zero : idealCount K 0 = 1 := by
  rw [idealCount, Nat.card_eq_one_iff_unique]
  refine ⟨⟨fun a b => Subtype.ext ?_⟩, ⟨⟨⊥, by simp⟩⟩⟩
  rw [absNorm_eq_zero_iff.mp a.2, absNorm_eq_zero_iff.mp b.2]

/-- ★★★★★**イデアル計数関数は乗法的**。 -/
theorem idealCount_mul {m n : ℕ} (hmn : Nat.Coprime m n) (hm : m ≠ 0) (hn : n ≠ 0) :
    idealCount K (m * n) = idealCount K m * idealCount K n :=
  card_absNorm_eq_mul hmn hm hn

theorem dedekindZeta_eq_lseries (s : ℂ) :
    dedekindZeta K s = LSeries (fun n ↦ (idealCount K n : ℂ)) s := rfl

/-! ## ★3. 部分和の増大度と絶対収束 -/

/-- ★`Icc 1 n` 上の計数の和は「ノルム `≤ n` の非零イデアルの個数」。 -/
theorem sum_card_absNorm_eq (n : ℕ) :
    Nat.card {I : (Ideal (𝓞 K))⁰ // absNorm (I : Ideal (𝓞 K)) ≤ n}
      = ∑ k ∈ Finset.Icc 1 n, Nat.card {I : Ideal (𝓞 K) // absNorm I = k} := by
  rw [← add_left_inj 1, ← card_norm_le_eq_card_norm_le_add_one,
    show Finset.Icc 1 n = Finset.Ioc 0 n from Finset.Icc_succ_left_eq_Ioc _ _,
    show 1 = Nat.card {I : Ideal (𝓞 K) // absNorm I = 0} by simp [Ideal.absNorm_eq_zero_iff],
    Finset.sum_Ioc_add_eq_sum_Icc (n.zero_le),
    ← Finset.card_preimage_eq_sum_card_image_eq (fun k _ ↦ finite_setOf_absNorm_eq k)]
  simp [Set.coe_eq_subtype]

/-- ★★部分和は `O(n)` —— mathlib の類数公式の漸近から。 -/
theorem sum_card_absNorm_isBigO :
    (fun n : ℕ ↦ ∑ k ∈ Finset.Icc 1 n, (idealCount K k : ℝ))
      =O[atTop] fun n : ℕ ↦ (n : ℝ) ^ (1 : ℝ) := by
  have hten : Tendsto (fun n : ℕ ↦
      (Nat.card {I : (Ideal (𝓞 K))⁰ // absNorm (I : Ideal (𝓞 K)) ≤ (n : ℝ)} : ℝ) / (n : ℝ))
      atTop (𝓝 _) :=
    (Ideal.tendsto_norm_le_div_atTop₀ K).comp tendsto_natCast_atTop_atTop
  have hmul := (hten.isBigO_one ℝ).mul (isBigO_refl (fun n : ℕ ↦ (n : ℝ)) atTop)
  refine hmul.congr (fun n ↦ ?_) (fun n ↦ ?_)
  · rcases Nat.eq_zero_or_pos n with rfl | hn
    · simp
    · have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn.ne'
      rw [div_mul_cancel₀ _ hn']
      have h : Nat.card {I : (Ideal (𝓞 K))⁰ // absNorm (I : Ideal (𝓞 K)) ≤ n}
          = ∑ k ∈ Finset.Icc 1 n, idealCount K k := sum_card_absNorm_eq K n
      simp only [Nat.cast_le] at *
      rw [h, Nat.cast_sum]
  · rw [one_mul, Real.rpow_one]

/-! ## ★4. Euler 積 -/

variable {s : ℂ}

/-- ★`n ↦ a_K(n)·n^{-s}`。 -/
noncomputable def zetaSummand (K : Type*) [Field K] [NumberField K] (s : ℂ) (n : ℕ) : ℂ :=
  (idealCount K n : ℂ) * (n : ℂ) ^ (-s)

theorem ne_zero_of_one_lt_re' (hs : 1 < s.re) : s ≠ 0 := by
  intro h; rw [h] at hs; norm_num at hs

theorem zetaSummand_eq_term (hs : s ≠ 0) (n : ℕ) :
    zetaSummand K s n = LSeries.term (fun n ↦ (idealCount K n : ℂ)) s n := by
  rcases eq_or_ne n 0 with rfl | hn
  · rw [LSeries.term_zero, zetaSummand,
      show ((0 : ℕ) : ℂ) = 0 from Nat.cast_zero, Complex.zero_cpow (neg_ne_zero.mpr hs),
      mul_zero]
  · rw [LSeries.term_of_ne_zero hn, zetaSummand, Complex.cpow_neg, div_eq_mul_inv]

theorem lSeriesSummable_idealCount (hs : 1 < s.re) :
    LSeriesSummable (fun n ↦ (idealCount K n : ℂ)) s := by
  have h := LSeriesSummable_of_sum_norm_bigO_and_nonneg (r := 1)
    (f := fun n : ℕ ↦ (idealCount K n : ℝ)) (sum_card_absNorm_isBigO K)
    (fun _ ↦ Nat.cast_nonneg _) zero_le_one hs
  simpa using h

theorem summable_zetaSummand (hs : 1 < s.re) :
    Summable (fun n ↦ ‖zetaSummand K s n‖) := by
  have h := (lSeriesSummable_idealCount K hs).norm
  exact h.congr fun n ↦ by rw [zetaSummand_eq_term K (ne_zero_of_one_lt_re' hs) n]

theorem tsum_zetaSummand (hs : 1 < s.re) :
    ∑' n : ℕ, zetaSummand K s n = dedekindZeta K s := by
  rw [dedekindZeta_eq_lseries, LSeries]
  exact tsum_congr fun n ↦ zetaSummand_eq_term K (ne_zero_of_one_lt_re' hs) n

@[simp] theorem zetaSummand_zero (hs : s ≠ 0) : zetaSummand K s 0 = 0 := by
  rw [zetaSummand, show ((0 : ℕ) : ℂ) = 0 from Nat.cast_zero,
    Complex.zero_cpow (neg_ne_zero.mpr hs), mul_zero]

@[simp] theorem zetaSummand_one : zetaSummand K s 1 = 1 := by
  rw [zetaSummand, idealCount_one]
  simp

/-- ★★★**加数は互いに素な引数で乗法的**。 -/
theorem zetaSummand_mul (hs : s ≠ 0) {m n : ℕ} (h : Nat.Coprime m n) :
    zetaSummand K s (m * n) = zetaSummand K s m * zetaSummand K s n := by
  rcases eq_or_ne m 0 with rfl | hm
  · rw [Nat.coprime_zero_left] at h
    subst h
    rw [zero_mul, zetaSummand_zero K hs, zero_mul]
  rcases eq_or_ne n 0 with rfl | hn
  · rw [Nat.coprime_zero_right] at h
    subst h
    rw [one_mul, zetaSummand_zero K hs, mul_zero]
  · rw [zetaSummand, zetaSummand, zetaSummand, idealCount_mul K h hm hn]
    have hcp : ((m * n : ℕ) : ℂ) ^ (-s) = (m : ℂ) ^ (-s) * (n : ℂ) ^ (-s) := by
      simpa only [Nat.cast_mul, ofReal_natCast]
        using Complex.mul_cpow_ofReal_nonneg m.cast_nonneg n.cast_nonneg (-s)
    rw [hcp]
    push_cast
    ring

/-- ★★★★★★**[cheb-euler] Dedekind ζ の Euler 積**(有理素数ごとの局所因子)。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★中身はイデアル計数関数の乗法性 1 本(`idealCount_mul`)と、
mathlib の `EulerProduct.eulerProduct_hasProd` である。 -/
theorem dedekindZeta_eulerProduct_hasProd (hs : 1 < s.re) :
    HasProd (fun p : Nat.Primes ↦ ∑' e : ℕ, zetaSummand K s (p ^ e)) (dedekindZeta K s) := by
  rw [← tsum_zetaSummand K hs]
  exact EulerProduct.eulerProduct_hasProd (zetaSummand_one K)
    (fun {_ _} h ↦ zetaSummand_mul K (ne_zero_of_one_lt_re' hs) h)
    (summable_zetaSummand K hs) (zetaSummand_zero K (ne_zero_of_one_lt_re' hs))

/-- ★★★★★**Euler 積**(`tprod` 版)。 -/
theorem dedekindZeta_eulerProduct_tprod (hs : 1 < s.re) :
    ∏' p : Nat.Primes, (∑' e : ℕ, zetaSummand K s (p ^ e)) = dedekindZeta K s :=
  (dedekindZeta_eulerProduct_hasProd K hs).tprod_eq

/-- ★★★★★**Euler 積**(有限部分積の収束の形)。 -/
theorem dedekindZeta_eulerProduct (hs : 1 < s.re) :
    Tendsto (fun n : ℕ ↦ ∏ p ∈ Nat.primesBelow n, ∑' e : ℕ, zetaSummand K s (p ^ e))
      atTop (𝓝 (dedekindZeta K s)) := by
  rw [← tsum_zetaSummand K hs]
  exact EulerProduct.eulerProduct (zetaSummand_one K)
    (fun {_ _} h ↦ zetaSummand_mul K (ne_zero_of_one_lt_re' hs) h)
    (summable_zetaSummand K hs) (zetaSummand_zero K (ne_zero_of_one_lt_re' hs))

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Theorem 6.4, (iv)` が使う Tchebotarev の一段目。 -/
def dedekindZeta_eulerProduct_hasProd.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Dedekind ζ の Euler 積(密度論法の一段目)",
    sectionId := "frdi-thm-6-4" }

/-- ★★★★★locator —— イデアル計数関数の乗法性(mathlib に無いので証明した)。 -/
def idealCount_mul.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — イデアル計数関数の乗法性",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
