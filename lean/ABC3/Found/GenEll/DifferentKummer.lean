import ABC3.Meta.Claim
import Mathlib.RingTheory.DedekindDomain.Different

/-!
# [GenEll] Proposition 1.7 の "elementary claim" —— **Kummer 拡大の different**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★原文が畳んでいる "elementary claim"

原文 p.10 は `Proposition 1.7, (i)` の証明を、最後に**局所体の主張 1 つ**へ落としている:

> Fix a prime number `p` and a positive integer `d`. Then there exists a positive
> integer `n` such that for any finite Galois extension `L/K` of finite extensions
> of `ℚ_p` with `[L : K] ≤ d`, the different ideal of `L/K` contains `p^n · O_L`.

★原文の証明は 6 段である:

1. `L[ζ]/K` を野生分岐と馴分岐に分ける(馴分岐なら `n = 1` で足りる)
2. `ζ_p ∈ K` の野生分岐の場合に帰着
3. `Gal(L/K)` は位数 `≤ d` の(**可解な**)`p`-群なので `[L:K] = p` の場合に帰着
4. Kummer 理論で `L = K(λ)`、`κ ≝ λ^p ∈ K`
5. `κ` を `(K^×)^p` の元倍して `κ = 単元 × π_K^a`(`0 ≤ a < p`)に正規化。
   `κ ∈ O_K` かつ `κ ∉ p^p·O_K` なので **`O_L ⊇ λ·O_L ⊇ p·O_L`**
6. `O_K[X]/(X^p − κ) ↪ O_L` から **different は `p·λ^{p−1}·O_L ⊇ p^p·O_L` を含む**

## ★★本ファイルが取るのは 5 と 6 である

★6 の中身は mathlib の `aeval_derivative_mem_differentIdeal`
(`f'(λ) ∈ different`)を `f = X^n − κ` に当てるだけである
——`f' = n·X^{n−1}` なので `n·λ^{n−1} ∈ different` が出る。

★★最小多項式を `A = O_K` の上で読むところに 1 段いる:
`minpoly K` の側で `X^n − κ` と分かっていても、`different` は `minpoly A` を使うので、
**`A` が整閉であること**(`minpoly.isIntegrallyClosed_eq_field_fractions'`)で降ろす。

## ★残り

| 段 | 状態 |
|---|---|
| 1(野生/馴の分離) | 未着手 |
| 2(`ζ_p ∈ K` への帰着) | 未着手 |
| 3(`p`-群の可解性で `[L:K] = p` へ) | 未着手 |
| 4(Kummer 理論) | 未着手 |
| **5(`p ∈ λ·O_L`)** | ★**本ファイル**——付値環の全順序性 1 行 |
| **6(different の下界)** | ★**本ファイル** |
| **5+6 を繋いだ結論 `p^p·O_L ⊆ different`** | ★★**本ファイル** |
-/

namespace ABC3.Found.GenEll

open Polynomial

/-! ## ★最小多項式を整閉な底へ降ろす -/

/-- ★★**`minpoly K` が `X^n − κ` なら `minpoly A` も `X^n − κ`**(`A` は整閉、`K = Frac A`)。

★`different` は `minpoly A` を使うので、この降ろしが要る。 -/
theorem minpoly_eq_X_pow_sub_C_of_map
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [IsIntegrallyClosed A] [IsDomain B] [IsFractionRing B L]
    (n : ℕ) (lam : B) (kappa : A) (hint : IsIntegral A lam)
    (hmapK : minpoly K (algebraMap B L lam) = X ^ n - C (algebraMap A K kappa)) :
    minpoly A lam = X ^ n - C kappa := by
  have hintL : IsIntegral A (algebraMap B L lam) := hint.map (IsScalarTower.toAlgHom A B L)
  have h1 := minpoly.isIntegrallyClosed_eq_field_fractions' (R := A) (S := L) K hintL
  have h2 : minpoly A (algebraMap B L lam) = minpoly A lam :=
    minpoly.algebraMap_eq (IsFractionRing.injective B L) lam
  rw [h2] at h1
  rw [hmapK] at h1
  refine Polynomial.map_injective (algebraMap A K) (IsFractionRing.injective A K) ?_
  rw [← h1]
  simp [Polynomial.map_sub, Polynomial.map_pow]

/-! ## ★★★★different の下界 -/

/-- ★★★★**原文 p.10 の最終計算** —— `minpoly_A(λ) = X^n − κ` なら
`n·λ^{n−1} ∈ different(L/K)`。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★原文は『one computes easily that the different ideal of `L/K` contains
`p · λ^{p−1} · O_L`』と 1 文で済ませている。
★★mathlib の `aeval_derivative_mem_differentIdeal` に `f = X^n − κ` を当てるだけである
——`f' = n·X^{n−1}`。 -/
theorem natCast_mul_pow_mem_differentIdeal
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    (n : ℕ) (lam : B) (kappa : A)
    (hmin : minpoly A lam = X ^ n - C kappa)
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤) :
    (n : B) * lam ^ (n - 1) ∈ differentIdeal A B := by
  have h := aeval_derivative_mem_differentIdeal A K L lam hgen
  rw [hmin] at h
  simpa [derivative_sub, derivative_X_pow, derivative_C] using h

/-- ★★★★★**上の 2 つを繋いだ形** —— `minpoly K` の側だけ分かっていればよい。 -/
theorem natCast_mul_pow_mem_differentIdeal_of_map
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    [IsFractionRing B L]
    (n : ℕ) (lam : B) (kappa : A) (hint : IsIntegral A lam)
    (hmapK : minpoly K (algebraMap B L lam) = X ^ n - C (algebraMap A K kappa))
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤) :
    (n : B) * lam ^ (n - 1) ∈ differentIdeal A B :=
  natCast_mul_pow_mem_differentIdeal A K L n lam kappa
    (minpoly_eq_X_pow_sub_C_of_map A K L n lam kappa hint hmapK) hgen

/-! ## ★★★★★段 5 —— `λ` は `p` を割る -/

/-- ★★★**付値環では `¬(p^n ∣ λ^n)` から `λ ∣ p` が出る**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★原文の『`κ ∈ O_K` but `κ ∉ p^p·O_K`, hence `O_L ⊇ λ·O_L ⊇ p·O_L`』の段である。
★★付値環なので**割り切れどうしが全順序**であることしか使わない。 -/
theorem dvd_of_not_pow_dvd_pow (B : Type*) [CommRing B] [IsDomain B] [ValuationRing B]
    (n : ℕ) (lam pp : B) (h : ¬ (pp ^ n ∣ lam ^ n)) : lam ∣ pp := by
  rcases ValuationRing.cond lam pp with ⟨c, hc | hc⟩
  · exact ⟨c, hc.symm⟩
  · exact absurd (hc ▸ (pow_dvd_pow_of_dvd ⟨c, rfl⟩ n)) h

/-- ★★**`λ ∣ p` なら `p·λ^{n−1} ∣ p^n`**——原文の `p·λ^{p−1}·O_L ⊇ p^{1+(p−1)}·O_L`。 -/
theorem mul_pow_dvd_pow (B : Type*) [CommRing B] (n : ℕ) (hn : 0 < n) (lam pp : B)
    (h : lam ∣ pp) : pp * lam ^ (n - 1) ∣ pp ^ n := by
  have h1 : lam ^ (n - 1) ∣ pp ^ (n - 1) := pow_dvd_pow_of_dvd h _
  have h2 : pp ^ n = pp * pp ^ (n - 1) := by
    conv_lhs => rw [show n = 1 + (n - 1) by omega]
    rw [pow_add, pow_one]
  rw [h2]
  exact mul_dvd_mul_left pp h1

/-- ★★★★★★**原文 p.10 の結論** —— `p^p·O_L ⊆ different(L/K)`。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★段 5(`λ ∣ p`)と段 6(`p·λ^{p−1} ∈ different`)を繋いだものである。
★★`hcoef` は「係数の `n` が剰余体の標数 `p` に一致する」ことである。 -/
theorem pow_mem_differentIdeal
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    (n : ℕ) (hn : 0 < n) (lam : B) (kappa : A) (pp : B)
    (hdvd : lam ∣ pp) (hcoef : (n : B) = pp)
    (hmin : minpoly A lam = X ^ n - C kappa)
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤) :
    pp ^ n ∈ differentIdeal A B := by
  have h1 := natCast_mul_pow_mem_differentIdeal A K L n lam kappa hmin hgen
  rw [hcoef] at h1
  obtain ⟨c, hc⟩ := mul_pow_dvd_pow B n hn lam pp hdvd
  rw [hc]
  exact Ideal.mul_mem_right c _ h1

/-! ## ★出典の紐付け(`.src`) -/

def minpoly_eq_X_pow_sub_C_of_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(最小多項式を O_K の上で読む段)",
    sectionId := "genell-prop-1-7" }

def natCast_mul_pow_mem_differentIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(different は p·λ^{p−1}·O_L を含む)",
    sectionId := "genell-prop-1-7" }

def natCast_mul_pow_mem_differentIdeal_of_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(minpoly_K の側から different の下界へ)",
    sectionId := "genell-prop-1-7" }

def dvd_of_not_pow_dvd_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(κ ∉ p^p·O_K から λ·O_L ⊇ p·O_L)",
    sectionId := "genell-prop-1-7" }

def pow_mem_differentIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(p^p·O_L ⊆ different(L/K))",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
