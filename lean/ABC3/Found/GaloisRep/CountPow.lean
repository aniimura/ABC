import ABC3.Found.GaloisRep.PullbackAlg
import Mathlib.RingTheory.DedekindDomain.Factorization

/-!
# Galois (G5) 第 142 ブロック —— **★★★★★★指数が `n` で割れれば `n` 乗**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★D1 の最終組み立てが取れた

D1 は `∃ J, J^n = (μ f_P)` である。★本ブロックはその**最後の一歩**を出す:

    ∀ v, n ∣ count_v(I)     ⟹     ∃ J, J^n = I

★★これで D1 に残るのは「**各素点での指数が `n` で割れること**」だけになる。

## ★★★★mathlib の在庫が効いた(2026-08-20 実測)

`FractionalIdeal.count : FractionalIdeal R⁰ K → ℤ` は**整数値**であり、
`count_pow`・`count_finprod`・`finite_factors`・
`finprod_heightOneSpectrum_factorization'` が揃っている。
★付値(`WithZero (Multiplicative ℤ)`、乗法的)を経由せずに**加法的に扱える**。

★★証明は `exps v := count_v(I)/n` と置いて `J := ∏ᶠ v, v^(exps v)` を取るだけ。
`finprod_pow` で `J^n = ∏ᶠ v, v^(exps v · n) = ∏ᶠ v, v^(count_v I) = I`。

## ★★`Ψ₂Sq` の次数は 3(奇数)

D1 の場合 B(`[n]Q = O`、無限遠の上)で使う。
★`z² = Ψ₂Sq(x)` と `deg Ψ₂Sq = 3` から `2ν(z) = 3ν(x)` となり、
`gcd(2,3) = 1` より **`ν(x)` は偶数**である。これが「無限遠で分岐する」ことの初等的な形。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_pow_of_dvd_count` | ★★★★★★**指数が `n` で割れれば `n` 乗になる** |
| `fP_ne_zero` | ★`f_P ≠ 0` |
| `spanSingleton_mu_ne_zero` | ★★`(μ f_P) ≠ 0` |
| `degree_Psi2Sq` / `natDegree_Psi2Sq` | ★★`deg Ψ₂Sq = 3` |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain nonZeroDivisors

/-! ## ★★★★★★指数が `n` で割れれば `n` 乗 -/

/-- ★★★★★★**各素点での指数が `n` で割れる分数イデアルは `n` 乗である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`exps v := count_v(I)/n` と置き、`J := ∏ᶠ v, v^(exps v)` を取る。
★★`finprod_pow` と `finprod_heightOneSpectrum_factorization'`(どちらも mathlib)で閉じる。 -/
theorem exists_pow_of_dvd_count {R K : Type} [CommRing R] [IsDedekindDomain R]
    [Field K] [Algebra R K] [IsFractionRing R K]
    {I : FractionalIdeal R⁰ K} (hI : I ≠ 0) {n : ℕ}
    (h : ∀ v : HeightOneSpectrum R, (n : ℤ) ∣ FractionalIdeal.count K v I) :
    ∃ J : FractionalIdeal R⁰ K, J ^ n = I := by
  classical
  set exps : HeightOneSpectrum R → ℤ := fun v => FractionalIdeal.count K v I / n with hexps
  have hcof : ∀ᶠ v : HeightOneSpectrum R in Filter.cofinite, exps v = 0 := by
    filter_upwards [FractionalIdeal.finite_factors (K := K) I] with v hv
    simp [hexps, hv]
  have hfin : Function.HasFiniteMulSupport
      (fun v : HeightOneSpectrum R => (v.asIdeal : FractionalIdeal R⁰ K) ^ exps v) := by
    refine Set.Finite.subset (Filter.eventually_cofinite.mp hcof) ?_
    intro v hv
    simp only [Function.mem_mulSupport] at hv
    simp only [Set.mem_setOf_eq]
    intro h0
    exact hv (by rw [h0, zpow_zero])
  refine ⟨∏ᶠ v : HeightOneSpectrum R, (v.asIdeal : FractionalIdeal R⁰ K) ^ exps v, ?_⟩
  rw [finprod_pow hfin, ← FractionalIdeal.finprod_heightOneSpectrum_factorization' K hI]
  refine finprod_congr fun v => ?_
  rw [← zpow_natCast ((v.asIdeal : FractionalIdeal R⁰ K) ^ exps v) n, ← zpow_mul]
  congr 1
  exact Int.ediv_mul_cancel (h v)

/-! ## ★`f_P ≠ 0` と `(μ f_P) ≠ 0` -/

variable {F : Type} [Field F]

/-- ★`f_P ≠ 0`——`XClass^n` が `(f_P)` に入ることから。 -/
theorem fP_ne_zero (W : WeierstrassCurve.Affine F) {x y : F} (n : ℕ) (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP}) : fP ≠ 0 := by
  intro h0
  have hmem : (CoordinateRing.XClass W x) ^ n ∈ Ideal.span ({fP} : Set W.CoordinateRing) := by
    rw [← hfP]; exact Ideal.pow_mem_pow (Ideal.subset_span (by simp)) n
  rw [h0, Ideal.span_singleton_eq_bot.2 rfl, Ideal.mem_bot] at hmem
  exact pow_ne_zero n (CoordinateRing.XClass_ne_zero x) hmem

/-- ★★`(μ f_P) ≠ 0`——`μ` が単射だから。 -/
theorem spanSingleton_mu_ne_zero (W : WeierstrassCurve.Affine F) {x y : F} (n : ℕ)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ) :
    FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) ≠ 0 := by
  rw [Ne, FractionalIdeal.spanSingleton_eq_zero_iff]
  intro hz
  exact fP_ne_zero W n fP hfP (hμinj (by rw [hz, map_zero]))

/-! ## ★★`Ψ₂Sq` の次数は 3 -/

/-- ★★**`deg Ψ₂Sq = 3`**(標数 ≠ 2)——`twoTorsionPolynomial` の主係数が `4` だから。 -/
theorem degree_Psi2Sq (W : WeierstrassCurve F) (h2 : IsUnit (2 : F)) :
    W.Ψ₂Sq.degree = 3 := by
  have h4 : (W.twoTorsionPolynomial).a ≠ 0 := by
    show (4 : F) ≠ 0
    have h44 : (4 : F) = 2 * 2 := by norm_num
    rw [h44]; exact mul_ne_zero h2.ne_zero h2.ne_zero
  rw [Psi2Sq_eq_twoTorsionPolynomial]
  exact Cubic.degree_of_a_ne_zero h4

/-- ★★`natDegree Ψ₂Sq = 3`。 -/
theorem natDegree_Psi2Sq (W : WeierstrassCurve F) (h2 : IsUnit (2 : F)) :
    W.Ψ₂Sq.natDegree = 3 := by
  have h4 : (W.twoTorsionPolynomial).a ≠ 0 := by
    show (4 : F) ≠ 0
    have h44 : (4 : F) = 2 * 2 := by norm_num
    rw [h44]; exact mul_ne_zero h2.ne_zero h2.ne_zero
  rw [Psi2Sq_eq_twoTorsionPolynomial]
  exact Cubic.natDegree_of_a_ne_zero h4

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_of_dvd_count.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——指数が n で割れる分数イデアルが n 乗であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
