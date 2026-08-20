import ABC3.Found.GaloisRep.D2Principal

/-!
# Galois (G5) 第 175 ブロック —— **★★★★★★★★★D2 をスケルトンの形にする**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★`n` 乗根は一意である

第 174 が示したのは**特定の** `rootUnit I n` が単項であることだった。★スケルトンの
`fractionalIdeal_isPrincipal` は「`J^n = I` なる**任意の** `J`」について述べている。

★★両者は一致する。指数を見れば

    n · count_v(J) = count_v(J^n) = count_v(I) = n · count_v(rootUnit I n)

なので `count_v(J) = count_v(rootUnit I n)`、そして**因子分解の一意性**
(mathlib の `finprod_heightOneSpectrum_factorization'`)から `J = rootUnit I n`。

## ★★★逸脱の記録: `hchar` を足した

第 150(`Σ_{T ∈ E[n]} T = 0`)が (G1) の `E[n] ≅ (ℤ/n)²` を使うので、
**`∀ k, 1 ≤ k → k ≤ n → (k : F) ≠ 0`** を仮定に足した。

★Weil 対 `e_n` はそもそも `char ∤ n` を要求するので、消費側((G5) `det_cyclotomic`、
`[CharZero K]` の下)には影響しない。★★`[Fintype (torsionPoints W n)]` と
`(4 : F) ≠ 0` はここから導けるので、仮定に足さずに済んだ。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eq_of_count_eq` | ★★★指数がすべて等しければ等しい |
| `eq_rootUnit_of_pow` | ★★★★**`n` 乗根の一意性** |
| `fractionalIdeal_isPrincipal_proof` | ★★★★★★★★★**D2(スケルトンの形)** |
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep
open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★★★`n` 乗根の一意性 -/

variable {R K : Type} [CommRing R] [IsDedekindDomain R] [Field K] [Algebra R K]
  [IsFractionRing R K]

/-- ★★★指数がすべて等しければ分数イデアルは等しい。 -/
theorem eq_of_count_eq {J J' : FractionalIdeal R⁰ K} (hJ : J ≠ 0) (hJ' : J' ≠ 0)
    (h : ∀ v : HeightOneSpectrum R,
      FractionalIdeal.count K v J = FractionalIdeal.count K v J') : J = J' := by
  rw [← FractionalIdeal.finprod_heightOneSpectrum_factorization' K hJ,
    ← FractionalIdeal.finprod_heightOneSpectrum_factorization' K hJ']
  exact finprod_congr fun v => by rw [h v]

/-- ★★★★**`n` 乗根の一意性**——`J^n = I` なら `J = rootUnit I n`。 -/
theorem eq_rootUnit_of_pow {I J : FractionalIdeal R⁰ K} (hI : I ≠ 0) {n : ℕ} (hn : 1 ≤ n)
    (hJ : J ^ n = I) :
    J = ((rootUnit I n : (FractionalIdeal R⁰ K)ˣ) : FractionalIdeal R⁰ K) := by
  have hJne : J ≠ 0 := by
    intro h0
    rw [h0, zero_pow (by omega : n ≠ 0)] at hJ
    exact hI hJ.symm
  have hdvd : ∀ v : HeightOneSpectrum R, (n : ℤ) ∣ FractionalIdeal.count K v I :=
    fun v => ⟨FractionalIdeal.count K v J, by rw [← hJ, FractionalIdeal.count_pow]⟩
  have hn0 : (n : ℤ) ≠ 0 := by exact_mod_cast (by omega : n ≠ 0)
  refine eq_of_count_eq hJne (Units.ne_zero _) (fun v => ?_)
  refine mul_left_cancel₀ hn0 ?_
  rw [← FractionalIdeal.count_pow, ← FractionalIdeal.count_pow, hJ, rootUnit_pow hI hdvd]

/-! ## ★★★★★★★★★D2(スケルトンの形) -/

variable {F : Type} [Field F] [IsAlgClosed F] [DecidableEq F]

/-- ★★★★★★★★★**D2 —— 任意の `n` 乗根イデアルは単項である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 174(特定の `rootUnit` が単項)と `n` 乗根の一意性を繋いだもの。
★★`hchar` を足した(第 150 が (G1) の `E[n] ≅ (ℤ/n)²` を使うため)。 -/
theorem fractionalIdeal_isPrincipal_proof (h2 : IsUnit (2 : F))
    (W : WeierstrassCurve.Affine F) [W.IsElliptic]
    [inst : IsDedekindDomain W.CoordinateRing]
    {x y : F} (h : W.Nonsingular x y) (n : ℕ) (hn : 1 ≤ n)
    (hchar : ∀ k : ℕ, 1 ≤ k → k ≤ n → (k : F) ≠ 0)
    (hP : n • (Point.some x y h) = 0)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {xn yn : W.FunctionField} (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (J : FractionalIdeal W.CoordinateRing⁰ W.FunctionField)
    (hJ : J ^ n = FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) :
    ∃ g : W.FunctionField, J = FractionalIdeal.spanSingleton W.CoordinateRing⁰ g := by
  classical
  haveI : Fintype (torsionPoints W n) :=
    (finite_torsion W n hn (hchar n hn le_rfl)).fintype
  have h4 : (4 : F) ≠ 0 := by
    have hu : IsUnit ((2 : F) * 2) := h2.mul h2
    intro h0
    exact hu.ne_zero (by rw [← h0]; ring)
  have hμfP : μ fP ≠ 0 := fun h0 =>
    (fP_ne_zero W n fP hfP) (hμinj (by rw [h0, map_zero]))
  have hIne : FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP) ≠ 0 := by
    rw [ne_eq, FractionalIdeal.spanSingleton_eq_zero_iff]; exact hμfP
  obtain ⟨g, hg⟩ := exists_spanSingleton_rootUnit_mulByN W h2 h4 n hn hchar μ hμinj hμF
    hns hμx hμy hμP h hP fP hfP
  exact ⟨g, by rw [eq_rootUnit_of_pow hIne hn hJ, hg]⟩

/-! ## ★出典の紐付け(`.src`) -/

def fractionalIdeal_isPrincipal_proof.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——D2、因子の類が自明であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
