import ABC3.Meta.Claim
import Mathlib.RingTheory.LocalRing.ResidueField.Defs
import Mathlib.RingTheory.LocalRing.Basic
import Mathlib.Algebra.CharP.Basic
import Mathlib.RingTheory.Polynomial.Eisenstein.Basic
import Mathlib.RingTheory.PrincipalIdealDomain
import Mathlib.RingTheory.Polynomial.Eisenstein.IsIntegral

/-!
# [GenEll] Proposition 1.7 の "elementary claim" —— **馴分岐**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★★★mathlib に馴分岐は無い(2026-08-26 実測)

REPL で名前を引いた結果:

| 探したもの | mathlib |
|---|---|
| `Algebra.IsUnramifiedAt` | ✅ ある |
| `IsTamelyRamified` / `Algebra.IsTamelyRamified` / `Ideal.IsTamelyRamified` / `IsTame` | ❌ **どれも無い** |

★したがって**建てる**しかない。本ファイルはその最初の一段である。

## ★★★何を取るか —— 第 379 の仮説は馴分岐そのものである

`DifferentKummer.lean` の `mem_differentIdeal_of_isUnit_natCast`(第 379)は
**`IsUnit (n : B)`** を仮定していた。★これが**馴分岐の言い換え**であることを示す:

> 局所環 `B` で `IsUnit (n : B)` ⟺ **剰余標数が `n` を割らない**

★★これで「野生か馴か」が **`n` が単元か**という 1 つの述語に集約される。
★★★原文の段 1(野生と馴に分ける)は、形式化の側では**この 1 行の場合分け**になる。
-/

namespace ABC3.Found.GenEll

open IsLocalRing

/-! ## ★★★★馴分岐の定義 -/

/-- ★★★★**馴分岐の次数**(局所環の言葉で)—— 剰余標数が `n` を割らない。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★原文が『tamely ramified』と呼ぶものを、**次数の側だけ**取り出した形である
(剰余体の分離性は `Algebra.IsSeparable` として別に持てる)。 -/
def IsTameDegree (B : Type*) [CommRing B] [IsLocalRing B] (n : ℕ) : Prop :=
  ¬ (ringChar (ResidueField B) ∣ n)

/-! ## ★★★★★`IsUnit (n : B)` は馴分岐そのもの -/

/-- ★★★★★**局所環では `IsUnit (n : B)` ⟺ 剰余標数が `n` を割らない**。

★機構は 3 行——単元であることは極大イデアルに入らないこと、
極大イデアルに入ることは剰余体で `0` になること、
剰余体で `0` になることは標数が割ることである。 -/
theorem isUnit_natCast_iff (B : Type*) [CommRing B] [IsLocalRing B] (ell n : ℕ)
    (hchar : ringChar (ResidueField B) = ell) :
    IsUnit (n : B) ↔ ¬ (ell ∣ n) := by
  have hres : ∀ x : B, residue B x = 0 ↔ x ∈ maximalIdeal B := fun x =>
    ⟨fun h => (Ideal.Quotient.eq_zero_iff_mem (I := maximalIdeal B) (a := x)).mp h,
     fun h => (Ideal.Quotient.eq_zero_iff_mem (I := maximalIdeal B) (a := x)).mpr h⟩
  have hcast : ((n : ℕ) : ResidueField B) = residue B (n : B) := by simp
  rw [← notMem_maximalIdeal, ← hres, ← hcast, ringChar.spec, hchar]

/-- ★★★★★★**第 379 の仮説は馴分岐である**(言い換え)。

★これで `DifferentKummer.lean` の `mem_differentIdeal_of_isUnit_natCast` が
「**馴分岐なら `p·O_L ⊆ 𝔡`**」を言っていることが型で見える。 -/
theorem isUnit_natCast_iff_isTameDegree (B : Type*) [CommRing B] [IsLocalRing B] (n : ℕ) :
    IsUnit (n : B) ↔ IsTameDegree B n :=
  isUnit_natCast_iff B _ n rfl

/-- ★不分岐(次数 `1`)は馴分岐である。 -/
theorem isTameDegree_one (B : Type*) [CommRing B] [IsLocalRing B]
    (hchar : ringChar (ResidueField B) ≠ 1) : IsTameDegree B 1 := by
  intro h
  exact hchar (Nat.dvd_one.mp h)

/-! ## ★★★★★★Eisenstein の導関数 —— 馴なら指数はちょうど `e−1` -/

/-- ★単元に極大イデアルの元を足しても単元である。 -/
theorem isUnit_add_mem_maximalIdeal {B : Type*} [CommRing B] [IsLocalRing B] {x y : B}
    (hx : IsUnit x) (hy : y ∈ maximalIdeal B) : IsUnit (x + y) := by
  rw [← notMem_maximalIdeal]
  intro hmem
  rw [← notMem_maximalIdeal] at hx
  exact hx (by simpa using Ideal.sub_mem _ hmem hy)

open Polynomial Finset in
/-- ★★★★★**Eisenstein 多項式の導関数の値**——
`f′(λ) = e·λ^{e−1} + λ^e·(…)`。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★Eisenstein の仮説は「係数が `π` で割れる」であり、全分岐なら `λ^e ∣ π` なので
**`λ^e ∣ a_i`** として受けている。 -/
theorem aeval_derivative_eisenstein {B : Type*} [CommRing B] (e : ℕ) (a : ℕ → B) (lam : B)
    (hdvd : ∀ i ∈ range e, lam ^ e ∣ a i) :
    ∃ c : B, (aeval lam) (derivative (X ^ e + ∑ i ∈ range e, C (a i) * X ^ i))
      = (e : B) * lam ^ (e - 1) + lam ^ e * c := by
  have hder : derivative (X ^ e + ∑ i ∈ range e, C (a i) * X ^ i)
      = C (e : B) * X ^ (e - 1) + ∑ i ∈ range e, C (a i) * (C (i : B) * X ^ (i - 1)) := by
    rw [derivative_add, derivative_X_pow, derivative_sum]
    congr 1
    refine Finset.sum_congr rfl (fun i _ => ?_)
    rw [derivative_C_mul, derivative_X_pow]
  have hsum : lam ^ e ∣ ∑ i ∈ range e, a i * ((i : B) * lam ^ (i - 1)) := by
    refine Finset.dvd_sum (fun i hi => ?_)
    exact Dvd.dvd.mul_right (hdvd i hi) _
  obtain ⟨c, hc⟩ := hsum
  refine ⟨c, ?_⟩
  rw [hder]
  simp only [map_add, map_mul, map_sum, aeval_C, aeval_X, map_pow]
  rw [← hc]
  simp

open Polynomial Finset in
/-- ★★★★★★★**馴分岐なら `f′(λ)` は `λ^{e−1}` の単元倍**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★単生成なら `𝔡 = (f′(λ))` なので(`DifferentKummer.lean` の
`differentIdeal_eq_span_of_adjoin_eq_top`、第 378)、これは
**different の指数がちょうど `e−1`** であることを意味する——
mathlib の `pow_sub_one_dvd_differentIdeal` は**下界しか**与えないので、
こちらが**上界**にあたる。

★★機構: `f′(λ) = λ^{e−1}·(e + λ·c)` で、馴(`e` が単元)なら
`e + λ·c` も単元である(`λ` は極大イデアルの元だから)。 -/
theorem aeval_derivative_eisenstein_tame {B : Type*} [CommRing B] [IsLocalRing B]
    (e : ℕ) (he : 0 < e) (a : ℕ → B) (lam : B)
    (hdvd : ∀ i ∈ range e, lam ^ e ∣ a i)
    (hlam : lam ∈ maximalIdeal B) (hunit : IsUnit (e : B)) :
    ∃ v : B, IsUnit v ∧
      (aeval lam) (derivative (X ^ e + ∑ i ∈ range e, C (a i) * X ^ i)) = lam ^ (e - 1) * v := by
  obtain ⟨c, hc⟩ := aeval_derivative_eisenstein e a lam hdvd
  refine ⟨(e : B) + lam * c, isUnit_add_mem_maximalIdeal hunit (Ideal.mul_mem_right c _ hlam), ?_⟩
  rw [hc]
  have hsplit : lam ^ e = lam ^ (e - 1) * lam := by
    conv_lhs => rw [show e = (e - 1) + 1 by omega]
    rw [pow_succ]
  rw [hsplit]
  ring

/-! ## ★★★★全分岐のとき Eisenstein の条件は `λ^e ∣ ·` になる -/

/-- ★★**`λ^e` と `π` が単元倍の違いなら、割り切れどうしが一致する**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★**全分岐**なら `π_K` は `B` で `λ^e` の単元倍になる。
★★したがって Eisenstein の標準的な条件『係数が `π_K` で割れる』は、
`DifferentKummer.lean` の `mem_differentIdeal_of_eisenstein_tame` が要求する
**`λ^e ∣ ·`** にそのままなる。 -/
theorem pow_dvd_iff_of_isUnit_mul {B : Type*} [CommRing B] {e : ℕ} {lam pi u x : B}
    (hu : IsUnit u) (hpi : lam ^ e = pi * u) :
    lam ^ e ∣ x ↔ pi ∣ x := by
  constructor
  · intro h
    exact dvd_trans (Dvd.intro u hpi.symm) h
  · intro h
    rw [hpi, mul_comm]
    exact (hu.mul_left_dvd).mpr h

/-- ★★★**Eisenstein の条件を `λ^e ∣ ·` へ移す**——全分岐の場合。 -/
theorem pow_dvd_of_dvd_of_isUnit_mul {B : Type*} [CommRing B] {e : ℕ} {lam pi u : B}
    (hu : IsUnit u) (hpi : lam ^ e = pi * u) {x : B} (hx : pi ∣ x) :
    lam ^ e ∣ x :=
  (pow_dvd_iff_of_isUnit_mul hu hpi).mpr hx

/-! ## ★★★★★全分岐 ⇒ `π` は `λ^e` の単元倍 -/

/-- ★★**イデアルが一致すれば単元倍である**。 -/
theorem associated_pow_of_span_eq {B : Type*} [CommRing B] [IsDomain B] {lam x : B} {e : ℕ}
    (h : Ideal.span {x} = Ideal.span {lam} ^ e) : Associated x (lam ^ e) := by
  rw [Ideal.span_singleton_pow] at h
  exact (Ideal.span_singleton_eq_span_singleton).mp h

/-- ★★★★**全分岐なら `π` は `λ^e` の単元倍**——
`pow_dvd_iff_of_isUnit_mul` が要求する形で出す。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★仮説の `Ideal.span {π} = Ideal.span {λ}^e` が**全分岐の定義そのもの**である
——`π·B = 𝔪_B^e`。 -/
theorem exists_isUnit_pow_eq_of_span_eq {B : Type*} [CommRing B] [IsDomain B] {lam x : B} {e : ℕ}
    (h : Ideal.span {x} = Ideal.span {lam} ^ e) : ∃ u : B, IsUnit u ∧ lam ^ e = x * u := by
  obtain ⟨u, hu⟩ := associated_pow_of_span_eq h
  exact ⟨(u : B), u.isUnit, hu.symm⟩

/-! ## ★★★★★★Eisenstein なら `p^k` を割って降りられる -/

open Polynomial in
/-- ★★★★★★**Eisenstein なら `p^k·z ∈ A[λ]` から `z ∈ A[λ]` が出る**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★これが構造定理の 3 つ目(**`B = A[λ]`**)への降下の段である——
任意の `z ∈ O_L` に対して `p^k·z ∈ O_K[λ]` となる `k` が取れれば、
この補題で `k` を `0` まで落とせる。

★★**1 段分は mathlib にあった**——
`mem_adjoin_of_smul_prime_smul_of_minpoly_isEisensteinAt`。
本定理はそれを `k` 回繰り返す帰納である。 -/
theorem mem_adjoin_of_pow_smul_of_isEisensteinAt {R : Type*} {K : Type*} {L : Type*} {pi : R}
    [CommRing R] [Field K] [Field L] [Algebra K L] [Algebra R L] [Algebra R K]
    [IsScalarTower R K L] [IsDomain R] [IsFractionRing R K] [IsIntegrallyClosed R]
    {PB : PowerBasis K L} (hp : Prime pi) (hBint : IsIntegral R PB.gen)
    (hEis : (minpoly R PB.gen).IsEisensteinAt (Ideal.span {pi})) :
    ∀ (k : ℕ) (z : L), IsIntegral R z → (pi ^ k) • z ∈ Algebra.adjoin R {PB.gen} →
      z ∈ Algebra.adjoin R {PB.gen} := by
  intro k
  induction k with
  | zero => intro z _ h; simpa using h
  | succ n ih =>
      intro z hz h
      have hpz : IsIntegral R (pi • z) := by
        rw [Algebra.smul_def]
        exact (isIntegral_algebraMap).mul hz
      have h' : (pi ^ n) • (pi • z) ∈ Algebra.adjoin R {PB.gen} := by
        rw [smul_smul]
        convert h using 2
        ring
      exact mem_adjoin_of_smul_prime_smul_of_minpoly_isEisensteinAt hp hBint hz
        (ih (pi • z) hpz h') hEis

/-! ## ★★★★★分母を払う——`d·z ∈ A[λ]` なる `d` は常に取れる -/

/-- ★★★★★**どの `z ∈ L` にも、`d·z ∈ A[λ]` となる `0 ≠ d ∈ A` がある**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★冪基底の係数は `K = Frac A` の元なので、**共通の分母を取る**だけである
(`IsLocalization.exist_integer_multiples`)。
★★`A` が DVR なら `d = π^k × 単元` なので、第 389 の降下の段に渡せる形になる
——この 2 つで **`B = A[λ]`** が出る。 -/
theorem exists_smul_mem_adjoin_powerBasis {R : Type*} {K : Type*} {L : Type*}
    [CommRing R] [Field K] [Field L] [Algebra K L] [Algebra R L] [Algebra R K]
    [IsScalarTower R K L] [IsDomain R] [IsFractionRing R K]
    (PB : PowerBasis K L) (z : L) :
    ∃ d : R, d ≠ 0 ∧ d • z ∈ Algebra.adjoin R {PB.gen} := by
  classical
  obtain ⟨b, hb⟩ := IsLocalization.exist_integer_multiples (nonZeroDivisors R)
    (Finset.univ : Finset (Fin PB.dim)) (fun i => PB.basis.repr z i)
  refine ⟨(b : R), nonZeroDivisors.coe_ne_zero b, ?_⟩
  have hz : ∑ i, (PB.basis.repr z) i • PB.basis i = z := PB.basis.sum_repr z
  rw [← hz, Finset.smul_sum]
  refine Subalgebra.sum_mem _ (fun i _ => ?_)
  obtain ⟨c, hc⟩ := hb i (Finset.mem_univ i)
  rw [PB.basis_eq_pow, ← smul_assoc, ← hc, IsScalarTower.algebraMap_smul]
  exact Subalgebra.smul_mem _
    (Subalgebra.pow_mem _ (Algebra.self_mem_adjoin_singleton R PB.gen) _) c

/-! ## ★★★★★★Eisenstein の定数項——`π²` では割れない -/

open Finset in
/-- ★★★★★★**定数項は `λ^e` の単元倍**——Eisenstein の『定数項が `π²` で割れない』。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★付値を使わず**割り切れだけ**で出る:
`a₀ = −λ^e − ∑_{i≥1} a_i λ^i` で、`i ≥ 1` の項はすべて `λ^{e+1}` で割れるので、
`a₀ = λ^e·(−1 − λ·s)` となり、`−1 − λ·s` は単元である。

★★これで **`λ^{e+1} ∤ a₀`**、とくに `e ≥ 1` なら `λ^{2e} ∤ a₀`、
すなわち `π² ∤ a₀` が出る。 -/
theorem constCoeff_eq_pow_mul_isUnit {B : Type*} [CommRing B] [IsLocalRing B]
    (e : ℕ) (he : 0 < e) (a : ℕ → B) (lam : B)
    (hlam : lam ∈ maximalIdeal B)
    (hdvd : ∀ i ∈ range e, lam ^ e ∣ a i)
    (hroot : lam ^ e + ∑ i ∈ range e, a i * lam ^ i = 0) :
    ∃ w : B, IsUnit w ∧ a 0 = lam ^ e * w := by
  classical
  -- ★`i ≥ 1` の項は `λ^{e+1}` で割れる
  have hIco : lam ^ (e + 1) ∣ ∑ i ∈ Ico 1 e, a i * lam ^ i := by
    refine Finset.dvd_sum (fun i hi => ?_)
    have hi1 : 1 ≤ i := (Finset.mem_Ico.mp hi).1
    have hie : i < e := (Finset.mem_Ico.mp hi).2
    obtain ⟨c, hc⟩ := hdvd i (Finset.mem_range.mpr hie)
    refine ⟨c * lam ^ (i - 1), ?_⟩
    rw [hc]
    have : lam ^ i = lam * lam ^ (i - 1) := by
      conv_lhs => rw [show i = 1 + (i - 1) by omega]
      rw [pow_add, pow_one]
    rw [this, pow_add, pow_one]
    ring
  obtain ⟨s, hs⟩ := hIco
  -- ★和を `a₀` と残りに分ける
  have hsplit : ∑ i ∈ range e, a i * lam ^ i
      = a 0 + ∑ i ∈ Ico 1 e, a i * lam ^ i := by
    have h0 : range e = insert 0 (Ico 1 e) := by
      ext x
      simp only [Finset.mem_range, Finset.mem_insert, Finset.mem_Ico]
      omega
    rw [h0, Finset.sum_insert (by simp)]
    simp
  rw [hsplit, hs] at hroot
  refine ⟨-1 - lam * s, ?_, ?_⟩
  · have hneg : IsUnit (-1 : B) := isUnit_one.neg
    have hmem : -(lam * s) ∈ maximalIdeal B :=
      Submodule.neg_mem _ (Ideal.mul_mem_right s _ hlam)
    have h2 := isUnit_add_mem_maximalIdeal hneg hmem
    simpa [sub_eq_add_neg] using h2
  · have hsplit2 : lam ^ (e + 1) = lam ^ e * lam := by rw [pow_succ]
    rw [hsplit2] at hroot
    have hz : a 0 = -(lam ^ e) - lam ^ e * lam * s := by
      have := hroot
      linear_combination this
    rw [hz]
    ring

/-! ## ★★★★★★★Eisenstein の係数——すべて `π` で割れる -/

open Finset in
/-- ★★★★★★★**全分岐の生成元の最小多項式の係数はすべて非単元**——
Eisenstein の『`i < e` の係数が `π` で割れる』。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★これも**付値なし**で出る——最小の反例 `i` を取ると:

* `λ^e` は `e > i` なので `λ^{i+1}` で割れる
* `j > i` の項は `λ^j` で割れ、`j ≥ i+1`
* `j < i` の項は最小性から `a_j` が非単元、すなわち `λ^e ∣ a_j` なので `λ^{e+j}` で割れ、`e+j > i`

★★したがって `a_i·λ^i` も `λ^{i+1}` で割れ、`λ^i` を約すと `λ ∣ a_i`、
つまり `a_i` は単元でない——矛盾である。

★★★仮説 `hpi` が「`A` の非単元は `λ^e` で割れる」であり、
これが**全分岐**(`π ~ λ^e`、第 388)の使い所である。 -/
theorem not_isUnit_coeff_of_root {A : Type*} {B : Type*} [CommRing A] [CommRing B]
    [IsDomain B] [Algebra A B]
    (e : ℕ) (a : ℕ → A) (lam : B) (hlam : lam ≠ 0) (hlamnu : ¬ IsUnit lam)
    (hpi : ∀ x : A, ¬ IsUnit x → lam ^ e ∣ algebraMap A B x)
    (hloc : ∀ x : A, ¬ IsUnit (algebraMap A B x) → ¬ IsUnit x)
    (hroot : lam ^ e + ∑ i ∈ range e, algebraMap A B (a i) * lam ^ i = 0) :
    ∀ i, i < e → ¬ IsUnit (a i) := by
  classical
  intro i
  induction i using Nat.strong_induction_on with
  | _ i ih =>
    intro hie hunit
    -- ★他の項はすべて `λ^{i+1}` で割れる
    have hrest : lam ^ (i + 1) ∣ ∑ j ∈ (range e).erase i, algebraMap A B (a j) * lam ^ j := by
      refine Finset.dvd_sum (fun j hj => ?_)
      have hjne : j ≠ i := (Finset.mem_erase.mp hj).1
      have hje : j < e := Finset.mem_range.mp (Finset.mem_erase.mp hj).2
      rcases lt_or_gt_of_ne hjne with hlt | hgt
      · obtain ⟨c, hc⟩ := hpi (a j) (ih j hlt hje)
        rw [hc]
        have h1 : lam ^ (i + 1) ∣ lam ^ (e + j) := pow_dvd_pow lam (by omega)
        have h2 : lam ^ e * c * lam ^ j = lam ^ (e + j) * c := by rw [pow_add]; ring
        rw [h2]
        exact h1.mul_right c
      · exact Dvd.dvd.mul_left (pow_dvd_pow lam (by omega)) _
    have hlead : lam ^ (i + 1) ∣ lam ^ e := pow_dvd_pow lam (by omega)
    have hsplit : ∑ j ∈ range e, algebraMap A B (a j) * lam ^ j
        = algebraMap A B (a i) * lam ^ i
          + ∑ j ∈ (range e).erase i, algebraMap A B (a j) * lam ^ j :=
      (Finset.add_sum_erase _ _ (Finset.mem_range.mpr hie)).symm
    rw [hsplit] at hroot
    have hkey : lam ^ (i + 1) ∣ algebraMap A B (a i) * lam ^ i := by
      have hneg : algebraMap A B (a i) * lam ^ i
          = -(lam ^ e + ∑ j ∈ (range e).erase i, algebraMap A B (a j) * lam ^ j) := by
        linear_combination hroot
      rw [hneg]
      exact dvd_neg.mpr (Dvd.dvd.add hlead hrest)
    obtain ⟨c, hc⟩ := hkey
    have hcancel : algebraMap A B (a i) = lam * c := by
      have hne : (lam : B) ^ i ≠ 0 := pow_ne_zero i hlam
      have heq : lam ^ i * algebraMap A B (a i) = lam ^ i * (lam * c) := by
        rw [pow_succ] at hc
        linear_combination hc
      exact mul_left_cancel₀ hne heq
    have hnu : ¬ IsUnit (algebraMap A B (a i)) := by
      rw [hcancel]
      intro hu
      exact hlamnu (isUnit_of_mul_isUnit_left hu)
    exact hloc (a i) hnu hunit

/-! ## ★★★★★★★★`minpoly` への配線 -/

open Polynomial Finset in
/-- ★★★★★★**モニックな多項式の根の関係式**——
`f` がモニックで次数 `e`、`λ` が根なら
`λ^e + ∑_{i<e} f.coeff i · λ^i = 0`。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★これが、本ファイルの Eisenstein の 2 本
(`constCoeff_eq_pow_mul_isUnit`･`not_isUnit_coeff_of_root`)への**入口**である。
★★mathlib の `Polynomial.Monic.as_sum` を使うだけである。 -/
theorem root_relation_of_monic {A : Type*} {B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    (f : A[X]) (hf : f.Monic) (e : ℕ) (hdeg : f.natDegree = e) (lam : B)
    (hroot : (aeval lam) f = 0) :
    lam ^ e + ∑ i ∈ range e, algebraMap A B (f.coeff i) * lam ^ i = 0 := by
  have h := hf.as_sum
  rw [hdeg] at h
  rw [h] at hroot
  simpa using hroot

open Polynomial Finset in
/-- ★★★★★★★**`minpoly` に当てた形**——`λ` が `A` 上整で最小多項式の次数が `e` なら
`λ^e + ∑_{i<e} (minpoly A λ).coeff i · λ^i = 0`。 -/
theorem root_relation_minpoly {A : Type*} {B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [Nontrivial B] (lam : B) (hint : IsIntegral A lam) (e : ℕ)
    (hdeg : (minpoly A lam).natDegree = e) :
    lam ^ e + ∑ i ∈ range e, algebraMap A B ((minpoly A lam).coeff i) * lam ^ i = 0 :=
  root_relation_of_monic (minpoly A lam) (minpoly.monic hint) e hdeg lam (minpoly.aeval A lam)

/-! ## ★★★★★★全分岐の仮説を局所環のデータから出す -/

/-- ★局所準同型の逆向きは常に成り立つ——単元は単元に移るから。 -/
theorem not_isUnit_of_not_isUnit_algebraMap {A : Type*} {B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] (x : A) (h : ¬ IsUnit (algebraMap A B x)) : ¬ IsUnit x :=
  fun hu => h (hu.map (algebraMap A B))

/-- ★★★★★★**全分岐の仮説を局所環のデータから出す**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★`A` が局所環で極大イデアルが `(π)`、かつ `λ^e` が `π` の単元倍なら、
`A` の非単元は `λ^e` で割れる——
`DifferentKummer.lean` の `mem_differentIdeal_of_totallyRamified_tame` が要求する形である。

★★仮説 `hmax` は「`A` は `π` を素元とする離散付値環」、
`hpi` は「全分岐(`π·B = 𝔪_B^e`)」である(第 388 が与える)。 -/
theorem pow_dvd_algebraMap_of_not_isUnit {A : Type*} {B : Type*} [CommRing A] [IsLocalRing A]
    [CommRing B] [Algebra A B] (e : ℕ) (pi : A) (lam u : B) (hu : IsUnit u)
    (hmax : IsLocalRing.maximalIdeal A = Ideal.span {pi})
    (hpi : lam ^ e = algebraMap A B pi * u) :
    ∀ x : A, ¬ IsUnit x → lam ^ e ∣ algebraMap A B x := by
  intro x hx
  have hmem : x ∈ IsLocalRing.maximalIdeal A := by
    by_contra hcon
    exact hx (IsLocalRing.notMem_maximalIdeal.mp hcon)
  rw [hmax, Ideal.mem_span_singleton] at hmem
  exact (pow_dvd_iff_of_isUnit_mul hu hpi).mpr (map_dvd (algebraMap A B) hmem)

/-! ## ★出典の紐付け(`.src`) -/

def IsTameDegree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1——馴分岐の定義。次数の側)",
    sectionId := "genell-prop-1-7" }

def isUnit_natCast_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1——IsUnit (n : B) が馴分岐であること)",
    sectionId := "genell-prop-1-7" }

def aeval_derivative_eisenstein.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(Eisenstein の導関数の値)",
    sectionId := "genell-prop-1-7" }

def aeval_derivative_eisenstein_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 1——馴なら different の指数はちょうど e−1)",
    sectionId := "genell-prop-1-7" }

def pow_dvd_iff_of_isUnit_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(全分岐で Eisenstein の条件が λ^e ∣ · になる)",
    sectionId := "genell-prop-1-7" }

def exists_isUnit_pow_eq_of_span_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(全分岐なら π は λ^e の単元倍)",
    sectionId := "genell-prop-1-7" }

def mem_adjoin_of_pow_smul_of_isEisensteinAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(B = A[λ] への降下の段)",
    sectionId := "genell-prop-1-7" }

def exists_smul_mem_adjoin_powerBasis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(B = A[λ] —— 分母を払う段)",
    sectionId := "genell-prop-1-7" }

def constCoeff_eq_pow_mul_isUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(Eisenstein の定数項は π² で割れない)",
    sectionId := "genell-prop-1-7" }

def not_isUnit_coeff_of_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(Eisenstein の係数はすべて π で割れる)",
    sectionId := "genell-prop-1-7" }

def root_relation_of_monic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(minpoly への配線——根の関係式)",
    sectionId := "genell-prop-1-7" }

def pow_dvd_algebraMap_of_not_isUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(全分岐の仮説を局所環のデータから)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
