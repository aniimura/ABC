import ABC3.Meta.Claim
import Mathlib.RingTheory.LocalRing.ResidueField.Defs
import Mathlib.RingTheory.LocalRing.Basic
import Mathlib.Algebra.CharP.Basic
import Mathlib.RingTheory.Polynomial.Eisenstein.Basic
import Mathlib.RingTheory.PrincipalIdealDomain

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

end ABC3.Found.GenEll
