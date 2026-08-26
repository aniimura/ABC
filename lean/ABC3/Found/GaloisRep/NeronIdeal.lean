import ABC3.Found.GaloisRep.NeronValuation
import Mathlib.RingTheory.DedekindDomain.Factorization

/-!
# Galois (G7) 第 316 ブロック —— **★★★★★★★★Néron 微分の分数イデアル**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `ω_E := ∏ᶠ p, p^{Néron 指数}` を**分数イデアルとして構成**し、
> **零でないこと**と**各素点での指数がちょうど Néron 指数であること**を示した。

★★★これが界面(第 311)の `omegaFrac` の中身である。

## ★★★★★★★★橋——`count` と付値

mathlib の分数イデアル論は `FractionalIdeal.count`(素因子の重複度)で書かれ、
私の Néron 指数は **adic 付値**で書かれている。両者はどちらも
`Associates.count` で定義されているので、**主分数イデアルの上で一致する**:

    count K p (spanSingleton x) = valAdd p x        (`count_spanSingleton_eq_valAdd`)

★★★★証明は `x = a/b` と書いて `count_well_defined` と `intValuationDef` を突き合わせるだけ。
★★これがあるので、以降は **mathlib の `count` の代数(`count_mul`・`count_zpow`・
`count_finprod`)がそのまま使える**。

## ★★★★★`log` で書き直しておく

`valAdd p x = -(v(x)).log`(`valAdd_eq_neg_log`)。
★`WithZero.log`(`exp` の逆)を使うと `unzero` の取り回しが消える。

## ★★★★構成と基本性質

* `omegaFracIdeal W := ∏ᶠ p, p^{neronExp p W}`(有限台なので実質有限積)
* `count_omegaFracIdeal`:`count p (ω_E) = neronExp p W`(`count_finprod`)
* `omegaFracIdeal_ne_zero`:有限個の非零因子の積だから非零
* 外延性(「`count` がすべて一致すれば分数イデアルは一致」)は
  **既に `Found/GaloisRep/D2Bridge.lean` にあった**(第 2 重複を這った——在庫の確認漏れ)

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valAdd_eq_neg_log` | ★★`valAdd` を `log` で書く |
| `count_spanSingleton_eq_valAdd` | ★★★★★★**`count` と付値の橋** |
| `omegaFracIdeal`・`count_omegaFracIdeal` | ★★★★★★★**分数イデアルの構成** |
| `omegaFracIdeal_ne_zero` | ★★★★非零性(外延性は既存の `D2Bridge.eq_of_count_eq`) |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve
open scoped nonZeroDivisors

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★`valAdd` を `log` で書く -/

theorem valAdd_eq_neg_log (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    valAdd p x = -((p.valuation L (x : L)).log) := by
  rw [valAdd]
  congr 1
  have h : ((WithZero.unzero (valuationP_ne_zero p x) : Multiplicative ℤ)
      : WithZero (Multiplicative ℤ)) = (p.valuation L) (x : L) := WithZero.coe_unzero _
  conv_rhs => rw [← h]
  rfl

/-! ## ★★★★★★`count` と付値の橋 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**主分数イデアルの `count` は加法付値そのもの**。

★どちらも `Associates.count` で定義されているので、`x = a/b` と書いて突き合わせる。 -/
theorem count_spanSingleton_eq_valAdd (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    FractionalIdeal.count L p (FractionalIdeal.spanSingleton (𝓞 L)⁰ (x : L)) = valAdd p x := by
  obtain ⟨a, b, hb, hab⟩ := IsFractionRing.div_surjective (𝓞 L) (x : L)
  have hb0 : b ≠ 0 := nonZeroDivisors.ne_zero hb
  have hbL : (algebraMap (𝓞 L) L) b ≠ 0 :=
    (map_ne_zero_iff _ (IsFractionRing.injective (𝓞 L) L)).2 hb0
  have hxne : (x : L) ≠ 0 := x.ne_zero
  have ha0 : a ≠ 0 := by
    intro h
    rw [h, map_zero, zero_div] at hab
    exact hxne hab.symm
  have hdec : FractionalIdeal.spanSingleton (𝓞 L)⁰ (x : L)
      = FractionalIdeal.spanSingleton (𝓞 L)⁰ ((algebraMap (𝓞 L) L) b)⁻¹
        * ↑(Ideal.span {a}) := by
    rw [FractionalIdeal.coeIdeal_span_singleton, FractionalIdeal.spanSingleton_mul_spanSingleton,
      ← hab]
    congr 1
    field_simp
  have hne : FractionalIdeal.spanSingleton (𝓞 L)⁰ (x : L) ≠ 0 :=
    FractionalIdeal.spanSingleton_ne_zero_iff.2 hxne
  rw [FractionalIdeal.count_well_defined L p hne hdec]
  have hval : (p.valuation L) (x : L)
      = WithZero.exp (((Associates.mk p.asIdeal).count (Associates.mk (Ideal.span {b})).factors : ℤ)
        - ((Associates.mk p.asIdeal).count (Associates.mk (Ideal.span {a})).factors : ℤ)) := by
    rw [← hab, map_div₀, IsDedekindDomain.HeightOneSpectrum.valuation_of_algebraMap,
      IsDedekindDomain.HeightOneSpectrum.valuation_of_algebraMap,
      IsDedekindDomain.HeightOneSpectrum.intValuation_apply,
      IsDedekindDomain.HeightOneSpectrum.intValuation_apply,
      IsDedekindDomain.HeightOneSpectrum.intValuationDef_if_neg _ ha0,
      IsDedekindDomain.HeightOneSpectrum.intValuationDef_if_neg _ hb0,
      WithZero.exp, WithZero.exp, WithZero.exp, ← WithZero.coe_div]
    congr 1
    rw [← ofAdd_sub]
    congr 1
    omega
  rw [valAdd_eq_neg_log, hval, WithZero.log_exp]
  omega

/-! ## ★★★★★★★分数イデアルの構成 -/

/-- ★★★★★★★**Néron 微分の分数イデアル** `ω_E = ∏ᶠ p^{Néron 指数}`。 -/
noncomputable def omegaFracIdeal (W : WeierstrassCurve L) : FractionalIdeal (𝓞 L)⁰ L :=
  ∏ᶠ p : HeightOneSpectrum (𝓞 L), (p.asIdeal : FractionalIdeal (𝓞 L)⁰ L) ^ (neronExp p W)

theorem neronExp_eventually_zero (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    ∀ᶠ p : HeightOneSpectrum (𝓞 L) in Filter.cofinite, neronExp p W = 0 := by
  rw [Filter.eventually_cofinite]
  exact finite_bad_primes' W hΔ

/-- ★★★★★各素点での重複度はちょうど Néron 指数。 -/
theorem count_omegaFracIdeal (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (p : HeightOneSpectrum (𝓞 L)) :
    FractionalIdeal.count L p (omegaFracIdeal W) = neronExp p W :=
  FractionalIdeal.count_finprod L p _ (neronExp_eventually_zero W hΔ)

/-- ★★★★零でない(有限個の非零因子の積だから)。 -/
theorem omegaFracIdeal_ne_zero (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    omegaFracIdeal W ≠ 0 := by
  rw [omegaFracIdeal]
  rw [finprod_eq_prod_of_mulSupport_subset _
    (s := (finite_bad_primes' W hΔ).toFinset) ?_]
  · rw [Finset.prod_ne_zero_iff]
    intro p _
    refine zpow_ne_zero _ ?_
    simpa using p.ne_bot
  · intro p hp
    simp only [Set.Finite.coe_toFinset, Set.mem_setOf_eq]
    intro hc
    apply hp
    show (p.asIdeal : FractionalIdeal (𝓞 L)⁰ L) ^ (neronExp p W) = 1
    rw [hc, zpow_zero]

/-! ## ★出典の紐付け(`.src`) -/

def count_spanSingleton_eq_valAdd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(主分数イデアルの count と加法付値の一致)",
    sectionId := "genell-def-3-3" }

def omegaFracIdeal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Néron 微分の分数イデアル)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
