/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInf
import Mathlib.NumberTheory.RamificationInertia.Inertia
import ABC3.Meta.Claim

/-!
# 素点の剰余標数（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.22。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

## ★★★★★★★★これは何か

原文 p.22 は

> write `S∘` for the union of `S`, the primes of **`ℚ`** that lie under primes of
> potentially multiplicative reduction of `E_L`, …

と、**`ℚ` の素数**で帳簿を作る。`Interface/GenEll/EllModuli.lean` の `multPrime` 欄も
`Nat.Prime` を要求する。★そこで `P : HeightOneSpectrum (𝓞 L)` から
その下の有理素数を取り出す道具を作る。

## ★★★機構

* `P ∩ ℤ` は `ℤ` の非零素イデアル（`Ideal.under_ne_bot`）
* `ℤ` は PID なので生成元があり、その `natAbs` が求める素数
* `absNorm P = p^f`（mathlib の `Ideal.absNorm_eq_pow_inertiaDeg'`）
  なので `log N(P) = f·log p`

★★これで `d·deg∞ = ∑_p v_p(Δ_min)·f_p·log(剰余標数)` と書け、
`sum_localHt_eq` 欄（`∑ h_j·log(p_j) = 23040·d·deg∞`）が組める。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★剰余標数 -/

/-- ★★★★**素点 `P` の剰余標数**——`P ∩ ℤ` の生成元の絶対値。 -/
noncomputable def resChar (P : HeightOneSpectrum (𝓞 L)) : ℕ :=
  (Submodule.IsPrincipal.generator (P.asIdeal.under ℤ)).natAbs

theorem under_ne_bot (P : HeightOneSpectrum (𝓞 L)) : P.asIdeal.under ℤ ≠ ⊥ :=
  Ideal.under_ne_bot (A := ℤ) P.ne_bot

instance under_isPrime (P : HeightOneSpectrum (𝓞 L)) : (P.asIdeal.under ℤ).IsPrime := by
  haveI := P.isPrime
  exact Ideal.comap_isPrime (algebraMap ℤ (𝓞 L)) P.asIdeal

theorem span_generator_under (P : HeightOneSpectrum (𝓞 L)) :
    Ideal.span {Submodule.IsPrincipal.generator (P.asIdeal.under ℤ)} = P.asIdeal.under ℤ :=
  Ideal.span_singleton_generator _

theorem generator_ne_zero (P : HeightOneSpectrum (𝓞 L)) :
    Submodule.IsPrincipal.generator (P.asIdeal.under ℤ) ≠ 0 := by
  intro h
  apply under_ne_bot P
  rw [← span_generator_under P, h]
  simp

/-- ★★★★★**剰余標数は素数である**。 -/
theorem resChar_prime (P : HeightOneSpectrum (𝓞 L)) : Nat.Prime (resChar P) := by
  have hprime : Prime (Submodule.IsPrincipal.generator (P.asIdeal.under ℤ)) := by
    rw [← Ideal.span_singleton_prime (generator_ne_zero P), span_generator_under P]
    exact under_isPrime P
  exact Int.prime_iff_natAbs_prime.1 hprime

theorem span_resChar (P : HeightOneSpectrum (𝓞 L)) :
    Ideal.span {(resChar P : ℤ)} = P.asIdeal.under ℤ := by
  rw [resChar, Int.span_natAbs, span_generator_under]

/-- ★★★★`P` は剰余標数の上にある。 -/
instance liesOver_resChar (P : HeightOneSpectrum (𝓞 L)) :
    P.asIdeal.LiesOver (Ideal.span {(resChar P : ℤ)}) :=
  ⟨span_resChar P⟩

/-! ## ★★★★★★`N(P) = p^f` -/

/-- ★★★★★★**`N(P) = p^{f(P)}`**（`p` は剰余標数）。 -/
theorem absNorm_eq_resChar_pow (P : HeightOneSpectrum (𝓞 L)) :
    Ideal.absNorm P.asIdeal
      = resChar P ^ ((Ideal.span {(resChar P : ℤ)}).inertiaDeg P.asIdeal) :=
  Ideal.absNorm_eq_pow_inertiaDeg' P.asIdeal (resChar_prime P)

/-- ★★★★★★★**`log N(P) = f(P)·log p`**——帳簿の要。 -/
theorem log_absNorm_eq (P : HeightOneSpectrum (𝓞 L)) :
    Real.log (Ideal.absNorm P.asIdeal)
      = ((Ideal.span {(resChar P : ℤ)}).inertiaDeg P.asIdeal : ℝ) * Real.log (resChar P) := by
  rw [absNorm_eq_resChar_pow P]
  push_cast
  rw [Real.log_pow]

/-- ★★★★慣性次数は正である（`N(P) ≥ 2` だから）。 -/
theorem inertiaDeg_pos (P : HeightOneSpectrum (𝓞 L)) :
    0 < (Ideal.span {(resChar P : ℤ)}).inertiaDeg P.asIdeal := by
  rcases Nat.eq_zero_or_pos ((Ideal.span {(resChar P : ℤ)}).inertiaDeg P.asIdeal) with h0 | h0
  · exfalso
    have h1 := absNorm_eq_resChar_pow P
    rw [h0, pow_zero] at h1
    exact P.isPrime.ne_top (Ideal.absNorm_eq_one_iff.1 h1)
  · exact h0

/-! ## ★出典の紐付け(`.src`) -/

def resChar.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(素点の剰余標数——原文の「the primes of ℚ that lie under…」)",
    sectionId := "genell-cor-4-3" }

def log_absNorm_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(log N(P) = f(P)·log p。★無条件)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GaloisRep
