/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInf
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.NumberTheory.RamificationInertia.Basic
import ABC3.Meta.Claim

/-!
# 素点の「下」—— `HeightOneSpectrum` の contraction（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★これは何か

`ResearchPaper/ellmoduli-witness-status.json` の `simplification20260830` で
「難所は `finsum` の scaling 補題 1 本」と絞り込んだ、その**最初の部品**である。

☆mathlib の `IsDedekindDomain.HeightOneSpectrum.comap` は**全射**な環準同型に対する
ものなので（`Ideal/Lemmas.lean:588`）、`𝓞L → 𝓞L′`（単射だが全射でない）には使えない。
★そこで「`P` の下の `p`」を自分で作る。

## ★機構

* 素性は `Ideal.comap` で保たれる
* 非零性は `Ideal.under_ne_bot`（`Mathlib/RingTheory/Ideal/GoingUp.lean:422`）
  ——整拡大なら `P ≠ ⊥` から `P ∩ 𝓞L ≠ ⊥`
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField

variable {A B : Type*} [CommRing A] [IsDomain A] [IsDedekindDomain A]
  [CommRing B] [IsDomain B] [IsDedekindDomain B] [Algebra A B]

/-- ★★★★★★**素点の「下」**——`P ∩ A`。

★`𝓞L → 𝓞L′` は全射でないので mathlib の `HeightOneSpectrum.comap` は使えない。
本定義がその代わりである。 -/
def HeightOneSpectrumUnder [Algebra.IsIntegral A B]
    (P : HeightOneSpectrum B) : HeightOneSpectrum A where
  asIdeal := P.asIdeal.under A
  isPrime := by
    haveI := P.isPrime
    exact Ideal.comap_isPrime (algebraMap A B) P.asIdeal
  ne_bot := Ideal.under_ne_bot (A := A) P.ne_bot

@[simp] theorem HeightOneSpectrumUnder_asIdeal [Algebra.IsIntegral A B]
    (P : HeightOneSpectrum B) :
    (HeightOneSpectrumUnder P).asIdeal = P.asIdeal.under A := rfl

/-- ★★★★`P` はその「下」の上にある。 -/
instance liesOver_under [Algebra.IsIntegral A B]
    (P : HeightOneSpectrum B) :
    P.asIdeal.LiesOver (HeightOneSpectrumUnder (A := A) P).asIdeal :=
  ⟨rfl⟩

/-! ## ★★★★★★★★★★素点の上での和 -/

section NumberField

open scoped Classical

variable (L L' : Type) [Field L] [NumberField L] [Field L'] [NumberField L']
  [Algebra L L'] [Algebra (𝓞 L) L'] [IsScalarTower (𝓞 L) L L']
  [IsScalarTower (𝓞 L) (𝓞 L') L'] [Module.Finite (𝓞 L) (𝓞 L')]

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★
**`p` の上の素点にわたる `e·log N(P)` の和は `[L′ : L]·log N(p)`**

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★これが「`finsum` の scaling 補題」（`§9-1152`、第 730）の**数学の核**である。

★機構は 2 つの mathlib 補題:

* `log N(P) = f(P|p)·log N(p)`（`Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver`）
* `Σ_{P|p} e·f = [L′ : L]`（`Ideal.sum_ramification_inertia`） -/
theorem sum_primesOver_ramificationIdx_log (p : HeightOneSpectrum (𝓞 L)) :
    ∑ P ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L'),
        (p.asIdeal.ramificationIdx P : ℝ) * Real.log (Ideal.absNorm P)
      = (Module.finrank L L' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  haveI hmax : p.asIdeal.IsMaximal := p.isPrime.isMaximal p.ne_bot
  have hlog : ∀ P ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L'),
      (p.asIdeal.ramificationIdx P : ℝ) * Real.log (Ideal.absNorm P)
        = ((p.asIdeal.ramificationIdx P * p.asIdeal.inertiaDeg P : ℕ) : ℝ)
          * Real.log (Ideal.absNorm p.asIdeal) := by
    intro P hP
    haveI : P.LiesOver p.asIdeal :=
      ((IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1 hP).2
    rw [Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver P p.asIdeal p.isPrime p.ne_bot]
    push_cast
    rw [Real.log_pow]
    ring
  rw [Finset.sum_congr rfl hlog, ← Finset.sum_mul, ← Nat.cast_sum,
    Ideal.sum_ramification_inertia (R := 𝓞 L) (S := 𝓞 L') L L' p.ne_bot]

def sum_primesOver_ramificationIdx_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(p の上の素点にわたる e·log N(P) の和は [L′:L]·log N(p)。★無条件)",
    sectionId := "genell-prop-3-4" }

/-! ## ★★★★★★★★`HeightOneSpectrum` としての繊維 -/

open scoped Classical in
/-- ★★★★★★**`p` の上の素点を `HeightOneSpectrum` の `Finset` として**。

☆mathlib の `IsDedekindDomain.primesOverFinset` は `Finset (Ideal B)` なので、
`degInfOf` の添字（`HeightOneSpectrum`）に合わせるために像を取る。 -/
noncomputable def primesOverH (p : HeightOneSpectrum (𝓞 L)) :
    Finset (HeightOneSpectrum (𝓞 L')) :=
  haveI hmax : p.asIdeal.IsMaximal := p.isPrime.isMaximal p.ne_bot
  (IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L')).attach.image
    (fun P : {x : Ideal (𝓞 L') //
        x ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L')} =>
      { asIdeal := P.1
        isPrime :=
          ((IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1 P.2).1
        ne_bot := by
          haveI : P.1.IsPrime :=
            ((IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1 P.2).1
          haveI : P.1.LiesOver p.asIdeal :=
            ((IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1 P.2).2
          exact Ideal.ne_bot_of_liesOver_of_ne_bot p.ne_bot P.1 })

open scoped Classical in
/-- ★★★★★`primesOverH` の元は `p` の上にある。 -/
theorem mem_primesOverH (p : HeightOneSpectrum (𝓞 L))
    {P : HeightOneSpectrum (𝓞 L')} (hP : P ∈ primesOverH L L' p) :
    P.asIdeal ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L') := by
  obtain ⟨Q, -, hQ⟩ := Finset.mem_image.1 hP
  rw [← hQ]
  exact Q.2

def primesOverH.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(p の上の素点を HeightOneSpectrum の Finset として)",
    sectionId := "genell-prop-3-4" }

end NumberField

/-! ## ★出典の紐付け(`.src`) -/

def HeightOneSpectrumUnder.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(素点の「下」——HeightOneSpectrum の contraction。★無条件)",
    sectionId := "genell-prop-3-4" }

def HeightOneSpectrumUnder.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.under_ne_bot(整拡大なら P ≠ ⊥ から P ∩ A ≠ ⊥)"
      (.inMathlib "Ideal.under_ne_bot") 1,
    .implicitStep
      ("☆mathlib の IsDedekindDomain.HeightOneSpectrum.comap は全射な環準同型に" ++
       "対するものなので 𝓞L → 𝓞L′ には使えない——本定義がその代わりである") 1 ]

end ABC3.Found.GaloisRep
