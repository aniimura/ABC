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
abbrev HeightOneSpectrumUnder [Algebra.IsIntegral A B]
    (P : HeightOneSpectrum B) : HeightOneSpectrum A :=
  HeightOneSpectrum.under A P

theorem HeightOneSpectrumUnder_asIdeal [Algebra.IsIntegral A B]
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

/-! ## ★★★★★★★★★★★★繊維分解 -/

open scoped Classical in
/-- ★★★★★`primesOverH` に入れるための十分条件。 -/
theorem mem_primesOverH_of_mem (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L'))
    (hP : P.asIdeal ∈ IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L')) :
    P ∈ primesOverH L L' p :=
  Finset.mem_image.2 ⟨⟨P.asIdeal, hP⟩, Finset.mem_attach _ _, HeightOneSpectrum.ext rfl⟩

open scoped Classical in
/-- ★★★★★★★★**`primesOverH L L' p` に入ることは「下が `p`」と同値**。

★これが繊維分解の要——`primesOverH` たちは互いに素で、全体を覆う。 -/
theorem mem_primesOverH_iff (p : HeightOneSpectrum (𝓞 L))
    (P : HeightOneSpectrum (𝓞 L')) :
    P ∈ primesOverH L L' p ↔ p = HeightOneSpectrumUnder (A := 𝓞 L) P := by
  constructor
  · intro hP
    haveI hmax : p.asIdeal.IsMaximal := p.isPrime.isMaximal p.ne_bot
    have h := (IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1
      (mem_primesOverH L L' p hP)
    exact HeightOneSpectrum.ext h.2.over
  · rintro rfl
    haveI hmax : (HeightOneSpectrumUnder (A := 𝓞 L) P).asIdeal.IsMaximal :=
      (HeightOneSpectrumUnder (A := 𝓞 L) P).isPrime.isMaximal
        (HeightOneSpectrumUnder (A := 𝓞 L) P).ne_bot
    refine mem_primesOverH_of_mem L L' _ P ?_
    exact (IsDedekindDomain.mem_primesOverFinset_iff
      (HeightOneSpectrumUnder (A := 𝓞 L) P).ne_bot (𝓞 L')).2 ⟨P.isPrime, ⟨rfl⟩⟩

open scoped Classical in
/-- ★★★★★★異なる `p` の繊維は交わらない。 -/
theorem primesOverH_disjoint {p q : HeightOneSpectrum (𝓞 L)} (hpq : p ≠ q) :
    Disjoint (primesOverH L L' p) (primesOverH L L' q) := by
  rw [Finset.disjoint_left]
  intro P hp hq
  exact hpq (((mem_primesOverH_iff L L' p P).1 hp).trans
    ((mem_primesOverH_iff L L' q P).1 hq).symm)

open scoped Classical in
/-- ★★★★★★`primesOverH` の `asIdeal` 像は `primesOverFinset` そのもの。 -/
theorem image_asIdeal_primesOverH (p : HeightOneSpectrum (𝓞 L)) :
    (primesOverH L L' p).image HeightOneSpectrum.asIdeal
      = IsDedekindDomain.primesOverFinset p.asIdeal (𝓞 L') := by
  haveI hmax : p.asIdeal.IsMaximal := p.isPrime.isMaximal p.ne_bot
  ext I
  constructor
  · intro hI
    obtain ⟨P, hP, rfl⟩ := Finset.mem_image.1 hI
    exact mem_primesOverH L L' p hP
  · intro hI
    have hmem := (IsDedekindDomain.mem_primesOverFinset_iff p.ne_bot (𝓞 L')).1 hI
    haveI : I.IsPrime := hmem.1
    haveI : I.LiesOver p.asIdeal := hmem.2
    exact Finset.mem_image.2
      ⟨⟨I, hmem.1, Ideal.ne_bot_of_liesOver_of_ne_bot p.ne_bot I⟩,
        mem_primesOverH_of_mem L L' p _ hI, rfl⟩

open scoped Classical in
/-- ★★★★★★★★★★★★第 734 を `HeightOneSpectrum` の言葉で。 -/
theorem sum_primesOverH_ramificationIdx_log (p : HeightOneSpectrum (𝓞 L)) :
    ∑ P ∈ primesOverH L L' p,
        (p.asIdeal.ramificationIdx P.asIdeal : ℝ) * Real.log (Ideal.absNorm P.asIdeal)
      = (Module.finrank L L' : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  rw [← sum_primesOver_ramificationIdx_log L L' p, ← image_asIdeal_primesOverH L L' p,
    Finset.sum_image (fun P _ Q _ h => HeightOneSpectrum.ext h)]

open scoped Classical in
/-- ★★★★★★★★★★**`finsum` を `p` ごとの繊維に束ねる**。 -/
theorem finsum_eq_sum_primesOverH (f : HeightOneSpectrum (𝓞 L') → ℝ)
    (S : Finset (HeightOneSpectrum (𝓞 L)))
    (hsupp : ∀ P, f P ≠ 0 → HeightOneSpectrumUnder (A := 𝓞 L) P ∈ S) :
    ∑ᶠ P, f P = ∑ p ∈ S, ∑ P ∈ primesOverH L L' p, f P := by
  rw [← Finset.sum_biUnion (fun p _ q _ hpq => primesOverH_disjoint L L' hpq)]
  refine finsum_eq_finset_sum_of_support_subset f ?_
  intro P hP
  simp only [Function.mem_support] at hP
  exact Finset.mem_coe.2 (Finset.mem_biUnion.2
    ⟨_, hsupp P hP, (mem_primesOverH_iff L L' _ P).2 rfl⟩)

open scoped Classical in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`finsum` の scaling 補題**——`f′(P) = e(P|p)·f(p)` なら

    ∑ᶠ_P f′(P)·log N(P) = [L′ : L] · ∑ᶠ_p f(p)·log N(p)

★★これが `degInfOf` の基底変換不変性（したがって `ht^Falt` の不変性）の核である。
`ResearchPaper/ellmoduli-witness-status.json` の `simplification20260830` で
「難所は `finsum` の scaling 補題 1 本」と絞り込んだ、その本体。 -/
theorem finsum_scaling (f : HeightOneSpectrum (𝓞 L) → ℝ)
    (f' : HeightOneSpectrum (𝓞 L') → ℝ)
    (hf' : ∀ P : HeightOneSpectrum (𝓞 L'),
      f' P = ((HeightOneSpectrumUnder (A := 𝓞 L) P).asIdeal.ramificationIdx P.asIdeal : ℝ)
        * f (HeightOneSpectrumUnder (A := 𝓞 L) P))
    (S : Finset (HeightOneSpectrum (𝓞 L))) (hS : ∀ p, f p ≠ 0 → p ∈ S) :
    ∑ᶠ P, f' P * Real.log (Ideal.absNorm P.asIdeal)
      = (Module.finrank L L' : ℝ)
        * ∑ᶠ p, f p * Real.log (Ideal.absNorm p.asIdeal) := by
  have hsub : ∀ P : HeightOneSpectrum (𝓞 L'),
      f' P * Real.log (Ideal.absNorm P.asIdeal) ≠ 0 →
      HeightOneSpectrumUnder (A := 𝓞 L) P ∈ S := by
    intro P hP
    refine hS _ (fun h0 => hP ?_)
    rw [hf' P, h0, mul_zero, zero_mul]
  have hR : ∑ᶠ p, f p * Real.log (Ideal.absNorm p.asIdeal)
      = ∑ p ∈ S, f p * Real.log (Ideal.absNorm p.asIdeal) := by
    refine finsum_eq_finset_sum_of_support_subset _ ?_
    intro p hp
    simp only [Function.mem_support] at hp
    exact Finset.mem_coe.2 (hS p (fun h0 => hp (by rw [h0, zero_mul])))
  rw [finsum_eq_sum_primesOverH L L' _ S hsub, hR, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun p _ => ?_)
  have hin : ∀ P ∈ primesOverH L L' p,
      f' P * Real.log (Ideal.absNorm P.asIdeal)
        = f p * ((p.asIdeal.ramificationIdx P.asIdeal : ℝ)
            * Real.log (Ideal.absNorm P.asIdeal)) := by
    intro P hP
    have hu : HeightOneSpectrumUnder (A := 𝓞 L) P = p :=
      ((mem_primesOverH_iff L L' p P).1 hP).symm
    rw [hf' P, hu]; ring
  rw [Finset.sum_congr rfl hin, ← Finset.mul_sum,
    sum_primesOverH_ramificationIdx_log L L' p]
  ring

def finsum_scaling.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(finsum の scaling 補題——f′ = e·f なら和は [L′:L] 倍。★無条件)",
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
