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
