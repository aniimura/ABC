/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.UnramQuad
import ABC3.Found.GaloisRep.RamifiedLocalData
import ABC3.Meta.Claim

/-!
# 第 1382 ブロック —— **非分岐 2 次拡大の段を `e` 倍版にする**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1029 の `e` 倍版

分裂性の二者択一（第 1035-1036）の**非分裂側**は
`R[X]/(f)`（`f` は分裂性の 2 次式の monic 化）へ上げる段である。
☆この拡大は**非分岐**（`e = 1`）なので `hp` はそのまま保たれる（第 1024）。

★したがって、下から来る `hpe`（`e` 倍）を合成すれば
`L → Frac(R[X]/(f))` でも `hpe`（同じ `e`）が成り立つ。

☆本ブロックは第 1029（`valuation_algebraMap_ext`・`isMinimal_baseChange_ext`・
`hasMultiplicativeReduction_ext`）の `e` 倍版と、
`R[X]/(f)` 自身の `m`-進完備性（第 1013 ＋ 第 1022）を作る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Polynomial

open scoped Classical

section AdicCompleteAdjoinRoot

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]

/-- ★★★★★★★★★★★★
**`R[X]/(f)` はそれ自身の極大イデアルに沿って完備**——★**無条件**（第 1382）。

☆第 1013（`m_R`-進完備）と第 1022（`m_{R′} = m_R R′`）を繋ぐだけである。 -/
theorem isAdicComplete_maximalIdeal_adjoinRoot
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsLocalRing (AdjoinRoot f)] :
    IsAdicComplete (IsLocalRing.maximalIdeal (AdjoinRoot f)) (AdjoinRoot f) := by
  haveI : IsAdicComplete (IsLocalRing.maximalIdeal R) (AdjoinRoot f) :=
    isAdicComplete_adjoinRoot (IsLocalRing.maximalIdeal R) hf hdeg
  rw [maximalIdeal_adjoinRoot_eq_map hf hdeg hirr]
  exact isAdicComplete_map_algebraMap (IsLocalRing.maximalIdeal R)

end AdicCompleteAdjoinRoot

section ExtRam

variable {L : Type} [Field L] [NumberField L]
variable {Lv : Type} [Field Lv] [Algebra L Lv]
variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
variable [Algebra R Lv] [IsFractionRing R Lv] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]

/-- ★★★★★★★★★★★★**拡大の上でも `hpe` は成り立つ**（第 1382）。

☆非分岐なので `e` は変わらない。 -/
theorem valuation_algebraMap_ext_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (x : L) :
    (HeightOneSpectrum.valuation (FractionRing (AdjoinRoot f))
        (IsDiscreteValuationRing.maximalIdeal (AdjoinRoot f)))
      (algebraMap L (FractionRing (AdjoinRoot f)) x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e := by
  rw [halg, valuation_quadFieldHom hf hdeg hirr, hpe]

/-- ★★★★★★★★**極小性は拡大にも移る**——`e` 倍版（第 1382）。 -/
theorem isMinimal_baseChange_ext_ram {e : ℕ} (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0) :
    WeierstrassCurve.IsMinimal (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f))) :=
  isMinimal_baseChange_at_bad_prime_ram p
    (valuation_algebraMap_ext_ram p hpe hf hdeg hirr halg) E C hC hc4ne hc4

/-- ★★★★★★★★**乗法還元も拡大に移る**——`e` 倍版（第 1382）。 -/
theorem hasMultiplicativeReduction_ext_ram {e : ℕ} (he : 1 ≤ e)
    (p : HeightOneSpectrum (𝓞 L))
    (hpe : ∀ x : L, (HeightOneSpectrum.valuation Lv
        (IsDiscreteValuationRing.maximalIdeal R)) (algebraMap L Lv x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e)
    {f : R[X]} (hf : f.Monic) (hdeg : f.natDegree = 2)
    (hirr : Irreducible (f.map (Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))))
    [IsDomain (AdjoinRoot f)] [IsDiscreteValuationRing (AdjoinRoot f)]
    [Algebra L (FractionRing (AdjoinRoot f))]
    (halg : ∀ x : L, algebraMap L (FractionRing (AdjoinRoot f)) x
      = quadFieldHom hf hdeg (algebraMap L Lv x))
    (E : WeierstrassCurve L) [E.IsElliptic]
    (C : WeierstrassCurve.VariableChange L)
    (hC : WeierstrassCurve.IsMinimal (primeSubring p) (C • E))
    (hc4ne : (C • E).c₄ ≠ 0) (hc4 : valAdd p (Units.mk0 ((C • E).c₄) hc4ne) = 0)
    (hj : jExp p E < 0)
    [WeierstrassCurve.IsMinimal (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f)))] :
    WeierstrassCurve.HasMultiplicativeReduction (AdjoinRoot f)
      ((C • E).baseChange (FractionRing (AdjoinRoot f))) :=
  hasMultiplicativeReduction_at_bad_prime_ram he p
    (valuation_algebraMap_ext_ram p hpe hf hdeg hirr halg) E C hC hc4ne hc4 hj

end ExtRam

/-! ## ★出典の紐付け(`.src`) -/

def isAdicComplete_maximalIdeal_adjoinRoot.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(R[X]/(f) はそれ自身の極大イデアルに沿って完備。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valuation_algebraMap_ext_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(拡大の上でも hpe は成り立つ。非分岐なので e は変わらない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isMinimal_baseChange_ext_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極小性は拡大にも移る。e 倍版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def hasMultiplicativeReduction_ext_ram.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(乗法還元も拡大に移る。e 倍版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valuation_algebraMap_ext_ram.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "valuation_quadFieldHom(第 1024、証明済み。非分岐なので等号)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.valuation_quadFieldHom") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1382）**——第 1029 の `e` 倍版である。" ++
       "☆非分岐 2 次拡大は `e = 1` なので、下から来る `e` がそのまま伝わる。") 17 ]

end ABC3.Found.GaloisRep
