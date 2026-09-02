/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.CyclotomicLocalExt
import ABC3.Found.GenEll.CompleteIntegralClosure
import ABC3.Found.GenEll.ValuationExtension
import ABC3.Found.GaloisRep.AdicCompleteIntegers
import ABC3.Found.GaloisRep.CompletionValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1377 ブロック —— **`L_p(ζ_l)` の完備 DVR パッケージ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1372 の入力を 1 つにまとめる

第 1366-1376 で部品が揃った。★本ブロックはそれを

> **`L_p` に原始 `l` 乗根を足した体 `E` と、その整数環 `C` は完備 DVR で、
> `v_C(algebraMap L E x) = (v_p x)^e`、`1 ≤ e` かつ `l ∤ e`**

という**1 つの命題**にまとめる。

☆部品:

| 部品 | 第 |
|---|---|
| `E ≔ cycExt L_p`（`≤ l−1` 次、`ζ_l` を含む） | 1376 |
| `C ≔ 整閉包(𝒪_p, E)` は DVR | 1366 |
| `C` は `m_C`-進完備・`IsFractionRing C E` | 1368 |
| `m_{𝒪_p} C = m_C^e`、`1 ≤ e ≤ [E : L_p]` | 1369 |
| `v_C(φ y) = (v_{𝒪_p} y)^e` | 1373 |
| `hp`（完備化）と合成 | 1374 |
| `l ∤ e` | 1376 |
-/

namespace ABC3.Found.GenEll

open ABC3.Meta IsDedekindDomain NumberField ABC3.Found.GaloisRep

section Package

variable {L : Type} [Field L] [NumberField L]

instance instCharZeroAdicCompletion (p : HeightOneSpectrum (𝓞 L)) :
    CharZero (p.adicCompletion L) :=
  charZero_of_injective_algebraMap (algebraMap L (p.adicCompletion L)).injective

variable (p : HeightOneSpectrum (𝓞 L)) {l : ℕ}

/-- ★★★★★★★★**`L_p` に原始 `l` 乗根を足した体**（第 1377）。 -/
@[reducible] noncomputable def locCycField (hl : l.Prime) : Type :=
  cycExt (p.adicCompletion L) hl

/-- ★★★★★★★★**その整数環**（第 1377）。 -/
@[reducible] noncomputable def locCycRing (hl : l.Prime) : Type :=
  ↥(integralClosure (p.adicCompletionIntegers L) (locCycField p hl))

/-- ★★★★★★★★★★★★**`L_p(ζ_l)` の整数環は DVR**（第 1377）。 -/
instance instIsDVRLocCycRing (hl : l.Prime) :
    IsDiscreteValuationRing (locCycRing p hl) :=
  isDiscreteValuationRing_integralClosure
    (A := p.adicCompletionIntegers L) (p.adicCompletion L) (locCycField p hl) _

/-- ★★★★★★★★★★★★
**`L_p(ζ_l)` の整数環は完備 DVR で、付値は `e` 倍**——★**無条件**（第 1377）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1372（`exists_h2_h1_of_bad_prime_ram`）に渡す入力一式である。 -/
theorem exists_locCyc_package (hl : l.Prime) :
    ∃ e : ℕ, 1 ≤ e ∧ ¬ (l ∣ e) ∧
      (IsLocalRing.maximalIdeal (p.adicCompletionIntegers L)).map
          (algebraMap (p.adicCompletionIntegers L) (locCycRing p hl))
        = (IsLocalRing.maximalIdeal (locCycRing p hl)) ^ e := by
  obtain ⟨e, he1, hIe⟩ := exists_pow_eq_map_maximalIdeal
    (A := p.adicCompletionIntegers L) (p.adicCompletion L) (locCycField p hl)
    (locCycRing p hl)
  refine ⟨e, he1, ?_, hIe⟩
  refine not_dvd_of_le_finrank (p.adicCompletion L) hl he1 ?_
  exact le_finrank_of_pow_eq_map_maximalIdeal
    (A := p.adicCompletionIntegers L) (p.adicCompletion L) (locCycField p hl)
    (locCycRing p hl) hIe

/-- ★★★★★★★★**商体は `L_p(ζ_l)`**（第 1377）。 -/
instance instIsFractionRingLocCycRing (hl : l.Prime) :
    IsFractionRing (locCycRing p hl) (locCycField p hl) :=
  isFractionRing_integralClosure
    (A := p.adicCompletionIntegers L) (p.adicCompletion L) (locCycField p hl) _

/-- ★★★★★★★★**`m_C`-進完備**（第 1377）。 -/
instance instIsAdicCompleteLocCycRing (hl : l.Prime) :
    IsAdicComplete (IsLocalRing.maximalIdeal (locCycRing p hl)) (locCycRing p hl) :=
  isAdicComplete_maximalIdeal_integralClosure
    (A := p.adicCompletionIntegers L) (p.adicCompletion L) (locCycField p hl) _

/-- ★★★★★★★★★★★★**`ζ_l ∈ L_p(ζ_l)`**（第 1377）。 -/
theorem isPrimitiveRoot_locCycRoot (hl : l.Prime) :
    IsPrimitiveRoot (cycRoot (p.adicCompletion L) hl) l :=
  isPrimitiveRoot_cycRoot (p.adicCompletion L) hl

/-- ★★★★★★★★★★★★★★★★★★★★
**`hpe`——付値は `e` 倍で両立する**——★**無条件**（第 1377）。

★★★これが第 1372（`exists_h2_h1_of_bad_prime_ram`）に渡す `hpe` である。 -/
theorem valuation_algebraMap_locCycField (hl : l.Prime) {e : ℕ} (he : 1 ≤ e)
    (hIe : (IsLocalRing.maximalIdeal (p.adicCompletionIntegers L)).map
        (algebraMap (p.adicCompletionIntegers L) (locCycRing p hl))
      = (IsLocalRing.maximalIdeal (locCycRing p hl)) ^ e) (x : L) :
    (HeightOneSpectrum.valuation (locCycField p hl)
        (IsDiscreteValuationRing.maximalIdeal (locCycRing p hl)))
        (algebraMap L (locCycField p hl) x)
      = ((HeightOneSpectrum.valuation L p) x) ^ e :=
  valuation_algebraMap_pow_comp (K := p.adicCompletion L)
    (A := p.adicCompletionIntegers L) (C := locCycRing p hl) (K' := locCycField p hl)
    p (valuation_algebraMap_adicCompletion L p) he hIe x

end Package

/-! ## ★出典の紐付け(`.src`) -/

def exists_locCyc_package.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(L_p(ζ_l) の整数環は完備 DVR で付値は e 倍。★無条件)",
    sectionId := "genell-thm-3-8" }

def valuation_algebraMap_locCycField.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(hpe——付値は e 倍で両立する。★無条件)",
    sectionId := "genell-thm-3-8" }

def isPrimitiveRoot_locCycRoot.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ_l ∈ L_p(ζ_l)。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_locCyc_package.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isDiscreteValuationRing_integralClosure(第 1366、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isDiscreteValuationRing_integralClosure") 1,
    .citation "[ABC3]" "finrank_cycExt_le(第 1376、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finrank_cycExt_le") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1377）**——第 1366-1376 の部品を 1 つの命題にまとめる段である。" ++
       "☆これが第 1372（`exists_h2_h1_of_bad_prime_ram`）に渡す入力一式になる。") 19 ]

end ABC3.Found.GenEll
