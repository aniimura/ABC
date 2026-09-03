/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.BadPrimeData
import ABC3.Found.GaloisRep.CompletionValuationBridge
import ABC3.Meta.Claim

/-!
# 第 1350 ブロック —— **乗法還元から悪い素点の局所データを出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`α` の側の入口

`EllModuliWitness` の `imageContainsSL2J_torsionExt`（#4）は
`HasMultRed` から出発する。★その第一歩は

1. 乗法還元の素点 `p` を取る（`v_p(Δ_min) ≠ 0`）
2. 半安定なので `v_p(j) < 0` である
3. `SemistableAt` の右の選択肢が**極小モデル `C • E` と `v_p(c₄) = 0`** を与える（第 954）

であり、これが `isMinimal_baseChange_at_bad_prime`・
`hasMultiplicativeReduction_at_bad_prime`（第 976、証明済み）の入力そのものである。

☆本ブロックはその 3 段を `SSCurve` の語彙で 1 本にまとめる。
★残るのは**分裂性**と **`ζ_l ∈ L_v`** の 2 枚だけになる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★**乗法還元の素点では `v_p(j) < 0`**——★**無条件**（第 1350）。 -/
theorem SSCurve.exists_jExp_neg (E : SSCurve) (h : E.HasMultRed) :
    ∃ p : HeightOneSpectrum (𝓞 E.fld), jExp p E.W < 0 := by
  obtain ⟨p, hp⟩ := h
  refine ⟨p, ?_⟩
  have hmd : minDeltaExp p E.W = max 0 (- jExp p E.W) :=
    minDeltaExp_eq_maxJ_of_semistable p E.W (E.ss p)
  have hlh : E.localHeightAt p = minDeltaExp p E.W := rfl
  rw [hlh] at hp
  by_cases hc : jExp p E.W < 0
  · exact hc
  · exact absurd (by rw [hmd]; omega) hp

/-- ★★★★★★★★★★★★★★★★★★★★
**乗法還元から悪い素点の局所データが出る**——★**無条件**（第 1350）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`HasMultRed` → `v_p(j) < 0`（半安定だから）→ 第 954 の 4 つのデータ。

★★★これが `isMinimal_baseChange_at_bad_prime`・
`hasMultiplicativeReduction_at_bad_prime`（第 976）の入力そのものである。 -/
theorem SSCurve.exists_bad_prime_data (E : SSCurve) (h : E.HasMultRed) :
    ∃ (p : HeightOneSpectrum (𝓞 E.fld)) (C : WeierstrassCurve.VariableChange E.fld),
      jExp p E.W < 0 ∧ IsMinimal (primeSubring p) (C • E.W) ∧
      ∃ hc : (C • E.W).c₄ ≠ 0, valAdd p (Units.mk0 ((C • E.W).c₄) hc) = 0 := by
  obtain ⟨p, hj⟩ := E.exists_jExp_neg h
  obtain ⟨C, hC, hc4⟩ := exists_minimal_c4_unit_of_jExp_neg p E.W (E.ss p) hj
  exact ⟨p, C, hj, hC, hc4⟩

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**乗法還元から完備化の上の乗法還元が出る**——★**無条件**（第 1351）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`hp` は `valuation_algebraMap_adicCompletion`（在庫、証明済み）が与える。
★極小性と乗法還元は第 976（証明済み）が与える。

★★★これで #4 に残るのは**分裂性**と **`ζ_l ∈ L_v`** の 2 枚だけである。 -/
theorem SSCurve.exists_local_multRed (E : SSCurve) (h : E.HasMultRed) :
    ∃ (p : HeightOneSpectrum (𝓞 E.fld)) (C : WeierstrassCurve.VariableChange E.fld),
      jExp p E.W < 0 ∧
      ∃ _ : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers E.fld)
          ((C • E.W).baseChange (p.adicCompletion E.fld)),
        WeierstrassCurve.HasMultiplicativeReduction (p.adicCompletionIntegers E.fld)
          ((C • E.W).baseChange (p.adicCompletion E.fld)) := by
  obtain ⟨p, C, hj, hC, hc4ne, hc4⟩ := E.exists_bad_prime_data h
  have hp := valuation_algebraMap_adicCompletion E.fld p
  haveI hmin : WeierstrassCurve.IsMinimal (p.adicCompletionIntegers E.fld)
      ((C • E.W).baseChange (p.adicCompletion E.fld)) :=
    isMinimal_baseChange_at_bad_prime (Lv := p.adicCompletion E.fld)
      (R := p.adicCompletionIntegers E.fld) p hp E.W C hC hc4ne hc4
  exact ⟨p, C, hj, hmin,
    hasMultiplicativeReduction_at_bad_prime (Lv := p.adicCompletion E.fld)
      (R := p.adicCompletionIntegers E.fld) p hp E.W C hC hc4ne hc4 hj⟩

/-! ## ★出典の紐付け(`.src`) -/

def SSCurve.exists_local_multRed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法還元から完備化の上の乗法還元が出る。★無条件)",
    sectionId := "genell-thm-3-8" }

def SSCurve.exists_local_multRed.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isMinimal_baseChange_at_bad_prime / hasMultiplicativeReduction_at_bad_prime(第 976、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.hasMultiplicativeReduction_at_bad_prime") 1,
    .citation "[ABC3]" "valuation_algebraMap_adicCompletion(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.valuation_algebraMap_adicCompletion") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1351）**——#4 に残るのは" ++
       "**分裂性**と **`ζ_l ∈ L_v`** の 2 枚だけになった。") 2 ]

def SSCurve.exists_jExp_neg.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法還元の素点では v_p(j) < 0。★無条件)",
    sectionId := "genell-thm-3-8" }

def SSCurve.exists_bad_prime_data.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法還元から悪い素点の局所データが出る。★無条件)",
    sectionId := "genell-thm-3-8" }

def SSCurve.exists_bad_prime_data.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_minimal_c4_unit_of_jExp_neg(第 954、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_minimal_c4_unit_of_jExp_neg") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1350）**——`EllModuliWitness` の " ++
       "`imageContainsSL2J_torsionExt`（#4）の**入口**である。" ++
       "☆これで `isMinimal_baseChange_at_bad_prime`・" ++
       "`hasMultiplicativeReduction_at_bad_prime`（第 976）が回り、" ++
       "残るのは**分裂性**と **`ζ_l ∈ L_v`** の 2 枚だけになる。") 2 ]

end ABC3.Found.GenEll
