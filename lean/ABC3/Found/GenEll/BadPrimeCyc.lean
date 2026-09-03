/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.BadPrimeFromMultRed
import ABC3.Found.GenEll.LocalCycPackage
import ABC3.Found.GaloisRep.RamifiedLocalData
import ABC3.Meta.Claim

/-!
# 第 1378 ブロック —— **円分拡大の上の悪い素点データ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1351 の**円分拡大版**

第 1351（`SSCurve.exists_local_multRed`）は完備化 `L_p` の上で
極小性と乗法還元を出した。★本ブロックはそれを

> **`L_p(ζ_l)` の上で**——`hpe`（付値は `e` 倍）つきで

出す。☆`ζ_l` が入っているので、第 1372（`exists_h2_h1_of_bad_prime_ram`）の
「原始 `l` 乗根が局所体にある」という要求が満たされる。

★残るのは**分裂性**（剰余体で 2 次式が分裂すること）だけである
——それは第 1035-1036 の二者択一（分裂するか、非分岐 2 次拡大で分裂させるか）が扱う。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField
open ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★★★★★★★★★
**乗法還元から `L_p(ζ_l)` の上の乗法還元が出る**——★**無条件**（第 1378）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`hpe`（付値は `e` 倍）と `l ∤ e` が一緒に出る。
★★★これで第 1372（`exists_h2_h1_of_bad_prime_ram`）の入力のうち
**分裂性以外がすべて揃う**。 -/
theorem SSCurve.exists_local_multRed_cyc (E : SSCurve) (h : E.HasMultRed)
    {l : ℕ} (hl : l.Prime) :
    ∃ (p : HeightOneSpectrum (𝓞 E.fld)) (C : WeierstrassCurve.VariableChange E.fld) (e : ℕ),
      jExp p E.W < 0 ∧ 1 ≤ e ∧ ¬ (l ∣ e) ∧
      (∀ x : E.fld, (HeightOneSpectrum.valuation (locCycField p hl)
            (IsDiscreteValuationRing.maximalIdeal (locCycRing p hl)))
            (algebraMap E.fld (locCycField p hl) x)
          = ((HeightOneSpectrum.valuation E.fld p) x) ^ e) ∧
      ∃ _ : WeierstrassCurve.IsMinimal (locCycRing p hl)
          ((C • E.W).baseChange (locCycField p hl)),
        WeierstrassCurve.HasMultiplicativeReduction (locCycRing p hl)
          ((C • E.W).baseChange (locCycField p hl)) := by
  obtain ⟨p, C, hj, hC, hc4ne, hc4⟩ := E.exists_bad_prime_data h
  obtain ⟨e, he1, hle, hIe⟩ := exists_locCyc_package p hl
  have hpe := valuation_algebraMap_locCycField p hl he1 hIe
  haveI hmin : WeierstrassCurve.IsMinimal (locCycRing p hl)
      ((C • E.W).baseChange (locCycField p hl)) :=
    isMinimal_baseChange_at_bad_prime_ram (Lv := locCycField p hl)
      (R := locCycRing p hl) p hpe E.W C hC hc4ne hc4
  exact ⟨p, C, e, hj, he1, hle, hpe, hmin,
    hasMultiplicativeReduction_at_bad_prime_ram (Lv := locCycField p hl)
      (R := locCycRing p hl) he1 p hpe E.W C hC hc4ne hc4 hj⟩

/-! ## ★出典の紐付け(`.src`) -/

def SSCurve.exists_local_multRed_cyc.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法還元から L_p(ζ_l) の上の乗法還元が出る。★無条件)",
    sectionId := "genell-thm-3-8" }

def SSCurve.exists_local_multRed_cyc.needs : List ProofObligation :=
  [ .citation "[ABC3]" "SSCurve.exists_bad_prime_data(第 1350、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.SSCurve.exists_bad_prime_data") 1,
    .citation "[ABC3]" "exists_locCyc_package(第 1377、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_locCyc_package") 1,
    .citation "[ABC3]" "isMinimal_baseChange_at_bad_prime_ram(第 1375、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isMinimal_baseChange_at_bad_prime_ram") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1378）**——第 1351 の**円分拡大版**である。" ++
       "☆`ζ_l` が局所体に入っているので、第 1372 の要求が満たされる。" ++
       "★残るのは**分裂性**だけ——第 1035-1036 の二者択一が扱う。") 19 ]

end ABC3.Found.GenEll
