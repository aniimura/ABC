/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaSplitDichotomy
import ABC3.Found.GenEll.BadPrimeCyc
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.Lemma37Full
import ABC3.Meta.Claim

/-!
# 第 1384 ブロック —— **悪い素点の惰性は幂単かつ非自明**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——α 側の到達点

第 1352 で節点にした `exists_h2_h1_unipotent` を、第 1365-1383 の道で閉じる。

☆組み立て:

| 段 | 第 |
|---|---|
| 悪い素点 `p`・極小化 `C`・`v_p(c₄) = 0` | 1350 |
| `L_p(ζ_l)` の完備 DVR パッケージ（`e`・`l ∤ e`・`hpe`） | 1377 |
| 拡大の上の極小性・乗法還元 | 1375 |
| 分裂性の二者択一 → `h2`・`h1` | 1383 |

★`hcop`（`¬ l ∣ v_p(j)`）は `PrimeToLocalHeights` から出る
——半安定なので `localHeight = −jExp` である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine IsDedekindDomain NumberField
open ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 悪い素点の惰性は `T_l E` の上で `mod l` で幂単かつ非自明**
——★**無条件**（第 1384）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `imageContainsSL2J_torsionExt`（#4）に残っていた節点である。 -/
theorem SSCurve.exists_h2_h1_unipotent_of_multRed (E : SSCurve) (l : ℕ) [hlf : Fact l.Prime]
    (hm : E.HasMultRed) (hpr : E.PrimeToLocalHeights l) :
    ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
      (∀ x : E.tate l, ∃ u : E.tate l,
          galTate E.W l σ (galTate E.W l σ x) + x
            = galTate E.W l σ x + galTate E.W l σ x + l • u) ∧
        (∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u) := by
  have hl : l.Prime := hlf.out
  obtain ⟨p, C, hj, hC, hc4ne, hc4⟩ := E.exists_bad_prime_data hm
  have hcop := E.not_dvd_jExp_of_primeToLocalHeights hl hpr p hj
  obtain ⟨e, he1, hle, hIe⟩ := exists_locCyc_package p hl
  have hpe := valuation_algebraMap_locCycField p hl he1 hIe
  haveI hCE : (C • E.W).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.W.isUnit_Δ
  haveI hCZ : CharZero (locCycRing p hl) :=
    RingHom.charZero (algebraMap (locCycRing p hl) (locCycField p hl))
  haveI hmin : WeierstrassCurve.IsMinimal (locCycRing p hl)
      ((C • E.W).baseChange (locCycField p hl)) :=
    isMinimal_baseChange_at_bad_prime_ram (Lv := locCycField p hl)
      (R := locCycRing p hl) p hpe E.W C hC hc4ne hc4
  haveI hmult : WeierstrassCurve.HasMultiplicativeReduction (locCycRing p hl)
      ((C • E.W).baseChange (locCycField p hl)) :=
    hasMultiplicativeReduction_at_bad_prime_ram (Lv := locCycField p hl)
      (R := locCycRing p hl) he1 p hpe E.W C hC hc4ne hc4 hj
  have hζ0 := isPrimitiveRoot_locCycRoot p hl
  have hu : IsUnit (cycRoot (p.adicCompletion E.fld) hl) := hζ0.isUnit hl.ne_zero
  have hζv : IsPrimitiveRoot ((hu.unit : (locCycField p hl)ˣ) : locCycField p hl) l := by
    rw [IsUnit.unit_spec]
    exact hζ0
  exact exists_h2_h1_of_multRed_local E.fld he1 p hpe E.W C hC hc4ne hc4 hj hle hcop hζv

/-! ## ★出典の紐付け(`.src`) -/

def SSCurve.exists_h2_h1_unipotent_of_multRed.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(悪い素点の惰性は T_l E の上で mod l で幂単かつ非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def SSCurve.exists_h2_h1_unipotent_of_multRed.needs : List ProofObligation :=
  [ .citation "[ABC3]" "SSCurve.exists_bad_prime_data(第 1350、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.SSCurve.exists_bad_prime_data") 1,
    .citation "[ABC3]" "exists_locCyc_package(第 1377、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_locCyc_package") 1,
    .citation "[ABC3]" "exists_h2_h1_of_multRed_local(第 1383、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_of_multRed_local") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1384）**——第 1352 で節点にした `exists_h2_h1_unipotent` が" ++
       "第 1365-1383 の道で閉じた。" ++
       "☆`hcop` は `PrimeToLocalHeights` から出る——半安定なので `localHeight = −jExp`。") 19 ]

end ABC3.Found.GenEll
