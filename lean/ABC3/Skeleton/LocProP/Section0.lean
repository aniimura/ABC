import ABC3.Meta.Claim
import ABC3.Found.LocProP.PAdicField
import ABC3.Found.LocProP.PDerivate
import ABC3.Interface.LocProP.EtaleSetup

/-!
# [LocProP] §0 —— Preliminaries and Notations(4 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.14-15。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

| 項目 | 状態 | 場所 |
|---|---|---|
| Definition 0.1(`p`-adic (local) field) | ★済(sorry 無し) | `Found/LocProP/PAdicField.lean` |
| Definition 0.2(`p`-derivate) | ★済(sorry 無し、逸脱記録つき) | `Found/LocProP/PDerivate.lean` |
| Definition 0.3(arithmetic first Chern class) | Interface 待ち | 本ファイル・`Interface/LocProP/EtaleSetup.lean` |
| Lemma 0.4(étale cohomology ≅ 群 cohomology) | 未着手 | — |

★§0 は 2/4(2026-09-04 時点)。
-/

namespace ABC3.Skeleton.LocProP

open ABC3.Interface.LocProP

/-- **[LocProP] Definition 0.3** —— 直線束 `L` の arithmetic first Chern class
`c₁(L) ∈ H²(X_K, ℤ_p(1))`。

原文 (LocProP p.15):
> Definition 0.3. We shall refer to c1(L) as the arithmetic first Chern class of L. -/
def arithmeticFirstChernClass (E : EtaleSetup) (L : E.Pic) : E.H2Zp1 := E.c1 L

def arithmeticFirstChernClass.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 15, item := "Definition 0.3",
    sectionId := "locprop-def-0-3" }

def arithmeticFirstChernClass.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "Kummer 完全列そのもの(0 → ℤ/p^nℤ(1) → 𝔾_m → 𝔾_m → 0)と、それが定める
       H¹(X_K,𝔾_m) → H²(X_K,(ℤ/p^nℤ)(1)) の連結射、および n についての逆極限。
       mathlib にスキームの étale cohomology が無いため
       Interface/LocProP/EtaleSetup.lean の c1 として posit した(証明していない)。" 15 ]

end ABC3.Skeleton.LocProP
