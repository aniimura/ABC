import ABC3.Meta.Claim
import ABC3.Found.LocProP.PAdicField
import ABC3.Found.LocProP.PDerivate
import ABC3.Interface.LocProP.EtaleSetup
import ABC3.Found.LocProP.CohomologyComparison

/-!
# [LocProP] §0 —— Preliminaries and Notations(4 項目)

原典: S. Mochizuki, *The Local Pro-p Anabelian Geometry of Curves* [LocProP]、
物理 p.14-15。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

| 項目 | 状態 | 場所 |
|---|---|---|
| Definition 0.1(`p`-adic (local) field) | ★済(sorry 無し) | `Found/LocProP/PAdicField.lean` |
| Definition 0.2(`p`-derivate) | ★済(sorry 無し、逸脱記録つき) | `Found/LocProP/PDerivate.lean` |
| Definition 0.3(arithmetic first Chern class) | ★済(Interface posit 経由) | 本ファイル・`Interface/LocProP/EtaleSetup.lean` |
| Lemma 0.4(étale cohomology ≅ 群 cohomology) | ★済(Interface posit 経由、★逸脱記録) | 本ファイル・`Interface/LocProP/EtaleSetup.lean` |

★§0 は 4/4(2026-09-04 時点)。Definition 0.3・Lemma 0.4 はどちらも
étale cohomology 自体が mathlib に無いための Interface posit(証明していない、
`.needs` に明記)。
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

/-- **[LocProP] Lemma 0.4** —— `Δ_X`(・`Π_X`)の群 cohomology と `X_K̄`(・`X_K`)の
étale cohomology を比べる自然な射は同型である。

原文 (LocProP p.15):
> Lemma 0.4. Assume that K is a field. For all integers i, r, the natural morphisms

実装は `Found/LocProP/CohomologyComparison.lean` へ委譲する。 -/
theorem etaleGroupCohomologyComparison_bijective (E : CohomologyComparisonSetup)
    (i r : ℤ) : Function.Bijective (E.compare i r) :=
  ABC3.Found.LocProP.cohomologyComparison_bijective E i r

def etaleGroupCohomologyComparison_bijective.src : ABC3.Meta.Source :=
  { paper := "LocProP", pdfPage := 15, item := "Lemma 0.4",
    sectionId := "locprop-lemma-0-4" }

/-- G9 の非空虚性対照(`Found/LocProP/CohomologyComparison.lean` の具体例へ委譲)。 -/
theorem etaleGroupCohomologyComparison_bijective.nonvacuous :
    Function.Bijective (CohomologyComparisonSetup.nonvacuous.some.compare 0 0) :=
  ABC3.Found.LocProP.cohomologyComparison_bijective.nonvacuous

def etaleGroupCohomologyComparison_bijective.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      "★★★逸脱: 原文の証明は Leray-Serre spectral sequence(Δ_X ⊴ Π_X、商 Γ_K)・
       Künneth の公式・étale 被覆と H¹ の対応の三段を使うが、mathlib には
       profinite 群の Lyndon–Hochschild–Serre spectral sequence も
       スキームの étale cohomology も無い(2026-09-04 実測)。
       比較射の全単射性そのものを Interface/LocProP/EtaleSetup.lean の
       compare_bijective として posit した(証明していない)。
       これは Definition 0.3 の c1 と同じ扱い方針である。" 16 ]

end ABC3.Skeleton.LocProP
