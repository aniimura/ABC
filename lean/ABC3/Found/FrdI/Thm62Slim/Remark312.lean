/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Profinite
import ABC3.Found.FrdI.Remark312
import ABC3.Found.FrdI.Prop113
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Thm62Slim.Theorem62

/-!
# Thm62Slim —— `[FrdI] Remark 3.1.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory ABC3.Meta
universe v u w w2 v2 u2

/-- ★★**slim なら Frobenius-slim** の、この節点での確認 ——
`Z_G(H_A)` がすべて自明なら、副有限性を経由しなくても Frobenius-slim。 -/
theorem isFrobeniusSlim_of_isSlimCat {E : Type u} [Category.{v} E]
    (hslim : IsSlimCat E) : IsFrobeniusSlim E :=
  isFrobeniusSlim_of_subsingleton E (fun A => ⟨fun η₁ η₂ => by
    rw [hslim A η₁, hslim A η₂]⟩)

def isFrobeniusSlim_of_isSlimCat.src : Source :=
  { paper := "FrdI", pdfPage := 58,
    item := "Remark 3.1.2 — slim なら Frobenius-slim",
    sectionId := "frdi-remark-3-1-2" }

def isFrobeniusSlim_of_isSlimCat.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isFrobeniusSlim_of_subsingleton"
      (.inProject "ABC3" "ABC3.Found.FrdI.isFrobeniusSlim_of_subsingleton") 58 ]

end ABC3.Found.FrdI
