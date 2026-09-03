/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ElementaryFrobenioid
import Mathlib.CategoryTheory.Endomorphism
import Mathlib.CategoryTheory.Widesubcategory
import Mathlib.CategoryTheory.MorphismProperty.Composition
import Mathlib.Data.Nat.Prime.Basic
import ABC3.Meta.Claim
import ABC3.Found.FrdI.MorphismTypes.Definition12

/-!
# MorphismTypes —— `[FrdI] Remark 1.2.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite
universe v u w u2 v2
variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### ★Remark 1.2.1 —— 原文が「定義から形式的に従う」と述べる含意

原文 (FrdI p.24):
> Remark 1.2.1. The following implications follow formally from the definitions:

原文 (FrdI p.24):
> pull-back morphism which is a base-isomorphism ⇐⇒ isomorphism

原文 (FrdI p.24):
> base-trivial =⇒metrically trivial

原文 (FrdI p.24):
> base-identity =⇒ Div-identity

原文 (FrdI p.24):
> universally Div-Frobenius-trivial =⇒ Div-Frobenius-trivial

★原文は「follow formally」と書くだけで証明を置かない。**それを実際に示す。**

★**含意は4本である**(2026-08-15 訂正)。一度この docstring は3本しか引用しておらず、
2本目(`base-trivial ⟹ metrically trivial`)を写し落としていた。
下の `isMetricallyTrivial_of_isBaseTrivial` がその4本目である。
-/

/-- ★**`base-identity ⟹ Div-identity`**(Remark 1.2.1 の第3行)。

`Base(φ) = Base(𝟙)` なら、`Φ` を当てた `Φ(φ) = Φ(𝟙)` も従う。 -/
theorem isDivIdentity_of_isBaseIdentity {A : C} {φ : A ⟶ A}
    (h : IsBaseIdentity P φ) : IsDivIdentity P φ := by
  show Φ.map (P.Base φ) = Φ.map (P.Base (𝟙 A))
  rw [show P.Base φ = P.Base (𝟙 A) from h]

/-! ### ★出典の紐付け(`.src`) -/

def isDivIdentity_of_isBaseIdentity.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Remark 1.2.1",
    sectionId := "frdi-remark-1-2-1" }

end ABC3.Found.FrdI
