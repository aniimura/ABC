/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Equiv

/-!
# [FrdI] `Proposition 4.4, (ii)` —— `𝒞^birat` は **group-like 型**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.83。

原文 (FrdI p.83):
> (ii) The functor Cbirat →F0D of (i) determines a structure of Frobenioid

★原文は「**of group-like type**」と言う。★★その中身は
**`0_𝒟` の値がすべて 1 元単系である**ことに尽きる ——
`Div` が情報を持たないので、`Φ^birat` の値は群(実は自明群)になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-- ★★★**[FrdI] Proposition 4.4, (ii)** —— `𝒞^birat` は **group-like 型**。

★`biratPre` の単系は `0_𝒟`(値がすべて 1 元)なので、
すべての元が `0` に等しく、したがって加法的可逆である。 -/
theorem birat_isOfGroupLikeType : IsOfGroupLikeType (biratPre P G) := by
  intro A
  refine (isGroupLike_iff _).mpr (fun a => ?_)
  rw [Subsingleton.elim a 0]
  exact isAddUnit_zero

/-- ★locator —— `Proposition 4.4, (ii)` の「group-like 型」。 -/
def birat_isOfGroupLikeType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (ii) — group-like 型",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
