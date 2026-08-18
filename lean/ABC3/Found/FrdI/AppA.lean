/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop113

/-!
# [FrdI] Appendix: Slim Exponentiation

原文 (FrdI p.118):
> In the present Appendix, we discuss some elementary general nonsense con-

★★`Theorem 3.4, (v)` の `Ψ_Base : 𝒟₁ → 𝒟₂` の構成が
`Definition A.1` と `Proposition A.2` を使うので、ここを先に取る。

## ★測ったこと(2026-08-19)

★★★**`Definition A.1, (i)` の `2-slim` は、我々の `IsRigidFunctor` そのもの**だった ——

原文 (FrdI p.118):
> in the 2-category has no nontrivial automorphisms.

つまり「2-圏の 1-射(＝関手)が非自明な自己同型を持たない」。
`IsRigidFunctor F := ∀ η : F ≅ F, η = Iso.refl F`(`Prop113.lean`)と同じである。

★`Proposition A.2` が扱う 2-圏は**具体的**で、
対象はスライス `𝒞_A`、1-射は `f : A ⟶ B` が誘導する `f_! : 𝒞_A → 𝒞_B`、
2-射はそれらの間の同型である。★mathlib では `f_! = Over.map f`。
-/

universe v u

namespace ABC3.Found.FrdI

open CategoryTheory

/-! ## ★`Definition A.1, (i)` —— スライスの 2-圏が `2-slim` であること -/

/-- ★★★★★**[FrdI] Proposition A.2 の第 1 主張** —— スライスの 2-圏は `2-slim`。

原文 (FrdI p.118):
> isomorphisms between these functors [cf. §0]. Then D is 2-slim. Moreover, the

★★`𝒞` が slim なら、`f : A ⟶ B` が誘導する `f_! = Over.map f` は rigid。

★手筋: `Over.map f ⋙ Over.forget B = Over.forget A` は**定義的に等しい**
(どちらも `X ↦ X.left`、`u ↦ u.left`)。
`η` を `Over.forget B` で whisker すれば `Over.forget A` の自己同型になり、
slim 性でそれが恒等。★`Over` の射は `.left` で決まるので、そこから `η = 𝟙`。 -/
theorem isRigidFunctor_overMap {C₁ : Type u} [Category.{v} C₁] (hslim : IsSlimCat C₁)
    {A B : C₁} (f : A ⟶ B) : IsRigidFunctor (Over.map f) := by
  intro η
  -- ★whisker して `Over.forget A` の自己同型にし、slim 性を当てる。
  have hw : (Functor.isoWhiskerRight η (Over.forget B) : Over.forget A ≅ Over.forget A)
      = Iso.refl _ := hslim A _
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have hleft : (η.hom.app X).left = 𝟙 X.left :=
    congrArg (fun t : Over.forget A ≅ Over.forget A => t.hom.app X) hw
  exact Over.OverMorphism.ext hleft

def isRigidFunctor_overMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 118,
    item := "Proposition A.2 — スライスの 2-圏は 2-slim",
    sectionId := "frdi-app-a" }

/-- ★★★★★**[FrdI] Proposition A.2 の第 2 主張** ——
`E : 𝒞 → |D|`(`A ↦ 𝒞_A`、`f ↦ f_!`)は**忠実**。

原文 (FrdI p.118):
> determines an equivalence of categories C

★★`E` は対象でも射でも構成から全射なので、圏同値であることの中身は
**「`f_! ≅ g_!` なら `f = g`」**、すなわちこの補題である。

★手筋: whisker と slim 性で成分の `.left` が恒等になり、
`Over` の三角形 (`Over.w`) がそのまま `f = g` を与える。 -/
theorem overMap_injective {C₁ : Type u} [Category.{v} C₁] (hslim : IsSlimCat C₁)
    {A B : C₁} {f g : A ⟶ B} (ε : Over.map f ≅ Over.map g) : f = g := by
  have hw : (Functor.isoWhiskerRight ε (Over.forget B) : Over.forget A ≅ Over.forget A)
      = Iso.refl _ := hslim A _
  have hleft : (ε.hom.app (Over.mk (𝟙 A))).left = 𝟙 A :=
    congrArg (fun t : Over.forget A ≅ Over.forget A => t.hom.app (Over.mk (𝟙 A))) hw
  have htri := Over.w (ε.hom.app (Over.mk (𝟙 A)))
  rw [hleft] at htri
  simpa using htri.symm

def overMap_injective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 118,
    item := "Proposition A.2 — slim exponentiation は忠実",
    sectionId := "frdi-app-a" }

end ABC3.Found.FrdI
