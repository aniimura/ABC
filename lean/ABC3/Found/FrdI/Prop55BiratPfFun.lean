/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55BiratPf

/-!
# [FrdI] Proposition 5.5, (ii) の birat の側を**関手**として組む

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★★`Prop55BiratPf.lean` は **射の全単射**

  `Hom_{(𝒞^birat)^pf}(A,B) ≃ Hom_{(𝒞^pf)^birat}(A,B)`

(`biratPfHomEquiv`)と、**恒等射の保存**(`biratPfHom_id`)・
**合成との両立**(`biratPfHom_comp`)まで作ってある。
★本ファイルはそれを**関手**に束ね、**充満忠実**であることを言う。

## ★★対象について(記録)

`(𝒞^birat)^pf` の**根 1 の部分**(`PfCat (biratPre P G) F'`)から
`(𝒞^pf)^birat`(`BiratCat (pfRootPre P F) Gpf`、対象は対 `(A,n)`)への関手になる。
★像は根 1 の対象だけなので**本質的全射ではない** ——
原文の「equivalence of categories」に届かせるには
根 `n` を一般にした版(`scaleRootEquiv` で根を揃える段)が要る。
★un-tr の側(`Prop55UntrFun.lean`)では両側の対象が一致していたので
そこは無料だったが、birat の側はそうではない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BiratPfFun

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  {G : Frobenioid P}

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)の関手** ——

  `(𝒞^birat)^pf` の根 1 の部分 ⥤ `(𝒞^pf)^birat`

★対象は `A ↦ ⟨A, 1⟩`、射は `biratPfHom`。
★恒等・合成は `biratPfHom_id` / `biratPfHom_comp` そのもの。 -/
noncomputable def biratPfFunctor (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    PfCat (biratPre P G) F' ⥤ BiratCat (pfRootPre P F) Gpf where
  obj A := (⟨biratDown P G A, 1⟩ : PfRootObj P F)
  map {A B} f := biratPfHom hfi Gpf F' (biratDown P G A) (biratDown P G B) f
  map_id A := biratPfHom_id hfi Gpf F' (biratDown P G A)
  map_comp {A B E} f g := biratPfHom_comp hfi Gpf F'
    (biratDown P G A) (biratDown P G B) (biratDown P G E) f g

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** —— 関手は**忠実**。 -/
theorem biratPfFunctor_faithful (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    (biratPfFunctor hfi Gpf F').Faithful where
  map_injective {A B} {f g} h :=
    (biratPfHom_bijective hfi hiso Gpf F'
      (biratDown P G A) (biratDown P G B)).1 h

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** —— 関手は**充満**。 -/
theorem biratPfFunctor_full (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) :
    (biratPfFunctor hfi Gpf F').Full where
  map_surjective {A B} h :=
    (biratPfHom_bijective hfi hiso Gpf F'
      (biratDown P G A) (biratDown P G B)).2 h

/-- ★★★★★**[FrdI] Proposition 5.5, (ii)(birat の側)** ——
関手が根 1 の対象のあいだで与える全単射は `biratPfHomEquiv` そのもの。 -/
theorem biratPfFunctor_map_eq (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    {A B : PfCat (biratPre P G) F'} (f : A ⟶ B) :
    (biratPfFunctor hfi Gpf F').map f
      = biratPfHomEquiv hfi hiso Gpf F' (biratDown P G A) (biratDown P G B) f :=
  rfl

end BiratPfFun

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Proposition 5.5, (ii)` の birat の側を関手として組んだもの
(充満忠実まで。本質的全射は根を一般にする段が残る)。 -/
def biratPfFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (ii) — (𝒞^birat)^pf ⥤ (𝒞^pf)^birat(充満忠実)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
