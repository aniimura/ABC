/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Frob
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 5.5, (ii) の `birat` の側 —— 添字圏の出発点

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★測って分かった構造

`(𝒞^pf)^birat` と `(𝒞^birat)^pf` はどちらも**二重の余極限**である:

* 左 = `colim_{Z ∈ IdxBirat(𝒞^pf)(A)} Hom^pf(Z, B)`
* 右 = `colim_{W ∈ IdxPf(𝒞^birat)(A,B)} Hom^birat(W)`

★内側の添字圏が外側の対象に依るので**単純な `colim_{I×J}` にはならない**。

## ★★★本ファイルの一歩 —— 左の添字圏が簡単になること

★`𝒞^pf` では**すべての射が co-angular**(`pfRoot_isCoAngular`、`Proposition 1.4, (i)`)
なので、`IdxBirat(𝒞^pf)(A)` の定義にある「co-angular pre-step」は
**「pre-step」だけ**になる。
★さらに `Proposition 3.2, (ii)` の判定(`isPreStep_mk_iff`)により、
それは `𝒞` の pre-step で代表される。

★★**対象は両側とも `𝒞` の対象である**(`PfRootObj` は根つきだが `BiratCat _ _ = C`)ので、
残るのは射の集合の全単射だけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★**`𝒞^pf` では「co-angular pre-step」＝「pre-step」**
(すべての射が co-angular だから)。

★★これで `(𝒞^pf)^birat` の添字圏(`IdxBirat`)は
**`𝒞^pf` の pre-step のスライスそのもの**になる。 -/
theorem coaPreProp_pfRoot_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) : coaPreProp (pfRootPre P F) f ↔ IsPreStep (pfRootPre P F) f :=
  ⟨fun h => h.2, fun h => ⟨pfRoot_isCoAngular hfi f, h⟩⟩

/-- ★★★**`𝒞^pf` の co-angular pre-step は `𝒞` の pre-step で代表される**。 -/
theorem coaPreProp_pfRoot_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    coaPreProp (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ) ↔ IsPreStep P φ :=
  (coaPreProp_pfRoot_iff hfi _).trans (isPreStep_mk_iff Z φ)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Proposition 5.5, (ii)` の `birat` の側の添字圏。 -/
def coaPreProp_pfRoot_mk_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — (𝒞^pf)^birat の添字圏は 𝒞^pf の pre-step のスライス",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
