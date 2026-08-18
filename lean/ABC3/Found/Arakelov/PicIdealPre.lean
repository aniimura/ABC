import ABC3.Found.Arakelov.PicIdealSec

/-!
# Arakelov (B2) 第 149 ブロック —— **イデアル層の前層加群**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★第 148 の切断を `PresheafOfModules` に組む

`PresheafOfModules R` の `map` は

    obj U ⟶ (ModuleCat.restrictScalars (R.map f).hom).obj (obj V)

という**係数制限つき**の射である。★これが素直に書けない
——`ModuleCat.ofHom` を書くと制限前の対象で instance を探してしまう。

## ★★逃げ道——`letI` + `Module.compHom`

    letI : Module Γ(X, U) (idealSections D V) := Module.compHom _ (X.presheaf.map f).hom

★これで `ModuleCat.ofHom` が通り、結果は `restrictScalars` の対象と**定義から一致**する。
★★[[ring-instance-two-paths]] の 8 例目。

## ★★★`map_id` / `map_comp` は `cat_disch` が自動で片付けた

★構造の既定値がそのまま通る——**書かなくてよかった**。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

/-- ★イデアル層の前層加群。 -/
noncomputable def idealPresheaf : X.PresheafOfModules where
  obj U := ModuleCat.of (Γ(X, U.unop) : Type u) (idealSections D U.unop)
  map {U V} f :=
    letI : Module (Γ(X, U.unop) : Type u) (idealSections D V.unop) :=
      Module.compHom _ (X.presheaf.map f).hom
    ModuleCat.ofHom
    { toFun := fun s => ⟨(X.presheaf.map f).hom s.1,
        idealSections_res D (leOfHom f.unop) s.1 s.2⟩
      map_add' := fun a b => Subtype.ext (map_add _ _ _)
      map_smul' := fun c a => Subtype.ext (map_mul _ _ _) }


/-! ## ★出典の紐付け(`.src`) -/

def idealPresheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——イデアル層の前層加群)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
