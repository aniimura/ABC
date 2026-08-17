import ABC3.Found.Arakelov.PicYonedaInf

/-!
# Arakelov (B1) 第 27 ブロック —— **生成元の上で引き戻しはテンソルを保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★4 本を繋ぐ

第 23–26 ブロックで揃った 4 本を、生成元の上で 1 本に繋ぐ:

    f^*(free(yV) ⊗ free(yW))
      ≅ f^*(free(y(V ⊓ W)))              ★第 26(米田と交わり)
      ≅ free(y(f⁻¹(V ⊓ W)))              ★第 24(生成元の引き戻し)
      = free(y(f⁻¹V ⊓ f⁻¹W))             ★第 23(逆像は ⊓ を保つ、`rfl`)
      ≅ free(y f⁻¹V) ⊗ free(y f⁻¹W)      ★第 26(逆向き)
      ≅ f^*free(yV) ⊗ f^*free(yW)        ★第 24(逆向き、両側)

★★★**これが `δ` の同型性の「生成元の段」である。**

## ★★残り

★余極限による持ち上げ——第 23 ブロックの余極限保存 3 本と
mathlib の `isColimitFreeYonedaCoproductsCokernelCofork` を使う。

★★★`PicSheaf` は**同型による商**なので、`pullback_mul` に要るのは
**対象の同型**だけである——`δ` そのものである必要はない。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★★★★★生成元の上での対象の同型 -/

/-- ★★★★★★**生成元の上で引き戻しはテンソルを保つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 23–26 ブロックの 4 本を繋いだものである。
★★★★これが `Pic` の引き戻しが群準同型になることの**核**である。 -/
noncomputable def pullbackTensorFreeIso (V W : Y.Opens) :
    (pullbackPre f).obj
        (PresheafOfModules.freeObj
            (R := Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (yoneda.obj V)
          ⊗ PresheafOfModules.freeObj
            (R := Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (yoneda.obj W))
      ≅ (pullbackPre f).obj (PresheafOfModules.freeObj
            (R := Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (yoneda.obj V))
        ⊗ (pullbackPre f).obj (PresheafOfModules.freeObj
            (R := Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}) (yoneda.obj W)) :=
  (pullbackPre f).mapIso (freeYonedaInfIso (R := Y.presheaf) V W)
    ≪≫ pullbackFreeYonedaIso f (V ⊓ W)
    ≪≫ (freeYonedaInfIso (R := X.presheaf)
          ((Opens.map f.base).obj V) ((Opens.map f.base).obj W)).symm
    ≪≫ ((pullbackFreeYonedaIso f V).symm ⊗ᵢ (pullbackFreeYonedaIso f W).symm)

/-! ## ★出典の紐付け(`.src`) -/

def pullbackTensorFreeIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の上で引き戻しがテンソルを保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
