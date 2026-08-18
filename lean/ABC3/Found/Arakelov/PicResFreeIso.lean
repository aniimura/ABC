import ABC3.Found.Arakelov.PicOnFree

/-!
# Arakelov (B1) 第 49 ブロック —— **生成元の制限は生成元(加群の段)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★制限と `free` は `rfl` で交換する

    (free F)|_V  =  free ((Over.forget V).op ⋙ F)

★★★どちらも**切断ごと**に `ModuleCat.free` を当てるだけなので、
**定義から等しい**(2026-08-18 実測、`rfl`)。

★これに第 47 ブロック(`(yoneda W)|_V ≅ yoneda (W ⊓ V)`)を継ぐと

    (free (yoneda W))|_V  ≅  free (yoneda (W ⊓ V))

が出る。★★これが Beck–Chevalley の mate の**左端**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `restrict_freeObj` | ★★制限と `free` は `rfl` で交換 |
| `restrictFreeYonedaIso` | ★★★★**`(free (yoneda W))|_V ≅ free (yoneda (W ⊓ V))`** |

## ★★★これで mate の両端が揃った

    左端: (f|)^*_pre ((free (yoneda W))|_V)
            ≅ (f|)^*_pre (free (yoneda (W ⊓ V)))      ★本ブロック
            ≅ free (yoneda (overPost (W ⊓ V)))        ★第 48
    右端: (f^*_pre (free (yoneda W)))|_{f⁻¹V}
            ≅ (free (yoneda f⁻¹W))|_{f⁻¹V}            ★第 24
            ≅ free (yoneda (f⁻¹W ⊓ f⁻¹V))             ★本ブロック

★`f⁻¹(W ⊓ V) = f⁻¹W ⊓ f⁻¹V`(第 23)なので**同じ対象**である。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {Y : Scheme.{u}} (V : Y.Opens)

/-- ★★**制限と `free` は `rfl` で交換する**——どちらも切断ごとだからである。 -/
theorem restrict_freeObj (F : (Y.Opens)ᵒᵖ ⥤ Type u) :
    (restrictPresheafFunctor Y V).obj
        (PresheafOfModules.freeObj
          (R := Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}) F)
      = PresheafOfModules.freeObj
          (R := ((Over.forget V).op ⋙ Y.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})
          ((Over.forget V).op ⋙ F) := rfl

/-- ★★★★**生成元の制限は生成元である**(加群の段)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 47 ブロック(型前層の段)を `free` で送っただけである。
★これが Beck–Chevalley の mate の左端である。 -/
noncomputable def restrictFreeYonedaIso (W : Y.Opens) :
    (restrictPresheafFunctor Y V).obj
        ((PresheafOfModules.free
          (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj (yoneda.obj W))
      ≅ (PresheafOfModules.free
          (((Over.forget V).op ⋙ Y.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj (objOn V W)) :=
  eqToIso (restrict_freeObj V (yoneda.obj W))
    ≪≫ (PresheafOfModules.free
        (((Over.forget V).op ⋙ Y.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})).mapIso
      (yonedaRestrictIso V W)

/-! ## ★出典の紐付け(`.src`) -/

def restrict_freeObj.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限と free が交換すること)",
    sectionId := "genell-def-1-1-i" }

def restrictFreeYonedaIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の制限が生成元であること(加群の段))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
