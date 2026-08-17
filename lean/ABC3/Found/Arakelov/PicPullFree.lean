import ABC3.Found.Arakelov.PicDeltaColimit

/-!
# Arakelov (B1) 第 24 ブロック —— **引き戻しの生成元上の具体形**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★前層の引き戻しは「抽象」だが、生成元の上では書ける

mathlib の `PresheafOfModules.pullback` は**随伴関手定理で作られている**——
対象の具体的な記述は無い。★★しかし生成元 `free (yoneda V)` の上では話が違う:

    pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ V
      : (pushforward φ ⋙ coyoneda.obj (op (free (yoneda V)))).CorepresentableBy
          (free (yoneda (F.obj V)))

★★★**同じ関手を引き戻しも余表現する**(随伴だから)。
`CorepresentableBy.uniqueUpToIso` で**同型が出る**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackCorepresentableBy` | ★★随伴から余表現可能性 |
| `pullbackFreeYonedaIso` | ★★★★★**`f^*(free (yoneda V)) ≅ free (yoneda (f⁻¹V))`** |

★★★これが `δ` の同型性の**第 3 段**(生成元の上で同型)の半分である。
残る半分は `free (yoneda V) ⊗ free (yoneda W) ≅ free (yoneda (V ⊓ W))` であり、
第 23 ブロックの `opensMap_inf`(逆像は `⊓` を保つ)と噛み合う。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★スキームの射から定まる環前層の射(第 20 ブロックの `alphaR`)。 -/
noncomputable abbrev pullbackPhi := alphaR (Opens.map f.base) X.presheaf Y.presheaf f.c

/-! ## ★★随伴からの余表現可能性 -/

/-- ★★**引き戻しは押し出しとの hom 集合を余表現する**——随伴そのものである。

★これを mathlib 側の余表現と突き合わせるのが本ブロックの手口である。 -/
noncomputable def pullbackCorepresentableBy (M : Y.PresheafOfModules) :
    (PresheafOfModules.pushforward (pullbackPhi f) ⋙ coyoneda.obj (op M)).CorepresentableBy
      ((pullbackPre f).obj M) where
  homEquiv := (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv _ _
  homEquiv_comp g h :=
    (PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).homEquiv_naturality_right h g

/-! ## ★★★★★生成元の上での具体形 -/

/-- ★★★★★**生成元の引き戻しは生成元である** `f^*(free (yoneda V)) ≅ free (yoneda (f⁻¹V))`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**抽象に作られた関手の、具体的な値が 1 つ確定した。**
★機構は「同じ関手を余表現する対象は同型」——`CorepresentableBy.uniqueUpToIso` である。

★★これが `δ` の同型性の第 3 段(生成元の上で同型)の半分を与える。 -/
noncomputable def pullbackFreeYonedaIso (V : Y.Opens) :
    (pullbackPre f).obj ((PresheafOfModules.free
        (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj (yoneda.obj V))
      ≅ (PresheafOfModules.free
          (X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
          (yoneda.obj ((Opens.map f.base).obj V)) :=
  (pullbackCorepresentableBy f _).uniqueUpToIso
    (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy (pullbackPhi f) V)

/-! ## ★出典の紐付け(`.src`) -/

def pullbackFreeYonedaIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元の引き戻しが生成元になること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
