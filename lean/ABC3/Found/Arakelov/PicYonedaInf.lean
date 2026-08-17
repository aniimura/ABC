import ABC3.Found.Arakelov.PicFreeTensor

/-!
# Arakelov (B1) 第 26 ブロック —— **前順序では米田がテンソルを交わりに送る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★前順序だから効く

一般の圏では `yoneda V ⊗ yoneda W`(点ごとの直積)は表現可能とは限らない。
★★**前順序では表現可能である**——

    Hom(U, V) × Hom(U, W) ≃ Hom(U, V ⊓ W)

は「`U ≤ V` かつ `U ≤ W` ⟺ `U ≤ V ⊓ W`」そのものだからである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `yonedaInfIso` | ★★★★**`yoneda V ⊗ yoneda W ≅ yoneda (V ⊓ W)`** |
| `freeYonedaInfIso` | ★★★★★**`free (yoneda V) ⊗ free (yoneda W) ≅ free (yoneda (V ⊓ W))`** |

★★★これで生成元の段の 3 本目が閉じる:

    第 23: 逆像は ⊓ を保つ                      `opensMap_inf`
    第 24: f^*(free (yoneda V)) ≅ free (yoneda (f⁻¹V))
    第 25: free (F ⊗ G) ≅ free F ⊗ free G
    第 26: yoneda V ⊗ yoneda W ≅ yoneda (V ⊓ W)  ← 本ブロック
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite Limits

variable {α : Type u} [Lattice α] {R : αᵒᵖ ⥤ CommRingCat.{u}}

/-! ## ★★★米田と交わり -/

/-- ★★★★**前順序では `yoneda V ⊗ yoneda W ≅ yoneda (V ⊓ W)`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★各点で `Hom(U,V) × Hom(U,W) ≃ Hom(U, V ⊓ W)` である——
射は高々 1 本なので、成り立つかどうかだけを見ればよい。 -/
def yonedaInfIso (V W : α) :
    (yoneda.obj V ⊗ yoneda.obj W : αᵒᵖ ⥤ Type u) ≅ yoneda.obj (V ⊓ W) :=
  NatIso.ofComponents
    (fun U => Equiv.toIso
      { toFun := fun p => homOfLE (le_inf (leOfHom p.1) (leOfHom p.2))
        invFun := fun p => (homOfLE (le_trans (leOfHom p) inf_le_left),
          homOfLE (le_trans (leOfHom p) inf_le_right))
        left_inv := fun _ => Prod.ext
          (Subsingleton.elim (α := ((unop U) ⟶ V)) _ _)
          (Subsingleton.elim (α := ((unop U) ⟶ W)) _ _)
        right_inv := fun _ => Subsingleton.elim (α := ((unop U) ⟶ V ⊓ W)) _ _ })
    (by
      intro U U' φ
      ext p
      exact Subsingleton.elim (α := ((unop U') ⟶ V ⊓ W)) _ _)

/-- ★★★★★**自由前層加群の言葉で** `free (yoneda V) ⊗ free (yoneda W) ≅ free (yoneda (V ⊓ W))`。

★★★これが `δ` の同型性の第 3 段(生成元の上で同型)の**対象の側**を閉じる。 -/
noncomputable def freeYonedaInfIso (V W : α) :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj V)
      ⊗ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj W)
      ≅ PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj (V ⊓ W)) :=
  (freeTensorIso (R := R) (yoneda.obj V) (yoneda.obj W)).symm
    ≪≫ (PresheafOfModules.free (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})).mapIso
        (yonedaInfIso V W)

/-! ## ★出典の紐付け(`.src`) -/

def yonedaInfIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——前順序で米田がテンソルを交わりに送ること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
