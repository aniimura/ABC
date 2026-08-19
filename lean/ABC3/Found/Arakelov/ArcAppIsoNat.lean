import ABC3.Found.Arakelov.ArcBridgeSmul

/-!
# Arakelov (C3) 第 267 ブロック —— **`appIso` の `hom` 側の自然性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★橋の自然性に要る最後の材料

`bridgeApp` が前層の射になるには、`e.hom` の自然性(在庫)に加えて
**`appIso` の自然性**が要る。★mathlib が持っているのは `inv` 側だけである:

    X.presheaf.map i ≫ (f.appIso V).inv = (f.appIso U).inv ≫ Y.presheaf.map (opensFunctor i)

★★`hom` 側は同型の相殺 2 回で出る(`appIso_hom_naturality`)。

## ★★制限射は `rfl` で一致する

`Over V` 側の制限射と mathlib の `restrict` の制限射は、どちらも
`F.val.map (opensFunctor.map i).op` で**定義的に等しい**(`restrictPre_map`)。
★★★第 255(切断が `rfl`)に続き、**射も `rfl`** である
——構造だけが違い、データは同じという状況が確定した。

| 定理 | 内容 |
|---|---|
| `overMap` | ★`V.Opens` の射 → `Over V` の射 |
| `restrictPre_map` | ★★制限射は `rfl` で一致 |
| `appIso_hom_naturality` | ★★★`appIso` の `hom` 側の自然性 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★`Over V` の射。 -/
def overMap {W W' : V.toScheme.Opens} (i : W' ⟶ W) : overObj V W' ⟶ overObj V W :=
  Over.homMk (homOfLE (Scheme.Hom.image_mono _ (leOfHom i)))

/-- ★★`Over V` 側の前層の制限射。 -/
theorem restrictPre_map {W W' : V.toScheme.Opens} (i : W' ⟶ W)
    (x : ((((restrictPresheafFunctor X V).obj F.val).obj (op (overObj V W))) : Type)) :
    (((restrictPresheafFunctor X V).obj F.val).map (overMap V i).op).hom x
      = (F.val.map (Scheme.Hom.opensFunctor V.ι |>.map i).op).hom x :=
  rfl


/-- ★★`appIso` の `hom` 版の自然性。 -/
theorem appIso_hom_naturality {W W' : V.toScheme.Opens} (i : op W ⟶ op W') :
    (V.ι.appIso W).hom ≫ V.toScheme.presheaf.map i
      = X.presheaf.map ((Scheme.Hom.opensFunctor V.ι).op.map i) ≫ (V.ι.appIso W').hom := by
  have h := Scheme.Hom.appIso_inv_naturality (f := V.ι) i
  calc (V.ι.appIso W).hom ≫ (V.toScheme).presheaf.map i
      = (V.ι.appIso W).hom ≫ ((V.toScheme).presheaf.map i ≫ (V.ι.appIso W').inv)
          ≫ (V.ι.appIso W').hom := by
        rw [Category.assoc, Iso.inv_hom_id, Category.comp_id]
    _ = (V.ι.appIso W).hom ≫ ((V.ι.appIso W).inv ≫
          X.presheaf.map ((Scheme.Hom.opensFunctor V.ι).op.map i)) ≫ (V.ι.appIso W').hom := by
        rw [h]
        rfl
    _ = X.presheaf.map ((Scheme.Hom.opensFunctor V.ι).op.map i) ≫ (V.ι.appIso W').hom := by
        rw [← Category.assoc, ← Category.assoc, Iso.hom_inv_id, Category.id_comp]
        rfl



/-! ## ★出典の紐付け(`.src`) -/

def appIso_hom_naturality.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——appIso の hom 側の自然性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
