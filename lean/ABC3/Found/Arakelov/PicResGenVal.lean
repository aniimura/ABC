import ABC3.Found.Arakelov.PicGenConsistent

/-!
# Arakelov (B1) 第 52 ブロック —— **制限の同型の生成元での値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 37 ブロックと同じ発想

`δ` のときは、最終段の前に**「同型の生成元での値」を先に切った**
(第 37 ブロック `freeYonedaEquiv_freeYonedaInfIso_inv`)。
★★それで最終段(第 40)が通った。

★★★Beck–Chevalley の mate でも同じ手が要る:

    freeYonedaEquiv ((restrictFreeYonedaIso V W).inv) = freeMk (W ⊓ V ⟶ W の包含)

## ★★機構

`restrictFreeYonedaIso` は `eqToIso ≪≫ free.mapIso (yonedaRestrictIso)` で作った。
★`eqToIso` の元になる等式(`restrict_freeObj`)は **`rfl`** なので、
その `inv` は**恒等射**である(実測)。
★★残るのは `free.map` を生成元に当てるだけ(第 37 ブロックの `free_map_app_freeMk`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `eqToIso_restrict_freeObj_inv` | ★`eqToIso` の `inv` は恒等射 |
| `freeYonedaEquiv_restrictFreeYonedaIso_inv` | ★★★★**制限の同型の生成元での値** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {Y : Scheme.{u}} (V W : Y.Opens)

/-- ★**`eqToIso` の `inv` は恒等射である**——元の等式が `rfl` だから。 -/
theorem eqToIso_restrict_freeObj_inv :
    (eqToIso (restrict_freeObj V (yoneda.obj W))).inv
      = 𝟙 ((restrictPresheafFunctor Y V).obj
          ((PresheafOfModules.free
            (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj (yoneda.obj W))) := rfl

/-- ★★★★**制限の同型の生成元での値**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Beck–Chevalley の mate の最終段に要る——
`δ` のときの第 37 ブロックと**同じ発想**である。 -/
theorem freeYonedaEquiv_restrictFreeYonedaIso_inv :
    PresheafOfModules.freeYonedaEquiv ((restrictFreeYonedaIso V W).inv)
      = (ModuleCat.freeMk (homOfLE (inf_le_left : W ⊓ V ≤ W))
          : ((restrictPresheafFunctor Y V).obj
              ((PresheafOfModules.free
                (Y.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u})).obj
                (yoneda.obj W))).obj (op (objOn V W))) := by
  rw [freeYonedaEquiv_apply]
  show ((PresheafOfModules.free _).map (yonedaRestrictIso V W).inv
      ≫ (eqToIso (restrict_freeObj V (yoneda.obj W))).inv).app _ _ = _
  erw [CategoryTheory.comp_apply, free_map_app_freeMk]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_restrictFreeYonedaIso_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限の同型の生成元での値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
