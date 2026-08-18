import ABC3.Found.Arakelov.PicUnitBij

/-!
# Arakelov (B1) 第 108 ブロック —— **`unitHomOfSection` の全単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★合成して閉じる

    unitHomOfSection s = (freeYonedaTermIso).inv ≫ freeYonedaEquiv.symm s   (`rfl`)

★**前者は同型**(`PresheafOfModules.evaluation` で `app` へ移す)、
★★**後者は第 107 ブロック**で全単射である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `unitHomOfSection_app_eq` | ★`app` は合成(`rfl`) |
| `bijective_unitHomOfSection_app` | ★★★★★**`unitHomOfSection` の `app` の全単射性** |

## ★★★次

これで第 102 ブロック(被覆で全単射なら同型)が当たり、
`(tilde M)|_{D g} ≅ 𝟙_` が出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens) (P : PresheafModulesOn X V)

/-- ★**`unitHomOfSection` の `app` は合成である**(`rfl`)。 -/
theorem unitHomOfSection_app_eq (s : P.obj (op (Over.mk (𝟙 V)))) (W : Over V) :
    (unitHomOfSection V P s).app (op W)
      = (freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
          (overTerminalUnique V)).inv.app (op W)
        ≫ (PresheafOfModules.freeYonedaEquiv.symm s).app (op W) := rfl

/-- ★★★★★**`unitHomOfSection` の `app` は、生成元の乗法が全単射なら全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 102 ブロック(被覆で全単射なら同型)が当たる。 -/
theorem bijective_unitHomOfSection_app (s : P.obj (op (Over.mk (𝟙 V)))) (W : Over V)
    [Unique ((yoneda.obj (Over.mk (𝟙 V))).obj (op W))]
    (hb : Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
        ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
      c • restrictSec V P s W)) :
    Function.Bijective ((unitHomOfSection V P s).app (op W)) := by
  haveI : IsIso ((freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
      (overTerminalUnique V)).inv.app (op W)) :=
    inferInstanceAs (IsIso ((PresheafOfModules.evaluation _ (op W)).map
      (freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
        (overTerminalUnique V)).inv))
  have h1 := ConcreteCategory.bijective_of_isIso
    ((freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
      (overTerminalUnique V)).inv.app (op W))
  have h2 := bijective_freeYonedaEquiv_symm_app V P s W hb
  rw [unitHomOfSection_app_eq]
  exact Function.Bijective.comp (f := ⇑(ConcreteCategory.hom
      ((freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
        (overTerminalUnique V)).inv.app (op W))))
    (g := ⇑(ConcreteCategory.hom
      ((PresheafOfModules.freeYonedaEquiv.symm s).app (op W)))) h2 h1

/-! ## ★出典の紐付け(`.src`) -/

def bijective_unitHomOfSection_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unitHomOfSection の全単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
