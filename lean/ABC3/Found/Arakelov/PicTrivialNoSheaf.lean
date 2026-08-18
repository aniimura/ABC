import ABC3.Found.Arakelov.PicTrivialIso

/-!
# Arakelov (B2) 第 170 ブロック —— **層の仮定なしの同型判定**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★双対には層の仮定が使えない

第 110・115 の同型判定は `P` が**層である**ことを要求する。
★ところが双対 `F^∨`(第 169)は層であることをまだ示していない。

★★しかし双対の場合、**すべての `W` で全単射**が言える
(自明化 `e : F|_V ≅ 𝟙_` が `Hom(F|_W, 𝟙_) ≅ End(𝟙_) ≅ Γ(X,W)` を与えるため)。
★★★そのときは**層の仮定は要らない**——成分ごとに同型なら射も同型である。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `isIso_unitHomOfSection'` | ★★層の仮定なしの `IsIso` |
| `trivialIsoOfSection'` | ★★★層の仮定なしの `𝟙_ ≅ P` |

★★★★これで双対の局所自明性が**層であることを示す前に**言える。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★★すべての `W` で全単射なら層の仮定なしで同型。 -/
theorem isIso_unitHomOfSection' (V : X.Opens) (P : PresheafModulesOn X V)
    (s : P.obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, Function.Bijective
      (fun c : ((((Over.forget V).op ⋙ X.presheaf)
          ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
        c • restrictSec V P s W)) :
    IsIso (unitHomOfSection V P s) := by
  haveI happ : ∀ W, IsIso ((unitHomOfSection V P s).app W) := by
    intro W
    exact (ConcreteCategory.isIso_iff_bijective _).2
      (bijective_unitHomOfSection_app V P s W.unop (h W.unop))
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (unitHomOfSection V P s)) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro W
    haveI := happ W
    exact inferInstanceAs
      (IsIso ((forget₂ _ AddCommGrpCat).map ((unitHomOfSection V P s).app W)))
  exact isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)

/-- ★★★層の仮定なしの `𝟙_ ≅ P`。 -/
noncomputable def trivialIsoOfSection' (V : X.Opens) (P : PresheafModulesOn X V)
    (s : P.obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, Function.Bijective
      (fun c : ((((Over.forget V).op ⋙ X.presheaf)
          ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op W) : Type u) =>
        c • restrictSec V P s W)) :
    𝟙_ (PresheafModulesOn X V) ≅ P :=
  @asIso _ _ _ _ (unitHomOfSection V P s) (isIso_unitHomOfSection' V P s h)


/-! ## ★出典の紐付け(`.src`) -/

def trivialIsoOfSection'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層の仮定なしの同型判定)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
