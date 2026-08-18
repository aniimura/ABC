import ABC3.Found.Arakelov.PicUnitSecBij

/-!
# Arakelov (B1) 第 109 ブロック —— **制限した前層加群でも「被覆で全単射なら同型」**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`IsLocallyTrivial` が要求する形そのもの

`IsLocallyTrivial` は `(restrict V).obj P ≅ 𝟙_`——
**`Over V` 上の前層加群の同型**を要求する。

★★第 102 ブロックは `Sheaf J Ab` の版だった。本ブロックはそれを
**`PresheafModulesOn X V`** に持ち込む。

## ★★機構 —— 4 段の梱包と反射

| 段 | 内容 |
|---|---|
| 1 | `⟨P.presheaf, hP⟩` を `Sheaf ((Opens.grothendieckTopology X).over V) Ab` に梱包 |
| 2 | 射は `⟨(toPresheaf _).map φ⟩` |
| 3 | 第 102 ブロックで `IsIso`(層の射として) |
| 4 | `sheafToPresheaf` と `toPresheaf` で 2 段反射 |

★**一発で通った**——§9-73 で詰まった `SheafOfModules.toSheaf` を経由しないのが要点である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_of_bijective_on_cover_mod` | ★★★★★★**被覆で全単射なら前層加群の射は同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {V : X.Opens} {P Q : PresheafModulesOn X V}

/-- ★★★★★★**被覆の上で全単射なら、制限した前層加群の射は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `IsLocallyTrivial` の同型を作るための最後の器具である。 -/
theorem isIso_of_bijective_on_cover_mod
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) P.presheaf)
    (hQ : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) Q.presheaf)
    (φ : P ⟶ Q)
    (h : ∀ W : Over V, ∃ S : Sieve W, S ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), S.arrows i →
        Function.Bijective ((PresheafOfModules.toPresheaf _).map φ |>.app (op Z))) :
    IsIso φ := by
  haveI : IsIso ((⟨(PresheafOfModules.toPresheaf _).map φ⟩ :
      (⟨P.presheaf, hP⟩ : Sheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u})
        ⟶ ⟨Q.presheaf, hQ⟩)) :=
    isIso_of_bijective_on_cover _ _ h
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map φ) := by
    have := (sheafToPresheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u}).map_isIso
      (⟨(PresheafOfModules.toPresheaf _).map φ⟩ :
      (⟨P.presheaf, hP⟩ : Sheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u})
        ⟶ ⟨Q.presheaf, hQ⟩)
    exact this
  exact isIso_of_reflects_iso φ (PresheafOfModules.toPresheaf _)

/-! ## ★出典の紐付け(`.src`) -/

def isIso_of_bijective_on_cover_mod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した前層加群でも被覆で全単射なら同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
