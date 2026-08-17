import ABC3.Found.Arakelov.PicSheafTensor
import Mathlib.Algebra.Category.ModuleCat.Sheaf.PushforwardContinuous

/-!
# Arakelov (B1) 第 4 ブロック —— **層加群の開集合への制限**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★なぜ制限が要るか

残る 2 公理(**結合律**と**逆元**)はどちらも

    「可逆層は局所的に構造層である」

を使う局所論法に帰着する。★★そのためには**開集合への制限**が要る。

★★★2026-08-17 の実測で、mathlib が
`Algebra/Category/ModuleCat/Sheaf/PushforwardContinuous.lean` に
**`SheafOfModules.over M X`(スライス圏 `Over X` への制限)を持っている**
ことが分かった。★探りで、スキームの開集合 `U : X.Opens` に対して
そのまま使えることを確認した。

## ★本ブロックで取れるもの

| 名前 | 内容 |
|---|---|
| `restrictOpen` | 層加群を開集合へ制限する |
| `restrictOpenFunctor` | その関手性 |
| `restrictOpenIso` | 同型は同型に移る(★局所論法の要) |
| `restrictOpen_unit` | 構造層の制限は構造層 |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

variable (X : Scheme.{u})

/-! ## ★★開集合への制限 -/

/-- ★**開集合 `U` 上の層加群の圏**(係数は `𝒪_X|_U`)。 -/
noncomputable abbrev ModulesOn (U : X.Opens) : Type _ :=
  SheafOfModules (X.ringCatSheaf.over U)

/-- ★★**層加群を開集合へ制限する関手**。

★★★mathlib の `SheafOfModules.overFunctor` をスキームの開集合に当てる。 -/
noncomputable abbrev restrictOpenFunctor (U : X.Opens) : X.Modules ⥤ ModulesOn X U :=
  SheafOfModules.overFunctor X.ringCatSheaf U

variable {X}

/-- ★**層加群の開集合への制限**。 -/
noncomputable abbrev restrictOpen (M : X.Modules) (U : X.Opens) : ModulesOn X U :=
  M.over U

/-- ★★**同型は制限しても同型である**。

★★★これが局所論法の要である——「可逆層は局所的に構造層」を使うとき、
局所自明化の同型を制限側へ運ぶのに要る。 -/
noncomputable def restrictOpenIso {M N : X.Modules} (e : M ≅ N) (U : X.Opens) :
    restrictOpen M U ≅ restrictOpen N U :=
  (restrictOpenFunctor X U).mapIso e

/-- ★**構造層の制限は構造層である**。 -/
noncomputable def restrictOpen_unit (U : X.Opens) :
    restrictOpen (unitModules X) U ≅ SheafOfModules.unit (X.ringCatSheaf.over U) :=
  Iso.refl _

/-- ★制限は射に沿って関手的である。 -/
noncomputable abbrev restrictOpenMap {M N : X.Modules} (f : M ⟶ N) (U : X.Opens) :
    restrictOpen M U ⟶ restrictOpen N U :=
  (restrictOpenFunctor X U).map f

@[simp] theorem restrictOpenMap_id (M : X.Modules) (U : X.Opens) :
    restrictOpenMap (𝟙 M) U = 𝟙 _ :=
  (restrictOpenFunctor X U).map_id M

@[simp] theorem restrictOpenMap_comp {M N P : X.Modules} (f : M ⟶ N) (g : N ⟶ P)
    (U : X.Opens) :
    restrictOpenMap (f ≫ g) U = restrictOpenMap f U ≫ restrictOpenMap g U :=
  (restrictOpenFunctor X U).map_comp f g

/-! ## ★出典の紐付け(`.src`) -/

def restrictOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層加群の開集合への制限)",
    sectionId := "genell-def-1-1-i" }

def restrictOpenIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限が同型を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
