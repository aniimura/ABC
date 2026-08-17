import ABC3.Found.Arakelov.PicRestrict
import Mathlib.Algebra.Category.ModuleCat.Presheaf.PushforwardZeroMonoidal

/-!
# Arakelov (B1) 第 5 ブロック —— **制限はテンソル積と両立する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★実測: mathlib が既に持っていた

「制限がテンソル積と可換」は自作せねばならないと見積もっていたが、**誤りだった**。
★★★mathlib の `Presheaf/PushforwardZeroMonoidal.lean` に

    instance : (pushforward₀OfCommRingCat F R).Monoidal

が在る(2026-08-17 実測)。★これは「前層加群の押し出し関手はモノイダル関手である」
——すなわち **`(P ⊗ Q)|_U ≅ P|_U ⊗ Q|_U`** そのものである。

★★★C1 のとき「複素解析空間と GAGA が要る」と見積もったのが誤りだったのと同じで、
**「無い」と決める前に測る**のが効いた。これで 3 度目である。

## ★本ブロックが渡す橋

`X.PresheafOfModules` から開集合 `U` 上の前層加群の圏への制限を、
mathlib の `pushforward₀OfCommRingCat` に当てて取り、
そのモノイダル性(`μIso`)を使える形にする。

## ★★局所論法の残り

| 部品 | 状態 |
|---|---|
| 開集合への制限(層) | ✅ 第 4 ブロック |
| 制限がテンソル積と両立(前層) | ★★★★**本ブロック** |
| 局所全単射が局所的に確かめられること | ★残 |
| 可逆層の局所自明化 | ★残 |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory

variable (X : Scheme.{u}) (U : X.Opens)

/-! ## ★★前層加群の制限 -/

/-- ★**開集合 `U` の上の前層加群の圏**(係数は `𝒪_X` を `Over U` へ引いたもの)。 -/
noncomputable abbrev PresheafModulesOn : Type _ :=
  PresheafOfModules.{u}
    (((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})

/-- ★★**前層加群を開集合へ制限する関手**。

★mathlib の `pushforward₀OfCommRingCat` を `Over.forget U` に当てる。 -/
noncomputable abbrev restrictPresheafFunctor :
    X.PresheafOfModules ⥤ PresheafModulesOn X U :=
  PresheafOfModules.pushforward₀OfCommRingCat (Over.forget U) X.presheaf

/-- ★制限先の圏もモノイダル圏である。 -/
noncomputable instance : MonoidalCategory (PresheafModulesOn X U) :=
  inferInstanceAs (MonoidalCategory (PresheafOfModules.{u}
    (((Over.forget U).op ⋙ X.presheaf) ⋙ forget₂ CommRingCat.{u} RingCat.{u})))

/-- ★★★**制限はモノイダル関手である**——mathlib の在庫。 -/
noncomputable instance : (restrictPresheafFunctor X U).Monoidal :=
  inferInstanceAs (PresheafOfModules.pushforward₀OfCommRingCat
    (Over.forget U) X.presheaf).Monoidal

variable {X U}

/-- ★★★★**制限はテンソル積と両立する**: `P|_U ⊗ Q|_U ≅ (P ⊗ Q)|_U`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが局所論法の要である——「可逆層は局所的に構造層」を使うとき、
テンソル積を制限の側へ運ぶのに要る。 -/
noncomputable def restrictPresheafTensor (P Q : X.PresheafOfModules) :
    (restrictPresheafFunctor X U).obj P ⊗ (restrictPresheafFunctor X U).obj Q
      ≅ (restrictPresheafFunctor X U).obj (P ⊗ Q) :=
  Functor.Monoidal.μIso (restrictPresheafFunctor X U) P Q

/-- ★**制限は単位を単位へ送る**。 -/
noncomputable def restrictPresheafUnit :
    𝟙_ (PresheafModulesOn X U) ≅ (restrictPresheafFunctor X U).obj (𝟙_ X.PresheafOfModules) :=
  Functor.Monoidal.εIso (restrictPresheafFunctor X U)

/-! ## ★出典の紐付け(`.src`) -/

def restrictPresheafFunctor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——前層加群の開集合への制限)",
    sectionId := "genell-def-1-1-i" }

def restrictPresheafTensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限がテンソル積と両立すること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
