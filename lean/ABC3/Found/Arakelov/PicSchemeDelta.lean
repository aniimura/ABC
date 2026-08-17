import ABC3.Found.Arakelov.PicSchemePullback
import ABC3.Found.Arakelov.PicPullOplax

/-!
# Arakelov (B1) 第 22 ブロック —— **スキームの射に沿った比較射 `δ`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 20 ブロックが層の段に直結した

第 20 ブロックは**抽象の**前層引き戻しについて oplax monoidal 性を作った。
★★スキームの射 `f : X ⟶ Y` について、mathlib の環層射は

    Scheme.Hom.toRingCatSheafHom f |>.hom = Functor.whiskerRight f.c (forget₂ CommRingCat RingCat)

であり、これは第 20 ブロックの `alphaR` **そのもの**である(`rfl` で一致、2026-08-18 実測)。
★★★したがって **`δ` がスキームの射について得られる**。

## ★★★★★層と前層を繋ぐ橋

mathlib の `SheafOfModules.pullbackIso` は

    層の引き戻し = 前層に落とす ⋙ 前層の引き戻し ⋙ 層化

を与える。さらに `sheafificationCompPullback` は

    層化 ⋙ 層の引き戻し ≅ 前層の引き戻し ⋙ 層化

を与える。★★**これで前層の段の `δ` を層の段へ運べる。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pullbackPre` | ★スキームの射に沿った**前層**の引き戻し |
| `pullbackPreOplax` | ★★★それが oplax monoidal |
| `pullbackDelta` | ★★★★**比較射** `f^*(P ⊗ Q) ⟶ f^*P ⊗ f^*Q` |
| `pullbackEta` | ★★`f^*𝒪 ⟶ 𝒪`(前層の段) |
| `sheafifyPullbackApp` | ★★★★★**層化と引き戻しの交換**(対象ごと) |

★★★残るのは `pullbackDelta` が(可逆層について)**同型であること**だけである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-! ## ★スキームの射に沿った前層の引き戻し -/

/-- ★**スキームの射に沿った前層加群の引き戻し**。

★★第 20 ブロックの `alphaR` に `f.c`(構造層の射)を入れたものである。 -/
noncomputable abbrev pullbackPre : Y.PresheafOfModules ⥤ X.PresheafOfModules :=
  PresheafOfModules.pullback (alphaR (Opens.map f.base) X.presheaf Y.presheaf f.c)

/-- ★★★**mathlib の環層射は第 20 ブロックの `alphaR` そのものである**。

★★これが「抽象で作ったものがスキームに効く」ことの根拠である。 -/
theorem toRingCatSheafHom_hom_eq :
    f.toRingCatSheafHom.hom = alphaR (Opens.map f.base) X.presheaf Y.presheaf f.c := rfl

/-- ★★★**前層の引き戻しは oplax monoidal である**(第 20 ブロックの適用)。 -/
noncomputable instance pullbackPreOplax : (pullbackPre f).OplaxMonoidal :=
  pullbackOplax (Opens.map f.base) X.presheaf Y.presheaf f.c

/-! ## ★★★★比較射 -/

/-- ★★★★★**比較射** `f^*(P ⊗ Q) ⟶ f^*P ⊗ f^*Q`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Pic` の**引き戻しが群準同型になる**ための射である。
★残るのは、可逆層についてこれが**同型**であることだけである。 -/
noncomputable def pullbackDelta (P Q : Y.PresheafOfModules) :
    (pullbackPre f).obj (P ⊗ Q) ⟶ (pullbackPre f).obj P ⊗ (pullbackPre f).obj Q :=
  Functor.OplaxMonoidal.δ (pullbackPre f) P Q

/-- ★★**構造層の比較射**(前層の段) `f^*𝒪_Y ⟶ 𝒪_X`。

★層の段では第 21 ブロックの `pullbackUnitIso` が**同型**を与える。 -/
noncomputable def pullbackEta : (pullbackPre f).obj (𝟙_ Y.PresheafOfModules) ⟶ 𝟙_ X.PresheafOfModules :=
  Functor.OplaxMonoidal.η (pullbackPre f)

/-! ## ★★★★★層化と引き戻しの交換 -/

/-- ★層化と引き戻しの交換(mathlib の自然な型で)。

★★★**2 段構えにする理由**は第 21 ブロックの `pullbackUnitIso` と同じ——
`Scheme.Modules.pullback` は `def` なので、期待型を先に与えると
インスタンス探索がメタ変数に阻まれる。 -/
noncomputable def sheafifyPullbackIsoAux :
    PresheafOfModules.sheafification (𝟙 Y.ringCatSheaf.obj)
        ⋙ SheafOfModules.pullback f.toRingCatSheafHom
      ≅ PresheafOfModules.pullback f.toRingCatSheafHom.hom
        ⋙ PresheafOfModules.sheafification (R₀ := X.ringCatSheaf.obj) (𝟙 X.ringCatSheaf.obj) :=
  SheafOfModules.sheafificationCompPullback f.toRingCatSheafHom

/-- ★★★★★**層化と引き戻しは交換する**(対象ごと)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが**前層の段の `δ` を層の段へ運ぶ橋**である。 -/
noncomputable def sheafifyPullbackApp (P : Y.PresheafOfModules) :
    (Scheme.Modules.pullback f).obj ((sheafifyFunctor Y).obj P)
      ≅ (sheafifyFunctor X).obj ((pullbackPre f).obj P) :=
  (sheafifyPullbackIsoAux f).app P

/-! ## ★出典の紐付け(`.src`) -/

def pullbackDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しの比較射)",
    sectionId := "genell-def-1-1-i" }

def sheafifyPullbackApp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層化と引き戻しの交換)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
