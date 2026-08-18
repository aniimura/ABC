import ABC3.Found.Arakelov.PicSieveTransport

/-!
# Arakelov (B1) 第 59 ブロック —— ★★★★★★★★★★**局所自明性は引き戻しで保たれる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★取れたもの

    IsLocallyTrivial Y L  ⟹  IsLocallyTrivial X (f^*_pre L)

★★★これが §9-49 で旗を立て、§9-52 で「**mathlib に無いので我々が作るしかない**」
と測った Beck–Chevalley の**目的**であった。

## ★★★★4 つの部品を繋いだだけ(一発で通った)

| 部品 | ブロック |
|---|---|
| 被覆の輸送 | ★第 58 |
| Beck–Chevalley | ★★★★第 54 |
| `(f|)^*_pre 𝒪 ≅ 𝒪` | ★第 56 |
| 制限の推移律(`rfl`) | ★第 57 |

★★証明は `calc` 5 段——**すべて既存の同型を繋ぐだけ**である。

## ★★★これで第 41 ブロックの仮定が落とせる

第 41 ブロック `pullbackSheafTensorIso` の仮定は
`IsLocallyRankOne X ((pullbackPre f).obj M.val)` であった。
★★局所自明なら局所階数 1 である(第 15 ブロック `IsLocallyTrivial.isLocallyRankOne`)ので、
**可逆層について仮定が満たされる**。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★★★★★★★**局所自明性は引き戻しで保たれる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★機構は 4 つの部品を繋ぐだけである:
第 58(被覆の輸送)・第 54(Beck–Chevalley)・第 56(`(f|)^* 𝒪 ≅ 𝒪`)・
第 57(制限の推移律)。 -/
theorem isLocallyTrivial_pullbackPre (L : Y.PresheafOfModules) (hL : IsLocallyTrivial Y L) :
    IsLocallyTrivial X ((pullbackPre f).obj L) := by
  intro U
  obtain ⟨S, hS, hiso⟩ := hL ⊤
  refine ⟨transportSieve f S U, transportSieve_mem f S hS U, ?_⟩
  rintro V i ⟨W, hW, hVW⟩
  obtain ⟨e⟩ := hiso (homOfLE le_top : W ⟶ ⊤) hW
  refine ⟨?_⟩
  calc (restrictPresheafFunctor X V).obj ((pullbackPre f).obj L)
      ≅ (restrictOnFunctor hVW).obj
          ((restrictPresheafFunctor X ((Opens.map f.base).obj W)).obj
            ((pullbackPre f).obj L)) := Iso.refl _
    _ ≅ (restrictOnFunctor hVW).obj
          ((pullbackPreOn f W).obj ((restrictPresheafFunctor Y W).obj L)) :=
        (restrictOnFunctor hVW).mapIso (bcIso f W L).symm
    _ ≅ (restrictOnFunctor hVW).obj ((pullbackPreOn f W).obj (𝟙_ (PresheafModulesOn Y W))) :=
        (restrictOnFunctor hVW).mapIso ((pullbackPreOn f W).mapIso e)
    _ ≅ (restrictOnFunctor hVW).obj (𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj W))) :=
        (restrictOnFunctor hVW).mapIso (pullbackOnUnitIso f W)
    _ ≅ 𝟙_ (PresheafModulesOn X V) := Iso.refl _


/-! ## ★出典の紐付け(`.src`) -/

def isLocallyTrivial_pullbackPre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所自明性が引き戻しで保たれること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
