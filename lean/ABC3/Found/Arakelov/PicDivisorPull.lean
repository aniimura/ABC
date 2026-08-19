import ABC3.Found.Arakelov.PicGoodOpen

/-!
# Arakelov (B2) 第 240 ブロック —— ★★★★★★★★**`𝒪(D)` の引き戻し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★B2 の 14 欄目——`CartierPicData` が埋まった

    f^* 𝒪_Y(D) = 𝒪_X(f^{-1}D)      (f は平坦、D は Cartier)

★★★これが `CartierPicData.ofDivisor_pullback` であり、
**B2 の最後の欄**である(2026-08-19)。

## ★同型の連鎖

| 段 | 同型 |
|---|---|
| 1 | `f^*(𝒪(D)) ≅ 層化((f^*_pre)(idealPresheaf D))` (第 41 `pullbackValIso`) |
| 2 | `層化((f^*_pre)(idealPresheaf D)) ≅ 層化(idealPresheaf (D.comap f))` (★第 239) |
| 3 | `層化(idealPresheaf (D.comap f)) ≅ 𝒪(f^{-1}D)` (第 12 `sheafifyValIso`) |

★★★第 213–239 の 27 ブロックが、この 3 行のために積み上げられていた。

## ★平坦性が要る理由(Interface に記録済み)

`Y = Spec k[x]`、`D = (x)`、`f : Spec k ⟶ Spec k[x]` を原点とすると
`D.comap f = ⊥` で可逆でない。★イデアル層の引き戻しは可逆とは限らない
——可逆**層** `𝒪(D)` の引き戻しは常に可逆だが、両者が一致するのは平坦なときである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) [AlgebraicGeometry.Flat f] (D : Y.IdealSheafData)

/-- ★★★★★★★★**`𝒪(D)` は平坦な引き戻しと両立する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `CartierPicData.ofDivisor_pullback`——**B2 の 14 欄目**である。 -/
theorem ofDivisorSheaf_pullback (hD : IsCartier Y D) :
    picPullback f (ofDivisorSheaf D) = ofDivisorSheaf (D.comap f) := by
  classical
  show PicSheaf.mk (InvSheaf.pullback f (divisorInvSheaf D))
    = PicSheaf.mk (divisorInvSheaf (D.comap f))
  refine PicSheaf.mk_eq_mk ?_
  show (Scheme.Modules.pullback f).obj (divisorInvSheaf D).carrier
    ≅ (divisorInvSheaf (D.comap f)).carrier
  rw [divisorInvSheaf_carrier D hD,
    divisorInvSheaf_carrier (D.comap f) (isCartier_comap f D hD)]
  exact pullbackValIso f (idealSheaf D)
    ≪≫ @asIso _ _ _ _ ((sheafifyFunctor X).map (pullIdealHom f D))
        (isIso_sheafify_pullIdealHom f D hD)
    ≪≫ sheafifyValIso (idealSheaf (D.comap f))

/-! ## ★出典の紐付け(`.src`) -/

def ofDivisorSheaf_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——𝒪(D) は平坦な引き戻しと両立する)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
