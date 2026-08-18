import ABC3.Found.Arakelov.PicLTPull
import ABC3.Found.Arakelov.PicSheafPullTensor

/-!
# Arakelov (B1) 第 60 ブロック —— ★★★★★★★★★★**第 41 の仮定が落ちた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★仮定なしの形になった

第 41 ブロック `pullbackSheafTensorIso` は 2 つの `IsLocallyRankOne` 仮定を持っていた。
★§9-49 で「**放置できない**」と旗を立て、§9-52 で
「**Beck–Chevalley を作るしかない**」と測り、第 43–59 で作った。

★★★本ブロックで、その仮定を**局所自明性だけ**から出す:

| 仮定 | 出どころ |
|---|---|
| `IsLocallyRankOne X (f^*_pre M.val)` | ★第 59(輸送)+ 第 15(自明 ⟹ 階数 1) |
| `IsLocallyRankOne X (層化 (f^*_pre L.val)).val` | ★★第 59 + 第 16(層化は自明性を保つ)+ 第 15 |

## ★★取れたもの

    f^*(L ⊗ M) ≅ f^*L ⊗ f^*M     (層の段、**局所自明性だけを仮定**)

★★★これで `PicardData.pullback_mul` が書ける——
可逆層は定義から局所自明だからである。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★**引き戻した前層は局所階数 1 である**(局所自明性から)。 -/
theorem isLocallyRankOne_pullbackPre (L : Y.PresheafOfModules) (hL : IsLocallyTrivial Y L) :
    IsLocallyRankOne X ((pullbackPre f).obj L) :=
  (isLocallyTrivial_pullbackPre f L hL).isLocallyRankOne

/-- ★★★**引き戻して層化したものも局所階数 1 である**。 -/
theorem isLocallyRankOne_sheafify_pullbackPre (L : Y.PresheafOfModules)
    (hL : IsLocallyTrivial Y L) :
    IsLocallyRankOne X ((sheafifyFunctor X).obj ((pullbackPre f).obj L)).val :=
  (isLocallyTrivial_sheafify X ((pullbackPre f).obj L)
    (isLocallyTrivial_pullbackPre f L hL)).isLocallyRankOne

/-- ★★★★★★★★★★**層の段で引き戻しはテンソル積を保つ**——**局所自明性だけを仮定**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★§9-49 で旗を立てた仮定が、第 43–59 の 17 ブロックで**落ちた**。
★これで `PicardData.pullback_mul` が書ける。 -/
noncomputable def pullbackSheafTensorIso' (L M : Y.Modules)
    (hL : IsLocallyTrivial Y L.val) (hM : IsLocallyTrivial Y M.val) :
    (Scheme.Modules.pullback f).obj (tensorModules L M)
      ≅ tensorModules ((Scheme.Modules.pullback f).obj L)
          ((Scheme.Modules.pullback f).obj M) :=
  pullbackSheafTensorIso f L M
    (isLocallyRankOne_pullbackPre f M.val hM)
    (isLocallyRankOne_sheafify_pullbackPre f L.val hL)

/-! ## ★出典の紐付け(`.src`) -/

def pullbackSheafTensorIso'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層の段で引き戻しがテンソル積を保つこと(仮定なし))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
