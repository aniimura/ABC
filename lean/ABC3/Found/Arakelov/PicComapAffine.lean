import ABC3.Found.Arakelov.PicPushoutKer
import ABC3.Found.Arakelov.PicPrincipalAffine
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial
import Mathlib.AlgebraicGeometry.Morphisms.Affine

/-!
# Arakelov (B2) 第 195 ブロック —— **アフィンでの `comap` の記述**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`comap` が**イデアルの拡大**であることが出た(アフィンの場合)

    X, Y アフィン ⟹ (D.comap f).ideal ⊤ = (D.ideal ⊤).map f.appTop

## ★★4 つの在庫を繋いだだけ

| 段 | 使ったもの | 出どころ |
|---|---|---|
| 1 | `comap = (pullback.fst f ι).ker` | mathlib(定義) |
| 2 | `Hom.ker_apply`(準コンパクトなら核は成分ごと) | ★mathlib |
| 3 | `isPushout_appTop_of_isPullback`(アフィンの引き戻しは環の押し出し) | ★mathlib |
| 4 | 押し出しの核はイデアルの拡大 | ★第 194 |
| 5 | `ker_subschemeι_app`、`subschemeι_app_surjective` | ★mathlib |

★★★準コンパクト性は**自動**である(`subschemeι` は閉埋め込みで、
準コンパクト性は底変換で保たれる)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_ideal_top` | ★★★★**アフィンでの `comap` の記述** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

/-- ★★★★**アフィンでは `comap` はイデアルの拡大である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `isCartierDivisor_comap` / `ofDivisor_pullback` の核である。 -/
theorem comap_ideal_top [IsAffine X] [IsAffine Y] :
    (D.comap f).ideal ⟨⊤, isAffineOpen_top X⟩
      = (D.ideal ⟨⊤, isAffineOpen_top Y⟩).map f.appTop.hom := by
  haveI : IsAffine D.subscheme := isAffine_of_isAffineHom D.subschemeι
  have hpb : IsPullback (Limits.pullback.fst f D.subschemeι)
      (Limits.pullback.snd f D.subschemeι) f D.subschemeι := IsPullback.of_hasPullback _ _
  have hpo := isPushout_appTop_of_isPullback hpb
  have hker := ker_of_isPushout hpo (D.subschemeι_app_surjective ⟨⊤, isAffineOpen_top Y⟩)
  have h1 : RingHom.ker ((D.subschemeι).appTop).hom = D.ideal ⟨⊤, isAffineOpen_top Y⟩ :=
    Scheme.IdealSheafData.ker_subschemeι_app _ ⟨⊤, isAffineOpen_top Y⟩
  rw [Scheme.IdealSheafData.comap, Scheme.Hom.ker_apply _ ⟨⊤, isAffineOpen_top X⟩]
  show RingHom.ker ((Limits.pullback.fst f D.subschemeι).appTop).hom = _
  rw [hker, h1]

/-! ## ★出典の紐付け(`.src`) -/

def comap_ideal_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィンでの comap の記述)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
