import ABC3.Found.Arakelov.PicLocalTrivial
import Mathlib.CategoryTheory.Sites.PreservesLocallyBijective

/-!
# Arakelov (B1) 第 16 ブロック(再) —— **制限は局所全単射を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★壁を破る 3 部品

残る壁は「**層化が局所自明性を保つ**」であり、その証明は

    η_P|_V : P|_V ⟶ (層化 P).val|_V     -- 局所全単射
    P|_V は層(仮定)、(層化 P).val|_V も層
    ⟹ η_P|_V は同型

という 3 段である。★★★2026-08-17 の実測で、**3 部品すべてが mathlib に在った**:

| 部品 | 在庫 |
|---|---|
| `Over.forget` が cocontinuous | ★`Sites/Over.lean` の instance |
| 制限が局所全射を保つ | ★`Sites/PreservesLocallyBijective.lean` |
| 制限が局所単射を保つ | 同上 |

★★★★**前 turn に「mathlib に該当補題なし」と判定したのは誤りだった**(7 度目)。

## ★★本ブロックが渡す形

スキームの開集合に当てた形で名前を付けておく。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable (X : Scheme.{u}) (V : X.Opens)

/-- ★**`Over.forget V` は cocontinuous である**(mathlib の instance の確認)。 -/
theorem isCocontinuous_overForget :
    (Over.forget V).IsCocontinuous
      ((Opens.grothendieckTopology X).over V) (Opens.grothendieckTopology X) :=
  inferInstance

/-- ★★★**制限は局所全射性を保つ**。 -/
theorem isLocallySurjective_restrict {F G : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u}} (f : F ⟶ G)
    [Presheaf.IsLocallySurjective (Opens.grothendieckTopology X) f] :
    Presheaf.IsLocallySurjective ((Opens.grothendieckTopology X).over V)
      (Functor.whiskerLeft (Over.forget V).op f) :=
  Presheaf.isLocallySurjective_whisker
    ((Opens.grothendieckTopology X).over V) (Opens.grothendieckTopology X) (Over.forget V) f

/-- ★★★**制限は局所単射性を保つ**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これと前定理で「層化の単位を制限しても局所全単射」が言え、
**層どうしの局所全単射は同型**なので壁が破れる。 -/
theorem isLocallyInjective_restrict {F G : (X.Opens)ᵒᵖ ⥤ AddCommGrpCat.{u}} (f : F ⟶ G)
    [Presheaf.IsLocallyInjective (Opens.grothendieckTopology X) f] :
    Presheaf.IsLocallyInjective ((Opens.grothendieckTopology X).over V)
      (Functor.whiskerLeft (Over.forget V).op f) :=
  Presheaf.isLocallyInjective_whisker
    ((Opens.grothendieckTopology X).over V) (Opens.grothendieckTopology X) (Over.forget V) f

/-! ## ★出典の紐付け(`.src`) -/

def isLocallySurjective_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限が局所全射性を保つこと)",
    sectionId := "genell-def-1-1-i" }

def isLocallyInjective_restrict.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限が局所単射性を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
