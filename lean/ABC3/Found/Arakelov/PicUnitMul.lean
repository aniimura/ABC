import ABC3.Found.Arakelov.PicIdealLT
import ABC3.Found.Arakelov.PicUnitBij

/-!
# Arakelov (B2) 第 163 ブロック —— **単位対象の「`c` 倍」自己射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層の双対へ向けた第 1 歩

`ofDivisor` には**逆層**が要り、それには双対

    F^∨(U) := ((restrictPresheafFunctor X U).obj F ⟶ 𝟙_)

を作る必要がある(mathlib に `SheafOfModules` の内部 Hom は無い、2026-08-18 実測)。

★その `Γ(X,U)` 加群構造は `c • φ := φ ≫ unitMul c` で入れる。

## ★★逃げ道——**既にあった**

「`c` 倍」自己射は、第 108 の `unitHomOfSection` を `P := 𝟙_` に当てたものである:

    unitMul c = unitHomOfSection U (𝟙_) c

★★★手で `app` と `naturality` を書こうとして 4 回失敗した後に気づいた
——**単位対象の切断は `Γ(X,U)` そのもの**だから、`unitHomOfSection` がそのまま使える。

★★[[ring-instance-two-paths]] の 9 例目でもある:
`𝟙_` の切断への `Γ(X,U)` 作用は、`CommRingCat` 綴りでは見つからず
`RingCat` 綴りでのみ見つかる——`unitHomOfSection` はその差を吸収している。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

/-- ★★**単位対象の「`c` 倍」自己射**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★第 108 の `unitHomOfSection` を `P := 𝟙_` に当てたものである。 -/
noncomputable def unitMul (c : ((𝟙_ (PresheafModulesOn X U)).obj (op (Over.mk (𝟙 U))) : Type u)) :
    𝟙_ (PresheafModulesOn X U) ⟶ 𝟙_ (PresheafModulesOn X U) :=
  unitHomOfSection U (𝟙_ (PresheafModulesOn X U)) c

/-! ## ★出典の紐付け(`.src`) -/

def unitMul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——単位対象の c 倍自己射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
