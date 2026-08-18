import ABC3.Found.Arakelov.PicTildeStock

/-!
# Arakelov (B1) 第 64 ブロック —— **`tilde` は左随伴である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これで生成元の議論がまた使える

mathlib は

    tilde.adjunction : tilde.functor R ⊣ moduleSpecΓFunctor

を持つ(2026-08-18 実測)。★★したがって **`tilde` は余極限を保つ**。

★★★`ModuleCat R` は `R` 自身で生成されるので、
「`tilde` がテンソル積を保つ」は**生成元の上で見れば済む**——
`tilde R = 𝒪` だからである。

★★★★これは `δ`(第 40)と Beck–Chevalley(第 54)で**2 度使った手口**の 3 度目である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tildeAdjunction` | ★`tilde ⊣ Γ`(mathlib の在庫、機械確認) |
| `tildeFunctor_preservesColimits` | ★★★**`tilde` は余極限を保つ** |

## ★★★残り(B1)

| # | 主張 |
|---|---|
| 1 | `tilde (M ⊗ N) ≅ tensorModules (tilde M) (tilde N)`(★生成元 + 本ブロック) |
| 2 | 可逆層は `tilde` の本質像に入る |
| 3 | `PicSheaf (Spec R) ≃* CommRing.Pic R` |
| 4 | `PicardData` の組み立て |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★**`tilde ⊣ Γ`**(mathlib の在庫、機械確認)。 -/
noncomputable abbrev tildeAdjunction :
    tildeFunctor R ⊣ AlgebraicGeometry.moduleSpecΓFunctor :=
  AlgebraicGeometry.tilde.adjunction

/-- ★★★**`tilde` は余極限を保つ**——左随伴だからである。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで「`tilde` がテンソル積を保つ」を**生成元の上で見れば済む**
——`δ`(第 40)・Beck–Chevalley(第 54)と**同じ手口の 3 度目**である。 -/
@[reducible] noncomputable def tildeFunctor_preservesColimits :
    PreservesColimits (tildeFunctor R) :=
  (tildeAdjunction R).leftAdjoint_preservesColimits

/-! ## ★出典の紐付け(`.src`) -/

def tildeFunctor_preservesColimits.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde が左随伴で余極限を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
