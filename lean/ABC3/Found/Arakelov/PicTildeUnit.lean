import ABC3.Found.Arakelov.PicTildeAdj

/-!
# Arakelov (B1) 第 65 ブロック —— **`tilde R = 𝒪` は定義から**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★生成元の上での値が `rfl` で出た

    tilde (R として見た R) = 構造層

★★これは **`Iso.refl` で通る**(2026-08-18 実測)——`tilde` の定義がそうなっている。

★★★これが「`tilde` がテンソル積を保つ」を生成元の上で見るときの**出発点**である:

    tilde (R ⊗ R) = tilde R = 𝒪 = 𝒪 ⊗ 𝒪 = tilde R ⊗ tilde R

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tildeUnit` | ★★★★**`tilde R ≅ 𝒪`**(`Iso.refl`) |

## ★★★残り(B1)

| # | 主張 | 材料 |
|---|---|---|
| 1 | `tilde` がテンソル積を保つ | ★第 64(余極限保存)+ 本ブロック(生成元での値) |
| 2 | 可逆層は `tilde` の本質像に入る | ★`isIso_fromTildeΓ_iff` |
| 3 | `PicSheaf (Spec R) ≃* CommRing.Pic R` | ★★1・2 から |
| 4 | `PicardData` の組み立て | ★構成から出る |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★★★★**`tilde R` は構造層である**——`Iso.refl` で通る。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「`tilde` がテンソル積を保つ」を生成元の上で見るときの出発点である。 -/
noncomputable def tildeUnit :
    (tildeFunctor R).obj (ModuleCat.of R R) ≅ unitModules (Spec R) :=
  Iso.refl _

/-! ## ★出典の紐付け(`.src`) -/

def tildeUnit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde R が構造層であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
