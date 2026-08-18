import ABC3.Found.Arakelov.PicPullGroup
import Mathlib.AlgebraicGeometry.Modules.Tilde

/-!
# Arakelov (B1) 第 63 ブロック —— **アフィン比較の在庫**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★B1 の残りは `equivPicRing` だけ

    PicSheaf (Spec R) ≃* CommRing.Pic R

★第 62 ブロックで `pullback` 系 4 欄が埋まったので、残るのはこれである。

## ★★★mathlib の在庫(2026-08-18 実測、本ブロックで機械確認)

| 在庫 | 内容 |
|---|---|
| `AlgebraicGeometry.tilde.functor R` | ★`ModuleCat R ⥤ (Spec R).Modules` |
| `tilde.toTildeΓNatIso` | ★★**`𝟭 ≅ tilde ⋙ Γ`**——`tilde` は充満忠実 |
| `Scheme.Modules.fromTildeΓ` | ★逆向きの比較射 |
| `isIso_fromTildeΓ_iff` | ★★`IsIso (fromTildeΓ M) ↔ M ∈ tilde の本質像` |
| `Module.Invertible R M` / `CommRing.Pic R` | ★`RingTheory/PicardGroup.lean` |

★★★**`tilde` が充満忠実で、本質像が特徴づけられている**のが効く。

## ★★残る 2 点(第 64–)

| # | 主張 | 見通し |
|---|---|---|
| 1 | `tilde` はテンソル積を保つ(`tilde (M ⊗ N) ≅ tilde M ⊗ tilde N`) | ★★局所的には自明。★我々の `tensorModules` は層化を通すので**測ること** |
| 2 | 可逆層は `tilde` の本質像に入る | ★★局所自由 ⟹ 準連接 ⟹ `isIso_fromTildeΓ_iff` |

★★★1 は「決める前に測る」の対象である——
`tensorModules` は前層テンソル + 層化なので、`tilde` との交換は自明ではない。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★**`tilde` は関手である**(mathlib の在庫、機械確認)。 -/
noncomputable abbrev tildeFunctor : ModuleCat R ⥤ (Spec R).Modules :=
  AlgebraicGeometry.tilde.functor R

/-- ★★**`Γ ∘ tilde ≅ 恒等`**——`tilde` は充満忠実である(mathlib の在庫、機械確認)。

★★★これが `equivPicRing` の骨である。 -/
noncomputable abbrev tildeGammaIso :
    𝟭 (ModuleCat R) ≅ tildeFunctor R ⋙ AlgebraicGeometry.moduleSpecΓFunctor :=
  AlgebraicGeometry.tilde.toTildeΓNatIso

/-! ## ★出典の紐付け(`.src`) -/

def tildeGammaIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——アフィン比較の在庫)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
