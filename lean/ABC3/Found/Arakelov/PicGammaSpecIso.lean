import ABC3.Found.Arakelov.PicGammaSections

/-!
# Arakelov (B1) 第 67 ブロック —— **`Γ(Spec R, ⊤) ≅ R`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★比較射の係数環の読み替え

第 66 ブロックで比較射の形が

    F.val(⊤) ⊗_R G.val(⊤)  →  F.val(⊤) ⊗_{𝒪(⊤)} G.val(⊤)  →  (層化 ..)(⊤)

と分かった。★**前段は係数環の読み替え**であり、
mathlib の `Scheme.ΓSpecIso` がそれを与える(2026-08-18 実測):

    Scheme.ΓSpecIso R : Γ(Spec R, ⊤) ≅ R

★★**同型**なので、読み替えは情報を失わない
——`M ⊗_R N ≅ M ⊗_{𝒪(⊤)} N` になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaSpecIso` | ★`Γ(Spec R, ⊤) ≅ R`(mathlib の在庫、機械確認) |
| `specPresheafTopIso` | ★構造層の `⊤` での切断は `R` |

## ★★★残り(B1)

| # | 主張 |
|---|---|
| 1 | 比較射 `Γ(F) ⊗ Γ(G) → Γ(F ⊗ G)`(★第 66・67 で足場は揃った) |
| 2 | `tilde` がテンソルを保つ(★★第 64・65 + 第 29) |
| 3 | 可逆層は `tilde` の本質像に入る |
| 4 | `equivPicRing` / 5 `PicardData` の組み立て |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★**`Γ(Spec R, ⊤) ≅ R`**(mathlib の在庫、機械確認)。 -/
noncomputable abbrev gammaSpecIso : Γ(Spec R, ⊤) ≅ R := Scheme.ΓSpecIso R

/-- ★**構造層の `⊤` での切断は `R` である**。

★★これが比較射の係数環の読み替えを与える。 -/
noncomputable abbrev specPresheafTopIso :
    ((Spec R).presheaf.obj (op ⊤) : CommRingCat.{u}) ≅ R := Scheme.ΓSpecIso R

/-! ## ★出典の紐付け(`.src`) -/

def gammaSpecIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Γ(Spec R, ⊤) ≅ R)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
