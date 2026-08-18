import ABC3.Found.Arakelov.PicTildeUnit

/-!
# Arakelov (B1) 第 66 ブロック —— **`Γ` とテンソルの切断は `rfl` で読める**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★比較射を作るための足場

`equivPicRing` に要る比較射 `Γ(F) ⊗ Γ(G) → Γ(F ⊗ G)` を作る前に、
**両端の切断が何か**を確かめる(2026-08-18 実測、どちらも `rfl`):

| 主張 | 結果 |
|---|---|
| `Γ(F) = F.val(⊤)` | ★`rfl` |
| `(F ⊗ G).val(⊤) = (層化 (F.val ⊗ G.val)).val(⊤)` | ★★`rfl` |

★★★したがって比較射は

    F.val(⊤) ⊗_R G.val(⊤)  →  (F.val ⊗ G.val)(⊤)  →  (層化 (F.val ⊗ G.val))(⊤)

の 2 段である。★前段は係数環の読み替え(`R ≅ 𝒪(⊤)`)、後段は層化の unit。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `gammaSections` | ★`Γ(F)` は `F.val(⊤)`(`rfl`) |
| `tensorModulesSections` | ★★テンソルの切断は層化の切断(`rfl`) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (R : CommRingCat.{u})

/-- ★**`Γ(F)` は `F.val(⊤)` である**(`rfl`)。 -/
theorem gammaSections (F : (Spec R).Modules) :
    (AlgebraicGeometry.moduleSpecΓFunctor.obj F : Type u)
      = ((F.val.obj (op ⊤) : ModuleCat _) : Type u) := rfl

/-- ★★**テンソル積の切断は層化の切断である**(`rfl`)。 -/
theorem tensorModulesSections (F G : (Spec R).Modules) :
    ((tensorModules F G).val.obj (op ⊤) : Type u)
      = (((PresheafOfModules.sheafification (R := (Spec R).ringCatSheaf)
          (𝟙 (Spec R).ringCatSheaf.obj)).obj (F.val ⊗ G.val)).val.obj (op ⊤) : Type u) := rfl

/-! ## ★出典の紐付け(`.src`) -/

def gammaSections.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Γ の切断)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
