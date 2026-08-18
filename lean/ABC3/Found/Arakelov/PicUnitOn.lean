import ABC3.Found.Arakelov.PicFreeTerm

/-!
# Arakelov (B1) 第 56 ブロック —— **制限した site でも `(f|)^* 𝒪 ≅ 𝒪`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の輸送に要る最後の部品

第 54 ブロック(Beck–Chevalley)と合わせると

    (f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V) ≅ (f|)^*_pre 𝒪_V ≅ 𝒪_{f⁻¹V}

となり、★★**`f^*_pre L` が局所自明**であることが出る。

## ★★機構

| 部品 | 在庫 |
|---|---|
| `overPost` は終対象を終対象に送る | ★**`rfl`**(実測) |
| `free (yoneda 終対象) ≅ 𝟙_` | ★第 55 ブロック(一般形) |
| 制限した site の生成元の引き戻し | ★第 51 ブロック |

★★★3 つを繋ぐだけで**一発で通った**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `overPost_terminal` | ★`overPost` は終対象を終対象に送る(`rfl`) |
| `pullbackOnUnitIso` | ★★★★**`(f|)^*_pre 𝒪_V ≅ 𝒪_{f⁻¹V}`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-- ★**`overPost` は終対象を終対象に送る**(`rfl`)。 -/
theorem overPost_terminal :
    (overPost f V).obj (Over.mk (𝟙 V))
      = Over.mk (𝟙 ((Opens.map f.base).obj V)) := rfl

/-- ★★★★**制限した site でも `(f|)^*_pre 𝒪_V ≅ 𝒪_{f⁻¹V}`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「局所自明性が引き戻しで保たれる」ことの最後の部品である
——第 54 ブロック(Beck–Chevalley)と合わせて

    (f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V) ≅ (f|)^*_pre 𝒪_V ≅ 𝒪_{f⁻¹V}

となる。 -/
noncomputable def pullbackOnUnitIso :
    (pullbackPreOn f V).obj (𝟙_ (PresheafModulesOn Y V))
      ≅ 𝟙_ (PresheafModulesOn X ((Opens.map f.base).obj V)) :=
  (pullbackPreOn f V).mapIso
      (freeYonedaTermIso (R := (Over.forget V).op ⋙ Y.presheaf) (Over.mk (𝟙 V))
        (overTerminalUnique V)).symm
    ≪≫ pullbackOnFreeIsoGen f V (Over.mk (𝟙 V))
    ≪≫ freeYonedaTermIso
        (R := (Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
        (Over.mk (𝟙 ((Opens.map f.base).obj V)))
        (overTerminalUnique ((Opens.map f.base).obj V))

/-! ## ★出典の紐付け(`.src`) -/

def pullbackOnUnitIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した site でも f^* 𝒪 ≅ 𝒪)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
