import ABC3.Found.Arakelov.PicPullTensorFinal

/-!
# Arakelov (B1) 第 61 ブロック —— ★★★★★★★★**可逆層の引き戻し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★可逆層は引き戻せる

    InvSheaf Y  →  InvSheaf X

★★4 つの欄すべてが揃った:

| 欄 | 出どころ |
|---|---|
| `carrier` / `inv` | ★mathlib の `Scheme.Modules.pullback` |
| `isInv`(`f^*L ⊗ f^*L' ≅ 𝒪`) | ★★★★第 60(仮定なしのテンソル保存)+ 第 21(`f^* 𝒪 ≅ 𝒪`) |
| `trivial` / `invTrivial` | ★★★第 59(局所自明性の輸送)+ 第 41(`pullbackValIso`) |

## ★★本ブロックで取れるもの

| 定理・定義 | 内容 |
|---|---|
| `isLocallyTrivial_of_iso` | ★局所自明性は同型で移る |
| `isLocallyTrivial_pullbackModules` | ★★層の引き戻しの局所自明性 |
| `InvSheaf.pullback` | ★★★★★★**可逆層の引き戻し** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★**局所自明性は同型で移る**。 -/
theorem isLocallyTrivial_of_iso {M N : X.PresheafOfModules} (e : M ≅ N)
    (h : IsLocallyTrivial X M) : IsLocallyTrivial X N := by
  intro U
  obtain ⟨S, hS, hV⟩ := h U
  refine ⟨S, hS, ?_⟩
  intro V i hi
  obtain ⟨g⟩ := hV i hi
  exact ⟨((restrictPresheafFunctor X V).mapIso e).symm ≪≫ g⟩

/-- ★★**層の引き戻しは局所自明性を保つ**。

★機構は第 41 ブロックの `pullbackValIso`(層の引き戻し = 前層で引き戻して層化)
+ 第 59(前層の輸送)+ 第 16(層化は自明性を保つ)。 -/
theorem isLocallyTrivial_pullbackModules (L : Y.Modules)
    (hL : IsLocallyTrivial Y L.val) :
    IsLocallyTrivial X ((Scheme.Modules.pullback f).obj L).val :=
  isLocallyTrivial_of_iso
    ((SheafOfModules.forget X.ringCatSheaf).mapIso (pullbackValIso f L)).symm
    (isLocallyTrivial_sheafify X ((pullbackPre f).obj L.val)
      (isLocallyTrivial_pullbackPre f L.val hL))

/-- ★★★★★★**可逆層の引き戻し**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `PicSheaf` の `pullback` が定義できる。 -/
noncomputable def InvSheaf.pullback (L : InvSheaf Y) : InvSheaf X where
  carrier := (Scheme.Modules.pullback f).obj L.carrier
  inv := (Scheme.Modules.pullback f).obj L.inv
  isInv := ⟨(pullbackSheafTensorIso' f L.carrier L.inv L.trivial L.invTrivial).symm
    ≪≫ (Scheme.Modules.pullback f).mapIso L.isInv.some
    ≪≫ pullbackUnitIso f⟩
  trivial := isLocallyTrivial_pullbackModules f L.carrier L.trivial
  invTrivial := isLocallyTrivial_pullbackModules f L.inv L.invTrivial

/-! ## ★出典の紐付け(`.src`) -/

def InvSheaf.pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層の引き戻し)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
