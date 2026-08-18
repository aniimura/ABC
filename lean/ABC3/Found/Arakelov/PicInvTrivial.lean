import ABC3.Found.Arakelov.PicSectionTrivial

/-!
# Arakelov (B1) 第 72 ブロック —— **逆もまた局所自明である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`sheafOf_surjective` に要る最後の 1 点

`Interface` の `IsInvertibleSheaf F` が言うのは

| 条 | 内容 |
|---|---|
| (a) | `∃ G, F ⊗ G ≅ 𝒪_X` |
| (b) | `F` は局所自明 |

★★一方 `Found` の `InvSheaf` は**逆の側の局所自明性** `invTrivial` も持つ。
★★★したがって `IsInvertibleSheaf F` から `InvSheaf X` を作るには

    F ⊗ G ≅ 𝒪_X  かつ  F 局所自明   ⟹   G 局所自明

を**証明せねばならない**。本ブロックがそれである。

## ★★機構 —— 「G は層である」が効く

`F` が自明な `V` の上で、**前層の**テンソルは

    (F.val ⊗ G.val)|_V ≅ F.val|_V ⊗ G.val|_V ≅ 𝟙_ ⊗ G.val|_V ≅ G.val|_V

となる(第 14 ブロック `restrictPresheafTensor`)。

★★★**右辺は層である**(`G : X.Modules` なので `G.val` は層、
層の制限は層——第 17 ブロック `isSheaf_restrict`)。
★したがって左辺も層であり、**層化の単位が `V` の上で同型になる**
(第 17 ブロック `isIso_restrictMap_sheafifyUnit`)。

★★これで層化を跨げるので

    G.val|_V ≅ (F.val ⊗ G.val)|_V ≅ (F ⊗ G).val|_V ≅ 𝒪|_V ≅ 𝟙_

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isLocallyTrivial_of_tensor_unit` | ★★★★**逆も局所自明** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (X : Scheme.{u})

/-- ★★★★**テンソル積が構造層になる相手も局所自明である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `Interface` の `IsInvertibleSheaf` から `Found` の `InvSheaf` を
作れるようになる——`sheafOf_surjective` に要る最後の 1 点である。 -/
theorem isLocallyTrivial_of_tensor_unit (F G : X.Modules)
    (e : tensorModules F G ≅ unitModules X)
    (hF : IsLocallyTrivial X F.val) :
    IsLocallyTrivial X G.val := by
  intro U
  obtain ⟨S, hS, hV⟩ := hF U
  refine ⟨S, hS, fun V i hi => ?_⟩
  obtain ⟨eF⟩ := hV i hi
  -- ★前層のテンソルは制限を跨ぐ——`(F ⊗ G)|_V ≅ G|_V`
  have eT : (restrictPresheafFunctor X V).obj (F.val ⊗ G.val)
      ≅ (restrictPresheafFunctor X V).obj G.val :=
    (restrictPresheafTensor F.val G.val).symm
      ≪≫ (eF ⊗ᵢ Iso.refl ((restrictPresheafFunctor X V).obj G.val)) ≪≫ λ_ _
  -- ★★右辺は層なので、左辺も層である
  have hG : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ G.val.presheaf) := isSheaf_restrict X V _ G.isSheaf
  have hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((Over.forget V).op ⋙ (F.val ⊗ G.val).presheaf) :=
    (Presheaf.isSheaf_of_iso_iff ((PresheafOfModules.toPresheaf _).mapIso eT)).2 hG
  -- ★★★したがって層化の単位が `V` の上で同型になる
  haveI := isIso_restrictMap_sheafifyUnit X V (F.val ⊗ G.val) hP
  exact ⟨eT.symm
    ≪≫ asIso ((restrictPresheafFunctor X V).map (sheafifyUnit X (F.val ⊗ G.val)))
    ≪≫ (restrictPresheafFunctor X V).mapIso
        ((SheafOfModules.forget X.ringCatSheaf).mapIso e)
    ≪≫ (restrictPresheafUnit (X := X) (U := V)).symm⟩

/-! ## ★出典の紐付け(`.src`) -/

def isLocallyTrivial_of_tensor_unit.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——逆もまた局所自明であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
