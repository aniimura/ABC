import ABC3.Found.Arakelov.PicQuotient
import ABC3.Interface.Arakelov.LineBundle

/-!
# Arakelov (B1) 第 19 ブロック —— **`Found` と `Interface` の橋**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★第 18 ブロックで語彙を揃えた結果

`Interface/Arakelov/LineBundle.lean` の `IsInvertiblePre` / `preTensor` / `preUnit` /
`preRestrict` は、`Found/Arakelov/` の
`InvertiblePresheaf` / `⊗` / `𝟙_` / `restrictPresheafFunctor` と**同じもの**である。

★★★**したがって橋は定義の一致だけで架かる。**

| Interface | Found |
|---|---|
| `preTensor X F G` | `F ⊗ G`(第 1 ブロックのモノイダル構造) |
| `preUnit X` | `𝟙_ X.PresheafOfModules` |
| `preRestrict X V` | `restrictPresheafFunctor X V`(第 5 ブロック) |
| `IsInvertiblePre F` | `InvertiblePresheaf` の `isInv` + `trivial` |
-/

universe u

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.Arakelov

/-- ★★★★**可逆前層は `Interface` の `IsInvertiblePre` を満たす**。

★★★これが B1 の witness の中核部分である。
★語彙を揃えたので、証明は構造体のフィールドを取り出すだけである。 -/
theorem isInvertiblePre_of_invertiblePresheaf {X : Scheme.{0}}
    (L : InvertiblePresheaf X) : IsInvertiblePre L.carrier := by
  refine ⟨⟨L.inv, L.isInv⟩, ?_⟩
  intro U
  obtain ⟨S, hS, h⟩ := L.trivial U
  exact ⟨S, hS, h⟩

/-- ★★**逆に、`IsInvertiblePre` を満たす前層は可逆前層である**。

★★★これが `carrier_surjective`(`Pic` に足りない元が無いこと)の根拠になる。 -/
theorem exists_invertiblePresheaf_of_isInvertiblePre {X : Scheme.{0}}
    {F : X.PresheafOfModules} (h : IsInvertiblePre F) :
    ∃ L : InvertiblePresheaf X, Nonempty (L.carrier ≅ F) := by
  obtain ⟨⟨G, hG⟩, htriv⟩ := h
  refine ⟨{ carrier := F, inv := G, isInv := hG, trivial := ?_, invTrivial := ?_ },
    ⟨Iso.refl _⟩⟩
  · intro U
    obtain ⟨S, hS, hV⟩ := htriv U
    exact ⟨S, hS, hV⟩
  · -- ★逆 `G` の局所自明性は `F ⊗ G ≅ 𝟙_` と `F` の局所自明性から出る
    intro U
    obtain ⟨S, hS, hV⟩ := htriv U
    refine ⟨S, hS, ?_⟩
    intro V i hi
    obtain ⟨eF⟩ := hV i hi
    obtain ⟨e⟩ := hG
    -- G|_V ≅ 𝟙_ ⊗ G|_V ≅ (F|_V) ⊗ G|_V ≅ (F ⊗ G)|_V ≅ 𝟙_|_V ≅ 𝟙_
    exact ⟨(λ_ _).symm ≪≫ (eF.symm ⊗ᵢ Iso.refl _)
      ≪≫ restrictPresheafTensor F G
      ≪≫ (restrictPresheafFunctor X V).mapIso e
      ≪≫ (restrictPresheafUnit (X := X) (U := V)).symm⟩

/-! ## ★出典の紐付け(`.src`) -/

def isInvertiblePre_of_invertiblePresheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆前層が Interface の条件を満たすこと)",
    sectionId := "genell-def-1-1-i" }

def exists_invertiblePresheaf_of_isInvertiblePre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Interface の条件を満たす前層は可逆前層)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
