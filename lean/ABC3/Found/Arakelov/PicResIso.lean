import ABC3.Found.Arakelov.PicGammaAppIso

/-!
# Arakelov (B2) 第 223 ブロック —— ★★★**等しい開集合の間の制限は同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-252 の逃げ道(b)が効いた

§9-252 で「転送が各段に伝播する」と測り、逃げ道を 3 つ挙げた:

| 候補 | 判定 |
|---|---|
| (a) 座標系を取り直す | 本命と書いた |
| ★(b) **制限が `IsIso`** | ★**これが効いた** |
| (c) `Iso` として組む | — |

★★**制限が同型**だと言えれば、第 222 の等式を**そのまま逆に解ける**
——転送を追う必要がない。

## ★★機構は `subst` 1 回

`W = V` を `subst` すると `homOfLE h : V ⟶ V` になり、
`Opens` は poset なので**`𝟙` と等しい**(`Subsingleton`)。あとは `simp`。

★★★§9-252 で「本命は (a)」と書いたが、**(b) のほうが 1 ブロックで済んだ**。
★逃げ道を**複数挙げておく**と、当たりを引ける。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `image_preimage_self` | ★★`fromSpec ''ᵁ (fromSpec ⁻¹ᵁ U) = U` |
| `isIso_res_of_eq` | ★★★**等しい開集合の間の制限は同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U : X.Opens} (hU : IsAffineOpen U)

/-- ★★**`fromSpec ''ᵁ (fromSpec ⁻¹ᵁ U) = U`**。 -/
theorem image_preimage_self :
    hU.fromSpec ''ᵁ (hU.fromSpec ⁻¹ᵁ U) = U := by
  rw [hU.fromSpec_preimage_self]
  exact fromSpec_image_top hU

/-- ★★★**等しい開集合の間の制限は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Opens` は poset なので、`subst` すれば `homOfLE h = 𝟙` である。 -/
theorem isIso_res_of_eq {V W : X.Opens} (h : V ≤ W) (heq : W = V) :
    IsIso (X.presheaf.map (homOfLE h).op) := by
  subst heq
  simp
  infer_instance

/-! ## ★出典の紐付け(`.src`) -/

def isIso_res_of_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——等しい開集合の間の制限は同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
