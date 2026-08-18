import ABC3.Found.Arakelov.PicAppLESec
import ABC3.Found.Arakelov.PicSectionPull

/-!
# Arakelov (B2) 第 214 ブロック —— **可換正方形を元に当てる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 213 を**元のレベル**に落とす

    (Spec.map appLE).appTop (hB.fromSpec.appTop v) = hA.fromSpec.appTop (f.appTop v)

★これで「切断を `Spec` 側へ運んでから `appLE` を当てる」と
「先に引き戻してから `Spec` 側へ運ぶ」が、**元ごとに一致する**と言える。

## ★★測って分かった —— `appTop` は `Γ(Y,⊤)` から始まる

    Scheme.Hom.appTop (f : X ⟶ Y) : Γ(Y, ⊤) ⟶ Γ(X, ⊤)

★したがって本ブロックの `v` は **`Γ(Y, ⊤)` の元**である。
★★生成元の場合で要るのは `Γ(Y, B')`(`B' = fromSpec ''ᵁ ⊤`)から始まる形なので、
`appTop` を `制限 ∘ appIso` に分解する 1 段が**まだ要る**——
`hB.fromSpec.appTop = (Y の制限 ⊤ → B') ≫ (hB.fromSpec.appIso ⊤).hom` である。

★★★これが B2 の**本当に最後の 1 段**である(見積もり 1–2 ブロック)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appLE_fromSpec_apply` | ★★★**可換正方形を元に当てた形** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★**可換正方形を元に当てる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 213 の等式に元 `v` を当てただけである。 -/
theorem appLE_fromSpec_apply {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (v : (Γ(Y, ⊤) : Type u)) :
    ((Spec.map (f.appLE B A i)).appTop).hom ((hB.fromSpec.appTop).hom v)
      = (hA.fromSpec.appTop).hom ((f.appTop).hom v) := by
  have h := appLE_fromSpec_sections f hA hB i
  rw [lhs_decomp f hB i, rhs_decomp f hA] at h
  exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) v) h

/-! ## ★出典の紐付け(`.src`) -/

def appLE_fromSpec_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可換正方形を元に当てる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
