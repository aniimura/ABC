import ABC3.Found.Arakelov.ArcCoverPou
import ABC3.Found.Arakelov.PicLocalTrivial

/-!
# Arakelov (C3) 第 285 ブロック —— **★★★★★★★連続な計量は存在する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★(C3) の `metric_nonempty` が出た

局所自明な `F` に対し、`X^arc` がコンパクト Hausdorff なら
**連続な弧計量が存在する**。

## ★★機構——4 段

| 段 | 使うもの |
|---|---|
| 1 | `IsLocallyTrivial` を `⊤` に適用して**自明化篩**を取る |
| 2 | 篩の添字型 `trivIndex` を作り、選択公理で自明化を選ぶ |
| 3 | 第 284 `exists_pou_of_cover` —— 1 の分割を取る |
| 4 | 第 283 `gluedArcMetric` —— 貼り合わせる |

★★★`Opens` の Grothendieck 位相が**点ごと**であることが、段 3 の被覆条件を直接与える。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `trivIndex` | ★自明化篩の添字型 |
| `exists_contArcMetric` | ★★★★★★★**連続な計量の存在** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}}

/-- ★**自明化篩の添字型**——篩に属する開集合の全体。 -/
def trivIndex (S : Sieve (⊤ : X.Opens)) : Type :=
  { V : X.Opens // S.arrows (homOfLE (le_top : V ≤ ⊤)) }

/-- ★★★★★★★**連続な計量は存在する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★局所自明な可逆層の上には、`X^arc` がコンパクト Hausdorff であるかぎり
**4 法則すべてを満たす計量が存在する**。 -/
theorem exists_contArcMetric (F : X.Modules) (hF : IsLocallyTrivial X F.val)
    (hc : @CompactSpace (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X))
    (ht : @T2Space (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)) :
    Nonempty (ContArcMetric X F) := by
  classical
  obtain ⟨S, hSmem, htriv⟩ := hF ⊤
  set U : trivIndex S → X.Opens := fun i => i.1 with hUdef
  have e : ∀ i : trivIndex S,
      (restrictPresheafFunctor X (U i)).obj F.val ≅ 𝟙_ (PresheafModulesOn X (U i)) :=
    fun i => Classical.choice (htriv (homOfLE le_top) i.2)
  have hU : ∀ x : X, ∃ i, x ∈ U i := by
    intro x
    obtain ⟨V, g, hg, hx⟩ := hSmem x trivial
    exact ⟨⟨V, (Subsingleton.elim g (homOfLE (le_top : V ≤ ⊤))) ▸ hg⟩, hx⟩
  obtain ⟨f, hsub⟩ := exists_pou_of_cover U hU hc ht
  exact ⟨gluedArcMetric U f hsub F e (fun i => isOpen_arcOpenSet (U i))⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_contArcMetric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——連続な計量の存在)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
