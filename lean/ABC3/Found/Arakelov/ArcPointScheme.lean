import ABC3.Found.Arakelov.ArcMetric

/-!
# Arakelov (C3) 第 247 ブロック —— **一点スキームでは局所自明 = 自明**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★ファイバーが 1 次元であることの根拠

`arcFiber p L` にノルムを入れるには、それが **1 次元 `ℂ` 加群**であることが要る。
★根拠は「`Spec ℂ` は**一点**である」ことだけである:

    S が ⊤ の被覆篩  ⟹  S は 𝟙 を含む      (点が入る開集合は ⊤ しか無い)

★★したがって `IsLocallyTrivial` は **⊤ での自明化**を直接与える
——`Spec ℂ` の上では「局所」と「大域」が一致する。

★★★mathlib の `instance {K} [Field K] : Unique (Spec (.of K))` を使う。

| 定理 | 内容 |
|---|---|
| `opens_eq_top` | ★空でない開集合は `⊤` |
| `coverSieve_top` | ★★被覆篩は `𝟙` を含む |
| `trivial_of_locallyTrivial` | ★★★**局所自明 ⟹ `⊤` で自明** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {Z : Scheme.{u}} [Unique Z]

/-- ★**一点スキームでは、点を含む開集合は `⊤` である**。 -/
theorem opens_eq_top {V : Z.Opens} (h : (default : Z) ∈ V) : V = ⊤ := by
  ext x
  simp only [Opens.coe_top, Set.mem_univ, iff_true]
  exact (Subsingleton.elim x default) ▸ h

/-- ★★**一点スキームでは、`⊤` の被覆篩は `𝟙` を含む**。 -/
theorem coverSieve_top {S : Sieve (⊤ : Z.Opens)}
    (hS : S ∈ (Opens.grothendieckTopology Z) ⊤) : S.arrows (𝟙 (⊤ : Z.Opens)) := by
  obtain ⟨V, i, hi, hxV⟩ := hS default trivial
  have hV : V = ⊤ := opens_eq_top hxV
  subst hV
  exact (Subsingleton.elim i (𝟙 _)) ▸ hi

/-- ★★★**一点スキームでは局所自明性は `⊤` での自明化を与える**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `arcFiber p L`(`p^* L` の `Spec ℂ` 上の大域切断)が
`ℂ` と同型であることが出る——ノルムを入れる場所が確定する。 -/
theorem trivial_of_locallyTrivial {M : Z.PresheafOfModules} (hM : IsLocallyTrivial Z M) :
    Nonempty ((restrictPresheafFunctor Z ⊤).obj M ≅ 𝟙_ (PresheafModulesOn Z ⊤)) := by
  obtain ⟨S, hS, hSt⟩ := hM ⊤
  exact hSt (𝟙 _) (coverSieve_top hS)

/-! ## ★出典の紐付け(`.src`) -/

def trivial_of_locallyTrivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——一点スキームでは局所自明は自明)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
