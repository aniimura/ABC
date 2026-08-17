import ABC3.Found.Arakelov.PicType

/-!
# Arakelov (B1) 第 15 ブロック —— **強い局所自明性とそのテンソル閉性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★なぜ「強い」形が要るか

第 9 ブロック(局所単射性)には **`V` の上での基底 1 本**で足りた
(`IsLocallyRankOne`)。★★しかし**テンソル積で閉じること**を示すには、
それでは足りない——`(A ⊗ B)(V)` の基底を作るのに `A(V)`, `B(V)` の基底だけでは
「制限と両立する」ことが言えないからである。

★★★そこで **`Over V` 上の前層同型** `A|_V ≅ 𝒪|_V` を要求する形(`IsLocallyTrivial`)を置く。
これは可逆層(直線束)の**標準の定義**そのものである。

## ★★機構——第 5 ブロックがそのまま効く

    (A ⊗ B)|_V ≅ A|_V ⊗ B|_V ≅ 𝟙_ ⊗ 𝟙_ ≅ 𝟙_

★1 つ目は `restrictPresheafTensor`(第 5 ブロック、mathlib の
`pushforward₀OfCommRingCat.Monoidal` から出たもの)、
★2 つ目は仮定、3 つ目は単位律。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite

variable (X : Scheme.{u})

/-- ★★★**強い局所自明性** —— 各開集合が「`𝒪` と前層同型になる開集合」で覆われる。

★★これが可逆層(直線束)の標準の定義である。 -/
def IsLocallyTrivial (M : X.PresheafOfModules) : Prop :=
  ∀ U : X.Opens, ∃ S : Sieve U, S ∈ (Opens.grothendieckTopology X) U ∧
    ∀ ⦃V : X.Opens⦄ (i : V ⟶ U), S i →
      Nonempty ((restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))

variable {X}

/-- ★★★★**強い局所自明性はテンソル積で閉じる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `Pic X` の**乗法が閉じている**ことの根拠である。

★機構は第 5 ブロック(制限はモノイダル関手)+ 被覆篩の交わり。 -/
theorem IsLocallyTrivial.tensor {A B : X.PresheafOfModules}
    (hA : IsLocallyTrivial X A) (hB : IsLocallyTrivial X B) :
    IsLocallyTrivial X (A ⊗ B) := by
  intro U
  obtain ⟨SA, hSA, hA'⟩ := hA U
  obtain ⟨SB, hSB, hB'⟩ := hB U
  refine ⟨SA ⊓ SB, (Opens.grothendieckTopology X).intersection_covering hSA hSB, ?_⟩
  intro V i hi
  obtain ⟨eA⟩ := hA' i hi.1
  obtain ⟨eB⟩ := hB' i hi.2
  exact ⟨(restrictPresheafTensor A B).symm ≪≫
    (eA ⊗ᵢ eB) ≪≫ λ_ (𝟙_ (PresheafModulesOn X V))⟩

/-- ★**強い局所自明性は弱い形(各点での基底)を含意する**。

★★`Over V` の終対象 `⟨V, 𝟙⟩` で評価するだけである。 -/
theorem IsLocallyTrivial.isLocallyRankOne {M : X.PresheafOfModules}
    (h : IsLocallyTrivial X M) : IsLocallyRankOne X M := by
  intro U
  obtain ⟨S, hS, hV⟩ := h U
  refine ⟨S, hS, ?_⟩
  intro V i hi
  obtain ⟨e⟩ := hV i hi
  exact ⟨(((PresheafOfModules.evaluation _ (op (Over.mk (𝟙 V)))).mapIso e).toLinearEquiv).symm⟩

/-! ## ★出典の紐付け(`.src`) -/

def IsLocallyTrivial.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆層の強い局所自明性)",
    sectionId := "genell-def-1-1-i" }

def IsLocallyTrivial.tensor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所自明性がテンソル積で閉じること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
