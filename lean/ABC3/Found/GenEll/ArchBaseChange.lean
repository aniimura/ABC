import ABC3.Found.GenEll.ArchConj

/-!
# [GenEll] Definition 1.2, (i) —— **アルキメデス側の底変換**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★共役の場合分けを閉じる

`ArchConj.lean` は「埋め込みが制限で一致すれば ℂ-点も一致する」を取った。
★実際には **一致するとは限らない**——`embedding_mk_eq` により
`v.embedding` は `w.embedding|_F` **または その複素共役**である。

★★★**本ファイルは共役の場合を閉じる。** そこで `IsConjInvariant`
(原文の「計量は `ι_X` と両立する」)が効く。

## ★★機構

複素共役でずれた埋め込みの `Spec` は、`conjSpec` との合成である:

    `Spec.map (conj ∘ ψ) = conjSpec ≫ Spec.map ψ`

★`Spec` は反変なので合成の順が入れ替わる(`Spec.map_comp`)。
あとは `IsConjInvariant` で値が変わらないことを言うだけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]

/-! ## ★★制限した素点の埋め込みは共役を除いて決まる -/

/-- ★**制限した素点は、埋め込みの合成が定める素点である**。 -/
theorem comap_eq_mk {k : Type*} [Field k] (w : InfinitePlace K) (f : k →+* K) :
    w.comap f = InfinitePlace.mk (w.embedding.comp f) := by
  conv_lhs => rw [← InfinitePlace.mk_embedding w]
  exact InfinitePlace.comap_mk _ _

/-- ★★**制限した素点の埋め込みは、合成した埋め込み**または**その共役**である。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これが「計量が `ι_X` と両立する」ことを要求する理由である。**
共役の分だけ ℂ-点がずれるので、Green 関数が共役不変でないと
高さが定義体の取り方に依ってしまう。 -/
theorem embedding_comap_dichotomy {k : Type*} [Field k] (w : InfinitePlace K)
    (f : k →+* K) :
    (w.comap f).embedding = w.embedding.comp f
      ∨ (w.comap f).embedding
          = NumberField.ComplexEmbedding.conjugate (w.embedding.comp f) := by
  rw [comap_eq_mk]
  exact InfinitePlace.embedding_mk_eq _

/-! ## ★★共役でずれた ℂ-点 -/

variable {F K}

/-- ★★**共役でずれた埋め込みの `Spec` は `conjSpec` との合成である**。

★`Spec` は反変なので合成の順が入れ替わる。 -/
theorem specMap_conj_comp {R : Type} [CommRing R] (ψ : R →+* ℂ) :
    Spec.map (CommRingCat.ofHom ((starRingEnd ℂ).comp ψ))
      = conjSpec ≫ Spec.map (CommRingCat.ofHom ψ) := by
  rw [conjSpec, ← Spec.map_comp]
  rfl

/-- ★★★**共役でずれても Green 関数の値は変わらない**。

原文 (GenEll p.3):
> There is an evident notion of tensor product of arithmetic line bundles on X. The

★★★これが `IsConjInvariant`(原文の `ι_X` 両立)の効き所である。 -/
theorem green_specMap_conj {X : Scheme.{0}} {R : Type} [CommRing R] (ψ : R →+* ℂ)
    (p : Spec (CommRingCat.of R) ⟶ X) (g : GreenFn X) (hg : IsConjInvariant g) :
    g (Spec.map (CommRingCat.ofHom ((starRingEnd ℂ).comp ψ)) ≫ p)
      = g (Spec.map (CommRingCat.ofHom ψ) ≫ p) := by
  rw [specMap_conj_comp, Category.assoc]
  exact hg _

/-! ## ★★★アルキメデス側の底変換 -/

variable (F K)

/-- ★★★**Green 関数が共役不変なら、ℂ-点での値は底変換で不変**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

★★★**これがアルキメデス側の底変換不変性である。**
仮定は 2 つだけ:
- `IsConjInvariant g`(原文の `ι_X` 両立)
- 環準同型の両立(`archRingHom` が制限と合う——共役を除いて)

★★共役でずれる場合を `IsConjInvariant` が吸収する。 -/
theorem green_archPoint_baseChange {X : Scheme.{0}} (v : InfinitePlace F)
    (w : InfinitePlace K) (xF : specRingOfIntegers F ⟶ X)
    (φ : CommRingCat.of (𝓞 F) ⟶ CommRingCat.of (𝓞 K))
    (g : GreenFn X) (hg : IsConjInvariant g)
    (hφ : (archRingHom K w).comp φ.hom = archRingHom F v
      ∨ (archRingHom K w).comp φ.hom = (starRingEnd ℂ).comp (archRingHom F v)) :
    g (archPoint (Spec.map φ ≫ xF) w) = g (archPoint xF v) := by
  have hcomp : archPoint (Spec.map φ ≫ xF) w
      = Spec.map (CommRingCat.ofHom ((archRingHom K w).comp φ.hom)) ≫ xF := by
    rw [archPoint, ← Category.assoc, archSpecMap_comp]
    rfl
  rcases hφ with h | h
  · rw [hcomp, h, archPoint, archSpecMap]
  · rw [hcomp, h]
    exact green_specMap_conj (archRingHom F v) xF g hg

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
`degNormalized_baseChange` との接続と `X(ℚ̄)` の型は `AlgPointClass.lean`(§9-744)で入った。 -/

def green_archPoint_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(アルキメデス側の底変換不変性——ℂ-点での値のみ)",
    sectionId := "genell-def-1-2-i" }

def embedding_comap_dichotomy.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1(計量と ι_X の両立を要求する理由——共役の場合分け)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.GenEll
