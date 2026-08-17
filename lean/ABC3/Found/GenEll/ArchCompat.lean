import ABC3.Found.GenEll.PullbackNatural

/-!
# [GenEll] Definition 1.2, (i) —— **アルキメデス側の仮定を discharge する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

## ★★★最後の非設計要素

`ArchADivBase.lean` の `archADiv_baseChange` は
「埋め込みが制限と両立する(共役を除く)」を仮定として受けていた。
★★★**本ファイルでそれを定理にする。**

## ★★機構は 2 つ

1. **タワー** —— `𝓞_F → 𝓞_K → K` と `𝓞_F → F → K` は同じ射である
   (`IsScalarTower`。★mathlib に instance がある——実測)
2. **共役の場合分け** —— `ArchBaseChange.lean` の `embedding_comap_dichotomy`

★★★これで `Definition 1.2, (i)` に残るのは **`X(ℚ̄)` の型の設計だけ**になる。
-/

namespace ABC3.Found.GenEll

open NumberField

variable (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K]
  [Algebra F K]

/-! ## ★★★埋め込みは制限と両立する(共役を除く) -/

/-- ★★★**`archRingHom` は底変換と両立する**(共役を除く)。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle M.

    `(archRingHom K w) ∘ (𝓞_F → 𝓞_K) = archRingHom F (w|_F)`  または その共役

★★★**`ArchADivBase.lean` の仮定が、これで定理になる。**

★機構は 2 つ:
`IsScalarTower`(`𝓞_F → 𝓞_K → K` と `𝓞_F → F → K` は同じ射)と
`embedding_comap_dichotomy`(制限した素点の埋め込みは合成かその共役)。 -/
theorem archRingHom_compat (w : InfinitePlace K) :
    (archRingHom K w).comp (algebraMap (𝓞 F) (𝓞 K))
        = archRingHom F (w.comap (algebraMap F K))
      ∨ (archRingHom K w).comp (algebraMap (𝓞 F) (𝓞 K))
        = (starRingEnd ℂ).comp (archRingHom F (w.comap (algebraMap F K))) := by
  -- ★タワー: `𝓞_F → 𝓞_K → K` = `𝓞_F → F → K`
  have htower : ∀ x : 𝓞 F,
      algebraMap (𝓞 K) K (algebraMap (𝓞 F) (𝓞 K) x)
        = algebraMap F K (algebraMap (𝓞 F) F x) := by
    intro x
    rw [← IsScalarTower.algebraMap_apply, ← IsScalarTower.algebraMap_apply]
  rcases embedding_comap_dichotomy K w (algebraMap F K) with h | h
  · left
    ext x
    show w.embedding (algebraMap (𝓞 K) K (algebraMap (𝓞 F) (𝓞 K) x))
      = (w.comap (algebraMap F K)).embedding (algebraMap (𝓞 F) F x)
    rw [htower, h]
    rfl
  · right
    ext x
    show w.embedding (algebraMap (𝓞 K) K (algebraMap (𝓞 F) (𝓞 K) x))
      = (starRingEnd ℂ) ((w.comap (algebraMap F K)).embedding (algebraMap (𝓞 F) F x))
    rw [htower, h]
    simp [NumberField.ComplexEmbedding.conjugate]

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Definition 1.2, (i)` 全体には
`X(ℚ̄)` の型そのものの構成(数体についての colimit)が残っている——
★それは**設計の仕事であり、証明の仕事ではない**。 -/

def archRingHom_compat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(アルキメデス側の底変換——仮定でなく定理として)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.GenEll
