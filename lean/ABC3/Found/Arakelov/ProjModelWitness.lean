import ABC3.Found.Arakelov.UltraCompact
import ABC3.Found.Arakelov.ArcSpaceInterface

/-!
# Arakelov (C2) 第 91 ブロック —— **★★★★★★★★★C2 達成**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★層 C の律速が外れた

    X が ℤ 上固有・平坦 ⟹ X^arc はコンパクト

★これは `Definition 1.1` の「コンパクト正規複素解析空間 `X^arc`」の
**コンパクト性そのもの**である。

★★**Chow の補題は使っていない**——mathlib に無いからである。
代わりに **超フィルター + 付値判定法**(第 80–90)で作った。

## ★★★★2026-08-20 の訂正——`projectiveCase` は偽だった

以前は「連続かつ単射で像が閉有界」からコンパクト性を要求していたが、
★**連続単射の像がコンパクトでも定義域はコンパクトとは限らない**
(例: `[0,1) → 円周`)。★★**埋め込み**(`IsInducing`)が要る
——`Found/GenEll/ArcModel.lean` が実際に持っているのもそれである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `projectiveModelDataWitness` | ★`ProjectiveModelData` の実装 |
| `ProjectiveModelData.nonvacuous` | ★★★★★★★★★**C2 達成** |
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.Arakelov

/-- ★★★★★★★★★**`ProjectiveModelData` の実装**。 -/
noncomputable def projectiveModelDataWitness : ProjectiveModelData where
  toArcSpaceData := arcSpaceDataImpl
  ProperFlatOverZ := fun X => IsProper (specZIsTerminal.from X) ∧ Flat (specZIsTerminal.from X)
  properFlatOverZ_iff := fun _ => Iff.rfl
  properFlat_compact := fun X h => by
    have h1 := h.1
    rw [IsProper.eq_valuativeCriterion] at h1
    haveI : CompactSpace X :=
      (HasAffineProperty.iff_of_isAffine (P := @QuasiCompact)).1 h1.1.1.2
    exact compactSpace_arc h1.1.1.1
  projectiveCase := fun X n emb hind hclosed hbdd => by
    letI : TopologicalSpace (arcSpaceDataImpl.Arc X) := arcSpaceDataImpl.topology X
    rw [← isCompact_univ_iff, hind.isCompact_iff, Set.image_univ]
    exact Metric.isCompact_of_isClosed_isBounded hclosed hbdd

/-- ★★★★★★★★★**`ProjectiveModelData` は非空虚である**——C2 達成。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが Arakelov 理論の 9 件のうち **C2** である——**層 C の律速**であった。 -/
theorem ProjectiveModelData.nonvacuous : Nonempty ProjectiveModelData :=
  ⟨projectiveModelDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def projectiveModelDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——ℤ-固有からのコンパクト性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.Arakelov
