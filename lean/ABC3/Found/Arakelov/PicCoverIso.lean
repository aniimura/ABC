import ABC3.Found.Arakelov.PicLocalBij

/-!
# Arakelov (B1) 第 102 ブロック —— **被覆で全単射なら同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★任意の site で使える器具

    各対象に被覆篩があって、その上で `f.app` が全単射  ⟹  `f` は同型

★★**空間の上とは限らない**——`Over V` の `over` 位相でも使える。
これが §9-105 で特定した構造上の 1 点を解く。

## ★★機構 —— 3 段

| 段 | 出典 |
|---|---|
| 被覆 → 局所全射 | ★第 101 ブロック |
| 被覆 → 局所単射 | ★第 101 ブロック(自然性で移す) |
| 局所全単射 → 同型 | ★mathlib `Sheaf.isLocallyBijective_iff_isIso` |

★★★§9-73 で「使えない」と止めた器具が、**`Sheaf J Ab` の版なら通る**
(§9-107 で切り分けた)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_of_bijective_on_cover` | ★★★★★★**被覆で全単射なら同型** |
-/

universe v u

namespace ABC3.Found.Arakelov

open CategoryTheory Opposite

variable {C : Type u} [Category.{v} C] (J : GrothendieckTopology C)
  [HasWeakSheafify J AddCommGrpCat.{v}] [J.WEqualsLocallyBijective AddCommGrpCat.{v}]

/-- ★★★★★★**被覆の上で全単射なら層の射は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★任意の site で使える——`Over V` の `over` 位相でも当たる。 -/
theorem isIso_of_bijective_on_cover {F G : Sheaf J AddCommGrpCat.{v}} (f : F ⟶ G)
    (h : ∀ U : C, ∃ S : Sieve U, S ∈ J U ∧
      ∀ ⦃W : C⦄ (i : W ⟶ U), S.arrows i → Function.Bijective (f.hom.app (op W))) :
    IsIso f := by
  haveI : Presheaf.IsLocallySurjective J f.hom := by
    refine isLocallySurjective_of_cover J f.hom (fun U s => ?_)
    obtain ⟨S, hS, hb⟩ := h U
    exact ⟨S, hS, fun {W} i hi => (hb i hi).2 _⟩
  haveI : Presheaf.IsLocallyInjective J f.hom := by
    refine isLocallyInjective_of_coverSieve J f.hom (fun U x y hxy => ?_)
    obtain ⟨S, hS, hb⟩ := h U
    refine ⟨S, hS, fun {W} i hi => ?_⟩
    refine (hb i hi).1 ?_
    rw [NatTrans.naturality_apply f.hom i.op x, NatTrans.naturality_apply f.hom i.op y, hxy]
  exact (Sheaf.isLocallyBijective_iff_isIso f).1 ⟨inferInstance, inferInstance⟩

/-! ## ★出典の紐付け(`.src`) -/

def isIso_of_bijective_on_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——被覆で全単射なら同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
