import ABC3.Found.Arakelov.PicInvPull

/-!
# Arakelov (B1) 第 62 ブロック —— ★★★★★★★★★★**`Pic` の引き戻しと 3 公理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★取れたもの(sorry 0)

    picPullback     : PicSheaf Y → PicSheaf X
    picPullback_id  : picPullback (𝟙 X) L = L
    picPullback_comp: picPullback (f ≫ g) L = picPullback f (picPullback g L)
    picPullback_mul : picPullback f (L * M) = picPullback f L * picPullback f M

★★★**`PicardData` の `pullback` 系 4 欄がすべて埋まった。**

## ★★★★なぜ短いか —— `PicSheaf` は同型による商だから

`PicSheaf` は可逆層の**同型類**である。したがって
**対象の同型がそのまま等式になる**(`mk_eq_mk`)。

| 公理 | 使う同型 | 出どころ |
|---|---|---|
| `id` | `pullback (𝟙 X) ≅ 𝟭` | ★mathlib `Scheme.Modules.pullbackId` |
| `comp` | `pullback g ⋙ pullback f ≅ pullback (f ≫ g)` | ★mathlib `pullbackComp` |
| `mul` | `f^*(L ⊗ M) ≅ f^*L ⊗ f^*M` | ★★★★**第 60 ブロック**(我々が 40 ブロックで作った) |

★★`id` と `comp` は mathlib から**ただで**出る。
★★★`mul` だけが山であり、それが第 18–60 の 43 ブロックであった。

## ★★残り(B1)

    `equivPicRing`(アフィンでの `CommRing.Pic` との比較)と `PicardData` の組み立て。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y Z : Scheme.{u}} (f : X ⟶ Y) (g : Y ⟶ Z)

/-- ★★★★**`Pic` の引き戻し**——可逆層の引き戻し(第 61)を同型類へ降ろす。 -/
noncomputable def picPullback : PicSheaf Y → PicSheaf X :=
  Quotient.map (InvSheaf.pullback f)
    (by rintro L L' ⟨e⟩; exact ⟨(Scheme.Modules.pullback f).mapIso e⟩)

/-- ★★**恒等射の引き戻しは何もしない**。 -/
theorem picPullback_id (L : PicSheaf X) : picPullback (𝟙 X) L = L := by
  refine Quotient.inductionOn L ?_
  intro M
  exact PicSheaf.mk_eq_mk ((Scheme.Modules.pullbackId X).app M.carrier)

/-- ★★**合成の引き戻しは引き戻しの合成**。 -/
theorem picPullback_comp (L : PicSheaf Z) :
    picPullback (f ≫ g) L = picPullback f (picPullback g L) := by
  refine Quotient.inductionOn L ?_
  intro M
  exact PicSheaf.mk_eq_mk (((Scheme.Modules.pullbackComp f g).app M.carrier).symm)

/-- ★★★★★★★★★★**引き戻しは群準同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが第 18–60 の 43 ブロックの**目的**であった——
`f^*(L ⊗ M) ≅ f^*L ⊗ f^*M`(第 60)がそのまま等式になる。 -/
theorem picPullback_mul (L M : PicSheaf Y) :
    picPullback f (L * M) = picPullback f L * picPullback f M := by
  refine Quotient.inductionOn₂ L M ?_
  intro A B
  exact PicSheaf.mk_eq_mk
    (pullbackSheafTensorIso' f A.carrier B.carrier A.trivial B.trivial)

/-! ## ★出典の紐付け(`.src`) -/

def picPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Pic の引き戻し)",
    sectionId := "genell-def-1-1-i" }

def picPullback_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——引き戻しが群準同型であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
