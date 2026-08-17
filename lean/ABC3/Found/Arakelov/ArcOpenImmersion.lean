import ABC3.Found.Arakelov.ArcTopologyAffineEq

/-!
# Arakelov (C1) の第七段 —— **開埋め込みに沿った連続性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★C1 に残る最後の条件

`Interface/Arakelov/ArcSpace.lean` の C1 は

    topology X = induced (map f) (topology Y)      (f : X ⟶ Y が開埋め込み)

を要求する。★これは 2 向きある:

| 向き | 意味 | 状態 |
|---|---|---|
| `topology X ≤ induced` | ★**`(· ≫ f)` が連続** | ★★★**本ファイル** |
| `induced ≤ topology X` | 開埋め込みの像で部分空間位相になること | ★未取得 |

## ★★★機構: 開埋め込みはアフィン開をアフィン開に送る

`f` が開埋め込みで `U ⊆ X` がアフィン開なら `f ''ᵁ U ⊆ Y` もアフィン開
(`IsAffineOpen.image_of_isOpenImmersion`、mathlib)。★★さらに

    U.ι ≫ f = e.hom ≫ (f ''ᵁ U).ι

という分解がある(`IsOpenImmersion.isoOfRangeEq_hom_fac`)。
★★★**したがって `X` の chart は `Y` の chart 1 つに収まる**——ここが要点である。
一般の射だと像が複数の chart にまたがるので、この議論は開埋め込みに固有である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace

/-! ## ★★アフィン開の像への同型 -/

/-- ★**開埋め込みに沿ったアフィン開の像**。 -/
def imageAffineOpen {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) : Y.affineOpens :=
  ⟨f ''ᵁ U.1, U.2.image_of_isOpenImmersion f⟩

/-- ★★**`U` とその像は同型**であり、包含と両立する。 -/
noncomputable def imageAffineOpenIso {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) : U.1.toScheme ≅ (imageAffineOpen f U).1.toScheme :=
  IsOpenImmersion.isoOfRangeEq (U.1.ι ≫ f) (imageAffineOpen f U).1.ι (by
    show Set.range ((U.1.ι ≫ f).base) = _
    rw [Scheme.Opens.range_ι]
    have : ((U.1.ι ≫ f).base : _ → _) = (f.base : _ → _) ∘ (U.1.ι.base : _ → _) := rfl
    rw [this, Set.range_comp, Scheme.Opens.range_ι]
    exact (Scheme.Hom.coe_image (f := f) (U := U.1)).symm)

/-- ★★**分解**: `U.ι ≫ f = iso.hom ≫ (f ''ᵁ U).ι`。 -/
theorem imageAffineOpenIso_fac {X Y : Scheme.{0}} (f : X ⟶ Y) [IsOpenImmersion f]
    (U : X.affineOpens) :
    (imageAffineOpenIso f U).hom ≫ (imageAffineOpen f U).1.ι = U.1.ι ≫ f :=
  IsOpenImmersion.isoOfRangeEq_hom_fac _ _ _

/-! ## ★出典の紐付け(`.src`) -/

def imageAffineOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みがアフィン開を保つこと)",
    sectionId := "genell-def-1-1-i" }

def imageAffineOpenIso_fac.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みの chart の分解)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
