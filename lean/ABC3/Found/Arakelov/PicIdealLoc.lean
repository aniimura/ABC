import ABC3.Found.Arakelov.PicIdealBO

/-!
# Arakelov (B2) 第 153 ブロック —— **切断は元のイデアルの局所化**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★mathlib に `Submodule.toLocalized'` があった

    Submodule.localized' S p f M' : Submodule S N
    instance : IsLocalizedModule p (M'.toLocalized' S p f)

★これと `IsAffineOpen.isLocalization_basicOpen`(制限は `Away h` の局所化)を繋ぐと

    idealSections D (X.basicOpen h) = (D.ideal A).localized' Γ(X, D(h)) (powers h) …

が出る。★★`localized'_eq_span` と `Ideal.map = span (f '' I)` が**同じもの**だからである。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `resAlg` | ★制限による代数構造(`@[reducible]`) |
| `isLocalizationAway_bo` | ★★mathlib の `Away` 局所化 |
| `idealSections_eq_localized` | ★★★**切断は元のイデアルの局所化** |

★★★これで「`D.ideal A` が可逆 ⟹ 局所自由」が第 92・130 の再利用で書ける。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

/-- ★制限による代数構造。 -/
@[reducible] noncomputable def resAlg (h : (Γ(X, A.1) : Type u)) :
    Algebra (Γ(X, A.1) : Type u) (Γ(X, X.basicOpen h) : Type u) :=
  (X.presheaf.map (homOfLE (X.basicOpen_le h)).op).hom.toAlgebra

/-- ★★制限は `Away h` の局所化である(mathlib)。 -/
theorem isLocalizationAway_bo (h : (Γ(X, A.1) : Type u)) :
    @IsLocalization.Away _ _ h (Γ(X, X.basicOpen h) : Type u) _ (resAlg A h) :=
  A.2.isLocalization_basicOpen h

/-- ★★★切断は元のイデアルの局所化である。 -/
theorem idealSections_eq_localized (h : (Γ(X, A.1) : Type u)) :
    letI := resAlg A h
    letI := isLocalizationAway_bo A h
    idealSections D (X.basicOpen h)
      = (D.ideal A).localized' (Γ(X, X.basicOpen h) : Type u) (Submonoid.powers h)
          (Algebra.linearMap (Γ(X, A.1) : Type u) (Γ(X, X.basicOpen h) : Type u)) := by
  letI := resAlg A h
  haveI := isLocalizationAway_bo A h
  rw [idealSections_basicOpen D h, Submodule.localized'_eq_span]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def idealSections_eq_localized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断は元のイデアルの局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
