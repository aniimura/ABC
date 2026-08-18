import ABC3.Found.Arakelov.PicAffCover

/-!
# Arakelov (B2) 第 156 ブロック —— **切断と局所化加群の同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 124–127 を `idealSections` で敷き直す

`tilde M` のときは第 118–127 で `Γ(tilde M, D f) ≅ M_f` を作った。
★イデアル層では **mathlib の `Submodule.toLocalized'` が同じ役割**を果たすので、
第 153 と合わせて**一気に**同型が出る。

| 宣言 | 内容 | `tilde` 版の対応 |
|---|---|---|
| `idealAwayEquiv` | ★`I_f ≅ idealSections D (D f)` | 第 118 `tildeAwayEquiv` |
| `awayRingEquivX` | ★`Γ(X, D f) ≃ₐ Localization (powers f)` | 第 124 `awayRingEquiv` |
| `modOnLocalizedX` | ★`I_f` への `Γ(X, D f)` 作用 | 第 124 `modOnLocalized` |

★★`awayRingEquivX` は `IsLocalization.algEquiv` 一発——
アフィン開の制限が `Away f` の局所化であること(mathlib)が効く。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

/-- ★★局所化加群と切断の同型。 -/
noncomputable def idealAwayEquiv (f : (Γ(X, A.1) : Type u)) :
    LocalizedModule (Submonoid.powers f) (D.ideal A)
      ≃ₗ[(Γ(X, A.1) : Type u)] (idealSections D (X.basicOpen f)) :=
  letI := resAlg A f
  haveI := isLocalizationAway_bo A f
  (IsLocalizedModule.iso (Submonoid.powers f)
      ((D.ideal A).toLocalized' (Γ(X, X.basicOpen f) : Type u) (Submonoid.powers f)
        (Algebra.linearMap (Γ(X, A.1) : Type u) (Γ(X, X.basicOpen f) : Type u))))
    ≪≫ₗ (LinearEquiv.ofEq _ _ (idealSections_eq_localized D A f).symm).restrictScalars
        (Γ(X, A.1) : Type u)

/-- ★`mk m 1` での値——`m` の制限そのもの。 -/
theorem idealAwayEquiv_mk_one (f : (Γ(X, A.1) : Type u)) (m : (D.ideal A)) :
    ((idealAwayEquiv D A f (LocalizedModule.mk m 1)) : (Γ(X, X.basicOpen f) : Type u))
      = (X.presheaf.map (homOfLE (X.basicOpen_le f)).op).hom (m : (Γ(X, A.1) : Type u)) := by
  letI := resAlg A f
  haveI := isLocalizationAway_bo A f
  show (((LinearEquiv.ofEq _ _ (idealSections_eq_localized D A f).symm).restrictScalars
      (Γ(X, A.1) : Type u))
      (IsLocalizedModule.iso (Submonoid.powers f)
        ((D.ideal A).toLocalized' (Γ(X, X.basicOpen f) : Type u) (Submonoid.powers f)
          (Algebra.linearMap (Γ(X, A.1) : Type u) (Γ(X, X.basicOpen f) : Type u)))
        (LocalizedModule.mk m 1)) : (Γ(X, X.basicOpen f) : Type u)) = _
  rw [IsLocalizedModule.iso_mk_one]
  rfl

/-- ★`Γ(X, D(f)) ≃ₐ Localization (powers f)`。 -/
noncomputable def awayRingEquivX (f : (Γ(X, A.1) : Type u)) :
    letI := resAlg A f
    (Γ(X, X.basicOpen f) : Type u) ≃ₐ[(Γ(X, A.1) : Type u)] Localization (Submonoid.powers f) :=
  letI := resAlg A f
  haveI := isLocalizationAway_bo A f
  IsLocalization.algEquiv (Submonoid.powers f) _ _

/-- ★★局所化加群への `Γ(X, D(f))` 作用。 -/
@[reducible] noncomputable def modOnLocalizedX (f : (Γ(X, A.1) : Type u)) :
    Module (Γ(X, X.basicOpen f) : Type u)
      (LocalizedModule (Submonoid.powers f) (D.ideal A)) :=
  letI := resAlg A f
  Module.compHom _ (awayRingEquivX A f).toRingHom


/-! ## ★出典の紐付け(`.src`) -/

def idealAwayEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断と局所化加群の同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
