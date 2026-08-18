import ABC3.Found.Arakelov.PicFiniteCover
import ABC3.Found.Arakelov.PicResSquare

/-!
# Arakelov (B1) 第 136 ブロック —— **構造層の制限は局所化である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★mathlib に在庫があった

`IsAffineOpen.isLocalization_of_eq_basicOpen` は

    U アフィン開、V = X.basicOpen f  ⟹  Γ(X,V) は Γ(X,U) の f での局所化

を与える。★★これを `U := D(g)`、`f := algebraMap R Γ(D g) t`、`V := D(t·g)` に当てる。

## ★★2 つの橋

| 橋 | 中身 |
|---|---|
| `specD_eq_basicOpen` | ★`D(t·g) = (Spec R).basicOpen (algebraMap t)` |
| `isLocalizedModule_specRes` | ★★`IsLocalization` を **`R` 上の `IsLocalizedModule`** へ |

★後者は `isLocalizedModule_iff_isLocalization`(`Algebra.algebraMapSubmonoid`)を通す
——`(powers t).map (algebraMap R A) = powers (algebraMap R A t)` が要点である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `specD` / `specDle` | ★`(Spec R).Opens` として固定した基本開集合 |
| `specD_eq_basicOpen` | ★`D(t·g)` は制限した `t` の basicOpen |
| `specResAlgHom` | ★★制限は `R` 代数射 |
| `isLocalizationAway_specRes` | ★★★`Away` 局所化 |
| `isLocalizedModule_specRes` | ★★★★**`powers t` の局所化加群** |

## ★★★型の 2 経路(記録)

`PrimeSpectrum.basicOpen g : Opens (PrimeSpectrum R)` と
`(Spec R).Opens` は defeq だが**推論されない**。
★`abbrev specD` で `(Spec R).Opens` 側に固定してしまうのが一番安い。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (g t : (R : Type u))


/-- ★D(r) を `(Spec R).Opens` として固定する。 -/
abbrev specD (r : (R : Type u)) : (Spec R).Opens := PrimeSpectrum.basicOpen r

/-- ★D(t·g) ≤ D(g)。 -/
theorem specDle : specD R (t * g) ≤ specD R g := basicOpenMul_le R g t

-- ★CT2: t の D(g) への制限の basicOpen は D(t*g)
theorem specD_eq_basicOpen : (specD R (t * g))
    = (Spec R).basicOpen (algebraMap (R : Type u) (Γ(Spec R, specD R g) : Type u) t) := by

  show specD R (t * g) = (Spec R).basicOpen
    (((Spec R).presheaf.map (homOfLE (le_top (a := specD R g))).op)
      ((Scheme.ΓSpecIso R).inv t))
  rw [Scheme.basicOpen_res, basicOpen_eq_of_affine]
  show specD R (t * g) = specD R g ⊓ specD R t
  rw [specD, PrimeSpectrum.basicOpen_mul]
  exact inf_comm _ _

/-- ★★制限は `R` 代数射である。 -/
noncomputable def specResAlgHom :
    (Γ(Spec R, specD R g) : Type u) →ₐ[(R : Type u)] (Γ(Spec R, specD R (t * g)) : Type u) where
  toRingHom := ((Spec R).presheaf.map (homOfLE (specDle R g t)).op).hom
  commutes' := fun r => by
    show ((Spec R).presheaf.map (homOfLE (specDle R g t)).op).hom
        (algebraMap (R : Type u) (Γ(Spec R, specD R g) : Type u) r) = _

    show (((Scheme.ΓSpecIso R).inv ≫ (Spec R).presheaf.map (homOfLE le_top).op)
        ≫ (Spec R).presheaf.map (homOfLE (specDle R g t)).op).hom r = _
    rw [Category.assoc, ← Functor.map_comp, ← op_comp]
    rfl

/-- ★★★制限は Away 局所化である。 -/
theorem isLocalizationAway_specRes :
    @IsLocalization.Away _ _ (algebraMap (R : Type u) (Γ(Spec R, specD R g) : Type u) t)
      (Γ(Spec R, specD R (t * g)) : Type u) _ (specResAlgHom R g t).toRingHom.toAlgebra :=
  (IsAffineOpen.Spec_basicOpen g).isLocalization_of_eq_basicOpen _ _ (specD_eq_basicOpen R g t)

/-- ★★★★制限は powers t の局所化加群である。 -/
theorem isLocalizedModule_specRes :
    IsLocalizedModule (Submonoid.powers t) (specResAlgHom R g t).toLinearMap := by
  letI : Algebra (Γ(Spec R, specD R g) : Type u) (Γ(Spec R, specD R (t * g)) : Type u) :=
    (specResAlgHom R g t).toRingHom.toAlgebra
  haveI : IsScalarTower (R : Type u) (Γ(Spec R, specD R g) : Type u)
      (Γ(Spec R, specD R (t * g)) : Type u) :=
    IsScalarTower.of_algebraMap_eq (fun r => ((specResAlgHom R g t).commutes r).symm)
  haveI : IsLocalization (Algebra.algebraMapSubmonoid (Γ(Spec R, specD R g) : Type u)
      (Submonoid.powers t)) (Γ(Spec R, specD R (t * g)) : Type u) := by
    show IsLocalization ((Submonoid.powers t).map
      (algebraMap (R : Type u) (Γ(Spec R, specD R g) : Type u))) _
    rw [Submonoid.map_powers]
    exact isLocalizationAway_specRes R g t
  exact isLocalizedModule_iff_isLocalization.2 inferInstance


/-! ## ★出典の紐付け(`.src`) -/

def isLocalizedModule_specRes.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——構造層の制限は局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
