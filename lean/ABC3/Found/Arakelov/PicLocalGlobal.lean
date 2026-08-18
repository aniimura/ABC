import ABC3.Found.Arakelov.PicIdealLTPt
import ABC3.Found.Arakelov.PicEvalIso
import ABC3.Found.Arakelov.PicGammaInv
import ABC3.Found.Arakelov.PicPrincipalAffine
import ABC3.Found.Arakelov.PicBaseChange

/-!
# Arakelov (B2) 第 198 ブロック —— ★★★★★★**可逆性は Zariski 局所的**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★mathlib の TODO を 3 ブロックで埋めた

`Mathlib/RingTheory/PicardGroup.lean` の TODO:

> Establish other characterizations of invertible modules, e.g. they are modules that
> become free of rank one when localized at every prime ideal.

その一部——**局所的に可逆なら可逆**——が本ブロックである:

    ∀ x ∈ Spec R, ∃ アフィン開 A ∋ x, D.ideal A が可逆
      ⟹  D.ideal ⊤ が可逆

## ★★★環の側から作ると長い、層の側から作ると短い(§9-215 の実測)

| 道 | 欠落 | 見積 |
|---|---|---|
| 環(`contractLeft` の局所性) | `(Mᵛ)_S ≅ (M_S)ᵛ`(有限表示) | 8–15 |
| ★層(第 196 → 182 → 132 系) | 係数環の橋のみ(★第 197) | **3** |

★★★**層の側を先に作っておいたことが効いた。**B1 で積んだ 146 ブロックが、
ここで「環の定理を 3 段で出す」形で返ってきている。

## ★★機構

    点ごとに可逆
      → IsLocallyTrivial            (第 196)
      → InvSheaf (Spec R)           (第 182)
      → Γ が R 上可逆               (第 132 系 invertible_gammaCarrier)
      → Γ が Γ(Spec R,⊤) 上可逆     (★第 197、algebraMap が全射)

★係数環の綴りの橋(`R` と `Γ(Spec R, ⊤)`)は instance を 1 つ置いて渡す。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `invertible_gamma_of_pointwise` | ★★大域切断は `R` 上可逆 |
| `moduleGammaCarrier` | ★綴りの橋(instance) |
| `invertible_ideal_top_of_pointwise` | ★★★★★★**可逆性は Zariski 局所的** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (D : (Spec R).IdealSheafData)

/-- ★★★局所的に可逆なら大域切断は `R` 上可逆。 -/
theorem invertible_gamma_of_pointwise
    (hcart : ∀ x : Spec R, ∃ A : (Spec R).affineOpens, ∃ _ : x ∈ A.1,
      Module.Invertible (Γ(Spec R, A.1) : Type u) (D.ideal A)) :
    Module.Invertible (R : Type u)
      (AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D) : Type u) :=
  invertible_gammaCarrier R
    (InvSheaf.ofLocallyTrivial (idealSheaf D)
      (isLocallyTrivial_idealSheaf_of_pointwise D hcart))




/-- ★大域切断加群の `Γ(⊤)` 加群構造(綴りの橋)。 -/
noncomputable instance moduleGammaCarrier : Module (Γ(Spec R, ⊤) : Type u)
    ((AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D)) : Type u) :=
  inferInstanceAs (Module (Γ(Spec R, ⊤) : Type u) ((idealSections D ⊤) : Type u))

/-- ★★★★★★局所的に可逆なら ⊤ でも可逆(可逆性の Zariski 局所性)。 -/
theorem invertible_ideal_top_of_pointwise
    (hcart : ∀ x : Spec R, ∃ A : (Spec R).affineOpens, ∃ _ : x ∈ A.1,
      Module.Invertible (Γ(Spec R, A.1) : Type u) (D.ideal A)) :
    Module.Invertible (Γ(Spec R, ⊤) : Type u)
      (D.ideal ⟨⊤, isAffineOpen_top (Spec R)⟩) := by
  have hs : Function.Surjective (algebraMap (R : Type u) (Γ(Spec R, ⊤) : Type u)) := by
    show Function.Surjective (Scheme.ΓSpecIso R).inv.hom
    exact (Scheme.ΓSpecIso R).commRingCatIsoToRingEquiv.symm.surjective
  haveI : IsScalarTower (R : Type u) (Γ(Spec R, ⊤) : Type u)
      ((AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D)) : Type u) :=
    ⟨fun r s x => by rw [Algebra.smul_def, mul_smul]; rfl⟩
  haveI := invertible_gamma_of_pointwise R D hcart
  have key := invertible_of_surjective_algebraMap
    (M := ((AlgebraicGeometry.moduleSpecΓFunctor.obj (idealSheaf D)) : Type u)) hs
  exact cast (congrArg (fun (J : Ideal (Γ(Spec R, ⊤) : Type u)) =>
    Module.Invertible (Γ(Spec R, ⊤) : Type u) J)
    (idealSections_affine D ⟨⊤, isAffineOpen_top (Spec R)⟩)) key

/-! ## ★出典の紐付け(`.src`) -/

def invertible_ideal_top_of_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆性は Zariski 局所的)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
