import ABC3.Found.Arakelov.PicIdealMulHom
import ABC3.Found.Arakelov.PicAffineSieve
import ABC3.Found.Arakelov.PicEvalIso

/-!
# Arakelov (B2) 第 190 ブロック —— ★★★★★★**因子の積は直線束のテンソル積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`ofDivisor_mul` が出た

    𝒪_X(D + E) = 𝒪_X(D) ⊗ 𝒪_X(E)

すなわち `ofDivisor` は**準同型**である。これで `CartierPicData` は 14 欄中 9 欄。

## ★★筋は第 182 と同じ 3 段

| 段 | 内容 | ブロック |
|---|---|---|
| 1 | アフィン開の上で全単射 | ★第 188 |
| 2 | ゆえに**局所全単射** | ★第 189 |
| 3 | 局所全単射は層化で同型 | mathlib(第 182 で実測) |

★★★第 189 を**再利用できる形**(「アフィン開で全単射なら局所全単射」)で
書いておいたので、本ブロックは 3 行で済んだ。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `isIso_sheafify_mulHom` | ★★★★積の射の層化は同型 |
| `mulModulesIso` | ★★★★★`𝒪(D) ⊗ 𝒪(E) ≅ 𝒪(D+E)` |
| `ofDivisorSheaf_mul` | ★★★★★★**`ofDivisor` は準同型** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
variable {X : Scheme.{u}} (D E : X.IdealSheafData)


/-- ★★★★**積の射の層化は同型である**。 -/
theorem isIso_sheafify_mulHom (hD : IsCartier X D) :
    IsIso ((sheafifyFunctor X).map (mulHom D E)) := by
  obtain ⟨hi, hs⟩ := locallyBijective_of_bijective_on_affine
    ((PresheafOfModules.toPresheaf _).map (mulHom D E))
    (fun A => bijective_mulHom_app D E A hD)
  haveI hi' : Presheaf.IsLocallyInjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (mulHom D E)) := hi
  haveI hs' : Presheaf.IsLocallySurjective (Opens.grothendieckTopology X)
      ((PresheafOfModules.toPresheaf X.ringCatSheaf.obj).map (mulHom D E)) := hs
  have hW : (MorphismProperty.inverseImage (Opens.grothendieckTopology X).W
      (PresheafOfModules.toPresheaf X.ringCatSheaf.obj)) (mulHom D E) :=
    GrothendieckTopology.W_of_isLocallyBijective _ _
  have heq := PresheafOfModules.inverseImage_W_toPresheaf_eq_inverseImage_isomorphisms
    (R := X.ringCatSheaf) (𝟙 X.ringCatSheaf.obj)
  have h2 : (MorphismProperty.isomorphisms (SheafOfModules X.ringCatSheaf)).inverseImage
      (PresheafOfModules.sheafification (𝟙 X.ringCatSheaf.obj)) (mulHom D E) := heq ▸ hW
  exact h2

/-- ★★★★★★因子の積は層のテンソル積。 -/
noncomputable def mulModulesIso (hD : IsCartier X D) :
    tensorModules (idealSheaf D) (idealSheaf E) ≅ idealSheaf (D * E) :=
  @asIso _ _ _ _ ((sheafifyFunctor X).map (mulHom D E)) (isIso_sheafify_mulHom D E hD)
    ≪≫ sheafifyValIso (idealSheaf (D * E))

theorem ofDivisorSheaf_mul (hD : IsCartier X D) (hE : IsCartier X E) :
    ofDivisorSheaf (D * E) = ofDivisorSheaf D * ofDivisorSheaf E := by
  classical
  refine PicSheaf.mk_eq_mk ?_
  rw [divisorInvSheaf_carrier (D * E) (isCartier_mul D E hD hE)]
  show idealSheaf (D * E) ≅
    tensorModules (divisorInvSheaf D).carrier (divisorInvSheaf E).carrier
  rw [divisorInvSheaf_carrier D hD, divisorInvSheaf_carrier E hE]
  exact (mulModulesIso D E hD).symm


/-! ## ★出典の紐付け(`.src`) -/

def ofDivisorSheaf_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——因子の積は直線束のテンソル積)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
