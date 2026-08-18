import ABC3.Found.Arakelov.PicStalkBasis
import Mathlib.Topology.Sheaves.Sheafify

/-!
# Arakelov (B1) 第 91 ブロック —— ★★★★★★**`tilde` はテンソルを保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★組み立て

    tensorModules (tilde M) (tilde N)  ≅  tilde (M ⊗_R N)

★★★**`M`・`N` の可逆性は要らない**——任意の加群で成り立つ。

## ★★機構 —— 5 段の合成

| 段 | 出典 |
|---|---|
| 比較射は基本開集合で全単射 | ★第 89 ブロック |
| 基本開集合は基底 | ★mathlib `PrimeSpectrum.isBasis_basic_opens` |
| 基底で全単射 ⟹ 茎で同型 | ★第 90 ブロック |
| 層化の unit は茎で同型 | ★mathlib `stalkFunctor_map_unit_toSheafify_isIso` |
| 茎で同型 ⟹ 同型 | ★第 77 ブロック |

★★★`desc` と `pre` は層化の unit で繋がっている
(mathlib `toPresheaf_map_sheafificationHomEquiv_def`)ので、
2 つの茎同型から 3 つ目が出る。

## ★★本ブロックで取れるもの

| 定理・定義 | 内容 |
|---|---|
| `isIso_tildeTensorDesc` | ★★★★★★**比較射は同型** |
| `tildeTensorIso` | ★★★★★★**`tilde M ⊗ tilde N ≅ tilde (M ⊗ N)`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace TensorProduct
open StructureSheaf

variable (R : CommRingCat.{u}) (M N : ModuleCat.{u} (R : Type u))

/-- ★`desc` と `pre` は層化の unit で繋がっている。 -/
theorem toPresheaf_tildeTensorPre :
    (PresheafOfModules.toPresheaf (Spec R).ringCatSheaf.obj).map (tildeTensorPre R M N)
      = CategoryTheory.toSheafify (Opens.grothendieckTopology (Spec R))
          ((tilde M).val ⊗ (tilde N).val).presheaf
        ≫ (PresheafOfModules.toPresheaf (Spec R).ringCatSheaf.obj).map
          (tildeTensorDesc R M N).val := by
  have h : (PresheafOfModules.sheafificationHomEquiv (𝟙 (Spec R).ringCatSheaf.obj))
      (tildeTensorDesc R M N) = tildeTensorPre R M N := Equiv.apply_symm_apply _ _
  rw [← h]
  exact PresheafOfModules.toPresheaf_map_sheafificationHomEquiv_def _ _

/-- ★★前層射は茎で同型である(基本開集合が基底であることから)。 -/
theorem isIso_stalk_tildeTensorPre (x : (Spec R)) :
    IsIso ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      ((PresheafOfModules.toPresheaf (Spec R).ringCatSheaf.obj).map (tildeTensorPre R M N))) := by
  refine isIso_stalkFunctor_map_of_isBasis
    (PrimeSpectrum.isBasis_basic_opens (R := (R : Type u))) ?_ x
  rintro U ⟨f, rfl⟩
  exact tensorSectionMap_bijective (R : Type u) M N f

/-- ★★★★★★**比較射は同型である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが「**`tilde` はテンソルを保つ**」である。 -/
theorem isIso_tildeTensorDesc : IsIso (tildeTensorDesc R M N) := by
  refine isIso_of_stalkIso _ (fun x => ?_)
  haveI hsh : IsIso ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      (CategoryTheory.toSheafify (Opens.grothendieckTopology (Spec R))
        ((tilde M).val ⊗ (tilde N).val).presheaf)) :=
    TopCat.Presheaf.stalkFunctor_map_unit_toSheafify_isIso x AddCommGrpCat.{u} _
  haveI hc : IsIso ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      (CategoryTheory.toSheafify (Opens.grothendieckTopology (Spec R))
        ((tilde M).val ⊗ (tilde N).val).presheaf)
      ≫ (TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
        ((SheafOfModules.toSheaf (Spec R).ringCatSheaf).map (tildeTensorDesc R M N)).hom) := by
    rw [← Functor.map_comp]
    have h := isIso_stalk_tildeTensorPre R M N x
    rwa [toPresheaf_tildeTensorPre] at h
  exact @IsIso.of_isIso_comp_left _ _ _ _ _
    ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      (CategoryTheory.toSheafify (Opens.grothendieckTopology (Spec R))
        ((tilde M).val ⊗ (tilde N).val).presheaf))
    ((TopCat.Presheaf.stalkFunctor AddCommGrpCat.{u} x).map
      ((SheafOfModules.toSheaf (Spec R).ringCatSheaf).map (tildeTensorDesc R M N)).hom)
    hsh hc

/-- ★★★★★★**`tilde M ⊗ tilde N ≅ tilde (M ⊗_R N)`**。 -/
noncomputable def tildeTensorIso :
    tensorModules (tilde M) (tilde N)
      ≅ tilde (ModuleCat.of (R : Type u) (M ⊗[(R : Type u)] N)) :=
  haveI := isIso_tildeTensorDesc R M N
  asIso (tildeTensorDesc R M N)

/-! ## ★出典の紐付け(`.src`) -/

def tildeTensorIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde はテンソルを保つ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
