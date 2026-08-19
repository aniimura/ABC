import ABC3.Found.Arakelov.PicWitness
import Mathlib.RingTheory.PicardGroup

/-!
# Arakelov (B3) 第 303 ブロック —— **★★★★★★★`PicSpecData` の非空虚 witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★★★★★★(B3) は (B1) の系である

Interface 自身がそう書いている:

> mathlib には `ClassGroup.equivPic : ClassGroup R ≃* CommRing.Pic R` があるので、
> 本条件は `equivClassGroup` を**完全に決めてしまう**
> ——すなわち **(B3) は独立の難所ではなく、(B1) の系である**。

★実際そのとおりで、

    equivClassGroup F := equivPicRing (𝓞 F) ≫ (ClassGroup.equivPic (𝓞 F))⁻¹

と取れば `equivClassGroup_compat` は **`apply_symm_apply`** で出る。

## ★★これで Arakelov 理論の 9 件がすべて閉じる

| 義務 | 状態 |
|---|---|
| B1 `PicardData` | ★閉 |
| B2 `CartierPicData` | ★閉 |
| **B3 `PicSpecData`** | ★★★**本ブロックで閉** |
| C1 `ArcSpaceData` | ★閉 |
| C2 `ProjectiveModelData` | ★閉 |
| C3 `HermitianMetricData` | ★閉 |
| D1 `APicData` | ★閉 |
| D2 `APicSpecData` | ★閉 |
| D3 `ArakelovHeightData` | ★閉 |
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

/-- ★**(B1) の `equivPicRing` と mathlib の `ClassGroup.equivPic` を繋いだ同型**。 -/
noncomputable def picSpecEquiv (F : Type) [Field F] [NumberField F] :
    picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F))) ≃ ClassGroup (𝓞 F) :=
  letI := picardDataWitness.group (Spec (CommRingCat.of (𝓞 F)))
  (picardDataWitness.equivPicRing (CommRingCat.of (𝓞 F))).toEquiv.trans
    (ClassGroup.equivPic (𝓞 F)).toEquiv.symm

theorem picSpecEquiv_compat (F : Type) [Field F] [NumberField F]
    (L : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F)))) :
    ClassGroup.equivPic (𝓞 F) (picSpecEquiv F L)
      = picardDataWitness.equivPicRing (CommRingCat.of (𝓞 F)) L :=
  (ClassGroup.equivPic (𝓞 F)).toEquiv.apply_symm_apply _

/-- ★★★★★★★**`PicSpecData` の実装**。 -/
noncomputable def picSpecDataWitness : PicSpecData where
  toPicardData := picardDataWitness
  equivClassGroup := fun F _ _ => picSpecEquiv F
  equivClassGroup_compat := fun F _ _ L => picSpecEquiv_compat F L

/-- ★★★★★★★**`PicSpecData` は非空虚である**——B3 達成。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★これが Arakelov 理論の 9 件のうち **B3** である。 -/
theorem PicSpecData.nonvacuous : Nonempty PicSpecData :=
  ⟨picSpecDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def picSpecDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(層 B——Spec 𝓞_F 上の Pic)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Interface.Arakelov
