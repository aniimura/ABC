import ABC3.Found.Arakelov.PicWitness

/-!
# Arakelov (B3) 第 147 ブロック —— ★★★★★★**`Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> (i) We shall write ADiv(F) for the group of arithmetic divisors on F

## ★★★★★★B1 が入ったら **1 ブロック**で出た

`PicSpecData` は 3 欄しかない:

| 欄 | witness |
|---|---|
| `toPicardData` | ★第 146(B1) |
| `equivClassGroup` | ★第 145 と mathlib の `ClassGroup.equivPic` の合成 |
| `equivClassGroup_compat` | ★`apply_symm_apply` 一発 |

## ★★在庫の実測(2026-08-18)

    ClassGroup.equivPic (R) [CommRing R] [IsDomain R] : ClassGroup R ≃* CommRing.Pic R

★mathlib の `RingTheory/PicardGroup.lean` 849 行。**有った。**

★★したがって

    Pic(Spec 𝓞_F) ≃* CommRing.Pic 𝓞_F ≃* ClassGroup 𝓞_F
        ↑ 第 145                ↑ mathlib

★★★Interface 自身が「(B3) は独立の難所ではなく (B1) の系である」と
書いていた通りであった——**測ってから作れば見積もりは効く**の 2 例目。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory NumberField

/-- ★★★★★★**B3 の witness**——`Pic(Spec 𝓞_F) ≅ ClassGroup 𝓞_F`。

原文 (GenEll p.4):
> (i) We shall write ADiv(F) for the group of arithmetic divisors on F

★★★B1(第 146)と mathlib の `ClassGroup.equivPic` の合成である。 -/
noncomputable def picSpecDataWitness : ABC3.Interface.Arakelov.PicSpecData where
  toPicardData := picardDataWitness
  equivClassGroup := fun F =>
    (equivPicRingSheaf (CommRingCat.of (𝓞 F))).toEquiv.trans
      (ClassGroup.equivPic (𝓞 F)).toEquiv.symm
  equivClassGroup_compat := by
    intro F _ _ L
    exact (ClassGroup.equivPic (𝓞 F)).apply_symm_apply _

/-! ## ★出典の紐付け(`.src`) -/

def picSpecDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.2, (i)(層 B——Spec 𝓞_F 上の Pic)",
    sectionId := "genell-def-1-2-i" }

end ABC3.Found.Arakelov
