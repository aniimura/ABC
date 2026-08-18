import ABC3.Found.Arakelov.PicPullSurj
import ABC3.Found.Arakelov.PicComapImage

/-!
# Arakelov (B2) 第 219 ブロック —— ★★★**同型で挟んだ等式を外す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★第 211 と第 218 を繋ぐ

第 211(`comap_ideal_image`)は**転送つき**の等式:

    I.comap eA⁻¹ = (J.comap eB⁻¹).map ψ

第 218 が要求するのは**転送なし**の等式:

    I = J.map (eB ≫ ψ ≫ eA⁻¹)

★本ブロックがその変換である。

## ★★機構は `map_comap` 2 回

| 段 | 補題 |
|---|---|
| `I = (I.comap eA⁻¹).map eA⁻¹` | ★`Ideal.map_comap_of_surjective`(`eA⁻¹` は全射) |
| `J.comap eB⁻¹ = J.map eB` | ★第 201 の `comap_inv_eq_map` |
| 合成の collapse | ★`Ideal.map_map` |

★★★§9-232 で「転送を外した形は `whnf` で落ちる」と測ったが、
**環の同型を抽象化した形**(スキームの綴りを持ち込まない)なら**通る**。
★これが逃げ道であった——[[ring-instance-two-paths]] の
「抽象の側で書いて具体に当てる」形である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ideal_of_comap_eq` | ★★★**同型で挟んだ等式を外す** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★★**同型で挟んだイデアルの等式を外す**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが第 211(転送つき)と第 218(転送なし)を繋ぐ。 -/
theorem ideal_of_comap_eq {R S T U : CommRingCat.{u}} (eA : R ≅ S) (eB : T ≅ U)
    (I : Ideal (R : Type u)) (J : Ideal (T : Type u)) (ψ : U ⟶ S)
    (h : I.comap eA.inv.hom = (J.comap eB.inv.hom).map ψ.hom) :
    I = J.map (eB.hom ≫ ψ ≫ eA.inv).hom := by
  have hsurj : Function.Surjective eA.inv.hom :=
    eA.commRingCatIsoToRingEquiv.symm.surjective
  have h1 : I = (I.comap eA.inv.hom).map eA.inv.hom :=
    (Ideal.map_comap_of_surjective _ hsurj I).symm
  rw [h1, h, Ideal.map_map, comap_inv_eq_map eB J, Ideal.map_map]
  rfl


/-! ## ★出典の紐付け(`.src`) -/

def ideal_of_comap_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——同型で挟んだ等式を外す)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
