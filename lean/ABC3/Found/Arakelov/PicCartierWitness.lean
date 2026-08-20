import ABC3.Found.Arakelov.PicDivisorPull
import ABC3.Found.Arakelov.PicPrincipalAffine
import ABC3.Found.Arakelov.PicWitness

/-!
# Arakelov (B2) 第 241 ブロック —— ★★★★★★★★★★**`CartierPicData` の witness**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★★★B2 が閉じた(2026-08-19)

`Interface` の `CartierPicData` は 14 欄。★**すべて埋まった。**

| 欄 | witness | 出所 |
|---|---|---|
| `toPicardData` | `picardDataWitness` | ★第 146(B1) |
| `IsCartierDivisor` | `IsCartier` | 第 184 |
| `isCartierDivisor_top` | `isCartier_top` | 第 186 |
| `isCartierDivisor_mul` | `isCartier_mul` | 第 186 |
| `isCartierDivisor_comap` | `isCartier_comap` | 第 202(平坦) |
| `isCartierDivisor_affine` | `isCartier_iff_top` | 第 184 |
| `ofDivisor` | `ofDivisorSheaf` | 第 186 |
| `ofDivisor_mul` | `ofDivisorSheaf_mul` | 第 190 |
| `ofDivisor_top` | `ofDivisorSheaf_top` | 第 187 |
| **`ofDivisor_pullback`** | **`ofDivisorSheaf_pullback`** | ★★**第 240** |
| `IsPrincipalDivisor` | `IsPrincipalDivisor` | 第 191 |
| `ofDivisor_eq_one_iff` | `ofDivisorSheaf_eq_one_iff` | 第 191 |
| `isPrincipalDivisor_mul` | `isPrincipalDivisor_mul` | 第 191 |
| `isPrincipalDivisor_affine` | `isPrincipalDivisor_affine` | 第 193 |

★★★第 176 ブロック(B2 の開始)から 66 ブロックである。

## ★★Interface を 3 回直した——**満たそうとして初めて見つかった誤り**

| 日付 | 直した欄 | 反例 |
|---|---|---|
| 2026-08-18 | `ofDivisor` を Cartier に限定 | ★`ℚ[x,y]` の `(x,y)` |
| 2026-08-19 | `isCartierDivisor_comap` に平坦性 | ★`Spec k ⟶ Spec k[x]` の原点 |
| 2026-08-19 | `isPrincipalDivisor_affine` に Cartier | ★`k[x]/(x²)` の `(x)` |

★★★どれも「証明を書こうとしたら通らなかった」ことで見つかった
——**Interface は消費するだけでは検証できない**。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-- ★★★★★★★★★★**B2 の witness**——`CartierPicData` の 14 欄すべて。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 176 から第 241 までの 66 ブロックの到達点である。 -/
noncomputable def cartierPicDataWitness : ABC3.Interface.Arakelov.CartierPicData where
  toPicardData := picardDataWitness
  IsCartierDivisor := fun X D => IsCartier X D
  isCartierDivisor_top := fun X => isCartier_top X
  isCartierDivisor_mul := fun _ D E hD hE => isCartier_mul D E hD hE
  isCartierDivisor_comap := fun f hf D hD => @isCartier_comap _ _ f hf D hD
  isCartierDivisor_affine := fun _ D => isCartier_iff_top D
  ofDivisor := fun _ D => ofDivisorSheaf D
  ofDivisor_mul := fun _ D E hD hE => ofDivisorSheaf_mul D E hD hE
  ofDivisor_top := fun X => ofDivisorSheaf_top X
  ofDivisor_pullback := fun f hf D hD => @ofDivisorSheaf_pullback _ _ f hf D hD
  IsPrincipalDivisor := fun X D => IsPrincipalDivisor X D
  ofDivisor_eq_one_iff := fun _ D hD => ofDivisorSheaf_eq_one_iff D hD
  isPrincipalDivisor_mul := fun _ D E hD hE => isPrincipalDivisor_mul D E hD hE
  isPrincipalDivisor_affine := fun R D hD => isPrincipalDivisor_affine R D hD

/-- ★★★★★★★★★★**B2 は非空虚である**。 -/
theorem cartierPicData_nonvacuous : Nonempty ABC3.Interface.Arakelov.CartierPicData :=
  ⟨cartierPicDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def cartierPicDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——CartierPicData の witness)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
