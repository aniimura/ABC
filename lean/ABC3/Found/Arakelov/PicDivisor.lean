import ABC3.Found.Arakelov.PicIdealMul
import ABC3.Found.Arakelov.PicEvalIso
import ABC3.Found.Arakelov.PicIdealLT

/-!
# Arakelov (B2) 第 186 ブロック —— ★★★★★★**因子から直線束へ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`ofDivisor` が書けた

    D  ↦  idealSheaf D  ↦  InvSheaf X  ↦  PicSheaf X

3 本の合流である:

| 部品 | ブロック |
|---|---|
| `idealSheaf D`(イデアル層は層加群) | 第 151 |
| `IsCartier ⟹ IsLocallyTrivial` | 第 162 |
| **局所自明 ⟹ 可逆層** | ★第 182 |

★Cartier でないときは `1` を返す(`dite`)——`ofDivisor` は全域でなければならず、
かつ Interface の法則はすべて Cartier 条件つきだから、これで整合する。

## ★★Cartier の閉性は環レベルの補題がそのまま上がる

`(D * E).ideal A = D.ideal A * E.ideal A` は mathlib で **`rfl`**、
`(⊤).ideal A = ⊤` も **`rfl`** なので、第 185 の
`invertible_mul` / `invertible_top` が**そのまま**当たる。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `isCartier_top` | ★空因子は Cartier |
| `isCartier_mul` | ★★Cartier は積で閉じる |
| `divisorInvSheaf` | ★★★★因子から可逆層 |
| `ofDivisorSheaf` | ★★★★★★**因子から `Pic`** |
| `divisorInvSheaf_carrier` | ★Cartier なら台はイデアル層 |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★**空因子は Cartier である**。 -/
theorem isCartier_top (X : Scheme.{u}) : IsCartier X ⊤ := fun _ => invertible_top

/-- ★★**Cartier は積で閉じる**。 -/
theorem isCartier_mul (D E : X.IdealSheafData)
    (hD : IsCartier X D) (hE : IsCartier X E) : IsCartier X (D * E) := fun A => by
  haveI := hD A
  haveI := hE A
  exact invertible_mul

open Classical in
/-- ★★★★**因子から可逆層**——Cartier でなければ単位元を返す。 -/
noncomputable def divisorInvSheaf (D : X.IdealSheafData) : InvSheaf X :=
  if h : IsCartier X D then
    InvSheaf.ofLocallyTrivial (idealSheaf D) (isLocallyTrivial_idealSheaf D h)
  else InvSheaf.one X

/-- ★★★★★★**因子から `Pic`**——`D ↦ 𝒪_X(D)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `CartierPicData.ofDivisor` である。 -/
noncomputable def ofDivisorSheaf (D : X.IdealSheafData) : PicSheaf X :=
  PicSheaf.mk (divisorInvSheaf D)

/-- ★**Cartier なら台はイデアル層そのものである**。 -/
theorem divisorInvSheaf_carrier (D : X.IdealSheafData) (h : IsCartier X D) :
    (divisorInvSheaf D).carrier = idealSheaf D := by
  classical
  rw [divisorInvSheaf, dif_pos h]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def ofDivisorSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——因子から直線束)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
