import ABC3.Found.Arakelov.PicIdealBase
import Mathlib.AlgebraicGeometry.Morphisms.Flat
import Mathlib.AlgebraicGeometry.IdealSheaf.Basic

/-!
# Arakelov (B2) 第 184 ブロック —— **Cartier 因子はアフィン開全体で可逆**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`isCartierDivisor_affine` の ⟸ 向きが出た

    アフィン `X` で `D.ideal ⊤` が可逆 ⟹ すべてのアフィン開 `A` で `D.ideal A` が可逆

★★機構は 3 つの合流:

| 部品 | 出どころ |
|---|---|
| `D.ideal A = (D.ideal B).map (制限)` | mathlib `IdealSheafData.map_ideal` |
| **アフィン開の制限は平坦** | mathlib `Scheme.Hom.flat_appLE`(`𝟙 X` は開埋め込み) |
| 平坦底変換で可逆イデアルは可逆 | ★第 183 |

## ★★これで `IsCartier` の定義が「アフィン開すべて」で書ける

    IsCartier X D := ∀ A : X.affineOpens, Module.Invertible Γ(X,A) (D.ideal A)

★この形なら `isLocallyTrivial_idealSheaf`(第 162)がそのまま当たり、
かつアフィンでは `⊤` 1 点に落ちる(`isCartier_iff_top`)。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `flat_restrict` | ★★アフィン開の制限は平坦 |
| `invertible_ideal_of_le` | ★★★大きい方で可逆なら小さい方でも可逆 |
| `IsCartier` | ★★★**Cartier 因子** |
| `isCartier_iff_top` | ★★★★**アフィンでは `⊤` だけ見ればよい** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★**アフィン開の制限は平坦である**——`𝟙 X` は開埋め込みだから。 -/
theorem flat_restrict {X : Scheme.{u}} {A B : X.Opens} (hA : IsAffineOpen A) (hB : IsAffineOpen B)
    (h : A ≤ B) :
    letI : Algebra (Γ(X, B) : Type u) (Γ(X, A) : Type u) :=
      (X.presheaf.map (homOfLE h).op).hom.toAlgebra
    Module.Flat (Γ(X, B) : Type u) (Γ(X, A) : Type u) := by
  have e : A ≤ (𝟙 X) ⁻¹ᵁ B := by simpa using h
  have hf := Scheme.Hom.flat_appLE (𝟙 X) hB hA e
  have happ : Scheme.Hom.appLE (𝟙 X) B A e = X.presheaf.map (homOfLE h).op := by
    show 𝟙 _ ≫ X.presheaf.map (homOfLE e).op = _
    exact Category.id_comp _
  rw [happ] at hf
  exact hf

/-- ★★★**大きいアフィン開で可逆なら小さいアフィン開でも可逆である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`D.ideal A` は `D.ideal B` の**拡大**(`map_ideal`)であり、制限は**平坦**だから
第 183 が当たる。 -/
theorem invertible_ideal_of_le {X : Scheme.{u}} (D : X.IdealSheafData) {A B : X.affineOpens}
    (h : A ≤ B) (hB : Module.Invertible (Γ(X, B.1) : Type u) (D.ideal B)) :
    Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A) := by
  letI : Algebra (Γ(X, B.1) : Type u) (Γ(X, A.1) : Type u) :=
    (X.presheaf.map (homOfLE (show A.1 ≤ B.1 from h)).op).hom.toAlgebra
  haveI : Module.Flat (Γ(X, B.1) : Type u) (Γ(X, A.1) : Type u) :=
    flat_restrict A.2 B.2 (show A.1 ≤ B.1 from h)
  haveI := hB
  have key := invertible_map_of_flat (R := (Γ(X, B.1) : Type u)) (S := (Γ(X, A.1) : Type u))
    (I := D.ideal B)
  rw [RingHom.algebraMap_toAlgebra] at key
  exact cast (congrArg (fun (J : Ideal (Γ(X, A.1) : Type u)) =>
    Module.Invertible (Γ(X, A.1) : Type u) J) (D.map_ideal h)) key

/-- ★★★**Cartier 因子**——すべてのアフィン開で可逆イデアル。 -/
def IsCartier (X : Scheme.{u}) (D : X.IdealSheafData) : Prop :=
  ∀ A : X.affineOpens, Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)

/-- ★★★★**アフィンスキームでは `⊤` だけ見ればよい**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `CartierPicData.isCartierDivisor_affine` の中身である。 -/
theorem isCartier_iff_top {X : Scheme.{u}} [IsAffine X] (D : X.IdealSheafData) :
    IsCartier X D ↔
      Module.Invertible (Γ(X, ⊤) : Type u) (D.ideal ⟨⊤, isAffineOpen_top X⟩) := by
  refine ⟨fun h => h ⟨⊤, isAffineOpen_top X⟩, fun h A => ?_⟩
  exact invertible_ideal_of_le D (B := ⟨⊤, isAffineOpen_top X⟩)
    (show A.1 ≤ (⊤ : X.Opens) from le_top) h

/-! ## ★出典の紐付け(`.src`) -/

def isCartier_iff_top.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Cartier 因子はアフィンでは ⊤ だけで決まる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
