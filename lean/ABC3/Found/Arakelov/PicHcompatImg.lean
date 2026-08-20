import ABC3.Found.Arakelov.PicSquareGen

/-!
# Arakelov (B2) 第 235 ブロック —— ★★★★★★★**`hcompat` が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★B2 の最後の 1 点

第 221 ブロック `surjective_pair` は仮定として

    f.appLE (hB.fromSpec ''ᵁ ⊤) (hA.fromSpec ''ᵁ ⊤) hle
      = (hB.fromSpec.appIso ⊤).hom ≫ (ΓSpecIso.hom ≫ f.appLE B A i ≫ ΓSpecIso.inv)
          ≫ (hA.fromSpec.appIso ⊤).inv

を取っていた。★本ブロックでこれが**証明された**——仮定が消える。

## ★★両辺を `appLE` の合成に直すと、第 234 の正方形そのものになる

| 段 | 使うもの |
|---|---|
| `(appIso ⊤).hom = appLE (fromSpec ''ᵁ ⊤) ⊤` | ★第 232 |
| `appLE ≫ appLE = (合成).appLE` | ★mathlib `appLE_comp_appLE` |
| `ΓSpecIso.hom ≫ φ ≫ ΓSpecIso.inv = (Spec.map φ).appLE ⊤ ⊤` | ★第 228 |
| 2 つの合成が等しい | ★★★★第 234(両側自由な正方形) |

★★★左辺は `(hA.fromSpec ≫ f).appLE (hB.fromSpec ''ᵁ ⊤) ⊤`、
右辺は `(Spec.map (f.appLE B A i) ≫ hB.fromSpec).appLE (hB.fromSpec ''ᵁ ⊤) ⊤` になり、
**第 234 の正方形を `V := hB.fromSpec ''ᵁ ⊤`、`W := ⊤` で読む**だけで一致する。

| 定理 | 内容 |
|---|---|
| `hcompat_key` | ★★★★★`.hom` を掛けた形(同型を外す前) |
| `hcompat_image` | ★★★★★★★**`hcompat` 本体** |
| `surjective_affine` | ★★★★★★**比較射はアフィンの対で全射**(仮定なし) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★★**`hcompat` を `.hom` を掛けた形で**——同型を外す前の等式。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★両辺とも `appLE` 2 つの合成に直り、第 234 の正方形に落ちる。 -/
theorem hcompat_key {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (hle : (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens))
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))) :
    f.appLE (hB.fromSpec ''ᵁ ⊤) (hA.fromSpec ''ᵁ ⊤) hle ≫ (hA.fromSpec.appIso ⊤).hom
      = (hB.fromSpec.appIso ⊤).hom ≫ (Scheme.ΓSpecIso (Γ(Y, B))).hom
          ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv := by
  rw [appIso_top_as_app hA, appIso_top_as_app hB, Scheme.Hom.appLE_comp_appLE]
  have e3 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (Spec.map (f.appLE B A i)) ⁻¹ᵁ ⊤ := le_top
  have hid : (Spec Γ(X, A)).presheaf.map (homOfLE (le_top : (⊤ : (Spec Γ(X, A)).Opens) ≤ ⊤)).op
      = 𝟙 _ := by
    rw [show (homOfLE (le_top : (⊤ : (Spec Γ(X, A)).Opens) ≤ ⊤)) = 𝟙 _ from rfl]
    simp
  rw [show (Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv
      = (Spec.map (f.appLE B A i)).appLE ⊤ ⊤ e3 by
    rw [specMap_appLE (f.appLE B A i) ⊤ e3, hid]; simp]
  rw [Scheme.Hom.appLE_comp_appLE]
  exact (square_appLE_gen f hA hB i _ _ _ _).symm

/-- ★★★★★★★**`hcompat`**——第 221 ブロックが仮定していた等式。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが B2(`ofDivisor_pullback`)の最後の 1 点であった。 -/
theorem hcompat_image {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (hle : (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens))
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))) :
    f.appLE (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))
        (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)) hle
      = (hB.fromSpec.appIso ⊤).hom ≫ ((Scheme.ΓSpecIso (Γ(Y, B))).hom
          ≫ f.appLE B A i ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv)
          ≫ (hA.fromSpec.appIso ⊤).inv := by
  conv_rhs => rw [← Category.assoc]
  rw [← hcompat_key f hA hB i hle, Category.assoc, Iso.hom_inv_id, Category.comp_id]

variable (D : Y.IdealSheafData)

/-- ★★★★★★**比較射はアフィンの対で全射である**(仮定なし)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 221 に第 235(本ブロック)を渡しただけで、仮定が消えた。 -/
theorem surjective_affine {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (hle : (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens))
      ≤ f ⁻¹ᵁ (hB.fromSpec ''ᵁ (⊤ : (Spec Γ(Y, B)).Opens))) :
    Function.Surjective (((pullIdealHom f D).app
      (op (hA.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A)).Opens)))).hom) :=
  surjective_pair f D hA hB i hle (hcompat_image f hA hB i hle)

/-! ## ★出典の紐付け(`.src`) -/

def hcompat_image.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——hcompat が閉じた)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
