import ABC3.Found.Arakelov.PicFromSpecTransport
import ABC3.Found.Arakelov.PicIdealFlat

/-!
# Arakelov (B2) 第 202 ブロック —— ★★★★★★**Cartier は平坦な引き戻しで保たれる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★`CartierPicData` の 13 欄目

    Flat f → IsCartier Y D → IsCartier X (D.comap f)

★§9-206 で「平坦性なしでは**偽**」と反例つきで記録した欄である。

## ★★筋は「アフィンの対に落として平坦底変換」

`x` ごとに `f(x)` の近傍アフィン `B`、`f⁻¹B` の中のアフィン `A ∋ x` を取る。すると

    (D.comap f)|_A  ≅  (D|_B) を f.appLE B A で拡大したもの

であり、`f.appLE` は**平坦**(`Scheme.Hom.flat_appLE`)だから第 183 が当たる。
★あとは第 200(点ごとで足りる)で大域へ。

## ★★可換性は mathlib に在った

    Spec.map (f.appLE B A i) ≫ hB.fromSpec = hA.fromSpec ≫ f

これ(`IsAffineOpen.SpecMap_appLE_fromSpec`)と `comap_comp` で

    (D.comap f).comap hA.fromSpec = (D.comap hB.fromSpec).comap (Spec.map (f.appLE B A i))

が出る。★`(Spec.map φ).appTop` は `ΓSpecIso_naturality` で
`e₁ ≫ φ ≫ e₂` に分解でき、`Ideal.map_map` で 3 段に割れる
——両端は同型(第 201)、中央は平坦(第 183)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appTop_decomp` | ★`(Spec.map φ).appTop` の分解 |
| `comap_decomp` | ★`comap` をアフィンの対へ落とす |
| `invertible_comap_pair` | ★★★★アフィンの対での可逆性 |
| `isCartier_comap` | ★★★★★★**Cartier は平坦な引き戻しで保たれる** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (D : Y.IdealSheafData)

theorem appTop_decomp {A : X.Opens} {B : Y.Opens} (i : A ≤ f ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i)).appTop
      = (Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ (f.appLE B A i)
        ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv := by
  rw [← Scheme.ΓSpecIso_naturality]
  simp

theorem comap_decomp {A : X.Opens} {B : Y.Opens} (hA : IsAffineOpen A) (hB : IsAffineOpen B)
    (i : A ≤ f ⁻¹ᵁ B) :
    (D.comap f).comap hA.fromSpec
      = (D.comap hB.fromSpec).comap (Spec.map (f.appLE B A i)) := by
  rw [← Scheme.IdealSheafData.comap_comp, ← Scheme.IdealSheafData.comap_comp,
    IsAffineOpen.SpecMap_appLE_fromSpec]




theorem invertible_comap_pair {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    [AlgebraicGeometry.Flat f]
    (hDB : Module.Invertible (Γ(Y, B) : Type u) (D.ideal ⟨B, hB⟩)) :
    Module.Invertible (Γ(X, A) : Type u) ((D.comap f).ideal ⟨A, hA⟩) := by
  refine invertible_of_comap_fromSpec (D.comap f) ⟨A, hA⟩ ?_
  haveI hI := invertible_comap_fromSpec D ⟨B, hB⟩ hDB
  rw [comap_decomp f D hA hB i, comap_ideal_top, appTop_decomp f i]
  letI : Algebra (Γ(Y, B) : Type u) (Γ(X, A) : Type u) := (f.appLE B A i).hom.toAlgebra
  haveI : Module.Flat (Γ(Y, B) : Type u) (Γ(X, A) : Type u) :=
    Scheme.Hom.flat_appLE f hB hA i
  haveI h1 : Module.Invertible (Γ(Y, B) : Type u)
      (((D.comap hB.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩).map
        (Scheme.ΓSpecIso (Γ(Y, B))).hom.hom) :=
    invertible_map_of_isIso _ _
  haveI h2 := invertible_map_of_flat (R := (Γ(Y, B) : Type u)) (S := (Γ(X, A) : Type u))
    (I := ((D.comap hB.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩).map
      (Scheme.ΓSpecIso (Γ(Y, B))).hom.hom)
  rw [RingHom.algebraMap_toAlgebra] at h2
  haveI := h2
  have h3 := invertible_map_of_isIso (Scheme.ΓSpecIso (Γ(X, A))).symm
    ((((D.comap hB.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩).map
      (Scheme.ΓSpecIso (Γ(Y, B))).hom.hom).map (f.appLE B A i).hom)
  rw [CommRingCat.hom_comp, CommRingCat.hom_comp, ← Ideal.map_map, ← Ideal.map_map]
  exact h3



theorem isCartier_comap {X Y : Scheme.{u}} (f : X ⟶ Y) [AlgebraicGeometry.Flat f]
    (D : Y.IdealSheafData) (hD : IsCartier Y D) : IsCartier X (D.comap f) := by
  refine isCartier_of_pointwise _ (fun x => ?_)
  obtain ⟨B, hB, hfxB, -⟩ := Opens.isBasis_iff_nbhd.1 Y.isBasis_affineOpens
    (U := ⊤) (x := f.base x) trivial
  obtain ⟨A, hA, hxA, hAsub⟩ := Opens.isBasis_iff_nbhd.1 X.isBasis_affineOpens
    (U := (Opens.map f.base).obj B) (x := x) hfxB
  exact ⟨⟨A, hA⟩, hxA, invertible_comap_pair f D hA hB hAsub (hD ⟨B, hB⟩)⟩

/-! ## ★出典の紐付け(`.src`) -/

def isCartier_comap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Cartier は平坦な引き戻しで保たれる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
