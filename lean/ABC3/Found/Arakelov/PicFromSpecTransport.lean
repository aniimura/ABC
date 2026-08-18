import ABC3.Found.Arakelov.PicCartierPt
import ABC3.Found.Arakelov.PicComapAffine

/-!
# Arakelov (B2) 第 201 ブロック —— **`fromSpec` に沿った可逆性の往復**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★最後の 2 欄に要る往復

第 200 の中で使った「`Spec Γ(X,A)` と `A` の間で可逆性を運ぶ」を、
**両向き**に切り出す。これがあれば残りは

    (D.comap f).ideal A  →(往)→  Spec Γ(X,A) の ⊤  →(第 195)→  平坦な拡大  →(復)→

と繋がる。

## ★★併せてイデアルの `map` 側も取る

第 199 は `comap` 側であった。`Ideal.map_comap_of_equiv`(mathlib)で
`I.map e = I.comap e⁻¹` なので、`map` 側も同じ補題から出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `comap_inv_eq_map` | ★同型なら `comap` と `map` は入れ替わる |
| `invertible_map_of_isIso` | ★★同型に沿った `map` でも可逆 |
| `invertible_of_comap_fromSpec` | ★★★`Spec Γ(X,A)` → `A`(復) |
| `invertible_comap_fromSpec` | ★★★`A` → `Spec Γ(X,A)`(往) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

section RingIso

variable {R S : CommRingCat.{u}} (e : R ≅ S) (I : Ideal (R : Type u))

/-- ★**同型なら `comap` と `map` は入れ替わる**。 -/
theorem comap_inv_eq_map : I.comap e.inv.hom = I.map e.hom.hom :=
  (Ideal.map_comap_of_equiv (I := I) (e.commRingCatIsoToRingEquiv)).symm

/-- ★★**同型に沿った `map` でも可逆である**。 -/
theorem invertible_map_of_isIso [Module.Invertible (R : Type u) I] :
    Module.Invertible (S : Type u) (I.map e.hom.hom) := by
  have h := invertible_comap_of_isIso e.symm I
  rwa [show (e.symm : S ≅ R).hom = e.inv from rfl, comap_inv_eq_map e I] at h

end RingIso

variable {X : Scheme.{u}} (E : X.IdealSheafData)

/-- ★★★**`Spec Γ(X,A)` から `A` へ可逆性を戻す**。 -/
theorem invertible_of_comap_fromSpec (A : X.affineOpens)
    (h : Module.Invertible (Γ(Spec Γ(X, A.1), ⊤) : Type u)
      ((E.comap A.2.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩)) :
    Module.Invertible (Γ(X, A.1) : Type u) (E.ideal A) := by
  have hAeq : A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens) = A.1 := by
    simp [Scheme.Hom.image_top_eq_opensRange, A.2.opensRange_fromSpec]
  have key := h
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion] at key
  have hB : (⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
      (isAffineOpen_top _).image_of_isOpenImmersion _⟩ : X.affineOpens) = A := Subtype.ext hAeq
  haveI : Module.Invertible (Γ(Spec Γ(X, A.1), ⊤) : Type u)
      ((E.ideal ⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
        ((A.2.fromSpec.appIso ⊤).symm.hom.hom)) := key
  exact cast (congrArg (fun C : X.affineOpens =>
      Module.Invertible (Γ(X, C.1) : Type u) (E.ideal C)) hB)
    (invertible_of_comap_isIso (A.2.fromSpec.appIso ⊤).symm _)

/-- ★逆向き——`A` から `Spec Γ(X,A)` へ運ぶ。 -/
theorem invertible_comap_fromSpec (A : X.affineOpens)
    (h : Module.Invertible (Γ(X, A.1) : Type u) (E.ideal A)) :
    Module.Invertible (Γ(Spec Γ(X, A.1), ⊤) : Type u)
      ((E.comap A.2.fromSpec).ideal ⟨⊤, isAffineOpen_top _⟩) := by
  have hAeq : A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens) = A.1 := by
    simp [Scheme.Hom.image_top_eq_opensRange, A.2.opensRange_fromSpec]
  have hB : (⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
      (isAffineOpen_top _).image_of_isOpenImmersion _⟩ : X.affineOpens) = A := Subtype.ext hAeq
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion]
  haveI : Module.Invertible (Γ(X, (A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens))) : Type u)
      (E.ideal ⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
        (isAffineOpen_top _).image_of_isOpenImmersion _⟩) :=
    cast (congrArg (fun C : X.affineOpens =>
      Module.Invertible (Γ(X, C.1) : Type u) (E.ideal C)) hB.symm) h
  exact invertible_comap_of_isIso _ _


/-! ## ★出典の紐付け(`.src`) -/

def invertible_of_comap_fromSpec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——fromSpec に沿った可逆性の往復)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
