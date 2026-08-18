import ABC3.Found.Arakelov.PicLocalGlobal
import ABC3.Found.Arakelov.PicComapIso

/-!
# Arakelov (B2) 第 200 ブロック —— ★★★★★★**Cartier 性は点ごとに確かめれば足りる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★一般のスキームで局所 ⟹ 大域が出た

    ∀ x ∈ X, ∃ アフィン開 A ∋ x, D.ideal A が可逆
      ⟹  IsCartier X D    （= ∀ アフィン開で可逆）

★★★第 198(可逆性の Zariski 局所性)は `Spec R` の話であった。
本ブロックはそれを**任意のスキームの任意のアフィン開**へ運ぶ。

## ★★運び方は `fromSpec` と「同型に沿った引き戻し」

アフィン開 `A` に対し `j := A.2.fromSpec : Spec Γ(X,A) ⟶ X` は**開埋め込み**で、
mathlib の `ideal_comap_of_isOpenImmersion` が

    (D.comap j).ideal V = (D.ideal (j ''ᵁ V)).comap ((j.appIso V).inv.hom)

を与える——★**環の同型に沿った引き戻し**である。したがって第 199 で可逆性が往復でき、
`Spec Γ(X,A)` の上で第 198 を使ってから `A` に戻せばよい。

★`j ''ᵁ ⊤ = A.1` は `opensRange_fromSpec` で出るが、**依存位置**にあるので
`rw` では動かない——`X.affineOpens` 上の motive で `cast` する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pointwise_comap` | ★点ごとの仮定は `Spec Γ(X,A)` へ運べる |
| `comap_comap_of_iso` | ★同型で往復すると元に戻る |
| `invertible_of_comap_isIso` | ★★第 199 の逆向き |
| `isCartier_of_pointwise` | ★★★★★★**Cartier 性は点ごとで足りる** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

/-- ★**点ごとの仮定は `Spec Γ(X,A)` へ運べる**。 -/
theorem pointwise_comap (A : X.affineOpens)
    (h : ∀ x : X, ∃ A : X.affineOpens, ∃ _ : x ∈ A.1,
      Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)) :
    ∀ y : Spec Γ(X, A.1), ∃ V : (Spec Γ(X, A.1)).affineOpens, ∃ _ : y ∈ V.1,
      Module.Invertible (Γ(Spec Γ(X, A.1), V.1) : Type u)
        ((D.comap A.2.fromSpec).ideal V) := by
  intro y
  obtain ⟨A', hxA', hinv'⟩ := h (A.2.fromSpec.base y)
  obtain ⟨V, hV, hyV, hVsub⟩ := Opens.isBasis_iff_nbhd.1
    (Spec Γ(X, A.1)).isBasis_affineOpens
    (U := (Opens.map A.2.fromSpec.base).obj A'.1) (x := y) hxA'
  refine ⟨⟨V, hV⟩, hyV, ?_⟩
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion]
  haveI : Module.Invertible
      (Γ(X, (A.2.fromSpec ''ᵁ V)) : Type u)
      (D.ideal ⟨A.2.fromSpec ''ᵁ V, hV.image_of_isOpenImmersion _⟩) := by
    refine invertible_ideal_of_le D
      (A := ⟨A.2.fromSpec ''ᵁ V, hV.image_of_isOpenImmersion _⟩) (B := A') ?_ hinv'
    rintro z ⟨w, hw, rfl⟩
    exact hVsub hw
  exact invertible_comap_of_isIso _ _




theorem comap_comap_of_iso {R S : CommRingCat.{u}} (e : R ≅ S) (I : Ideal (S : Type u)) :
    (I.comap e.hom.hom).comap e.inv.hom = I := by
  ext x
  simp only [Ideal.mem_comap]
  constructor
  · intro hx
    have : e.hom.hom (e.inv.hom x) = x :=
      congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) x) e.inv_hom_id
    rwa [this] at hx
  · intro hx
    have : e.hom.hom (e.inv.hom x) = x :=
      congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) x) e.inv_hom_id
    rwa [this]

theorem invertible_of_comap_isIso {R S : CommRingCat.{u}} (e : R ≅ S) (I : Ideal (S : Type u))
    [Module.Invertible (R : Type u) (I.comap e.hom.hom)] :
    Module.Invertible (S : Type u) I := by
  have h := invertible_comap_of_isIso e.symm (I.comap e.hom.hom)
  rwa [show (e.symm : S ≅ R).hom = e.inv from rfl, comap_comap_of_iso e I] at h




theorem isCartier_of_pointwise
    (h : ∀ x : X, ∃ A : X.affineOpens, ∃ _ : x ∈ A.1,
      Module.Invertible (Γ(X, A.1) : Type u) (D.ideal A)) :
    IsCartier X D := by
  intro A
  have hAeq : A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens) = A.1 := by
    simp [Scheme.Hom.image_top_eq_opensRange, A.2.opensRange_fromSpec]
  have key := invertible_ideal_top_of_pointwise (Γ(X, A.1)) (D.comap A.2.fromSpec)
    (pointwise_comap D A h)
  rw [Scheme.IdealSheafData.ideal_comap_of_isOpenImmersion] at key
  have hB : (⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
      (isAffineOpen_top _).image_of_isOpenImmersion _⟩ : X.affineOpens) = A := Subtype.ext hAeq
  haveI : Module.Invertible (Γ(Spec Γ(X, A.1), ⊤) : Type u)
      ((D.ideal ⟨A.2.fromSpec ''ᵁ (⊤ : (Spec Γ(X, A.1)).Opens),
          (isAffineOpen_top _).image_of_isOpenImmersion _⟩).comap
        ((A.2.fromSpec.appIso ⊤).symm.hom.hom)) := key
  exact cast (congrArg (fun B : X.affineOpens =>
      Module.Invertible (Γ(X, B.1) : Type u) (D.ideal B)) hB)
    (invertible_of_comap_isIso (A.2.fromSpec.appIso ⊤).symm _)

/-! ## ★出典の紐付け(`.src`) -/

def isCartier_of_pointwise.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——Cartier 性は点ごとで足りる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
