import ABC3.Found.Arakelov.PicDivisorTop

/-!
# Arakelov (B2) 第 188 ブロック —— **因子の積は切断の積**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`ofDivisor_mul` に向けた射

    idealPresheaf D ⊗ idealPresheaf E ⟶ idealPresheaf (D * E)

を `s ⊗ t ↦ s·t` で作り、★**アフィン開の上では全単射**であることを示す。

## ★★アフィン開では第 185 がそのまま当たる

アフィン開 `A` では `idealSections D A = D.ideal A`(第 148)なので、この射は
環レベルの掛け算写像 `mulTensorMap`(第 185)そのものである:

| 向き | 根拠 |
|---|---|
| 単射 | ★`D.ideal A` の**平坦性**(第 185) |
| 全射 | ★像が `D.ideal A * E.ideal A = (D*E).ideal A`(第 185) |

★★★`(D*E).ideal A = D.ideal A * E.ideal A` は mathlib で **`rfl`** なので繋がる。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `mul_mem_idealSections` | ★切断の積は積因子の切断 |
| `mulBil` / `mulHom` | ★★★積の射 |
| `idealSections_mul_affine` | ★アフィンでは切断も積 |
| `invertible_idealSections` | ★Cartier なら切断加群は可逆 |
| `sVal` / `mulHom_app_val` | ★環レベルの写像との一致 |
| `bijective_mulHom_app` | ★★★★**アフィン開の上では全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
variable {X : Scheme.{u}} (D E : X.IdealSheafData)

/-- ★切断の積は積因子の切断である。 -/
theorem mul_mem_idealSections {U : X.Opens} {s t : (Γ(X, U) : Type u)}
    (hs : s ∈ idealSections D U) (ht : t ∈ idealSections E U) :
    s * t ∈ idealSections (D * E) U := by
  intro V h
  have := Ideal.mul_mem_mul (hs V h) (ht V h)
  simpa [map_mul] using this

/-- ★★積の双線型写像。 -/
noncomputable def mulBil (U : (X.Opens)ᵒᵖ) :
    ((idealPresheaf D).obj U : Type u) →ₗ[(Γ(X, U.unop) : Type u)]
      ((idealPresheaf E).obj U : Type u) →ₗ[(Γ(X, U.unop) : Type u)]
        ((idealPresheaf (D * E)).obj U : Type u) where
  toFun s :=
    { toFun := fun t => ⟨s.1 * t.1, mul_mem_idealSections D E s.2 t.2⟩
      map_add' := fun a b => Subtype.ext (mul_add _ _ _)
      map_smul' := fun c t => Subtype.ext (by
        show s.1 * (c * t.1) = c * (s.1 * t.1)
        ring) }
  map_add' := fun a b => LinearMap.ext fun t => Subtype.ext (add_mul _ _ _)
  map_smul' := fun c s => LinearMap.ext fun t => Subtype.ext (by
    show (c * s.1) * t.1 = c * (s.1 * t.1)
    ring)



variable {X : Scheme.{u}} (D E : X.IdealSheafData)

set_option maxHeartbeats 1000000 in
/-- ★★★積の射。 -/
noncomputable def mulHom :
    (idealPresheaf D) ⊗ (idealPresheaf E) ⟶ idealPresheaf (D * E) where
  app U :=
    letI : Module (((X.presheaf ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U) : Type u)
        (↥(idealSections (D * E) U.unop)) := ((idealPresheaf (D * E)).obj U).isModule
    ModuleCat.ofHom (TensorProduct.lift (mulBil D E U))
  naturality := by
    intro U V f
    refine ModuleCat.MonoidalCategory.tensor_ext ?_
    intro s t
    erw [ConcreteCategory.comp_apply, ConcreteCategory.comp_apply]
    exact Subtype.ext (map_mul (X.presheaf.map f).hom s.1 t.1).symm



variable {X : Scheme.{u}} (D E : X.IdealSheafData)

theorem idealSections_mul_affine (A : X.affineOpens) :
    idealSections (D * E) A.1 = idealSections D A.1 * idealSections E A.1 := by
  rw [idealSections_affine, idealSections_affine, idealSections_affine]
  rfl

theorem invertible_idealSections (A : X.affineOpens) (hD : IsCartier X D) :
    Module.Invertible (Γ(X, A.1) : Type u) (idealSections D A.1) :=
  cast (congrArg (fun (J : Ideal (Γ(X, A.1) : Type u)) =>
    Module.Invertible (Γ(X, A.1) : Type u) J) (idealSections_affine D A).symm) (hD A)



variable {X : Scheme.{u}} (D E : X.IdealSheafData)
/-- ★イデアル層の切断を環の元と見る橋。 -/
def sVal {X : Scheme.{u}} {F : X.IdealSheafData} {U : X.Opens}
    (s : ((idealPresheaf F).obj (op U) : Type u)) : (Γ(X, U) : Type u) := s.1

/-- ★積の射の値は環レベルの掛け算写像と一致する。 -/
theorem mulHom_app_val (A : X.affineOpens)
    (x : ((((idealPresheaf D) ⊗ (idealPresheaf E)).obj (op A.1)) : Type u)) :
    sVal (((mulHom D E).app (op A.1)).hom x)
      = mulTensorMap (idealSections D A.1) (idealSections E A.1) x := by
  induction x using TensorProduct.induction_on with
  | zero =>
      exact (congrArg sVal (map_zero ((mulHom D E).app (op A.1)).hom)).trans
        (map_zero (mulTensorMap (idealSections D A.1) (idealSections E A.1))).symm
  | tmul s t =>
      exact (mulTensorMap_tmul (idealSections D A.1) (idealSections E A.1) s t).symm
  | add a b ha hb =>
      refine (congrArg sVal (map_add ((mulHom D E).app (op A.1)).hom a b)).trans ?_
      refine Eq.trans (congrArg₂ (· + ·) ha hb) ?_
      exact (map_add (mulTensorMap (idealSections D A.1) (idealSections E A.1)) a b).symm

theorem bijective_mulHom_app (A : X.affineOpens) (hD : IsCartier X D) :
    Function.Bijective (((mulHom D E).app (op A.1)).hom) := by
  haveI := invertible_idealSections D A hD
  constructor
  · intro x y hxy
    refine mulTensorMap_injective (I := idealSections D A.1) (J := idealSections E A.1) ?_
    rw [← mulHom_app_val, ← mulHom_app_val, hxy]
  · intro u
    have hu : sVal u ∈ idealSections D A.1 * idealSections E A.1 := by
      rw [← idealSections_mul_affine]
      exact u.2
    rw [← mulTensorMap_range] at hu
    obtain ⟨x, hx⟩ := hu
    exact ⟨x, Subtype.ext ((mulHom_app_val D E A x).trans hx)⟩


/-! ## ★出典の紐付け(`.src`) -/

def bijective_mulHom_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——因子の積は切断の積)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
