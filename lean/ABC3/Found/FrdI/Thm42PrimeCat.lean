/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42Prime
import ABC3.Found.FrdI.Prop41Cat

/-!
# [FrdI] Theorem 4.2, (ii) —— `Prime(Φ(A))` を **primary step の同値類**として復元する

原文 (FrdI p.79):
> class of primary steps to or from Ai [where the correspondence between elements of

★★`Thm42Prime.lean` の `suppElt_eq_iff_exists_common` は
「primary な 2 元の台が一致する ⟺ 0 でない共通の下界がある」を与えた。
★★★本ファイルではその「共通の下界」を **圏論的な条件**に翻訳する ——
`Definition 1.3, (iii), (d)` の圏同値(在庫の `exists_factor_through` /
`exists_factor_of_mle'`)で、`MLe` は「共通の商 step を通る分解」になる。

## ★結論

**2 つの primary な co-angular pre-step `ϵ`、`ϵ'` が同じ素点を定める
⟺ 両方が通る co-angular **step** `ζ` が存在する。**

★右辺は `IsCoAngular`・`IsStep`・`IsPreStep`・合成だけで書かれているので、
`Theorem 3.4, (ii)` により **`Ψ` でそのまま移る**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w u2 v2 w'

section Cat

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★**「同じ素点」の圏論的条件** —— 両方が通る co-angular step があること。 -/
def SamePrimeCat (P : PreFrobenioid C Φ) {A B B' : C} (ϵ : B ⟶ A) (ϵ' : B' ⟶ A) : Prop :=
  ∃ (Z : C) (ζ : Z ⟶ A), IsCoAngular P ζ ∧ IsStep P ζ ∧
    (∃ ψ : B ⟶ Z, IsPreStep P ψ ∧ ϵ = ψ ≫ ζ) ∧
    (∃ ψ' : B' ⟶ Z, IsPreStep P ψ' ∧ ϵ' = ψ' ≫ ζ)

/-- ★★★★★★**[FrdI] Theorem 4.2, (ii) の中核** ——
2 つの primary な co-angular pre-step が同じ素点を定めることの**圏論的判定**。 -/
theorem suppElt_preStepVal_eq_iff_samePrimeCat (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A : C}
    {ι₀ : Prime (Φ.val (P.toElem.obj A).base) → Pf (Φ.val (P.toElem.obj A).base) → ℝ≥0}
    (H : IsPerfFactorialWith (Φ.val (P.toElem.obj A).base) ι₀)
    (hperf : IsPerfectMonoid (Φ.val (P.toElem.obj A).base))
    (hdivM : IsDivisorial (Φ.val (P.toElem.obj A).base))
    {B B' : C} (ϵ : B ⟶ A) (hcϵ : IsCoAngular P ϵ) (hsϵ : IsPreStep P ϵ)
    (hpϵ : IsPrimaryElt (preStepVal P ϵ hsϵ))
    (ϵ' : B' ⟶ A) (hcϵ' : IsCoAngular P ϵ') (hsϵ' : IsPreStep P ϵ')
    (hpϵ' : IsPrimaryElt (preStepVal P ϵ' hsϵ')) :
    SuppElt ι₀ (preStepVal P ϵ hsϵ) = SuppElt ι₀ (preStepVal P ϵ' hsϵ')
      ↔ SamePrimeCat P ϵ ϵ' := by
  rw [suppElt_eq_iff_exists_common H hperf hdivM hpϵ hpϵ']
  constructor
  · rintro ⟨d, hd0, hdϵ, hdϵ'⟩
    -- ★`d` を実現する商 step `χ` を `ϵ` の側から取り出す
    obtain ⟨Z, ψ, χ, hψs, hχs, hψc, hχc, hfac, hχval⟩ :=
      exists_factor_of_mle' G ϵ hcϵ hsϵ hdϵ
    -- ★`ϵ'` は `χ` を通る
    have hle' : MLe (preStepVal P χ hχs) (preStepVal P ϵ' hsϵ') := by
      rw [hχval]; exact hdϵ'
    obtain ⟨ψ', hψ'c, hψ's, hfac'⟩ :=
      exists_factor_through G χ ϵ' hχc hχs hcϵ' hsϵ' hle'
    refine ⟨Z, χ, hχc, ?_, ⟨ψ, hψs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩
    refine (isStep_iff_preStepVal_ne_zero hiso χ hχs).mpr ?_
    rw [hχval]; exact hd0
  · rintro ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩
    refine ⟨preStepVal P ζ hζst.1, (isStep_iff_preStepVal_ne_zero hiso ζ hζst.1).mp hζst, ?_, ?_⟩
    · subst hfac
      exact ⟨Φ.map (@inv _ _ _ _ (P.Base (ψ ≫ ζ)) hsϵ.2) (P.Div ψ),
        (preStepVal_comp ψ ζ hψs hζst.1 hsϵ).symm⟩
    · subst hfac'
      exact ⟨Φ.map (@inv _ _ _ _ (P.Base (ψ' ≫ ζ)) hsϵ'.2) (P.Div ψ'),
        (preStepVal_comp ψ' ζ hψ's hζst.1 hsϵ').symm⟩

def SamePrimeCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 77,
    item := "Theorem 4.2, (ii) — 同素性の圏論的判定",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★同素性は `Ψ` で移る -/

section Transport

variable {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★**`SamePrimeCat` は `Ψ` で移る** —— 右辺が純粋に圏論的だから。 -/
theorem samePrimeCat_map (Ψ : C ⥤ C₂)
    (hCA : ∀ {X Y : C} (f : X ⟶ Y), IsCoAngular P f → IsCoAngular P₂ (Ψ.map f))
    (hPS : ∀ {X Y : C} (f : X ⟶ Y), IsPreStep P f → IsPreStep P₂ (Ψ.map f))
    (hST : ∀ {X Y : C} (f : X ⟶ Y), IsStep P f → IsStep P₂ (Ψ.map f))
    {A B B' : C} {ϵ : B ⟶ A} {ϵ' : B' ⟶ A} (h : SamePrimeCat P ϵ ϵ') :
    SamePrimeCat P₂ (Ψ.map ϵ) (Ψ.map ϵ') := by
  obtain ⟨Z, ζ, hζc, hζst, ⟨ψ, hψs, hfac⟩, ⟨ψ', hψ's, hfac'⟩⟩ := h
  exact ⟨Ψ.obj Z, Ψ.map ζ, hCA ζ hζc, hST ζ hζst,
    ⟨Ψ.map ψ, hPS ψ hψs, by rw [hfac, Ψ.map_comp]⟩,
    ⟨Ψ.map ψ', hPS ψ' hψ's, by rw [hfac', Ψ.map_comp]⟩⟩

def samePrimeCat_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 79,
    item := "Theorem 4.2, (ii) — 同素性は Ψ で移る",
    sectionId := "frdi-thm-4-2" }

end Transport

end Cat

end ABC3.Found.FrdI
