/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm42PsiPrime

/-!
# [FrdI] Theorem 4.2, (i) —— `Div-identity` 自己射の一般形へ

原文 (FrdI p.81):
> is a Div-identity endomorphism if and only if it admits a factorization

★★測って立てた**圏論的特徴づけ**(原文の 2 段可換図式の中身):

`α = γ ≫ β`(`γ` は base-isomorphism、`β` は pull-back)と割ると
`Φ.map (Base α) = Φ.map (Base γ) ∘ Φ.map (Base β)` である。
★そして
- `β` に沿った四角形 `β ≫ v = w' ≫ β` は `Φ.map (Base β) (Div v) = Div w'` を与える
  (`β` は pull-back なので linear)
- `γ` に沿った四角形 `w ≫ γ = γ ≫ w'` は
  `Φ.map (Base γ) (Div w') = n • Div w`(`n = degFr γ`)を与える

★★したがって **`IsDivIdentity α` ⟺ すべての `v ∈ 𝒪^▷(A)` について `n • Div w = Div v`**
であり、`n • Div w = Div (w^n)` と `exists_iso_of_div_eq`(コスライス)により
これは**圏論的**な条件になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w u2 v2

section PullBackOtri

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**pull-back に沿って `𝒪^▷` を引き戻す** —— `Definition 1.2, (ii)` の普遍性から。 -/
theorem otriPullBack (Q : PreFrobenioid C Φ) {A B : C} (β : B ⟶ A) (hβ : IsPullBack Q β)
    {v : End A} (hv : v ∈ OTri Q A) :
    ∃ w : End B, w ∈ OTri Q B ∧ β ≫ (((v : End A)) : A ⟶ A) = (((w : End B)) : B ⟶ B) ≫ β := by
  obtain ⟨w, ⟨hwb, hwsq⟩, -⟩ := prop_1_11_iii Q β hβ v (𝟙 _) (by
    rw [show Q.Base (((v : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hv.1, Q.Base_id,
      Category.comp_id, Category.id_comp])
  refine ⟨w, ⟨?_, ?_⟩, hwsq⟩
  · show Q.Base (((w : End B)) : B ⟶ B) = Q.Base (𝟙 B)
    rw [hwb, Q.Base_id]
  · have h := congrArg Q.degFr hwsq
    rw [Q.degFr_comp, Q.degFr_comp,
      show Q.degFr (((v : End A)) : A ⟶ A) = 1 from hv.2, one_mul] at h
    show Q.degFr (((w : End B)) : B ⟶ B) = 1
    exact (mul_left_cancel (a := Q.degFr β) (by rw [mul_one]; exact h)).symm

/-- ★★**pull-back の四角形の `Div`** —— `β` は linear なので `Φ.map (Base β) (Div v) = Div w`。 -/
theorem div_square_pullBack (Q : PreFrobenioid C Φ) {A B : C} (β : B ⟶ A)
    (hβl : Q.degFr β = 1)
    {w : End B} (hw : w ∈ OTri Q B) {v : End A} (hv : v ∈ OTri Q A)
    (hsq : β ≫ (((v : End A)) : A ⟶ A) = (((w : End B)) : B ⟶ B) ≫ β) :
    Φ.map (Q.Base β) (Q.Div (((v : End A)) : A ⟶ A))
      = Q.Div (((w : End B)) : B ⟶ B) := by
  have h := div_square_frob Q β hw hv hsq.symm
  rw [hβl] at h
  rw [h, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]

def otriPullBack.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 81,
    item := "Theorem 4.2, (i) — pull-back に沿った 𝒪^▷ の引き戻し",
    sectionId := "frdi-thm-4-2" }

/-! ## ★★★★★★`Div-identity` の圏論的特徴づけ -/

/-- ★★★★**`Div-identity` の圏論版** ——
`Definition 1.3, (iv), (a)` の 3 分解に沿って `𝒪^▷` を運び、
最後に `n = degFr γ` 乗と同型で一致することを要求する。 -/
def DivIdCat (P : PreFrobenioid C Φ) {A : C} (α : A ⟶ A) : Prop :=
  ∃ (X Y : C) (γ : A ⟶ X) (pre : X ⟶ Y) (plb : Y ⟶ A),
    α = γ ≫ pre ≫ plb ∧ IsFrobeniusType P γ ∧ IsPreStep P pre ∧ IsPullBack P plb ∧
    ∀ (v : End A) (v₁ : End Y) (v₂ : End X) (v₃ : End A) (ε : End X),
      v ∈ OTri P A → v₁ ∈ OTri P Y → v₂ ∈ OTri P X → v₃ ∈ OTri P A → ε ∈ OTimes P X →
      plb ≫ (((v : End A)) : A ⟶ A) = (((v₁ : End Y)) : Y ⟶ Y) ≫ plb →
      pre ≫ (((v₁ : End Y)) : Y ⟶ Y) = (((v₂ : End X)) : X ⟶ X) ≫ pre →
      (((v₃ : End A)) : A ⟶ A) ≫ γ
        = γ ≫ ((((v₂ : End X)) : X ⟶ X) ≫ (((ε : End X)) : X ⟶ X)) →
        ∃ θ : A ⟶ A, IsIso θ ∧
          ((((v₃ ^ (((P.degFr γ : ℕ+) : ℕ)) : End A)) : A ⟶ A)) ≫ θ
            = (((v : End A)) : A ⟶ A)

/-- ★★`𝒪^▷` の冪の `Div` は倍数。 -/
theorem div_pow_otri (Q : PreFrobenioid C Φ) {A : C} {v : End A} (hv : v ∈ OTri Q A) (n : ℕ) :
    Q.Div ((((v ^ n : End A)) : A ⟶ A)) = n • Q.Div ((((v : End A)) : A ⟶ A)) := by
  induction n with
  | zero =>
    show Q.Div ((1 : End A) : A ⟶ A) = (0 : ℕ) • _
    rw [zero_nsmul]
    exact Q.Div_id A
  | succ m ih =>
    have hvm : (v ^ m) ∈ OTri Q A := pow_mem hv m
    show Q.Div ((((v ^ (m + 1) : End A)) : A ⟶ A)) = _
    rw [pow_succ]
    show Q.Div ((((v : End A)) : A ⟶ A) ≫ (((v ^ m : End A)) : A ⟶ A)) = _
    rw [Q.Div_comp, show Q.Base ((((v : End A)) : A ⟶ A)) = Q.Base (𝟙 A) from hv.1,
      Q.Base_id, show Q.degFr ((((v ^ m : End A)) : A ⟶ A)) = 1 from hvm.2]
    have h3 : Φ.map (𝟙 ((Q.toElem.obj A).base)) (Q.Div ((((v ^ m : End A)) : A ⟶ A)))
      = Q.Div ((((v ^ m : End A)) : A ⟶ A)) := Φ.map_id _ _
    rw [h3, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul, ih, succ_nsmul]

/-- ★`𝒪^×` を後置しても `Div` は変わらない。 -/
theorem div_otimes_comp (Q : PreFrobenioid C Φ) {A : C} {v : End A} (hv : v ∈ OTri Q A)
    {ε : End A} (hε : ε ∈ OTimes Q A) :
    Q.Div ((((v : End A)) : A ⟶ A) ≫ (((ε : End A)) : A ⟶ A))
      = Q.Div ((((v : End A)) : A ⟶ A)) := by
  haveI : IsIso ((((ε : End A)) : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso ε).mp hε.2
  exact div_comp_iso ((((v : End A)) : A ⟶ A)) ((((ε : End A)) : A ⟶ A))

def DivIdCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 81,
    item := "Theorem 4.2, (i) — Div-identity の圏論的特徴づけ",
    sectionId := "frdi-thm-4-2" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`IsDivIdentity` ⟹ `DivIdCat`**。 -/
theorem divIdCat_of_isDivIdentity (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    {A : C} (α : A ⟶ A) (h : IsDivIdentity P α) : DivIdCat P α := by
  obtain ⟨X, Y, γ, pre, plb, hfac, hγ, hpre, hplb⟩ := F.arbFactor α
  refine ⟨X, Y, γ, pre, plb, hfac, hγ, hpre, hplb, ?_⟩
  intro v v₁ v₂ v₃ ε hv hv₁ hv₂ hv₃ hε hsq1 hsq2 hsq3
  have hplbl : P.degFr plb = 1 := (F.pullBackLB plb hplb).2
  have e1 : Φ.map (P.Base plb) (P.Div ((((v : End A)) : A ⟶ A)))
      = P.Div ((((v₁ : End Y)) : Y ⟶ Y)) :=
    div_square_pullBack P plb hplbl hv₁ hv hsq1
  have e2 : Φ.map (P.Base pre) (P.Div ((((v₁ : End Y)) : Y ⟶ Y)))
      = P.Div ((((v₂ : End X)) : X ⟶ X)) :=
    map_base_div_otri P pre hpre.1 hv₁ hv₂ hsq2
  set v₂' : End X := (ε : End X) * (v₂ : End X) with hv₂'def
  have hv₂'m : v₂' ∈ OTri P X := (OTri P X).mul_mem (OTimes_le_OTri P X hε) hv₂
  have hsq3' : (((v₃ : End A)) : A ⟶ A) ≫ γ = γ ≫ (((v₂' : End X)) : X ⟶ X) := hsq3
  have e3 : Φ.map (P.Base γ) (P.Div ((((v₂' : End X)) : X ⟶ X)))
      = ((P.degFr γ : ℕ+) : ℕ) • P.Div ((((v₃ : End A)) : A ⟶ A)) :=
    div_square_frob P γ hv₃ hv₂'m hsq3'
  have e4 : P.Div ((((v₂' : End X)) : X ⟶ X)) = P.Div ((((v₂ : End X)) : X ⟶ X)) :=
    div_otimes_comp P hv₂ hε
  -- ★合成して `Φ.map (Base α) (Div v) = n • Div v₃`
  have hbase : P.Base α = P.Base γ ≫ P.Base pre ≫ P.Base plb := by
    rw [hfac, P.Base_comp, P.Base_comp]
  have ekey : Φ.map (P.Base α) (P.Div ((((v : End A)) : A ⟶ A)))
      = ((P.degFr γ : ℕ+) : ℕ) • P.Div ((((v₃ : End A)) : A ⟶ A)) := by
    rw [hbase, ← e3, e4, ← e2, ← e1]
    rw [Φ.map_comp, Φ.map_comp]
  -- ★`IsDivIdentity` から `= Div v`
  have hid : Φ.map (P.Base α) (P.Div ((((v : End A)) : A ⟶ A)))
      = P.Div ((((v : End A)) : A ⟶ A)) := by
    have h1 : Φ.map (P.Base α) = Φ.map (P.Base (𝟙 A)) := h
    rw [h1, P.Base_id]
    exact Φ.map_id _ _
  have hdvv : P.Div ((((v₃ ^ (((P.degFr γ : ℕ+) : ℕ)) : End A)) : A ⟶ A))
      = P.Div ((((v : End A)) : A ⟶ A)) := by
    rw [div_pow_otri P hv₃, ← ekey, hid]
  obtain ⟨θ, hθiso, hθ⟩ := exists_iso_of_div_eq G
    ((((v₃ ^ (((P.degFr γ : ℕ+) : ℕ)) : End A)) : A ⟶ A)) ((((v : End A)) : A ⟶ A))
    (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ (pow_mem hv₃ _))
    (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ hv) hdvv
  exact ⟨θ, hθiso, hθ⟩

def divIdCat_of_isDivIdentity.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 81,
    item := "Theorem 4.2, (i) — IsDivIdentity ⟹ 圏論的条件",
    sectionId := "frdi-thm-4-2" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`DivIdCat` ⟹ `IsDivIdentity`**。 -/
theorem isDivIdentity_of_divIdCat (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hperfM : ∀ Y : C, IsPerfectMonoid (Φ.val (P.toElem.obj Y).base))
    {A : C} (α : A ⟶ A) (h : DivIdCat P α) : IsDivIdentity P α := by
  obtain ⟨X, Y, γ, pre, plb, hfac, hγ, hpre, hplb, hcond⟩ := h
  have hplbl : P.degFr plb = 1 := (F.pullBackLB plb hplb).2
  haveI hγi : IsIso (P.Base γ) := hγ.2
  have hbase : P.Base α = P.Base γ ≫ P.Base pre ≫ P.Base plb := by
    rw [hfac, P.Base_comp, P.Base_comp]
  show Φ.map (P.Base α) = Φ.map (P.Base (𝟙 A))
  rw [P.Base_id]
  refine AddMonoidHom.ext (fun x => ?_)
  rw [show Φ.map (𝟙 ((P.toElem.obj A).base)) x = x from Φ.map_id _ _]
  obtain ⟨v, hvx⟩ := hdivS A x
  have hv : (v : End A) ∈ OTri P A := v.2
  subst hvx
  obtain ⟨v₁, hv₁, hsq1⟩ := otriPullBack P plb hplb hv
  have hprec : IsCoAngular P pre := prop_1_4_i P _ (fun Z _ => hiso Z)
  set v₂ := otriPull P F pre hprec hpre.1 ⟨v₁, hv₁⟩ with hv₂def
  have hv₂ : (v₂ : End X) ∈ OTri P X := v₂.2
  have hsq2 : pre ≫ (((v₁ : End Y)) : Y ⟶ Y) = (((v₂ : End X)) : X ⟶ X) ≫ pre :=
    otriPull_spec P F pre hprec hpre.1 ⟨v₁, hv₁⟩
  obtain ⟨d, hdd⟩ := (hperfM A (P.degFr γ)).2
    (Φ.map (P.Base γ) (P.Div ((((v₂ : End X)) : X ⟶ X))))
  have hddb : ((P.degFr γ : ℕ+) : ℕ) • d
      = Φ.map (P.Base γ) (P.Div ((((v₂ : End X)) : X ⟶ X))) := hdd
  obtain ⟨v₃, hv₃d⟩ := hdivS A d
  have hv₃ : (v₃ : End A) ∈ OTri P A := v₃.2
  obtain ⟨u₀, hsq0, -⟩ :=
    prop_1_10_i_exists_given P F ((((v₃ : End A)) : A ⟶ A)) γ hγ γ hγ rfl
  have hb0 : P.Base ((((v₃ : End A)) : A ⟶ A)) ≫ P.Base γ = P.Base γ ≫ P.Base u₀ := by
    rw [← P.Base_comp, ← P.Base_comp, hsq0]
  have hu₀b : IsBaseIdentity P u₀ := by
    show P.Base u₀ = P.Base (𝟙 X)
    rw [P.Base_id]
    rw [show P.Base ((((v₃ : End A)) : A ⟶ A)) = P.Base (𝟙 A) from hv₃.1, P.Base_id,
      Category.id_comp] at hb0
    exact ((cancel_epi (P.Base γ)).mp (by rw [Category.comp_id]; exact hb0)).symm
  have hu₀l : IsLinear P u₀ := by
    have hdd2 : P.degFr ((((v₃ : End A)) : A ⟶ A) ≫ γ) = P.degFr (γ ≫ u₀) := by rw [hsq0]
    rw [P.degFr_comp, P.degFr_comp,
      show P.degFr ((((v₃ : End A)) : A ⟶ A)) = 1 from hv₃.2, mul_one] at hdd2
    show P.degFr u₀ = 1
    exact (mul_right_cancel (b := P.degFr γ) (by rw [one_mul]; exact hdd2)).symm
  have hu₀m : u₀ ∈ OTri P X := ⟨hu₀b, hu₀l⟩
  have hdiv0 : Φ.map (P.Base γ) (P.Div u₀)
      = ((P.degFr γ : ℕ+) : ℕ) • P.Div ((((v₃ : End A)) : A ⟶ A)) :=
    div_square_frob P γ hv₃ hu₀m hsq0
  have hdu : P.Div u₀ = P.Div ((((v₂ : End X)) : X ⟶ X)) := by
    refine (Φ.map_bijective_of_iso (@asIso _ _ _ _ (P.Base γ) hγi)).1 ?_
    show Φ.map (P.Base γ) (P.Div u₀)
      = Φ.map (P.Base γ) (P.Div ((((v₂ : End X)) : X ⟶ X)))
    rw [hdiv0, hv₃d]
    exact hdd
  obtain ⟨ε, hεm, hεeq⟩ := F.faithfulUpToUnits u₀ ((((v₂ : End X)) : X ⟶ X))
    (show P.Base u₀ = P.Base ((((v₂ : End X)) : X ⟶ X)) by
      rw [show P.Base u₀ = P.Base (𝟙 X) from hu₀b,
        show P.Base ((((v₂ : End X)) : X ⟶ X)) = P.Base (𝟙 X) from hv₂.1])
    hdu (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ hu₀m)
    (prop_1_4_i P _ (fun Z _ => hiso Z)) (isPreStep_of_otri _ hv₂)
  have hsq3 : ((((v₃ : End A)) : A ⟶ A)) ≫ γ
      = γ ≫ (((((v₂ : End X)) : X ⟶ X)) ≫ ((((ε : End X)) : X ⟶ X))) := by
    rw [hsq0, ← hεeq]
  obtain ⟨θ, hθiso, hθ⟩ := hcond (v : End A) v₁ (v₂ : End X) (v₃ : End A) ε
    hv hv₁ hv₂ hv₃ hεm hsq1 hsq2 hsq3
  haveI := hθiso
  have hdvv : P.Div ((((v : End A)) : A ⟶ A))
      = ((P.degFr γ : ℕ+) : ℕ) • P.Div ((((v₃ : End A)) : A ⟶ A)) := by
    rw [← div_pow_otri P hv₃, ← hθ]
    exact div_comp_iso (P₂ := P) ((((v₃ ^ (((P.degFr γ : ℕ+) : ℕ)) : End A)) : A ⟶ A)) θ
  have e1 : Φ.map (P.Base plb) (P.Div ((((v : End A)) : A ⟶ A)))
      = P.Div ((((v₁ : End Y)) : Y ⟶ Y)) :=
    div_square_pullBack P plb hplbl hv₁ hv hsq1
  have e2 : Φ.map (P.Base pre) (P.Div ((((v₁ : End Y)) : Y ⟶ Y)))
      = P.Div ((((v₂ : End X)) : X ⟶ X)) :=
    map_base_div_otri P pre hpre.1 hv₁ hv₂ hsq2
  rw [hbase, Φ.map_comp, Φ.map_comp, e1, e2, ← hddb, ← hv₃d, ← hdvv]

def isDivIdentity_of_divIdCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 81,
    item := "Theorem 4.2, (i) — 圏論的条件 ⟹ IsDivIdentity",
    sectionId := "frdi-thm-4-2" }

end PullBackOtri

end ABC3.Found.FrdI
