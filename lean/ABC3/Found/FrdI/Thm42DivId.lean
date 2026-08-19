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
    rfl
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

end PullBackOtri

end ABC3.Found.FrdI
