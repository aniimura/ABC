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

end PullBackOtri

end ABC3.Found.FrdI
