/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Ref

/-!
# [FrdI] Remark 2.7.2 の 3 分解 —— 等長な場合は skeleton の仮定が要らない

原文 (FrdI p.52):
> is a base-identity pre-step endomorphism;

`Def27.lean` の `rem_2_7_2_factor` は **`𝒞` が base-trivial 型かつ skeleton**
を仮定する。原文も「WLOG `𝒞` は skeleton としてよい」と書く。

★★**しかし `Theorem 5.2, (iv)` が実際に使うのは `𝒞^Fr-tr` と `𝒞^birat` であり、
そこでは「すべての射が等長」**である。そしてその場合、
**skeleton も base-trivial 型も要らない**:

`φ : A → B`(`A, B ∈ Ob(𝒫)`)に対し
1. `α := 𝒫` の一意な持ち上げ(底 `Base φ`)。これは pull-back 射。
2. その普遍性で `g ≫ α = φ`、`Base g = 𝟙` なる `g : A → A` が一意に取れる。
3. ★**`g` は等長なので Frobenius 型**(isotropic 型なので co-angular は自動)。
4. `Definition 1.3, (ii)` の一意性で `F_n ≫ β = g` なる同型 `β` が取れる。

★一意性(`rem_2_7_2_uniq`)は**もともと skeleton を使っていない**ので、そのまま使える。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D] {S : BaseSection P}

/-- ★★★★**`Remark 2.7.2` の 3 分解** —— **すべての射が等長**なら
`skeleton` と `base-trivial 型`の仮定なしで出る。

原文 (FrdI p.52):
> where α is P-distinguished; β is a base-identity pre-step endomorphism; γ is F distinguished.
-/
theorem rem_2_7_2_factor_isometric (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {Fs : ℕ+ →* SectionEnd S} (hFs : IsFrobeniusSection S Fs)
    {A B : C} (hA : S.objP A) (hB : S.objP B) (φ : A ⟶ B) :
    ∃ (n : ℕ+) (β : A ⟶ A) (α : A ⟶ B),
      φ = (((Fs n).app ⟨A, hA⟩ : End A) : A ⟶ A) ≫ β ≫ α ∧
      IsBaseIdentity P β ∧ IsPreStep P β ∧ S.homP α := by
  have hαP := S.lift_homP hA hB (P.Base φ)
  have hαpb := S.isPullBack hαP
  have hdegα : P.degFr (S.lift hA hB (P.Base φ)) = 1 := (Fc.pullBackLB _ hαpb).2
  obtain ⟨g, hg1, hg2⟩ := pullBack_lift hαpb φ (P.Base (𝟙 A)) (by
    rw [P.Base_id, Category.id_comp, S.lift_base])
  have hgb : P.Base g = 𝟙 _ := by rw [hg2, P.Base_id]
  have hdeg : P.degFr g = P.degFr φ := by
    have h := congrArg P.degFr hg1
    rw [P.degFr_comp, hdegα, one_mul] at h
    exact h
  have hgF : IsFrobeniusType P g :=
    ⟨⟨prop_1_4_i P g (fun Y _ => hiso Y), hmet g⟩, by
      show IsIso (P.Base g); rw [hgb]; infer_instance⟩
  have hFn : IsFrobeniusType P (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) :=
    hFs.frobType _ _
  have hFnb : P.Base (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = 𝟙 _ := by
    have hh := hFs.baseIdentity (P.degFr φ) ⟨A, hA⟩
    show P.Base _ = _
    rw [show P.Base (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = P.Base (𝟙 A) from hh,
      P.Base_id]
  have hdegFn : P.degFr (((Fs (P.degFr φ)).app ⟨A, hA⟩ : End A) : A ⟶ A) = P.degFr φ :=
    (SectionEnd.deg_eq (Fs (P.degFr φ)) ⟨A, hA⟩).symm.trans (hFs.degSection _)
  obtain ⟨β, hβiso, hβ⟩ := Fc.frobDegUniq A A A _ g hFn hgF (hdegFn.trans hdeg.symm)
  haveI := hβiso
  have hbβ : P.Base β = 𝟙 _ := by
    have h := congrArg P.Base hβ
    rw [P.Base_comp, hFnb, Category.id_comp, hgb] at h
    exact h
  refine ⟨P.degFr φ, β, S.lift hA hB (P.Base φ), ?_, ?_, ?_, hαP⟩
  · rw [← Category.assoc, hβ, hg1]
  · show P.Base β = P.Base (𝟙 A)
    rw [P.Base_id]; exact hbβ
  · exact ⟨degFr_of_isIso P β, by show IsIso (P.Base β); rw [hbβ]; infer_instance⟩

/-- ★★★**3 分解の `β` は `𝒪^×(A)` の元**(等長な場合)。

★`β` は base-identity な pre-step で、等長なので isotropic 型のもとで**同型**、
したがって `𝒪^×(A)` に入る。★これが `Theorem 5.2, (iv)` の `u` 成分の正体である。 -/
theorem rem_2_7_2_beta_otimes (hiso : ∀ X : C, IsIsotropic P X)
    (hmet : ∀ {A B : C} (f : A ⟶ B), IsIsometric P f)
    {A : C} (β : A ⟶ A) (hb : IsBaseIdentity P β) (hs : IsPreStep P β) :
    (β : End A) ∈ OTimes P A := by
  haveI : IsIso β := hiso A A β (hmet β) hs
  exact ⟨⟨hb, hs.1⟩, (isUnit_iff_isIso (β : End A)).mpr inferInstance⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Remark 2.7.2` の 3 分解(等長版)。 -/
def rem_2_7_2_factor_isometric.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52,
    item := "Remark 2.7.2 — 3 分解の存在(等長な場合、skeleton 不要)",
    sectionId := "frdi-def-2-7" }

end ABC3.Found.FrdI
