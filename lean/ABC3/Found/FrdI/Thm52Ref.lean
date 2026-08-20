/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Path

/-!
# [FrdI] Theorem 5.2, (iv) の第 2 段 —— base-Frobenius pair が定める参照射

原文 (FrdI p.101):
> Let (P, F ) be a base-Frobenius pair of C [cf. our assumption that C is of model type].

★`𝒫` は `𝒟` と圏同値なので、`𝒫` の 2 対象の間の射は**底で一意に決まる**
(`BaseSection.lift`)。★そこへ Frobenius-section `F` の `n` 次の元を掛けると、
**底 `θ`・次数 `n`・零因子 0** の射 `r(θ, n)` が得られる。

★★**この `r` は関手的**である:
```
r(θ ≫ θ', m·n) = r(θ, n) ≫ r(θ', m)
```
—— `F` の自然性で `F_n` を右へ通し、`F` が単系準同型であることを使う。
★これが `𝒞 → 𝒞^model` の「基準」になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `𝒫` の射は底で決まる -/

/-- ★`𝒫` の 2 対象と底の射から定まる `𝒫` の射(充満性)。 -/
noncomputable def BaseSection.lift (S : BaseSection P) {A B : C} (hA : S.objP A)
    (hB : S.objP B) (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) : A ⟶ B :=
  (S.fullP hA hB θ).choose

theorem BaseSection.lift_homP (S : BaseSection P) {A B : C} (hA : S.objP A) (hB : S.objP B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) : S.homP (S.lift hA hB θ) :=
  (S.fullP hA hB θ).choose_spec.1

@[simp] theorem BaseSection.lift_base (S : BaseSection P) {A B : C} (hA : S.objP A)
    (hB : S.objP B) (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) :
    P.Base (S.lift hA hB θ) = θ :=
  (S.fullP hA hB θ).choose_spec.2

/-- ★`𝒫` の射は底で決まる(忠実性)。 -/
theorem BaseSection.eq_lift (S : BaseSection P) {A B : C} (hA : S.objP A) (hB : S.objP B)
    (f : A ⟶ B) (hf : S.homP f) : f = S.lift hA hB (P.Base f) :=
  S.faithfulP hf (S.lift_homP hA hB _) (S.lift_base hA hB _).symm

theorem BaseSection.lift_id (S : BaseSection P) {A : C} (hA : S.objP A) :
    S.lift hA hA (P.Base (𝟙 A)) = 𝟙 A :=
  (S.eq_lift hA hA (𝟙 A) (S.id_mem hA)).symm

theorem BaseSection.lift_comp (S : BaseSection P) {A B E : C} (hA : S.objP A) (hB : S.objP B)
    (hE : S.objP E) (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base)
    (θ' : (P.toElem.obj B).base ⟶ (P.toElem.obj E).base) :
    S.lift hA hE (θ ≫ θ') = S.lift hA hB θ ≫ S.lift hB hE θ' := by
  have hb : P.Base (S.lift hA hB θ ≫ S.lift hB hE θ') = θ ≫ θ' := by
    rw [P.Base_comp, S.lift_base, S.lift_base]
  rw [← hb]
  exact (S.eq_lift hA hE _ (S.comp_mem (S.lift_homP hA hB θ) (S.lift_homP hB hE θ'))).symm

/-! ## ★2. 参照射 -/

/-- ★★**参照射** `r(θ, n) = lift(θ) ≫ F_n` —— 底 `θ`、次数 `n`、零因子 0。 -/
noncomputable def refHom (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S)
    {A B : C} (hA : S.objP A) (hB : S.objP B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) (n : ℕ+) : A ⟶ B :=
  S.lift hA hB θ ≫ (((Fs n).app ⟨B, hB⟩ : End B) : B ⟶ B)

variable [IsConnected D]

theorem refHom_base (S : BaseSection P) {Fs : ℕ+ →* SectionEnd S}
    (hF : IsFrobeniusSection S Fs) {A B : C} (hA : S.objP A) (hB : S.objP B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) (n : ℕ+) :
    P.Base (refHom S Fs hA hB θ n) = θ := by
  rw [refHom, P.Base_comp, S.lift_base,
    show P.Base (((Fs n).app ⟨B, hB⟩ : End B) : B ⟶ B) = P.Base (𝟙 B) from hF.baseIdentity n _,
    P.Base_id, Category.comp_id]

theorem refHom_deg (F : FrobenioidCore P) (S : BaseSection P) {Fs : ℕ+ →* SectionEnd S}
    (hF : IsFrobeniusSection S Fs) {A B : C} (hA : S.objP A) (hB : S.objP B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) (n : ℕ+) :
    P.degFr (refHom S Fs hA hB θ n) = n := by
  have hl : P.degFr (S.lift hA hB θ) = 1 :=
    (F.pullBackLB _ (S.isPullBack (S.lift_homP hA hB θ))).2
  have hn : P.degFr (((Fs n).app ⟨B, hB⟩ : End B) : B ⟶ B) = n :=
    (SectionEnd.deg_eq (Fs n) ⟨B, hB⟩).symm.trans (hF.degSection n)
  rw [refHom, P.degFr_comp, hl, hn, mul_one]

theorem refHom_div (F : FrobenioidCore P) (S : BaseSection P) {Fs : ℕ+ →* SectionEnd S}
    (hF : IsFrobeniusSection S Fs) {A B : C} (hA : S.objP A) (hB : S.objP B)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) (n : ℕ+) :
    P.Div (refHom S Fs hA hB θ n) = 0 := by
  have hl : P.Div (S.lift hA hB θ) = 0 :=
    (F.pullBackLB _ (S.isPullBack (S.lift_homP hA hB θ))).1.2
  have hn : P.Div (((Fs n).app ⟨B, hB⟩ : End B) : B ⟶ B) = 0 := (hF.frobType n ⟨B, hB⟩).1.2
  rw [refHom, P.Div_comp, hl, hn, map_zero, smul_zero, add_zero]

/-- ★★★**参照射は関手的**。

★`F` の自然性で `F_n` を `lift θ'` の右へ通し、`F` が単系準同型であることを使う。 -/
theorem refHom_comp (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S)
    {A B E : C} (hA : S.objP A) (hB : S.objP B) (hE : S.objP E)
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base)
    (θ' : (P.toElem.obj B).base ⟶ (P.toElem.obj E).base) (n m : ℕ+) :
    refHom S Fs hA hE (θ ≫ θ') (m * n)
      = refHom S Fs hA hB θ n ≫ refHom S Fs hB hE θ' m := by
  have hnat : (S.lift hB hE θ') ≫ (((Fs n).app ⟨E, hE⟩ : End E) : E ⟶ E)
      = (((Fs n).app ⟨B, hB⟩ : End B) : B ⟶ B) ≫ (S.lift hB hE θ') :=
    (Fs n).naturality (A := (⟨B, hB⟩ : S.Obj)) (B := (⟨E, hE⟩ : S.Obj))
      ⟨S.lift hB hE θ', S.lift_homP hB hE θ'⟩
  have hmul : (((Fs (m * n)).app ⟨E, hE⟩ : End E) : E ⟶ E)
      = (((Fs n).app ⟨E, hE⟩ : End E) : E ⟶ E) ≫ (((Fs m).app ⟨E, hE⟩ : End E) : E ⟶ E) := by
    rw [map_mul]; rfl
  rw [refHom, refHom, refHom, S.lift_comp hA hB hE, hmul]
  rw [Category.assoc, Category.assoc, ← Category.assoc (S.lift hB hE θ'), hnat,
    Category.assoc]

omit [IsConnected D] in
theorem refHom_id (S : BaseSection P) (Fs : ℕ+ →* SectionEnd S) {A : C} (hA : S.objP A) :
    refHom S Fs hA hA (P.Base (𝟙 A)) 1 = 𝟙 A := by
  rw [refHom, S.lift_id hA, show (((Fs 1).app ⟨A, hA⟩ : End A) : A ⟶ A) = 𝟙 A from by
    rw [map_one]; rfl, Category.comp_id]

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.2, (iv)` の参照射(base-Frobenius pair の使い方)。 -/
def refHom_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — base-Frobenius pair が定める参照射",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
