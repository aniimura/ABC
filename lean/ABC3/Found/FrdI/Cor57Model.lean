/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor57Base
import ABC3.Found.FrdI.Def45

/-!
# [FrdI] Corollary 5.7, (i) の後半 —— base-Frobenius 対と model 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.108。

原文 (FrdI p.108):
> to base-sections (respectively, quasi-base-Frobenius pairs) of C2. In particular, C1 is

## ★★何が残っていたか

`Cor57Base.lean` は `BaseSection P₁ → BaseSection P₂` を作った。
残るのは **`𝒫`-Frobenius-section `F` の輸送**である ——
これが済むと `BaseFrobeniusPair` が移り、
原文の「In particular, `C₁` is of model type if and only if `C₂` is」が出る。

## ★★原像は一意である

輸送した `objP₂` は `∃ A₁, objP₁ A₁ ∧ Ψ A₁ = A₂` なので、
原像 `A₁` を `Classical.choose` で取る。★**一意性は `objP_eq_of_map_eq`** が与えるので、
どの証人を取っても同じである(`preObj_eq`)。

## ★共役は 3 本の補題で消す

`app X := eqToHom ≫ Ψ.map (ε.app A₁) ≫ eqToHom` の形なので、
`eqToHom` による共役が base-identity・Frobenius 型・次数を変えないことを
`subst` 1 行ずつで用意しておく。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section FrobNormIso

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★同型による共役。 -/
def conjIso {X Y : C} (e : X ≅ Y) (α : End Y) : End X :=
  (e.hom ≫ (α : Y ⟶ Y) ≫ e.inv : X ⟶ X)

@[simp] theorem conjIso_hom {X Y : C} (e : X ≅ Y) (α : End Y) :
    ((conjIso e α : End X) : X ⟶ X) = e.hom ≫ (α : Y ⟶ Y) ≫ e.inv := rfl

/-- ★共役は冪を保つ。 -/
theorem conjIso_pow {X Y : C} (e : X ≅ Y) (α : End Y) (n : ℕ) :
    (conjIso e α) ^ n = conjIso e (α ^ n) := by
  induction n with
  | zero =>
      show ((𝟙 X : X ⟶ X) : End X) = conjIso e 1
      rw [conjIso]
      show (𝟙 X : X ⟶ X) = e.hom ≫ (𝟙 Y : Y ⟶ Y) ≫ e.inv
      simp
  | succ n ih =>
      show ((conjIso e α : End X) : X ⟶ X) ≫ (((conjIso e α) ^ n : End X) : X ⟶ X)
        = ((conjIso e (α ^ (n + 1)) : End X) : X ⟶ X)
      rw [ih]
      show (e.hom ≫ (α : Y ⟶ Y) ≫ e.inv) ≫ (e.hom ≫ ((α ^ n : End Y) : Y ⟶ Y) ≫ e.inv)
        = e.hom ≫ (((α : Y ⟶ Y) ≫ ((α ^ n : End Y) : Y ⟶ Y))) ≫ e.inv
      simp

/-- ★共役は base-identity を保つ。 -/
theorem isBaseIdentity_conjIso {X Y : C} (e : X ≅ Y) (φ : End Y)
    (hφ : IsBaseIdentity P φ) : IsBaseIdentity P (conjIso e φ) := by
  have hb : P.Base ((φ : Y ⟶ Y)) = 𝟙 _ := by
    rw [show P.Base ((φ : Y ⟶ Y)) = P.Base (𝟙 Y) from hφ, P.Base_id]
  show P.Base (e.hom ≫ (φ : Y ⟶ Y) ≫ e.inv) = P.Base (𝟙 X)
  rw [P.Base_comp, P.Base_comp, hb, Category.id_comp, ← P.Base_comp, e.hom_inv_id]

/-- ★共役は次数を保つ。 -/
theorem degFr_conjIso {X Y : C} (e : X ≅ Y) (φ : End Y) :
    P.degFr ((conjIso e φ : End X) : X ⟶ X) = P.degFr ((φ : Y ⟶ Y)) := by
  show P.degFr (e.hom ≫ (φ : Y ⟶ Y) ≫ e.inv) = _
  rw [P.degFr_comp, P.degFr_comp, isLinear_of_isIso P e.hom, isLinear_of_isIso P e.inv]
  simp

/-- ★共役は `𝒪^▷` を保つ。 -/
theorem oTri_conjIso {X Y : C} (e : X ≅ Y) (α : End Y) (hα : α ∈ OTri P Y) :
    conjIso e α ∈ OTri P X :=
  ⟨isBaseIdentity_conjIso e α hα.1, by
    show P.degFr ((conjIso e α : End X) : X ⟶ X) = 1
    rw [degFr_conjIso]
    exact hα.2⟩

/-- ★★★**Frobenius-normalized は対象の同型で保たれる**。 -/
theorem isFrobeniusNormalized_of_iso {X Y : C} (e : X ≅ Y)
    (h : IsFrobeniusNormalized P X) : IsFrobeniusNormalized P Y := by
  intro φ hφ α hα
  have key := h (conjIso e φ) (isBaseIdentity_conjIso e φ hφ)
    (conjIso e α) (oTri_conjIso e α hα)
  rw [degFr_conjIso, conjIso_pow] at key
  have h2 := congrArg (fun t : X ⟶ X => e.inv ≫ t ≫ e.hom) key
  simp only [conjIso_hom] at h2
  simpa using h2

end FrobNormIso

section Cor57Model

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  {Ψ : C₁ ⥤ C₂} {ΨBase : D₁ ⥤ D₂}

/-! ## ★1. `eqToHom` による共役 -/

theorem isBaseIdentity_eqToHom_conj {A B : C₂} (h : A = B) (φ : A ⟶ A)
    (hφ : IsBaseIdentity P₂ φ) :
    IsBaseIdentity P₂ ((eqToHom h.symm ≫ φ ≫ eqToHom h : B ⟶ B)) := by
  subst h
  simpa using hφ

theorem isFrobeniusType_eqToHom_conj {A B : C₂} (h : A = B) (φ : A ⟶ A)
    (hφ : IsFrobeniusType P₂ φ) :
    IsFrobeniusType P₂ ((eqToHom h.symm ≫ φ ≫ eqToHom h : B ⟶ B)) := by
  subst h
  simpa using hφ

theorem degFr_eqToHom_conj {A B : C₂} (h : A = B) (φ : A ⟶ A) :
    P₂.degFr ((eqToHom h.symm ≫ φ ≫ eqToHom h : B ⟶ B)) = P₂.degFr φ := by
  subst h
  simp

/-! ## ★2. 原像 -/

/-- ★輸送した対象の**原像**（`Classical.choose`）。 -/
noncomputable def preObj (S : BaseSection P₁) (F : C₁ ⥤ C₂) {A₂ : C₂}
    (h : objPmap S F A₂) : C₁ := h.choose

theorem preObj_objP (S : BaseSection P₁) (F : C₁ ⥤ C₂) {A₂ : C₂} (h : objPmap S F A₂) :
    S.objP (preObj S F h) := h.choose_spec.1

theorem map_preObj (S : BaseSection P₁) (F : C₁ ⥤ C₂) {A₂ : C₂} (h : objPmap S F A₂) :
    F.obj (preObj S F h) = A₂ := h.choose_spec.2

/-- ★★**原像は一意** —— `objP_eq_of_map_eq` による。 -/
theorem preObj_eq {S : BaseSection P₁} [Ψ.Full] [Ψ.Faithful] {A₂ : C₂}
    (h : objPmap S Ψ A₂) {A₁ : C₁} (hA₁ : S.objP A₁) (e : Ψ.obj A₁ = A₂) :
    preObj S Ψ h = A₁ :=
  objP_eq_of_map_eq S (preObj_objP S Ψ h) hA₁ ((map_preObj S Ψ h).trans e.symm)

/-- ★`SectionEnd` の成分は対象の等式で移る。 -/
theorem sectionEnd_app_congr {S : BaseSection P₁} (ε : SectionEnd S) {A B : C₁} (h : A = B)
    (hA : S.objP A) (hB : S.objP B) :
    ((ε.app ⟨B, hB⟩ : End B) : B ⟶ B)
      = eqToHom h.symm ≫ ((ε.app ⟨A, hA⟩ : End A) : A ⟶ A) ≫ eqToHom h := by
  subst h
  simp

/-! ## ★3. `SectionEnd` の輸送 -/

section Trans

variable {S : BaseSection P₁} [Ψ.Full] [Ψ.Faithful]
  [ΨBase.Full] [ΨBase.Faithful] [ΨBase.EssSurj]
  (bs : BaseSquare Ψ ΨBase P₁ P₂)
  (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
  (hft : ∀ {A : C₁}, S.objP A → IsFrobeniusTrivial P₂ (Ψ.obj A))

/-- ★輸送した成分の値。 -/
noncomputable def transApp (ε : SectionEnd S) (X : (S.map bs hpb hft).Obj) : End X.1 :=
  (eqToHom (map_preObj S Ψ X.2).symm
    ≫ Ψ.map ((ε.app ⟨preObj S Ψ X.2, preObj_objP S Ψ X.2⟩ : End _) : _ ⟶ _)
    ≫ eqToHom (map_preObj S Ψ X.2) : X.1 ⟶ X.1)

/-- ★★**輸送した成分は証人の取り方に依らない**。 -/
theorem transApp_eq (ε : SectionEnd S) (X : (S.map bs hpb hft).Obj)
    {A₁ : C₁} (hA₁ : S.objP A₁) (hA : Ψ.obj A₁ = X.1) :
    ((transApp bs hpb hft ε X : End X.1) : X.1 ⟶ X.1)
      = eqToHom hA.symm ≫ Ψ.map ((ε.app ⟨A₁, hA₁⟩ : End A₁) : A₁ ⟶ A₁) ≫ eqToHom hA := by
  have h : preObj S Ψ X.2 = A₁ := preObj_eq X.2 hA₁ hA
  rw [sectionEnd_app_congr ε h (preObj_objP S Ψ X.2) hA₁]
  show eqToHom (map_preObj S Ψ X.2).symm
      ≫ Ψ.map ((ε.app ⟨preObj S Ψ X.2, preObj_objP S Ψ X.2⟩ : End _) : _ ⟶ _)
      ≫ eqToHom (map_preObj S Ψ X.2)
    = eqToHom hA.symm
      ≫ Ψ.map (eqToHom h.symm
          ≫ ((ε.app ⟨preObj S Ψ X.2, preObj_objP S Ψ X.2⟩ : End _) : _ ⟶ _)
          ≫ eqToHom h)
      ≫ eqToHom hA
  simp [eqToHom_map]

/-- ★★★★**`SectionEnd` の輸送**。 -/
noncomputable def SectionEnd.mapTo (ε : SectionEnd S) : SectionEnd (S.map bs hpb hft) where
  app X := transApp bs hpb hft ε X
  naturality := fun {X Y} f => by
    obtain ⟨A₁, B₁, f₁, hA₁, hB₁, hf₁, hA, hB, hfe⟩ := f.2
    have hX := transApp_eq bs hpb hft ε X hA₁ hA
    have hY := transApp_eq bs hpb hft ε Y hB₁ hB
    show (f.1 : X.1 ⟶ Y.1) ≫ ((transApp bs hpb hft ε Y : End Y.1) : Y.1 ⟶ Y.1)
      = ((transApp bs hpb hft ε X : End X.1) : X.1 ⟶ X.1) ≫ (f.1 : X.1 ⟶ Y.1)
    rw [hX, hY]
    let XA : S.Obj := ⟨A₁, hA₁⟩
    let YB : S.Obj := ⟨B₁, hB₁⟩
    have hnat0 : (f₁ : A₁ ⟶ B₁) ≫ ((ε.app YB : End B₁) : B₁ ⟶ B₁)
        = ((ε.app XA : End A₁) : A₁ ⟶ A₁) ≫ (f₁ : A₁ ⟶ B₁) :=
      ε.naturality (⟨f₁, hf₁⟩ : XA ⟶ YB)
    have hnat := congrArg Ψ.map hnat0
    rw [Ψ.map_comp, Ψ.map_comp, hfe] at hnat
    have h3 := congrArg (fun t : Ψ.obj A₁ ⟶ Ψ.obj B₁ => eqToHom hA.symm ≫ t ≫ eqToHom hB) hnat
    simpa using h3

theorem transApp_def (ε : SectionEnd S) (X : (S.map bs hpb hft).Obj) :
    ((transApp bs hpb hft ε X : End X.1) : X.1 ⟶ X.1)
      = eqToHom (map_preObj S Ψ X.2).symm
        ≫ Ψ.map ((ε.app ⟨preObj S Ψ X.2, preObj_objP S Ψ X.2⟩ : End _) : _ ⟶ _)
        ≫ eqToHom (map_preObj S Ψ X.2) := rfl

theorem transApp_one (X : (S.map bs hpb hft).Obj) :
    (transApp bs hpb hft (1 : SectionEnd S) X : End X.1) = 1 := by
  refine (?_ : ((transApp bs hpb hft (1 : SectionEnd S) X : End X.1) : X.1 ⟶ X.1)
    = ((1 : End X.1) : X.1 ⟶ X.1))
  rw [transApp_def]
  show eqToHom (map_preObj S Ψ X.2).symm ≫ Ψ.map (𝟙 _) ≫ eqToHom (map_preObj S Ψ X.2)
    = 𝟙 X.1
  simp

theorem transApp_mul (ε ε' : SectionEnd S) (X : (S.map bs hpb hft).Obj) :
    (transApp bs hpb hft (ε * ε') X : End X.1)
      = transApp bs hpb hft ε X * transApp bs hpb hft ε' X := by
  refine (?_ : ((transApp bs hpb hft (ε * ε') X : End X.1) : X.1 ⟶ X.1)
    = ((transApp bs hpb hft ε X * transApp bs hpb hft ε' X : End X.1) : X.1 ⟶ X.1))
  rw [transApp_def, transApp_def, transApp_def]
  simp [SectionEnd.mul_app]

/-- ★★★★**`𝒫`-Frobenius-section の輸送**。 -/
noncomputable def frobMapTo (Fs : ℕ+ →* SectionEnd S) :
    ℕ+ →* SectionEnd (S.map bs hpb hft) where
  toFun n := SectionEnd.mapTo bs hpb hft (Fs n)
  map_one' := by
    refine SectionEnd.ext ?_
    funext X
    show transApp bs hpb hft (Fs 1) X = _
    rw [map_one]
    exact transApp_one bs hpb hft X
  map_mul' m n := by
    refine SectionEnd.ext ?_
    funext X
    show transApp bs hpb hft (Fs (m * n)) X = _
    rw [map_mul]
    exact transApp_mul bs hpb hft (Fs m) (Fs n) X


/-! ## ★4. Frobenius-section の性質 -/

/-- ★★★★★**輸送した Frobenius-section は Frobenius-section である**。 -/
theorem isFrobeniusSection_frobMapTo [IsConnected D₁] [IsConnected D₂] [Ψ.IsEquivalence]
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂)
    (Fs : ℕ+ →* SectionEnd S) (hFs : IsFrobeniusSection S Fs) :
    IsFrobeniusSection (S.map bs hpb hft) (frobMapTo bs hpb hft Fs) where
  degSection n := by
    haveI := (S.map bs hpb hft).isConnected_obj
    have hX : ∀ X : (S.map bs hpb hft).Obj,
        (frobMapTo bs hpb hft Fs n).deg = n := by
      intro X
      rw [SectionEnd.deg_eq _ X]
      show P₂.degFr ((transApp bs hpb hft (Fs n) X : End X.1) : X.1 ⟶ X.1) = n
      rw [transApp_def, degFr_eqToHom_conj, hdeg]
      exact (SectionEnd.deg_eq (Fs n) _).symm.trans (hFs.degSection n)
    exact hX (Classical.arbitrary _)
  baseIdentity n X := by
    show IsBaseIdentity P₂ ((transApp bs hpb hft (Fs n) X : End X.1) : X.1 ⟶ X.1)
    rw [transApp_def]
    exact isBaseIdentity_eqToHom_conj _ _ (bs.isBaseIdentity_map _ (hFs.baseIdentity n _))
  frobType n X := by
    show IsFrobeniusType P₂ ((transApp bs hpb hft (Fs n) X : End X.1) : X.1 ⟶ X.1)
    rw [transApp_def]
    exact isFrobeniusType_eqToHom_conj _ _
      (isFrobeniusType_map Ψ P₁ P₂ F₁ F₂ hdeg (hFs.frobType n _))

end Trans

/-! ## ★5. base-Frobenius 対と pre-model 型 -/

section Model

variable [Ψ.IsEquivalence] [ΨBase.Full] [ΨBase.Faithful] [ΨBase.EssSurj]
  [IsConnected D₁] [IsConnected D₂]

/-- ★★★★★★**[FrdI] Corollary 5.7, (i) の後半** —— base-Frobenius 対が移る。 -/
noncomputable def BaseFrobeniusPair.mapTo (p : BaseFrobeniusPair P₁)
    (bs : BaseSquare Ψ ΨBase P₁ P₂)
    (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
    (hft : ∀ {A : C₁}, IsFrobeniusTrivial P₁ A → IsFrobeniusTrivial P₂ (Ψ.obj A))
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂) :
    BaseFrobeniusPair P₂ where
  sec := p.sec.map bs hpb (fun h => hft (p.sec.frobTrivial h))
  frob := frobMapTo bs hpb (fun h => hft (p.sec.frobTrivial h)) p.frob
  isFrobSection :=
    isFrobeniusSection_frobMapTo bs hpb (fun h => hft (p.sec.frobTrivial h))
      F₁ F₂ hdeg p.frob p.isFrobSection

/-- ★★★★★**[FrdI] Corollary 5.7, (i)** —— **pre-model 型は圏同値で移る**。

★★原文の「In particular, `C₁` is of model type if and only if `C₂` is」のうち、
`Definition 4.5, (i)` の第 1 条（pre-model 型）の分である。
★第 2 条（birationally Frobenius-normalized 型）の輸送は
birationalization の関手 `Ψ^birat` を経由するので別の段である。

原文 (FrdI p.108):
> to base-sections (respectively, quasi-base-Frobenius pairs) of C2. In particular, C1 is -/
theorem isPreModelType_map
    (bs : BaseSquare Ψ ΨBase P₁ P₂)
    (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
    (hft : ∀ {A : C₁}, IsFrobeniusTrivial P₁ A → IsFrobeniusTrivial P₂ (Ψ.obj A))
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂)
    (h : IsPreModelType P₁) : IsPreModelType P₂ := by
  obtain ⟨p⟩ := h
  exact ⟨p.mapTo bs hpb hft F₁ F₂ hdeg⟩

/-- ★★★★★★**[FrdI] Corollary 5.7, (i)** —— **model 型は圏同値で移る**。

★第 1 条（pre-model 型）は `isPreModelType_map`、
第 2 条（birationally Frobenius-normalized 型）は仮定 `hbfn`（`Corollary 4.11, (ii)` が与える）と
`isFrobeniusNormalized_of_iso`（対象の同型で保たれる）ですべての対象に広げる。

原文 (FrdI p.108):
> to base-sections (respectively, quasi-base-Frobenius pairs) of C2. In particular, C1 is -/
theorem isOfModelType_map [Ψ.EssSurj] (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
    (bs : BaseSquare Ψ ΨBase P₁ P₂)
    (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
    (hft : ∀ {A : C₁}, IsFrobeniusTrivial P₁ A → IsFrobeniusTrivial P₂ (Ψ.obj A))
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hdeg : PreservesDegFr Ψ P₁ P₂)
    (hbfn : ∀ A : C₁, IsBirationallyFrobeniusNormalized P₁ G₁ A →
      IsBirationallyFrobeniusNormalized P₂ G₂ (Ψ.obj A))
    (h : IsOfModelType C₁ P₁ G₁) : IsOfModelType C₂ P₂ G₂ := by
  refine ⟨isPreModelType_map bs hpb hft F₁ F₂ hdeg h.1, ?_⟩
  intro A₂
  obtain ⟨A₁, ⟨e⟩⟩ : ∃ A₁ : C₁, Nonempty (Ψ.obj A₁ ≅ A₂) :=
    ⟨Ψ.objPreimage A₂, ⟨Ψ.objObjPreimageIso A₂⟩⟩
  exact isFrobeniusNormalized_of_iso ((toBiratCat P₂ G₂).mapIso e) (hbfn A₁ (h.2 A₁))

/-! ### ★出典の紐付け -/

/-- ★★★★★★**[FrdI] Corollary 5.7**（条なしの locator）。

| 条 | 実装 |
|---|---|
| (i) base-section を移す | `BaseSection.map`（`Cor57Base.lean`） |
| (i) base-Frobenius 対を移す | `BaseFrobeniusPair.mapTo` |
| (i) model 型は移る | `isOfModelType_map` |
| (ii) unit-profinite 型は移る | `isOfUnitProfiniteType_map`（`Cor57Unit.lean`） |
| (iii)(iv) quasi-base-Frobenius 対を移す | `BaseFrobPair.map`（`Cor57Pair.lean`） |

★★**仮定は原文が挙げるものそのまま**である —— 原文は
「定義を辿れば、`Ψ` が isotropic 対象・prime-Frobenius 射・pull-back 射・
birationalization・射影関手（したがって単元）を保つことを示せば足りる。
それは `Theorem 3.4, (i), (iii)`; `Corollary 4.10`; `Corollary 4.11, (ii)` から出る」
と書くので、底の 1-可換図式（`BaseSquare`）・pull-back 射の保存・
Frobenius-trivial 性の保存・次数の保存・birationally Frobenius-normalized 性の保存を
仮定として受け取る。

★逸脱の記録: 底の四角形を**素の型**で持つ `BaseSquare` を置いた
（`Cor57Pair.lean` の冒頭）。合成関手版からは `BaseSquare.ofNatIso` で作れる。 -/
def cor_5_7.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 108,
    item := "Corollary 5.7",
    sectionId := "frdi-cor-5-7" }

/-- ★★★★locator —— `Corollary 5.7, (i)` の後半（base-Frobenius 対と pre-model 型）。 -/
def isPreModelType_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 108,
    item := "Corollary 5.7, (i) — base-Frobenius 対と pre-model 型が移る",
    sectionId := "frdi-cor-5-7" }

end Model

end Cor57Model

end ABC3.Found.FrdI
