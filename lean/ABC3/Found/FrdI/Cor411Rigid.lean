/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] Corollary 4.11, (i) の rigidity

原文 (FrdI p.93):
> is a Div-identity automorphism. [Indeed, this follows by applying the functoriality

## ★★★★★測って分かった段取り(2026-08-19)

原文は 3 段で書く:
1. `α ∈ Aut(𝒞ᵢ → 𝒞ᵢ^un-tr)` が誘導する自己同型は **Div-identity**
   (co-angular pre-step への関手性 ＋ `Definition 1.3, (iii), (d)` の第 2 圏同値)
2. `𝒟ᵢ` が **Div-slim** なので **base-identity**
3. `𝒞ᵢ^un-tr` は **unit-trivial** 型なので **恒等**

★★1 の中身を測ると、実は**圏同値そのものは要らず**、
「どの `x ∈ Φ(A)` も `A` から出る co-angular pre-step の零因子である」
(＝圏同値の**本質的全射性**)だけでよい。★あとは `Div_comp` の計算である:

- `Div (α ≫ ϵ) = Φ.map (Base α) (Div ϵ) + (degFr ϵ) • Div α = Φ.map (Base α) x`
- `Div (ϵ ≫ β) = Φ.map (Base ϵ) (Div β) + (degFr β) • Div ϵ = x`

(`α`・`β` は同型なので零因子 0、`ϵ` は pre-step なので次数 1。)
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section DivIdOfAut

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★★**`Definition 1.3, (iii), (d)` の本質的全射性** ——
`Φ(A)` のどの元も `A` から出る co-angular pre-step の零因子である。

★圏同値の像は `Order(Φ(A))` の中で**同型**を与えるだけだが、
`Φ(A)` は integral かつ sharp なので `mle_antisymm` で**等式**に上がる。 -/
theorem exists_coaPreStep_div (G : Frobenioid P) (A : C)
    (x : Φ.val (P.toElem.obj A).base) :
    ∃ (B : C) (ϵ : A ⟶ B), IsCoAngular P ϵ ∧ IsPreStep P ϵ ∧ P.Div ϵ = x := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreUnderEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := coaPreUnderFunctor P A) (toOrderCat x)
  refine ⟨Z.right.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
  have h1 : MLe (P.Div Z.hom.hom) x := leOfHom e.hom
  have h2 : MLe x (P.Div Z.hom.hom) := leOfHom e.inv
  exact mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 h1 h2

/-- ★★★★**自己同型が co-angular pre-step すべてに沿って自然なら Div-identity**。 -/
theorem isDivIdentity_of_preStep_natural (G : Frobenioid P) {A : C} (α : A ⟶ A)
    [IsIso α]
    (hnat : ∀ (B : C) (ϵ : A ⟶ B), IsCoAngular P ϵ → IsPreStep P ϵ →
      ∃ β : B ⟶ B, IsIso β ∧ α ≫ ϵ = ϵ ≫ β) :
    IsDivIdentity P α := by
  show Φ.map (P.Base α) = Φ.map (P.Base (𝟙 A))
  refine AddMonoidHom.ext (fun x => ?_)
  obtain ⟨B, ϵ, hϵc, hϵs, hdiv⟩ := exists_coaPreStep_div G A x
  obtain ⟨β, hβ, hsq⟩ := hnat B ϵ hϵc hϵs
  haveI := hβ
  have hda : P.Div α = 0 := isIsometric_of_isIso P α
  have hdb : P.Div β = 0 := isIsometric_of_isIso P β
  have hL : P.Div (α ≫ ϵ) = Φ.map (P.Base α) x := by
    rw [P.Div_comp, hda, hdiv, smul_zero, add_zero]
  have hR : P.Div (ϵ ≫ β) = x := by
    rw [P.Div_comp, hdb, map_zero, zero_add, degFr_of_isIso P β, hdiv]
    show ((1 : ℕ+) : ℕ) • x = x
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hkey : Φ.map (P.Base α) x = x := by
    rw [← hL, hsq, hR]
  refine hkey.trans ?_
  show x = Φ.map (P.Base (𝟙 A)) x
  rw [P.Base_id, Φ.map_id]

def isDivIdentity_of_preStep_natural.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — 誘導される自己同型は Div-identity",
    sectionId := "frdi-cor-4-11" }

end DivIdOfAut

/-! ## ★2. `𝒞 → 𝒞^un-tr` の rigidity -/

section RigidUnTr

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★`Aut(𝒞 → 𝒞^un-tr)` の自然性は `𝒞^un-tr` の**すべての射**に及ぶ
(`istrToUnTr` が充満だから)。 -/
theorem untrAut_naturality (η : (istrToUnTr P) ≅ (istrToUnTr P))
    {A B : UnTr P} (f : A ⟶ B) :
    f ≫ η.hom.app B = η.hom.app A ≫ f := by
  obtain ⟨f₀, hf₀⟩ := (istrToUnTr P).map_surjective f
  rw [← hf₀]
  exact η.hom.naturality f₀

/-- ★成分は同型。 -/
theorem untrAut_app_isIso (η : (istrToUnTr P) ≅ (istrToUnTr P)) (A : UnTr P) :
    IsIso (η.hom.app A) :=
  ⟨η.inv.app A, η.hom_inv_id_app A, η.inv_hom_id_app A⟩

/-- ★★★**`η` が定めるスライスの自己同型** —— 成分は `Base (η_Z)`。 -/
noncomputable def untrAutSlice (Fc : FrobenioidCore P)
    (η : (istrToUnTr P) ≅ (istrToUnTr P)) (A : UnTr P) :
    (Over.forget A ⋙ (unTrPre P Fc).proj) ≅ (Over.forget A ⋙ (unTrPre P Fc).proj) :=
  NatIso.ofComponents
    (fun Z => @asIso _ _ _ _ ((unTrPre P Fc).Base (η.hom.app Z.left))
      ⟨(unTrPre P Fc).Base (η.inv.app Z.left),
        ((unTrPre P Fc).Base_comp (η.hom.app Z.left) (η.inv.app Z.left)).symm.trans
          ((congrArg (fun t => (unTrPre P Fc).Base t) (η.hom_inv_id_app Z.left)).trans
            ((unTrPre P Fc).Base_id _)),
        ((unTrPre P Fc).Base_comp (η.inv.app Z.left) (η.hom.app Z.left)).symm.trans
          ((congrArg (fun t => (unTrPre P Fc).Base t) (η.inv_hom_id_app Z.left)).trans
            ((unTrPre P Fc).Base_id _))⟩)
    (fun {Z W} g =>
      ((unTrPre P Fc).Base_comp g.left (η.hom.app W.left)).symm.trans
        ((congrArg (fun t => (unTrPre P Fc).Base t) (untrAut_naturality η g.left)).trans
          ((unTrPre P Fc).Base_comp (η.hom.app Z.left) g.left)))

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**成分は恒等** —— 原文 p.93 の 3 段。

原文 (FrdI p.93):
> hence trivial [since Ciun-tr is of unit-trivial type]. This completes the proof of

★★1. 誘導される自己同型は **Div-identity**(`isDivIdentity_of_preStep_natural`)
★★2. `𝒟` が **Div-slim** なので **base-identity**(`divSlim_over_aut_eq_id`)
★★3. `𝒞^un-tr` は **unit-trivial** なので **恒等**(`unTr_unitTrivial`) -/
theorem untrAut_app_eq_id (Fc : FrobenioidCore P)
    (G' : Frobenioid (unTrPre P Fc)) (hds : IsDivSlim Φ)
    (η : (istrToUnTr P) ≅ (istrToUnTr P)) (A : UnTr P) :
    η.hom.app A = 𝟙 A := by
  -- ★段 1
  have hdivid : ∀ B : UnTr P, IsDivIdentity (unTrPre P Fc) (η.hom.app B) := by
    intro B
    haveI := untrAut_app_isIso η B
    refine isDivIdentity_of_preStep_natural G' (η.hom.app B) (fun E ϵ _ _ => ?_)
    exact ⟨η.hom.app E, untrAut_app_isIso η E, (untrAut_naturality η ϵ).symm⟩
  -- ★段 2
  have hdiv : ∀ (Z : Over A) (x : Φ.val ((unTrPre P Fc).proj.obj Z.left)),
      Φ.map ((untrAutSlice Fc η A).hom.app Z) x = x := by
    intro Z x
    have h2 : Φ.map ((unTrPre P Fc).Base (η.hom.app Z.left))
        = Φ.map ((unTrPre P Fc).Base (𝟙 Z.left)) := hdivid Z.left
    have h3 := congrArg (fun t : Φ.val _ →+ Φ.val _ => t x) h2
    refine h3.trans ?_
    show Φ.map ((unTrPre P Fc).Base (𝟙 Z.left)) x = x
    rw [(unTrPre P Fc).Base_id, Φ.map_id]
  have hrefl : untrAutSlice Fc η A = Iso.refl _ :=
    divSlim_over_aut_eq_id (unTrPre P Fc) (unTr_frobenioidCore P Fc) hds A
      (untrAutSlice Fc η A) hdiv
  have hZ₀ := congrArg (fun t : (Over.forget A ⋙ (unTrPre P Fc).proj)
      ≅ (Over.forget A ⋙ (unTrPre P Fc).proj) => t.hom.app (Over.mk (𝟙 A))) hrefl
  -- ★段 3
  have hbi : IsBaseIdentity (unTrPre P Fc) (η.hom.app A) := by
    show (unTrPre P Fc).Base (η.hom.app A) = (unTrPre P Fc).Base (𝟙 A)
    rw [(unTrPre P Fc).Base_id]
    exact hZ₀
  haveI := untrAut_app_isIso η A
  obtain ⟨α, hα⟩ : ∃ α : End A, ((α : A ⟶ A)) = η.hom.app A := ⟨η.hom.app A, rfl⟩
  have hbi' : IsBaseIdentity (unTrPre P Fc) ((α : A ⟶ A)) := by rw [hα]; exact hbi
  have hu' : IsUnit α := by
    haveI : IsIso ((α : A ⟶ A)) := by rw [hα]; exact untrAut_app_isIso η A
    exact (CategoryTheory.isUnit_iff_isIso α).mpr inferInstance
  have hmem : α ∈ OTimes (unTrPre P Fc) A :=
    (mem_otimes_iff (unTrPre P Fc) α).mpr ⟨hu', hbi'⟩
  rw [unTr_unitTrivial P Fc A] at hmem
  have hh : α = 1 := Submonoid.mem_bot.mp hmem
  rw [← hα, hh]
  rfl

/-- ★★★★★★**[FrdI] Corollary 4.11, (i) の rigidity** —— `𝒞 → 𝒞^un-tr` は rigid。 -/
theorem isRigidFunctor_istrToUnTr_of_divSlim (Fc : FrobenioidCore P)
    (G' : Frobenioid (unTrPre P Fc)) (hds : IsDivSlim Φ) :
    IsRigidFunctor (istrToUnTr P) := by
  intro η
  apply Iso.ext
  apply NatTrans.ext
  funext A
  exact untrAut_app_eq_id Fc G' hds η A

def isRigidFunctor_istrToUnTr_of_divSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — 𝒞 → 𝒞^un-tr の rigidity",
    sectionId := "frdi-cor-4-11" }

end RigidUnTr

/-! ## ★3. 合成関手の rigidity -/

section CompRigid

/-- ★★★★**右から充満忠実な関手を合成しても rigid 性は保たれる**。

★在庫の `isRigidFunctor_of_equivalence_comp` は**左から**圏同値を合成する版。
ここは**右から**の版で、`Corollary 4.11` の「合成関手はどれも rigid」に要る。 -/
theorem isRigidFunctor_comp_fullyFaithful {X : Type*} [Category X] {Y : Type*} [Category Y]
    {Z : Type*} [Category Z] (F : X ⥤ Y) (E : Y ⥤ Z) [E.Full] [E.Faithful]
    (h : IsRigidFunctor F) : IsRigidFunctor (F ⋙ E) := by
  intro η
  have hex : ∀ W : X, ∃ t : F.obj W ⟶ F.obj W, E.map t = η.hom.app W :=
    fun W => E.map_surjective _
  choose t ht using hex
  have hnat : ∀ {W V : X} (f : W ⟶ V), F.map f ≫ t V = t W ≫ F.map f := by
    intro W V f
    refine E.map_injective ?_
    rw [E.map_comp, E.map_comp, ht, ht]
    exact η.hom.naturality f
  have hiso : ∀ W : X, IsIso (t W) := by
    intro W
    haveI : IsIso (E.map (t W)) := by rw [ht]; exact ⟨η.inv.app W,
      η.hom_inv_id_app W, η.inv_hom_id_app W⟩
    exact isIso_of_reflects_iso (t W) E
  have hrefl : (NatIso.ofComponents (fun W => @asIso _ _ _ _ (t W) (hiso W))
      (fun {W V} f => hnat f)) = Iso.refl F := h _
  apply Iso.ext
  apply NatTrans.ext
  funext W
  have htW : t W = 𝟙 (F.obj W) :=
    congrArg (fun u : F ≅ F => u.hom.app W) hrefl
  show η.hom.app W = 𝟙 ((F ⋙ E).obj W)
  rw [← ht W, htW]
  exact E.map_id _

def isRigidFunctor_comp_fullyFaithful.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — 合成関手の rigidity",
    sectionId := "frdi-cor-4-11" }

end CompRigid

/-! ## ★4. `Corollary 4.11, (i)` の組み立て -/

section Assemble

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**[FrdI] Corollary 4.11, (i) の rigidity(合成関手)**。

原文 (FrdI p.91):
> where the vertical arrows are the natural projection functors; the horizontal

★`𝒞₁^istr → 𝒞₁^un-tr` が rigid(`isRigidFunctor_istrToUnTr_of_divSlim`)で、
`Ψ^un-tr` は圏同値(充満忠実)なので、合成も rigid。 -/
theorem cor_4_11_i_comp_rigid (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂)
    (Ψ : C ≌ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.functor.map f))
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.functor.map f) → IsPullBack P f)
    (hds₂ : IsDivSlim Φ₂)
    (hDivId : ∀ {X : C} (α : X ⟶ X), IsDivIdentity P α →
      IsDivIdentity P₂ (Ψ.functor.map α))
    (h₁ : IsOfQuasiIsotropicType C P) (h₂ : IsOfQuasiIsotropicType C₂ P₂)
    (hUE' : ∀ {A B : C} (α₁ α₂ : A ⟶ B),
      P₂.toElem.map (Ψ.functor.map α₁) = P₂.toElem.map (Ψ.functor.map α₂) →
        P.toElem.map α₁ = P.toElem.map α₂)
    (G' : Frobenioid (unTrPre P F)) (hds : IsDivSlim Φ) :
    IsRigidFunctor (istrToUnTr P
      ⋙ psiUnTrOfDivSlim P F hiso P₂ F₂ Ψ hPB hPB' hds₂ hDivId h₁ h₂) := by
  haveI := Ψ.isEquivalence_functor
  haveI : (psiUnTrOfDivSlim P F hiso P₂ F₂ Ψ hPB hPB' hds₂ hDivId h₁ h₂).IsEquivalence :=
    psiUnTr_isEquivalence Ψ.functor h₁ h₂
      (fun α₁ α₂ hh => toElem_map_congr_of_otimes Ψ.functor F hiso
        (fun _X δ hδ => otimes_map_of_divSlim P F hiso P₂ F₂ Ψ hPB hPB' hds₂ hDivId δ hδ)
        α₁ α₂ hh) hUE'
  exact isRigidFunctor_comp_fullyFaithful _ _
    (isRigidFunctor_istrToUnTr_of_divSlim F G' hds)

def cor_4_11_i_comp_rigid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.11, (i) — 合成関手はどれも rigid",
    sectionId := "frdi-cor-4-11" }

end Assemble

/-! ## ★5. Div-slim 版の「`𝒞 → 𝔽_Φ` は rigid」 -/

section RigidToElem

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**Div-slim なら `𝒞 → 𝔽_Φ` は rigid**。

原文 (FrdI p.94):
> the rigidity assertion follows via the same argument as was applied to prove the

★★(i) の rigidity と同じ 3 段:
1. 成分は `𝔽_Φ` の自己同型なので **`deg = 1`・`div = 0`**(`isIso_iff` / sharp)
2. co-angular pre-step への自然性で **底の自己同型は `Φ` を恒等へ送る**
3. **Div-slim** で底が恒等 ⟹ `⟨𝟙, 0, 1⟩ = 𝟙` -/
theorem isRigidFunctor_toElem_of_divSlim (F : FrobenioidCore P) (G : Frobenioid P)
    (hds : IsDivSlim Φ) : IsRigidFunctor P.toElem := by
  intro η
  have hI : ∀ A : C, IsIso (η.hom.app A) := fun A =>
    ⟨η.inv.app A, η.hom_inv_id_app A, η.inv_hom_id_app A⟩
  have hdg : ∀ A : C, (η.hom.app A).deg = 1 := fun A =>
    ((ElemFrobCat.isIso_iff _).mp (hI A)).2.2
  have hdv : ∀ A : C, (η.hom.app A).div = 0 := by
    intro A
    haveI := hI A
    exact ElemFrobCat.div_eq_zero_of_isIso (fun X => (P.divisorial X).2) _
  have hbI : ∀ A : C, IsIso ((η.hom.app A).base) := fun A =>
    ((ElemFrobCat.isIso_iff _).mp (hI A)).1
  -- ★段 2: 底の自己同型は `Φ` を恒等へ送る
  have hDivId : ∀ (A : C) (x : Φ.val ((P.toElem.obj A).base)),
      Φ.map ((η.hom.app A).base) x = x := by
    intro A x
    obtain ⟨B, ϵ, hϵc, hϵs, hdivϵ⟩ := exists_coaPreStep_div G A x
    have h := congrArg ElemFrobCat.Hom.div (η.hom.naturality ϵ)
    rw [ElemFrobCat.comp_div, ElemFrobCat.comp_div, hdv B, hdv A, hdg B,
      map_zero, zero_add, smul_zero, add_zero] at h
    have h2 : ((1 : ℕ+) : ℕ) • P.Div ϵ = Φ.map ((η.hom.app A).base) (P.Div ϵ) := h
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul] at h2
    rw [← hdivϵ]
    exact h2.symm
  -- ★段 3: Div-slim で底が恒等
  have hbase : ∀ A : C, (η.hom.app A).base = 𝟙 ((P.toElem.obj A).base) := by
    intro A
    set ζ : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj) :=
      NatIso.ofComponents
        (fun Z => @asIso _ _ _ _ ((η.hom.app Z.left).base) (hbI Z.left))
        (fun {Z W} g => congrArg ElemFrobCat.Hom.base (η.hom.naturality g.left)) with hζ
    have hdiv2 : ∀ (Z : Over A) (x : Φ.val (P.proj.obj Z.left)),
        Φ.map (ζ.hom.app Z) x = x := fun Z x => hDivId Z.left x
    have hrefl : ζ = Iso.refl _ := divSlim_over_aut_eq_id P F hds A ζ hdiv2
    exact congrArg (fun t : (Over.forget A ⋙ P.proj) ≅ (Over.forget A ⋙ P.proj) =>
      t.hom.app (Over.mk (𝟙 A))) hrefl
  apply Iso.ext
  apply NatTrans.ext
  funext A
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · rw [hbase A]; rfl
  · rw [hdv A]; rfl
  · rw [hdg A]; rfl

def isRigidFunctor_toElem_of_divSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (iv) — Div-slim なら 𝒞 → 𝔽_Φ は rigid",
    sectionId := "frdi-cor-4-11" }

end RigidToElem

end ABC3.Found.FrdI
