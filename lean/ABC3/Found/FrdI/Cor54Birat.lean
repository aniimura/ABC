/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54
import ABC3.Found.FrdI.Cor411Birat
import ABC3.Found.FrdI.Cor411Phi

/-!
# [FrdI] Corollary 5.4 の `hbirat` —— `η` は `Φ^birat` を `Φ^birat` へ移す

原文 (FrdI p.104):
> 5.4 follows immediately from Corollaries 4.10; 4.11, (iii), (iv). [Here, we note that

★★`Cor54.lean` は `Ψ^rlf` を作るとき、因子の単系の同型
`η : Φ₁ ≅ Ψ_𝒟^* Φ₂` が **`Φ^birat` を `Φ^birat` へ移すこと**(`hbirat`)を
仮定として受け取っていた。本ファイルはそれを `Corollary 4.10`/`4.11` から**導く**。

## ★筋

1. `y ∈ Φ^birat(d)` を `biratDivGp δ`(`δ ∈ 𝒪^×(A^birat)`)の形にする
   (`mem_phiBiratOn_iff`)。
2. `Ψ^birat δ ∈ 𝒪^×(Ψ(A)^birat)`(`biratPsi_otimes_map`、`Corollary 4.11, (ii)`)。
3. ★**残っていた 1 段**は `η ∘ Div^gp = Div^gp ∘ Ψ^birat`。これを
   `sliceDivGpOf`(`Div^gp` の代表元の形)の層で証明する(`sliceDivGpOf_transport`)。
   材料は 3 つだけ:
   * `η` の自然性(`gpMap_phiIsoApp_nat`)
   * `η ∘ Div = Div ∘ Ψ`(`cor_4_11_iv_hdivc`、`Corollary 4.11, (iv)`)
   * `sq` の自然性(底の 1-可換図式)
4. 代表元は `HomBirat.exists_rep` で取り、`Ψ^birat` の代表元での値は
   `biratPsiMap_mk`(`Corollary 4.10`)が与える。

★★したがって `hbirat` は原文どおり **`Corollary 4.10`; `4.11, (ii)(iii)(iv)`** から出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)

/-! ## ★1. `η` の自然性を `Gp` の層で -/

/-- ★`η` の自然性を `Gp` の層で。 -/
theorem gpMap_phiIsoApp_nat (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    {A B : D₁} (f : B ⟶ A) (x : Gp (Φ₁.val A)) :
    gpMap _ (phiIsoApp ΨB η B) (gpMap _ (Φ₁.map f) x)
      = gpMap _ (Φ₂.map (ΨB.map f)) (gpMap _ (phiIsoApp ΨB η A) x) := by
  have hc : (phiIsoApp ΨB η B).comp (Φ₁.map f)
      = (Φ₂.map (ΨB.map f)).comp (phiIsoApp ΨB η A) := by
    ext y
    exact phiIsoApp_nat ΨB η f y
  rw [← AddMonoidHom.comp_apply, ← gpMap_comp, ← AddMonoidHom.comp_apply, ← gpMap_comp, hc]

/-! ## ★2. `Div^gp` の代表元の形が移ること

★★ここが `cor54-hbirat` の**残っていた 1 段**である。 -/

/-- ★★★★**`Div^gp`(代表元の形)は `Ψ` と `η` で移る** —— `Corollary 4.11, (iii)(iv)` の合流点。

★`sliceDivGpOf a ha φ = (Base a)⁻¹ ·(Div φ - deg φ · Div a)` なので、
必要なのは (1) `η` の自然性、(2) `η ∘ Div = Div ∘ Ψ`、(3) `sq` の自然性の 3 つだけ。 -/
theorem sliceDivGpOf_transport (Ψ : C₁ ⥤ C₂) (ΨB : D₁ ⥤ D₂)
    (sqm : ∀ X : C₁, ΨB.obj ((P₁.toElem.obj X).base) ⟶ (P₂.toElem.obj (Ψ.obj X)).base)
    (hsq : ∀ {X Y : C₁} (f : X ⟶ Y), ΨB.map (P₁.Base f) ≫ sqm Y = sqm X ≫ P₂.Base (Ψ.map f))
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (hdiv : ∀ {X Y : C₁} (f : X ⟶ Y),
      phiIsoApp ΨB η ((P₁.toElem.obj X).base) (P₁.Div f)
        = Φ₂.map (sqm X) (P₂.Div (Ψ.map f)))
    (hdeg : ∀ {X Y : C₁} (f : X ⟶ Y), P₂.degFr (Ψ.map f) = P₁.degFr f)
    {A B A' : C₁} (a : A' ⟶ A) (ha : IsIso (P₁.Base a))
    (ha₂ : IsIso (P₂.Base (Ψ.map a))) (φ : A' ⟶ B) :
    gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj A).base)) (sliceDivGpOf (P := P₁) a ha φ)
      = gpMap _ (Φ₂.map (sqm A))
          (sliceDivGpOf (P := P₂) (Ψ.map a) ha₂ (Ψ.map φ)) := by
  haveI := ha
  haveI := ha₂
  set Y : Gp (Φ₂.val ((P₂.toElem.obj (Ψ.obj A')).base)) :=
    toGp _ (P₂.Div (Ψ.map φ)) - ((P₁.degFr φ : ℕ+) : ℕ) • toGp _ (P₂.Div (Ψ.map a)) with hY
  have hmor : ΨB.map (inv (P₁.Base a)) ≫ sqm A' = sqm A ≫ inv (P₂.Base (Ψ.map a)) := by
    have h0 : ΨB.map (inv (P₁.Base a)) ≫ sqm A' ≫ P₂.Base (Ψ.map a) = sqm A := by
      rw [← hsq a, ← Category.assoc, ← ΨB.map_comp, IsIso.inv_hom_id, ΨB.map_id,
        Category.id_comp]
    have h1 := congrArg (fun t => t ≫ inv (P₂.Base (Ψ.map a))) h0
    refine Eq.trans ?_ h1
    show ΨB.map (inv (P₁.Base a)) ≫ sqm A'
      = (ΨB.map (inv (P₁.Base a)) ≫ sqm A' ≫ P₂.Base (Ψ.map a)) ≫ inv (P₂.Base (Ψ.map a))
    rw [Category.assoc, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  have hstep1 : gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj A').base))
      (toGp _ (P₁.Div φ) - ((P₁.degFr φ : ℕ+) : ℕ) • toGp _ (P₁.Div a))
      = gpMap _ (Φ₂.map (sqm A')) Y := by
    rw [hY, map_sub, map_sub, map_nsmul, map_nsmul, gpMap_toGp, gpMap_toGp,
      gpMap_toGp, gpMap_toGp, hdiv φ, hdiv a]
  have hL : gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj A).base)) (sliceDivGpOf (P := P₁) a ha φ)
      = gpMap _ (Φ₂.map (ΨB.map (inv (P₁.Base a)) ≫ sqm A')) Y := by
    rw [sliceDivGpOf_eq, gpMap_phiIsoApp_nat, hstep1, gpMap_phi_comp]
  have hR : gpMap _ (Φ₂.map (sqm A))
      (sliceDivGpOf (P := P₂) (Ψ.map a) ha₂ (Ψ.map φ))
      = gpMap _ (Φ₂.map (sqm A ≫ inv (P₂.Base (Ψ.map a)))) Y := by
    rw [sliceDivGpOf_eq, gpMap_phi_comp, hY, hdeg φ]
  rw [hL, hR, hmor]

variable (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))

/-- ★★★★★**`Div^gp` は `Ψ^birat` と両立する** —— `η ∘ Div^gp = Div^gp ∘ Ψ^birat`。 -/
theorem biratDivGp_transport (ΨB : D₁ ⥤ D₂)
    (sqm : ∀ X : C₁, ΨB.obj ((P₁.toElem.obj X).base) ⟶ (P₂.toElem.obj (Ψ.obj X)).base)
    (hsq : ∀ {X Y : C₁} (f : X ⟶ Y), ΨB.map (P₁.Base f) ≫ sqm Y = sqm X ≫ P₂.Base (Ψ.map f))
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (hdiv : ∀ {X Y : C₁} (f : X ⟶ Y),
      phiIsoApp ΨB η ((P₁.toElem.obj X).base) (P₁.Div f)
        = Φ₂.map (sqm X) (P₂.Div (Ψ.map f)))
    (hdeg : ∀ {X Y : C₁} (f : X ⟶ Y), P₂.degFr (Ψ.map f) = P₁.degFr f)
    {A B : C₁} (f : HomBirat P₁ G₁ A B) :
    gpMap _ (phiIsoApp ΨB η ((P₁.toElem.obj A).base)) (biratDivGp (P := P₁) (G := G₁) f)
      = gpMap _ (Φ₂.map (sqm A))
          (biratDivGp (P := P₂) (G := G₂) (biratPsiMap G₁ G₂ Ψ hfwd A B f)) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  rw [biratDivGp_mk, biratPsiMap_mk, biratDivGp_mk]
  exact sliceDivGpOf_transport Ψ ΨB sqm hsq η hdiv hdeg _ _ _ φ

/-! ## ★3. `Φ^birat` が移ること -/

/-- ★★★★★**`η` は `Φ^birat` を `Φ^birat` へ移す** —— `Corollary 5.4` の `hbirat`。

★`d : 𝒟₁` に対して底が `d` と同型な `A`(`biratBaseObj`)を取り、
`Φ^birat(d)` の元を `𝒪^×(A^birat)` の元の `Div^gp` として書く。
行き先の同型は `(sq A)⁻¹ ≫ Ψ_𝒟(e)` である。 -/
theorem phiBiratOn_transport [Ψ.IsEquivalence]
    (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (ΨB : D₁ ⥤ D₂)
    (sqm : ∀ X : C₁, ΨB.obj ((P₁.toElem.obj X).base) ⟶ (P₂.toElem.obj (Ψ.obj X)).base)
    (hsqiso : ∀ X : C₁, IsIso (sqm X))
    (hsq : ∀ {X Y : C₁} (f : X ⟶ Y), ΨB.map (P₁.Base f) ≫ sqm Y = sqm X ≫ P₂.Base (Ψ.map f))
    (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (hdiv : ∀ {X Y : C₁} (f : X ⟶ Y),
      phiIsoApp ΨB η ((P₁.toElem.obj X).base) (P₁.Div f)
        = Φ₂.map (sqm X) (P₂.Div (Ψ.map f)))
    (hdeg : ∀ {X Y : C₁} (f : X ⟶ Y), P₂.degFr (Ψ.map f) = P₁.degFr f)
    (Fc₂ : FrobenioidCore (biratPre P₂ G₂)) (hds₂ : IsDivSlim Φ₂)
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f)
    (hPBb : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y), IsPullBack (biratPre P₁ G₁) f →
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f))
    (hPBb' : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y),
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f) →
        IsPullBack (biratPre P₁ G₁) f)
    (hDivEq : ∀ {X Y : C₁} (f g : X ⟶ Y), DivEquivalent P₁ f g →
      DivEquivalent P₂ (Ψ.map f) (Ψ.map g))
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hBIso : ∀ {X Y : C₁} (f : X ⟶ Y), IsBaseIsomorphism P₁ f →
      IsBaseIsomorphism P₂ (Ψ.map f))
    (d : D₁) (y : Gp (Φ₁.val d)) (hy : y ∈ phiBiratOn G₁ d) :
    gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d) := by
  set A : C₁ := biratBaseObj G₁ d with hA
  set e : (P₁.toElem.obj A).base ≅ d := biratBaseIso G₁ d with he
  haveI := hsqiso A
  obtain ⟨δ, hδ, hδeq⟩ := (mem_phiBiratOn_iff G₁ hiso₁ e y).mp hy
  have hδ₂ := biratPsi_otimes_map G₁ G₂ Ψ hfwd Fc₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso _ δ hδ
  have hδeq2 : biratDivGp (P := P₁) (G := G₁) (A := A) (B := A) δ
      = gpMap _ (Φ₁.map e.hom) y := hδeq
  have key := biratDivGp_transport G₁ G₂ Ψ hfwd ΨB sqm hsq η hdiv hdeg (A := A) (B := A) δ
  rw [hδeq2, gpMap_phiIsoApp_nat] at key
  have h2 := congrArg (fun z => gpMap _ (Φ₂.map (inv (sqm A))) z) key
  rw [gpMap_phi_comp, gpMap_phi_inv_right] at h2
  refine (mem_phiBiratOn_iff G₂ hiso₂
    ((asIso (sqm A)).symm ≪≫ ΨB.mapIso e) (gpMap _ (phiIsoApp ΨB η d) y)).mpr ?_
  exact ⟨_, hδ₂, h2.symm⟩

/-! ## ★4. `Corollary 4.11` からの供給 -/

/-- ★★★★★★**`hbirat` を `Corollary 4.11` から導く** ——
`η := Ψ_Φ`(`Corollary 4.11, (iii)` の `psiPhiOnBase`)、
`η ∘ Div = Div ∘ Ψ`(`Corollary 4.11, (iv)` の `cor_4_11_iv_hdivc`)を入れるだけ。

★★これで `Cor54.lean` の `hbirat` は**仮定でなくなった**。 -/
theorem phiBiratOn_transport_of_cor411 (Ψe : C₁ ≌ C₂)
    (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψe.functor.map f))
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hiso : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    (hdivS : ∀ (Y : C₁) (a : Φ₁.val (P₁.toElem.obj Y).base),
      ∃ u : OTri P₁ Y, P₁.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hdivS₂ : ∀ (Y : C₂) (a : Φ₂.val (P₂.toElem.obj Y).base),
      ∃ u : OTri P₂ Y, P₂.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hOTri : ∀ (Z : C₁) (δ : End Z), δ ∈ OTri P₁ Z →
      ((Ψe.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψe.functor.obj Z))
        ∈ OTri P₂ (Ψe.functor.obj Z))
    (hOTri' : ∀ (Z : C₁) (δ : End Z),
      ((Ψe.functor.map ((((δ : End Z)) : Z ⟶ Z))) : End (Ψe.functor.obj Z))
        ∈ OTri P₂ (Ψe.functor.obj Z) → δ ∈ OTri P₁ Z)
    (hperfM : ∀ Y : C₁, IsPerfectMonoid (Φ₁.val (P₁.toElem.obj Y).base))
    (hdeg : ∀ {X Y : C₁} (g : X ⟶ Y), P₂.degFr (Ψe.functor.map g) = P₁.degFr g)
    (hFT : ∀ {X Y : C₁} (g : X ⟶ Y), IsFrobeniusType P₁ g →
      IsFrobeniusType P₂ (Ψe.functor.map g))
    (hPS : ∀ {X Y : C₁} (g : X ⟶ Y), IsPreStep P₁ g → IsPreStep P₂ (Ψe.functor.map g))
    (hPB : ∀ {X Y : C₁} (g : X ⟶ Y), IsPullBack P₁ g → IsPullBack P₂ (Ψe.functor.map g))
    (ΨB : D₁ ⥤ D₂) (sq : P₁.proj ⋙ ΨB ≅ Ψe.functor ⋙ P₂.proj)
    (Fc₂ : FrobenioidCore (biratPre P₂ G₂)) (hds₂ : IsDivSlim Φ₂)
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψe.functor.map f) → coaPreProp P₁ f)
    (hPBb : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y), IsPullBack (biratPre P₁ G₁) f →
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψe.functor hfwd).map f))
    (hPBb' : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y),
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψe.functor hfwd).map f) →
        IsPullBack (biratPre P₁ G₁) f)
    (hDivEq : ∀ {X Y : C₁} (f g : X ⟶ Y), DivEquivalent P₁ f g →
      DivEquivalent P₂ (Ψe.functor.map f) (Ψe.functor.map g))
    (hBIso : ∀ {X Y : C₁} (f : X ⟶ Y), IsBaseIsomorphism P₁ f →
      IsBaseIsomorphism P₂ (Ψe.functor.map f))
    (d : D₁) (y : Gp (Φ₁.val d)) (hy : y ∈ phiBiratOn G₁ d) :
    gpMap _ (phiIsoApp ΨB
        (psiPhiOnBase Ψe.functor F₁ ΨB sq
          (psiPhi Ψe G₁ G₂ F₁ hiso hiso₂ hdivS hdivS₂ hOTri hOTri' hperfM hdeg)) d) y
      ∈ phiBiratOn G₂ (ΨB.obj d) :=
  phiBiratOn_transport G₁ G₂ Ψe.functor hfwd hiso hiso₂ ΨB
    (fun X => sq.hom.app X) (fun X => (sq.app X).isIso_hom)
    (fun f => sq.hom.naturality f) _
    (fun f => cor_4_11_iv_hdivc Ψe G₁ G₂ F₁ F₂ hiso hiso₂ hdivS hdivS₂ hOTri hOTri'
      hperfM hdeg hFT hPS hPB ΨB sq f)
    (fun f => hdeg f) Fc₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso d y hy

/-! ### ★出典の紐付け -/

/-- ★locator —— `Corollary 5.4` の `η ∘ Div^gp = Div^gp ∘ Ψ^birat`。 -/
def sliceDivGpOf_transport.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — η ∘ Div^gp = Div^gp ∘ Ψ^birat(代表元の層)",
    sectionId := "frdi-cor-5-4" }

/-- ★locator —— `Corollary 5.4` の `hbirat`(`Corollary 4.11` から導いたもの)。 -/
def phiBiratOn_transport_of_cor411.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ_Φ は Φ^birat を Φ^birat へ移す",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
