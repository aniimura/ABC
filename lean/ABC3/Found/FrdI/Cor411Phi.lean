/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411BaseFn

/-!
# [FrdI] Corollary 4.11, (iii) —— `Ψ_Φ` を `𝒟` へ降ろす

原文 (FrdI p.92):
> (iii) Suppose further that C1, C2 are of rationally standard type. Then there

## ★★★★★測って分かったこと(2026-08-19)

`Theorem 4.9` の `psiPhi` は **`𝒞ᵢ` の上の関手**の同型
(`phiOnC P₁ ≅ Ψ.op ⋙ phiOnC P₂`)である。
★(iii) が要求するのは **`𝒟ᵢ` の上**、すなわち `Φ₁.functor ≅ Ψ_Base.op ⋙ Φ₂.functor`。

★★(ii) の 1-可換図式 `P₁.proj ⋙ Ψ_Base ≅ Ψ ⋙ P₂.proj` で書き換えると
`psiPhi` は `P₁.proj.op ⋙ Φ₁.functor ≅ P₁.proj.op ⋙ (Ψ_Base.op ⋙ Φ₂.functor)` になるので、
★**`P.proj.op` との前合成が同型を反射する**ことを示せば降りる。

★その中身は `Theorem 3.4, (v)` の 3 分解(`base_three_factor`)——
`𝒟` のどの射も `Base` の像と(`Base` の像である)同型の逆の合成で書けるので、
`Dᵒᵖ` でも同じ(向きが逆になるだけ)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section ProjOpPrecomp

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

universe v₄ u₄

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`P.proj.op` との前合成は関手の同型を反射する**。 -/
noncomputable def projOpPrecompIsoGen (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    {Ee : Type u₄} [Category.{v₄} Ee] {G G' : Dᵒᵖ ⥤ Ee}
    (h : P.proj.op ⋙ G ≅ P.proj.op ⋙ G') : G ≅ G' := by
  let h₀ : ∀ Z : Cᵒᵖ, G.obj (P.proj.op.obj Z) ≅ G'.obj (P.proj.op.obj Z) :=
    fun Z => h.app Z
  have hh₀ : ∀ Z, (h₀ Z).hom = h.hom.app Z := fun _ => rfl
  have hnBase : ∀ {Z W : Cᵒᵖ} (g : Z ⟶ W),
      G.map (P.proj.op.map g) ≫ (h₀ W).hom = (h₀ Z).hom ≫ G'.map (P.proj.op.map g) := by
    intro Z W g
    rw [hh₀ W, hh₀ Z]
    exact h.hom.naturality g
  have hnComp : ∀ {Z W V : Cᵒᵖ} (t : P.proj.op.obj Z ⟶ P.proj.op.obj W)
      (s : P.proj.op.obj W ⟶ P.proj.op.obj V),
      (G.map t ≫ (h₀ W).hom = (h₀ Z).hom ≫ G'.map t) →
      (G.map s ≫ (h₀ V).hom = (h₀ W).hom ≫ G'.map s) →
      G.map (t ≫ s) ≫ (h₀ V).hom = (h₀ Z).hom ≫ G'.map (t ≫ s) := by
    intro Z W V t s ht hs
    rw [G.map_comp, G'.map_comp, Category.assoc, hs, ← Category.assoc, ht, Category.assoc]
  have hnInv : ∀ {Z W : Cᵒᵖ} (t : P.proj.op.obj Z ⟶ P.proj.op.obj W) (ht : IsIso t),
      (G.map t ≫ (h₀ W).hom = (h₀ Z).hom ≫ G'.map t) →
      G.map (@inv _ _ _ _ t ht) ≫ (h₀ Z).hom
        = (h₀ W).hom ≫ G'.map (@inv _ _ _ _ t ht) := by
    intro Z W t ht hnat
    haveI := ht
    have hGi : G.map (inv t) ≫ G.map t = 𝟙 _ := by
      rw [← Functor.map_comp, IsIso.inv_hom_id, G.map_id]
    have hGi' : G'.map t ≫ G'.map (inv t) = 𝟙 _ := by
      rw [← Functor.map_comp, IsIso.hom_inv_id, G'.map_id]
    calc G.map (inv t) ≫ (h₀ Z).hom
        = G.map (inv t) ≫ (h₀ Z).hom ≫ (G'.map t ≫ G'.map (inv t)) := by
          rw [hGi', Category.comp_id]
      _ = G.map (inv t) ≫ ((h₀ Z).hom ≫ G'.map t) ≫ G'.map (inv t) := by
          simp only [Category.assoc]
      _ = G.map (inv t) ≫ (G.map t ≫ (h₀ W).hom) ≫ G'.map (inv t) := by rw [hnat]
      _ = (G.map (inv t) ≫ G.map t) ≫ (h₀ W).hom ≫ G'.map (inv t) := by
          simp only [Category.assoc]
      _ = (h₀ W).hom ≫ G'.map (inv t) := by rw [hGi, Category.id_comp]
  have hnAll : ∀ {Z W : Cᵒᵖ} (t : P.proj.op.obj Z ⟶ P.proj.op.obj W),
      G.map t ≫ (h₀ W).hom = (h₀ Z).hom ≫ G'.map t := by
    intro Z W t
    obtain ⟨X, B, a, c, ps, ha, hc, hps, hfac⟩ := base_three_factor P F t.unop
    haveI : IsIso (P.Base a) := ha.2
    have ht : t = (P.proj.op.map ps.op) ≫ (P.proj.op.map c.op)
        ≫ (@inv _ _ _ _ (P.proj.op.map a.op) (by
            show IsIso (P.Base a).op
            infer_instance)) := by
      refine Quiver.Hom.unop_inj ?_
      show t.unop = _
      rw [hfac]
      simp
      rfl
    rw [ht]
    refine hnComp _ _ (hnBase ps.op) (hnComp _ _ (hnBase c.op) ?_)
    exact hnInv (P.proj.op.map a.op) (by show IsIso (P.Base a).op; infer_instance)
      (hnBase a.op)
  haveI : P.proj.op.EssSurj := ⟨fun Y => by
    obtain ⟨A, -, ⟨e⟩⟩ := F.baseSurj Y.unop
    exact ⟨Opposite.op A, ⟨e.op.symm⟩⟩⟩
  refine NatIso.ofComponents
    (fun X => (G.mapIso (P.proj.op.objObjPreimageIso X)).symm
      ≪≫ h₀ (P.proj.op.objPreimage X)
      ≪≫ G'.mapIso (P.proj.op.objObjPreimageIso X)) (fun {X Y} f => ?_)
  set pX := P.proj.op.objPreimage X with hpX
  set pY := P.proj.op.objPreimage Y with hpY
  set eX := P.proj.op.objObjPreimageIso X with heX
  set eY := P.proj.op.objObjPreimageIso Y with heY
  have hnat : G.map (eX.hom ≫ f ≫ eY.inv) ≫ (h₀ pY).hom
      = (h₀ pX).hom ≫ G'.map (eX.hom ≫ f ≫ eY.inv) := hnAll _
  rw [G.map_comp, G.map_comp, G'.map_comp, G'.map_comp] at hnat
  show G.map f ≫ (G.map eY.inv ≫ (h₀ pY).hom ≫ G'.map eY.hom)
    = (G.map eX.inv ≫ (h₀ pX).hom ≫ G'.map eX.hom) ≫ G'.map f
  have hEY : G'.map eY.inv ≫ G'.map eY.hom = 𝟙 (G'.obj Y) := by
    rw [← Functor.map_comp, eY.inv_hom_id, G'.map_id]
  have hEX : G.map eX.hom ≫ G.map eX.inv = 𝟙 (G.obj (P.proj.op.obj pX)) := by
    rw [← Functor.map_comp, eX.hom_inv_id, G.map_id]
  refine (cancel_epi (G.map eX.hom)).mp ?_
  simp only [← Category.assoc] at hnat ⊢
  rw [hnat]
  simp only [Category.assoc, hEY, hEX, Category.comp_id]
  exact (Category.id_comp _).symm

def projOpPrecompIsoGen.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iii) — P.proj.op との前合成は関手の同型を反射する",
    sectionId := "frdi-cor-4-11" }

end ProjOpPrecomp

/-! ## ★2. `Ψ_Φ` を `𝒟` へ降ろす -/

section PsiPhiDescend

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  (Ψ : C₁ ⥤ C₂)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (iii)** —— `Ψ_Φ : Φ₁ ≅ Ψ_Base.op ⋙ Φ₂`。

原文 (FrdI p.92):
> (iii) Suppose further that C1, C2 are of rationally standard type. Then there

★★`Theorem 4.9` の `psiPhi`(`𝒞ᵢ` の上)を、(ii) の 1-可換図式で書き換えてから
`projOpPrecompIsoGen` で **`𝒟ᵢ` の上へ降ろす**。 -/
noncomputable def psiPhiOnBase (F₁ : FrobenioidCore P₁) (ΨB : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hPhi : phiOnC P₁ ≅ Ψ.op ⋙ phiOnC P₂) :
    Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor := by
  refine projOpPrecompIsoGen P₁ F₁ (hPhi ≪≫ ?_)
  refine (Functor.associator Ψ.op P₂.proj.op Φ₂.functor) ≪≫ ?_
  refine Functor.isoWhiskerRight (NatIso.op sq) Φ₂.functor ≪≫ ?_
  exact (Functor.associator P₁.proj.op ΨB.op Φ₂.functor).symm

def psiPhiOnBase.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iii) — Ψ_Φ を 𝒟 の上へ降ろす",
    sectionId := "frdi-cor-4-11" }

end PsiPhiDescend

end ABC3.Found.FrdI
