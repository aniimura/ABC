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

variable (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
  {Ee : Type u₄} [Category.{v₄} Ee] {G G' : Dᵒᵖ ⥤ Ee}
  (h : P.proj.op ⋙ G ≅ P.proj.op ⋙ G')

include F in
set_option maxHeartbeats 1000000 in
/-- ★★★★★**`Base` の像の間のどの射でも自然** ——
`Theorem 3.4, (v)` の 3 分解(`base_three_factor`)を `Dᵒᵖ` で使う。 -/
theorem projOp_natAll {Z W : Cᵒᵖ} (t : P.proj.op.obj Z ⟶ P.proj.op.obj W) :
    G.map t ≫ h.hom.app W = h.hom.app Z ≫ G'.map t := by
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
  rw [← hh₀ W, ← hh₀ Z]
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

include F in
/-- ★`P.proj.op` は本質的全射(`baseSurj`)。 -/
theorem projOp_essSurj : P.proj.op.EssSurj := ⟨fun Y => by
  obtain ⟨A, -, ⟨e⟩⟩ := F.baseSurj Y.unop
  exact ⟨Opposite.op A, ⟨e.op.symm⟩⟩⟩

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`P.proj.op` との前合成は関手の同型を反射する**。 -/
noncomputable def projOpPrecompIsoGen : G ≅ G' := by
  haveI := projOp_essSurj P F
  -- ★綴りを固定する(`(F ⋙ G).obj` と `G.obj (F.obj _)` の食い違いを避ける)
  let h₀ : ∀ Z : Cᵒᵖ, G.obj (P.proj.op.obj Z) ≅ G'.obj (P.proj.op.obj Z) :=
    fun Z => h.app Z
  have hh₀ : ∀ Z, (h₀ Z).hom = h.hom.app Z := fun _ => rfl
  refine NatIso.ofComponents
    (fun X => (G.mapIso (P.proj.op.objObjPreimageIso X)).symm
      ≪≫ h₀ (P.proj.op.objPreimage X)
      ≪≫ G'.mapIso (P.proj.op.objObjPreimageIso X)) (fun {X Y} f => ?_)
  set pX := P.proj.op.objPreimage X with hpX
  set pY := P.proj.op.objPreimage Y with hpY
  set eX := P.proj.op.objObjPreimageIso X with heX
  set eY := P.proj.op.objObjPreimageIso Y with heY
  have hnat : G.map (eX.hom ≫ f ≫ eY.inv) ≫ (h₀ pY).hom
      = (h₀ pX).hom ≫ G'.map (eX.hom ≫ f ≫ eY.inv) := by
    rw [hh₀ pY, hh₀ pX]
    exact projOp_natAll P F h _
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

set_option maxHeartbeats 1000000 in
/-- ★★★★★**反射した同型の `Base` の像での成分は `h` そのもの**。

★★これが `Corollary 4.11, (iv)` の `div` 成分の一致に要る。 -/
theorem projOpPrecompIsoGen_app (Z : Cᵒᵖ) :
    (projOpPrecompIsoGen P F h).hom.app (P.proj.op.obj Z) = h.hom.app Z := by
  haveI := projOp_essSurj P F
  set Y := P.proj.op.obj Z with hY
  set pY := P.proj.op.objPreimage Y with hpY
  set eY := P.proj.op.objObjPreimageIso Y with heY
  show G.map eY.inv ≫ h.hom.app pY ≫ G'.map eY.hom = h.hom.app Z
  have hnat : G.map eY.hom ≫ h.hom.app Z = h.hom.app pY ≫ G'.map eY.hom :=
    projOp_natAll P F h eY.hom
  have hEX : G.map eY.inv ≫ G.map eY.hom = 𝟙 (G.obj Y) := by
    rw [← Functor.map_comp, eY.inv_hom_id, G.map_id]
  have s1 : G.map eY.inv ≫ h.hom.app pY ≫ G'.map eY.hom
      = G.map eY.inv ≫ G.map eY.hom ≫ h.hom.app Z := by rw [hnat]; rfl
  rw [s1, ← Category.assoc, hEX, Category.id_comp]

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
    Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor :=
  projOpPrecompIsoGen P₁ F₁
    (hPhi ≪≫ (Functor.associator Ψ.op P₂.proj.op Φ₂.functor).symm
      ≪≫ Functor.isoWhiskerRight (NatIso.op sq) Φ₂.functor
      ≪≫ Functor.associator P₁.proj.op ΨB.op Φ₂.functor)

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`Ψ_Φ` の `Base A` での成分** —— `Ψ_Φ` は `hPhi` を `sq` でずらしたもの。 -/
theorem psiPhiOnBase_app_apply (F₁ : FrobenioidCore P₁) (ΨB : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (hPhi : phiOnC P₁ ≅ Ψ.op ⋙ phiOnC P₂) (A : C₁)
    (x : Φ₁.val ((P₁.toElem.obj A).base)) :
    ((psiPhiOnBase Ψ F₁ ΨB sq hPhi).hom.app
        (Opposite.op ((P₁.toElem.obj A).base))).hom x
      = Φ₂.map (sq.hom.app A) ((hPhi.hom.app (Opposite.op A)).hom x) := by
  have h1 := projOpPrecompIsoGen_app P₁ F₁
    (hPhi ≪≫ (Functor.associator Ψ.op P₂.proj.op Φ₂.functor).symm
      ≪≫ Functor.isoWhiskerRight (NatIso.op sq) Φ₂.functor
      ≪≫ Functor.associator P₁.proj.op ΨB.op Φ₂.functor) (Opposite.op A)
  have h2 := congrArg (fun t : (Φ₁.functor).obj (Opposite.op ((P₁.toElem.obj A).base)) ⟶
      _ => (AddCommMonCat.Hom.hom t) x) h1
  exact h2

def psiPhiOnBase.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iii) — Ψ_Φ を 𝒟 の上へ降ろす",
    sectionId := "frdi-cor-4-11" }

end PsiPhiDescend

/-! ## ★3. `Ψ_𝔽 : 𝔽_{Φ₁} ⥤ 𝔽_{Φ₂}` -/

section ElemFrobOver

variable {D₁ : Type u} [Category.{v} D₁] {Φ₁ : MonoidOn.{v, u, w} D₁}
  {D₂ : Type u} [Category.{v} D₂] {Φ₂ : MonoidOn.{v, u, w} D₂}

/-- ★★★★★**`Ψ_Base` に沿った `𝔽_Φ` の関手** ——
在庫の `elemFrobMap` は**同じ `𝒟` の上**の 2 つの単系のためのものなので、
底の関手 `Ψ_Base : 𝒟₁ ⥤ 𝒟₂` に沿った版を作る。

原文 (FrdI p.92):
> In particular, ΨBase, ΨΦ induce an equivalence of categories -/
def elemFrobMapOver (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ⟶ ΨB.op ⋙ Φ₂.functor) :
    ElemFrobCat Φ₁ ⥤ ElemFrobCat Φ₂ where
  obj A := ⟨ΨB.obj A.base⟩
  map {A B} φ := ⟨ΨB.map φ.base, (η.app (Opposite.op A.base)).hom φ.div, φ.deg⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ ?_ rfl
    · exact ΨB.map_id A.base
    · show (η.app (Opposite.op A.base)).hom 0 = 0
      exact map_zero _
  map_comp {A B E} φ ψ := by
    refine ElemFrobCat.Hom.ext ?_ ?_ rfl
    · exact ΨB.map_comp φ.base ψ.base
    · show (η.app (Opposite.op A.base)).hom
        (Φ₁.map φ.base ψ.div + ((ψ.deg : ℕ+) : ℕ) • φ.div) = _
      rw [map_add, map_nsmul]
      congr 1
      exact congrArg (fun t => (AddCommMonCat.Hom.hom t) ψ.div) (η.naturality φ.base.op)

/-- ★★`Ψ_𝔽` は底の射影と可換(`rfl`)。 -/
theorem elemFrobMapOver_proj (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ⟶ ΨB.op ⋙ Φ₂.functor) :
    elemFrobMapOver ΨB η ⋙ ElemFrobCat.proj = ElemFrobCat.proj ⋙ ΨB := rfl

def elemFrobMapOver.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iii) — Ψ_Base と Ψ_Φ が誘導する 𝔽_Φ の関手",
    sectionId := "frdi-cor-4-11" }

end ElemFrobOver

/-! ## ★4. `Corollary 4.11, (iv)` の 1-可換図式 -/

section CorIV

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  (Ψ : C₁ ⥤ C₂)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (iv) の 1-可換図式** ——
`P₁.toElem ⋙ Ψ_𝔽 ≅ Ψ ⋙ P₂.toElem`。

原文 (FrdI p.92):
> (iv) Suppose further that C1, C2 are of rationally standard type. Then there exists

★★3 成分に分かれる:
- `base` 成分 …… (ii) の 1-可換図式の**自然性**
- `deg` 成分 …… `Ψ` が Frobenius 次数を保つこと(`Theorem 3.4, (iii)`)
- `div` 成分 …… ★`Ψ_Φ` が `Div` を `Div` へ送ること(`divMap_div_all`) -/
noncomputable def cor_4_11_iv_square (ΨB : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨB ≅ Ψ ⋙ P₂.proj)
    (η : Φ₁.functor ⟶ ΨB.op ⋙ Φ₂.functor)
    (hdivc : ∀ {A B : C₁} (φ : A ⟶ B),
      (η.app (Opposite.op ((P₁.toElem.obj A).base))).hom (P₁.Div φ)
        = Φ₂.map (sq.hom.app A) (P₂.Div (Ψ.map φ)))
    (hdeg : ∀ {A B : C₁} (φ : A ⟶ B), P₂.degFr (Ψ.map φ) = P₁.degFr φ) :
    P₁.toElem ⋙ elemFrobMapOver ΨB η ≅ Ψ ⋙ P₂.toElem :=
  NatIso.ofComponents
    (fun A =>
      { hom := ⟨sq.hom.app A, 0, 1⟩
        inv := ⟨sq.inv.app A, 0, 1⟩
        hom_inv_id := by
          refine ElemFrobCat.Hom.ext ?_ ?_ ?_
          · exact sq.hom_inv_id_app A
          · simp [ElemFrobCat.Hom.comp]
          · simp [ElemFrobCat.Hom.comp]
        inv_hom_id := by
          refine ElemFrobCat.Hom.ext ?_ ?_ ?_
          · exact sq.inv_hom_id_app A
          · simp [ElemFrobCat.Hom.comp]
          · simp [ElemFrobCat.Hom.comp] })
    (fun {A B} φ => by
      refine ElemFrobCat.Hom.ext ?_ ?_ ?_
      · exact sq.hom.naturality φ
      · simp only [Functor.comp_map, ElemFrobCat.comp_div, map_zero, zero_add, smul_zero,
          add_zero, PNat.one_coe, one_smul]
        exact hdivc φ
      · simp only [Functor.comp_map, ElemFrobCat.comp_deg, one_mul, mul_one]
        exact (hdeg φ).symm)

def cor_4_11_iv_square.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iv) — 𝒞ᵢ → 𝔽_{Φᵢ} の 1-可換図式",
    sectionId := "frdi-cor-4-11" }

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**(iv) の `div` 成分の一致** ——
`Ψ_Φ` の成分は `divMap` を `sq` でずらしたもの(`psiPhiOnBase_app_apply`)で、
`divMap` は `Div` を `Div` へ送る(`divMap_div_all`)。 -/
theorem cor_4_11_iv_hdivc (Ψe : C₁ ≌ C₂) (G₁ : Frobenioid P₁) (G₂ : Frobenioid P₂)
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
    {A B : C₁} (φ : A ⟶ B) :
    ((psiPhiOnBase Ψe.functor F₁ ΨB sq
        (psiPhi Ψe G₁ G₂ F₁ hiso hiso₂ hdivS hdivS₂ hOTri hOTri' hperfM hdeg)).hom.app
      (Opposite.op ((P₁.toElem.obj A).base))).hom (P₁.Div φ)
      = Φ₂.map (sq.hom.app A) (P₂.Div (Ψe.functor.map φ)) := by
  rw [psiPhiOnBase_app_apply, psiPhi_app_apply,
    divMap_div_all Ψe G₁ F₁ F₂ hiso hdivS hOTri hperfM hdeg hFT hPS hPB φ]

def cor_4_11_iv_hdivc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (iv) — div 成分の一致",
    sectionId := "frdi-cor-4-11" }

end CorIV

end ABC3.Found.FrdI
