/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop16.Dict

/-!
# Prop16 —— Frobenioid core と (iv) 以降

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']

section Dict

variable {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

variable (P : PreFrobenioid C Φ) (G : D' ⥤ D)
  (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
  (hD' : IsTotallyEpimorphic D')
  (hcC : IsConnected (CfpCat P G)) (hcD' : IsConnected D')

section Core

variable (F : FrobenioidCore P)

include F in
/-- **(iii)(a)** の移送 —— co-angular は合成で閉じる。 -/
theorem cfp_coAngularComp {X Y Z : CfpCat P G} (ψ : X ⟶ Y) (φ : Y ⟶ Z) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') ψ →
      IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ →
      IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') (ψ ≫ φ) := by
  intro hψ hφ
  refine cfp_coAngular_of P G hG hD' hcC hcD' _ ?_
  exact F.coAngularComp (CfpCat.fst ψ) (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' ψ).mp hψ) ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hφ)

include F in
/-- **(iii)(b)** の移送。 -/
theorem cfp_coAngularOfPreStep {X Y : CfpCat P G} (α : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (φ : X ⟶ Y) : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ :=
  cfp_coAngular_of P G hG hD' hcC hcD' φ
    (F.coAngularOfPreStep (CfpCat.fst α) ((cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hca)
      ((cfp_preStep_iff P G hG hD' hcC hcD' α hps.2).mp hps) (CfpCat.fst φ))

include F in
/-- **(v)(a)** の移送 —— pre-step は mono。

★**両成分がそれぞれ mono** であればよい: `𝒞` 側は `F.preStepMono`、
`𝒟'` 側は **pre-step の定義から `snd` が同型**。 -/
theorem cfp_preStepMono {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) : Mono φ := by
  haveI hm : Mono (CfpCat.fst φ) :=
    F.preStepMono (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hi : IsIso (CfpCat.snd φ) := hφ.2
  constructor
  intro Z g h hgh
  have e1 : CfpCat.fst g ≫ CfpCat.fst φ = CfpCat.fst h ≫ CfpCat.fst φ :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hgh
  have e2 : CfpCat.snd g ≫ CfpCat.snd φ = CfpCat.snd h ≫ CfpCat.snd φ :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hgh
  exact InducedCategory.hom_ext
    (CommaMorphism.ext ((cancel_mono (CfpCat.fst φ)).mp e1)
      ((cancel_mono (CfpCat.snd φ)).mp e2))

include F in
/-- **(vii)(b)** の移送 —— isotropic な対象から出る射の終域は isotropic。 -/
theorem cfp_isotropicClosed {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') X) :
    IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') Y :=
  (cfp_isotropic_iff P G hG hD' hcC hcD' Y).mpr
    (F.isotropicClosed (CfpCat.fst φ) ((cfp_isotropic_iff P G hG hD' hcC hcD' X).mp h))

include F in
/-- **(ii)** の移送 —— 各次数の Frobenius 型射が存在する。

★Frobenius 型は base-isomorphism なので**鎖**があり、
新しい対象 `B` の `𝒟'` 成分は `A` のものを流用して `snd φ = 𝟙` に取れる。 -/
theorem cfp_frobDegSurj (A : CfpCat P G) (n : ℕ+) :
    ∃ (B : CfpCat P G) (φ : A ⟶ B),
      IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') φ ∧
        (cfpPreFrobenioid P G hG hD' hcC hcD').degFr φ = n := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨B₀, φ₀, hft, hdeg⟩ := F.frobDegSurj A.obj.left n
  haveI hφb : IsIso (P.proj.map φ₀) := hft.2
  have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
  refine ⟨⟨⟨B₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, hdeg⟩
  exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft

include F in
/-- **(v)(b)** の移送 —— pre-step は「co-angular pre-step ≫ isometric pre-step」に分解する。

★中間対象は**両側とも pre-step に挟まれる**ので鎖がある。 -/
theorem cfp_preStepFactor {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβc, hβs, hαi, hαs⟩ :=
    F.preStepFactor (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hβb : IsIso (P.proj.map β₀) := hβs.2
  have hzi : IsIso (inv (P.proj.map β₀) ≫ X.obj.hom) := inferInstance
  refine ⟨⟨⟨Z₀, X.obj.right, inv (P.proj.map β₀) ≫ X.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩, ?_, ?_,
    ⟨hβs.1, by show IsIso (𝟙 _); infer_instance⟩, ?_, ?_⟩
  · show P.proj.map α₀ ≫ Y.obj.hom
      = (inv (P.proj.map β₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, ← cfp_square φ,
      show P.proj.map (CfpCat.fst φ) = P.proj.map β₀ ≫ P.proj.map α₀ from by
        rw [hfac, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hβc
  · exact (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hαi
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(v)(c)** の移送 —— pre-step は「isometric pre-step ≫ co-angular pre-step」に分解する。 -/
theorem cfp_preStepFactor' {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    ∃ (Z : CfpCat P G) (β : X ⟶ Z) (α : Z ⟶ Y),
      φ = β ≫ α ∧ IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨Z₀, β₀, α₀, hfac, hβi, hβs, hαc, hαs⟩ :=
    F.preStepFactor' (CfpCat.fst φ) ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
  haveI hβb : IsIso (P.proj.map β₀) := hβs.2
  have hzi : IsIso (inv (P.proj.map β₀) ≫ X.obj.hom) := inferInstance
  refine ⟨⟨⟨Z₀, X.obj.right, inv (P.proj.map β₀) ≫ X.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩, ?_, ?_,
    ⟨hβs.1, by show IsIso (𝟙 _); infer_instance⟩, ?_, ?_⟩
  · show P.proj.map α₀ ≫ Y.obj.hom
      = (inv (P.proj.map β₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, ← cfp_square φ,
      show P.proj.map (CfpCat.fst φ) = P.proj.map β₀ ≫ P.proj.map α₀ from by
        rw [hfac, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβi
  · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hαc
  · exact ⟨hαs.1, hφ.2⟩

include F in
/-- **(i)(a)** の移送 —— `𝒟'` のどの対象の上にも Frobenius-trivial な対象がある。

★`𝒞` の `baseSurj` を `G Y` に当てて得た同型を、そのまま CFP の三つ組の第3成分にする。
★**新しい対象を作るのに鎖は要らない** —— 同型が入力として与えられるから。 -/
theorem cfp_baseSurj (Y : D') :
    ∃ A : CfpCat P G, IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ∧
      Nonempty (((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj A).base ≅ Y) := by
  obtain ⟨A₀, hft, ⟨e⟩⟩ := F.baseSurj (G.obj Y)
  haveI : IsIso e.hom := e.isIso_hom
  refine ⟨⟨⟨A₀, Y, e.hom⟩, inferInstanceAs (IsIso e.hom)⟩, ?_, ⟨Iso.refl _⟩⟩
  exact (cfp_frobTrivial_iff P G hG hD' hcC hcD' _).mpr hft

include F in
/-- **(ii)** の移送 —— 同じ次数の Frobenius 型射の本質的一意性。

★★**`𝒟'` 成分は `(snd φ)⁻¹ ≫ snd ψ` に取れる** —— Frobenius 型は base-isomorphism なので
両方の `𝒟'` 成分が同型であり、**`G` の充満性は要らない**。 -/
theorem cfp_frobDegUniq (A B E : CfpCat P G) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hψ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') ψ)
    (hd : (cfpPreFrobenioid P G hG hD' hcC hcD').degFr φ = (cfpPreFrobenioid P G hG hD' hcC hcD').degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  haveI hE : IsIso E.obj.hom := E.property
  haveI hsφ : IsIso (CfpCat.snd φ) := hφ.2
  haveI hsψ : IsIso (CfpCat.snd ψ) := hψ.2
  haveI hpφ : IsIso (P.proj.map (CfpCat.fst φ)) := cfp_baseIso_fst P G hG hD' hcC hcD' φ hφ.2
  obtain ⟨β₀, hβiso, hβ⟩ := F.frobDegUniq A.obj.left B.obj.left E.obj.left
    (CfpCat.fst φ) (CfpCat.fst ψ)
    ((cfp_frobType_iff P G hG hD' hcC hcD' φ hφ.2).mp hφ)
    ((cfp_frobType_iff P G hG hD' hcC hcD' ψ hψ.2).mp hψ) hd
  have hsq : P.proj.map β₀ ≫ E.obj.hom
      = B.obj.hom ≫ G.map (inv (CfpCat.snd φ) ≫ CfpCat.snd ψ) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst φ))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hβ, cfp_square ψ, ← Category.assoc,
      cfp_square φ, Category.assoc, ← G.map_comp, ← Category.assoc,
      IsIso.hom_inv_id, Category.id_comp]
  have hsnd : IsIso (inv (CfpCat.snd φ) ≫ CfpCat.snd ψ) := inferInstance
  refine ⟨InducedCategory.homMk ⟨β₀, inv (CfpCat.snd φ) ≫ CfpCat.snd ψ, hsq⟩,
    cfp_isIso_of P G _ hβiso hsnd, ?_⟩
  refine InducedCategory.hom_ext (CommaMorphism.ext hβ ?_)
  show CfpCat.snd φ ≫ inv (CfpCat.snd φ) ≫ CfpCat.snd ψ = CfpCat.snd ψ
  rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(vii)(a)** の移送 —— isotropic hull の存在(普遍性つき)。

★hull は isometric pre-step なので**鎖**があり、新しい対象が持ち上がる。
★普遍性の `∃!` は、**`𝒟'` 成分が `snd γ` に一意に決まる**ので `𝒞` の一意性から出る。 -/
theorem cfp_isotropicHullExists (A : CfpCat P G) :
    ∃ (B : CfpCat P G) (α : A ⟶ B), IsIsotropicHull (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨B₀, α₀, hαi, hαs, hBiso, huniv⟩ := F.isotropicHullExists A.obj.left
  haveI hαb : IsIso (P.proj.map α₀) := hαs.2
  have hzi : IsIso (inv (P.proj.map α₀) ≫ A.obj.hom) := inferInstance
  refine ⟨⟨⟨B₀, A.obj.right, inv (P.proj.map α₀) ≫ A.obj.hom⟩, hzi⟩,
    InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩,
    (cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hαi,
    ⟨hαs.1, by show IsIso (𝟙 _); infer_instance⟩,
    (cfp_isotropic_iff P G hG hD' hcC hcD' _).mpr hBiso, ?_⟩
  intro Cc hCc γ
  haveI hC : IsIso Cc.obj.hom := Cc.property
  obtain ⟨β₀, hβ₀, hβ₀u⟩ := huniv Cc.obj.left ((cfp_isotropic_iff P G hG hD' hcC hcD' Cc).mp hCc)
    (CfpCat.fst γ)
  have hsq : P.proj.map β₀ ≫ Cc.obj.hom
      = (inv (P.proj.map α₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd γ) := by
    rw [Category.assoc, ← cfp_square γ,
      show P.proj.map (CfpCat.fst γ) = P.proj.map α₀ ≫ P.proj.map β₀ from by
        rw [hβ₀, P.proj.map_comp],
      ← Category.assoc, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  refine ⟨InducedCategory.homMk ⟨β₀, (CfpCat.snd γ : A.obj.right ⟶ Cc.obj.right), hsq⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext hβ₀ ?_)
    show CfpCat.snd γ = 𝟙 _ ≫ CfpCat.snd γ
    simp
  · intro β hβ
    have hf : CfpCat.fst γ = α₀ ≫ CfpCat.fst β :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hβ
    have hs : CfpCat.snd γ = 𝟙 _ ≫ CfpCat.snd β :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hβ
    exact InducedCategory.hom_ext
      (CommaMorphism.ext (hβ₀u _ hf) (hs.trans (Category.id_comp _)).symm)

include F in
/-- ★`𝒪^▷` の元は `𝒞'` へそのまま持ち上がる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_otri_lift (A : CfpCat P G) (x : End A.obj.left) (hx : x ∈ OTri P A.obj.left) :
    ∃ y : End A, y ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') A ∧ CfpCat.fst y = x := by
  haveI hA : IsIso A.obj.hom := A.property
  have hsq : P.proj.map (x : A.obj.left ⟶ A.obj.left) ≫ A.obj.hom
      = A.obj.hom ≫ G.map (𝟙 A.obj.right) := by
    rw [show P.proj.map (x : A.obj.left ⟶ A.obj.left) = 𝟙 _ from hx.1.trans (P.Base_id _),
      G.map_id, Category.comp_id, Category.id_comp]
  exact ⟨InducedCategory.homMk ⟨(x : A.obj.left ⟶ A.obj.left), 𝟙 _, hsq⟩,
    ⟨by show CfpCat.snd _ = 𝟙 _; rfl, hx.2⟩, rfl⟩

include F in
/-- **(iii)(c)** 順方向の移送。

★`𝒪^▷` の元は base-identity なので **`𝒟'` 成分が `𝟙` に固定**され、
`𝒞` の `∃!` がそのまま上がる。 -/
theorem cfp_otriFwd {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (α : End X) (hα : α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X) :
    ∃! β : End Y, β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y ∧
      (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ := by
  obtain ⟨β₀, ⟨hβ₀m, hβ₀e⟩, hβ₀u⟩ := F.otriFwd (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
    (CfpCat.fst (α : End X))
    ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hα.1, hα.2⟩
  obtain ⟨β, hβm, hβf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F Y β₀ hβ₀m
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hα.1
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβm.1
  refine ⟨β, ⟨hβm, ?_⟩, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst φ ≫ CfpCat.fst (β : End Y) = CfpCat.fst (α : End X) ≫ CfpCat.fst φ
      rw [hβf]
      exact hβ₀e
    · show CfpCat.snd φ ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ
      rw [hβs, hαs, Category.comp_id, Category.id_comp]
  · rintro γ ⟨hγm, hγe⟩
    have hγf : CfpCat.fst φ ≫ CfpCat.fst (γ : End Y)
        = CfpCat.fst (α : End X) ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hγe
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst (γ : End Y) = CfpCat.fst (β : End Y)
      rw [hβf]
      exact hβ₀u _ ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hγm.1, hγm.2⟩, hγf⟩
    · show CfpCat.snd (γ : End Y) = CfpCat.snd (β : End Y)
      rw [hβs]
      exact hγm.1

include F in
/-- **(iii)(c)** 逆方向の移送。 -/
theorem cfp_otriBwd {X Y : CfpCat P G} (φ : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (β : End Y) (hβ : β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y) :
    ∃! α : End X, α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X ∧
      (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ := by
  obtain ⟨α₀, ⟨hα₀m, hα₀e⟩, hα₀u⟩ := F.otriBwd (CfpCat.fst φ)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
    (CfpCat.fst (β : End Y))
    ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hβ.1, hβ.2⟩
  obtain ⟨α, hαm, hαf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F X α₀ hα₀m
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβ.1
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hαm.1
  refine ⟨α, ⟨hαm, ?_⟩, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst φ ≫ CfpCat.fst (β : End Y) = CfpCat.fst (α : End X) ≫ CfpCat.fst φ
      rw [hαf]
      exact hα₀e
    · show CfpCat.snd φ ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ
      rw [hβs, hαs, Category.comp_id, Category.id_comp]
  · rintro γ ⟨hγm, hγe⟩
    have hγf : CfpCat.fst φ ≫ CfpCat.fst (β : End Y)
        = CfpCat.fst (γ : End X) ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hγe
    refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
    · show CfpCat.fst (γ : End X) = CfpCat.fst (α : End X)
      rw [hαf]
      exact hα₀u _ ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hγm.1, hγm.2⟩, hγf⟩
    · show CfpCat.snd (γ : End X) = CfpCat.snd (α : End X)
      rw [hαs]
      exact hγm.1

include F in
/-- **(vi)** の移送 —— `𝒪^×` を除いた忠実性。

★`Div` の移送には **`Φ(α⁻¹)` の単射性**（`Definition 1.1, (ii), (a)`）を使う。 -/
theorem cfp_faithfulUpToUnits {X Y : CfpCat P G} (φ ψ : X ⟶ Y)
    (hb : BaseEquivalent (cfpPreFrobenioid P G hG hD' hcC hcD') φ ψ)
    (hm : MetricallyEquivalent (cfpPreFrobenioid P G hG hD' hcC hcD') φ ψ)
    (hφc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hφs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hψc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') ψ)
    (hψs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') ψ) :
    ∃ α : End Y, α ∈ OTimes (cfpPreFrobenioid P G hG hD' hcC hcD') Y ∧
      φ = ψ ≫ (α : Y ⟶ Y) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  have hb' : CfpCat.snd φ = CfpCat.snd ψ := hb
  have hmC : MetricallyEquivalent P (CfpCat.fst φ) (CfpCat.fst ψ) :=
    Φ.map_injective (@inv _ _ _ _ X.obj.hom hX) hm
  have hbC : BaseEquivalent P (CfpCat.fst φ) (CfpCat.fst ψ) := by
    show P.proj.map (CfpCat.fst φ) = P.proj.map (CfpCat.fst ψ)
    rw [cfp_base_fst P G φ, cfp_base_fst P G ψ, hb']
  obtain ⟨α₀, hα₀, hφψ⟩ := F.faithfulUpToUnits (CfpCat.fst φ) (CfpCat.fst ψ) hbC hmC
    ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hφc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' φ hφs.2).mp hφs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' ψ).mp hψc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' ψ hψs.2).mp hψs)
  obtain ⟨α, hαm, hαf⟩ := cfp_otri_lift P G hG hD' hcC hcD' F Y α₀ hα₀.1
  haveI hα₀i : IsIso (α₀ : Y.obj.left ⟶ Y.obj.left) := (isUnit_iff_isIso _).mp hα₀.2
  have hsnd1 : CfpCat.snd (α : End Y) = 𝟙 _ := hαm.1
  have hsndi : IsIso (CfpCat.snd (α : End Y)) := by rw [hsnd1]; infer_instance
  have hfsti : IsIso (CfpCat.fst (α : End Y)) := by rw [hαf]; infer_instance
  refine ⟨α, ⟨hαm, (isUnit_iff_isIso _).mpr (cfp_isIso_of P G _ hfsti hsndi)⟩, ?_⟩
  refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
  · show CfpCat.fst φ = CfpCat.fst ψ ≫ CfpCat.fst (α : End Y)
    rw [hαf]
    exact hφψ
  · show CfpCat.snd φ = CfpCat.snd ψ ≫ CfpCat.snd (α : End Y)
    rw [hsnd1, Category.comp_id]
    exact hb'

include F in
/-- **(iii)(c)** 全単射が `Base` にしか依らないことの移送。 -/
theorem cfp_otriBase {X Y : CfpCat P G} (φ φ' : X ⟶ Y)
    (hca : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hps : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ)
    (hca' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ')
    (hps' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ')
    (hbase : (cfpPreFrobenioid P G hG hD' hcC hcD').Base φ = (cfpPreFrobenioid P G hG hD' hcC hcD').Base φ')
    (α : End X) (hα : α ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') X)
    (β : End Y) (hβ : β ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') Y)
    (h : (φ ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ) :
    (φ' ≫ β : X ⟶ Y) = (α : X ⟶ X) ≫ φ' := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  have hb' : CfpCat.snd φ = CfpCat.snd φ' := hbase
  have hαs : CfpCat.snd (α : End X) = 𝟙 _ := hα.1
  have hβs : CfpCat.snd (β : End Y) = 𝟙 _ := hβ.1
  have hf : CfpCat.fst φ ≫ CfpCat.fst (β : End Y)
      = CfpCat.fst (α : End X) ≫ CfpCat.fst φ :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h
  refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
  · refine F.otriBase (CfpCat.fst φ) (CfpCat.fst φ')
      ((cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hca)
      ((cfp_preStep_iff P G hG hD' hcC hcD' φ hps.2).mp hps)
      ((cfp_coAngular_iff P G hG hD' hcC hcD' φ').mp hca')
      ((cfp_preStep_iff P G hG hD' hcC hcD' φ' hps'.2).mp hps') ?_
      (CfpCat.fst (α : End X)) ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hα.1, hα.2⟩
      (CfpCat.fst (β : End Y)) ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hβ.1, hβ.2⟩ hf
    show P.proj.map (CfpCat.fst φ) = P.proj.map (CfpCat.fst φ')
    rw [cfp_base_fst P G φ, cfp_base_fst P G φ', hb']
  · show CfpCat.snd φ' ≫ CfpCat.snd (β : End Y) = CfpCat.snd (α : End X) ≫ CfpCat.snd φ'
    rw [hβs, hαs, Category.comp_id, Category.id_comp]

include F in
/-- **(iv)(a)** の移送 —— 任意の射の3分解。

★`γ`(Frobenius 型)と `β`(pre-step)は base-isomorphism なので**鎖**があり、
中間対象 2 つが持ち上がる。`α`(pull-back)は `𝒟'` 成分として `snd φ` を受け取る。 -/
theorem cfp_arbFactor {A B : CfpCat P G} (φ : A ⟶ B) :
    ∃ (X Y : CfpCat P G) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ, hβ, hα⟩ := F.arbFactor (CfpCat.fst φ)
  haveI hγb : IsIso (P.proj.map γ₀) := hγ.2
  haveI hβb : IsIso (P.proj.map β₀) := hβ.2
  have hxi : IsIso (inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hyi : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  refine ⟨⟨⟨X₀, A.obj.right, inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hxi⟩,
    ⟨⟨Y₀, A.obj.right, inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hyi⟩,
    InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : A.obj.right ⟶ B.obj.right), ?_⟩, ?_, ?_, ?_, ?_⟩
  · show P.proj.map α₀ ≫ B.obj.hom
      = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, Category.assoc, ← cfp_square φ, hfacp]
    simp only [Category.assoc]
    rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
      ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hγ
  · exact ⟨hβ.1, by show IsIso (𝟙 _); infer_instance⟩
  · exact cfp_isPullBack_of P G hG hD' hcC hcD' _ hα

include F in
/-- **(v)(b)** の一意性の移送。

★**同型を作るのに `𝒟'` 成分は `(snd β)⁻¹ ≫ snd β'` に取る** ——
`β`, `β'` はどちらも pre-step なので `𝒟'` 成分が同型であり、
★**逆射も成分ごとに書き下せる**ので `inv` を CFP の射に対して使わずに済む(表 #2)。 -/
theorem cfp_preStepFactorUniq {A B : CfpCat P G} (X X' : CfpCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hβs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hαi : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hαs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hβc' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hβs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hαi' : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') α')
    (hαs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hsβ : IsIso (CfpCat.snd β) := hβs.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβs'.2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) := cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) := cfp_baseIso_fst P G hG hD' hcC hcD' β' hβs'.2
  obtain ⟨γ₀, hγ1, hγ2⟩ := F.preStepFactorUniq X.obj.left X'.obj.left
    (CfpCat.fst β) (CfpCat.fst α) (CfpCat.fst β') (CfpCat.fst α')
    (congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' β).mp hβc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβs.2).mp hβs)
    ((cfp_isometric_iff P G hG hD' hcC hcD' α).mp hαi)
    ((cfp_preStep_iff P G hG hD' hcC hcD' α hαs.2).mp hαs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' β').mp hβc')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβs'.2).mp hβs')
    ((cfp_isometric_iff P G hG hD' hcC hcD' α').mp hαi')
    ((cfp_preStep_iff P G hG hD' hcC hcD' α' hαs'.2).mp hαs')
  -- ★四角形: `𝒞` 成分の底射は `𝒟'` 成分から決まる
  have hsq : P.proj.map γ₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← hγ2, cfp_square β', ← Category.assoc,
      cfp_square β, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hγ3 : CfpCat.fst β' ≫ γ₀.inv = CfpCat.fst β := by
    rw [hγ2, Category.assoc, γ₀.hom_inv_id, Category.comp_id]
  have hsq' : P.proj.map γ₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hγ3, cfp_square β, ← Category.assoc,
      cfp_square β', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  refine ⟨⟨InducedCategory.homMk ⟨γ₀.hom, inv (CfpCat.snd β) ≫ CfpCat.snd β', hsq⟩,
    InducedCategory.homMk ⟨γ₀.inv, inv (CfpCat.snd β') ≫ CfpCat.snd β, hsq'⟩, ?_, ?_⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ CfpCat.snd β') ≫ inv (CfpCat.snd β') ≫ CfpCat.snd β = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β' = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ1 ?_)
    show CfpCat.snd α' = (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ CfpCat.snd α
    have hs : CfpCat.snd β ≫ CfpCat.snd α = CfpCat.snd β' ≫ CfpCat.snd α' :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
    rw [Category.assoc, hs, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ2 ?_)
    show CfpCat.snd β' = CfpCat.snd β ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(v)(c)** の一意性の移送(★`preStepFactorUniq` と同じ形)。 -/
theorem cfp_preStepFactorUniq' {A B : CfpCat P G} (X X' : CfpCat P G)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβi : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hβs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hαc : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hαs : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hβi' : IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hβs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hαc' : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') α')
    (hαs' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hsβ : IsIso (CfpCat.snd β) := hβs.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβs'.2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) := cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) := cfp_baseIso_fst P G hG hD' hcC hcD' β' hβs'.2
  obtain ⟨γ₀, hγ1, hγ2⟩ := F.preStepFactorUniq' X.obj.left X'.obj.left
    (CfpCat.fst β) (CfpCat.fst α) (CfpCat.fst β') (CfpCat.fst α')
    (congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq)
    ((cfp_isometric_iff P G hG hD' hcC hcD' β).mp hβi)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβs.2).mp hβs)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hαc)
    ((cfp_preStep_iff P G hG hD' hcC hcD' α hαs.2).mp hαs)
    ((cfp_isometric_iff P G hG hD' hcC hcD' β').mp hβi')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβs'.2).mp hβs')
    ((cfp_coAngular_iff P G hG hD' hcC hcD' α').mp hαc')
    ((cfp_preStep_iff P G hG hD' hcC hcD' α' hαs'.2).mp hαs')
  have hsq : P.proj.map γ₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← hγ2, cfp_square β', ← Category.assoc,
      cfp_square β, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hγ3 : CfpCat.fst β' ≫ γ₀.inv = CfpCat.fst β := by
    rw [hγ2, Category.assoc, γ₀.hom_inv_id, Category.comp_id]
  have hsq' : P.proj.map γ₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, hγ3, cfp_square β, ← Category.assoc,
      cfp_square β', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  refine ⟨⟨InducedCategory.homMk ⟨γ₀.hom, inv (CfpCat.snd β) ≫ CfpCat.snd β', hsq⟩,
    InducedCategory.homMk ⟨γ₀.inv, inv (CfpCat.snd β') ≫ CfpCat.snd β, hsq'⟩, ?_, ?_⟩, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ CfpCat.snd β') ≫ inv (CfpCat.snd β') ≫ CfpCat.snd β = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext γ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β' = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ1 ?_)
    show CfpCat.snd α' = (inv (CfpCat.snd β') ≫ CfpCat.snd β) ≫ CfpCat.snd α
    have hs : CfpCat.snd β ≫ CfpCat.snd α = CfpCat.snd β' ≫ CfpCat.snd α' :=
      congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
    rw [Category.assoc, hs, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hγ2 ?_)
    show CfpCat.snd β' = CfpCat.snd β ≫ inv (CfpCat.snd β) ≫ CfpCat.snd β'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- **(i)(b)** の移送 —— `𝒟'` の同型を pre-step の span で実現する。

★`𝒟'` の同型 `α` を `G` で送り、両端の同型で共役して `𝒞` の `preStepSpan` に渡す。
★中間対象の `𝒟'` 成分は `A` のものを流用し、`φ` の `𝒟'` 成分は `𝟙`、
`ψ` の `𝒟'` 成分は `α` そのものに取れる。 -/
theorem cfp_preStepSpan (A B : CfpCat P G)
    (α : ((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj A).base ⟶
      ((cfpPreFrobenioid P G hG hD' hcC hcD').toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : CfpCat P G) (φ : X ⟶ A) (ψ : X ⟶ B)
      (hφ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') φ),
      IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') ψ ∧
        α = @inv _ _ _ _ ((cfpPreFrobenioid P G hG hD' hcC hcD').Base φ) hφ.2 ≫
          (cfpPreFrobenioid P G hG hD' hcC hcD').Base ψ := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  -- ★#3: 綴りの決まった変数を先に導入する
  obtain ⟨a, rfl⟩ : ∃ a : A.obj.right ⟶ B.obj.right, a = α := ⟨α, rfl⟩
  haveI hai : IsIso a := hα
  haveI hGa : IsIso (G.map a) := inferInstance
  obtain ⟨wB, hwB1, hwB2⟩ := hB.out
  haveI hwBi : IsIso wB := ⟨B.obj.hom, hwB2, hwB1⟩
  have hui : IsIso (A.obj.hom ≫ G.map a ≫ wB) := inferInstance
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ :=
    F.preStepSpan A.obj.left B.obj.left (A.obj.hom ≫ G.map a ≫ wB) hui
  haveI hφb : IsIso (P.proj.map φ₀) := hφ₀.2
  have hxi : IsIso (P.proj.map φ₀ ≫ A.obj.hom) := inferInstance
  have h1 : A.obj.hom ≫ G.map a ≫ wB
      = @inv _ _ _ _ (P.proj.map φ₀) hφ₀.2 ≫ P.proj.map ψ₀ := heq
  have h2 : A.obj.hom ≫ G.map a
      = @inv _ _ _ _ (P.proj.map φ₀) hφ₀.2 ≫ P.proj.map ψ₀ ≫ B.obj.hom := by
    have h3 := congrArg (fun t => t ≫ B.obj.hom) h1
    simp only [Category.assoc] at h3
    rw [hwB2, Category.comp_id] at h3
    exact h3
  have hkey : P.proj.map ψ₀ ≫ B.obj.hom
      = (P.proj.map φ₀ ≫ A.obj.hom) ≫ G.map a := by
    rw [Category.assoc, h2, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  refine ⟨⟨⟨X₀, A.obj.right, P.proj.map φ₀ ≫ A.obj.hom⟩, hxi⟩,
    InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨ψ₀, a, hkey⟩,
    ⟨hφ₀.1, by show IsIso (𝟙 _); infer_instance⟩,
    ⟨hψ₀.1, hai⟩, ?_⟩
  show a = @inv _ _ _ _ (𝟙 A.obj.right) _ ≫ a
  rw [IsIso.inv_id, Category.id_comp]

/-- ★★**pull-back は左から簡約できる** —— `f ≫ w` と `w` が pull-back なら `f` も。

★★**`Definition 1.2, (ii)` の全単射条件だけから出る**(Frobenioid の公理は要らない)。
`Proposition 1.7, (v)` の pull-back の段は `FrobenioidCore` を仮定するが、
**この向きだけなら仮定なしで示せる**。
★これが無いと `plBkEquiv` の充満性が循環する(`𝒞'` が Frobenioid であることを要してしまう)。 -/
theorem isPullBack_of_comp_left {Cc : Type u2} [Category.{v2} Cc] {Ψ : MonoidOn.{v, u, w} D}
    (Q : PreFrobenioid Cc Ψ) {A B E : Cc} (f : A ⟶ B) (wm : B ⟶ E)
    (hw : IsPullBack Q wm) (hq : IsPullBack Q (f ≫ wm)) : IsPullBack Q f := by
  intro T
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have e1 : (f₁ ≫ f) = f₂ ≫ f := congrArg Prod.fst hp
    have e2 : Q.Base f₁ = Q.Base f₂ := congrArg Prod.snd hp
    refine (hq T).1 (Subtype.ext (Prod.ext ?_ e2))
    show (f₁ ≫ f ≫ wm) = f₂ ≫ f ≫ wm
    rw [← Category.assoc, e1, Category.assoc]
  · rintro ⟨⟨g, u⟩, hcond⟩
    have hcond' : Q.Base (g ≫ wm) = u ≫ Q.Base (f ≫ wm) := by
      rw [Q.Base_comp, Q.Base_comp, hcond, Category.assoc]
    obtain ⟨h, hh⟩ := (hq T).2 ⟨(g ≫ wm, u), hcond'⟩
    have hp := Subtype.ext_iff.mp hh
    have h1 : (h ≫ f ≫ wm) = g ≫ wm := congrArg Prod.fst hp
    have h2 : Q.Base h = u := congrArg Prod.snd hp
    refine ⟨h, Subtype.ext (Prod.ext ?_ h2)⟩
    refine (hw T).1 (Subtype.ext (Prod.ext ?_ ?_))
    · show ((h ≫ f) ≫ wm) = g ≫ wm
      rw [Category.assoc]; exact h1
    · show Q.Base (h ≫ f) = Q.Base g
      rw [Q.Base_comp, h2, hcond]

include F in
/-- ★**(iv)(a) の移送に「pull-back 因子の射影がもとの pull-back である」ことを足した形**。

★これが (iii) の残っていた向き(`cfp_isPullBack_to`)の入口になる。 -/
theorem cfp_arbFactor' {A B : CfpCat P G} (φ : A ⟶ B) :
    ∃ (X Y : CfpCat P G) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ ∧
        IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β ∧
        IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
        IsPullBack P (CfpCat.fst α) := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ, hβ, hα⟩ := F.arbFactor (CfpCat.fst φ)
  haveI hγb : IsIso (P.proj.map γ₀) := hγ.2
  haveI hβb : IsIso (P.proj.map β₀) := hβ.2
  have hxi : IsIso (inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hyi : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) := inferInstance
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  refine ⟨⟨⟨X₀, A.obj.right, inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hxi⟩,
    ⟨⟨Y₀, A.obj.right, inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom⟩, hyi⟩,
    InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩,
    InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : A.obj.right ⟶ B.obj.right), ?_⟩,
    ?_, ?_, ?_, ?_, hα⟩
  · show P.proj.map α₀ ≫ B.obj.hom
      = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ A.obj.hom) ≫ G.map (CfpCat.snd φ)
    rw [Category.assoc, Category.assoc, ← cfp_square φ, hfacp]
    simp only [Category.assoc]
    rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
      ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
  · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
    show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
    simp
  · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hγ
  · exact ⟨hβ.1, by show IsIso (𝟙 _); infer_instance⟩
  · exact cfp_isPullBack_of P G hG hD' hcC hcD' _ hα

include F in
/-- ★★★★**(iii) の残っていた向き** —— **`𝒞'` の pull-back は `𝒞` の pull-back**。

原文 (FrdI p.28):
> (iii) A morphism of C is a(n) isometry (respectively, morphism of a given

★★**筋**: 任意射の 3 分解 `φ = γ ≫ β ≫ α`(`cfp_arbFactor'`)で
`γ ≫ β` は base-isomorphism、`α` は pull-back。
`φ` が pull-back なら **pull-back は左から簡約できる**(`isPullBack_of_comp_left`)ので
`γ ≫ β` も pull-back になり、**底が同型な pull-back は同型**
(`isIso_of_isPullBack_of_baseIso`)だから `γ ≫ β` は同型。
★あとは `fst φ = fst (γ ≫ β) ≫ fst α`(前者は同型、後者は `𝒞` の pull-back)である。

★★この向きが無かったために `pullBackLB` と `arbFactorUniq` の移送が止まっていた。 -/
theorem cfp_isPullBack_to {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ) :
    IsPullBack P (CfpCat.fst φ) := by
  obtain ⟨Z, W, γ, β, α, hfac, hγ, hβ, hα, hα₀⟩ :=
    cfp_arbFactor' P G hG hD' hcC hcD' F φ
  haveI hsγ : IsIso (CfpCat.snd γ) := hγ.2
  haveI hsβ : IsIso (CfpCat.snd β) := hβ.2
  have hbi : IsIso ((cfpPreFrobenioid P G hG hD' hcC hcD').Base (γ ≫ β)) := by
    show IsIso (CfpCat.snd γ ≫ CfpCat.snd β)
    infer_instance
  have hpb : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') ((γ ≫ β) ≫ α) := by
    rw [Category.assoc, ← hfac]; exact h
  have h1 : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') (γ ≫ β) :=
    isPullBack_of_comp_left _ (γ ≫ β) α hα hpb
  haveI : IsIso (γ ≫ β) :=
    isIso_of_isPullBack_of_baseIso (cfpPreFrobenioid P G hG hD' hcC hcD') h1 hbi
  haveI : IsIso (CfpCat.fst (γ ≫ β)) := cfp_isIso_fst P G _ inferInstance
  have hfst : CfpCat.fst φ = CfpCat.fst (γ ≫ β) ≫ CfpCat.fst α := by
    rw [hfac, ← Category.assoc]; rfl
  rw [hfst]
  exact IsPullBack.comp P (isPullBack_of_isIso P _) hα₀

/-- **(iii)** —— pull-back は射影で決まる(両向き)。 -/
theorem cfp_isPullBack_iff (F : FrobenioidCore P) {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsPullBack P (CfpCat.fst φ) :=
  ⟨cfp_isPullBack_to P G hG hD' hcC hcD' F φ, cfp_isPullBack_of P G hG hD' hcC hcD' φ⟩

include F in
/-- **(iv)(b)** の移送 —— pull-back は LB-invertible かつ linear。 -/
theorem cfp_pullBackLB {A B : CfpCat P G} (α : A ⟶ B)
    (h : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α) :
    IsLBInvertible (cfpPreFrobenioid P G hG hD' hcC hcD') α ∧
      IsLinear (cfpPreFrobenioid P G hG hD' hcC hcD') α := by
  obtain ⟨h1, h2⟩ := F.pullBackLB (CfpCat.fst α)
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α h)
  exact ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' α).mpr h1,
    (cfp_linear_iff P G hG hD' hcC hcD' α).mpr h2⟩

include F in
/-- **(iv)(a)** の一意性の移送。

★`𝒟'` 成分は `γ`(Frobenius 型)と `β`(pre-step)の `𝒟'` 成分がどちらも同型なので、
**逆射を成分ごとに書き下せる**(分類表 #2 の回避)。 -/
theorem cfp_arbFactorUniq {A B : CfpCat P G} (X Y X' Y' : CfpCat P G)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B)
    (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : γ ≫ β ≫ α = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ)
    (hβ : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β)
    (hα : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α)
    (hγ' : IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') γ')
    (hβ' : IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') β')
    (hα' : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hX : IsIso X.obj.hom := X.property
  haveI hX' : IsIso X'.obj.hom := X'.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI hY' : IsIso Y'.obj.hom := Y'.property
  haveI hsγ : IsIso (CfpCat.snd γ) := hγ.2
  haveI hsγ' : IsIso (CfpCat.snd γ') := hγ'.2
  haveI hsβ : IsIso (CfpCat.snd β) := hβ.2
  haveI hsβ' : IsIso (CfpCat.snd β') := hβ'.2
  have hfst : CfpCat.fst γ ≫ CfpCat.fst β ≫ CfpCat.fst α
      = CfpCat.fst γ' ≫ CfpCat.fst β' ≫ CfpCat.fst α' :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) heq
  have hsnd : CfpCat.snd γ ≫ CfpCat.snd β ≫ CfpCat.snd α
      = CfpCat.snd γ' ≫ CfpCat.snd β' ≫ CfpCat.snd α' :=
    congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) heq
  obtain ⟨δ₀, ε₀, e1, e2, e3⟩ := F.arbFactorUniq X.obj.left Y.obj.left X'.obj.left Y'.obj.left
    (CfpCat.fst γ) (CfpCat.fst β) (CfpCat.fst α)
    (CfpCat.fst γ') (CfpCat.fst β') (CfpCat.fst α') hfst
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ hγ.2).mp hγ)
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβ.2).mp hβ)
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α hα)
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ' hγ'.2).mp hγ')
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβ'.2).mp hβ')
    (cfp_isPullBack_to P G hG hD' hcC hcD' F α' hα')
  haveI hpγ : IsIso (P.proj.map (CfpCat.fst γ)) :=
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ hγ.2).mp hγ).2
  haveI hpγ' : IsIso (P.proj.map (CfpCat.fst γ')) :=
    ((cfp_frobType_iff P G hG hD' hcC hcD' γ' hγ'.2).mp hγ').2
  haveI hpβ : IsIso (P.proj.map (CfpCat.fst β)) :=
    ((cfp_preStep_iff P G hG hD' hcC hcD' β hβ.2).mp hβ).2
  haveI hpβ' : IsIso (P.proj.map (CfpCat.fst β')) :=
    ((cfp_preStep_iff P G hG hD' hcC hcD' β' hβ'.2).mp hβ').2
  have e3' : CfpCat.fst γ = CfpCat.fst γ' ≫ ε₀.inv := by
    rw [e3, Category.assoc, ε₀.hom_inv_id, Category.comp_id]
  have e2' : CfpCat.fst β ≫ δ₀.hom = ε₀.hom ≫ CfpCat.fst β' := by
    rw [e2, ← Category.assoc, ε₀.hom_inv_id, Category.id_comp]
  have e2'' : CfpCat.fst β' ≫ δ₀.inv = ε₀.inv ≫ CfpCat.fst β := by
    rw [e2, Category.assoc, Category.assoc, δ₀.hom_inv_id, Category.comp_id]
  have hsqε : P.proj.map ε₀.hom ≫ X'.obj.hom
      = X.obj.hom ≫ G.map (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst γ))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← e3, cfp_square γ', ← Category.assoc,
      cfp_square γ, Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hsqε' : P.proj.map ε₀.inv ≫ X.obj.hom
      = X'.obj.hom ≫ G.map (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst γ'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, ← e3', cfp_square γ, ← Category.assoc,
      cfp_square γ', Category.assoc, ← G.map_comp, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]
  have hsqδ : P.proj.map δ₀.hom ≫ Y'.obj.hom
      = Y.obj.hom ≫ G.map (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ')
          ≫ CfpCat.snd β') := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, e2', P.proj.map_comp, Category.assoc,
      cfp_square β', ← Category.assoc, hsqε, ← Category.assoc, cfp_square β]
    simp only [Category.assoc, ← G.map_comp]
    rw [← Category.assoc (CfpCat.snd β), IsIso.hom_inv_id, Category.id_comp]
  have hsqδ' : P.proj.map δ₀.inv ≫ Y.obj.hom
      = Y'.obj.hom ≫ G.map (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ)
          ≫ CfpCat.snd β) := by
    refine (cancel_epi (P.proj.map (CfpCat.fst β'))).mp ?_
    rw [← Category.assoc, ← P.proj.map_comp, e2'', P.proj.map_comp, Category.assoc,
      cfp_square β, ← Category.assoc, hsqε', ← Category.assoc, cfp_square β']
    simp only [Category.assoc, ← G.map_comp]
    rw [← Category.assoc (CfpCat.snd β'), IsIso.hom_inv_id, Category.id_comp]
  have hsα : CfpCat.snd α'
      = (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β)
        ≫ CfpCat.snd α := by
    simp only [Category.assoc]
    rw [hsnd]
    simp
  refine ⟨⟨InducedCategory.homMk ⟨δ₀.hom, inv (CfpCat.snd β) ≫
      (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β', hsqδ⟩,
    InducedCategory.homMk ⟨δ₀.inv, inv (CfpCat.snd β') ≫
      (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β, hsqδ'⟩, ?_, ?_⟩,
    ⟨InducedCategory.homMk ⟨ε₀.hom, inv (CfpCat.snd γ) ≫ CfpCat.snd γ', hsqε⟩,
    InducedCategory.homMk ⟨ε₀.inv, inv (CfpCat.snd γ') ≫ CfpCat.snd γ, hsqε'⟩, ?_, ?_⟩,
    ?_, ?_, ?_⟩
  · refine InducedCategory.hom_ext (CommaMorphism.ext δ₀.hom_inv_id ?_)
    show (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β')
      ≫ (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β) = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext δ₀.inv_hom_id ?_)
    show (inv (CfpCat.snd β') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β)
      ≫ (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β') = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext ε₀.hom_inv_id ?_)
    show (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) = 𝟙 _
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext ε₀.inv_hom_id ?_)
    show (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') = 𝟙 _
    simp
  · exact InducedCategory.hom_ext (CommaMorphism.ext e1 hsα)
  · refine InducedCategory.hom_ext (CommaMorphism.ext e2 ?_)
    show CfpCat.snd β' = (inv (CfpCat.snd γ') ≫ CfpCat.snd γ) ≫ CfpCat.snd β
      ≫ (inv (CfpCat.snd β) ≫ (inv (CfpCat.snd γ) ≫ CfpCat.snd γ') ≫ CfpCat.snd β')
    simp
  · refine InducedCategory.hom_ext (CommaMorphism.ext e3 ?_)
    show CfpCat.snd γ' = CfpCat.snd γ ≫ inv (CfpCat.snd γ) ≫ CfpCat.snd γ'
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

include F in
/-- ★★★★★★**(i)(c) の移送** —— `(𝒞'^pl-bk)_A → 𝒟'_{A_𝒟'}` は圏同値。

★★**忠実性・充満性は `𝒞'` の pull-back 性(対象の定義そのもの)を直接使う**。
充満性で要る「pull-back の左簡約」は `isPullBack_of_comp_left`(仮定なし)なので
**循環しない**。本質的全射性は `𝒞` 側の基底変換(`plBk_baseChange`)だけで足りる。

★★★**2026-08-20: ここで止まっていた原因が分かった** ——
`P.proj.obj X` と `(P.toElem.obj X).base` は `def proj` を展開すれば等しいが、
★**`rw` / `simp` は書き換え後の項を `instances` 透明度で型検査するので、
この 2 つが混ざった目標では「型が正しくない」と言って止まる。**
(以前の記録が挙げた候補 A(`inv` のインスタンス)・候補 B(`rw` の位置)は
どちらも原因ではなかった。)
★**対処は `rw` を使わないこと** —— `calc` と `congrArg` / `Category.assoc` を
**項として**書けば、既定の透明度で defeq が通る。 -/
theorem cfp_plBkEquiv (A : CfpCat P G) :
    (plBkOverFunctor (cfpPreFrobenioid P G hG hD' hcC hcD') A).IsEquivalence := by
  haveI hfaith : (plBkOverFunctor (cfpPreFrobenioid P G hG hD' hcC hcD') A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : (cfpPreFrobenioid P G hG hD' hcC hcD').Base f.left.hom
        = (cfpPreFrobenioid P G hG hD' hcC hcD').Base g.left.hom :=
      congrArg CommaMorphism.left hfg
    have hwf : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have hwg : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      (IsPullBack.hom_ext _ W.hom.property _ _ (hwf.trans hwg.symm) hb))
  haveI hfull : (plBkOverFunctor (cfpPreFrobenioid P G hG hD' hcC hcD') A).Full := by
    constructor
    intro Z W h
    obtain ⟨f₀, hf₀, -⟩ := IsPullBack.lift _ W.hom.property Z.left.obj Z.hom.hom h.left
      (Over.w h).symm
    have hpb : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') f₀ :=
      isPullBack_of_comp_left _ f₀ W.hom.hom W.hom.property
        (by rw [hf₀.1]; exact Z.hom.property)
    exact ⟨Over.homMk (⟨f₀, hpb⟩ : Z.left ⟶ W.left)
        (InducedWideCategory.Hom.ext hf₀.1),
      Over.OverMorphism.ext hf₀.2⟩
  haveI hess : (plBkOverFunctor (cfpPreFrobenioid P G hG hD' hcC hcD') A).EssSurj := by
    constructor
    intro Y
    haveI hA : IsIso A.obj.hom := A.property
    obtain ⟨q, hq⟩ : ∃ q : Y.left ⟶ A.obj.right, q = Y.hom := ⟨Y.hom, rfl⟩
    obtain ⟨w, hwq⟩ : ∃ w : G.obj Y.left ⟶ (P.toElem.obj A.obj.left).base,
        w = G.map q ≫ inv A.obj.hom := ⟨_, rfl⟩
    have hcancel : w ≫ A.obj.hom = G.map q :=
      calc w ≫ A.obj.hom
          = (G.map q ≫ inv A.obj.hom) ≫ A.obj.hom :=
            congrArg (fun t => t ≫ A.obj.hom) hwq
        _ = G.map q ≫ (inv A.obj.hom ≫ A.obj.hom) := Category.assoc _ _ _
        _ = G.map q ≫ 𝟙 _ := congrArg (fun t => G.map q ≫ t) (IsIso.inv_hom_id A.obj.hom)
        _ = G.map q := Category.comp_id _
    obtain ⟨Yt, αt, k, hpb, hw⟩ := plBk_baseChange P F A.obj.left w
    haveI hki : IsIso k.hom := inferInstance
    have hsq : P.proj.map αt ≫ A.obj.hom = k.hom ≫ G.map q :=
      calc P.proj.map αt ≫ A.obj.hom
          = (k.hom ≫ w) ≫ A.obj.hom := congrArg (fun t => t ≫ A.obj.hom) hw
        _ = k.hom ≫ (w ≫ A.obj.hom) := Category.assoc _ _ _
        _ = k.hom ≫ G.map q := congrArg (fun t => k.hom ≫ t) hcancel
    refine ⟨Over.mk (show (⟨(⟨⟨Yt, Y.left, k.hom⟩, hki⟩ : CfpCat P G)⟩ :
        PlBk (cfpPreFrobenioid P G hG hD' hcC hcD'))
          ⟶ (⟨A⟩ : PlBk (cfpPreFrobenioid P G hG hD' hcC hcD')) from
      ⟨InducedCategory.homMk ⟨αt, q, hsq⟩,
        cfp_isPullBack_of P G hG hD' hcC hcD' _ hpb⟩), ⟨?_⟩⟩
    refine Over.isoMk (Iso.refl _) ?_
    show 𝟙 Y.left ≫ Y.hom = q
    rw [Category.id_comp, hq]
  exact Functor.IsEquivalence.mk

include F in
/-- ★★★★★★★**[FrdI] Proposition 1.6, (ii)** ——
**`𝒞' = 𝒞 ×_𝒟 𝒟'` は `Definition 1.3` の 21 条をすべて満たす**。

原文 (FrdI p.28):
> equivalences, the conditions of Definition 1.3 follow via a routine verification. Thus,

★原文の「routine verification」の中身がこれである。★**21 条のうち最後まで
残っていたのは `plBkEquiv` と、それが要求する pull-back の「𝒞' ⟹ 𝒞」向き
(`cfp_isPullBack_to`)であった。** -/
theorem cfpFrobenioidCore : FrobenioidCore (cfpPreFrobenioid P G hG hD' hcC hcD') where
  baseSurj := cfp_baseSurj P G hG hD' hcC hcD' F
  preStepSpan := cfp_preStepSpan P G hG hD' hcC hcD' F
  plBkEquiv := cfp_plBkEquiv P G hG hD' hcC hcD' F
  frobDegSurj := cfp_frobDegSurj P G hG hD' hcC hcD' F
  frobDegUniq := cfp_frobDegUniq P G hG hD' hcC hcD' F
  coAngularComp := cfp_coAngularComp P G hG hD' hcC hcD' F
  coAngularOfPreStep := fun α hca hst φ =>
    cfp_coAngularOfPreStep P G hG hD' hcC hcD' F α hca hst φ
  otriFwd := fun φ hca hst α hα => cfp_otriFwd P G hG hD' hcC hcD' F φ hca hst α hα
  otriBwd := fun φ hca hst β hβ => cfp_otriBwd P G hG hD' hcC hcD' F φ hca hst β hβ
  otriBase := fun φ φ' hca hst hca' hst' hbase α hα β hβ h =>
    cfp_otriBase P G hG hD' hcC hcD' F φ φ' hca hst hca' hst' hbase α hα β hβ h
  arbFactor := cfp_arbFactor P G hG hD' hcC hcD' F
  arbFactorUniq := fun X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα' =>
    cfp_arbFactorUniq P G hG hD' hcC hcD' F X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα'
  pullBackLB := cfp_pullBackLB P G hG hD' hcC hcD' F
  preStepMono := cfp_preStepMono P G hG hD' hcC hcD' F
  preStepFactor := cfp_preStepFactor P G hG hD' hcC hcD' F
  preStepFactorUniq := fun X X' β α β' α' heq hca hβ hαi hα hca' hβ' hαi' hα' =>
    cfp_preStepFactorUniq P G hG hD' hcC hcD' F X X' β α β' α' heq hca hβ hαi hα hca' hβ' hαi' hα'
  preStepFactor' := cfp_preStepFactor' P G hG hD' hcC hcD' F
  preStepFactorUniq' := fun X X' β α β' α' heq hβi hβ hca hα hβi' hβ' hca' hα' =>
    cfp_preStepFactorUniq' P G hG hD' hcC hcD' F X X' β α β' α' heq hβi hβ hca hα hβi' hβ' hca' hα'
  faithfulUpToUnits := fun φ ψ hb hm hφc hφs hψc hψs =>
    cfp_faithfulUpToUnits P G hG hD' hcC hcD' F φ ψ hb hm hφc hφs hψc hψs
  isotropicHullExists := cfp_isotropicHullExists P G hG hD' hcC hcD' F
  isotropicClosed := cfp_isotropicClosed P G hG hD' hcC hcD' F

/-! ### ★`plBkEquiv` —— ★解決済(2026-08-20)。以下は止まっていた頃の記録である

★★**数学**: 忠実性・充満性は **`𝒞'` の pull-back 性(対象の定義そのもの)を直接使う**。
充満性で要る「pull-back の左簡約」は **`isPullBack_of_comp_left`(仮定なし)** で証明済みなので
**循環しない**。本質的全射性は `𝒞` の `plBkEquiv` の本質的全射性 ＋ `cfp_isPullBack_of`
(構成の向き)だけで足りる。★**pull-back の「`𝒞' ⟹ 𝒞`」向きは要らない。**

★**Lean で止まる唯一の箇所**(本質的全射性の中):
```
hw  : e ≫ (G.map Y.hom ≫ inv A.obj.hom) = P.proj.map Z'.hom.hom
⊢    P.proj.map Z'.hom.hom ≫ A.obj.hom = e ≫ G.map Y.hom
```
`rw [← hw, Category.assoc, Category.assoc]` の後、目標は表示上
`e ≫ G.map Y.hom ≫ inv A.obj.hom ≫ A.obj.hom = e ≫ G.map Y.hom` になるが、
★**`rw [IsIso.inv_hom_id]` が「`inv ?f ≫ ?f` が見つからない」と言う。**

★試して駄目だったもの(6 通り): `simp` / `simp only [Category.assoc, hwA2]` /
`Category.assoc` を引数明示で 2 回 / `calc` で括弧を明示 /
`wA`(`hA.out` から取った逆射)版 / `@inv _ _ _ _ A.obj.hom hA` 版。
★**原因は特定できていない。** 分類表 #1 の一種と見ているが、
**「症状ではなく原因を特定する」を実行できていない**ので、そのまま記録する。

★★**親による原因の候補(2026-08-15。★これは答えではなく「判別法」である)**

第12段で我々は「もっともらしい原因を思いついても、それが原因であることは
別に確かめる必要がある」を学んだ(C1 が偽だった)。だから候補として書く。

**候補A —— `inv` のインスタンス引数違い**
  `inv f` は `[IsIso f]` を暗黙に取る。目標中の `inv A.obj.hom` が担いでいる
  インスタンスと、`rw [IsIso.inv_hom_id]` が単一化で合成するインスタンスが
  **別物**なら、項として異なるので構文的な `rw` は照合に失敗する。
  `A.obj.hom` の `IsIso` には `A.property`(`FullSubcategory` の `Prop` フィールド)由来と
  `haveI hA` 由来の 2 つがありうる。
  ★**第12段で偽と判定した C1 が、別の場所で真になっている**可能性である。

**候補B —— `rw [Category.assoc]` の適用位置**
  `rw` は**最初の出現だけ**を書き換える。2 回当てても狙った位置に行くとは限らない。

★**判別法**: `simp only [Category.assoc]` **だけ**を当てて右結合に正規化してから
`rw [IsIso.inv_hom_id]` を試す。
  - **通れば B**(位置の問題だった)
  - **なお落ちれば A**(インスタンスの問題)。そのとき `set_option pp.all true` で
    `inv` の第4引数を目視すれば確定する。

★**未試行の手が 2 つある**(子の 6 通りに含まれない):
  1. `simp only [Category.assoc, IsIso.inv_hom_id, Category.comp_id]`
     —— 子が試したのは `hwA2` 版で、`IsIso.inv_hom_id` を simp 補題として
     入れた版は試していない。
  2. ★**`slice_lhs 3 4 => rw [IsIso.inv_hom_id]`**
     —— `Mathlib/Tactic/CategoryTheory/Slice.lean` にある。
     ★**位置で指定するので、結合にも `rw` の適用順にも依存しない。**
     候補 A・B のどちらであっても B なら確実に抜ける。
-/

/-! ### ★(参考) `plBkEquiv` の構造

★★**数学は片付きました**: 忠実性と充満性は **`𝒞'` の pull-back 性を直接使う**だけでよく、
充満性で要る「pull-back の左簡約」は **`isPullBack_of_comp_left`(上)で証明済み** ——
`Definition 1.2, (ii)` の全単射条件だけから出るので**循環しません**。
本質的全射性も `cfp_isPullBack_of`(構成の向き)だけで足ります。

★Lean では本質的全射性の中の `wA ≫ A.obj.hom` の書き換えが通らず(分類表 #1 の一種)、
**規模を超えたのでここで切りました**。★**pull-back の「`𝒞'` ⟹ `𝒞`」向きは要らない**
という結論は変わりません。
-/

end Core

end Dict
end ABC3.Found.FrdI
