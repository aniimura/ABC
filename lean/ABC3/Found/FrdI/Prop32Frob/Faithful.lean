/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop32Frob.Transition

/-!
# Prop32Frob —— (iii)(c) の順逆・3 脚添字

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★26. (iii)(c) の順方向(根が等しい場合)

★★根が等しいと `compRoot α φ` と `compRoot φ β` が**同じ `rtRootIso`** を使う
(どちらも `r*r = r*r` の証明で、証明無関係により同一)。
★したがって両者は `compPf_mk` の 1 本で比較でき、
`𝒞` の `otriFwd` がそのまま渡る。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**(iii)(c) の順方向**(根が等しい場合)。 -/
theorem pfRoot_otriFwd_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ)
    (α : End (⟨A, r⟩ : PfRootObj P F)) (hα : α ∈ OTri (pfRootPre P F) ⟨A, r⟩) :
    ∃ β : End (⟨B, r⟩ : PfRootObj P F), β ∈ OTri (pfRootPre P F) ⟨B, r⟩ ∧
      (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
        = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
  obtain ⟨V, a₀, f₀, ha₀, hf₀⟩ := exists_rep3 (P := P) (F := F)
    ((rtRootIso P F A A (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv
      (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩))
    ((rtRootIso P F A B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).inv φ)
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag (F := F) hfi V
  obtain ⟨a₁, ha₁⟩ : ∃ t : Dd ⟶ Dd,
      t = idxTransport P F ((idx12 P F _ _ _).map u) a₀ := ⟨_, rfl⟩
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  have haW : (rtRootIso P F A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁ := by
    rw [ha₁]
    exact ha₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) a₀).symm
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀).symm
  have hαmk : (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) rfl) a₁ := by
    exact ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A A
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) haW)).trans
      (rtRootIso_hom_mk (F := F) A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ := by
    exact ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have ha₁O : P.Base a₁ = P.Base (𝟙 Dd) ∧ P.degFr a₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) a₁).mp ?_
    rw [← hαmk]; exact hα
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨b₁, ⟨hb₁O, hb₁e⟩, -⟩ :=
    F.otriFwd f₁ hf₁c hf₁s (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
  refine ⟨(rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    refine (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) _).mpr ⟨hb₁O.1, hb₁O.2⟩
  · show compRoot P F φ _ = compRoot P F (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) φ
    have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
      compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ _
    have e2 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
        (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
      compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁
    have e3 : HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃))
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
      congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e
    unfold compRoot
    rw [hfW, Iso.hom_inv_id_apply, haW, e1, e2, e3]

/-! ## ★27. 3 脚添字を「第 2・第 3 脚が一致し、第 1 脚の終域が isotropic」へ

★(iii)(c) の**逆方向**では `β` が第 2・第 3 脚にまたがる自己射なので、
今度はそちら側を対角にする。 -/

variable {P F} in
/-- ★★**3 脚添字を「第 2・第 3 脚が同じ射で、第 1 脚の終域が isotropic」へ**。 -/
theorem exists_idx3_diag23 (hfi : IsOfFrobeniusIsotropicType P) {A E : C}
    (V : IdxPf3 P F A E E) :
    ∃ (Dd E₃ : C) (l : A ⟶ Dd) (hl : IsFrobeniusType P l) (m : E ⟶ E₃)
      (hm : IsFrobeniusType P m) (hd : P.degFr l = P.degFr m),
      IsIsotropic P Dd ∧ Nonempty (V ⟶ idxMk3 (F := F) l m m hl hm hm hd rfl) := by
  obtain ⟨hv1, hv2, hv3, h12, h23⟩ := V.hom.property
  obtain ⟨X₃, s, t, hs, ht, hsd, htd, hst⟩ :=
    frob_common_upper P F V.hom.hom.2.1 hv2 V.hom.hom.2.2 hv3
  obtain ⟨A₁, a, ha, had⟩ := F.frobDegSurj V.right.obj.1 (P.degFr s)
  obtain ⟨Dd, δ, hδ, hDd⟩ := hfi A₁
  obtain ⟨E₃, ε, hε, hεd⟩ := F.frobDegSurj X₃ (P.degFr δ)
  have hstε : V.hom.hom.2.1 ≫ (s ≫ ε) = V.hom.hom.2.2 ≫ (t ≫ ε) := by
    rw [← Category.assoc, ← Category.assoc, hst]
    rfl
  have hsdt : P.degFr s = P.degFr t := by
    rw [hsd, htd, h23]
    rfl
  have hdeg1 : P.degFr (a ≫ δ) = P.degFr (s ≫ ε) := by
    rw [P.degFr_comp a δ, P.degFr_comp s ε, had, hεd]
  have hdeg2 : P.degFr (s ≫ ε) = P.degFr (t ≫ ε) := by
    rw [P.degFr_comp s ε, P.degFr_comp t ε, hsdt]
  have hlm : P.degFr (V.hom.hom.1 ≫ (a ≫ δ)) = P.degFr (V.hom.hom.2.1 ≫ (s ≫ ε)) := by
    rw [P.degFr_comp V.hom.hom.1 (a ≫ δ), P.degFr_comp V.hom.hom.2.1 (s ≫ ε),
      hdeg1, h12]
  refine ⟨Dd, E₃, V.hom.hom.1 ≫ (a ≫ δ),
    IsFrobeniusType.comp P F hv1 (IsFrobeniusType.comp P F ha hδ),
    V.hom.hom.2.1 ≫ (s ≫ ε), IsFrobeniusType.comp P F hv2
      (IsFrobeniusType.comp P F hs hε), hlm, hDd,
    ⟨Under.homMk (show V.right ⟶ (⟨(Dd, E₃, E₃)⟩ : TriFr P F) from
      ⟨(a ≫ δ, s ≫ ε, t ≫ ε), IsFrobeniusType.comp P F ha hδ,
        IsFrobeniusType.comp P F hs hε, IsFrobeniusType.comp P F ht hε,
        hdeg1, hdeg2⟩)
      (WideSubcategory.hom_ext _ (Prod.ext rfl (Prod.ext rfl hstε.symm)))⟩⟩

/-! ## ★28. (iii)(c) の逆方向(根が等しい場合) -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**(iii)(c) の逆方向**(根が等しい場合)。 -/
theorem pfRoot_otriBwd_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ)
    (β : End (⟨B, r⟩ : PfRootObj P F)) (hβ : β ∈ OTri (pfRootPre P F) ⟨B, r⟩) :
    ∃ α : End (⟨A, r⟩ : PfRootObj P F), α ∈ OTri (pfRootPre P F) ⟨A, r⟩ ∧
      (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
        = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
  obtain ⟨V, f₀, b₀, hf₀, hb₀⟩ := exists_rep3 (P := P) (F := F)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).inv (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag23 (F := F) hfi V
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx12 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  obtain ⟨b₁, hb₁⟩ : ∃ t : E₃ ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) b₀ := ⟨_, rfl⟩
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) f₀).symm
  have hbW : (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) b₁ := by
    rw [hb₁]
    exact hb₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) b₀).symm
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have hβmk : (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) rfl) b₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F B B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hbW)).trans
      (rtRootIso_hom_mk (F := F) B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) m m hm hm rfl) b₁)
  have hb₁O : P.Base b₁ = P.Base (𝟙 E₃) ∧ P.degFr b₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm) b₁).mp ?_
    rw [← hβmk]; exact hβ
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨a₁, ⟨ha₁O, ha₁e⟩, -⟩ :=
    F.otriBwd f₁ hf₁c hf₁s (b₁ : End E₃) ⟨hb₁O.1, hb₁O.2⟩
  refine ⟨(rtRootIso P F A A (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) ((a₁ : End Dd) : Dd ⟶ Dd)), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) _).mpr ⟨ha₁O.1, ha₁O.2⟩
  · have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) b₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (f₁ ≫ b₁) :=
      compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ b₁
    have e2 : compPf P F
        (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) ((a₁ : End Dd) : Dd ⟶ Dd))
        (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁) :=
      compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) _ f₁
    have e3 : HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (f₁ ≫ b₁)
        = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
          (((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁) :=
      congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) ha₁e
    show compRoot P F φ _ = compRoot P F _ φ
    unfold compRoot
    rw [hfW, hbW, Iso.hom_inv_id_apply, e1, e2, e3]

/-! ## ★29. 同じ添字の 2 射の底が一致すれば、代表元の底も一致する -/

variable {P F} in
/-- ★**同じ添字なら `rootBase` の一致から `Base` の一致が出る**。 -/
theorem base_eq_of_rootBase_eq {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (f f' : Z.right.obj.1 ⟶ Z.right.obj.2)
    (h : rootBase (show HomRoot P F X Y from HomPf.mk Z f)
      = rootBase (show HomRoot P F X Y from HomPf.mk Z f')) :
    P.Base f = P.Base f' := by
  haveI h1 : IsIso (P.Base Z.hom.hom.1) := Z.hom.property.1.2
  haveI he : IsIso (P.Base (rtExt P F X.obj Y.root)) :=
    (rtExt_frobType P F X.obj Y.root).2
  have s1 := rootBase_spec (show HomRoot P F X Y from HomPf.mk Z f)
  have s2 := rootBase_spec (show HomRoot P F X Y from HomPf.mk Z f')
  rw [pfBase_mk] at s1 s2
  have e1 : P.Base (rtExt P F X.obj Y.root) ≫ repBase Z f
      = P.Base (rtExt P F X.obj Y.root) ≫ repBase Z f' := by
    rw [← s1, ← s2, h]
  have e2 : repBase Z f = repBase Z f' :=
    (cancel_epi (P.Base (rtExt P F X.obj Y.root))).mp e1
  have r1 := repBase_spec Z f
  have r2 := repBase_spec Z f'
  rw [e2] at r1
  exact (cancel_epi (P.Base Z.hom.hom.1)).mp (r1.symm.trans r2)

/-! ## ★30. (iii)(c) の「底にしか依らない」条(根が等しい場合) -/

set_option maxHeartbeats 4000000 in
variable {P F} in
/-- ★★★**(iii)(c) の第 3 条**(根が等しい場合)。 -/
theorem pfRoot_otriBase_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ φ' : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hφs : IsPreStep (pfRootPre P F) φ) (hφ's : IsPreStep (pfRootPre P F) φ')
    (hbase : (pfRootPre P F).Base φ = (pfRootPre P F).Base φ')
    (α : End (⟨A, r⟩ : PfRootObj P F)) (hα : α ∈ OTri (pfRootPre P F) ⟨A, r⟩)
    (β : End (⟨B, r⟩ : PfRootObj P F)) (_hβ : β ∈ OTri (pfRootPre P F) ⟨B, r⟩)
    (heq : (φ ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ)) :
    (φ' ≫ (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ') := by
  obtain ⟨V, a₀, f₀, f₀', ha₀, hf₀, hf₀'⟩ := exists_rep3_pair (P := P) (F := F)
    ((rtRootIso P F A A (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv
      (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩))
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ')
  obtain ⟨Dd, E₃, l, hl, m, hm, hdlm, hDd, ⟨u⟩⟩ := exists_idx3_diag (F := F) hfi V
  obtain ⟨a₁, ha₁⟩ : ∃ t : Dd ⟶ Dd,
      t = idxTransport P F ((idx12 P F _ _ _).map u) a₀ := ⟨_, rfl⟩
  obtain ⟨f₁, hf₁⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀ := ⟨_, rfl⟩
  obtain ⟨f₁', hf₁'⟩ : ∃ t : Dd ⟶ E₃,
      t = idxTransport P F ((idx23 P F _ _ _).map u) f₀' := ⟨_, rfl⟩
  have haW : (rtRootIso P F A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁ := by
    rw [ha₁]
    exact ha₀.trans (HomPf.mk_map (P := P) (F := F) ((idx12 P F _ _ _).map u) a₀).symm
  have hfW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁ := by
    rw [hf₁]
    exact hf₀.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀).symm
  have hfW' : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ'
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁' := by
    rw [hf₁']
    exact hf₀'.trans (HomPf.mk_map (P := P) (F := F) ((idx23 P F _ _ _).map u) f₀').symm
  have hφmk : φ = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
  have hφ'mk : φ' = HomPf.mk (idxMk (P := P) (F := F)
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (rtLift P F B (show r * r = r * r from rfl) ≫ m)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
      (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁' :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hfW')).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
  have hαmk : (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩)
      = HomPf.mk (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) rfl) a₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A A
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) haW)).trans
      (rtRootIso_hom_mk (F := F) A A (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
  rw [hφmk, hφ'mk] at hbase
  have hbf : P.Base f₁ = P.Base f₁' :=
    base_eq_of_rootBase_eq (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm])) f₁ f₁' hbase
  have ha₁O : P.Base a₁ = P.Base (𝟙 Dd) ∧ P.degFr a₁ = 1 := by
    refine (oTri_mk_diag (X := (⟨A, r⟩ : PfRootObj P F))
      (rtLift P F A (show r * r = r * r from rfl) ≫ l)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl) a₁).mp ?_
    rw [← hαmk]; exact hα
  have hidx : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) :=
    idxMk (P := P) (F := F) l m hl hm hdlm
  have hf₁s : IsPreStep P f₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁).mp ?_
    rw [← hφmk]; exact hφs
  have hf₁'s : IsPreStep P f₁' := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F))
      (idxMk (P := P) (F := F)
        (rtLift P F A (show r * r = r * r from rfl) ≫ l)
        (rtLift P F B (show r * r = r * r from rfl) ≫ m)
        (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hm)
        (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr, hdlm]))
      f₁').mp ?_
    rw [← hφ'mk]; exact hφ's
  have hf₁c : IsCoAngular P f₁ :=
    prop_1_4_i P f₁ (fun _ g => F.isotropicClosed g hDd)
  have hf₁'c : IsCoAngular P f₁' :=
    prop_1_4_i P f₁' (fun _ g => F.isotropicClosed g hDd)
  obtain ⟨b₁, ⟨hb₁O, hb₁e⟩, -⟩ :=
    F.otriFwd f₁ hf₁c hf₁s (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
  have hb₁e' : (f₁' ≫ ((b₁ : End E₃) : E₃ ⟶ E₃))
      = ((a₁ : End Dd) : Dd ⟶ Dd) ≫ f₁' :=
    F.otriBase f₁ f₁' hf₁c hf₁s hf₁'c hf₁'s hbf (a₁ : End Dd) ⟨ha₁O.1, ha₁O.2⟩
      (b₁ : End E₃) hb₁O hb₁e
  have e1 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁ ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁ _
  have e1' : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃))
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)
        (f₁' ≫ ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    compPf_mk (idxMk3 (F := F) l m m hl hm hm hdlm rfl) f₁' _
  have e2 : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
      (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁)
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁) :=
    compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁
  have e2' : compPf P F (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a₁)
      (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) f₁')
      = HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm) (a₁ ≫ f₁') :=
    compPf_mk (idxMk3 (F := F) l l m hl hl hm rfl hdlm) a₁ f₁'
  have heq0 : (φ ≫ (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ) := by
    show compRoot P F φ _ = compRoot P F _ φ
    unfold compRoot
    rw [hfW, Iso.hom_inv_id_apply, haW, e1, e2]
    exact congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom)
      (congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e)
  have heq0' : (φ' ≫ (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)))
      = ((α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨A, r⟩) ≫ φ') := by
    show compRoot P F φ' _ = compRoot P F _ φ'
    unfold compRoot
    rw [hfW', Iso.hom_inv_id_apply, haW, e1', e2']
    exact congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom)
      (congrArg (HomPf.mk (idxMk (P := P) (F := F) l m hl hm hdlm)) hb₁e')
  haveI : Epi φ := pfRoot_totEpi P F _ _ φ
  have hββ : (β : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
      = (rtRootIso P F B B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) m m hm hm rfl) ((b₁ : End E₃) : E₃ ⟶ E₃)) :=
    (cancel_epi φ).mp (heq.trans heq0.symm)
  rw [hββ]
  exact heq0'

/-! ## ★31. 根を揃えて (iii)(c) を一般の場合へ

★★`X = (A,n)`、`Y = (B,m)` を `X' = (A^m, n*m)`、`Y' = (B^n, n*m)` に取り替える。
★どちらも根が `n*m` なので `sameRoot` 版が当たり、結果は共役で戻す。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の順方向、`𝒞^pf` 版。 -/
theorem pfRoot_otriFwd (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (α : End X) (hα : α ∈ OTri (pfRootPre P F) X) :
    ∃! β : End Y, β ∈ OTri (pfRootPre P F) Y ∧
      ((φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hα' := oTri_conj (F := F) eX (inv eX) (IsIso.hom_inv_id eX) (IsIso.inv_hom_id eX)
    α hα
  obtain ⟨β', hβ'O, hβ'e⟩ := pfRoot_otriFwd_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY) hφ'
    (inv eX ≫ (α : X ⟶ X) ≫ eX) hα'
  have hβO := oTri_conj (F := F) (inv eY) eY (IsIso.inv_hom_id eY) (IsIso.hom_inv_id eY)
    β' hβ'O
  have h1 : inv eX ≫ (φ ≫ eY ≫ (β' : _ ⟶ _))
      = inv eX ≫ ((α : X ⟶ X) ≫ φ ≫ eY) := by
    have h := hβ'e
    simp only [Category.assoc] at h ⊢
    rw [h]
    simp only [Category.assoc, IsIso.hom_inv_id_assoc]
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  have h2 : φ ≫ eY ≫ (β' : _ ⟶ _) = (α : X ⟶ X) ≫ φ ≫ eY :=
    (cancel_epi (inv eX)).mp h1
  have hβe : (φ ≫ (eY ≫ (β' : _ ⟶ _) ≫ inv eY)) = ((α : X ⟶ X) ≫ φ) := by
    calc φ ≫ (eY ≫ (β' : _ ⟶ _) ≫ inv eY)
        = (φ ≫ eY ≫ (β' : _ ⟶ _)) ≫ inv eY := by simp only [Category.assoc]
      _ = ((α : X ⟶ X) ≫ φ ≫ eY) ≫ inv eY := by rw [h2]
      _ = (α : X ⟶ X) ≫ φ := by
          simp only [Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  refine ⟨eY ≫ (β' : _ ⟶ _) ≫ inv eY, ⟨hβO, hβe⟩, ?_⟩
  rintro y ⟨-, hy⟩
  haveI : Epi φ := pfRoot_totEpi P F _ _ φ
  exact (cancel_epi φ).mp (hy.trans hβe.symm)

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の逆方向、`𝒞^pf` 版。 -/
theorem pfRoot_otriBwd (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (β : End Y) (hβ : β ∈ OTri (pfRootPre P F) Y) :
    ∃! α : End X, α ∈ OTri (pfRootPre P F) X ∧
      ((φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hβ' := oTri_conj (F := F) eY (inv eY) (IsIso.hom_inv_id eY) (IsIso.inv_hom_id eY)
    β hβ
  obtain ⟨α', hα'O, hα'e⟩ := pfRoot_otriBwd_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY) hφ'
    (inv eY ≫ (β : Y ⟶ Y) ≫ eY) hβ'
  have hαO := oTri_conj (F := F) (inv eX) eX (IsIso.inv_hom_id eX) (IsIso.hom_inv_id eX)
    α' hα'O
  have h1 : inv eX ≫ (φ ≫ (β : Y ⟶ Y) ≫ eY)
      = inv eX ≫ ((eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ ≫ eY) := by
    have h := hα'e
    simp only [Category.assoc, IsIso.hom_inv_id_assoc] at h
    simp only [Category.assoc, IsIso.inv_hom_id_assoc]
    exact h
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  have h2 : φ ≫ (β : Y ⟶ Y) ≫ eY = (eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ ≫ eY :=
    (cancel_epi (inv eX)).mp h1
  haveI : Mono eY := IsIso.mono_of_iso _
  have hαe : (φ ≫ (β : Y ⟶ Y)) = ((eX ≫ (α' : _ ⟶ _) ≫ inv eX) ≫ φ) := by
    refine (cancel_mono eY).mp ?_
    simp only [Category.assoc]
    simp only [Category.assoc] at h2
    exact h2
  refine ⟨eX ≫ (α' : _ ⟶ _) ≫ inv eX, ⟨hαO, hαe⟩, ?_⟩
  rintro y ⟨-, hy⟩
  haveI : Mono φ := pfRoot_preStepMono (F := F) φ hφs
  exact (cancel_mono φ).mp (hy.symm.trans hαe)

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iii)(c)** の第 3 条、`𝒞^pf` 版。 -/
theorem pfRoot_otriBase (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (φ φ' : X ⟶ Y) (_hφc : IsCoAngular (pfRootPre P F) φ)
    (hφs : IsPreStep (pfRootPre P F) φ) (_hφ'c : IsCoAngular (pfRootPre P F) φ')
    (hφ's : IsPreStep (pfRootPre P F) φ')
    (hbase : (pfRootPre P F).Base φ = (pfRootPre P F).Base φ')
    (α : End X) (hα : α ∈ OTri (pfRootPre P F) X)
    (β : End Y) (hβ : β ∈ OTri (pfRootPre P F) Y)
    (heq : (φ ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ)) :
    (φ' ≫ (β : Y ⟶ Y)) = ((α : X ⟶ X) ≫ φ') := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hψ : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hψ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ' ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφ's (isPreStep_of_isIso (pfRootPre P F) eY))
  have hbase' : (pfRootPre P F).Base (inv eX ≫ φ ≫ eY)
      = (pfRootPre P F).Base (inv eX ≫ φ' ≫ eY) := by
    rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp,
      (pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp, hbase]
  have hα' := oTri_conj (F := F) eX (inv eX) (IsIso.hom_inv_id eX) (IsIso.inv_hom_id eX)
    α hα
  have hβ' := oTri_conj (F := F) eY (inv eY) (IsIso.hom_inv_id eY) (IsIso.inv_hom_id eY)
    β hβ
  have heq' : ((inv eX ≫ φ ≫ eY) ≫ (inv eY ≫ (β : Y ⟶ Y) ≫ eY))
      = ((inv eX ≫ (α : X ⟶ X) ≫ eX) ≫ (inv eX ≫ φ ≫ eY)) := by
    simp only [Category.assoc, IsIso.hom_inv_id_assoc, IsIso.inv_hom_id_assoc]
    have h3 : (φ ≫ (β : Y ⟶ Y)) ≫ eY = ((α : X ⟶ X) ≫ φ) ≫ eY :=
      congrArg (fun t => t ≫ eY) heq
    simp only [Category.assoc] at h3
    exact congrArg (fun t => inv eX ≫ t) h3
  have hres := pfRoot_otriBase_sameRoot (F := F) hfi (inv eX ≫ φ ≫ eY)
    (inv eX ≫ φ' ≫ eY) hψ hψ' hbase' (inv eX ≫ (α : X ⟶ X) ≫ eX) hα'
    (inv eY ≫ (β : Y ⟶ Y) ≫ eY) hβ' heq'
  simp only [Category.assoc, IsIso.hom_inv_id_assoc, IsIso.inv_hom_id_assoc] at hres
  haveI : Epi (inv eX) := IsIso.epi_of_iso _
  haveI : Mono eY := IsIso.mono_of_iso _
  refine (cancel_mono eY).mp ?_
  refine (cancel_epi (inv eX)).mp ?_
  simp only [Category.assoc]
  exact hres

end ABC3.Found.FrdI
