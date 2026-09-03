/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop32Frob.Pullback

/-!
# Prop32Frob —— `𝒞` の pull-back・底の比較・出典

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★47. ★★★★`𝒞` の pull-back は `𝒞^pf` でも pull-back(一意性の側) -/

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★普遍性の**一意性**の側。 -/
theorem pfRoot_isPullBack_mk_hext (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (α : X ⟶ Y)
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (χ : Z.right.obj.1 ⟶ Z.right.obj.2) (hαdef : α = HomPf.mk Z χ)
    (hχ : IsPullBack P χ)
    (T : PfRootObj P F) (u v : T ⟶ X) (h1 : u ≫ α = v ≫ α)
    (h2 : (pfRootPre P F).Base u = (pfRootPre P F).Base v) : u = v := by
  obtain ⟨V₀, φ₀, φ₀', ψ₀, hu₀, hv₀, hα₀, hc₀, hd₀⟩ :=
    compRoot_rep_pairL (F := F) u v α
  obtain ⟨V, tr, hisoV⟩ := exists_idx3_isotropic2 (F := F) hfi V₀
  obtain ⟨φ, hφ⟩ : ∃ t : V.right.obj.1 ⟶ V.right.obj.2.1,
      t = idxTransport P F ((idx12 P F _ _ _).map tr) φ₀ := ⟨_, rfl⟩
  obtain ⟨φ', hφ'⟩ : ∃ t : V.right.obj.1 ⟶ V.right.obj.2.1,
      t = idxTransport P F ((idx12 P F _ _ _).map tr) φ₀' := ⟨_, rfl⟩
  obtain ⟨ψ, hψ⟩ : ∃ t : V.right.obj.2.1 ⟶ V.right.obj.2.2,
      t = idxTransport P F ((idx23 P F _ _ _).map tr) ψ₀ := ⟨_, rfl⟩
  have hu : u = (rtRootIso P F T.obj X.obj
      (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * T.root = Y.root * T.root from rfl)).hom
      (HomPf.mk ((idx12 P F _ _ _).obj V) φ) := by
    rw [hφ, HomPf.mk_map]; exact hu₀
  have hv : v = (rtRootIso P F T.obj X.obj
      (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * T.root = Y.root * T.root from rfl)).hom
      (HomPf.mk ((idx12 P F _ _ _).obj V) φ') := by
    rw [hφ', HomPf.mk_map]; exact hv₀
  have hα : α = (rtRootIso P F X.obj Y.obj
      (show Y.root * T.root = T.root * Y.root from mul_comm _ _)
      (show X.root * T.root = T.root * X.root from mul_comm _ _)).hom
      (HomPf.mk ((idx23 P F _ _ _).obj V) ψ) := by
    rw [hψ, HomPf.mk_map]; exact hα₀
  have hc : u ≫ α = (rtRootIso P F T.obj Y.obj
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
      (show X.root * T.root = X.root * T.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)) := by
    refine Eq.trans hc₀ (congrArg _ ?_)
    rw [hφ, hψ]
    exact (HomPf.mk_map ((idx13 P F _ _ _).map tr) (φ₀ ≫ ψ₀)).symm.trans
      (congrArg (HomPf.mk _) (idxTransport_comp_pair (F := F) tr φ₀ ψ₀).symm)
  have hd : v ≫ α = (rtRootIso P F T.obj Y.obj
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
      (show X.root * T.root = X.root * T.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ)) := by
    refine Eq.trans hd₀ (congrArg _ ?_)
    rw [hφ', hψ]
    exact (HomPf.mk_map ((idx13 P F _ _ _).map tr) (φ₀' ≫ ψ₀)).symm.trans
      (congrArg (HomPf.mk _) (idxTransport_comp_pair (F := F) tr φ₀' ψ₀).symm)
  have hψco : IsCoAngular P ψ := prop_1_4_i P ψ (fun _ g => F.isotropicClosed g hisoV)
  have heqrep := hαdef.symm.trans hα
  rw [rtRootIso_hom_mk] at heqrep
  have hψpb : IsPullBack P ψ := pfRoot_rep_isPullBack (F := F) hχ heqrep hψco
  have hcomp : HomPf.mk ((idx13 P F _ _ _).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F _ _ _).obj V) (φ' ≫ ψ) := by
    have h := hc.symm.trans (h1.trans hd)
    have h' := congrArg (fun t => (rtRootIso P F T.obj Y.obj
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
      (show X.root * T.root = X.root * T.root from rfl)).inv t) h
    rwa [Iso.hom_inv_id_apply, Iso.hom_inv_id_apply] at h'
  have hbase : P.Base φ = P.Base φ' := by
    have hu' := hu
    have hv' := hv
    rw [rtRootIso_hom_mk] at hu' hv'
    rw [hu', hv'] at h2
    exact base_eq_of_rootBase_eq (X := T) (Y := X)
      ((pushIdx (F := F)
        (rtLift P F T.obj (show Y.root * X.root = Y.root * X.root from rfl))
        (rtLift_frobType P F T.obj _)
        (rtLift P F X.obj (show Y.root * T.root = Y.root * T.root from rfl))
        (rtLift_frobType P F X.obj _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx12 P F _ _ _).obj V)) φ φ' h2
  have hfin := homPf_cancel_pullBack (F := F) V φ φ' ψ hψpb hbase hcomp
  rw [hu, hv]
  exact congrArg _ hfin

/-! ## ★48. 底の比較 —— 3 つの位置で**同じ同型**が現れる

★★`compRoot` の 3 つの `rtRootIso` はそれぞれ違う持ち上げを使うが、
`rtLift_ext`(`rtExt d ≫ rtLift = rtExt t`)で潰すと
**`Base(rtExt · t) ≫ Base(脚)` という同じ形**になる。★これが底の追跡を可能にする。 -/

variable {P F} in
/-- ★`rootBase_mk_spec` の `idxMk` 版。 -/
theorem rootBase_mk_spec' {X Y : PfRootObj P F} {A' B' : C}
    (a : rtObj P F X.obj Y.root ⟶ A') (b : rtObj P F Y.obj X.root ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A' ⟶ B') :
    rootBase (show HomRoot P F X Y from
        HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) φ)
        ≫ (P.Base (rtExt P F Y.obj X.root) ≫ P.Base b)
      = (P.Base (rtExt P F X.obj Y.root) ≫ P.Base a) ≫ P.Base φ :=
  rootBase_mk_spec (idxMk (P := P) (F := F) a b ha hb hd) φ

variable {P F} in
/-- ★★**根の取り替えを通した底の特徴づけ**。 -/
theorem rootBase_rtRootIso_mk_spec {X Y : PfRootObj P F} {e tA tB : ℕ+}
    (hA : tA = e * Y.root) (hB : tB = e * X.root)
    (W : IdxPf P F (rtObj P F X.obj tA) (rtObj P F Y.obj tB))
    (φ : W.right.obj.1 ⟶ W.right.obj.2) :
    rootBase (show HomRoot P F X Y from
        (rtRootIso P F X.obj Y.obj hA hB).hom (HomPf.mk W φ))
        ≫ (P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2)
      = (P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1) ≫ P.Base φ := by
  rw [rtRootIso_hom_mk]
  have hL : P.Base (rtExt P F Y.obj X.root) ≫ P.Base (rtLift P F Y.obj hB ≫ W.hom.hom.2)
      = P.Base (rtExt P F Y.obj tB) ≫ P.Base W.hom.hom.2 := by
    rw [P.Base_comp, ← Category.assoc, ← P.Base_comp, rtLift_ext]
  have hR : P.Base (rtExt P F X.obj Y.root) ≫ P.Base (rtLift P F X.obj hA ≫ W.hom.hom.1)
      = P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1 := by
    rw [P.Base_comp, ← Category.assoc, ← P.Base_comp, rtLift_ext]
  have h := rootBase_mk_spec' (X := X) (Y := Y)
    (rtLift P F X.obj hA ≫ W.hom.hom.1) (rtLift P F Y.obj hB ≫ W.hom.hom.2)
    (IsFrobeniusType.comp P F (rtLift_frobType P F X.obj hA) W.hom.property.1)
    (IsFrobeniusType.comp P F (rtLift_frobType P F Y.obj hB) W.hom.property.2.1)
    (by rw [P.degFr_comp, P.degFr_comp, rtLift_degFr, rtLift_degFr]
        exact congrArg (fun m => m * e) W.hom.property.2.2)
    φ
  rw [hL, hR] at h
  exact h

/-! ## ★49. ★★★★`𝒞` の pull-back は `𝒞^pf` でも pull-back(持ち上げの側) -/

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★普遍性の**持ち上げ**の側。 -/
theorem pfRoot_isPullBack_mk_lift (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (α : X ⟶ Y)
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (χ : Z.right.obj.1 ⟶ Z.right.obj.2) (hαdef : α = HomPf.mk Z χ)
    (hχ : IsPullBack P χ)
    (T : PfRootObj P F) (f : T ⟶ Y)
    (b : ((pfRootPre P F).toElem.obj T).base ⟶ ((pfRootPre P F).toElem.obj X).base)
    (hb : (pfRootPre P F).Base f = b ≫ (pfRootPre P F).Base α) :
    ∃ g : T ⟶ X, g ≫ α = f ∧ (pfRootPre P F).Base g = b := by
  obtain ⟨V₀, θ₀, ψ₀, hf₀, hα₀⟩ := exists_rep3_cospan (P := P) (F := F)
    ((rtRootIso P F T.obj Y.obj
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
      (show X.root * T.root = X.root * T.root from rfl)).inv f)
    ((rtRootIso P F X.obj Y.obj
      (show Y.root * T.root = T.root * Y.root from mul_comm _ _)
      (show X.root * T.root = T.root * X.root from mul_comm _ _)).inv α)
  obtain ⟨V, tr, hisoV⟩ := exists_idx3_isotropic2 (F := F) hfi V₀
  obtain ⟨hV1, hV2, hV3, hV12, hV23⟩ := V.hom.property
  obtain ⟨θ, hθ⟩ : ∃ t : V.right.obj.1 ⟶ V.right.obj.2.2,
      t = idxTransport P F ((idx13 P F _ _ _).map tr) θ₀ := ⟨_, rfl⟩
  obtain ⟨ψ, hψ⟩ : ∃ t : V.right.obj.2.1 ⟶ V.right.obj.2.2,
      t = idxTransport P F ((idx23 P F _ _ _).map tr) ψ₀ := ⟨_, rfl⟩
  have hf : f = (rtRootIso P F T.obj Y.obj
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)
      (show X.root * T.root = X.root * T.root from rfl)).hom
      (HomPf.mk ((idx13 P F _ _ _).obj V) θ) := by
    rw [hθ, HomPf.mk_map, ← hf₀, Iso.inv_hom_id_apply]
  have hα : α = (rtRootIso P F X.obj Y.obj
      (show Y.root * T.root = T.root * Y.root from mul_comm _ _)
      (show X.root * T.root = T.root * X.root from mul_comm _ _)).hom
      (HomPf.mk ((idx23 P F _ _ _).obj V) ψ) := by
    rw [hψ, HomPf.mk_map, ← hα₀, Iso.inv_hom_id_apply]
  -- ★`ψ` は pull-back
  have hψco : IsCoAngular P ψ := prop_1_4_i P ψ (fun _ g => F.isotropicClosed g hisoV)
  have heqrep := hαdef.symm.trans hα
  rw [rtRootIso_hom_mk] at heqrep
  have hψpb : IsPullBack P ψ := pfRoot_rep_isPullBack (F := F) hχ heqrep hψco
  -- ★底の 3 本
  have hfb : rootBase f
      ≫ (P.Base (rtExt P F Y.obj (X.root * T.root)) ≫ P.Base V.hom.hom.2.2)
      = (P.Base (rtExt P F T.obj (Y.root * X.root)) ≫ P.Base V.hom.hom.1) ≫ P.Base θ := by
    rw [hf]
    exact rootBase_rtRootIso_mk_spec (X := T) (Y := Y) _ _ ((idx13 P F _ _ _).obj V) θ
  have hαb : rootBase α
      ≫ (P.Base (rtExt P F Y.obj (X.root * T.root)) ≫ P.Base V.hom.hom.2.2)
      = (P.Base (rtExt P F X.obj (Y.root * T.root)) ≫ P.Base V.hom.hom.2.1) ≫ P.Base ψ := by
    rw [hα]
    exact rootBase_rtRootIso_mk_spec (X := X) (Y := Y) _ _ ((idx23 P F _ _ _).obj V) ψ
  haveI hiT : IsIso (P.Base (rtExt P F T.obj (Y.root * X.root))) :=
    (rtExt_frobType P F T.obj (Y.root * X.root)).2
  haveI hiX : IsIso (P.Base (rtExt P F X.obj (Y.root * T.root))) :=
    (rtExt_frobType P F X.obj (Y.root * T.root)).2
  haveI hiY : IsIso (P.Base (rtExt P F Y.obj (X.root * T.root))) :=
    (rtExt_frobType P F Y.obj (X.root * T.root)).2
  haveI hi1 : IsIso (P.Base V.hom.hom.1) := hV1.2
  haveI hi2 : IsIso (P.Base V.hom.hom.2.1) := hV2.2
  haveI hi3 : IsIso (P.Base V.hom.hom.2.2) := hV3.2
  haveI hiuT : IsIso (P.Base (rtExt P F T.obj (Y.root * X.root)) ≫ P.Base V.hom.hom.1) :=
    IsIso.comp_isIso' hiT hi1
  haveI hiuX : IsIso (P.Base (rtExt P F X.obj (Y.root * T.root)) ≫ P.Base V.hom.hom.2.1) :=
    IsIso.comp_isIso' hiX hi2
  -- ★底の持ち上げ条件
  obtain ⟨bb, hbb⟩ : ∃ t : (P.toElem.obj V.right.obj.1).base
        ⟶ (P.toElem.obj V.right.obj.2.1).base,
      t = inv (P.Base (rtExt P F T.obj (Y.root * X.root)) ≫ P.Base V.hom.hom.1)
        ≫ (b ≫ (P.Base (rtExt P F X.obj (Y.root * T.root)) ≫ P.Base V.hom.hom.2.1)) :=
    ⟨_, rfl⟩
  have hbθ : P.Base θ = bb ≫ P.Base ψ := by
    have e1 : (P.Base (rtExt P F T.obj (Y.root * X.root)) ≫ P.Base V.hom.hom.1)
          ≫ P.Base θ
        = (b ≫ (P.Base (rtExt P F X.obj (Y.root * T.root)) ≫ P.Base V.hom.hom.2.1))
          ≫ P.Base ψ :=
      hfb.symm.trans (((congrArg (fun t => t
          ≫ (P.Base (rtExt P F Y.obj (X.root * T.root)) ≫ P.Base V.hom.hom.2.2)) hb).trans
        ((Category.assoc _ _ _).trans (congrArg (fun t => b ≫ t) hαb))).trans
        (Category.assoc _ _ _).symm)
    rw [hbb]
    exact (eq_inv_comp_of _ hiuT _ _ e1).trans
      (Category.assoc (inv (P.Base (rtExt P F T.obj (Y.root * X.root))
        ≫ P.Base V.hom.hom.1))
        (b ≫ (P.Base (rtExt P F X.obj (Y.root * T.root)) ≫ P.Base V.hom.hom.2.1))
        (P.Base ψ)).symm
  obtain ⟨φ, ⟨hφ1, hφ2⟩, -⟩ := IsPullBack.lift P hψpb V.right.obj.1 θ bb hbθ
  refine ⟨(rtRootIso P F T.obj X.obj
      (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * T.root = Y.root * T.root from rfl)).hom
      (HomPf.mk ((idx12 P F _ _ _).obj V) φ), ?_, ?_⟩
  · rw [hα, hf]
    refine Eq.trans (compRoot_mk3 (X := T) (Y := X) (Z := Y)
      V.hom.hom.1 V.hom.hom.2.1 V.hom.hom.2.2 hV1 hV2 hV3 hV12 hV23 φ ψ) ?_
    exact congrArg _ (congrArg (HomPf.mk _) hφ1)
  · have hgb := rootBase_rtRootIso_mk_spec (X := T) (Y := X)
      (show Y.root * X.root = Y.root * X.root from rfl)
      (show Y.root * T.root = Y.root * T.root from rfl) ((idx12 P F _ _ _).obj V) φ
    refine (cancel_mono (P.Base (rtExt P F X.obj (Y.root * T.root))
      ≫ P.Base V.hom.hom.2.1)).mp ?_
    refine Eq.trans hgb ?_
    refine Eq.trans (congrArg (fun t => (P.Base (rtExt P F T.obj (Y.root * X.root))
      ≫ P.Base V.hom.hom.1) ≫ t) hφ2) ?_
    rw [hbb]
    exact IsIso.hom_inv_id_assoc _ _

/-! ## ★50. ★★★★**`𝒞` の pull-back は `𝒞^pf` でも pull-back** -/

variable {P F} in
/-- ★★★★**代表元が `𝒞` の pull-back なら `𝒞^pf` でも pull-back**。 -/
theorem pfRoot_isPullBack_mk (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (α : X ⟶ Y)
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (χ : Z.right.obj.1 ⟶ Z.right.obj.2) (hαdef : α = HomPf.mk Z χ)
    (hχ : IsPullBack P χ) : IsPullBack (pfRootPre P F) α :=
  isPullBack_of_lift (pfRootPre P F) α
    (fun T u v h1 h2 => pfRoot_isPullBack_mk_hext hfi α Z χ hαdef hχ T u v h1 h2)
    (fun T f b hb => pfRoot_isPullBack_mk_lift hfi α Z χ hαdef hχ T f b hb)

variable {P F} in
/-- ★`rtMap` は pull-back を保つ。 -/
theorem rtMap_isPullBack (k : ℕ+) {A B : C} (φ : A ⟶ B) (hφ : IsPullBack P φ) :
    IsPullBack P (rtMap (F := F) k φ) :=
  prop_1_10_i_pullBack_of P F (rtExt_frobType P F A k) (rtExt_frobType P F B k)
    (by rw [rtExt_degFr, rtExt_degFr]) (rtMap_spec (F := F) k φ) hφ

/-! ## ★51. ★★★★(iv)(a) —— 任意射の 3 分解 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★**Frobenius 型 ≫ 次数 1** への分解(根が等しい場合)。 -/
theorem pfRoot_frobSplit_sameRoot (hfi : IsOfFrobeniusIsotropicType P) {A B : C} {r : ℕ+}
    (f : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) :
    ∃ (G : C) (Gam : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨G, r⟩)
      (R : (⟨G, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩),
      IsFrobeniusType (pfRootPre P F) Gam ∧ IsLinear (pfRootPre P F) R ∧ f = Gam ≫ R := by
  obtain ⟨x, hx⟩ : ∃ x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)),
      (rtJ P F A B r).hom x = f := ⟨(rtJ P F A B r).inv f, Iso.inv_hom_id_apply _ _⟩
  obtain ⟨n, hn⟩ : ∃ n : ℕ+, pfDeg x = n := ⟨_, rfl⟩
  obtain ⟨z, hz1, hz2⟩ := homPf_frobSplit (P := P) (F := F) x hn
  have h1 : IsFrobeniusType P
      (rtExt P F A n ≫ rtExt P F (rtObj P F A n) (r * r)) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A n)
      (rtExt_frobType P F (rtObj P F A n) (r * r))
  have h2 : IsFrobeniusType P
      (rtExt P F A (r * r) ≫ rtExt P F (rtObj P F A (r * r)) n) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A (r * r))
      (rtExt_frobType P F (rtObj P F A (r * r)) n)
  have hdeg : P.degFr (rtExt P F A n ≫ rtExt P F (rtObj P F A n) (r * r))
      = P.degFr (rtExt P F A (r * r) ≫ rtExt P F (rtObj P F A (r * r)) n) := by
    rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr, rtExt_degFr, rtExt_degFr,
      mul_comm]
  obtain ⟨w, hw, -⟩ := F.frobDegUniq A _ _ _ _ h1 h2 hdeg
  haveI := hw
  obtain ⟨gam, hgam⟩ : ∃ t : rtObj P F A (r * r) ⟶ rtObj P F (rtObj P F A n) (r * r),
      t = rtExt P F (rtObj P F A (r * r)) n ≫ @inv _ _ _ _ w hw := ⟨_, rfl⟩
  have hgamF : IsFrobeniusType P gam := by
    haveI : IsIso (@inv _ _ _ _ w hw) := @IsIso.inv_isIso _ _ _ _ w hw
    rw [hgam]
    exact IsFrobeniusType.comp P F (rtExt_frobType P F (rtObj P F A (r * r)) n)
      (isFrobeniusType_of_isIso P (@inv _ _ _ _ w hw))
  refine ⟨rtObj P F A n,
    (rtJ P F A (rtObj P F A n) r).hom (toHomPf (F := F) gam),
    (rtJ P F (rtObj P F A n) B r).hom (compPf P F (toHomPf (F := F) w) z), ?_, ?_, ?_⟩
  · refine ⟨⟨pfRoot_isCoAngular hfi _, ?_⟩, ?_⟩
    · exact (rtJ_isIsometric_iff (F := F) (A := A) (B := rtObj P F A n) r
        (idxOne P F (rtObj P F A (r * r)) (rtObj P F (rtObj P F A n) (r * r))) gam).mpr
        hgamF.1.2
    · exact (rtJ_isBaseIso_iff (F := F) (A := A) (B := rtObj P F A n) r
        (idxOne P F (rtObj P F A (r * r)) (rtObj P F (rtObj P F A n) (r * r))) gam).mpr
        hgamF.2
  · show (pfRootPre P F).degFr _ = 1
    rw [rtJ_degFr', pfDeg_comp, hz1, pfDeg_toHomPf,
      show P.degFr w = 1 from isLinear_of_isIso P w, mul_one]
  · rw [rtJ_comp, ← compPf_assoc, ← toHomPf_comp, hgam, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id, hz2, hx]

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iv)(a)** の `𝒞^pf` 版(根が等しい場合)。 -/
theorem pfRoot_arbFactor_sameRoot (hfi : IsOfFrobeniusIsotropicType P) {A B : C} {r : ℕ+}
    (f : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) :
    ∃ (W₁ W₂ : PfRootObj P F) (γ : (⟨A, r⟩ : PfRootObj P F) ⟶ W₁) (β : W₁ ⟶ W₂)
      (a : W₂ ⟶ (⟨B, r⟩ : PfRootObj P F)),
      f = γ ≫ β ≫ a ∧ IsFrobeniusType (pfRootPre P F) γ ∧
        IsPreStep (pfRootPre P F) β ∧ IsPullBack (pfRootPre P F) a := by
  obtain ⟨G, Gam, R, hGam, hRlin, hfGR⟩ := pfRoot_frobSplit_sameRoot (F := F) hfi f
  obtain ⟨T, pp, ll, hpp, -, hReq, Z', χ', hlldef, hχ'⟩ :=
    pfRoot_linearSplit_sameRoot (F := F) R hRlin
  exact ⟨⟨G, r⟩, T, Gam, pp, ll, hfGR.trans (congrArg (fun t => Gam ≫ t) hReq),
    hGam, hpp, pfRoot_isPullBack_mk hfi ll Z' χ' hlldef hχ'⟩

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iv)(a)** —— `𝒞^pf` 版。 -/
theorem pfRoot_arbFactor (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) :
    ∃ (W₁ W₂ : PfRootObj P F) (γ : X ⟶ W₁) (β : W₁ ⟶ W₂) (a : W₂ ⟶ Y),
      f = γ ≫ β ≫ a ∧ IsFrobeniusType (pfRootPre P F) γ ∧
        IsPreStep (pfRootPre P F) β ∧ IsPullBack (pfRootPre P F) a := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  obtain ⟨W₁, W₂, γ, β, a, hfac, hγ, hβ, ha⟩ :=
    pfRoot_arbFactor_sameRoot (F := F) hfi (inv eX ≫ f ≫ eY)
  refine ⟨W₁, W₂, eX ≫ γ, β, a ≫ inv eY, ?_, ?_, hβ, ?_⟩
  · have h0 : eX ≫ (inv eX ≫ f ≫ eY) ≫ inv eY = f := by
      rw [Category.assoc, IsIso.hom_inv_id_assoc, Category.assoc, IsIso.hom_inv_id,
        Category.comp_id]
    conv_lhs => rw [← h0]
    rw [hfac]
    simp only [Category.assoc]
  · refine ⟨⟨pfRoot_isCoAngular hfi _, ?_⟩, ?_⟩
    · show (pfRootPre P F).Div (eX ≫ γ) = 0
      rw [(pfRootPre P F).Div_comp,
        show (pfRootPre P F).Div γ = 0 from hγ.1.2,
        show (pfRootPre P F).Div eX = 0 from isIsometric_of_isIso (pfRootPre P F) eX,
        map_zero, smul_zero, add_zero]
    · haveI : IsIso ((pfRootPre P F).Base γ) := hγ.2
      haveI : IsIso ((pfRootPre P F).Base eX) :=
        isBaseIsomorphism_of_isIso (pfRootPre P F) eX
      show IsIso ((pfRootPre P F).Base (eX ≫ γ))
      rw [(pfRootPre P F).Base_comp]
      infer_instance
  · exact IsPullBack.comp (pfRootPre P F) ha
      (isPullBack_of_isIso (pfRootPre P F) (inv eY))

/-! ## ★52. (iv)(a) の一意性 —— 3 分解は同型を除いて一意

★★`𝒞^birat` の `birat_arbFactorUniq` と**同じ手**である ——
`frobDegUniq` ＋ 全射性 ＋ `IsPullBack.lift` の 3 点だけで出るので、
どちらの Frobenioid でも同じように書ける。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iv)(a)** の一意性 —— `𝒞^pf` 版。 -/
theorem pfRoot_arbFactorUniq (hfi : IsOfFrobeniusIsotropicType P)
    {A B : PfRootObj P F} (X Y X' Y' : PfRootObj P F)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : γ ≫ β ≫ α = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (pfRootPre P F) γ) (hβ : IsPreStep (pfRootPre P F) β)
    (hα : IsPullBack (pfRootPre P F) α)
    (hγ' : IsFrobeniusType (pfRootPre P F) γ') (hβ' : IsPreStep (pfRootPre P F) β')
    (hα' : IsPullBack (pfRootPre P F) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  -- ★段 1: 次数はどちらも `degFr γ`
  have c3 : (pfRootPre P F).degFr (β ≫ α) = 1 :=
    ((pfRootPre P F).degFr_comp β α).trans (by
      rw [show (pfRootPre P F).degFr α = 1 from (pfRoot_pullBackLB hfi α hα).2,
        show (pfRootPre P F).degFr β = 1 from hβ.1]; simp)
  have c3' : (pfRootPre P F).degFr (β' ≫ α') = 1 :=
    ((pfRootPre P F).degFr_comp β' α').trans (by
      rw [show (pfRootPre P F).degFr α' = 1 from (pfRoot_pullBackLB hfi α' hα').2,
        show (pfRootPre P F).degFr β' = 1 from hβ'.1]; simp)
  have e1 : (pfRootPre P F).degFr (γ ≫ β ≫ α) = (pfRootPre P F).degFr γ :=
    ((pfRootPre P F).degFr_comp γ (β ≫ α)).trans (by rw [c3]; simp)
  have e1' : (pfRootPre P F).degFr (γ' ≫ β' ≫ α') = (pfRootPre P F).degFr γ' :=
    ((pfRootPre P F).degFr_comp γ' (β' ≫ α')).trans (by rw [c3']; simp)
  have hdγ : (pfRootPre P F).degFr γ = (pfRootPre P F).degFr γ' :=
    e1.symm.trans ((congrArg (pfRootPre P F).degFr heq).trans e1')
  -- ★段 2: `ε`
  obtain ⟨e, heiso, hee⟩ := pfRoot_frobDegUniq hfi A X X' γ γ' hγ hγ' hdγ
  haveI := heiso
  -- ★段 3: `γ` は全射なので消せる
  haveI hepiγ : Epi γ := pfRoot_totEpi P F _ _ γ
  have hfac : β ≫ α = (e ≫ β') ≫ α' := by
    have hcancel : β ≫ α = e ≫ (β' ≫ α') := by
      refine hepiγ.left_cancellation _ _ ?_
      have y1 : γ ≫ (e ≫ (β' ≫ α')) = (γ ≫ e) ≫ (β' ≫ α') :=
        (Category.assoc _ _ _).symm
      have y2 : (γ ≫ e) ≫ (β' ≫ α') = γ' ≫ (β' ≫ α') :=
        congrArg (fun t => t ≫ (β' ≫ α')) hee
      exact heq.trans (y1.trans y2).symm
    exact hcancel.trans (Category.assoc _ _ _).symm
  -- ★段 4: 底の同型
  haveI hbβ : IsIso ((pfRootPre P F).Base β) := hβ.2
  haveI hbβ' : IsIso ((pfRootPre P F).Base β') := hβ'.2
  haveI hbe : IsIso ((pfRootPre P F).Base e) := isBaseIsomorphism_of_isIso (pfRootPre P F) e
  have hbc : (pfRootPre P F).Base (e ≫ β')
      = (pfRootPre P F).Base e ≫ (pfRootPre P F).Base β' := (pfRootPre P F).Base_comp e β'
  haveI hbβ'' : IsIso ((pfRootPre P F).Base (e ≫ β')) := by rw [hbc]; infer_instance
  have hbfac : (pfRootPre P F).Base β ≫ (pfRootPre P F).Base α
      = (pfRootPre P F).Base (e ≫ β') ≫ (pfRootPre P F).Base α' :=
    ((pfRootPre P F).Base_comp β α).symm.trans
      ((congrArg (pfRootPre P F).Base hfac).trans ((pfRootPre P F).Base_comp _ _))
  -- ★段 5: `δ` とその逆
  have hb1 : (pfRootPre P F).Base α
      = (@inv _ _ _ _ ((pfRootPre P F).Base β) hbβ ≫ (pfRootPre P F).Base (e ≫ β'))
        ≫ (pfRootPre P F).Base α' := by
    rw [Category.assoc, ← hbfac, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  have hb2 : (pfRootPre P F).Base α'
      = (@inv _ _ _ _ ((pfRootPre P F).Base (e ≫ β')) hbβ'' ≫ (pfRootPre P F).Base β)
        ≫ (pfRootPre P F).Base α := by
    rw [Category.assoc, hbfac, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨d1, ⟨hd1a, hd1b⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα' Y α
    (@inv _ _ _ _ ((pfRootPre P F).Base β) hbβ ≫ (pfRootPre P F).Base (e ≫ β')) hb1
  obtain ⟨d2, ⟨hd2a, hd2b⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα Y' α'
    (@inv _ _ _ _ ((pfRootPre P F).Base (e ≫ β')) hbβ'' ≫ (pfRootPre P F).Base β) hb2
  -- ★段 6: `d1`・`d2` は互いに逆
  have hd12 : d1 ≫ d2 = 𝟙 Y := by
    refine IsPullBack.hom_ext (pfRootPre P F) hα _ _ ?_ ?_
    · have y1 : (d1 ≫ d2) ≫ α = d1 ≫ (d2 ≫ α) := Category.assoc _ _ _
      have y2 : d1 ≫ (d2 ≫ α) = d1 ≫ α' := congrArg (fun t => d1 ≫ t) hd2a
      exact ((y1.trans y2).trans hd1a).trans (Category.id_comp _).symm
    · exact ((pfRootPre P F).Base_comp d1 d2).trans
        (by rw [hd1b, hd2b, (pfRootPre P F).Base_id]; simp)
  have hd21 : d2 ≫ d1 = 𝟙 Y' := by
    refine IsPullBack.hom_ext (pfRootPre P F) hα' _ _ ?_ ?_
    · have y1 : (d2 ≫ d1) ≫ α' = d2 ≫ (d1 ≫ α') := Category.assoc _ _ _
      have y2 : d2 ≫ (d1 ≫ α') = d2 ≫ α := congrArg (fun t => d2 ≫ t) hd1a
      exact ((y1.trans y2).trans hd2a).trans (Category.id_comp _).symm
    · exact ((pfRootPre P F).Base_comp d2 d1).trans
        (by rw [hd2b, hd1b, (pfRootPre P F).Base_id]; simp)
  -- ★段 7: `β ≫ d1 = e ≫ β'`
  have hbd : β ≫ d1 = e ≫ β' := by
    refine IsPullBack.hom_ext (pfRootPre P F) hα' _ _ ?_ ?_
    · have y1 : (β ≫ d1) ≫ α' = β ≫ (d1 ≫ α') := Category.assoc _ _ _
      have y2 : β ≫ (d1 ≫ α') = β ≫ α := congrArg (fun t => β ≫ t) hd1a
      exact ((y1.trans y2).trans hfac)
    · exact ((pfRootPre P F).Base_comp β d1).trans
        (by rw [hd1b, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp])
  -- ★段 8: 組み立て
  refine ⟨⟨d1, d2, hd12, hd21⟩, asIso e, hd2a.symm, ?_, hee.symm⟩
  show β' = inv e ≫ (β ≫ d1)
  have y1 : inv e ≫ (β ≫ d1) = inv e ≫ (e ≫ β') :=
    congrArg (fun t => inv e ≫ t) hbd
  have y2 : inv e ≫ (e ≫ β') = (inv e ≫ e) ≫ β' := (Category.assoc _ _ _).symm
  have y3 : (inv e ≫ e) ≫ β' = 𝟙 _ ≫ β' :=
    congrArg (fun t => t ≫ β') (IsIso.inv_hom_id e)
  exact ((y1.trans (y2.trans y3)).trans (Category.id_comp _)).symm

/-! ## ★53. ★★★★(i)(c) —— `(𝒞^pf)^pl-bk_X → 𝒟_{X_𝒟}` は圏同値

★★3 つとも `𝒞^birat` と同じ手である:

| 条 | 手 |
|---|---|
| 忠実 | `W.hom` の `hom_ext` |
| 充満 | `W.hom` の `lift` ＋ pull-back の右キャンセル |
| 本質的全射 | `𝒞` の `plBk_baseChange` を `Λ_{X.root}` で押し出し、`pfRoot_isPullBack_mk` |
-/

variable {P F} in
/-- ★★`Λ_k` の底は元の底。 -/
@[simp] theorem lamHom_rootBase (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).Base (lamHom (F := F) k φ) = P.Base φ := by
  haveI : IsIso (P.Base (rtExt P F B k)) := (rtExt_frobType P F B k).2
  have h := rootBase_mk_spec (X := (⟨A, k⟩ : PfRootObj P F))
    (Y := (⟨B, k⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)
  have h' : (pfRootPre P F).Base (lamHom (F := F) k φ)
        ≫ (P.Base (rtExt P F B k) ≫ P.Base (𝟙 (rtObj P F B k)))
      = (P.Base (rtExt P F A k) ≫ P.Base (𝟙 (rtObj P F A k)))
        ≫ P.Base (rtMap (F := F) k φ) := h
  simp only [P.Base_id, Category.comp_id] at h'
  refine (cancel_mono (P.Base (rtExt P F B k))).mp ?_
  refine Eq.trans h' ?_
  rw [← P.Base_comp, ← rtMap_spec]
  exact P.Base_comp φ (rtExt P F B k)

variable {P F} in
/-- ★`(𝒞^pf)^pl-bk_X → 𝒟_{X_𝒟}` は忠実。 -/
theorem pfRoot_plBkOver_faithful (X : PfRootObj P F) :
    (plBkOverFunctor (pfRootPre P F) X).Faithful where
  map_injective {Z W} {f g} h := by
    refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
    refine IsPullBack.hom_ext (pfRootPre P F) W.hom.property _ _ ?_ ?_
    · exact (congrArg InducedWideCategory.Hom.hom (Over.w f)).trans
        (congrArg InducedWideCategory.Hom.hom (Over.w g)).symm
    · exact congrArg CommaMorphism.left h

variable {P F} in
/-- ★★`(𝒞^pf)^pl-bk_X → 𝒟_{X_𝒟}` は充満。 -/
theorem pfRoot_plBkOver_full (X : PfRootObj P F) :
    (plBkOverFunctor (pfRootPre P F) X).Full := by
  constructor
  intro Z W h
  obtain ⟨g, ⟨hg1, hg2⟩, -⟩ := IsPullBack.lift (pfRootPre P F) W.hom.property
    Z.left.obj Z.hom.hom h.left (Over.w h).symm
  have hgpb : IsPullBack (pfRootPre P F) g :=
    isPullBack_of_comp_right (pfRootPre P F) g W.hom.hom
      (by rw [hg1]; exact Z.hom.property) W.hom.property
  refine ⟨Over.homMk (show Z.left ⟶ W.left from ⟨g, hgpb⟩)
    (WideSubcategory.hom_ext _ hg1), ?_⟩
  exact Over.OverMorphism.ext hg2

variable {P F} in
/-- ★★`(𝒞^pf)^pl-bk_X → 𝒟_{X_𝒟}` は本質的全射。 -/
theorem pfRoot_plBkOver_essSurj (hfi : IsOfFrobeniusIsotropicType P) (X : PfRootObj P F) :
    (plBkOverFunctor (pfRootPre P F) X).EssSurj := by
  refine ⟨fun T => ?_⟩
  obtain ⟨Yt, αt, k, hαt, hbb⟩ := plBk_baseChange P F X.obj T.hom
  refine ⟨Over.mk (show (⟨(⟨Yt, X.root⟩ : PfRootObj P F)⟩ : PlBk (pfRootPre P F)) ⟶ ⟨X⟩ from
    ⟨lamHom (F := F) X.root αt, ?_⟩), ⟨?_⟩⟩
  · exact pfRoot_isPullBack_mk hfi _
      (idxOne P F (rtObj P F Yt X.root) (rtObj P F X.obj X.root))
      (rtMap (F := F) X.root αt) rfl (rtMap_isPullBack (F := F) X.root αt hαt)
  · refine Over.isoMk k ?_
    show k.hom ≫ T.hom = (pfRootPre P F).Base (lamHom (F := F) X.root αt)
    rw [lamHom_rootBase]
    exact hbb.symm

variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (i)(c)** の `𝒞^pf` 版。 -/
theorem pfRoot_plBkEquiv (hfi : IsOfFrobeniusIsotropicType P) (X : PfRootObj P F) :
    (plBkOverFunctor (pfRootPre P F) X).IsEquivalence :=
  ⟨pfRoot_plBkOver_faithful X, pfRoot_plBkOver_full X, pfRoot_plBkOver_essSurj hfi X⟩

/-! ## ★54. ★★★★★**組み立て** —— `𝒞^pf` は `Definition 1.3` の 21 条を満たす

原文 (FrdI p.59):
> (i), is a Frobenioid of perfect and isotropic type. Moreover, there is a natural
-/

variable {P F} in
/-- ★★★★★**[FrdI] Proposition 3.2, (iii)** の中核 ——
`𝒞` が Frobenius-isotropic 型の Frobenioid なら、
**`𝒞^pf` は `Definition 1.3` の 21 条をすべて満たす**。 -/
theorem pfRootCore (hfi : IsOfFrobeniusIsotropicType P) :
    FrobenioidCore (pfRootPre P F) where
  baseSurj := pfRoot_baseSurj hfi
  preStepSpan := pfRoot_preStepSpan
  plBkEquiv := pfRoot_plBkEquiv hfi
  frobDegSurj := pfRoot_frobDegSurj hfi
  frobDegUniq := pfRoot_frobDegUniq hfi
  coAngularComp := pfRoot_coAngularComp hfi
  coAngularOfPreStep := pfRoot_coAngularOfPreStep hfi
  otriFwd := pfRoot_otriFwd hfi
  otriBwd := pfRoot_otriBwd hfi
  otriBase := pfRoot_otriBase hfi
  arbFactor := pfRoot_arbFactor hfi
  arbFactorUniq := pfRoot_arbFactorUniq hfi
  pullBackLB := pfRoot_pullBackLB hfi
  preStepMono := pfRoot_preStepMono
  preStepFactor := pfRoot_preStepFactor hfi
  preStepFactorUniq := pfRoot_preStepFactorUniq hfi
  preStepFactor' := pfRoot_preStepFactor' hfi
  preStepFactorUniq' := pfRoot_preStepFactorUniq' hfi
  faithfulUpToUnits := pfRoot_faithfulUpToUnits hfi
  isotropicHullExists := pfRoot_isotropicHullExists hfi
  isotropicClosed := pfRoot_isotropicClosed hfi

/-! ## ★55. (iii)(d) への下ごしらえ —— `Λ_k` の零因子と、2 つの関手の忠実性 -/

/-- ★`Pf` の約分。 -/
theorem Pf.mk_smul_cancel {M : Type w} [AddCommMonoid M] (m : M) (k a : ℕ+) :
    Pf.mk (((k : ℕ+) : ℕ) • m) (k * a) = Pf.mk m a :=
  Pf.sound 1 (by
    push_cast
    rw [smul_smul, show (1 : ℕ) * ((a : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)
        = 1 * (((k : ℕ+) : ℕ) * ((a : ℕ+) : ℕ)) from by ring])

variable {P F} in
/-- ★★**`Λ_k` の零因子は `Div φ / k`**。 -/
theorem rootDiv_lamHom (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).Div (lamHom (F := F) k φ) = Pf.mk (P.Div φ) k := by
  show Pf.divBy (k * k) (Pf.map (Φ.map (P.Base (rtExt P F A k)))
    (pfDiv (toHomPf (F := F) (rtMap (F := F) k φ)))) = _
  rw [pfDiv_toHomPf, Pf.map_mk, Pf.divBy_mk, rtMap_Div, mul_one]
  exact Pf.mk_smul_cancel (P.Div φ) k k

/-- ★★`_A(𝒞^coa-pre) → Order(Φ(A))` は**つねに忠実**(全射性から)。 -/
theorem coaPreUnder_faithful {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3)
    [MorphismProperty.IsMultiplicative (coaPreProp Q)] (A : C3) :
    (coaPreUnderFunctor Q A).Faithful where
  map_injective {Z W} {f g} _ := by
    refine Under.UnderMorphism.ext (InducedWideCategory.Hom.ext ?_)
    haveI : Epi Z.hom.hom := Q.totEpiC _ _ _
    exact (cancel_epi Z.hom.hom).mp
      ((congrArg InducedWideCategory.Hom.hom (Under.w f)).trans
        (congrArg InducedWideCategory.Hom.hom (Under.w g)).symm)

/-- ★★`(𝒞^coa-pre)_A → Order(Φ(A))^opp` は**pre-step が mono なら忠実**。 -/
theorem coaPreOver_faithful {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3)
    [MorphismProperty.IsMultiplicative (coaPreProp Q)]
    (hmono : ∀ {X Y : C3} (φ : X ⟶ Y), IsPreStep Q φ → Mono φ) (A : C3) :
    (coaPreOverFunctor Q A).Faithful where
  map_injective {Z W} {f g} _ := by
    refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
    haveI : Mono W.hom.hom := hmono W.hom.hom W.hom.property.2
    exact (cancel_mono W.hom.hom).mp
      ((congrArg InducedWideCategory.Hom.hom (Over.w f)).trans
        (congrArg InducedWideCategory.Hom.hom (Over.w g)).symm)

/-! ## ★56. (iii)(d) の本質的全射性(前置の側)

★★**手**: `X = (A,r)` と `x = a/n` に対し、まず `X ≅ (A^{(n)}, r·n)` で根を上げ、
`𝒞` の `coaPreUnderEquiv` の本質的全射性を `A^{(n)}` で `r·a` に当て、
その co-angular pre-step を `Λ_{r·n}` で押し出す。
★`rootDiv (Λ_k ψ) = Div ψ / k` なので、出てくる零因子は `(r·a)/(r·n) = a/n` である。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**(iii)(d) 前置の本質的全射性** —— `𝒞^pf` 版。 -/
theorem pfRoot_coaPreUnder_essSurj (hfi : IsOfFrobeniusIsotropicType P) (G : Frobenioid P)
    (X : PfRootObj P F) :
    letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
    (coaPreUnderFunctor (pfRootPre P F) X).EssSurj := by
  letI := coaPreProp_isMultiplicative (pfRootPre P F) (pfRoot_coAngularComp hfi)
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  refine ⟨fun x => ?_⟩
  obtain ⟨a, n, rfl⟩ : ∃ (a : Φ.val (P.toElem.obj X.obj).base) (n : ℕ+),
      x = Pf.mk a n :=
    Pf.inductionOn (p := fun y => ∃ (a : Φ.val (P.toElem.obj X.obj).base) (n : ℕ+),
      y = Pf.mk a n) x (fun m b => ⟨m, b, rfl⟩)
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root n (X.root * n) rfl
  haveI := hXiso
  haveI hbiso : IsIso ((pfRootPre P F).Base eX) :=
    isBaseIsomorphism_of_isIso (pfRootPre P F) eX
  obtain ⟨Z₀, ⟨e₀⟩⟩ := (G.coaPreUnderEquiv (rtObj P F X.obj n)).essSurj.mem_essImage
    (toOrderCat (Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base eX) hbiso)
      (((X.root : ℕ+) : ℕ) • a)))
  have hdiv0 : P.Div Z₀.hom.hom
      = Φ.map (@inv _ _ _ _ ((pfRootPre P F).Base eX) hbiso) (((X.root : ℕ+) : ℕ) • a) :=
    mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2
      (leOfHom e₀.hom) (leOfHom e₀.inv)
  -- ★押し出した射
  obtain ⟨φ, hφdef⟩ : ∃ t : X ⟶ (⟨Z₀.right.obj, X.root * n⟩ : PfRootObj P F),
      t = eX ≫ lamHom (F := F) (X.root * n) Z₀.hom.hom := ⟨_, rfl⟩
  have hφstep : IsPreStep (pfRootPre P F) φ := by
    rw [hφdef]
    exact IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) eX)
      (lamHom_isPreStep (X.root * n) Z₀.hom.hom Z₀.hom.property.2)
  have hφdiv : (pfRootPre P F).Div φ = Pf.mk a n := by
    rw [hφdef, (pfRootPre P F).Div_comp,
      show (pfRootPre P F).Div eX = 0 from isIsometric_of_isIso (pfRootPre P F) eX,
      smul_zero, add_zero, rootDiv_lamHom]
    show Pf.map (Φ.map ((pfRootPre P F).Base eX)) (Pf.mk (P.Div Z₀.hom.hom) (X.root * n))
      = Pf.mk a n
    rw [Pf.map_mk, hdiv0, ← Φ.map_comp, IsIso.hom_inv_id, Φ.map_id]
    exact Pf.mk_smul_cancel a X.root n
  refine ⟨Under.mk (Y := (⟨(⟨Z₀.right.obj, X.root * n⟩ : PfRootObj P F)⟩ :
      WideSubcategory (coaPreProp (pfRootPre P F))))
    (show (⟨X⟩ : WideSubcategory (coaPreProp (pfRootPre P F))) ⟶ _ from
      ⟨φ, pfRoot_isCoAngular hfi φ, hφstep⟩), ⟨eqToIso ?_⟩⟩
  show toOrderCat ((pfRootPre P F).Div φ) = toOrderCat (Pf.mk a n)
  rw [hφdiv]
  rfl

/-! ## ★57. (iii)(d) の充満性への下ごしらえ -/

variable {P F} in
/-- ★**始域を共有する 2 射を同じ 3 脚添字へ**(`exists_rep3_cospan` の双対)。 -/
theorem exists_rep3_span {A B E : C} (f : HomPf P F A B) (g : HomPf P F A E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (ψ : V.right.obj.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx13 P F A B E).obj V) ψ := by
  obtain ⟨Zf, φ₀, rfl⟩ := HomPf.exists_rep f
  obtain ⟨Zg, ψ₀, rfl⟩ := HomPf.exists_rep g
  obtain ⟨hfa, hfb, hfab⟩ := Zf.hom.property
  obtain ⟨hga, hge, hgae⟩ := Zg.hom.property
  obtain ⟨E₁, e₁, he₁, he₁d⟩ := F.frobDegSurj E (P.degFr Zf.hom.hom.1)
  obtain ⟨B₂, b₂, hb₂, hb₂d⟩ := F.frobDegSurj B (P.degFr Zg.hom.hom.1)
  set Vf : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(Zf.right.obj.1, Zf.right.obj.2, E₁)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(Zf.hom.hom.1, Zf.hom.hom.2, e₁), hfa, hfb, he₁, hfab,
          hfab.symm.trans he₁d.symm⟩) with hVf
  set Vg : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(Zg.right.obj.1, B₂, Zg.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(Zg.hom.hom.1, b₂, Zg.hom.hom.2), hga, hb₂, hge, hb₂d.symm,
          hb₂d.trans hgae⟩) with hVg
  exact ⟨IsFiltered.max Vf Vg,
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.leftToMax Vf Vg)) φ₀,
    idxTransport P F ((idx13 P F A B E).map (IsFiltered.rightToMax Vf Vg)) ψ₀,
    (HomPf.mk_map _ φ₀).symm, (HomPf.mk_map _ ψ₀).symm⟩

/-- ★★**`Pf` の `≼` は、分子を `k` 倍すれば `M` の `≼` に落ちる**。 -/
theorem Pf.mle_num_of_mle {M : Type w} [AddCommMonoid M] {c₁ c₂ : M} {N : ℕ+}
    (h : MLe (Pf.mk c₁ N) (Pf.mk c₂ N)) :
    ∃ k : ℕ+, MLe (((k : ℕ+) : ℕ) • c₁) (((k : ℕ+) : ℕ) • c₂) := by
  obtain ⟨x, hx⟩ := h
  obtain ⟨f, k, rfl⟩ : ∃ (f : M) (k : ℕ+), x = Pf.mk f k :=
    Pf.inductionOn (p := fun y => ∃ (f : M) (k : ℕ+), y = Pf.mk f k) x
      (fun m b => ⟨m, b, rfl⟩)
  rw [Pf.mk_add_mk] at hx
  obtain ⟨j, hj⟩ := Quotient.exact hx
  refine ⟨j * N * k, ⟨(((j : ℕ+) : ℕ) * ((N : ℕ+) : ℕ)) • (((N : ℕ+) : ℕ) • f), ?_⟩⟩
  rw [smul_add, smul_smul] at hj
  push_cast at hj ⊢
  rw [← mul_assoc] at hj
  exact hj

/-! ## ★58. 根の取り替えを通した零因子の公式

★★`rootBase_rtRootIso_mk_spec` の零因子版。★★**同じ 3 脚添字の第 1 脚を共有する
2 本の射は、分母も分子の写像も一致する**——これが (iii)(d) の充満性の要である。 -/

variable {P F} in
/-- ★★**根の取り替えを通した零因子**。 -/
theorem rootDiv_rtRootIso_mk {X Y : PfRootObj P F} {e tA tB : ℕ+}
    (hA : tA = e * Y.root) (hB : tB = e * X.root)
    (W : IdxPf P F (rtObj P F X.obj tA) (rtObj P F Y.obj tB))
    (φ : W.right.obj.1 ⟶ W.right.obj.2) :
    (pfRootPre P F).Div (show HomRoot P F X Y from
        (rtRootIso P F X.obj Y.obj hA hB).hom (HomPf.mk W φ))
      = Pf.mk (Φ.map (P.Base (rtExt P F X.obj tA))
          (Φ.map (P.Base W.hom.hom.1) (P.Div φ)))
          (X.root * Y.root * (P.degFr W.hom.hom.1 * e)) := by
  refine Eq.trans (congrArg (fun t : X ⟶ Y => (pfRootPre P F).Div t)
    (rtRootIso_hom_mk (F := F) X.obj Y.obj hA hB W φ)) ?_
  have hnum : Φ.map (P.Base (rtExt P F X.obj Y.root))
        (Φ.map (P.Base (rtLift P F X.obj hA ≫ W.hom.hom.1)) (P.Div φ))
      = Φ.map (P.Base (rtExt P F X.obj tA)) (Φ.map (P.Base W.hom.hom.1) (P.Div φ)) := by
    rw [P.Base_comp, ← Φ.map_comp,
      show P.Base (rtExt P F X.obj Y.root)
          ≫ P.Base (rtLift P F X.obj hA) ≫ P.Base W.hom.hom.1
        = P.Base (rtExt P F X.obj tA) ≫ P.Base W.hom.hom.1 from by
          rw [← Category.assoc, ← P.Base_comp, rtLift_ext]]
    exact Φ.map_comp _ _ _
  have hden : P.degFr (rtLift P F X.obj hA ≫ W.hom.hom.1)
      = P.degFr W.hom.hom.1 * e := by
    rw [P.degFr_comp, rtLift_degFr]
    rfl
  show Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root)))
        (pfDiv (HomPf.mk ((pushIdx (F := F) (rtLift P F X.obj hA)
          (rtLift_frobType P F X.obj hA) (rtLift P F Y.obj hB)
          (rtLift_frobType P F Y.obj hB)
          (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ))) = _
  rw [pfDiv_mk]
  unfold repDiv
  rw [Pf.map_mk, Pf.divBy_mk]
  show Pf.mk (Φ.map (P.Base (rtExt P F X.obj Y.root))
      (Φ.map (P.Base (rtLift P F X.obj hA ≫ W.hom.hom.1)) (P.Div φ)))
      (X.root * Y.root * P.degFr (rtLift P F X.obj hA ≫ W.hom.hom.1)) = _
  rw [hnum, hden]

/-! ## ★32. ★★★★現在地と残り(2026-08-17 の測定)

### ★★★★★埋まった 21 条(`FrobenioidCore (pfRootPre P F)` の**全部**)

| 条 | 宣言 | 手 |
|---|---|---|
| (i)(a) `baseSurj` | `pfRoot_baseSurj` | `(A,1)` として送るだけ |
| (i)(b) `preStepSpan` | `pfRoot_preStepSpan` | ★span を `A^m`・`B^n` の間で取り `W := (V₀, n*m)` |
| (ii) `frobDegSurj` | `pfRoot_frobDegSurj` | `exists_rt_frob` |
| (ii) `frobDegUniq` | `pfRoot_frobDegUniq` | co-span を揃え、isotropic に押し上げる |
| (iii)(a) `coAngularComp` | `pfRoot_coAngularComp` | ★isotropic 型ゆえ全射が co-angular |
| (iii)(b) `coAngularOfPreStep` | `pfRoot_coAngularOfPreStep` | 同上 |
| (iii)(c) `otriFwd` | `pfRoot_otriFwd` | 根を揃える ＋ 対角 ＋ `𝒞` の `otriFwd` |
| (iii)(c) `otriBwd` | `pfRoot_otriBwd` | 同上 |
| (iii)(c) `otriBase` | `pfRoot_otriBase` | 同上 |
| (v)(a) `preStepMono` | `pfRoot_preStepMono` | `homPf_cancel_preStep` |
| (v)(b) `preStepFactor` | `pfRoot_preStepFactor` | `φ = φ ≫ 𝟙` |
| (v)(b) `preStepFactorUniq` | `pfRoot_preStepFactorUniq` | isometric pre-step は同型 |
| (v)(c) `preStepFactor'` | `pfRoot_preStepFactor'` | `φ = 𝟙 ≫ φ` |
| (v)(c) `preStepFactorUniq'` | `pfRoot_preStepFactorUniq'` | 同上 |
| (vi) `faithfulUpToUnits` | `pfRoot_faithfulUpToUnits` | ★根を揃える ＋ `rootDiv` を `Div` に降ろす(divisorial) |
| (iv)(b) `pullBackLB` | `pfRoot_pullBackLB` | ★★Frobenius 分解 ＋ linear 分解の 2 段で持ち上げる |
| (iv)(a) `arbFactor` | `pfRoot_arbFactor` | ★★Frobenius 分解 ＋ linear 分解 ＋ `pfRoot_isPullBack_mk` |
| (iv)(a) `arbFactorUniq` | `pfRoot_arbFactorUniq` | `frobDegUniq` ＋ 全射性 ＋ `IsPullBack.lift` |
| (i)(c) `plBkEquiv` | `pfRoot_plBkEquiv` | ★`plBk_baseChange` を `Λ_{X.root}` で押し出す |
| (vii)(a) `isotropicHullExists` | `pfRoot_isotropicHullExists` | `𝟙` |
| (vii)(b) `isotropicClosed` | `pfRoot_isotropicClosed` | 全対象 isotropic |

★土台は `pfRoot_isOfIsotropicType`(`𝒞` が Frobenius-isotropic 型 ⟹ `𝒞^pf` は isotropic 型)。

### ★★★★★2026-08-17 —— `FrobenioidCore` は**閉じた**(`pfRootCore`)

★律速は **`pfRoot_isPullBack_mk`**(代表元が `𝒞` の pull-back なら `𝒞^pf` でも
pull-back)の 1 本だった。これで `arbFactor` / `arbFactorUniq` / `plBkEquiv` が
順に落ちた。手は普遍性の 2 側:

| 側 | 宣言 | 手 |
|---|---|---|
| 一意性 | `pfRoot_isPullBack_mk_hext` | `compRoot_rep_pairL` で共通代表元 → 第 2 脚を isotropic に押し上げ → `homPf_cancel_pullBack` |
| 持ち上げ | `pfRoot_isPullBack_mk_lift` | `exists_rep3_cospan` で 3 脚へ → `𝒞` の `IsPullBack.lift` → `compRoot_mk3` で戻す |

★★**底の追跡が要点**だった。`compRoot` の 3 つの `rtRootIso` はそれぞれ違う
持ち上げを使うが、`rtLift_ext`(`rtExt d ≫ rtLift = rtExt t`)で潰すと
**3 つとも `Base(rtExt · t) ≫ Base(脚)` という同じ同型**になる
(`rootBase_rtRootIso_mk_spec`)。

### `Proposition 3.2` 全体として残るもの

1. ~~**`Frobenioid` の 2 本の圏同値** `coaPreUnderEquiv` / `coaPreOverEquiv`((iii)(d))~~
   ★★★**2026-08-18 に 6 条すべて閉じた**(`Prop32Equiv.lean` の
   `pfRoot_frobenioid`)。表は `Prop32Equiv.lean` の冒頭にある。
2. **(ii) の辞書の残り** —— Frobenius 型・pull-back・co-angular・
   base-identity 自己射・同型(`isPreStep` など 5 項は `Prop32.lean` で済み)
3. **(iii) の後半** —— 「perfect 型」と `𝒞^pf ≃ (𝒞^pf)^pf`
-/

/-- ★**locator** —— `Proposition 3.2, (iii)` の**中核**(`Definition 1.3` の 21 条)。

★★**条つき**(命題全体の完全実装の主張ではない)。残るのは
`Frobenioid` の 2 本の圏同値((iii)(d))、(ii) の辞書の残り、
「perfect 型」、`𝒞^pf ≃ (𝒞^pf)^pf` である。 -/
def pfRootCore.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 58, item := "Proposition 3.2, (iii) — Definition 1.3 の 21 条",
    sectionId := "frdi-prop-3-2" }

end ABC3.Found.FrdI
