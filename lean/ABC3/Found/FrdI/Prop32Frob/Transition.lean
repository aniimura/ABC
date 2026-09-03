/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop32Frob.Triple

/-!
# Prop32Frob —— 根を上げる道具・遷移は同型を保つ・pre-step

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★16. 根を上げる道具

★★`𝒞` の射を「両側の `k` 乗根の間」へ運ぶ。★添字圏の遷移として書けば
`idxTransport_isPreStep`(と `repDeg_map`)がそのまま使える。 -/

variable {P F} in
/-- ★`k` 乗根へ上げる添字。 -/
noncomputable def idxPow (U V : C) (k : ℕ+) : IdxPf P F U V :=
  idxMk (P := P) (F := F) (rtExt P F U k) (rtExt P F V k)
    (rtExt_frobType P F U k) (rtExt_frobType P F V k)
    (by rw [rtExt_degFr, rtExt_degFr])

variable {P F} in
/-- ★始対象からの遷移射。 -/
noncomputable def idxPowHom (U V : C) (k : ℕ+) : idxOne P F U V ⟶ idxPow (F := F) U V k :=
  Under.homMk (show (⟨(U, V)⟩ : BiFr P F) ⟶ (⟨(rtObj P F U k, rtObj P F V k)⟩ : BiFr P F) from
    ⟨(rtExt P F U k, rtExt P F V k), rtExt_frobType P F U k, rtExt_frobType P F V k,
      by rw [rtExt_degFr, rtExt_degFr]⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))

variable {P F} in
/-- ★★**`k` 乗へ運んだ射**。 -/
noncomputable def liftPow {U V : C} (p : U ⟶ V) (k : ℕ+) :
    rtObj P F U k ⟶ rtObj P F V k :=
  idxTransport P F (idxPowHom (F := F) U V k) p

variable {P F} in
/-- ★四角形。 -/
theorem liftPow_spec {U V : C} (p : U ⟶ V) (k : ℕ+) :
    p ≫ rtExt P F V k = rtExt P F U k ≫ liftPow (F := F) p k :=
  idxTransport_spec (F := F) (idxPowHom (F := F) U V k) p

variable {P F} in
/-- ★pre-step は保たれる。 -/
theorem liftPow_isPreStep {U V : C} (p : U ⟶ V) (k : ℕ+) (hp : IsPreStep P p) :
    IsPreStep P (liftPow (F := F) p k) :=
  idxTransport_isPreStep (F := F) (idxPowHom (F := F) U V k) p hp

variable {P F} in
/-- ★★**`rtObj (rtObj A m) n` と `rtObj A (n * m)` を同一視する**。 -/
theorem exists_rtObj_assoc (A : C) (m n t : ℕ+) (ht : t = n * m) :
    ∃ β : rtObj P F (rtObj P F A m) n ⟶ rtObj P F A t,
      IsIso β ∧ rtExt P F A m ≫ rtExt P F (rtObj P F A m) n ≫ β = rtExt P F A t := by
  have h1 : IsFrobeniusType P (rtExt P F A m ≫ rtExt P F (rtObj P F A m) n) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A m) (rtExt_frobType P F (rtObj P F A m) n)
  have h2 : IsFrobeniusType P (rtExt P F A t) := rtExt_frobType P F A t
  have hdeg : P.degFr (rtExt P F A m ≫ rtExt P F (rtObj P F A m) n)
      = P.degFr (rtExt P F A t) := by
    rw [P.degFr_comp, rtExt_degFr, rtExt_degFr, rtExt_degFr, ht]
  obtain ⟨β, hβ, hβe⟩ := F.frobDegUniq A _ _ _ _ h1 h2 hdeg
  exact ⟨β, hβ, (Category.assoc _ _ _).symm.trans hβe⟩

/-! ## ★17. (i)(b) —— 底の同型は pre-step の span で持ち上がる

★★根の取り方が要点である。`X = (A, n)`、`Y = (B, m)` に対し
**`𝒞` の span を `A^m` と `B^n` の間で取り**、その頂点 `V₀` を
`W := (V₀, n * m)` として使う。
★モデルで言えば `W = V₀/(nm) ≤ A/n ⟺ V₀ ≤ mA` であり、
これがちょうど `V₀ ⟶ A^m` という pre-step である。 -/

variable {P F} in
/-- ★始対象の添字での `repBase` は `Base` そのもの。 -/
theorem repBase_idxOne {A B : C} (φ : A ⟶ B) :
    repBase (idxOne P F A B) φ = P.Base φ := by
  have h := repBase_spec (F := F) (idxOne P F A B) φ
  show repBase (idxOne P F A B) φ = P.Base φ
  have h2 : repBase (idxOne P F A B) φ ≫ P.Base (𝟙 B) = P.Base (𝟙 A) ≫ P.Base φ := h
  rw [P.Base_id, P.Base_id, Category.comp_id, Category.id_comp] at h2
  exact h2

/-- ★`f ≫ α = g` から `α = f⁻¹ ≫ g`。 -/
theorem eq_inv_comp_of {W X Y : D} (f : W ⟶ X) (hf : IsIso f) (g : W ⟶ Y) (α : X ⟶ Y)
    (h : f ≫ α = g) : α = @inv _ _ _ _ f hf ≫ g := by
  haveI := hf
  rw [← h, IsIso.inv_hom_id_assoc]

variable {P F} in
/-- ★★**(i)(b)** —— `𝒞^pf` 版。 -/
theorem pfRoot_preStepSpan (X Y : PfRootObj P F)
    (α : ((pfRootPre P F).toElem.obj X).base ⟶ ((pfRootPre P F).toElem.obj Y).base)
    (hα : IsIso α) :
    ∃ (W : PfRootObj P F) (φ : W ⟶ X) (ψ : W ⟶ Y) (hφ : IsPreStep (pfRootPre P F) φ),
      IsPreStep (pfRootPre P F) ψ ∧
        α = @inv _ _ _ _ ((pfRootPre P F).Base φ) hφ.2 ≫ (pfRootPre P F).Base ψ := by
  haveI := hα
  haveI hieA : IsIso (P.Base (rtExt P F X.obj Y.root)) := (rtExt_frobType P F X.obj Y.root).2
  haveI hieB : IsIso (P.Base (rtExt P F Y.obj X.root)) := (rtExt_frobType P F Y.obj X.root).2
  haveI hinv : IsIso (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
      ≫ P.Base (rtExt P F Y.obj X.root)) :=
    IsIso.comp_isIso' (IsIso.inv_isIso) (IsIso.comp_isIso' hα hieB)
  obtain ⟨V₀, p, q, hp, hq, hspan⟩ := F.preStepSpan (rtObj P F X.obj Y.root)
    (rtObj P F Y.obj X.root)
    (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
      ≫ P.Base (rtExt P F Y.obj X.root)) hinv
  haveI hip : IsIso (P.Base p) := hp.2
  haveI hiq : IsIso (P.Base q) := hq.2
  obtain ⟨βA, hβA, hβAe⟩ := exists_rtObj_assoc (F := F) X.obj Y.root X.root (X.root * Y.root) rfl
  obtain ⟨βB, hβB, hβBe⟩ := exists_rtObj_assoc (F := F) Y.obj X.root Y.root (X.root * Y.root)
    (mul_comm _ _)
  haveI := hβA
  haveI := hβB
  have hφs : IsPreStep P (liftPow (F := F) p X.root ≫ βA) :=
    IsPreStep.comp P (liftPow_isPreStep (F := F) p X.root hp) (isPreStep_of_isIso P βA)
  have hψs : IsPreStep P (liftPow (F := F) q Y.root ≫ βB) :=
    IsPreStep.comp P (liftPow_isPreStep (F := F) q Y.root hq) (isPreStep_of_isIso P βB)
  have hbase : ∀ {U V : C} (r : V₀ ⟶ U) (k t : ℕ+)
      (L : rtObj P F V₀ k ⟶ rtObj P F U k)
      (_hL : r ≫ rtExt P F U k = rtExt P F V₀ k ≫ L)
      (γ : rtObj P F U k ⟶ rtObj P F V t) (eU : V ⟶ U) (heU : IsIso (P.Base eU))
      (_hcomp : eU ≫ rtExt P F U k ≫ γ = rtExt P F V t),
      rootBase (show HomRoot P F ⟨V₀, t⟩ ⟨V, k⟩ from
          HomPf.mk (idxOne P F (rtObj P F V₀ k) (rtObj P F V t)) (L ≫ γ))
        = P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU := by
    intro U V r k t L hL γ eU heU hcomp
    haveI := heU
    have hmor : rtExt P F V₀ k ≫ (L ≫ γ) = r ≫ (rtExt P F U k ≫ γ) := by
      rw [← Category.assoc, ← hL, Category.assoc]
    refine (rootBase_uniq _ _ ?_).symm
    show (P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU) ≫ P.Base (rtExt P F V t)
      = P.Base (rtExt P F V₀ k) ≫ pfBase (HomPf.mk (idxOne P F _ _) _)
    rw [pfBase_mk, repBase_idxOne, ← hcomp]
    have hRHS : P.Base (rtExt P F V₀ k) ≫ P.Base (L ≫ γ)
        = P.Base r ≫ P.Base (rtExt P F U k ≫ γ) :=
      (P.Base_comp _ _).symm.trans ((congrArg P.Base hmor).trans (P.Base_comp _ _))
    have hLHS : (P.Base r ≫ @inv _ _ _ _ (P.Base eU) heU)
          ≫ P.Base (eU ≫ (rtExt P F U k ≫ γ))
        = P.Base r ≫ P.Base (rtExt P F U k ≫ γ) := by
      rw [P.Base_comp eU (rtExt P F U k ≫ γ)]
      simp
    exact hLHS.trans hRHS.symm
  have hbφ := hbase p X.root (X.root * Y.root) (liftPow (F := F) p X.root)
    (liftPow_spec (F := F) p X.root) βA (rtExt P F X.obj Y.root) hieA hβAe
  have hbψ := hbase q Y.root (X.root * Y.root) (liftPow (F := F) q Y.root)
    (liftPow_spec (F := F) q Y.root) βB (rtExt P F Y.obj X.root) hieB hβBe
  refine ⟨⟨V₀, X.root * Y.root⟩,
    HomPf.mk (idxOne P F _ _) (liftPow (F := F) p X.root ≫ βA),
    HomPf.mk (idxOne P F _ _) (liftPow (F := F) q Y.root ≫ βB),
    (isPreStep_mk_iff (X := (⟨V₀, X.root * Y.root⟩ : PfRootObj P F)) (Y := X) _ _).mpr hφs,
    (isPreStep_mk_iff (X := (⟨V₀, X.root * Y.root⟩ : PfRootObj P F)) (Y := Y) _ _).mpr hψs, ?_⟩
  have hkey : rootBase (show HomRoot P F ⟨V₀, X.root * Y.root⟩ X from
        HomPf.mk (idxOne P F _ _) (liftPow (F := F) p X.root ≫ βA)) ≫ α
      = rootBase (show HomRoot P F ⟨V₀, X.root * Y.root⟩ Y from
        HomPf.mk (idxOne P F _ _) (liftPow (F := F) q Y.root ≫ βB)) := by
    have h1 : P.Base p ≫ (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ α
        ≫ P.Base (rtExt P F Y.obj X.root)) = P.Base q := by
      refine (congrArg (fun t => P.Base p ≫ t) hspan).trans ?_
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    have hcancel : P.Base (rtExt P F Y.obj X.root)
        ≫ @inv _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) hieB = 𝟙 _ :=
      @IsIso.hom_inv_id _ _ _ _ (P.Base (rtExt P F Y.obj X.root)) hieB
    rw [hbφ, hbψ, ← h1]
    simp only [Category.assoc]
    refine congrArg (fun t => P.Base p
      ≫ (@inv _ _ _ _ (P.Base (rtExt P F X.obj Y.root)) hieA ≫ t)) ?_
    exact (Category.comp_id α).symm.trans (congrArg (fun t => α ≫ t) hcancel.symm)
  exact eq_inv_comp_of _ _ _ _ hkey

/-! ## ★18. 遷移は同型を保つ

★★逆射の遷移が逆射になる。★示すのは `a ≫ (T₁ ≫ T₂) = a ≫ 𝟙` で、
`a` が epi(`𝒞` は totally epimorphic)なので消せる。 -/

variable {P F} in
/-- ★★**遷移は同型を保つ**。 -/
theorem frobTransport_isIso {A' B' A'' B'' : C}
    (a : A' ⟶ A'') (ha : IsFrobeniusType P a) (b : B' ⟶ B'') (hb : IsFrobeniusType P b)
    (hd : P.degFr a = P.degFr b) (θ : A' ⟶ B') (hθ : IsIso θ) :
    IsIso (frobTransport (F := F) a ha b hb hd θ) := by
  haveI := hθ
  haveI hea : Epi a := P.totEpiC _ _ _
  haveI heb : Epi b := P.totEpiC _ _ _
  have s1 : θ ≫ b = a ≫ frobTransport (F := F) a ha b hb hd θ :=
    frobTransport_spec _ _ _ _ _ θ
  have s2 : inv θ ≫ a = b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) :=
    frobTransport_spec _ _ _ _ _ (inv θ)
  refine ⟨frobTransport (F := F) b hb a ha hd.symm (inv θ), ?_, ?_⟩
  · refine (cancel_epi a).mp ?_
    have e1 : a ≫ (frobTransport (F := F) a ha b hb hd θ
        ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
        = (a ≫ frobTransport (F := F) a ha b hb hd θ)
          ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) := (Category.assoc _ _ _).symm
    have e2 : (a ≫ frobTransport (F := F) a ha b hb hd θ)
          ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)
        = (θ ≫ b) ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ) :=
      congrArg (fun t => t ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)) s1.symm
    have e3 : (θ ≫ b) ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)
        = θ ≫ (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ)) := Category.assoc _ _ _
    have e4 : θ ≫ (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
        = θ ≫ (inv θ ≫ a) := congrArg (fun t => θ ≫ t) s2.symm
    have e5 : θ ≫ (inv θ ≫ a) = a := by
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    exact ((((e1.trans e2).trans e3).trans e4).trans e5).trans (Category.comp_id a).symm
  · refine (cancel_epi b).mp ?_
    have e1 : b ≫ (frobTransport (F := F) b hb a ha hd.symm (inv θ)
        ≫ frobTransport (F := F) a ha b hb hd θ)
        = (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
          ≫ frobTransport (F := F) a ha b hb hd θ := (Category.assoc _ _ _).symm
    have e2 : (b ≫ frobTransport (F := F) b hb a ha hd.symm (inv θ))
          ≫ frobTransport (F := F) a ha b hb hd θ
        = (inv θ ≫ a) ≫ frobTransport (F := F) a ha b hb hd θ :=
      congrArg (fun t => t ≫ frobTransport (F := F) a ha b hb hd θ) s2.symm
    have e3 : (inv θ ≫ a) ≫ frobTransport (F := F) a ha b hb hd θ
        = inv θ ≫ (a ≫ frobTransport (F := F) a ha b hb hd θ) := Category.assoc _ _ _
    have e4 : inv θ ≫ (a ≫ frobTransport (F := F) a ha b hb hd θ)
        = inv θ ≫ (θ ≫ b) := congrArg (fun t => inv θ ≫ t) s1.symm
    have e5 : inv θ ≫ (θ ≫ b) = b := by
      rw [← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    exact ((((e1.trans e2).trans e3).trans e4).trans e5).trans (Category.comp_id b).symm

variable {P F} in
/-- ★添字の遷移は同型を保つ。 -/
theorem idxTransport_isIso {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (θ : Z.right.obj.1 ⟶ Z.right.obj.2) (hθ : IsIso θ) :
    IsIso (idxTransport P F u θ) :=
  frobTransport_isIso _ _ _ _ _ θ hθ

/-! ## ★19. ★★代表元が同型なら同型(添字を持ち上げない形)

★★`pfRoot_isIso_of_rep` は `compRoot` の定義に合わせた「持ち上げた添字」の形だが、
★実際に使うのは**素の添字**の形である。`pushIdx` の cofinal 性で移す。 -/

variable {P F} in
/-- ★★★**素の添字で代表元が同型なら `𝒞^pf` でも同型**。 -/
theorem pfRoot_isIso_mk {Y Z : PfRootObj P F}
    (W : IdxPf P F (rtObj P F Y.obj Z.root) (rtObj P F Z.obj Y.root))
    (β₀ : W.right.obj.1 ⟶ W.right.obj.2) (hβ₀ : IsIso β₀) :
    IsIso (show Y ⟶ Z from HomPf.mk W β₀) := by
  obtain ⟨V, ⟨u⟩⟩ := exists_hom_of_final
    (pushIdx (F := F)
      (rtLift P F Y.obj (show Y.root * Z.root = Y.root * Z.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (rtLift P F Z.obj (show Y.root * Y.root = Y.root * Y.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])) W
  set β₁ := idxTransport P F u β₀ with hβ₁
  have hiso : IsIso β₁ := by rw [hβ₁]; exact idxTransport_isIso (F := F) u β₀ hβ₀
  have h1 : (show Y ⟶ Z from HomPf.mk W β₀)
      = (rtRootIso P F Y.obj Z.obj (show Y.root * Z.root = Y.root * Z.root from rfl)
          (show Y.root * Y.root = Y.root * Y.root from rfl)).hom
        (HomPf.mk V β₁) := by
    rw [rtRootIso_hom_mk, hβ₁]
    exact (HomPf.mk_map (P := P) (F := F) u β₀).symm
  refine pfRoot_isIso_of_rep (show Y ⟶ Z from HomPf.mk W β₀) V β₁ (hφ := hiso) ?_
  rw [h1, Iso.hom_inv_id_apply]

/-! ## ★20. 同じ始域から出る 2 射を揃える(co-span)

★★`exists_rep3` は「合成できる 2 射」を揃えるが、
`frobDegUniq` に要るのは**同じ始域から出る 2 射**である。
★3 脚添字の `idx12` と `idx13` を使えば同じ形で作れる。 -/

variable {P F} in
/-- ★**co-span を共通の 3 脚添字へ**。 -/
theorem exists_rep_cospan {A B E : C} (f : HomPf P F A B) (g : HomPf P F A E) :
    ∃ (V : IdxPf3 P F A B E) (φ : V.right.obj.1 ⟶ V.right.obj.2.1)
      (χ : V.right.obj.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx12 P F A B E).obj V) φ ∧
      g = HomPf.mk ((idx13 P F A B E).obj V) χ := by
  obtain ⟨Zf, φ₀, hf⟩ := HomPf.exists_rep (P := P) (F := F) f
  obtain ⟨Zg, χ₀, hg⟩ := HomPf.exists_rep (P := P) (F := F) g
  obtain ⟨hfa, hfb, hfab⟩ := Zf.hom.property
  obtain ⟨hga, hge, hgae⟩ := Zg.hom.property
  obtain ⟨E₁, e₁, he₁, he₁d⟩ := F.frobDegSurj E (P.degFr Zf.hom.hom.1)
  obtain ⟨B₂, b₂, hb₂, hb₂d⟩ := F.frobDegSurj B (P.degFr Zg.hom.hom.1)
  refine ⟨IsFiltered.max
      (Under.mk (Y := (⟨(Zf.right.obj.1, Zf.right.obj.2, E₁)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(Zf.hom.hom.1, Zf.hom.hom.2, e₁), hfa, hfb, he₁, hfab,
            hfab.symm.trans he₁d.symm⟩))
      (Under.mk (Y := (⟨(Zg.right.obj.1, B₂, Zg.right.obj.2)⟩ : TriFr P F))
        (show triFrObj P F A B E ⟶ _ from
          ⟨(Zg.hom.hom.1, b₂, Zg.hom.hom.2), hga, hb₂, hge, hb₂d.symm,
            hb₂d.trans hgae⟩)),
    idxTransport P F ((idx12 P F A B E).map (IsFiltered.leftToMax _ _)) φ₀,
    idxTransport P F ((idx13 P F A B E).map (IsFiltered.rightToMax _ _)) χ₀, ?_, ?_⟩
  · rw [HomPf.mk_map]; exact hf.symm
  · rw [HomPf.mk_map]; exact hg.symm

variable {P F} in
/-- ★3 脚添字も「第 1 脚が isotropic」な所まで押し上げられる。 -/
theorem exists_idx3_isotropic (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E) :
    ∃ (W : IdxPf3 P F A B E) (u : V ⟶ W), IsIsotropic P W.right.obj.1 := by
  obtain ⟨Dd, a, ha, hDd⟩ := hfi V.right.obj.1
  obtain ⟨B₂, b, hb, hbd⟩ := F.frobDegSurj V.right.obj.2.1 (P.degFr a)
  obtain ⟨E₂, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr a)
  obtain ⟨hva, hvb, hve, hvab, hvbe⟩ := V.hom.property
  refine ⟨Under.mk (Y := (⟨(Dd, B₂, E₂)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(V.hom.hom.1 ≫ a, V.hom.hom.2.1 ≫ b, V.hom.hom.2.2 ≫ e),
          IsFrobeniusType.comp P F hva ha, IsFrobeniusType.comp P F hvb hb,
          IsFrobeniusType.comp P F hve he, ?_, ?_⟩),
    Under.homMk (show V.right ⟶ (⟨(Dd, B₂, E₂)⟩ : TriFr P F) from
      ⟨(a, b, e), ha, hb, he, hbd.symm, hbd.trans hed.symm⟩)
      (WideSubcategory.hom_ext _ rfl), hDd⟩
  · rw [P.degFr_comp, P.degFr_comp, hbd, hvab]
  · rw [P.degFr_comp, P.degFr_comp, hbd, hed, hvbe]

/-! ## ★21. (ii) の一意性 —— `frobDegUniq`

★★`φ`・`ψ` を **`compRoot` が使う根の高さ**へ持ち上げて co-span を揃え、
添字を isotropic まで押し上げてから `𝒞` の `frobDegUniq` を当てる。
★得た同型 `β₀` を `idx23` の添字で戻せば、合成則は `compPf_mk` そのもの。 -/

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★**(ii) の一意性** —— `𝒞^pf` 版。 -/
theorem pfRoot_frobDegUniq (hfi : IsOfFrobeniusIsotropicType P)
    (X Y Z : PfRootObj P F) (φ : X ⟶ Y) (ψ : X ⟶ Z)
    (hφ : IsFrobeniusType (pfRootPre P F) φ) (hψ : IsFrobeniusType (pfRootPre P F) ψ)
    (hdeg : (pfRootPre P F).degFr φ = (pfRootPre P F).degFr ψ) :
    ∃ β : Y ⟶ Z, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨V, φ₀, χ₀, hφ0, hχ0⟩ := exists_rep_cospan (P := P) (F := F)
    ((rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
      (show Z.root * X.root = Z.root * X.root from rfl)).inv φ)
    ((rtRootIso P F X.obj Z.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = Y.root * X.root from rfl)).inv ψ)
  obtain ⟨W, u, hW⟩ := exists_idx3_isotropic (F := F) hfi V
  set φ₁ := idxTransport P F ((idx12 P F _ _ _).map u) φ₀ with hφ₁
  set χ₁ := idxTransport P F ((idx13 P F _ _ _).map u) χ₀ with hχ₁
  have hφW : (rtRootIso P F X.obj Y.obj (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).inv φ
      = HomPf.mk ((idx12 P F _ _ _).obj W) φ₁ := by
    rw [hφ₁, HomPf.mk_map]; exact hφ0
  have hχW : (rtRootIso P F X.obj Z.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = Y.root * X.root from rfl)).inv ψ
      = HomPf.mk ((idx13 P F _ _ _).obj W) χ₁ := by
    rw [hχ₁, HomPf.mk_map]; exact hχ0
  -- ★`φ`・`ψ` を押し出した添字の代表元として書き直す
  have hφmk : φ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Z.root * Y.root from rfl))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Y.obj (show Z.root * X.root = Z.root * X.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx12 P F _ _ _).obj W)) φ₁ := by
    rw [← rtRootIso_hom_mk (F := F) X.obj Y.obj _ _ ((idx12 P F _ _ _).obj W) φ₁,
      ← hφW, Iso.inv_hom_id_apply]
  have hψmk : ψ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Z.obj (show Y.root * X.root = Y.root * X.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx13 P F _ _ _).obj W)) χ₁ := by
    rw [← rtRootIso_hom_mk (F := F) X.obj Z.obj _ _ ((idx13 P F _ _ _).obj W) χ₁,
      ← hχW, Iso.inv_hom_id_apply]
  set Wφ : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root) :=
    (pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Z.root * Y.root from rfl))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Y.obj (show Z.root * X.root = Z.root * X.root from rfl))
      (rtLift_frobType P F Y.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx12 P F _ _ _).obj W) with hWφ
  set Wψ : IdxPf P F (rtObj P F X.obj Z.root) (rtObj P F Z.obj X.root) :=
    (pushIdx (F := F)
      (rtLift P F X.obj (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _))
      (rtLift_frobType P F X.obj _)
      (rtLift P F Z.obj (show Y.root * X.root = Y.root * X.root from rfl))
      (rtLift_frobType P F Z.obj _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj ((idx13 P F _ _ _).obj W) with hWψ
  -- ★代表元は `𝒞` の Frobenius 型
  have hco : ∀ {U : C} (t : W.right.obj.1 ⟶ U), IsCoAngular P t :=
    fun t => prop_1_4_i P t (fun _ g => F.isotropicClosed g hW)
  have hφ₁F : IsFrobeniusType P φ₁ := by
    refine ⟨⟨hco φ₁, ?_⟩, ?_⟩
    · refine (isIsometric_mk_iff (X := X) (Y := Y) Wφ φ₁).mp ?_
      rw [← hφmk]; exact hφ.1.2
    · refine (isBaseIsomorphism_mk_iff (X := X) (Y := Y) Wφ φ₁).mp ?_
      rw [← hφmk]; exact hφ.2
  have hχ₁F : IsFrobeniusType P χ₁ := by
    refine ⟨⟨hco χ₁, ?_⟩, ?_⟩
    · refine (isIsometric_mk_iff (X := X) (Y := Z) Wψ χ₁).mp ?_
      rw [← hψmk]; exact hψ.1.2
    · refine (isBaseIsomorphism_mk_iff (X := X) (Y := Z) Wψ χ₁).mp ?_
      rw [← hψmk]; exact hψ.2
  have hd1 : P.degFr φ₁ = P.degFr χ₁ := by
    have e1 : (pfRootPre P F).degFr φ = P.degFr φ₁ := by
      rw [hφmk]
      exact (degFr_mk_iff (X := X) (Y := Y) Wφ φ₁ (P.degFr φ₁)).mpr rfl
    have e2 : (pfRootPre P F).degFr ψ = P.degFr χ₁ := by
      rw [hψmk]
      exact (degFr_mk_iff (X := X) (Y := Z) Wψ χ₁ (P.degFr χ₁)).mpr rfl
    rw [← e1, ← e2, hdeg]
  obtain ⟨β₀, hβ₀, hβ₀e⟩ := F.frobDegUniq _ _ _ φ₁ χ₁ hφ₁F hχ₁F hd1
  -- ★戻す
  refine ⟨(rtRootIso P F Y.obj Z.obj
      (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
      (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).hom
    (HomPf.mk ((idx23 P F _ _ _).obj W) β₀), ?_, ?_⟩
  · rw [rtRootIso_hom_mk]
    exact pfRoot_isIso_mk _ _ hβ₀
  · show compRoot P F φ _ = ψ
    unfold compRoot
    rw [hφW, Iso.hom_inv_id_apply, compPf_mk]
    refine Eq.trans (congrArg (ConcreteCategory.hom (rtRootIso P F X.obj Z.obj
        (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = Y.root * X.root from rfl)).hom)
      (congrArg (HomPf.mk ((idx13 P F _ _ _).obj W)) hβ₀e)) ?_
    rw [← hχW, Iso.inv_hom_id_apply]

/-! ## ★22. ★★★根の指数は上げられる —— `(A, n) ≅ (A^m, n*m)`

★★これで **2 つの対象を同じ根の指数に揃えられる**。
`X = (A, n)`、`Y = (B, m)` なら `X ≅ (A^m, n*m)`、`Y ≅ (B^n, n*m)` である。

★★★根が等しいと `compRoot` の 3 つの `rtRootIso` が**すべて同じ命題の証明**になり
(`r*r = r*r`)、証明無関係で**同一の射**になる。★これが (iii)(c) や (vi) の鍵。 -/

variable {P F} in
/-- ★★**根の指数を上げる同型** —— `(A, n) ≅ (A^m, t)`(`t = n * m`)。 -/
theorem pfRoot_exists_iso_root (A : C) (n m t : ℕ+) (ht : t = n * m) :
    ∃ e : (⟨A, n⟩ : PfRootObj P F) ⟶ ⟨rtObj P F A m, t⟩, IsIso e := by
  obtain ⟨β, hβ, -⟩ := exists_rtObj_assoc (F := F) A m n t ht
  haveI := hβ
  exact ⟨HomPf.mk (idxOne P F _ _) (@inv _ _ _ _ β hβ),
    pfRoot_isIso_mk _ _ (@IsIso.inv_isIso _ _ _ _ β hβ)⟩

/-! ## ★23. `𝒪^▷` と (iii)(c) は同型で移せる(圏の言葉だけ) -/

variable {P F} in
/-- ★同型で挟むと `𝒪^▷` の元は `𝒪^▷` の元に移る。 -/
theorem oTri_conj {X X' : PfRootObj P F} (e : X ⟶ X') (e' : X' ⟶ X)
    (he1 : e ≫ e' = 𝟙 X) (he2 : e' ≫ e = 𝟙 X') (α : End X)
    (hα : α ∈ OTri (pfRootPre P F) X) :
    (show End X' from e' ≫ (α : X ⟶ X) ≫ e) ∈ OTri (pfRootPre P F) X' := by
  haveI : IsIso e := ⟨⟨e', he1, he2⟩⟩
  haveI : IsIso e' := ⟨⟨e, he2, he1⟩⟩
  constructor
  · show (pfRootPre P F).Base (e' ≫ (α : X ⟶ X) ≫ e)
      = (pfRootPre P F).Base (𝟙 X')
    rw [(pfRootPre P F).Base_comp, (pfRootPre P F).Base_comp,
      show (pfRootPre P F).Base (α : X ⟶ X) = (pfRootPre P F).Base (𝟙 X) from hα.1,
      (pfRootPre P F).Base_id, Category.id_comp, (pfRootPre P F).Base_id,
      ← (pfRootPre P F).Base_comp, he2, (pfRootPre P F).Base_id]
  · show (pfRootPre P F).degFr (e' ≫ (α : X ⟶ X) ≫ e) = 1
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (α : X ⟶ X) = 1 from hα.2,
      show (pfRootPre P F).degFr e = 1 from isLinear_of_isIso (pfRootPre P F) e,
      show (pfRootPre P F).degFr e' = 1 from isLinear_of_isIso (pfRootPre P F) e']
    simp

/-! ## ★24. 3 脚添字を「対角かつ isotropic」へ

★★`𝒪^▷` の元を `𝒞` の `𝒪^▷` の元として読むには、
**第 1・第 2 脚が同じ射**でなければならない(そうでないと
`Base` が恒等にならず、共役が残る)。
★`frob_common_upper` で 2 脚を揃え、そのうえで `hfi` で isotropic まで押す。 -/

variable {P F} in
/-- ★★**3 脚添字を「第 1・第 2 脚が一致し、終域が isotropic」な所まで押し上げる**。 -/
theorem exists_idx3_diag (hfi : IsOfFrobeniusIsotropicType P) {A E : C}
    (V : IdxPf3 P F A A E) :
    ∃ (Dd E₃ : C) (l : A ⟶ Dd) (hl : IsFrobeniusType P l) (m : E ⟶ E₃)
      (hm : IsFrobeniusType P m) (hd : P.degFr l = P.degFr m),
      IsIsotropic P Dd ∧ Nonempty (V ⟶ idxMk3 (F := F) l l m hl hl hm rfl hd) := by
  obtain ⟨hv1, hv2, hv3, h12, h23⟩ := V.hom.property
  obtain ⟨X₃, a, c, ha, hc, had, hcd, hac⟩ :=
    frob_common_upper P F V.hom.hom.1 hv1 V.hom.hom.2.1 hv2
  obtain ⟨Dd, δ, hδ, hDd⟩ := hfi X₃
  obtain ⟨E₃, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr (a ≫ δ))
  have hacδ : V.hom.hom.1 ≫ (a ≫ δ) = V.hom.hom.2.1 ≫ (c ≫ δ) := by
    rw [← Category.assoc, ← Category.assoc, hac]
    rfl
  have hdaδ : P.degFr (a ≫ δ) = P.degFr (c ≫ δ) := by
    rw [P.degFr_comp a δ, P.degFr_comp c δ, had, hcd, h12]
    rfl
  have hlm : P.degFr (V.hom.hom.1 ≫ (a ≫ δ)) = P.degFr (V.hom.hom.2.2 ≫ e) := by
    rw [P.degFr_comp V.hom.hom.1 (a ≫ δ), P.degFr_comp V.hom.hom.2.2 e, hed, h12, h23]
  refine ⟨Dd, E₃, V.hom.hom.1 ≫ (a ≫ δ),
    IsFrobeniusType.comp P F hv1 (IsFrobeniusType.comp P F ha hδ),
    V.hom.hom.2.2 ≫ e, IsFrobeniusType.comp P F hv3 he, hlm, hDd,
    ⟨Under.homMk (show V.right ⟶ (⟨(Dd, Dd, E₃)⟩ : TriFr P F) from
      ⟨(a ≫ δ, c ≫ δ, e), IsFrobeniusType.comp P F ha hδ,
        IsFrobeniusType.comp P F hc hδ, he, hdaδ, ?_⟩)
      (WideSubcategory.hom_ext _ (Prod.ext rfl (Prod.ext hacδ.symm rfl)))⟩⟩
  rw [← hdaδ, hed]

/-! ## ★25. 対角の添字での `𝒪^▷` の判定

★★第 1・第 2 脚が**同じ射** `l` なら、`repBase = Base l ≫ Base a ≫ (Base l)⁻¹` が
共役なので、★`Base a = 𝟙` と `rootBase = 𝟙` は同値になる。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★**対角の添字での `𝒪^▷` の判定**。 -/
theorem oTri_mk_diag {X : PfRootObj P F} {D : C} (l : rtObj P F X.obj X.root ⟶ D)
    (hl : IsFrobeniusType P l) (a : D ⟶ D) :
    (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      ∈ OTri (pfRootPre P F) X ↔ (P.Base a = P.Base (𝟙 D) ∧ P.degFr a = 1) := by
  haveI hil : IsIso (P.Base l) := hl.2
  haveI hie : IsIso (P.Base (rtExt P F X.obj X.root)) := (rtExt_frobType P F X.obj X.root).2
  have hrep : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
      = P.Base l ≫ P.Base a :=
    repBase_spec (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hroot : rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root)
        ≫ pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a) :=
    rootBase_spec (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
  have hpf : pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
    pfBase_mk (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hidB : 𝟙 ((pfRootPre P F).toElem.obj X).base ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root) := Category.id_comp _
  have hidL : 𝟙 ((P.toElem.obj (rtObj P F X.obj X.root)).base) ≫ P.Base l
      = P.Base l := Category.id_comp _
  constructor
  · rintro ⟨hb, hd⟩
    refine ⟨?_, (degFr_mk_iff (X := X) (Y := X)
      (idxMk (P := P) (F := F) l l hl hl rfl) a 1).mp hd⟩
    have h0 : rootBase (show End X from
        HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        = (pfRootPre P F).Base (𝟙 X) := hb
    rw [(pfRootPre P F).Base_id] at h0
    rw [h0, hpf] at hroot
    -- hroot : 𝟙 ≫ Base (rtExt) = Base (rtExt) ≫ repBase
    have h1 : P.Base (rtExt P F X.obj X.root) ≫ 𝟙 _
        = P.Base (rtExt P F X.obj X.root)
          ≫ repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
      (Category.comp_id _).trans (hidB.symm.trans hroot)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      ((cancel_epi (P.Base (rtExt P F X.obj X.root))).mp h1).symm
    rw [h2] at hrep
    -- hrep : 𝟙 ≫ Base l = Base l ≫ Base a
    have h3 : P.Base l ≫ 𝟙 ((P.toElem.obj D).base) = P.Base l ≫ P.Base a :=
      (Category.comp_id _).trans (hidL.symm.trans hrep)
    rw [P.Base_id]
    exact ((cancel_epi (P.Base l)).mp h3).symm
  · rintro ⟨hb, hd⟩
    refine ⟨?_, (degFr_mk_iff (X := X) (Y := X)
      (idxMk (P := P) (F := F) l l hl hl rfl) a 1).mpr hd⟩
    rw [P.Base_id] at hb
    rw [hb] at hrep
    -- hrep : repBase ≫ Base l = Base l ≫ 𝟙
    have h1 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
        = 𝟙 _ ≫ P.Base l := hrep.trans ((Category.comp_id _).trans hidL.symm)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      (cancel_mono (P.Base l)).mp h1
    rw [hpf, h2] at hroot
    -- hroot : rootBase ≫ Base (rtExt) = Base (rtExt) ≫ 𝟙
    show rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = (pfRootPre P F).Base (𝟙 X)
    rw [(pfRootPre P F).Base_id]
    refine (cancel_mono (P.Base (rtExt P F X.obj X.root))).mp ?_
    exact hroot.trans ((Category.comp_id _).trans hidB.symm)

end ABC3.Found.FrdI
