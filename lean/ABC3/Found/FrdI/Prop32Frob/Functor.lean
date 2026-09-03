/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop32Frob.Faithful

/-!
# Prop32Frob —— 単系の補題・(vi) 単元を除く忠実性

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★33. 単系の補題 —— divisorial なら `n • a = n • b → a = b`

★★`Φ(A)` は divisorial(integral ＋ saturated ＋ sharp)。
`n • (a − b) = 0` は像に入るので **saturated** で `a − b` 自身が像に入り、
`b ≼ a` が出る。対称に `a ≼ b`。**integral の簡約**と **sharp** で `a = b`。

★これが `MetricallyEquivalent` を `Pf` から代表元へ降ろす鍵である。 -/

theorem eq_of_nsmul_eq_of_divisorial {M : Type w} [AddCommMonoid M]
    (hdiv : IsDivisorial M) {a b : M} {n : ℕ} (hn : 0 < n) (hab : n • a = n • b) : a = b := by
  obtain ⟨⟨hint, hsat, -⟩, hsharp⟩ := hdiv
  haveI : IsCancelAdd M := isCancelAdd_of_isIntegralMonoid (M := M) hint
  have key : ∀ x y : M, n • x = n • y → ∃ c : M, y + c = x := by
    intro x y hxy
    have h1 : n • (toGp M x - toGp M y) = 0 := by
      rw [smul_sub, ← toGp_nsmul, ← toGp_nsmul, hxy, sub_self]
    have h2 : n • (toGp M x - toGp M y) ∈ Set.range (toGp M) := by
      rw [h1]; exact ⟨0, toGp_zero M⟩
    obtain ⟨c, hc⟩ := hsat _ n hn h2
    refine ⟨c, hint ?_⟩
    rw [toGp_add, hc]
    abel
  obtain ⟨c, hc⟩ := key a b hab
  obtain ⟨d, hd⟩ := key b a hab.symm
  have h3 : b + (c + d) = b + 0 := by
    rw [← add_assoc, hc, hd, add_zero]
  have h4 : c + d = 0 := add_left_cancel h3
  have h5 : c = 0 := hsharp c ⟨⟨c, d, h4, by rw [add_comm]; exact h4⟩, rfl⟩
  rw [← hc, h5, add_zero]

/-- ★**`Pf` の同じ分母どうしは分子で決まる**(divisorial なら)。 -/
theorem Pf.mk_inj_of_divisorial {M : Type w} [AddCommMonoid M] (hdiv : IsDivisorial M)
    {a b : M} {n : ℕ+} (h : Pf.mk a n = Pf.mk b n) : a = b := by
  obtain ⟨k, hk⟩ := Quotient.exact h
  exact eq_of_nsmul_eq_of_divisorial hdiv (Nat.mul_pos k.pos n.pos) hk

variable {P F} in
/-- ★★**同じ添字なら `rootDiv` の一致から `Div` の一致が出る**。 -/
theorem div_eq_of_rootDiv_eq {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (f f' : Z.right.obj.1 ⟶ Z.right.obj.2)
    (h : rootDiv (show HomRoot P F X Y from HomPf.mk Z f)
      = rootDiv (show HomRoot P F X Y from HomPf.mk Z f')) :
    P.Div f = P.Div f' := by
  replace h : Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv (HomPf.mk Z f)))
    = Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root))) (pfDiv (HomPf.mk Z f'))) := h
  rw [pfDiv_mk, pfDiv_mk] at h
  replace h : Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root)))
        (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div f)) (repRoot Z)))
    = Pf.divBy (X.root * Y.root)
      (Pf.map (Φ.map (P.Base (rtExt P F X.obj Y.root)))
        (Pf.mk (Φ.map (P.Base Z.hom.hom.1) (P.Div f')) (repRoot Z))) := h
  rw [Pf.map_mk, Pf.map_mk, Pf.divBy_mk, Pf.divBy_mk] at h
  have h1 := Pf.mk_inj_of_divisorial (P.divisorial _) h
  exact Φ.map_injective (P.Base Z.hom.hom.1)
    (Φ.map_injective (P.Base (rtExt P F X.obj Y.root)) h1)

/-! ## ★34. (vi) 単元を除く忠実性(根が等しい場合) -/

variable {P F} in
/-- ★平行 2 射を共通の添字へ。 -/
theorem exists_rep_par {A B : C} (f g : HomPf P F A B) :
    ∃ (Z : IdxPf P F A B) (φ ψ : Z.right.obj.1 ⟶ Z.right.obj.2),
      f = HomPf.mk Z φ ∧ g = HomPf.mk Z ψ := by
  obtain ⟨Zf, φ₀, hf⟩ := HomPf.exists_rep (P := P) (F := F) f
  obtain ⟨Zg, ψ₀, hg⟩ := HomPf.exists_rep (P := P) (F := F) g
  exact ⟨IsFiltered.max Zf Zg,
    idxTransport P F (IsFiltered.leftToMax Zf Zg) φ₀,
    idxTransport P F (IsFiltered.rightToMax Zf Zg) ψ₀,
    (by rw [HomPf.mk_map]; exact hf.symm),
    (by rw [HomPf.mk_map]; exact hg.symm)⟩

set_option maxHeartbeats 2000000 in
variable {P F} in
/-- ★★★**(vi)**(根が等しい場合)。 -/
theorem pfRoot_faithfulUpToUnits_sameRoot (hfi : IsOfFrobeniusIsotropicType P)
    {A B : C} {r : ℕ+} (φ ψ : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hbe : BaseEquivalent (pfRootPre P F) φ ψ)
    (hme : MetricallyEquivalent (pfRootPre P F) φ ψ)
    (hφs : IsPreStep (pfRootPre P F) φ) (hψs : IsPreStep (pfRootPre P F) ψ) :
    ∃ α : End (⟨B, r⟩ : PfRootObj P F), α ∈ OTimes (pfRootPre P F) ⟨B, r⟩ ∧
      φ = (ψ ≫ ((α : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩))) := by
  obtain ⟨Z₀, φ₀, ψ₀, hφ0, hψ0⟩ := exists_rep_par (P := P) (F := F)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv φ)
    ((rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)).inv ψ)
  obtain ⟨W, u, hW⟩ := exists_idx_isotropic (F := F) hfi Z₀
  obtain ⟨φ₁, hφ₁⟩ : ∃ t : W.right.obj.1 ⟶ W.right.obj.2,
      t = idxTransport P F u φ₀ := ⟨_, rfl⟩
  obtain ⟨ψ₁, hψ₁⟩ : ∃ t : W.right.obj.1 ⟶ W.right.obj.2,
      t = idxTransport P F u ψ₀ := ⟨_, rfl⟩
  obtain ⟨hw1, hw2, hwd⟩ := W.hom.property
  have hφW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv φ = HomPf.mk W φ₁ := by
    rw [hφ₁, HomPf.mk_map]; exact hφ0
  have hψW : (rtRootIso P F A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl)).inv ψ = HomPf.mk W ψ₁ := by
    rw [hψ₁, HomPf.mk_map]; exact hψ0
  have hφmk : φ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
      (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hφW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) W φ₁)
  have hψmk : ψ = HomPf.mk ((pushIdx (F := F)
      (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
      (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj W) ψ₁ :=
    ((Iso.inv_hom_id_apply _ _).symm.trans
      (congrArg (ConcreteCategory.hom (rtRootIso P F A B
        (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) hψW)).trans
      (rtRootIso_hom_mk (F := F) A B (show r * r = r * r from rfl)
        (show r * r = r * r from rfl) W ψ₁)
  rw [hφmk, hψmk] at hbe hme
  have hbe' : BaseEquivalent P φ₁ ψ₁ :=
    base_eq_of_rootBase_eq (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ ψ₁ hbe
  have hme' : MetricallyEquivalent P φ₁ ψ₁ :=
    div_eq_of_rootDiv_eq (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
      ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁ ψ₁ hme
  have hφ₁s : IsPreStep P φ₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) φ₁).mp ?_
    rw [← hφmk]; exact hφs
  have hψ₁s : IsPreStep P ψ₁ := by
    refine (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushIdx (F := F)
        (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
        (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
        (by rw [rtLift_degFr, rtLift_degFr])).obj W) ψ₁).mp ?_
    rw [← hψmk]; exact hψs
  have hφ₁c : IsCoAngular P φ₁ :=
    prop_1_4_i P φ₁ (fun _ g => F.isotropicClosed g hW)
  have hψ₁c : IsCoAngular P ψ₁ :=
    prop_1_4_i P ψ₁ (fun _ g => F.isotropicClosed g hW)
  obtain ⟨α₀, hα₀U, hα₀e⟩ :=
    F.faithfulUpToUnits φ₁ ψ₁ hbe' hme' hφ₁c hφ₁s hψ₁c hψ₁s
  refine ⟨(rtRootIso P F B B (show r * r = r * r from rfl)
      (show r * r = r * r from rfl)).hom
    (HomPf.mk (idxMk (P := P) (F := F) W.hom.hom.2 W.hom.hom.2 hw2 hw2 rfl)
      ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2)), ?_, ?_⟩
  · refine ⟨?_, ?_⟩
    · rw [rtRootIso_hom_mk]
      exact (oTri_mk_diag (X := (⟨B, r⟩ : PfRootObj P F))
        (rtLift P F B (show r * r = r * r from rfl) ≫ W.hom.hom.2)
        (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hw2) _).mpr
        ⟨hα₀U.1.1, hα₀U.1.2⟩
    · refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
      rw [rtRootIso_hom_mk]
      exact pfRoot_isIso_mk _ _ ((CategoryTheory.isUnit_iff_isIso _).mp hα₀U.2)
  · have e1 : compPf P F (HomPf.mk W ψ₁)
        (HomPf.mk (idxMk (P := P) (F := F) W.hom.hom.2 W.hom.hom.2 hw2 hw2 rfl)
          ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2))
        = HomPf.mk W (ψ₁ ≫ ((α₀ : End W.right.obj.2) : W.right.obj.2 ⟶ W.right.obj.2)) :=
      compPf_mk (idxMk3 (F := F) W.hom.hom.1 W.hom.hom.2 W.hom.hom.2 hw1 hw2 hw2 hwd rfl)
        ψ₁ _
    show φ = compRoot P F ψ _
    unfold compRoot
    rw [hψW, Iso.hom_inv_id_apply]
    refine Eq.trans ?_ (congrArg (ConcreteCategory.hom (rtRootIso P F A B
      (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom) e1).symm
    refine Eq.trans ?_ (congrArg (fun t => ConcreteCategory.hom (rtRootIso P F A B
      (show r * r = r * r from rfl) (show r * r = r * r from rfl)).hom
      (HomPf.mk W t)) hα₀e)
    rw [← hφW, Iso.inv_hom_id_apply]

/-! ## ★35. (vi) を一般の根へ -/

/-- ★同型で挟んでも `Div` の一致は保たれる。 -/
theorem metricallyEquivalent_conj {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3) {X X' Y Y' : C3}
    (a : X' ⟶ X) (b : Y ⟶ Y') (ha : IsIso a) (hb : IsIso b) {f g : X ⟶ Y}
    (h : MetricallyEquivalent Q f g) :
    MetricallyEquivalent Q (a ≫ f ≫ b) (a ≫ g ≫ b) := by
  haveI := ha
  haveI := hb
  show Q.Div (a ≫ f ≫ b) = Q.Div (a ≫ g ≫ b)
  rw [Q.Div_comp, Q.Div_comp, Q.Div_comp, Q.Div_comp,
    show Q.Div a = 0 from isIsometric_of_isIso Q a,
    show Q.Div b = 0 from isIsometric_of_isIso Q b,
    show Q.Div f = Q.Div g from h]
  simp

/-- ★同型で挟んでも `Base` の一致は保たれる。 -/
theorem baseEquivalent_conj {C3 : Type u3} [Category.{v3} C3]
    {Φ3 : MonoidOn.{v, u, w} D} (Q : PreFrobenioid C3 Φ3) {X X' Y Y' : C3}
    (a : X' ⟶ X) (b : Y ⟶ Y') {f g : X ⟶ Y} (h : BaseEquivalent Q f g) :
    BaseEquivalent Q (a ≫ f ≫ b) (a ≫ g ≫ b) := by
  show Q.Base (a ≫ f ≫ b) = Q.Base (a ≫ g ≫ b)
  rw [Q.Base_comp, Q.Base_comp, Q.Base_comp, Q.Base_comp,
    show Q.Base f = Q.Base g from h]

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (vi)** —— `𝒞^pf` 版。 -/
theorem pfRoot_faithfulUpToUnits (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F} (φ ψ : X ⟶ Y)
    (hbe : BaseEquivalent (pfRootPre P F) φ ψ)
    (hme : MetricallyEquivalent (pfRootPre P F) φ ψ)
    (_hφc : IsCoAngular (pfRootPre P F) φ) (hφs : IsPreStep (pfRootPre P F) φ)
    (_hψc : IsCoAngular (pfRootPre P F) ψ) (hψs : IsPreStep (pfRootPre P F) ψ) :
    ∃ α : End Y, α ∈ OTimes (pfRootPre P F) Y ∧ φ = (ψ ≫ ((α : Y ⟶ Y))) := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hφ' : IsPreStep (pfRootPre P F) (inv eX ≫ φ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hφs (isPreStep_of_isIso (pfRootPre P F) eY))
  have hψ' : IsPreStep (pfRootPre P F) (inv eX ≫ ψ ≫ eY) :=
    IsPreStep.comp (pfRootPre P F) (isPreStep_of_isIso (pfRootPre P F) (inv eX))
      (IsPreStep.comp (pfRootPre P F) hψs (isPreStep_of_isIso (pfRootPre P F) eY))
  obtain ⟨α', hα'U, hα'e⟩ := pfRoot_faithfulUpToUnits_sameRoot (F := F) hfi
    (inv eX ≫ φ ≫ eY) (inv eX ≫ ψ ≫ eY)
    (baseEquivalent_conj (pfRootPre P F) (inv eX) eY hbe)
    (metricallyEquivalent_conj (pfRootPre P F) (inv eX) eY inferInstance inferInstance hme)
    hφ' hψ'
  refine ⟨eY ≫ ((α' : _ ⟶ _) ≫ inv eY), ⟨?_, ?_⟩, ?_⟩
  · exact (oTri_conj (F := F) (inv eY) eY (IsIso.inv_hom_id eY) (IsIso.hom_inv_id eY)
      α' hα'U.1)
  · refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
    haveI : IsIso ((α' : _ ⟶ _)) := (CategoryTheory.isUnit_iff_isIso _).mp hα'U.2
    show IsIso (eY ≫ ((α' : _ ⟶ _) ≫ inv eY))
    infer_instance
  · have h := hα'e
    simp only [Category.assoc] at h ⊢
    have h2 : eX ≫ inv eX ≫ φ ≫ eY = eX ≫ inv eX ≫ ψ ≫ eY ≫ (α' : _ ⟶ _) := by
      rw [h]
    simp only [IsIso.hom_inv_id_assoc] at h2
    have h3 : (φ ≫ eY) ≫ inv eY = (ψ ≫ eY ≫ (α' : _ ⟶ _)) ≫ inv eY :=
      congrArg (fun t => t ≫ inv eY) h2
    simpa only [Category.assoc, IsIso.hom_inv_id, Category.comp_id] using h3

/-! ## ★36. ★★★★`k` 乗根の高さへの持ち上げ `Λ_k`

★★**動機**(2026-08-17): 残る 4 条(`pullBackLB` / `arbFactor` / `arbFactorUniq` /
`plBkEquiv`)はどれも「`𝒞` の射を `𝒞^pf` へ**根つきで**送る」道具を要る。
★`toRootHom` は根 `1` にしか送れないが、必要なのは**任意の根 `k`** である。

★★**構成**: `rtExt A k : A ⟶ A^{(k)}` と `rtExt B k : B ⟶ B^{(k)}` は
同じ次数 `k` の Frobenius 型射なので、`Proposition 1.10, (i)` の遷移
(`frobTransport`)が `φ : A ⟶ B` を `A^{(k)} ⟶ B^{(k)}` へ一意に持ち上げる。 -/

variable {P F} in
/-- ★★`φ : A ⟶ B` を **`k` 乗根の高さへ持ち上げた射** `A^{(k)} ⟶ B^{(k)}`。 -/
noncomputable def rtMap (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    rtObj P F A k ⟶ rtObj P F B k :=
  frobTransport (F := F) (rtExt P F A k) (rtExt_frobType P F A k)
    (rtExt P F B k) (rtExt_frobType P F B k)
    (by rw [rtExt_degFr, rtExt_degFr]) φ

variable {P F} in
/-- ★持ち上げの特徴づけ(四角形)。 -/
theorem rtMap_spec (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    φ ≫ rtExt P F B k = rtExt P F A k ≫ rtMap (F := F) k φ :=
  frobTransport_spec _ _ _ _ _ φ

variable {P F} in
/-- ★持ち上げの一意性。 -/
theorem rtMap_eq (k : ℕ+) {A B : C} (φ : A ⟶ B) (ψ : rtObj P F A k ⟶ rtObj P F B k)
    (h : φ ≫ rtExt P F B k = rtExt P F A k ≫ ψ) : rtMap (F := F) k φ = ψ :=
  frobTransport_eq _ _ _ _ _ φ ψ h

variable {P F} in
/-- ★★**次数は変わらない**。 -/
theorem rtMap_degFr (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    P.degFr (rtMap (F := F) k φ) = P.degFr φ := by
  have h := congrArg P.degFr (rtMap_spec (F := F) k φ)
  rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr] at h
  have h2 : P.degFr φ * k = P.degFr (rtMap (F := F) k φ) * k := (mul_comm _ _).trans h
  exact (mul_right_cancel h2).symm

variable {P F} in
/-- ★★**零因子は `k` 倍で写る** —— `(Base ext)^*(Div (Λφ)) = k • Div φ`。 -/
theorem rtMap_Div (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    Φ.map (P.Base (rtExt P F A k)) (P.Div (rtMap (F := F) k φ))
      = ((k : ℕ+) : ℕ) • P.Div φ := by
  have h := congrArg P.Div (rtMap_spec (F := F) k φ)
  rw [P.Div_comp, P.Div_comp, show P.Div (rtExt P F B k) = 0 from (rtExt_frobType P F B k).1.2,
    show P.Div (rtExt P F A k) = 0 from (rtExt_frobType P F A k).1.2,
    rtExt_degFr, map_zero, zero_add, smul_zero, add_zero] at h
  exact h.symm

variable {P F} in
/-- ★★**等長性は行き来する**。 -/
theorem rtMap_isIsometric_iff (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    IsIsometric P (rtMap (F := F) k φ) ↔ IsIsometric P φ := by
  constructor
  · intro h
    have h1 : ((k : ℕ+) : ℕ) • P.Div φ = 0 := by
      rw [← rtMap_Div (F := F) k φ, show P.Div (rtMap (F := F) k φ) = 0 from h, map_zero]
    exact nsmul_eq_zero_of_isSharp (P.divisorial _).2 h1
  · intro h
    refine Φ.map_injective (P.Base (rtExt P F A k)) ?_
    rw [rtMap_Div, show P.Div φ = 0 from h, smul_zero, map_zero]

variable {P F} in
/-- ★★**底の同型性も行き来する**。 -/
theorem rtMap_isBaseIso_iff (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    IsBaseIsomorphism P (rtMap (F := F) k φ) ↔ IsBaseIsomorphism P φ := by
  haveI hA : IsIso (P.Base (rtExt P F A k)) := (rtExt_frobType P F A k).2
  haveI hB : IsIso (P.Base (rtExt P F B k)) := (rtExt_frobType P F B k).2
  refine (isIso_iff_of_sq (P.Base φ) (P.Base (rtExt P F B k))
    (P.Base (rtExt P F A k)) (P.Base (rtMap (F := F) k φ)) ?_).symm
  rw [← P.Base_comp, ← P.Base_comp, rtMap_spec]

variable {P F} in
/-- ★恒等射の持ち上げは恒等射。 -/
@[simp] theorem rtMap_id (k : ℕ+) (A : C) :
    rtMap (F := F) k (𝟙 A) = 𝟙 (rtObj P F A k) :=
  rtMap_eq k _ _ (by rw [Category.id_comp, Category.comp_id])

variable {P F} in
/-- ★★**持ち上げは合成を保つ**。 -/
theorem rtMap_comp (k : ℕ+) {A B E : C} (φ : A ⟶ B) (ψ : B ⟶ E) :
    rtMap (F := F) k (φ ≫ ψ) = rtMap (F := F) k φ ≫ rtMap (F := F) k ψ := by
  refine rtMap_eq k _ _ ?_
  rw [Category.assoc, rtMap_spec, ← Category.assoc, rtMap_spec, Category.assoc]

variable {P F} in
/-- ★同型の持ち上げは同型。 -/
theorem rtMap_isIso (k : ℕ+) {A B : C} (u : A ⟶ B) [IsIso u] :
    IsIso (rtMap (F := F) k u) :=
  ⟨rtMap (F := F) k (inv u),
    by rw [← rtMap_comp, IsIso.hom_inv_id, rtMap_id],
    by rw [← rtMap_comp, IsIso.inv_hom_id, rtMap_id]⟩

/-! ### ★`Λ_k` —— `𝒞` の射から `𝒞^pf` の根 `k` の射へ -/

variable {P F} in
/-- ★★★**`φ : A ⟶ B` の根 `k` での像** `(A,k) ⟶ (B,k)`。 -/
noncomputable def lamHom (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (⟨A, k⟩ : PfRootObj P F) ⟶ (⟨B, k⟩ : PfRootObj P F) :=
  HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)

variable {P F} in
@[simp] theorem lamHom_degFr (k : ℕ+) {A B : C} (φ : A ⟶ B) :
    (pfRootPre P F).degFr (lamHom (F := F) k φ) = P.degFr φ := by
  show rootDeg (show HomRoot P F ⟨A, k⟩ ⟨B, k⟩ from
    HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)) = _
  rw [rootDeg_mk]
  exact rtMap_degFr (F := F) k φ

variable {P F} in
theorem lamHom_isLinear (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsLinear P φ) :
    IsLinear (pfRootPre P F) (lamHom (F := F) k φ) := by
  show (pfRootPre P F).degFr (lamHom (F := F) k φ) = 1
  rw [lamHom_degFr]; exact h

variable {P F} in
theorem lamHom_isIsometric (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsIsometric P φ) :
    IsIsometric (pfRootPre P F) (lamHom (F := F) k φ) :=
  (isIsometric_mk_iff (X := (⟨A, k⟩ : PfRootObj P F)) (Y := (⟨B, k⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)).mpr
    ((rtMap_isIsometric_iff (F := F) k φ).mpr h)

variable {P F} in
theorem lamHom_isBaseIso (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism (pfRootPre P F) (lamHom (F := F) k φ) :=
  (isBaseIsomorphism_mk_iff (X := (⟨A, k⟩ : PfRootObj P F)) (Y := (⟨B, k⟩ : PfRootObj P F))
    (idxOne P F (rtObj P F A k) (rtObj P F B k)) (rtMap (F := F) k φ)).mpr
    ((rtMap_isBaseIso_iff (F := F) k φ).mpr h)

variable {P F} in
theorem lamHom_isPreStep (k : ℕ+) {A B : C} (φ : A ⟶ B) (h : IsPreStep P φ) :
    IsPreStep (pfRootPre P F) (lamHom (F := F) k φ) :=
  ⟨lamHom_isLinear k φ h.1, lamHom_isBaseIso k φ h.2⟩

variable {P F} in
/-- ★★同型の像は同型。 -/
theorem lamHom_isIso (k : ℕ+) {A B : C} (u : A ⟶ B) (hu : IsIso u) :
    IsIso (lamHom (F := F) k u) := by
  haveI := hu
  exact pfRoot_isIso_mk _ _ (rtMap_isIso (F := F) k u)

variable {P F} in
/-- ★★**Frobenius 型の像は Frobenius 型**(`𝒞^pf` が isotropic 型のとき)。 -/
theorem lamHom_isFrobeniusType (hfi : IsOfFrobeniusIsotropicType P) (k : ℕ+)
    {A B : C} (φ : A ⟶ B) (h : IsFrobeniusType P φ) :
    IsFrobeniusType (pfRootPre P F) (lamHom (F := F) k φ) :=
  ⟨⟨pfRoot_isCoAngular hfi _, lamHom_isIsometric k φ h.1.2⟩, lamHom_isBaseIso k φ h.2⟩

/-! ## ★37. 根を上げる同型の一般形と、根が等しいときの合成則 -/

variable {P F} in
/-- ★★**任意の Frobenius 型射で根を上げられる** —— `pfRoot_exists_iso_root` の一般形。

★`z : A ⟶ U` が次数 `d` の Frobenius 型射なら `(A, k) ≅ (U, k*d)`。
★`pfRoot_exists_iso_root` は `z = rtExt A d` の場合であった。 -/
theorem pfRoot_iso_of_frobType {A U : C} (z : A ⟶ U) (hz : IsFrobeniusType P z)
    {d k t : ℕ+} (hzd : P.degFr z = d) (ht : t = k * d) :
    ∃ e : (⟨A, k⟩ : PfRootObj P F) ⟶ ⟨U, t⟩, IsIso e := by
  obtain ⟨e₀, he₀⟩ := pfRoot_exists_iso_root (F := F) A k d t ht
  obtain ⟨θ, hθ, -⟩ := F.frobDegUniq A U (rtObj P F A d) z (rtExt P F A d) hz
    (rtExt_frobType P F A d) (by rw [hzd, rtExt_degFr])
  haveI := hθ
  haveI := he₀
  haveI : IsIso (lamHom (F := F) t (@inv _ _ _ _ θ hθ)) :=
    lamHom_isIso t _ (@IsIso.inv_isIso _ _ _ _ θ hθ)
  exact ⟨e₀ ≫ lamHom (F := F) t (@inv _ _ _ _ θ hθ), inferInstance⟩

/-- ★**根 `r` の対象どうしの `Hom` を `r·r` 乗の高さで書くときの添字の押し出し**。 -/
noncomputable def pushRt (A B : C) (r : ℕ+) :
    IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)) ⥤
      IdxPf P F (rtObj P F A r) (rtObj P F B r) :=
  pushIdx (F := F) (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
    (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
    (by rw [rtLift_degFr, rtLift_degFr])

/-- ★**根が等しいときの、`Hom` の同一視** ——
`Hom^pf(A^{(r·r)}, B^{(r·r)}) ≅ ((A,r) ⟶ (B,r))`。 -/
noncomputable def rtJ (A B : C) (r : ℕ+) :
    HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r))
      ≅ HomPf P F (rtObj P F A r) (rtObj P F B r) :=
  rtRootIso P F A B (show r * r = r * r from rfl) (show r * r = r * r from rfl)

variable {P F} in
/-- ★★**根が等しいときの合成則** —— `compRoot` は `rtJ` で共役した `compPf` である。

★★3 つの `rtRootIso` の添字はすべて `(dA, dB, e, tA, tB) = (r, r, r, r·r, r·r)` で、
命題の証明部分は証明無関係で潰れる。★これが (iii)(c) や (vi) で効いた仕掛けと同じもの。 -/
theorem compRoot_sameRoot {A B E : C} {r : ℕ+}
    (f : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) (g : (⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) :
    f ≫ g = (rtJ P F A E r).hom
      (compPf P F ((rtJ P F A B r).inv f) ((rtJ P F B E r).inv g)) := rfl

variable {P F} in
/-- ★★**`rtJ` は合成を移す**。 -/
theorem rtJ_comp {A B E : C} {r : ℕ+}
    (x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (y : HomPf P F (rtObj P F B (r * r)) (rtObj P F E (r * r))) :
    (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from (rtJ P F A B r).hom x)
        ≫ (show ((⟨B, r⟩ : PfRootObj P F) ⟶ ⟨E, r⟩) from (rtJ P F B E r).hom y)
      = (rtJ P F A E r).hom (compPf P F x y) := by
  rw [compRoot_sameRoot, Iso.hom_inv_id_apply, Iso.hom_inv_id_apply]

variable {P F} in
/-- ★`rtJ` の代表元での計算則。 -/
theorem rtJ_hom_mk {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (rtJ P F A B r).hom (HomPf.mk V φ)
      = HomPf.mk ((pushRt P F A B r).obj V) φ :=
  rtRootIso_hom_mk (F := F) A B _ _ V φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 次数。 -/
theorem rtJ_degFr {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    (pfRootPre P F).degFr
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      = P.degFr φ := by
  rw [rtJ_hom_mk]
  exact rootDeg_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
    ((pushRt P F A B r).obj V) φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 等長性。 -/
theorem rtJ_isIsometric_iff {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    IsIsometric (pfRootPre P F)
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      ↔ IsIsometric P φ := by
  rw [rtJ_hom_mk]
  exact isIsometric_mk_iff (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F))
    ((pushRt P F A B r).obj V) φ

variable {P F} in
/-- ★★**根が等しいときの辞書** —— 底の同型性。 -/
theorem rtJ_isBaseIso_iff {A B : C} (r : ℕ+)
    (V : IdxPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)))
    (φ : V.right.obj.1 ⟶ V.right.obj.2) :
    IsBaseIsomorphism (pfRootPre P F)
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from
          (rtJ P F A B r).hom (HomPf.mk V φ))
      ↔ IsBaseIsomorphism P φ := by
  rw [rtJ_hom_mk]
  exact isBaseIsomorphism_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
    (Y := (⟨B, r⟩ : PfRootObj P F)) ((pushRt P F A B r).obj V) φ

/-! ## ★38. ★★★★`Hom^pf` の Frobenius 分解

★★**これが `pullBackLB` / `arbFactor` の共通の足場**である。
`𝒞` の `arbFactor` を代表元に当て、**Frobenius 型の部分を添字の脚へ吸収する**。
★吸収先は **`A₀` の次数 `n` の標準拡大 `A₀^{(n)}`** に取る —— そうすると
`𝒞^pf` 側で中間対象を `(A^{(n)}, r)` の形に書ける。 -/

variable {P F} in
/-- ★添字の始対象から `idxMk` への遷移射。 -/
noncomputable def idxOneHom {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    idxOne P F A B ⟶ idxMk (P := P) (F := F) a b ha hb hd :=
  Under.homMk (show (⟨(A, B)⟩ : BiFr P F) ⟶ (⟨(A', B')⟩ : BiFr P F) from
    ⟨(a, b), ha, hb, hd⟩)
    (WideSubcategory.hom_ext _ (Prod.ext (Category.id_comp _) (Category.id_comp _)))

variable {P F} in
/-- ★★**`𝒞` の射の像を、脚を伸ばした代表元で書く**。 -/
theorem toHomPf_eq_mk {A B A' B' : C} (a : A ⟶ A') (b : B ⟶ B')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A ⟶ B) (ψ : A' ⟶ B') (h : φ ≫ b = a ≫ ψ) :
    toHomPf (F := F) φ = HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) ψ := by
  have h1 : idxTransport P F (idxOneHom (F := F) a b ha hb hd) φ = ψ :=
    frobTransport_eq a ha b hb hd φ ψ h
  rw [← h1]
  exact (HomPf.mk_map (idxOneHom (F := F) a b ha hb hd) φ).symm

variable {P F} in
/-- ★★★★**`Hom^pf` の Frobenius 分解** ——
`x : A₀ ⟶ B₀`(`Hom^pf` の元、次数 `n`)は
**標準の Frobenius 拡大 `A₀ → A₀^{(n)}` の像**と**次数 `1` の射**の合成になる。 -/
theorem homPf_frobSplit {A₀ B₀ : C} (x : HomPf P F A₀ B₀) {n : ℕ+} (hn : pfDeg x = n) :
    ∃ z : HomPf P F (rtObj P F A₀ n) B₀,
      pfDeg z = 1 ∧ compPf P F (toHomPf (F := F) (rtExt P F A₀ n)) z = x := by
  obtain ⟨V, χ, hx⟩ := HomPf.exists_rep (P := P) (F := F) x
  have hχ : P.degFr χ = n := by rw [← hn, ← hx, pfDeg_mk]; rfl
  obtain ⟨hv1, hv2, hvd⟩ := V.hom.property
  obtain ⟨G, H, γ, β, a, hfac, hγ, hβ, ha⟩ := F.arbFactor χ
  have haLin : P.degFr a = 1 := (F.pullBackLB a ha).2
  have hρ : P.degFr (β ≫ a) = 1 := by
    rw [P.degFr_comp, haLin, show P.degFr β = 1 from hβ.1, mul_one]
  have hγn : P.degFr γ = n := by
    have h0 := congrArg P.degFr hfac
    rw [P.degFr_comp, hρ, one_mul, hχ] at h0
    exact h0.symm
  -- ★`A₀^{(n)}` の `d` 次拡大と `G` を合わせる
  have h1 : IsFrobeniusType P (rtExt P F A₀ n ≫
      rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1)) :=
    IsFrobeniusType.comp P F (rtExt_frobType P F A₀ n)
      (rtExt_frobType P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
  have h2 : IsFrobeniusType P (V.hom.hom.1 ≫ γ) := IsFrobeniusType.comp P F hv1 hγ
  have hdeg : P.degFr (rtExt P F A₀ n ≫
        rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
      = P.degFr (V.hom.hom.1 ≫ γ) := by
    rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr, hγn, mul_comm]
  obtain ⟨θ, hθ, hθe⟩ := F.frobDegUniq A₀ _ G _ _ h1 h2 hdeg
  haveI := hθ
  obtain ⟨a₂, ha₂eq⟩ : ∃ t : rtObj P F A₀ n ⟶ G,
      t = rtExt P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1) ≫ θ := ⟨_, rfl⟩
  have ha₂F : IsFrobeniusType P a₂ := by
    rw [ha₂eq]
    exact IsFrobeniusType.comp P F
      (rtExt_frobType P F (rtObj P F A₀ n) (P.degFr V.hom.hom.1))
      (isFrobeniusType_of_isIso P θ)
  have ha₂d : P.degFr a₂ = P.degFr V.hom.hom.1 := by
    rw [ha₂eq, P.degFr_comp, show P.degFr θ = 1 from isLinear_of_isIso P θ, rtExt_degFr,
      one_mul]
  have hsq : rtExt P F A₀ n ≫ a₂ = V.hom.hom.1 ≫ γ := by
    rw [ha₂eq, ← Category.assoc]
    exact hθe
  -- ★2 つの代表元
  have hy : toHomPf (F := F) (rtExt P F A₀ n)
      = HomPf.mk (idxMk (P := P) (F := F) V.hom.hom.1 a₂ hv1 ha₂F ha₂d.symm) γ :=
    toHomPf_eq_mk V.hom.hom.1 a₂ hv1 ha₂F ha₂d.symm _ γ hsq
  refine ⟨HomPf.mk (idxMk (P := P) (F := F) a₂ V.hom.hom.2 ha₂F hv2 (ha₂d.trans hvd))
      (β ≫ a), ?_, ?_⟩
  · rw [pfDeg_mk]; exact hρ
  · rw [hy]
    exact (compPf_mk_pair (P := P) (F := F) V.hom.hom.1 a₂ V.hom.hom.2
      hv1 ha₂F hv2 ha₂d.symm (ha₂d.trans hvd) γ (β ≫ a)).trans
      ((congrArg (HomPf.mk V) hfac.symm).trans hx)

end ABC3.Found.FrdI
