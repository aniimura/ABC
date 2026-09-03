/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop32Frob.Functor

/-!
# Prop32Frob —— `Λ_k` は関手・(iv)(b)・根が違う合成の組み立て

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-! ## ★39. ★★★★`Λ_k` は関手である

★★`lamHom` を `rtJ` で書き直すと、合成則が `rtJ_comp` ＋ `toHomPf_comp` ＋
`rtMap_comp` の 3 本から出る。 -/

variable {P F} in
/-- ★添字の始対象から「押し出した始対象」への遷移射。 -/
noncomputable def idxPushOneHom {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b) :
    idxOne P F A B ⟶ (pushIdx (F := F) a ha b hb hd).obj (idxOne P F A' B') :=
  Under.homMk (show (⟨(A, B)⟩ : BiFr P F) ⟶ (⟨(A', B')⟩ : BiFr P F) from
    ⟨(a, b), ha, hb, hd⟩)
    (WideSubcategory.hom_ext _ (Prod.ext
      ((Category.id_comp _).trans (Category.comp_id _).symm)
      ((Category.id_comp _).trans (Category.comp_id _).symm)))

variable {P F} in
/-- ★★**根の不変性で `𝒞` の射の像がどう写るか**。 -/
theorem toHomPf_rootIso {A B A' B' : C} (a : A ⟶ A') (ha : IsFrobeniusType P a)
    (b : B ⟶ B') (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (φ : A ⟶ B) (ψ : A' ⟶ B') (h : φ ≫ b = a ≫ ψ) :
    (rootIso (F := F) a ha b hb hd).hom (toHomPf (F := F) ψ) = toHomPf (F := F) φ := by
  have h1 : idxTransport P F (idxPushOneHom (F := F) a ha b hb hd) φ = ψ :=
    frobTransport_eq a ha b hb hd φ ψ h
  refine Eq.trans (rootIso_hom_mk a ha b hb hd (idxOne P F A' B') ψ) ?_
  rw [← h1]
  exact HomPf.mk_map (idxPushOneHom (F := F) a ha b hb hd) φ

variable {P F} in
/-- ★★**`rtMap` は 1 段上げても四角形で繋がる**。 -/
theorem rtMap_rtLift (r : ℕ+) {A B : C} (φ : A ⟶ B) :
    rtMap (F := F) r φ ≫ rtLift P F B (show r * r = r * r from rfl)
      = rtLift P F A (show r * r = r * r from rfl) ≫ rtMap (F := F) (r * r) φ := by
  haveI : Epi (rtExt P F A r) := P.totEpiC _ _ _
  refine (cancel_epi (rtExt P F A r)).mp ?_
  rw [← Category.assoc, ← rtMap_spec, Category.assoc, rtLift_ext, rtMap_spec,
    ← Category.assoc, rtLift_ext]

variable {P F} in
/-- ★★★**`Λ_r` を `rtJ` で書く**。 -/
theorem lamHom_eq_rtJ (r : ℕ+) {A B : C} (φ : A ⟶ B) :
    lamHom (F := F) r φ
      = (rtJ P F A B r).hom (toHomPf (F := F) (rtMap (F := F) (r * r) φ)) :=
  (toHomPf_rootIso (rtLift P F A (show r * r = r * r from rfl)) (rtLift_frobType P F A _)
    (rtLift P F B (show r * r = r * r from rfl)) (rtLift_frobType P F B _)
    (by rw [rtLift_degFr, rtLift_degFr]) (rtMap (F := F) r φ)
    (rtMap (F := F) (r * r) φ) (rtMap_rtLift (F := F) r φ)).symm

variable {P F} in
/-- ★★★**`Λ_k` は合成を保つ**。 -/
theorem lamHom_comp (k : ℕ+) {A B E : C} (φ : A ⟶ B) (ψ : B ⟶ E) :
    lamHom (F := F) k φ ≫ lamHom (F := F) k ψ = lamHom (F := F) k (φ ≫ ψ) := by
  rw [lamHom_eq_rtJ, lamHom_eq_rtJ, lamHom_eq_rtJ, rtJ_comp, ← toHomPf_comp, rtMap_comp]

variable {P F} in
/-- ★★★**`Λ_k` は恒等射を保つ**。 -/
theorem lamHom_id (k : ℕ+) (A : C) :
    lamHom (F := F) k (𝟙 A) = 𝟙 (⟨A, k⟩ : PfRootObj P F) := by
  show HomPf.mk (idxOne P F (rtObj P F A k) (rtObj P F A k)) (rtMap (F := F) k (𝟙 A))
    = idRoot P F ⟨A, k⟩
  rw [rtMap_id]
  rfl

/-! ## ★40. (iv)(b) の第 1 歩 —— pull-back は linear(根が等しい場合)

★★手は `𝒞^birat` の `birat_pullBack_repr` と同じ:
**`α` を「Frobenius 型 ≫ 次数 1」に分解し、次数 1 の側を `α` に沿って持ち上げる**と
`1 = degFr α · degFr g` になる。
★`𝒞^pf` では分解の中間対象を **`(A^{(n)}, r)`** に取れる(`homPf_frobSplit` の
吸収先が標準拡大 `A₀^{(n)}` だから)——ここが根つきでも閉じる理由である。 -/

variable {P F} in
/-- ★`rtJ` の像の次数(代表元を経由しない形)。 -/
theorem rtJ_degFr' {A B : C} (r : ℕ+)
    (x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r))) :
    (pfRootPre P F).degFr
        (show ((⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) from (rtJ P F A B r).hom x)
      = pfDeg x := by
  obtain ⟨V, φ, hφ⟩ := HomPf.exists_rep (P := P) (F := F) x
  rw [← hφ, rtJ_degFr, pfDeg_mk]
  rfl

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**pull-back は linear**(根が等しい場合)。 -/
theorem pfRoot_isLinear_of_pullBack_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hα : IsPullBack (pfRootPre P F) α) : IsLinear (pfRootPre P F) α := by
  obtain ⟨x, hx⟩ : ∃ x : HomPf P F (rtObj P F A (r * r)) (rtObj P F B (r * r)),
      (rtJ P F A B r).hom x = α := ⟨(rtJ P F A B r).inv α, Iso.inv_hom_id_apply _ _⟩
  obtain ⟨n, hn⟩ : ∃ n : ℕ+, pfDeg x = n := ⟨_, rfl⟩
  obtain ⟨z, hz1, hz2⟩ := homPf_frobSplit (P := P) (F := F) x hn
  -- ★`A^{(n)}` を `r·r` 乗したものと、`A^{(r·r)}` を `n` 乗したものを合わせる
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
  -- ★分解の 2 本
  obtain ⟨gam, hgam⟩ : ∃ t : rtObj P F A (r * r) ⟶ rtObj P F (rtObj P F A n) (r * r),
      t = rtExt P F (rtObj P F A (r * r)) n ≫ @inv _ _ _ _ w hw := ⟨_, rfl⟩
  obtain ⟨Gam, hGam⟩ : ∃ t : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨rtObj P F A n, r⟩,
      t = (rtJ P F A (rtObj P F A n) r).hom (toHomPf (F := F) gam) := ⟨_, rfl⟩
  obtain ⟨R, hR⟩ : ∃ t : (⟨rtObj P F A n, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩,
      t = (rtJ P F (rtObj P F A n) B r).hom
        (compPf P F (toHomPf (F := F) w) z) := ⟨_, rfl⟩
  -- ★合成すると `α`
  have hcomp : Gam ≫ R = α := by
    rw [hGam, hR, rtJ_comp, ← compPf_assoc, ← toHomPf_comp, hgam, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id, hz2, hx]
  -- ★`R` は次数 1
  have hRdeg : (pfRootPre P F).degFr R = 1 := by
    rw [hR, rtJ_degFr', pfDeg_comp, hz1, pfDeg_toHomPf,
      show P.degFr w = 1 from isLinear_of_isIso P w, mul_one]
  -- ★`Gam` の底は同型
  have hGamBase : IsIso ((pfRootPre P F).Base Gam) := by
    have hb : IsBaseIsomorphism P gam := by
      haveI : IsIso (@inv _ _ _ _ w hw) := @IsIso.inv_isIso _ _ _ _ w hw
      rw [hgam]
      exact (IsFrobeniusType.comp P F (rtExt_frobType P F (rtObj P F A (r * r)) n)
        (isFrobeniusType_of_isIso P (@inv _ _ _ _ w hw))).2
    have := (rtJ_isBaseIso_iff (F := F) (A := A) (B := rtObj P F A n) r
      (idxOne P F (rtObj P F A (r * r)) (rtObj P F (rtObj P F A n) (r * r))) gam).mpr hb
    rw [hGam]
    exact this
  -- ★持ち上げ
  have hbase : (pfRootPre P F).Base R
      = @inv _ _ _ _ ((pfRootPre P F).Base Gam) hGamBase ≫ (pfRootPre P F).Base α := by
    refine eq_inv_comp_of ((pfRootPre P F).Base Gam) hGamBase _ _ ?_
    rw [← (pfRootPre P F).Base_comp, hcomp]
  obtain ⟨g, ⟨hg1, -⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα
    (⟨rtObj P F A n, r⟩ : PfRootObj P F) R _ hbase
  have hmul : (pfRootPre P F).degFr α * (pfRootPre P F).degFr g = 1 := by
    rw [← hRdeg, ← hg1, (pfRootPre P F).degFr_comp]
  have hcoe : (((pfRootPre P F).degFr α : ℕ+) : ℕ)
      * (((pfRootPre P F).degFr g : ℕ+) : ℕ) = 1 := by
    rw [← PNat.mul_coe, hmul]; rfl
  show (pfRootPre P F).degFr α = 1
  exact PNat.coe_injective (Nat.dvd_one.mp ⟨_, hcoe.symm⟩)

/-! ## ★41. ★★★★根が違ってもよい合成の**組み立て**規則

★`compRoot_rep` は「合成を分解する」向きだった。★ここでは**逆向き** ——
3 脚の添字と 2 本の `𝒞` の射から、`𝒞^pf` の 2 射とその合成を**作る**。
★★これがあれば `arbFactor` の中間対象を根の違う対象に取れる。 -/

variable {P F} in
/-- ★★★★**3 脚の添字から `𝒞^pf` の合成を組み立てる**。 -/
theorem compRoot_mk3 {X Y Z : PfRootObj P F} {A' B' E' : C}
    (a : rtObj P F X.obj (Z.root * Y.root) ⟶ A')
    (b : rtObj P F Y.obj (Z.root * X.root) ⟶ B')
    (e : rtObj P F Z.obj (Y.root * X.root) ⟶ E')
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (he : IsFrobeniusType P e)
    (hab : P.degFr a = P.degFr b) (hbe : P.degFr b = P.degFr e)
    (φ : A' ⟶ B') (ψ : B' ⟶ E') :
    (show (X ⟶ Y) from (rtRootIso P F X.obj Y.obj
        (show Z.root * Y.root = Z.root * Y.root from rfl)
        (show Z.root * X.root = Z.root * X.root from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) a b ha hb hab) φ))
      ≫ (show (Y ⟶ Z) from (rtRootIso P F Y.obj Z.obj
        (show Z.root * X.root = X.root * Z.root from mul_comm _ _)
        (show Y.root * X.root = X.root * Y.root from mul_comm _ _)).hom
      (HomPf.mk (idxMk (P := P) (F := F) b e hb he hbe) ψ))
      = (rtRootIso P F X.obj Z.obj
          (show Z.root * Y.root = Y.root * Z.root from mul_comm _ _)
          (show Y.root * X.root = Y.root * X.root from rfl)).hom
        (HomPf.mk (idxMk (P := P) (F := F) a e ha he (hab.trans hbe)) (φ ≫ ψ)) := by
  show compRoot P F _ _ = _
  unfold compRoot
  rw [Iso.hom_inv_id_apply, Iso.hom_inv_id_apply, compPf_mk_pair]

/-! ## ★42. 代表元の脚を `k` 乗ぶん伸ばす -/

variable {P F} in
/-- ★添字の脚を `k` 乗ぶん伸ばす遷移射。 -/
noncomputable def idxExtHom {A B : C} (V : IdxPf P F A B) (k : ℕ+) :
    V ⟶ idxMk (P := P) (F := F)
      (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k)
      (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k)
      (IsFrobeniusType.comp P F V.hom.property.1 (rtExt_frobType P F _ k))
      (IsFrobeniusType.comp P F V.hom.property.2.1 (rtExt_frobType P F _ k))
      (by rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr]
          exact congrArg (fun m => k * m) V.hom.property.2.2) :=
  Under.homMk (show V.right ⟶ (⟨(rtObj P F V.right.obj.1 k,
      rtObj P F V.right.obj.2 k)⟩ : BiFr P F) from
    ⟨(rtExt P F V.right.obj.1 k, rtExt P F V.right.obj.2 k),
      rtExt_frobType P F _ k, rtExt_frobType P F _ k,
      by rw [rtExt_degFr, rtExt_degFr]⟩)
    (WideSubcategory.hom_ext _ (Prod.ext rfl rfl))

variable {P F} in
/-- ★★**脚を伸ばした代表元** —— 射は `rtMap` で持ち上がる。 -/
theorem mk_rtExt {A B : C} (V : IdxPf P F A B) (k : ℕ+)
    (χ : V.right.obj.1 ⟶ V.right.obj.2) :
    HomPf.mk (idxMk (P := P) (F := F)
      (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k)
      (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k)
      (IsFrobeniusType.comp P F V.hom.property.1 (rtExt_frobType P F _ k))
      (IsFrobeniusType.comp P F V.hom.property.2.1 (rtExt_frobType P F _ k))
      (by rw [P.degFr_comp, P.degFr_comp, rtExt_degFr, rtExt_degFr]
          exact congrArg (fun m => k * m) V.hom.property.2.2))
      (rtMap (F := F) k χ)
      = HomPf.mk V χ := by
  have h1 : idxTransport P F (idxExtHom (F := F) V k) χ = rtMap (F := F) k χ :=
    frobTransport_eq _ _ _ _ _ χ _ (rtMap_spec (F := F) k χ)
  rw [← h1]
  exact HomPf.mk_map (idxExtHom (F := F) V k) χ

variable {P F} in
/-- ★★**脚を伸ばした代表元(伸ばし先を与える形)**。 -/
theorem mk_rtExt_gen {A B : C} (V : IdxPf P F A B) (k : ℕ+)
    (a : A ⟶ rtObj P F V.right.obj.1 k) (b : B ⟶ rtObj P F V.right.obj.2 k)
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b) (hd : P.degFr a = P.degFr b)
    (hva : V.hom.hom.1 ≫ rtExt P F V.right.obj.1 k = a)
    (hvb : V.hom.hom.2 ≫ rtExt P F V.right.obj.2 k = b)
    (χ : V.right.obj.1 ⟶ V.right.obj.2) :
    HomPf.mk (idxMk (P := P) (F := F) a b ha hb hd) (rtMap (F := F) k χ)
      = HomPf.mk V χ := by
  refine Eq.trans ?_ (HomPf.mk_map
    (Under.homMk (show V.right ⟶ (⟨(rtObj P F V.right.obj.1 k,
        rtObj P F V.right.obj.2 k)⟩ : BiFr P F) from
      ⟨(rtExt P F V.right.obj.1 k, rtExt P F V.right.obj.2 k),
        rtExt_frobType P F _ k, rtExt_frobType P F _ k,
        by rw [rtExt_degFr, rtExt_degFr]⟩)
      (WideSubcategory.hom_ext _ (Prod.ext hva hvb)) :
        V ⟶ idxMk (P := P) (F := F) a b ha hb hd) χ)
  exact congrArg (HomPf.mk _)
    (frobTransport_eq _ _ _ _ _ χ _ (rtMap_spec (F := F) k χ)).symm

/-! ## ★43. ★★★★linear な射は「pre-step ≫ 等長」に分解する

★★**中間対象の根が変わる**のが要点である。`X = (A,r)`、`Z = (B,r)` で
代表元の脚の次数が `c` のとき、中間対象は **`(Cc, r·r·c)`** に取る。
★`𝒞` の `Proposition 1.7, (iii)` の分解 `χ = α₁ ≫ γ₁` を
**`rtMap (r·r)` で丸ごと持ち上げ**、3 脚の添字の脚を**同型(次数 1)**に取る。 -/

variable {P F} in
/-- ★`idxMk` の脚が等しければ代表元も等しい。 -/
theorem mk_idxMk_congr {A B A' B' : C} {a a' : A ⟶ A'} {b b' : B ⟶ B'}
    (ha : a = a') (hb : b = b')
    {ha1 : IsFrobeniusType P a} {hb1 : IsFrobeniusType P b} {hd1 : P.degFr a = P.degFr b}
    {ha2 : IsFrobeniusType P a'} {hb2 : IsFrobeniusType P b'}
    {hd2 : P.degFr a' = P.degFr b'} (φ : A' ⟶ B') :
    HomPf.mk (idxMk (P := P) (F := F) a b ha1 hb1 hd1) φ
      = HomPf.mk (idxMk (P := P) (F := F) a' b' ha2 hb2 hd2) φ := by
  subst ha
  subst hb
  rfl

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★★**linear ⟹ 「pre-step ≫ 等長」**(根が等しい場合)。 -/
theorem pfRoot_linearSplit_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩) (hlin : IsLinear (pfRootPre P F) α) :
    ∃ (T : PfRootObj P F) (pp : (⟨A, r⟩ : PfRootObj P F) ⟶ T) (ll : T ⟶ ⟨B, r⟩),
      IsPreStep (pfRootPre P F) pp ∧ IsIsometric (pfRootPre P F) ll ∧ α = pp ≫ ll ∧
      ∃ (Z' : IdxPf P F (rtObj P F T.obj r) (rtObj P F B T.root))
        (χ' : Z'.right.obj.1 ⟶ Z'.right.obj.2),
        ll = HomPf.mk Z' χ' ∧ IsPullBack P χ' := by
  obtain ⟨V, χ, hV⟩ := HomPf.exists_rep (P := P) (F := F)
    (show HomPf P F (rtObj P F A r) (rtObj P F B r) from α)
  have hχlin : IsLinear P χ :=
    (rootDeg_mk (X := (⟨A, r⟩ : PfRootObj P F)) (Y := (⟨B, r⟩ : PfRootObj P F)) V χ).symm.trans
      (by rw [hV]; exact hlin)
  obtain ⟨Cc, α₁, γ₁, hfac, hα₁, hγ₁⟩ := (prop_1_7_iii_linear_factor P F χ).mp hχlin
  obtain ⟨hv1, hv2, hvd⟩ := V.hom.property
  obtain ⟨t, ht⟩ : ∃ t : ℕ+, t = r * r * P.degFr V.hom.hom.1 := ⟨_, rfl⟩
  -- ★2 本の脚(次数 1、すなわち同型)
  have hd1 : P.degFr (rtLift P F A (show r * t = t * r from mul_comm _ _))
      = P.degFr (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 (r * r)) := by
    rw [rtLift_degFr, P.degFr_comp, rtExt_degFr, ht]
  have hd3 : P.degFr (rtLift P F B (show t * r = t * r from rfl))
      = P.degFr (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 (r * r)) := by
    rw [rtLift_degFr, P.degFr_comp, rtExt_degFr, ht, hvd]
  obtain ⟨l₁, hl₁iso, hl₁⟩ := F.frobDegUniq (rtObj P F A r) _ _
    (rtLift P F A (show r * t = t * r from mul_comm _ _))
    (V.hom.hom.1 ≫ rtExt P F V.right.obj.1 (r * r))
    (rtLift_frobType P F A _)
    (IsFrobeniusType.comp P F hv1 (rtExt_frobType P F _ (r * r))) hd1
  obtain ⟨l₃, hl₃iso, hl₃⟩ := F.frobDegUniq (rtObj P F B r) _ _
    (rtLift P F B (show t * r = t * r from rfl))
    (V.hom.hom.2 ≫ rtExt P F V.right.obj.2 (r * r))
    (rtLift_frobType P F B _)
    (IsFrobeniusType.comp P F hv2 (rtExt_frobType P F _ (r * r))) hd3
  haveI := hl₁iso
  haveI := hl₃iso
  have hl₁F : IsFrobeniusType P l₁ := isFrobeniusType_of_isIso P l₁
  have hl₃F : IsFrobeniusType P l₃ := isFrobeniusType_of_isIso P l₃
  have hl₂F : IsFrobeniusType P (𝟙 (rtObj P F Cc (r * r))) :=
    isFrobeniusType_of_isIso P (𝟙 _)
  have hab : P.degFr l₁ = P.degFr (𝟙 (rtObj P F Cc (r * r))) := by
    rw [show P.degFr l₁ = 1 from isLinear_of_isIso P l₁, P.degFr_id]
  have hbe : P.degFr (𝟙 (rtObj P F Cc (r * r))) = P.degFr l₃ := by
    rw [show P.degFr l₃ = 1 from isLinear_of_isIso P l₃, P.degFr_id]
  refine ⟨⟨Cc, t⟩, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · exact (rtRootIso P F A Cc (show r * t = r * t from rfl)
      (show r * r = r * r from rfl)).hom
      (HomPf.mk (idxMk (P := P) (F := F) l₁ (𝟙 (rtObj P F Cc (r * r))) hl₁F hl₂F hab)
        (rtMap (F := F) (r * r) α₁))
  · exact (rtRootIso P F Cc B (show r * r = r * r from mul_comm _ _)
      (show t * r = r * t from mul_comm _ _)).hom
      (HomPf.mk (idxMk (P := P) (F := F) (𝟙 (rtObj P F Cc (r * r))) l₃ hl₂F hl₃F hbe)
        (rtMap (F := F) (r * r) γ₁))
  · rw [rtRootIso_hom_mk]
    exact (isPreStep_mk_iff (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨Cc, t⟩ : PfRootObj P F)) _ _).mpr
      ⟨(rtMap_degFr (F := F) (r * r) α₁).trans hα₁.1,
        (rtMap_isBaseIso_iff (F := F) (r * r) α₁).mpr hα₁.2⟩
  · rw [rtRootIso_hom_mk]
    exact (isIsometric_mk_iff (X := (⟨Cc, t⟩ : PfRootObj P F))
      (Y := (⟨B, r⟩ : PfRootObj P F)) _ _).mpr
      ((rtMap_isIsometric_iff (F := F) (r * r) γ₁).mpr (F.pullBackLB γ₁ hγ₁).1.2)
  · have hcc : rtMap (F := F) (r * r) α₁ ≫ rtMap (F := F) (r * r) γ₁
        = rtMap (F := F) (r * r) χ := by rw [← rtMap_comp, ← hfac]
    refine Eq.trans ?_ (compRoot_mk3 (X := (⟨A, r⟩ : PfRootObj P F))
      (Y := (⟨Cc, t⟩ : PfRootObj P F)) (Z := (⟨B, r⟩ : PfRootObj P F))
      l₁ (𝟙 (rtObj P F Cc (r * r))) l₃ hl₁F hl₂F hl₃F hab hbe
      (rtMap (F := F) (r * r) α₁) (rtMap (F := F) (r * r) γ₁)).symm
    refine Eq.trans ?_ (congrArg
      (fun u : rtObj P F V.right.obj.1 (r * r) ⟶ rtObj P F V.right.obj.2 (r * r) =>
        (rtRootIso P F A B (show r * t = t * r from mul_comm _ _)
          (show t * r = t * r from rfl)).hom
          (HomPf.mk (idxMk (P := P) (F := F) l₁ l₃ hl₁F hl₃F (hab.trans hbe)) u))
      hcc.symm)
    rw [rtRootIso_hom_mk]
    symm
    exact (mk_rtExt_gen (F := F) V (r * r)
      (rtLift P F A (show r * t = t * r from mul_comm _ _) ≫ l₁)
      (rtLift P F B (show t * r = t * r from rfl) ≫ l₃)
      (IsFrobeniusType.comp P F (rtLift_frobType P F A _) hl₁F)
      (IsFrobeniusType.comp P F (rtLift_frobType P F B _) hl₃F)
      (by rw [P.degFr_comp, P.degFr_comp,
          show P.degFr l₁ = 1 from isLinear_of_isIso P l₁,
          show P.degFr l₃ = 1 from isLinear_of_isIso P l₃,
          rtLift_degFr, rtLift_degFr])
      hl₁.symm hl₃.symm χ).trans hV
  · refine ⟨(pushIdx (F := F)
      (rtLift P F Cc (show r * r = r * r from mul_comm _ _)) (rtLift_frobType P F Cc _)
      (rtLift P F B (show t * r = r * t from mul_comm _ _)) (rtLift_frobType P F B _)
      (by rw [rtLift_degFr, rtLift_degFr])).obj
      (idxMk (P := P) (F := F) (𝟙 (rtObj P F Cc (r * r))) l₃ hl₂F hl₃F hbe),
      rtMap (F := F) (r * r) γ₁, ?_, ?_⟩
    · rw [rtRootIso_hom_mk]
    · exact prop_1_10_i_pullBack_of P F
        (rtExt_frobType P F Cc (r * r)) (rtExt_frobType P F V.right.obj.2 (r * r))
        ((rtExt_degFr P F Cc (r * r)).trans
          (rtExt_degFr P F V.right.obj.2 (r * r)).symm)
        (rtMap_spec (F := F) (r * r) γ₁) hγ₁

/-! ## ★44. (iv)(b) —— pull-back は LB-invertible かつ linear

★★等長性の手も `𝒞^birat` と同じ: `α = pp ≫ ll`(pre-step ≫ 等長)に分解し、
`ll` を `α` に沿って持ち上げると **`pp` が同型**になる。 -/

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**pull-back は等長**(根が等しい場合)。 -/
theorem pfRoot_isIsometric_of_pullBack_sameRoot {A B : C} {r : ℕ+}
    (α : (⟨A, r⟩ : PfRootObj P F) ⟶ ⟨B, r⟩)
    (hα : IsPullBack (pfRootPre P F) α) : IsIsometric (pfRootPre P F) α := by
  obtain ⟨T, pp, ll, hpp, hll, heq, -⟩ := pfRoot_linearSplit_sameRoot α
    (pfRoot_isLinear_of_pullBack_sameRoot α hα)
  haveI hppb : IsIso ((pfRootPre P F).Base pp) := hpp.2
  have hbase : (pfRootPre P F).Base ll
      = @inv _ _ _ _ ((pfRootPre P F).Base pp) hppb ≫ (pfRootPre P F).Base α := by
    refine eq_inv_comp_of _ hppb _ _ ?_
    rw [← (pfRootPre P F).Base_comp, ← heq]
  obtain ⟨u, ⟨hu1, hu2⟩, -⟩ := IsPullBack.lift (pfRootPre P F) hα T ll _ hbase
  have hkey : pp ≫ u = 𝟙 (⟨A, r⟩ : PfRootObj P F) := by
    refine IsPullBack.hom_ext (pfRootPre P F) hα _ _ ?_ ?_
    · rw [Category.assoc, hu1, ← heq, Category.id_comp]
    · rw [(pfRootPre P F).Base_comp, hu2, IsIso.hom_inv_id, (pfRootPre P F).Base_id]
  haveI : IsIso pp := isIso_of_comp_eq_id (pfRoot_totEpi P F) pp u hkey
  show (pfRootPre P F).Div α = 0
  rw [heq, (pfRootPre P F).Div_comp, show (pfRootPre P F).Div ll = 0 from hll,
    show (pfRootPre P F).Div pp = 0 from isIsometric_of_isIso (pfRootPre P F) pp,
    map_zero, smul_zero, add_zero]

variable {P F} in
/-- ★★★★**[FrdI] Definition 1.3, (iv)(b)** —— `𝒞^pf` 版。 -/
theorem pfRoot_pullBackLB (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (α : X ⟶ Y) (hα : IsPullBack (pfRootPre P F) α) :
    IsLBInvertible (pfRootPre P F) α ∧ IsLinear (pfRootPre P F) α := by
  obtain ⟨eX, hXiso⟩ := pfRoot_exists_iso_root (F := F) X.obj X.root Y.root
    (X.root * Y.root) rfl
  obtain ⟨eY, hYiso⟩ := pfRoot_exists_iso_root (F := F) Y.obj Y.root X.root
    (X.root * Y.root) (mul_comm _ _)
  haveI := hXiso
  haveI := hYiso
  have hpb' : IsPullBack (pfRootPre P F) (inv eX ≫ α ≫ eY) :=
    IsPullBack.comp (pfRootPre P F) (isPullBack_of_isIso (pfRootPre P F) (inv eX))
      (IsPullBack.comp (pfRootPre P F) hα (isPullBack_of_isIso (pfRootPre P F) eY))
  have hlin' := pfRoot_isLinear_of_pullBack_sameRoot (inv eX ≫ α ≫ eY) hpb'
  have hiso' := pfRoot_isIsometric_of_pullBack_sameRoot (inv eX ≫ α ≫ eY) hpb'
  have hlin : IsLinear (pfRootPre P F) α := by
    have h0 : (pfRootPre P F).degFr (inv eX ≫ α ≫ eY) = 1 := hlin'
    rw [(pfRootPre P F).degFr_comp, (pfRootPre P F).degFr_comp,
      show (pfRootPre P F).degFr (inv eX) = 1 from isLinear_of_isIso (pfRootPre P F) (inv eX),
      show (pfRootPre P F).degFr eY = 1 from isLinear_of_isIso (pfRootPre P F) eY,
      one_mul, mul_one] at h0
    exact h0
  refine ⟨⟨pfRoot_isCoAngular hfi α, ?_⟩, hlin⟩
  have h1 : (pfRootPre P F).Div (inv eX ≫ α ≫ eY) = 0 := hiso'
  rw [(pfRootPre P F).Div_comp, (pfRootPre P F).Div_comp,
    show (pfRootPre P F).Div (inv eX) = 0 from
      isIsometric_of_isIso (pfRootPre P F) (inv eX),
    show (pfRootPre P F).Div eY = 0 from isIsometric_of_isIso (pfRootPre P F) eY,
    smul_zero, add_zero, map_zero, zero_add,
    show (pfRootPre P F).degFr eY = 1 from isLinear_of_isIso (pfRootPre P F) eY] at h1
  simp only [PNat.one_coe, one_smul] at h1
  show (pfRootPre P F).Div α = 0
  exact (Φ.pfOn (phiSharp P)).map_injective _ (h1.trans (map_zero _).symm)

/-! ## ★45. pull-back を `Hom^pf` の中で消す

★★残る 3 条(`arbFactor` / `arbFactorUniq` / `plBkEquiv`)はすべて
**「`𝒞` の pull-back は `𝒞^pf` でも pull-back」**に帰する。
★その普遍性のうち**一意性の側**をここで用意する ——
`homPf_cancel_preStep` と同じ形で、`Mono` の代わりに
`IsPullBack.hom_ext`(底の一致も要る)を使う。 -/

variable {P F} in
/-- ★★**pull-back 性は Frobenius 遷移で戻せる**(co-angular が要る)。

★`prop_1_10_i_pullBack_of` の逆向き。★`Proposition 1.4, (ii)` を両向きに使う:
linear と等長は四角形から**両向きに**出るが、co-angular だけは仮定に置く
(`𝒞^pf` の側では始域を isotropic に取れるので、そこはただで手に入る)。 -/
theorem pullBack_of_transport_back (Fc : FrobenioidCore P) {A B A' B' : C}
    {φ : A ⟶ B} {a : A ⟶ A'} {b : B ⟶ B'} {φ' : A' ⟶ B'}
    (ha : IsFrobeniusType P a) (hb : IsFrobeniusType P b)
    (hd : P.degFr a = P.degFr b) (hsq : φ ≫ b = a ≫ φ')
    (hφ' : IsPullBack P φ') (hco : IsCoAngular P φ) : IsPullBack P φ := by
  obtain ⟨hlb', hlin'⟩ := Fc.pullBackLB φ' hφ'
  have hlin : IsLinear P φ := by
    have h0 := congrArg P.degFr hsq
    rw [P.degFr_comp, P.degFr_comp, show P.degFr φ' = 1 from hlin', one_mul, hd] at h0
    show P.degFr φ = 1
    exact mul_left_cancel (h0.trans (mul_one _).symm)
  refine prop_1_4_ii_mpr P Fc φ ⟨hco, ?_⟩ hlin
  have h0 := congrArg P.Div hsq
  rw [P.Div_comp, P.Div_comp, show P.Div b = 0 from hb.1.2,
    show P.Div a = 0 from ha.1.2, show P.Div φ' = 0 from hlb'.2,
    map_zero, map_zero, zero_add, smul_zero, add_zero] at h0
  exact nsmul_eq_zero_of_isSharp (P.divisorial _).2 h0

variable {P F} in
/-- ★遷移は底の一致を保つ。 -/
theorem idxTransport_Base_eq {A B : C} {Z W : IdxPf P F A B} (u : Z ⟶ W)
    (φ φ' : Z.right.obj.1 ⟶ Z.right.obj.2) (h : P.Base φ = P.Base φ') :
    P.Base (idxTransport P F u φ) = P.Base (idxTransport P F u φ') := by
  haveI : IsIso (P.Base u.right.hom.1) := u.right.property.1.2
  refine (cancel_epi (P.Base u.right.hom.1)).mp ?_
  rw [← P.Base_comp, ← P.Base_comp, ← idxTransport_spec, ← idxTransport_spec,
    P.Base_comp, P.Base_comp, h]

variable {P F} in
/-- ★3 脚添字を「**第 2 脚**が isotropic」な所まで押し上げる。 -/
theorem exists_idx3_isotropic2 (hfi : IsOfFrobeniusIsotropicType P) {A B E : C}
    (V : IdxPf3 P F A B E) :
    ∃ (W : IdxPf3 P F A B E) (u : V ⟶ W), IsIsotropic P W.right.obj.2.1 := by
  obtain ⟨Dd, a, ha, hDd⟩ := hfi V.right.obj.2.1
  obtain ⟨A₂, p, hp, hpd⟩ := F.frobDegSurj V.right.obj.1 (P.degFr a)
  obtain ⟨E₂, e, he, hed⟩ := F.frobDegSurj V.right.obj.2.2 (P.degFr a)
  obtain ⟨hva, hvb, hve, hvab, hvbe⟩ := V.hom.property
  refine ⟨Under.mk (Y := (⟨(A₂, Dd, E₂)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(V.hom.hom.1 ≫ p, V.hom.hom.2.1 ≫ a, V.hom.hom.2.2 ≫ e),
          IsFrobeniusType.comp P F hva hp, IsFrobeniusType.comp P F hvb ha,
          IsFrobeniusType.comp P F hve he, ?_, ?_⟩),
    Under.homMk (show V.right ⟶ (⟨(A₂, Dd, E₂)⟩ : TriFr P F) from
      ⟨(p, a, e), hp, ha, he, hpd, hed.symm⟩)
      (WideSubcategory.hom_ext _ rfl), hDd⟩
  · rw [P.degFr_comp, P.degFr_comp, hpd, hvab]
  · rw [P.degFr_comp, P.degFr_comp, hed, hvbe]

set_option maxHeartbeats 1000000 in
variable {P F} in
/-- ★★★**`Hom^pf` の中で pull-back は消せる**(底の一致つき)。 -/
theorem homPf_cancel_pullBack {A B E : C} (V : IdxPf3 P F A B E)
    (φ φ' : V.right.obj.1 ⟶ V.right.obj.2.1) (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2)
    (hψ : IsPullBack P ψ) (hb : P.Base φ = P.Base φ')
    (h : HomPf.mk ((idx13 P F A B E).obj V) (φ ≫ ψ)
      = HomPf.mk ((idx13 P F A B E).obj V) (φ' ≫ ψ)) :
    HomPf.mk ((idx12 P F A B E).obj V) φ = HomPf.mk ((idx12 P F A B E).obj V) φ' := by
  obtain ⟨V', t, t', ht⟩ := HomPf.eq_iff.mp h
  rw [idx_hom_ext t' t] at ht
  obtain ⟨V'', ⟨k⟩⟩ := exists_hom_of_final (idx13 P F A B E) V'
  set s : V ⟶ IsFiltered.max V V'' := IsFiltered.leftToMax V V'' with hs
  set r : V'' ⟶ IsFiltered.max V V'' := IsFiltered.rightToMax V V'' with hr
  have hm : t ≫ k ≫ (idx13 P F A B E).map r = (idx13 P F A B E).map s :=
    idx_hom_ext _ _
  have hA : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have hB : idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map r)
          (idxTransport P F k (idxTransport P F t (φ' ≫ ψ))) := by
    rw [← hm, idxTransport_comp, idxTransport_comp]
  have key : idxTransport P F ((idx13 P F A B E).map s) (φ ≫ ψ)
      = idxTransport P F ((idx13 P F A B E).map s) (φ' ≫ ψ) := by
    rw [hA, hB, ht]
  have e1 := idxTransport_comp_pair (F := F) s φ ψ
  have e2 := idxTransport_comp_pair (F := F) s φ' ψ
  have hcan : idxTransport P F ((idx12 P F A B E).map s) φ
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ
      = idxTransport P F ((idx12 P F A B E).map s) φ'
      ≫ idxTransport P F ((idx23 P F A B E).map s) ψ := by
    rw [e1, e2]
    exact key
  have hψ' : IsPullBack P (idxTransport P F ((idx23 P F A B E).map s) ψ) :=
    prop_1_10_i_pullBack_of P F
      ((idx23 P F A B E).map s).right.property.1
      ((idx23 P F A B E).map s).right.property.2.1
      ((idx23 P F A B E).map s).right.property.2.2
      (idxTransport_spec ((idx23 P F A B E).map s) ψ) hψ
  have hb' : P.Base (idxTransport P F ((idx12 P F A B E).map s) φ)
      = P.Base (idxTransport P F ((idx12 P F A B E).map s) φ') :=
    idxTransport_Base_eq ((idx12 P F A B E).map s) φ φ' hb
  have hφφ' : idxTransport P F ((idx12 P F A B E).map s) φ
      = idxTransport P F ((idx12 P F A B E).map s) φ' :=
    IsPullBack.hom_ext P hψ' _ _ hcan hb'
  rw [← HomPf.mk_map ((idx12 P F A B E).map s) φ, ← HomPf.mk_map ((idx12 P F A B E).map s) φ',
    hφφ']

/-! ## ★46. pull-back の普遍性の**持ち上げの側**への下ごしらえ -/

variable {P F} in
/-- ★★**代表元での底の特徴づけ**(`inv` を使わない形)。 -/
theorem rootBase_mk_spec {X Y : PfRootObj P F}
    (W : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (ψ : W.right.obj.1 ⟶ W.right.obj.2) :
    rootBase (show HomRoot P F X Y from HomPf.mk W ψ)
        ≫ P.Base (rtExt P F Y.obj X.root) ≫ P.Base W.hom.hom.2
      = (P.Base (rtExt P F X.obj Y.root) ≫ P.Base W.hom.hom.1) ≫ P.Base ψ := by
  rw [← Category.assoc, rootBase_spec, Category.assoc, pfBase_mk, repBase_spec,
    Category.assoc]
  rfl

variable {P F} in
/-- ★★**代表元が pull-back であることは代表元の取り方によらない**
(co-angular が言えている限り)。 -/
theorem pfRoot_rep_isPullBack {A B : C} {Z₀ Z : IdxPf P F A B}
    {χ₀ : Z₀.right.obj.1 ⟶ Z₀.right.obj.2} {χ : Z.right.obj.1 ⟶ Z.right.obj.2}
    (hχ₀ : IsPullBack P χ₀) (heq : HomPf.mk Z₀ χ₀ = HomPf.mk Z χ)
    (hco : IsCoAngular P χ) : IsPullBack P χ := by
  obtain ⟨V, u, v, huv⟩ := HomPf.eq_iff.mp heq
  have h1 : IsPullBack P (idxTransport P F u χ₀) :=
    prop_1_10_i_pullBack_of P F u.right.property.1 u.right.property.2.1
      u.right.property.2.2 (idxTransport_spec u χ₀) hχ₀
  rw [huv] at h1
  exact pullBack_of_transport_back F v.right.property.1 v.right.property.2.1
    v.right.property.2.2 (idxTransport_spec v χ) h1 hco

variable {P F} in
/-- ★**終域を共有する 2 射を同じ 3 脚添字へ**(`exists_rep3` の co-span 版)。 -/
theorem exists_rep3_cospan {A B E : C} (f : HomPf P F A E) (w : HomPf P F B E) :
    ∃ (V : IdxPf3 P F A B E) (θ : V.right.obj.1 ⟶ V.right.obj.2.2)
      (ψ : V.right.obj.2.1 ⟶ V.right.obj.2.2),
      f = HomPf.mk ((idx13 P F A B E).obj V) θ ∧
      w = HomPf.mk ((idx23 P F A B E).obj V) ψ := by
  obtain ⟨Zf, θ₀, rfl⟩ := HomPf.exists_rep f
  obtain ⟨Zw, ψ₀, rfl⟩ := HomPf.exists_rep w
  obtain ⟨hfa, hfe, hfae⟩ := Zf.hom.property
  obtain ⟨hwb, hwe, hwbe⟩ := Zw.hom.property
  obtain ⟨B₁, b₁, hb₁, hb₁d⟩ := F.frobDegSurj B (P.degFr Zf.hom.hom.1)
  obtain ⟨A₂, a₂, ha₂, ha₂d⟩ := F.frobDegSurj A (P.degFr Zw.hom.hom.1)
  set Vf : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(Zf.right.obj.1, B₁, Zf.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(Zf.hom.hom.1, b₁, Zf.hom.hom.2), hfa, hb₁, hfe, hb₁d.symm,
          hb₁d.trans hfae⟩) with hVf
  set Vw : IdxPf3 P F A B E :=
    Under.mk (Y := (⟨(A₂, Zw.right.obj.1, Zw.right.obj.2)⟩ : TriFr P F))
      (show triFrObj P F A B E ⟶ _ from
        ⟨(a₂, Zw.hom.hom.1, Zw.hom.hom.2), ha₂, hwb, hwe, ha₂d, hwbe⟩) with hVw
  exact ⟨IsFiltered.max Vf Vw,
    idxTransport P F ((idx13 P F A B E).map (IsFiltered.leftToMax Vf Vw)) θ₀,
    idxTransport P F ((idx23 P F A B E).map (IsFiltered.rightToMax Vf Vw)) ψ₀,
    (HomPf.mk_map _ θ₀).symm, (HomPf.mk_map _ ψ₀).symm⟩

end ABC3.Found.FrdI
