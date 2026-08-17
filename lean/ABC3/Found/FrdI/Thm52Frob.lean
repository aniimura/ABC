/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52
import ABC3.Found.FrdI.Prop15
import ABC3.Found.FrdI.Prop17
import ABC3.Found.FrdI.Remark311
import ABC3.Found.FrdI.Def31

/-!
# [FrdI] Theorem 5.2, (ii) —— model Frobenioid は Frobenioid である

原文 (FrdI p.101):
> (ii) The category C is a Frobenioid [with respect to the functor C

★原文は「(i)(ii) は routine な検証で従う。(ii) については
`Proposition 1.5, (i)` の『初等 Frobenioid は Frobenioid である』の検証に似ている」
とだけ言う。ここではその検証を実際に行う。

★★**`𝔽_Φ` との違いは 1 点**である —— 射に `u ∈ B(A_𝒟)` が乗っており、
関係式 `deg·α + Div = Base^*(β) + Div_B(u)` で `α` と結ばれている。
`B` は group-like なので **`u` は必ず可逆**で、たいていの条件では
「`u` を後から合わせる」ことができる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

namespace ModelData

/-! ## ★1. 同型射の判定

★★**`𝔽_Φ` と同じ**(`ElemFrobCat.isIso_iff`)—— 底が同型・零因子が 0・次数 1。
`u` 成分は条件に出てこない: `B` が group-like なので**逆元がいつでも取れる**からである。 -/

/-- ★`B` の元はつねに可逆(`B` は group-like)。 -/
theorem bmon_isAddUnit (h : Hyp M) (A : D) (x : M.bmon.val A) : IsAddUnit x :=
  (isGroupLike_iff _).mp (h.bmonGroupLike A) x

/-- ★★**同型射の判定**。 -/
theorem model_isIso_iff (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    IsIso φ ↔ IsIso φ.base ∧ φ.div = 0 ∧ φ.deg = 1 := by
  constructor
  · intro hiso
    have hd1 := congrArg Hom.deg (IsIso.hom_inv_id φ)
    rw [comp_deg, id_deg] at hd1
    have hd1' : ((Hom.deg (inv φ) : ℕ+) : ℕ) * ((Hom.deg φ : ℕ+) : ℕ) = 1 := by
      exact_mod_cast congrArg (fun n : ℕ+ => (n : ℕ)) hd1
    have hdegi : Hom.deg (inv φ) = 1 :=
      PNat.coe_eq_one_iff.mp (Nat.dvd_one.mp ⟨(Hom.deg φ : ℕ+), hd1'.symm⟩)
    have hdeg : Hom.deg φ = 1 :=
      PNat.coe_eq_one_iff.mp
        (Nat.dvd_one.mp ⟨(Hom.deg (inv φ) : ℕ+), by rw [mul_comm]; exact hd1'.symm⟩)
    have hdv := congrArg Hom.div (IsIso.hom_inv_id φ)
    rw [comp_div, id_div, hdegi] at hdv
    simp only [PNat.one_coe, one_smul] at hdv
    have hunit : IsAddUnit φ.div :=
      ⟨⟨φ.div, M.phi.map φ.base (Hom.div (inv φ)), by rw [add_comm]; exact hdv, hdv⟩, rfl⟩
    refine ⟨⟨Hom.base (inv φ), ?_, ?_⟩, (h.divisorial A.base).2 _ hunit, hdeg⟩
    · rw [← comp_base, IsIso.hom_inv_id, id_base]
    · rw [← comp_base, IsIso.inv_hom_id, id_base]
  · rintro ⟨hbase, hdiv, hdeg⟩
    -- ★`u` の逆元を取る
    obtain ⟨U, hU⟩ := bmon_isAddUnit h B.base (M.bmon.map (inv φ.base) φ.u)
    -- ★`Base(φ)^*` と `Base(φ)⁻¹^*` は互いに逆
    have hmapAB : ∀ x : M.bmon.val A.base,
        M.bmon.map φ.base (M.bmon.map (inv φ.base) x) = x := by
      intro x
      rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
    have hmapBA : ∀ x : Gp (M.phi.val B.base),
        M.phi.gpMapOn (inv φ.base) (M.phi.gpMapOn φ.base x) = x := by
      intro x
      rw [← MonoidOn.gpMapOn_comp, IsIso.inv_hom_id, MonoidOn.gpMapOn_id]
    -- ★逆射の関係式
    have hcond : ((1 : ℕ+) : ℕ) • B.cls + toGpHom _ (0 : M.phi.val B.base)
        = M.phi.gpMapOn (inv φ.base) A.cls + M.divB _ U.neg := by
      have hφ := φ.cond
      rw [hdiv, hdeg] at hφ
      have h2 := congrArg (M.phi.gpMapOn (inv φ.base)) hφ
      rw [map_add, map_add, map_nsmul, hmapBA] at h2
      have h3 : M.phi.gpMapOn (inv φ.base) (M.divB A.base φ.u)
          = M.divB B.base (M.bmon.map (inv φ.base) φ.u) := (M.divB_nat _ _).symm
      rw [h3, ← hU] at h2
      simp only [PNat.one_coe, one_smul, map_zero, add_zero] at h2 ⊢
      rw [h2, add_assoc, ← map_add, U.val_neg, map_zero, add_zero]
    refine ⟨⟨⟨inv φ.base, 0, 1, U.neg, hcond⟩, ?_, ?_⟩⟩
    · refine Hom.ext (IsIso.hom_inv_id _) ?_ ?_ ?_
      · show M.phi.map φ.base (0 : M.phi.val B.base) + ((1 : ℕ+) : ℕ) • φ.div
          = (0 : M.phi.val A.base)
        simp [hdiv]
      · show (1 : ℕ+) * φ.deg = 1
        simp [hdeg]
      · show M.bmon.map φ.base U.neg + ((1 : ℕ+) : ℕ) • φ.u = (0 : M.bmon.val A.base)
        have h4 : M.bmon.map φ.base (U : M.bmon.val B.base) = φ.u := by
          rw [hU, hmapAB]
        have h5 : M.bmon.map φ.base (U : M.bmon.val B.base)
            + M.bmon.map φ.base U.neg = 0 := by
          rw [← map_add, U.val_neg, map_zero]
        simp only [PNat.one_coe, one_smul]
        rw [add_comm, ← h4, h5]
    · refine Hom.ext (IsIso.inv_hom_id _) ?_ ?_ ?_
      · show M.phi.map (inv φ.base) φ.div + (φ.deg : ℕ) • (0 : M.phi.val B.base)
          = (0 : M.phi.val B.base)
        simp [hdiv]
      · show φ.deg * 1 = 1
        simp [hdeg]
      · show M.bmon.map (inv φ.base) φ.u + (φ.deg : ℕ) • U.neg = (0 : M.bmon.val B.base)
        rw [hdeg, ← hU]
        simp only [PNat.one_coe, one_smul]
        exact U.val_neg

/-! ## ★2. isotropic 型・すべての射は co-angular -/

/-- ★★**すべての対象は isotropic**。 -/
theorem model_isotropic (h : Hyp M) (A : Obj M) : IsIsotropic (modelPre h) A := by
  intro Dd φ hiso hstep
  exact (model_isIso_iff h φ).mpr ⟨hstep.2, hiso, hstep.1⟩

/-- ★★`𝒞` は **isotropic 型**。 -/
theorem model_isotropicType (h : Hyp M) : IsOfIsotropicType (modelPre h) :=
  model_isotropic h

/-- ★★**すべての射は co-angular**(`Proposition 1.4, (i)`)。 -/
theorem model_coAngular (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    IsCoAngular (modelPre h) φ :=
  prop_1_4_i _ φ (fun X _ => model_isotropic h X)

/-- ★**Frobenius 型 ⟺ isometric な base-isomorphism**。 -/
theorem model_frobeniusType_iff (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    IsFrobeniusType (modelPre h) φ ↔ (φ.div = 0 ∧ IsIso φ.base) :=
  ⟨fun hf => ⟨hf.1.2, hf.2⟩, fun hf => ⟨⟨model_coAngular h φ, hf.1⟩, hf.2⟩⟩

/-! ## ★3. `B` の逆元を道具にする

★`B` は group-like なので加法逆元がつねに取れる。**一意**なので
`B` の誘導射とも `Div_B` とも交換する。 -/

theorem exists_bneg (h : Hyp M) {A : D} (x : M.bmon.val A) : ∃ y, x + y = 0 := by
  obtain ⟨U, hU⟩ := bmon_isAddUnit h A x
  exact ⟨U.neg, by rw [← hU]; exact U.val_neg⟩

/-- ★`B(A)` の加法逆元。 -/
noncomputable def bneg (h : Hyp M) {A : D} (x : M.bmon.val A) : M.bmon.val A :=
  Classical.choose (exists_bneg h x)

@[simp] theorem add_bneg (h : Hyp M) {A : D} (x : M.bmon.val A) : x + bneg h x = 0 :=
  Classical.choose_spec (exists_bneg h x)

@[simp] theorem bneg_add (h : Hyp M) {A : D} (x : M.bmon.val A) : bneg h x + x = 0 := by
  rw [add_comm]; exact add_bneg h x

/-- ★可換モノイドでは逆元は一意。 -/
theorem eq_bneg (h : Hyp M) {A : D} {x y : M.bmon.val A} (hxy : x + y = 0) :
    y = bneg h x := by
  calc y = y + (x + bneg h x) := by rw [add_bneg, add_zero]
    _ = (x + y) + bneg h x := by abel
    _ = bneg h x := by rw [hxy, zero_add]

theorem map_bneg (h : Hyp M) {A B : D} (f : B ⟶ A) (x : M.bmon.val A) :
    M.bmon.map f (bneg h x) = bneg h (M.bmon.map f x) :=
  eq_bneg h (by rw [← map_add, add_bneg, map_zero])

theorem divB_bneg (h : Hyp M) {A : D} (x : M.bmon.val A) :
    M.divB A (bneg h x) = -(M.divB A x) :=
  eq_neg_of_add_eq_zero_right (by rw [← map_add, add_bneg, map_zero])

/-! ## ★4. pull-back 射の判定

原文 (FrdI p.27):
> pull-back morphism if and only if it is a linear isometry. The fact that FΦ satisfies

★★**`𝔽_Φ` と同じ**「線型な等長射」である。`u` 成分は
**`B` が group-like なので後から合わせられる**。

★★順向きの試験点は `𝔽_Φ` と同じには取れない —— `⟨Base φ, 0, 1⟩` に当たる射を
`A` から出すことは(`α` の縛りがあるので)一般にできない。
代わりに **`X := (A_𝒟, Base(φ)^*(β))`** を取ると `⟨Base φ, 0, 1, 0⟩` が作れる。 -/

theorem model_isPullBack_iff (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    IsPullBack (modelPre h) φ ↔ (φ.deg = 1 ∧ φ.div = 0) := by
  constructor
  · intro hpb
    have hcond0 : ((1 : ℕ+) : ℕ) • (M.phi.gpMapOn φ.base B.cls)
        + toGpHom _ (0 : M.phi.val A.base)
        = M.phi.gpMapOn φ.base B.cls + M.divB _ (0 : M.bmon.val A.base) := by simp
    obtain ⟨g0, hg0b, hg0d, hg0deg⟩ :
        ∃ g0 : (⟨A.base, M.phi.gpMapOn φ.base B.cls⟩ : Obj M) ⟶ B,
          g0.base = φ.base ∧ g0.div = 0 ∧ g0.deg = 1 :=
      ⟨⟨φ.base, 0, 1, 0, hcond0⟩, rfl, rfl, rfl⟩
    obtain ⟨f, hf⟩ := (hpb _).2
      ⟨(g0, 𝟙 A.base), by
        show g0.base = 𝟙 A.base ≫ φ.base
        rw [hg0b, Category.id_comp]⟩
    have hp := Subtype.ext_iff.mp hf
    have h1 : (f ≫ φ : (⟨A.base, M.phi.gpMapOn φ.base B.cls⟩ : Obj M) ⟶ B) = g0 :=
      congrArg Prod.fst hp
    have h2 : f.base = 𝟙 A.base := congrArg Prod.snd hp
    have hdeg : φ.deg = 1 := by
      have hd := congrArg Hom.deg h1
      rw [comp_deg, hg0deg] at hd
      exact pnat_left_eq_one hd
    refine ⟨hdeg, ?_⟩
    have hd := congrArg Hom.div h1
    rw [comp_div, h2, hdeg, hg0d] at hd
    simp only [PNat.one_coe, one_smul, MonoidOn.map_id] at hd
    exact (h.divisorial A.base).2 _
      ⟨⟨φ.div, f.div, hd, by rw [add_comm]; exact hd⟩, rfl⟩
  · rintro ⟨hdeg, hdiv⟩
    intro X
    constructor
    · intro f₁ f₂ hf
      have hp := Subtype.ext_iff.mp hf
      have hb : f₁.base = f₂.base := congrArg Prod.snd hp
      have hc : (f₁ ≫ φ : X ⟶ B) = f₂ ≫ φ := congrArg Prod.fst hp
      letI := isCancelAdd_of_isIntegralMonoid _ (h.divisorial X.base).1.1
      have hdg : f₁.deg = f₂.deg := by
        have hd := congrArg Hom.deg hc
        rw [comp_deg, comp_deg] at hd
        exact mul_left_cancel hd
      refine Hom.ext hb ?_ hdg ?_
      · have hd := congrArg Hom.div hc
        rw [comp_div, comp_div, hb, hdeg] at hd
        simp only [PNat.one_coe, one_smul] at hd
        exact add_left_cancel hd
      · have hd := congrArg Hom.u hc
        rw [comp_u, comp_u, hb, hdeg] at hd
        simp only [PNat.one_coe, one_smul] at hd
        have hd2 : f₁.u + M.bmon.map f₂.base φ.u = f₂.u + M.bmon.map f₂.base φ.u := by
          rw [add_comm f₁.u, add_comm f₂.u]; exact hd
        exact groupLike_add_right_cancel (h.bmonGroupLike X.base) hd2
    · rintro ⟨⟨g, b0⟩, hgb⟩
      simp only [modelPre_Base] at hgb
      -- ★`b0` の型は `toElem` 経由で書かれているので、素の型に取り直す
      obtain ⟨b, rfl⟩ : ∃ b : X.base ⟶ A.base, b = b0 := ⟨b0, rfl⟩
      have hgb' : g.base = b ≫ φ.base := hgb
      have hφ := φ.cond
      rw [hdeg, hdiv] at hφ
      simp only [PNat.one_coe, one_smul, map_zero, add_zero] at hφ
      have hb2 : M.phi.gpMapOn b A.cls
          = M.phi.gpMapOn g.base B.cls + M.divB X.base (M.bmon.map b φ.u) := by
        rw [hφ, (M.phi.gpMapOn b).map_add, ← MonoidOn.gpMapOn_comp, ← hgb',
          M.divB_nat b φ.u]
      refine ⟨⟨b, g.div, g.deg, bneg h (M.bmon.map b φ.u) + g.u, ?_⟩, ?_⟩
      · rw [g.cond, hb2, map_add, divB_bneg]
        abel
      · refine Subtype.ext (Prod.ext ?_ rfl)
        refine Hom.ext hgb'.symm ?_ ?_ ?_
        · show M.phi.map b φ.div + (φ.deg : ℕ) • g.div = g.div
          rw [hdiv, hdeg]
          simp
        · show φ.deg * g.deg = g.deg
          rw [hdeg, one_mul]
        · show M.bmon.map b φ.u + (φ.deg : ℕ) • (bneg h (M.bmon.map b φ.u) + g.u) = g.u
          rw [hdeg]
          simp only [PNat.one_coe, one_smul]
          rw [← add_assoc, add_bneg, zero_add]

/-- ★`𝔽_Φ` と同じく、pre-Frobenioid の段階で
「pull-back ⟺ LB-invertible かつ linear」が出る。 -/
theorem model_isPullBack_iff_lbInvertible (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    IsPullBack (modelPre h) φ ↔ (IsLBInvertible (modelPre h) φ ∧ IsLinear (modelPre h) φ) := by
  rw [model_isPullBack_iff h]
  exact ⟨fun hh => ⟨⟨model_coAngular h φ, hh.2⟩, hh.1⟩, fun hh => ⟨hh.2, hh.1.2⟩⟩

/-! ## ★5. `Definition 1.3` の 21 条 —— 易しい側

★`𝔽_Φ`(`Prop15.lean`)の証明をそのまま写し、`u` 成分を並走させる。
★★**対象を作るときに `α` を選べる**のが model Frobenioid の特徴で、
`𝔽_Φ` で「同じ対象」を使い回していた所は「`α` をずらした対象」に置き換わる。 -/

/-- ★`Gp` の元は `toGp a − toGp b` と書ける。 -/
theorem gp_sub_repr {N : Type w} [AddCommMonoid N] (x : Gp N) :
    ∃ a b : N, x = toGpHom N a - toGpHom N b := by
  induction x using AddLocalization.induction_on with | _ y =>
  refine ⟨y.1, (y.2 : N), eq_sub_of_add_eq ?_⟩
  exact mk_add_toGp N y.1 y.2

/-- ★`(d, 0)` の次数 `n` の base-identity Frobenius 自己射。 -/
def frobZeroEnd (d : D) (n : ℕ+) : End (⟨d, 0⟩ : Obj M) :=
  ⟨𝟙 d, 0, n, 0, by simp⟩

@[simp] theorem frobZeroEnd_base (d : D) (n : ℕ+) :
    (frobZeroEnd (M := M) d n).base = 𝟙 d := rfl
@[simp] theorem frobZeroEnd_div (d : D) (n : ℕ+) :
    (frobZeroEnd (M := M) d n).div = 0 := rfl
@[simp] theorem frobZeroEnd_deg (d : D) (n : ℕ+) :
    (frobZeroEnd (M := M) d n).deg = n := rfl
@[simp] theorem frobZeroEnd_u (d : D) (n : ℕ+) :
    (frobZeroEnd (M := M) d n).u = 0 := rfl

/-- ★★**`α = 0` の対象は Frobenius-trivial**。 -/
theorem model_frobTrivial_zero (h : Hyp M) (d : D) :
    IsFrobeniusTrivial (modelPre h) (⟨d, 0⟩ : Obj M) := by
  refine ⟨⟨⟨frobZeroEnd (M := M) d, ?_⟩, ?_⟩,
    fun n => rfl, fun n => ⟨rfl, ⟨⟨model_coAngular h _, rfl⟩, ?_⟩⟩⟩
  · refine Hom.ext ?_ ?_ ?_ ?_ <;> simp
  · intro m n
    have hmul : (frobZeroEnd (M := M) d m) * (frobZeroEnd (M := M) d n)
        = (frobZeroEnd (M := M) d n) ≫ (frobZeroEnd (M := M) d m) := rfl
    rw [hmul]
    refine Hom.ext ?_ ?_ ?_ ?_ <;> simp
  · show IsIso ((modelPre h).Base (frobZeroEnd (M := M) d _))
    show IsIso (𝟙 d)
    infer_instance

/-- **(i)(a)** —— `(Y, 0)` を取ればよい。 -/
theorem model_baseSurj (h : Hyp M) (Y : D) :
    ∃ A : Obj M, IsFrobeniusTrivial (modelPre h) A ∧
      Nonempty (((modelPre h).toElem.obj A).base ≅ Y) :=
  ⟨(⟨Y, 0⟩ : Obj M), model_frobTrivial_zero h Y, ⟨Iso.refl Y⟩⟩

/-- **(ii)** —— `(A_𝒟, n·α)` へ次数 `n` の Frobenius 型射が入る。 -/
theorem model_frobDegSurj (h : Hyp M) (A : Obj M) (n : ℕ+) :
    ∃ (B : Obj M) (φ : A ⟶ B), IsFrobeniusType (modelPre h) φ ∧
      (modelPre h).degFr φ = n := by
  refine ⟨(⟨A.base, (n : ℕ) • A.cls⟩ : Obj M),
    (⟨𝟙 A.base, 0, n, 0, by simp⟩ : A ⟶ ⟨A.base, (n : ℕ) • A.cls⟩),
    ⟨⟨model_coAngular h _, rfl⟩, ?_⟩, rfl⟩
  show IsIso (𝟙 A.base)
  infer_instance

/-- **(iii)(a)** —— すべての射が co-angular なので自明。 -/
theorem model_coAngularComp (h : Hyp M) {A B E : Obj M} (ψ : A ⟶ B) (φ : B ⟶ E) :
    IsCoAngular (modelPre h) ψ → IsCoAngular (modelPre h) φ →
      IsCoAngular (modelPre h) (ψ ≫ φ) :=
  fun _ _ => model_coAngular h _

/-- **(iii)(b)** —— 同上。 -/
theorem model_coAngularOfPreStep (h : Hyp M) {A' A : Obj M} (α : A' ⟶ A) :
    IsCoAngular (modelPre h) α → IsPreStep (modelPre h) α →
      ∀ φ : A' ⟶ A, IsCoAngular (modelPre h) φ :=
  fun _ _ φ => model_coAngular h φ

/-- **(iv)(b)** -/
theorem model_pullBackLB (h : Hyp M) {A B : Obj M} (α : A ⟶ B)
    (hα : IsPullBack (modelPre h) α) :
    IsLBInvertible (modelPre h) α ∧ IsLinear (modelPre h) α :=
  (model_isPullBack_iff_lbInvertible h α).mp hα

/-- **(iv)(a)** —— `Frobenius 型 → pre-step → pull-back` の 3 分解。

★`𝔽_Φ` と同じく **3 成分がそのまま 3 因子に分かれる**が、
対象の側では `α` が `deg·α`、`deg·α + Div` と 2 度ずれる。 -/
theorem model_arbFactor (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    ∃ (X Y : Obj M) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (modelPre h) γ ∧ IsPreStep (modelPre h) β ∧
        IsPullBack (modelPre h) α := by
  have hα : ((1 : ℕ+) : ℕ) • ((φ.deg : ℕ) • A.cls + toGpHom _ φ.div)
      + toGpHom _ (0 : M.phi.val A.base)
      = M.phi.gpMapOn φ.base B.cls + M.divB _ φ.u := by
    simpa using φ.cond
  refine ⟨(⟨A.base, (φ.deg : ℕ) • A.cls⟩ : Obj M),
    (⟨A.base, (φ.deg : ℕ) • A.cls + toGpHom _ φ.div⟩ : Obj M),
    (⟨𝟙 A.base, 0, φ.deg, 0, by simp⟩ : A ⟶ ⟨A.base, (φ.deg : ℕ) • A.cls⟩),
    (⟨𝟙 A.base, φ.div, 1, 0, by simp⟩ :
      (⟨A.base, (φ.deg : ℕ) • A.cls⟩ : Obj M) ⟶
        ⟨A.base, (φ.deg : ℕ) • A.cls + toGpHom _ φ.div⟩),
    (⟨φ.base, 0, 1, φ.u, hα⟩ :
      (⟨A.base, (φ.deg : ℕ) • A.cls + toGpHom _ φ.div⟩ : Obj M) ⟶ B), ?_,
    ⟨⟨model_coAngular h _, rfl⟩, ?_⟩, ⟨rfl, ?_⟩,
    (model_isPullBack_iff h _).mpr ⟨rfl, rfl⟩⟩
  · refine Hom.ext ?_ ?_ ?_ ?_ <;> simp
  · show IsIso (𝟙 A.base)
    infer_instance
  · show IsIso (𝟙 A.base)
    infer_instance

/-- **(v)(b)** —— `φ = φ ≫ 𝟙`。 -/
theorem model_preStepFactor (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    (hφ : IsPreStep (modelPre h) φ) :
    ∃ (X : Obj M) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular (modelPre h) β ∧ IsPreStep (modelPre h) β ∧
        IsIsometric (modelPre h) α ∧ IsPreStep (modelPre h) α :=
  ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, model_coAngular h φ, hφ,
    (modelPre h).Div_id B, isPreStep_id _ B⟩

/-- **(v)(c)** —— `φ = 𝟙 ≫ φ`。 -/
theorem model_preStepFactor' (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    (hφ : IsPreStep (modelPre h) φ) :
    ∃ (X : Obj M) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric (modelPre h) β ∧ IsPreStep (modelPre h) β ∧
        IsCoAngular (modelPre h) α ∧ IsPreStep (modelPre h) α :=
  ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm, (modelPre h).Div_id A,
    isPreStep_id _ A, model_coAngular h φ, hφ⟩

/-- **(vii)(a)** —— すべての対象が isotropic なので `𝟙` が isotropic hull。 -/
theorem model_isotropicHullExists (h : Hyp M) (A : Obj M) :
    ∃ (B : Obj M) (φ : A ⟶ B), IsIsotropicHull (modelPre h) φ :=
  ⟨A, 𝟙 A, (modelPre h).Div_id A, isPreStep_id _ A, model_isotropic h A,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

/-- **(vii)(b)** -/
theorem model_isotropicClosed (h : Hyp M) {A B : Obj M} (_φ : A ⟶ B) :
    IsIsotropic (modelPre h) A → IsIsotropic (modelPre h) B :=
  fun _ => model_isotropic h B

/-- **(v)(a)** —— pre-step は monomorphism。 -/
theorem model_preStepMono (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    (hφ : IsPreStep (modelPre h) φ) : Mono φ := by
  haveI : IsIso φ.base := hφ.2
  have hdeg : φ.deg = 1 := hφ.1
  refine ⟨fun {Z} f g hfg => ?_⟩
  letI := isCancelAdd_of_isIntegralMonoid _ (h.divisorial Z.base).1.1
  have hb : f.base = g.base := by
    have hh := congrArg Hom.base hfg
    rw [comp_base, comp_base] at hh
    exact (cancel_mono φ.base).mp hh
  refine Hom.ext hb ?_ ?_ ?_
  · have hh := congrArg Hom.div hfg
    rw [comp_div, comp_div, hb, hdeg] at hh
    simp only [PNat.one_coe, one_smul] at hh
    exact add_left_cancel hh
  · have hh := congrArg Hom.deg hfg
    rw [comp_deg, comp_deg] at hh
    exact mul_left_cancel hh
  · have hh := congrArg Hom.u hfg
    rw [comp_u, comp_u, hb, hdeg] at hh
    simp only [PNat.one_coe, one_smul] at hh
    refine groupLike_add_right_cancel (h.bmonGroupLike Z.base) (a := f.u) (b := g.u)
      (c := M.bmon.map g.base φ.u) ?_
    rw [add_comm f.u, add_comm g.u]
    exact hh

/-- **(ii)** の本質的一意性。 -/
theorem model_frobDegUniq (h : Hyp M) (A B E : Obj M) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (modelPre h) φ) (hψ : IsFrobeniusType (modelPre h) ψ)
    (hd : (modelPre h).degFr φ = (modelPre h).degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI hbφ : IsIso φ.base := hφ.2
  haveI hbψ : IsIso ψ.base := hψ.2
  have hφd : φ.div = 0 := hφ.1.2
  have hψd : ψ.div = 0 := hψ.1.2
  have hdeg : φ.deg = ψ.deg := hd
  -- ★`Base(φ)^*` と `Base(φ)⁻¹^*` は互いに逆
  have hmapAB : ∀ x : M.bmon.val A.base,
      M.bmon.map φ.base (M.bmon.map (inv φ.base) x) = x := by
    intro x
    rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]
  have hgpBB : ∀ x : Gp (M.phi.val B.base),
      M.phi.gpMapOn (inv φ.base) (M.phi.gpMapOn φ.base x) = x := by
    intro x
    rw [← MonoidOn.gpMapOn_comp, IsIso.inv_hom_id, MonoidOn.gpMapOn_id]
  have hφc := φ.cond
  have hψc := ψ.cond
  rw [hφd] at hφc
  rw [hψd, ← hdeg] at hψc
  simp only [map_zero, add_zero] at hφc hψc
  have hkey : M.phi.gpMapOn φ.base B.cls + M.divB _ φ.u
      = M.phi.gpMapOn ψ.base E.cls + M.divB _ ψ.u := by rw [← hφc, ← hψc]
  have hcond : ((1 : ℕ+) : ℕ) • B.cls + toGpHom _ (0 : M.phi.val B.base)
      = M.phi.gpMapOn (inv φ.base ≫ ψ.base) E.cls
        + M.divB _ (M.bmon.map (inv φ.base) ψ.u
            + bneg h (M.bmon.map (inv φ.base) φ.u)) := by
    have h2 := congrArg (M.phi.gpMapOn (inv φ.base)) hkey
    rw [(M.phi.gpMapOn (inv φ.base)).map_add, (M.phi.gpMapOn (inv φ.base)).map_add,
      hgpBB, ← MonoidOn.gpMapOn_comp, ← M.divB_nat (inv φ.base) φ.u,
      ← M.divB_nat (inv φ.base) ψ.u] at h2
    rw [map_add, divB_bneg]
    simp only [PNat.one_coe, one_smul, map_zero, add_zero]
    rw [eq_sub_of_add_eq h2]
    abel
  refine ⟨(⟨inv φ.base ≫ ψ.base, 0, 1,
    M.bmon.map (inv φ.base) ψ.u + bneg h (M.bmon.map (inv φ.base) φ.u), hcond⟩ : B ⟶ E),
    (model_isIso_iff h _).mpr ⟨inferInstance, rfl, rfl⟩, ?_⟩
  refine Hom.ext ?_ ?_ ?_ ?_
  · show φ.base ≫ (inv φ.base ≫ ψ.base) = ψ.base
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  · show M.phi.map φ.base (0 : M.phi.val B.base) + ((1 : ℕ+) : ℕ) • φ.div = ψ.div
    rw [hφd, hψd]
    simp
  · show (1 : ℕ+) * φ.deg = ψ.deg
    rw [one_mul, hdeg]
  · show M.bmon.map φ.base (M.bmon.map (inv φ.base) ψ.u
        + bneg h (M.bmon.map (inv φ.base) φ.u)) + ((1 : ℕ+) : ℕ) • φ.u = ψ.u
    rw [map_add, hmapAB, map_bneg, hmapAB]
    simp only [PNat.one_coe, one_smul]
    rw [add_assoc, bneg_add, add_zero]

/-! ## ★6. `𝒪^▷(A)` の記述

★★**`𝔽_Φ` との違いがはっきり出るところ** —— `𝔽_Φ` では `𝒪^▷(A) ≅ Φ(A)` だったが、
model Frobenioid では

  `𝒪^▷(A) ≅ {(a, u) ∈ Φ(A) × B(A) : toGp a = Div_B(u)}`

である(`α` を動かさないための縛りが入る)。 -/

/-- ★`(a, u)` から作る `𝒪^▷(A)` の元。 -/
def otriOf (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : End A :=
  ⟨𝟙 A.base, a, 1, u, by simp [hau]⟩

@[simp] theorem otriOf_base (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : (otriOf A a u hau).base = 𝟙 A.base := rfl
@[simp] theorem otriOf_div (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : (otriOf A a u hau).div = a := rfl
@[simp] theorem otriOf_deg (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : (otriOf A a u hau).deg = 1 := rfl
@[simp] theorem otriOf_u (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : (otriOf A a u hau).u = u := rfl

theorem otriOf_mem (h : Hyp M) (A : Obj M) (a : M.phi.val A.base) (u : M.bmon.val A.base)
    (hau : toGpHom _ a = M.divB _ u) : otriOf A a u hau ∈ OTri (modelPre h) A :=
  ⟨((modelPre h).Base_id A).symm, rfl⟩

/-- ★`𝒪^▷(A)` の元は「底が `𝟙`・次数 1」に他ならない。 -/
theorem mem_otri_iff (h : Hyp M) {A : Obj M} (x : End A) :
    x ∈ OTri (modelPre h) A ↔ (x.base = 𝟙 A.base ∧ x.deg = 1) := by
  constructor
  · rintro ⟨hb, hl⟩
    exact ⟨hb.trans ((modelPre h).Base_id A), hl⟩
  · rintro ⟨hb, hl⟩
    refine ⟨?_, hl⟩
    show (modelPre h).Base x = (modelPre h).Base (𝟙 A)
    rw [(modelPre h).Base_id A]
    exact hb

/-- ★`𝒪^▷(A)` の元は `Div` と `u` で決まる。 -/
theorem otri_ext (h : Hyp M) {A : Obj M} {x y : End A}
    (hx : x ∈ OTri (modelPre h) A) (hy : y ∈ OTri (modelPre h) A)
    (hd : x.div = y.div) (hu : x.u = y.u) : x = y := by
  obtain ⟨hxb, hxl⟩ := (mem_otri_iff h x).mp hx
  obtain ⟨hyb, hyl⟩ := (mem_otri_iff h y).mp hy
  exact Hom.ext (hxb.trans hyb.symm) hd (hxl.trans hyl.symm) hu

/-- ★`𝒪^▷(A)` の元は `toGp (Div) = Div_B(u)` を満たす。 -/
theorem otri_gp_eq (h : Hyp M) {A : Obj M} {x : End A} (hx : x ∈ OTri (modelPre h) A) :
    toGpHom _ x.div = M.divB _ x.u := by
  obtain ⟨hxb, hxl⟩ := (mem_otri_iff h x).mp hx
  have hc := x.cond
  rw [hxb, hxl] at hc
  simp only [PNat.one_coe, one_smul, MonoidOn.gpMapOn_id] at hc
  exact add_left_cancel hc

/-- ★★**可換条件の判定** —— `𝔽_Φ` と同じ形の式が `Div` と `u` の 2 本になる。 -/
theorem otri_comm_iff (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    {α : End A} (hα : α ∈ OTri (modelPre h) A)
    {β : End B} (hβ : β ∈ OTri (modelPre h) B) :
    (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ ↔
      (M.phi.map φ.base β.div = (φ.deg : ℕ) • α.div ∧
        M.bmon.map φ.base β.u = (φ.deg : ℕ) • α.u) := by
  obtain ⟨hαb, hαl⟩ := (mem_otri_iff h α).mp hα
  obtain ⟨hβb, hβl⟩ := (mem_otri_iff h β).mp hβ
  letI := isCancelAdd_of_isIntegralMonoid _ (h.divisorial A.base).1.1
  constructor
  · intro heq
    have hd := congrArg Hom.div heq
    have hu := congrArg Hom.u heq
    rw [comp_div, comp_div, hβl, hαb] at hd
    rw [comp_u, comp_u, hβl, hαb] at hu
    simp only [PNat.one_coe, one_smul, MonoidOn.map_id] at hd hu
    refine ⟨add_right_cancel (by rw [add_comm ((φ.deg : ℕ) • α.div)]; exact hd), ?_⟩
    exact groupLike_add_right_cancel (h.bmonGroupLike A.base) (c := φ.u)
      (by rw [add_comm ((φ.deg : ℕ) • α.u)]; exact hu)
  · rintro ⟨hd, hu⟩
    refine Hom.ext ?_ ?_ ?_ ?_
    · rw [comp_base, comp_base, hβb, hαb, Category.comp_id, Category.id_comp]
    · rw [comp_div, comp_div, hβl, hαb]
      simp only [PNat.one_coe, one_smul, MonoidOn.map_id]
      rw [hd, add_comm]
    · rw [comp_deg, comp_deg, hβl, hαl, one_mul, mul_one]
    · rw [comp_u, comp_u, hβl, hαb]
      simp only [PNat.one_coe, one_smul, MonoidOn.map_id]
      rw [hu, add_comm]

/-! ## ★7. (iii)(c) の 3 条 -/

theorem phi_map_roundtrip {A B : D} (f : A ⟶ B) [IsIso f] (x : M.phi.val A) :
    M.phi.map f (M.phi.map (inv f) x) = x := by
  rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]

theorem bmon_map_roundtrip {A B : D} (f : A ⟶ B) [IsIso f] (x : M.bmon.val A) :
    M.bmon.map f (M.bmon.map (inv f) x) = x := by
  rw [← MonoidOn.map_comp, IsIso.hom_inv_id, MonoidOn.map_id]

/-- **(iii)(c)** 順方向。 -/
theorem model_otriFwd (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    (_hca : IsCoAngular (modelPre h) φ) (hst : IsPreStep (modelPre h) φ)
    (α : End A) (hα : α ∈ OTri (modelPre h) A) :
    ∃! β : End B, β ∈ OTri (modelPre h) B ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  haveI hbi : IsIso φ.base := hst.2
  have hdeg : φ.deg = 1 := hst.1
  have hαgp := otri_gp_eq h hα
  have hb : toGpHom _ (M.phi.map (inv φ.base) α.div)
      = M.divB _ (M.bmon.map (inv φ.base) α.u) := by
    rw [← MonoidOn.gpMapOn_toGpHom, hαgp, ← M.divB_nat (inv φ.base) α.u]
  refine ⟨otriOf B (M.phi.map (inv φ.base) α.div) (M.bmon.map (inv φ.base) α.u) hb,
    ⟨otriOf_mem h B _ _ hb, ?_⟩, ?_⟩
  · refine (otri_comm_iff h φ hα (otriOf_mem h B _ _ hb)).mpr ?_
    rw [hdeg]
    simp only [PNat.one_coe, one_smul, otriOf_div, otriOf_u]
    exact ⟨phi_map_roundtrip φ.base α.div, bmon_map_roundtrip φ.base α.u⟩
  · rintro β ⟨hβ, hβe⟩
    obtain ⟨hd, hu⟩ := (otri_comm_iff h φ hα hβ).mp hβe
    rw [hdeg] at hd hu
    simp only [PNat.one_coe, one_smul] at hd hu
    refine otri_ext h hβ (otriOf_mem h B _ _ hb) ?_ ?_
    · rw [otriOf_div, ← hd, ← MonoidOn.map_comp, IsIso.inv_hom_id, MonoidOn.map_id]
    · rw [otriOf_u, ← hu, ← MonoidOn.map_comp, IsIso.inv_hom_id, MonoidOn.map_id]

/-- **(iii)(c)** 逆方向。 -/
theorem model_otriBwd (h : Hyp M) {A B : Obj M} (φ : A ⟶ B)
    (_hca : IsCoAngular (modelPre h) φ) (hst : IsPreStep (modelPre h) φ)
    (β : End B) (hβ : β ∈ OTri (modelPre h) B) :
    ∃! α : End A, α ∈ OTri (modelPre h) A ∧
      (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  haveI hbi : IsIso φ.base := hst.2
  have hdeg : φ.deg = 1 := hst.1
  have hβgp := otri_gp_eq h hβ
  have ha : toGpHom _ (M.phi.map φ.base β.div) = M.divB _ (M.bmon.map φ.base β.u) := by
    rw [← MonoidOn.gpMapOn_toGpHom, hβgp, ← M.divB_nat φ.base β.u]
  refine ⟨otriOf A (M.phi.map φ.base β.div) (M.bmon.map φ.base β.u) ha,
    ⟨otriOf_mem h A _ _ ha, ?_⟩, ?_⟩
  · refine (otri_comm_iff h φ (otriOf_mem h A _ _ ha) hβ).mpr ?_
    rw [hdeg]
    simp
  · rintro α ⟨hα, hαe⟩
    obtain ⟨hd, hu⟩ := (otri_comm_iff h φ hα hβ).mp hαe
    rw [hdeg] at hd hu
    simp only [PNat.one_coe, one_smul] at hd hu
    exact otri_ext h hα (otriOf_mem h A _ _ ha) hd.symm hu.symm

/-- **(iii)(c)** 全単射は `Base(φ)` にしか依らない。 -/
theorem model_otriBase (h : Hyp M) {A B : Obj M} (φ φ' : A ⟶ B)
    (_hca : IsCoAngular (modelPre h) φ) (hst : IsPreStep (modelPre h) φ)
    (_hca' : IsCoAngular (modelPre h) φ') (hst' : IsPreStep (modelPre h) φ')
    (hbase : (modelPre h).Base φ = (modelPre h).Base φ')
    (α : End A) (hα : α ∈ OTri (modelPre h) A)
    (β : End B) (hβ : β ∈ OTri (modelPre h) B)
    (heq : (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ) :
    (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ' := by
  obtain ⟨hd, hu⟩ := (otri_comm_iff h φ hα hβ).mp heq
  rw [show φ.deg = 1 from hst.1] at hd hu
  refine (otri_comm_iff h φ' hα hβ).mpr ?_
  rw [show φ'.deg = 1 from hst'.1, show φ'.base = φ.base from hbase.symm]
  exact ⟨hd, hu⟩

/-! ## ★8. (v)(b)(c) の一意性に共通する同型

★`𝔽_Φ` と同じく **1 つの構成が (b)(c) 両方を与える**。
★★`Φ` は sharp なので「可逆元だけ違う」は「**等しい**」になり、
`𝔽_Φ` の証明より短くなる。代わりに `u` 成分の合わせが要る。 -/

theorem model_preStepIso (h : Hyp M) {A B X X' : Obj M}
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B)
    (heq : (β ≫ α : A ⟶ B) = β' ≫ α')
    (hβ : IsPreStep (modelPre h) β) (hα : IsLinear (modelPre h) α)
    (hβ' : IsPreStep (modelPre h) β') (hα' : IsLinear (modelPre h) α')
    (hdiv : β'.div = β.div) :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  haveI hbβ : IsIso β.base := hβ.2
  haveI hbβ' : IsIso β'.base := hβ'.2
  have hdβ : β.deg = 1 := hβ.1
  have hdβ' : β'.deg = 1 := hβ'.1
  have hdα : α.deg = 1 := hα
  have hdα' : α'.deg = 1 := hα'
  have hbase : β.base ≫ α.base = β'.base ≫ α'.base := congrArg Hom.base heq
  have hueq : M.bmon.map β.base α.u + β.u = M.bmon.map β'.base α'.u + β'.u := by
    have hh := congrArg Hom.u heq
    rw [comp_u, comp_u, hdα, hdα'] at hh
    simpa using hh
  have hdeq : M.phi.map β.base α.div = M.phi.map β'.base α'.div := by
    have hh := congrArg Hom.div heq
    rw [comp_div, comp_div, hdα, hdα', hdiv] at hh
    simp only [PNat.one_coe, one_smul] at hh
    letI := isCancelAdd_of_isIntegralMonoid _ (h.divisorial A.base).1.1
    exact add_right_cancel hh
  -- ★`γ` の関係式
  have hcond : ((1 : ℕ+) : ℕ) • X.cls + toGpHom _ (0 : M.phi.val X.base)
      = M.phi.gpMapOn (inv β.base ≫ β'.base) X'.cls
        + M.divB _ (M.bmon.map (inv β.base) (β'.u + bneg h β.u)) := by
    have hb := β.cond
    have hb' := β'.cond
    rw [hdβ] at hb
    rw [hdβ', hdiv] at hb'
    simp only [PNat.one_coe, one_smul] at hb hb'
    have hkey : M.phi.gpMapOn β.base X.cls + M.divB _ β.u
        = M.phi.gpMapOn β'.base X'.cls + M.divB _ β'.u := by rw [← hb, ← hb']
    have h2 := congrArg (M.phi.gpMapOn (inv β.base)) hkey
    rw [(M.phi.gpMapOn (inv β.base)).map_add, (M.phi.gpMapOn (inv β.base)).map_add,
      ← MonoidOn.gpMapOn_comp, IsIso.inv_hom_id, MonoidOn.gpMapOn_id,
      ← MonoidOn.gpMapOn_comp, ← M.divB_nat (inv β.base) β.u,
      ← M.divB_nat (inv β.base) β'.u] at h2
    rw [map_add, map_bneg, (M.divB X.base).map_add, divB_bneg]
    simp only [PNat.one_coe, one_smul, map_zero, add_zero]
    rw [eq_sub_of_add_eq h2]
    abel
  obtain ⟨g, hgb, hgd, hgdeg, hgu⟩ :
      ∃ g : X ⟶ X', g.base = inv β.base ≫ β'.base ∧ g.div = 0 ∧ g.deg = 1 ∧
        g.u = M.bmon.map (inv β.base) (β'.u + bneg h β.u) :=
    ⟨⟨inv β.base ≫ β'.base, 0, 1,
      M.bmon.map (inv β.base) (β'.u + bneg h β.u), hcond⟩, rfl, rfl, rfl, rfl⟩
  haveI hgiso : IsIso g :=
    (model_isIso_iff h g).mpr ⟨by rw [hgb]; infer_instance, hgd, hgdeg⟩
  refine ⟨asIso g, ?_, ?_⟩
  · -- `α' = γ.inv ≫ α` ⟺ `γ.hom ≫ α' = α`
    rw [Iso.eq_inv_comp]
    show (g ≫ α' : X ⟶ B) = α
    refine Hom.ext ?_ ?_ ?_ ?_
    · show g.base ≫ α'.base = α.base
      rw [hgb, Category.assoc, ← hbase, ← Category.assoc, IsIso.inv_hom_id,
        Category.id_comp]
    · show M.phi.map g.base α'.div + (α'.deg : ℕ) • g.div = α.div
      rw [hgb, hgd, MonoidOn.map_comp, ← hdeq, ← MonoidOn.map_comp, IsIso.inv_hom_id,
        MonoidOn.map_id]
      simp
    · show α'.deg * g.deg = α.deg
      rw [hgdeg, mul_one, hdα, hdα']
    · show M.bmon.map g.base α'.u + (α'.deg : ℕ) • g.u = α.u
      rw [hdα', hgb, hgu]
      simp only [PNat.one_coe, one_smul]
      rw [MonoidOn.map_comp, ← map_add, ← add_assoc, ← hueq, add_assoc, add_bneg,
        add_zero, ← MonoidOn.map_comp, IsIso.inv_hom_id, MonoidOn.map_id]
  · show β' = β ≫ g
    refine Hom.ext ?_ ?_ ?_ ?_
    · show β'.base = β.base ≫ g.base
      rw [hgb, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    · show β'.div = M.phi.map β.base g.div + (g.deg : ℕ) • β.div
      rw [hgd, hgdeg, map_zero, zero_add]
      simp [hdiv]
    · show β'.deg = g.deg * β.deg
      rw [hgdeg, one_mul, hdβ, hdβ']
    · show β'.u = M.bmon.map β.base g.u + (g.deg : ℕ) • β.u
      rw [hgu, hgdeg, bmon_map_roundtrip]
      simp only [PNat.one_coe, one_smul]
      rw [add_assoc, bneg_add, add_zero]

/-- **(v)(b)** の一意性。 -/
theorem model_preStepFactorUniq (h : Hyp M) {A B : Obj M} (X X' : Obj M)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (_hcβ : IsCoAngular (modelPre h) β) (hβ : IsPreStep (modelPre h) β)
    (hiα : IsIsometric (modelPre h) α) (hα : IsPreStep (modelPre h) α)
    (_hcβ' : IsCoAngular (modelPre h) β') (hβ' : IsPreStep (modelPre h) β')
    (hiα' : IsIsometric (modelPre h) α') (hα' : IsPreStep (modelPre h) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  refine model_preStepIso h β α β' α' heq hβ hα.1 hβ' hα'.1 ?_
  have hh := congrArg Hom.div heq
  rw [comp_div, comp_div, show α.deg = 1 from hα.1, show α'.deg = 1 from hα'.1,
    show α.div = 0 from hiα, show α'.div = 0 from hiα'] at hh
  simp only [PNat.one_coe, one_smul, map_zero, zero_add] at hh
  exact hh.symm

/-- **(v)(c)** の一意性。 -/
theorem model_preStepFactorUniq' (h : Hyp M) {A B : Obj M} (X X' : Obj M)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (hiβ : IsIsometric (modelPre h) β) (hβ : IsPreStep (modelPre h) β)
    (_hcα : IsCoAngular (modelPre h) α) (hα : IsPreStep (modelPre h) α)
    (hiβ' : IsIsometric (modelPre h) β') (hβ' : IsPreStep (modelPre h) β')
    (_hcα' : IsCoAngular (modelPre h) α') (hα' : IsPreStep (modelPre h) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom :=
  model_preStepIso h β α β' α' heq hβ hα.1 hβ' hα'.1
    (by rw [show β.div = 0 from hiβ, show β'.div = 0 from hiβ'])

/-! ## ★9. (vi) 単元を除く忠実性

★`Φ` は sharp なので零因子の差は **0**。ずれは `u` にしか出ない。 -/

theorem model_faithfulUpToUnits (h : Hyp M) {A B : Obj M} (φ ψ : A ⟶ B)
    (hb : BaseEquivalent (modelPre h) φ ψ) (hm : MetricallyEquivalent (modelPre h) φ ψ)
    (hφ : IsPreStep (modelPre h) φ) (hψ : IsPreStep (modelPre h) ψ) :
    ∃ α : End B, α ∈ OTimes (modelPre h) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  haveI hbψ : IsIso ψ.base := hψ.2
  have hbe : φ.base = ψ.base := hb
  have hme : φ.div = ψ.div := hm
  have hdφ : φ.deg = 1 := hφ.1
  have hdψ : ψ.deg = 1 := hψ.1
  -- ★`Div_B(φ.u) = Div_B(ψ.u)`
  have hgp : M.divB _ φ.u = M.divB _ ψ.u := by
    have hc1 := φ.cond
    have hc2 := ψ.cond
    rw [hdφ, hme, hbe] at hc1
    rw [hdψ] at hc2
    exact add_left_cancel (hc1.symm.trans hc2)
  have hcond : toGpHom _ (0 : M.phi.val B.base)
      = M.divB _ (M.bmon.map (inv ψ.base) (φ.u + bneg h ψ.u)) := by
    rw [M.divB_nat (inv ψ.base) (φ.u + bneg h ψ.u), map_add, divB_bneg, hgp]
    simp
  refine ⟨otriOf B 0 (M.bmon.map (inv ψ.base) (φ.u + bneg h ψ.u)) hcond,
    ⟨otriOf_mem h B _ _ hcond, ?_⟩, ?_⟩
  · refine (isUnit_iff_isIso _).mpr ((model_isIso_iff h _).mpr ⟨?_, rfl, rfl⟩)
    show IsIso (𝟙 B.base)
    infer_instance
  · refine Hom.ext ?_ ?_ ?_ ?_
    · show φ.base = ψ.base ≫ 𝟙 B.base
      rw [Category.comp_id]
      exact hbe
    · show φ.div = M.phi.map ψ.base (0 : M.phi.val B.base) + ((1 : ℕ+) : ℕ) • ψ.div
      rw [map_zero, zero_add]
      simp [hme]
    · show φ.deg = 1 * ψ.deg
      rw [one_mul, hdφ, hdψ]
    · show φ.u = M.bmon.map ψ.base (M.bmon.map (inv ψ.base) (φ.u + bneg h ψ.u))
        + ((1 : ℕ+) : ℕ) • ψ.u
      rw [bmon_map_roundtrip]
      simp only [PNat.one_coe, one_smul]
      rw [add_assoc, bneg_add, add_zero]

/-! ## ★10. (i)(b) span

★★`𝔽_Φ` では `X := A` で済んだが、model では **`α` を下げた対象**を作らねばならない。
`Gp` の元は `toGp a − toGp b` と書けるので、
`X := (A_𝒟, −toGp(b₁ + b₂))` と取れば `A` にも `B` にも pre-step が入る。 -/

theorem model_preStepSpan (h : Hyp M) (A B : Obj M)
    (α : ((modelPre h).toElem.obj A).base ⟶ ((modelPre h).toElem.obj B).base)
    (hα : IsIso α) :
    ∃ (X : Obj M) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep (modelPre h) φ),
      IsPreStep (modelPre h) ψ ∧
        α = @inv _ _ _ _ ((modelPre h).Base φ) hφ.2 ≫ (modelPre h).Base ψ := by
  obtain ⟨a, rfl⟩ : ∃ a : A.base ⟶ B.base, a = α := ⟨α, rfl⟩
  obtain ⟨a₁, b₁, h1⟩ := gp_sub_repr A.cls
  obtain ⟨a₂, b₂, h2⟩ := gp_sub_repr (M.phi.gpMapOn a B.cls)
  have hc1 : ((1 : ℕ+) : ℕ) • (-(toGpHom (M.phi.val A.base) (b₁ + b₂)))
      + toGpHom _ (a₁ + b₂)
      = M.phi.gpMapOn (𝟙 A.base) A.cls + M.divB _ (0 : M.bmon.val A.base) := by
    rw [h1, MonoidOn.gpMapOn_id, map_zero, add_zero, map_add, map_add]
    simp only [PNat.one_coe, one_smul]
    abel
  have hc2 : ((1 : ℕ+) : ℕ) • (-(toGpHom (M.phi.val A.base) (b₁ + b₂)))
      + toGpHom _ (a₂ + b₁)
      = M.phi.gpMapOn a B.cls + M.divB _ (0 : M.bmon.val A.base) := by
    rw [h2, map_zero, add_zero, map_add, map_add]
    simp only [PNat.one_coe, one_smul]
    abel
  refine ⟨(⟨A.base, -(toGpHom _ (b₁ + b₂))⟩ : Obj M),
    (⟨𝟙 A.base, a₁ + b₂, 1, 0, hc1⟩ :
      (⟨A.base, -(toGpHom _ (b₁ + b₂))⟩ : Obj M) ⟶ A),
    (⟨a, a₂ + b₁, 1, 0, hc2⟩ :
      (⟨A.base, -(toGpHom _ (b₁ + b₂))⟩ : Obj M) ⟶ B),
    ⟨rfl, ?_⟩, ⟨rfl, hα⟩, ?_⟩
  · show IsIso (𝟙 A.base)
    infer_instance
  · haveI hii : IsIso ((modelPre h).Base
        (⟨𝟙 A.base, a₁ + b₂, 1, 0, hc1⟩ :
          (⟨A.base, -(toGpHom _ (b₁ + b₂))⟩ : Obj M) ⟶ A)) := by
      show IsIso (𝟙 A.base)
      infer_instance
    exact (@IsIso.eq_inv_comp _ _ _ _ _ _ hii _ _).mpr (by
      show 𝟙 A.base ≫ a = a
      rw [Category.id_comp])

/-! ## ★11. (iv)(a) の一意性 -/

/-- ★`α` に要るのは **linear + isometric** だけ(pre-step まで要らない)。 -/
theorem model_preStepFactorUniq_lin (h : Hyp M) {A B : Obj M} (X X' : Obj M)
    (β : A ⟶ X) (α : X ⟶ B) (β' : A ⟶ X') (α' : X' ⟶ B) (heq : β ≫ α = β' ≫ α')
    (hβ : IsPreStep (modelPre h) β) (hiα : IsIsometric (modelPre h) α)
    (hlα : IsLinear (modelPre h) α)
    (hβ' : IsPreStep (modelPre h) β') (hiα' : IsIsometric (modelPre h) α')
    (hlα' : IsLinear (modelPre h) α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  refine model_preStepIso h β α β' α' heq hβ hlα hβ' hlα' ?_
  have hh := congrArg Hom.div heq
  rw [comp_div, comp_div, show α.deg = 1 from hlα, show α'.deg = 1 from hlα',
    show α.div = 0 from hiα, show α'.div = 0 from hiα'] at hh
  simp only [PNat.one_coe, one_smul, map_zero, zero_add] at hh
  exact hh.symm

theorem model_arbFactorUniq (h : Hyp M) {A B : Obj M} (X Y X' Y' : Obj M)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (heq : (γ ≫ β ≫ α : A ⟶ B) = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType (modelPre h) γ) (hβ : IsPreStep (modelPre h) β)
    (hα : IsPullBack (modelPre h) α)
    (hγ' : IsFrobeniusType (modelPre h) γ') (hβ' : IsPreStep (modelPre h) β')
    (hα' : IsPullBack (modelPre h) α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  obtain ⟨hαl, hαi⟩ := (model_isPullBack_iff h α).mp hα
  obtain ⟨hαl', hαi'⟩ := (model_isPullBack_iff h α').mp hα'
  have hdegγ : (modelPre h).degFr γ = (modelPre h).degFr γ' := by
    have hh := congrArg Hom.deg heq
    rw [comp_deg, comp_deg, comp_deg, comp_deg,
      show β.deg = 1 from hβ.1, show α.deg = 1 from hαl,
      show β'.deg = 1 from hβ'.1, show α'.deg = 1 from hαl'] at hh
    simpa using hh
  obtain ⟨εm, hεiso, hγε⟩ := model_frobDegUniq h A X X' γ γ' hγ hγ' hdegγ
  obtain ⟨hεb, -, hεn⟩ := (model_isIso_iff h εm).mp hεiso
  haveI := hεiso
  haveI : Epi γ := model_totEpi h _ _ γ
  have hcancel : (β ≫ α : X ⟶ B) = (εm ≫ β') ≫ α' := by
    refine (cancel_epi γ).mp ?_
    rw [heq, ← hγε]
    simp
  have hβ'' : IsPreStep (modelPre h) (εm ≫ β') := by
    refine ⟨?_, ?_⟩
    · show (εm ≫ β').deg = 1
      rw [comp_deg, hεn, show β'.deg = 1 from hβ'.1, mul_one]
    · haveI hb1 : IsIso εm.base := hεb
      haveI hb2 : IsIso β'.base := hβ'.2
      show IsIso ((εm ≫ β').base)
      rw [comp_base]
      infer_instance
  obtain ⟨δ, hδ1, hδ2⟩ := model_preStepFactorUniq_lin h Y Y' β α (εm ≫ β') α'
    hcancel hβ hαi hαl hβ'' hαi' hαl'
  refine ⟨δ, asIso εm, hδ1, ?_, hγε.symm⟩
  rw [← hδ2]
  simp

/-! ## ★12. (i)(c) —— `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値

★★pull-back は「次数 1・零因子 0」なので、`𝒞^pl-bk` の射は
**底の射と `u` だけ**で決まる。その `u` は `B` が group-like なので一意に定まり、
**底の圏そのもの**が復元される。 -/

theorem model_plBkEquiv (h : Hyp M) (A : Obj M) :
    (plBkOverFunctor (modelPre h) A).IsEquivalence := by
  haveI hfaith : (plBkOverFunctor (modelPre h) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : f.left.hom.base = g.left.hom.base := congrArg CommaMorphism.left hfg
    obtain ⟨hWl, hWi⟩ := (model_isPullBack_iff h W.hom.hom).mp W.hom.property
    obtain ⟨hfl, hfi⟩ := (model_isPullBack_iff h f.left.hom).mp f.left.property
    obtain ⟨hgl, hgi⟩ := (model_isPullBack_iff h g.left.hom).mp g.left.property
    have hwf : (f.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have hwg : (g.left.hom ≫ W.hom.hom) = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    have hu : f.left.hom.u = g.left.hom.u := by
      have h1 := congrArg Hom.u (hwf.trans hwg.symm)
      rw [comp_u, comp_u, hb, hWl] at h1
      simp only [PNat.one_coe, one_smul] at h1
      exact groupLike_add_right_cancel (h.bmonGroupLike _)
        (c := M.bmon.map g.left.hom.base W.hom.hom.u)
        (by rw [add_comm f.left.hom.u, add_comm g.left.hom.u]; exact h1)
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      (Hom.ext hb (hfi.trans hgi.symm) (hfl.trans hgl.symm) hu))
  haveI hfull : (plBkOverFunctor (modelPre h) A).Full := by
    constructor
    intro Z W hh
    obtain ⟨hWl, hWi⟩ := (model_isPullBack_iff h W.hom.hom).mp W.hom.property
    obtain ⟨hZl, hZi⟩ := (model_isPullBack_iff h Z.hom.hom).mp Z.hom.property
    obtain ⟨f₀, hf₀⟩ := (W.hom.property Z.left.obj).2
      ⟨(Z.hom.hom, hh.left), (Over.w hh).symm⟩
    have hp := Subtype.ext_iff.mp hf₀
    have h1 : (f₀ ≫ W.hom.hom) = Z.hom.hom := congrArg Prod.fst hp
    have h2 : f₀.base = hh.left := congrArg Prod.snd hp
    have hdeg : f₀.deg = 1 := by
      have hx := congrArg Hom.deg h1
      rw [comp_deg, hWl, one_mul, hZl] at hx
      exact hx
    have hdiv : f₀.div = 0 := by
      have hx := congrArg Hom.div h1
      rw [comp_div, hWl, hWi, hZi] at hx
      simpa using hx
    exact ⟨Over.homMk (⟨f₀, (model_isPullBack_iff h f₀).mpr ⟨hdeg, hdiv⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext h1), Over.OverMorphism.ext h2⟩
  haveI hess : (plBkOverFunctor (modelPre h) A).EssSurj := by
    constructor
    intro Y
    obtain ⟨q, hq⟩ : ∃ q : Y.left ⟶ ((modelPre h).toElem.obj A).base, q = Y.hom :=
      ⟨Y.hom, rfl⟩
    have hcond : ((1 : ℕ+) : ℕ) • (M.phi.gpMapOn q A.cls)
        + toGpHom _ (0 : M.phi.val Y.left)
        = M.phi.gpMapOn q A.cls + M.divB _ (0 : M.bmon.val Y.left) := by simp
    refine ⟨Over.mk (show (⟨(⟨Y.left, M.phi.gpMapOn q A.cls⟩ : Obj M)⟩ : PlBk (modelPre h))
        ⟶ (⟨A⟩ : PlBk (modelPre h)) from
      ⟨(⟨q, 0, 1, 0, hcond⟩ : (⟨Y.left, M.phi.gpMapOn q A.cls⟩ : Obj M) ⟶ A),
        (model_isPullBack_iff h _).mpr ⟨rfl, rfl⟩⟩), ⟨?_⟩⟩
    refine Over.isoMk (Iso.refl _) ?_
    show 𝟙 Y.left ≫ Y.hom = q
    rw [Category.id_comp, hq]
  exact ⟨hfaith, hfull, hess⟩

/-! ## ★13. `Definition 1.3` の core 21 条

原文 (FrdI p.101):
> (ii) The category C is a Frobenioid [with respect to the functor C
-/

/-- ★★★★**model Frobenioid は `FrobenioidCore` の 21 条をすべて満たす**。 -/
theorem model_frobenioidCore (h : Hyp M) : FrobenioidCore (modelPre h) where
  baseSurj := model_baseSurj h
  preStepSpan := model_preStepSpan h
  plBkEquiv := model_plBkEquiv h
  frobDegSurj := model_frobDegSurj h
  frobDegUniq := model_frobDegUniq h
  coAngularComp := model_coAngularComp h
  coAngularOfPreStep := model_coAngularOfPreStep h
  otriFwd := fun φ hca hst => model_otriFwd h φ hca hst
  otriBwd := fun φ hca hst => model_otriBwd h φ hca hst
  otriBase := fun φ φ' hca hst hca' hst' => model_otriBase h φ φ' hca hst hca' hst'
  arbFactor := model_arbFactor h
  arbFactorUniq := model_arbFactorUniq h
  pullBackLB := model_pullBackLB h
  preStepMono := model_preStepMono h
  preStepFactor := model_preStepFactor h
  preStepFactorUniq := model_preStepFactorUniq h
  preStepFactor' := model_preStepFactor' h
  preStepFactorUniq' := model_preStepFactorUniq' h
  faithfulUpToUnits := fun φ ψ hb hm _hcφ hφ _hcψ hψ =>
    model_faithfulUpToUnits h φ ψ hb hm hφ hψ
  isotropicHullExists := model_isotropicHullExists h
  isotropicClosed := model_isotropicClosed h

/-! ## ★14. (iii)(d) の 2 本の圏同値

★★`𝔽_Φ` は `Φ^char` 上に載っていたので `toChar`・`CharRel` を経由したが、
model Frobenioid は **`Φ` そのもの**の上にあるので `MLe` がそのまま使える。
代わりに **`u` 成分と `α` の整合(`cond`)**を作る手間が増える。

★★**型の注意**: `WideSubcategory` を通ると対象が `{ obj := A }.obj` の形で現れ、
`A.base` と**構文上は違う**ものになる。そこで **射の構成は素の型を取る補題に出し**、
圏同値の証明ではそれを適用するだけにする。 -/

theorem gpMapOn_injective_of_iso (h : Hyp M) {X Y : D} (f : X ⟶ Y) [IsIso f] :
    Function.Injective (M.phi.gpMapOn f) := by
  intro x y hxy
  have h1 : ∀ z : Gp (M.phi.val Y), M.phi.gpMapOn (inv f) (M.phi.gpMapOn f z) = z := by
    intro z
    rw [← MonoidOn.gpMapOn_comp, IsIso.inv_hom_id, MonoidOn.gpMapOn_id]
  rw [← h1 x, ← h1 y, hxy]

/-- ★**前置** —— `z : A ⟶ Z`、`w : A ⟶ W` が pre-step で `Div z + x = Div w` なら、
`z ≫ t = w` なる pre-step `t : Z ⟶ W` がある。 -/
theorem exists_hom_under (h : Hyp M) {A Z W : Obj M} (z : A ⟶ Z) (w : A ⟶ W)
    (hz : IsIso z.base) (hw : IsIso w.base) (hzd : z.deg = 1) (hwd : w.deg = 1)
    (x : M.phi.val A.base) (hx : z.div + x = w.div) :
    ∃ t : Z ⟶ W, z ≫ t = w ∧ t.deg = 1 ∧ IsIso t.base := by
  haveI := hz
  haveI := hw
  have hZc := z.cond
  have hWc := w.cond
  rw [hzd] at hZc
  rw [hwd] at hWc
  simp only [PNat.one_coe, one_smul] at hZc hWc
  have hz2 : M.phi.gpMapOn z.base Z.cls
      = A.cls + toGpHom (M.phi.val A.base) z.div - M.divB A.base z.u :=
    eq_sub_of_add_eq hZc.symm
  have hw2 : M.phi.gpMapOn w.base W.cls
      = A.cls + toGpHom (M.phi.val A.base) w.div - M.divB A.base w.u :=
    eq_sub_of_add_eq hWc.symm
  have hxg : toGpHom (M.phi.val A.base) z.div + toGpHom (M.phi.val A.base) x
      = toGpHom (M.phi.val A.base) w.div := by
    rw [← map_add, hx]
  have hcond : ((1 : ℕ+) : ℕ) • Z.cls + toGpHom _ (M.phi.map (inv z.base) x)
      = M.phi.gpMapOn (inv z.base ≫ w.base) W.cls
        + M.divB _ (M.bmon.map (inv z.base) (w.u + bneg h z.u)) := by
    refine gpMapOn_injective_of_iso h z.base ?_
    simp only [PNat.one_coe, one_smul]
    rw [(M.phi.gpMapOn z.base).map_add, (M.phi.gpMapOn z.base).map_add,
      MonoidOn.gpMapOn_toGpHom, phi_map_roundtrip,
      ← MonoidOn.gpMapOn_comp, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp,
      ← M.divB_nat z.base _, bmon_map_roundtrip,
      hz2, hw2, (M.divB A.base).map_add, divB_bneg, ← hxg]
    abel
  refine ⟨⟨inv z.base ≫ w.base, M.phi.map (inv z.base) x, 1,
    M.bmon.map (inv z.base) (w.u + bneg h z.u), hcond⟩, ?_, rfl, ?_⟩
  · refine Hom.ext ?_ ?_ ?_ ?_
    · show z.base ≫ (inv z.base ≫ w.base) = w.base
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    · show M.phi.map z.base (M.phi.map (inv z.base) x) + ((1 : ℕ+) : ℕ) • z.div = w.div
      rw [phi_map_roundtrip]
      simp only [PNat.one_coe, one_smul]
      rw [add_comm]
      exact hx
    · show (1 : ℕ+) * z.deg = w.deg
      rw [one_mul, hzd, hwd]
    · show M.bmon.map z.base (M.bmon.map (inv z.base) (w.u + bneg h z.u))
        + ((1 : ℕ+) : ℕ) • z.u = w.u
      rw [bmon_map_roundtrip]
      simp only [PNat.one_coe, one_smul]
      rw [add_assoc, bneg_add, add_zero]
  · show IsIso (inv z.base ≫ w.base)
    infer_instance

/-- ★`Div` がちょうど `a` の pre-step を `A` から出す。 -/
theorem exists_hom_div (h : Hyp M) (A : Obj M) (a : M.phi.val A.base) :
    ∃ (Z : Obj M) (z : A ⟶ Z), z.div = a ∧ z.deg = 1 ∧ IsIso z.base := by
  refine ⟨⟨A.base, A.cls + toGpHom _ a⟩,
    ⟨𝟙 A.base, a, 1, 0, by simp⟩, rfl, rfl, ?_⟩
  show IsIso (𝟙 A.base)
  infer_instance

/-- ★**後置** —— `z : Z ⟶ A`、`w : W ⟶ A` が pre-step で
`Base(t)^*(Div w) + y = Div z` なら、`t ≫ w = z` なる pre-step `t : Z ⟶ W` がある。 -/
theorem exists_hom_over (h : Hyp M) {A Z W : Obj M} (z : Z ⟶ A) (w : W ⟶ A)
    (hz : IsIso z.base) (hw : IsIso w.base) (hzd : z.deg = 1) (hwd : w.deg = 1)
    (y : M.phi.val Z.base)
    (hy : M.phi.map (z.base ≫ inv w.base) w.div + y = z.div) :
    ∃ t : Z ⟶ W, t ≫ w = z ∧ t.deg = 1 ∧ IsIso t.base := by
  haveI := hz
  haveI := hw
  have hZc := z.cond
  have hWc := w.cond
  rw [hzd] at hZc
  rw [hwd] at hWc
  simp only [PNat.one_coe, one_smul] at hZc hWc
  have hyg : toGpHom (M.phi.val Z.base) (M.phi.map (z.base ≫ inv w.base) w.div)
      + toGpHom (M.phi.val Z.base) y = toGpHom (M.phi.val Z.base) z.div := by
    rw [← map_add, hy]
  have hcond : ((1 : ℕ+) : ℕ) • Z.cls + toGpHom _ y
      = M.phi.gpMapOn (z.base ≫ inv w.base) W.cls
        + M.divB _ (z.u + bneg h (M.bmon.map (z.base ≫ inv w.base) w.u)) := by
    have hWt := congrArg (M.phi.gpMapOn (z.base ≫ inv w.base)) hWc
    rw [(M.phi.gpMapOn (z.base ≫ inv w.base)).map_add,
      (M.phi.gpMapOn (z.base ≫ inv w.base)).map_add, MonoidOn.gpMapOn_toGpHom,
      ← MonoidOn.gpMapOn_comp, Category.assoc, IsIso.inv_hom_id, Category.comp_id,
      ← M.divB_nat (z.base ≫ inv w.base) w.u] at hWt
    have hW3 : M.phi.gpMapOn (z.base ≫ inv w.base) W.cls
        = M.phi.gpMapOn z.base A.cls
          + M.divB Z.base (M.bmon.map (z.base ≫ inv w.base) w.u)
          - toGpHom (M.phi.val Z.base) (M.phi.map (z.base ≫ inv w.base) w.div) :=
      eq_sub_of_add_eq hWt
    have hZ3 : M.phi.gpMapOn z.base A.cls
        = Z.cls + toGpHom (M.phi.val Z.base) z.div - M.divB Z.base z.u :=
      eq_sub_of_add_eq hZc.symm
    rw [hW3, hZ3, (M.divB Z.base).map_add, divB_bneg, ← hyg]
    simp only [PNat.one_coe, one_smul]
    abel
  refine ⟨⟨z.base ≫ inv w.base, y, 1,
    z.u + bneg h (M.bmon.map (z.base ≫ inv w.base) w.u), hcond⟩, ?_, rfl, ?_⟩
  · refine Hom.ext ?_ ?_ ?_ ?_
    · show (z.base ≫ inv w.base) ≫ w.base = z.base
      rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id]
    · show M.phi.map (z.base ≫ inv w.base) w.div + (w.deg : ℕ) • y = z.div
      rw [hwd]
      simp only [PNat.one_coe, one_smul]
      exact hy
    · show w.deg * 1 = z.deg
      rw [mul_one, hzd, hwd]
    · show M.bmon.map (z.base ≫ inv w.base) w.u
        + (w.deg : ℕ) • (z.u + bneg h (M.bmon.map (z.base ≫ inv w.base) w.u)) = z.u
      rw [hwd]
      simp only [PNat.one_coe, one_smul]
      rw [← add_assoc, add_comm (M.bmon.map (z.base ≫ inv w.base) w.u) z.u,
        add_assoc, add_bneg, add_zero]
  · show IsIso (z.base ≫ inv w.base)
    infer_instance

/-- ★co-angular pre-step は合成で閉じる。 -/
instance model_coaPre_mult (h : Hyp M) :
    MorphismProperty.IsMultiplicative (coaPreProp (modelPre h)) :=
  coaPreProp_isMultiplicative (modelPre h) (model_coAngularComp h)

/-- **(iii)(d)** コスライス側 `_A(𝒞^coa-pre) → Order(Φ(A))` は圏同値。 -/
theorem model_coaPreUnderEquiv (h : Hyp M) (A : Obj M) :
    (coaPreUnderFunctor (modelPre h) A).IsEquivalence := by
  haveI hfaith : (coaPreUnderFunctor (modelPre h) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    have h2 : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w g)
    haveI : Epi Z.hom.hom := model_totEpi h _ _ _
    exact Under.UnderMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_epi Z.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreUnderFunctor (modelPre h) A).Full := by
    constructor
    intro Z W hle
    obtain ⟨x, hx⟩ := leOfHom hle
    obtain ⟨t, ht, htd, htb⟩ := exists_hom_under h Z.hom.hom W.hom.hom
      Z.hom.property.2.2 W.hom.property.2.2 Z.hom.property.2.1 W.hom.property.2.1 x hx
    haveI := Preorder.subsingleton_hom
      ((coaPreUnderFunctor (modelPre h) A).obj Z) ((coaPreUnderFunctor (modelPre h) A).obj W)
    exact ⟨Under.homMk (⟨t, ⟨model_coAngular h _, htd, htb⟩⟩ : Z.right ⟶ W.right)
      (InducedWideCategory.Hom.ext (by simp only [WideSubcategory.comp_def]; exact ht)),
      Subsingleton.elim _ _⟩
  haveI hess : (coaPreUnderFunctor (modelPre h) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a, rfl⟩ : ∃ a : M.phi.val ((modelPre h).toElem.obj A).base, toOrderCat a = c :=
      ⟨c, rfl⟩
    obtain ⟨Z, z, hzdiv, hzdeg, hzb⟩ := exists_hom_div h A a
    exact ⟨Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp (modelPre h)))
        ⟶ (⟨Z⟩ : WideSubcategory (coaPreProp (modelPre h))) from
      ⟨z, ⟨model_coAngular h _, hzdeg, hzb⟩⟩),
      ⟨eqToIso (congrArg toOrderCat hzdiv)⟩⟩
  exact ⟨hfaith, hfull, hess⟩

/-- **(iii)(d)** スライス側 `(𝒞^coa-pre)_A → Order(Φ(A))^opp` は圏同値。 -/
theorem model_coaPreOverEquiv (h : Hyp M) (A : Obj M) :
    (coaPreOverFunctor (modelPre h) A).IsEquivalence := by
  haveI hfaith : (coaPreOverFunctor (modelPre h) A).Faithful := by
    constructor
    intro Z W f g _
    have h1 : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have h2 : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : Mono W.hom.hom :=
      model_preStepMono h W.hom.hom ⟨W.hom.property.2.1, W.hom.property.2.2⟩
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext
      ((cancel_mono W.hom.hom).mp (h1.trans h2.symm)))
  haveI hfull : (coaPreOverFunctor (modelPre h) A).Full := by
    constructor
    intro Z W hle
    haveI hz : IsIso Z.hom.hom.base := Z.hom.property.2.2
    haveI hw : IsIso W.hom.hom.base := W.hom.property.2.2
    obtain ⟨x0, hx0⟩ := leOfHom hle.unop
    obtain ⟨x, rfl⟩ : ∃ x : M.phi.val A.base, x = x0 := ⟨x0, rfl⟩
    have hx' : M.phi.map (inv W.hom.hom.base) W.hom.hom.div + x
        = M.phi.map (inv Z.hom.hom.base) Z.hom.hom.div := hx0
    have hy : M.phi.map (Z.hom.hom.base ≫ inv W.hom.hom.base) W.hom.hom.div
        + M.phi.map Z.hom.hom.base x = Z.hom.hom.div := by
      have h3 := congrArg (M.phi.map Z.hom.hom.base) hx'
      rw [map_add, phi_map_roundtrip, ← MonoidOn.map_comp] at h3
      exact h3
    obtain ⟨t, ht, htd, htb⟩ := exists_hom_over h Z.hom.hom W.hom.hom
      Z.hom.property.2.2 W.hom.property.2.2 Z.hom.property.2.1 W.hom.property.2.1 _ hy
    haveI := Preorder.subsingleton_hom
      ((coaPreOverFunctor (modelPre h) A).obj Z).unop
      ((coaPreOverFunctor (modelPre h) A).obj W).unop
    exact ⟨Over.homMk (⟨t, ⟨model_coAngular h _, htd, htb⟩⟩ : Z.left ⟶ W.left)
      (InducedWideCategory.Hom.ext (by simp only [WideSubcategory.comp_def]; exact ht)),
      Subsingleton.elim _ _⟩
  haveI hess : (coaPreOverFunctor (modelPre h) A).EssSurj := by
    constructor
    intro c
    obtain ⟨a0, rfl⟩ :
        ∃ a0 : M.phi.val ((modelPre h).toElem.obj A).base,
          Opposite.op (toOrderCat a0) = c :=
      ⟨c.unop, rfl⟩
    obtain ⟨a, rfl⟩ : ∃ a : M.phi.val A.base, a = a0 := ⟨a0, rfl⟩
    have hcond : ((1 : ℕ+) : ℕ) • (A.cls - toGpHom (M.phi.val A.base) a)
        + toGpHom _ a = M.phi.gpMapOn (𝟙 A.base) A.cls + M.divB _ (0 : M.bmon.val A.base) := by
      simp
    refine ⟨Over.mk (show (⟨(⟨A.base, A.cls - toGpHom _ a⟩ : Obj M)⟩ :
        WideSubcategory (coaPreProp (modelPre h)))
        ⟶ (⟨A⟩ : WideSubcategory (coaPreProp (modelPre h))) from
      ⟨(⟨𝟙 A.base, a, 1, 0, hcond⟩ : (⟨A.base, A.cls - toGpHom _ a⟩ : Obj M) ⟶ A),
        ⟨model_coAngular h _, rfl, ?_⟩⟩), ⟨eqToIso ?_⟩⟩
    · show IsIso (𝟙 A.base)
      infer_instance
    · refine congrArg Opposite.op (congrArg toOrderCat ?_)
      show M.phi.map (inv (𝟙 A.base)) a = a
      rw [IsIso.inv_id, MonoidOn.map_id]
  exact ⟨hfaith, hfull, hess⟩

/-- ★★★★★**[FrdI] Theorem 5.2, (ii)** —— model Frobenioid は **Frobenioid** である。

原文 (FrdI p.101):
> (ii) The category C is a Frobenioid [with respect to the functor C
-/
theorem model_frobenioid (h : Hyp M) : Frobenioid (modelPre h) where
  core := model_frobenioidCore h
  coaPreUnderEquiv := model_coaPreUnderEquiv h
  coaPreOverEquiv := model_coaPreOverEquiv h

/-! ## ★15. Frobenius-normalized 型

★`𝒪^▷(A)` の元の冪は `Div` と `u` を `n` 倍するだけなので、
`φ ≫ α^{deg φ} = α ≫ φ` は両辺とも `Div(φ) + n·Div(α)` になる。 -/

/-- ★`𝒪^▷(A)` の元の冪の 4 成分。 -/
theorem otri_pow (h : Hyp M) {A : Obj M} {α : End A} (hα : α ∈ OTri (modelPre h) A) (n : ℕ) :
    ((α ^ n : End A).base = 𝟙 A.base ∧ (α ^ n : End A).deg = 1) ∧
      (α ^ n : End A).div = n • α.div ∧ (α ^ n : End A).u = n • α.u := by
  obtain ⟨hb, hd⟩ := (mem_otri_iff h α).mp hα
  induction n with
  | zero => exact ⟨⟨rfl, rfl⟩, by simp, by simp⟩
  | succ k ih =>
      obtain ⟨⟨hbk, hdk⟩, hdivk, huk⟩ := ih
      have hmul : (α ^ (k + 1) : End A) = (α : A ⟶ A) ≫ (α ^ k : End A) := by
        rw [pow_succ]
        rfl
      refine ⟨⟨?_, ?_⟩, ?_, ?_⟩
      · rw [hmul]
        show α.base ≫ (α ^ k : End A).base = 𝟙 A.base
        rw [hb, hbk, Category.id_comp]
      · rw [hmul]
        show (α ^ k : End A).deg * α.deg = 1
        rw [hdk, hd, mul_one]
      · rw [hmul]
        show M.phi.map α.base (α ^ k : End A).div + ((α ^ k : End A).deg : ℕ) • α.div
          = (k + 1) • α.div
        rw [hb, hdk, hdivk, MonoidOn.map_id]
        simp [add_smul, add_comm]
      · rw [hmul]
        show M.bmon.map α.base (α ^ k : End A).u + ((α ^ k : End A).deg : ℕ) • α.u
          = (k + 1) • α.u
        rw [hb, hdk, huk, MonoidOn.map_id]
        simp [add_smul, add_comm]

/-- ★★`𝒞` は **Frobenius-normalized 型**。 -/
theorem model_frobNormalizedType (h : Hyp M) : IsOfFrobeniusNormalizedType (modelPre h) := by
  intro A φ hφb α hα
  have hφbase : φ.base = 𝟙 A.base := by
    have h0 : (modelPre h).Base φ = (modelPre h).Base (𝟙 A) := hφb
    rw [(modelPre h).Base_id A] at h0
    exact h0
  obtain ⟨hαb, hαd⟩ := (mem_otri_iff h α).mp hα
  obtain ⟨⟨hpb, hpd⟩, hpdiv, hpu⟩ := otri_pow h hα ((modelPre h).degFr φ : ℕ)
  refine Hom.ext ?_ ?_ ?_ ?_
  · show φ.base ≫ (α ^ ((modelPre h).degFr φ : ℕ) : End A).base = α.base ≫ φ.base
    rw [hpb, hαb, Category.comp_id, Category.id_comp]
  · show M.phi.map φ.base (α ^ ((modelPre h).degFr φ : ℕ) : End A).div
        + ((α ^ ((modelPre h).degFr φ : ℕ) : End A).deg : ℕ) • φ.div
      = M.phi.map α.base φ.div + (φ.deg : ℕ) • α.div
    rw [hpd, hpdiv, hαb, hφbase, MonoidOn.map_id, MonoidOn.map_id]
    simp [add_comm]
  · show (α ^ ((modelPre h).degFr φ : ℕ) : End A).deg * φ.deg = φ.deg * α.deg
    rw [hpd, hαd, one_mul, mul_one]
  · show M.bmon.map φ.base (α ^ ((modelPre h).degFr φ : ℕ) : End A).u
        + ((α ^ ((modelPre h).degFr φ : ℕ) : End A).deg : ℕ) • φ.u
      = M.bmon.map α.base φ.u + (φ.deg : ℕ) • α.u
    rw [hpd, hpu, hαb, hφbase, MonoidOn.map_id, MonoidOn.map_id]
    simp [add_comm]

/-- ★★`𝒞` は **quasi-isotropic 型**(isotropic 型だから)。 -/
theorem model_quasiIsotropicType (h : Hyp M) :
    IsOfQuasiIsotropicType (Obj M) (modelPre h) :=
  isOfQuasiIsotropicType_of_isOfIsotropicType (modelPre h) (model_frobenioidCore h)
    (model_isotropicType h)

/-- ★★`𝒞` は **Frobenius-isotropic 型**(すべての対象が isotropic)。 -/
theorem model_frobIsotropicType (h : Hyp M) : IsOfFrobeniusIsotropicType (modelPre h) :=
  fun A => ⟨A, 𝟙 A, isFrobeniusType_of_isIso (modelPre h) (𝟙 A), model_isotropic h A⟩

/-! ## ★16. (iii) standard 型の判定

原文 (FrdI p.101):
> (iii) C is of standard type if and only if the following conditions are satisfied:

★★**model Frobenioid は quasi-isotropic・Frobenius-isotropic・Frobenius-normalized を
無条件に満たす**ので、`standard 型`に残るのは原文の (a)(b)(c) の 3 条だけになる。

★原文 (a) の「`Φ` が零モノイド」は、`Φ` が divisorial(したがって sharp)なので
**`Φ` が group-like であること**と同値であり、我々の `groupLikeCompact` の前件そのもの。 -/

theorem model_standardType_iff (h : Hyp M) :
    IsOfStandardType D (Obj M) (modelPre h) (model_frobenioidCore h) ↔
      ((IsOfGroupLikeType (modelPre h) → ∃ A : Istr (modelPre h),
          IsFrobeniusCompact (istrPre (modelPre h) (model_frobenioidCore h)) A) ∧
        IsOfFSMFFType D ∧ MonoidOn.IsNonDilatingOn M.phi) := by
  constructor
  · intro hs
    exact ⟨hs.groupLikeCompact, hs.baseFSMFF, hs.phiNonDilating⟩
  · rintro ⟨h1, h2, h3⟩
    exact
      { quasiIsotropic := model_quasiIsotropicType h
        frobIsotropic := model_frobIsotropicType h
        groupLikeCompact := h1
        frobNormalized := model_frobNormalizedType h
        baseFSMFF := h2
        phiNonDilating := h3 }

/-! ## ★出典 -/

/-- ★★★★★**[FrdI] Theorem 5.2, (ii)** —— model Frobenioid は
**isotropic 型の Frobenioid** である。 -/
def model_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (ii) — Frobenioid かつ isotropic 型",
    sectionId := "frdi-thm-5-2" }

/-- ★★★**[FrdI] Theorem 5.2, (iii)** —— standard 型の判定(前半)。 -/
def model_standardType_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iii) — standard 型の判定",
    sectionId := "frdi-thm-5-2" }

end ModelData

end ABC3.Found.FrdI
