/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52
import ABC3.Found.FrdI.Prop15

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

end ModelData

end ABC3.Found.FrdI
