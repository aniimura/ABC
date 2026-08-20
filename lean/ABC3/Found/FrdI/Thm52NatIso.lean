/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52Model
import ABC3.Found.FrdI.Thm52Birat
import ABC3.Found.FrdI.Prop44Ker

/-!
# [FrdI] Theorem 5.2, (ii) の残り —— `𝒪^×(−^birat) ≅ B` と model 型

原文 (FrdI p.101):
> (ii) The category C is a Frobenioid [with respect to the functor C → FΦ] of isotropic and

★★`Theorem 5.2, (ii)` には `Thm52Frob.lean` で閉じていない主張が 2 つ残っていた:

1. **model 型**(= pre-model 型 ∧ **birationally Frobenius-normalized 型**)。
   前者は `model_isPreModelType` で済んでいたが、後者が未実装だった。
2. **`𝒞^birat` に伴う `𝒪^×(−)` と `B` の自然同型**(`Div_B` と両立)。

★★★どちらも **`Hom^birat` の第 4 の不変量 `u`** を作れば出る:
```
biratU ([a]⁻¹ ≫ [f]) := B(Base a)⁻¹ (u_f − deg(f)·u_a) ∈ B(Base A)^gp
```
これは `Prop44Gp.lean` の `biratDivGp`(`Div^gp`)と**まったく同じ形**である
——model Frobenioid では `u` の合成則が `Div` の合成則と同じ形だから。

★したがって `sliceDivGpOf` 系の機械をそのまま `B` に写して作った。
★★そのうえで **4 不変量 `(Base, deg, Div^gp, u)` による忠実性**
(`model_birat_faithful`)が出て、1 は `birat_frobNormalized_of_unitTrivial`
と同じ 4 行の計算になり、2 は `κ := gpEquiv ∘ biratU` の全単射性になる。

★`B` が **group-like**(`Hyp.bmonGroupLike`、原文の仮定)なので
`Gp (B(A)) ≅ B(A)`(`gpEquivOfGroupLike`)であり、`biratU` を `B` へ降ろせる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w

/-! ## ★1. group-like なモノイドと群化 -/

/-- ★group-like なモノイドでの「逆元」。 -/
noncomputable def negOfGroupLike {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) (a : N) : N :=
  ((-(((isGroupLike_iff N).mp hN a).choose) : AddUnits N) : N)

theorem negOfGroupLike_add {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) (a : N) :
    negOfGroupLike hN a + a = 0 := by
  set u : AddUnits N := ((isGroupLike_iff N).mp hN a).choose with hu
  have hspec : (u : N) = a := ((isGroupLike_iff N).mp hN a).choose_spec
  show ((-u : AddUnits N) : N) + a = 0
  rw [← hspec]
  exact AddUnits.neg_add u

theorem toGp_injective_of_groupLike {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) :
    Function.Injective (toGp N) := by
  intro a b hab
  obtain ⟨c, hc⟩ := toGp_eq_iff.mp hab
  exact ModelData.groupLike_add_right_cancel hN hc

theorem toGp_surjective_of_groupLike {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) :
    Function.Surjective (toGp N) := by
  intro x
  refine AddLocalization.induction_on x ?_
  rintro ⟨a, b⟩
  refine ⟨a + negOfGroupLike hN (b : N), ?_⟩
  show AddLocalization.mk (a + negOfGroupLike hN (b : N)) 0 = AddLocalization.mk a b
  rw [AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
  refine ⟨⟨0, trivial⟩, ?_⟩
  show (0 : N) + ((b : N) + (a + negOfGroupLike hN (b : N))) = (0 : N) + ((0 : N) + a)
  rw [zero_add, zero_add, zero_add, add_comm (a : N) (negOfGroupLike hN (b : N)),
    ← add_assoc, add_comm (b : N) (negOfGroupLike hN (b : N)), negOfGroupLike_add, zero_add]

/-- ★★group-like なモノイドは自分の群化と同型。 -/
noncomputable def gpEquivOfGroupLike {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) :
    Gp N ≃+ N :=
  (AddEquiv.ofBijective (toGpHom N)
    ⟨toGp_injective_of_groupLike hN, toGp_surjective_of_groupLike hN⟩).symm

theorem gpEquivOfGroupLike_toGp {N : Type w} [AddCommMonoid N] (hN : IsGroupLike N) (a : N) :
    gpEquivOfGroupLike hN (toGp N a) = a :=
  (AddEquiv.ofBijective (toGpHom N)
    ⟨toGp_injective_of_groupLike hN, toGp_surjective_of_groupLike hN⟩).symm_apply_apply a

theorem gpEquivOfGroupLike_map {N N' : Type w} [AddCommMonoid N] [AddCommMonoid N']
    (hN : IsGroupLike N) (hN' : IsGroupLike N') (φ : N →+ N') (z : Gp N) :
    gpEquivOfGroupLike hN' (gpMap _ φ z) = φ (gpEquivOfGroupLike hN z) := by
  obtain ⟨y, rfl⟩ := toGp_surjective_of_groupLike hN z
  rw [gpMap_toGp, gpEquivOfGroupLike_toGp, gpEquivOfGroupLike_toGp]

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

namespace ModelData

/-! ## ★2. `Hom^birat` の第 4 の不変量 `u`

★`Prop44Gp.lean` の `sliceDivGpOf` → `biratDivGp` をそのまま `B` に写す。 -/

/-- ★`u` 成分の切り出し(`sliceDivGpOf` の `B` 版)。 -/
noncomputable def sliceUGpOf {A B A' : Obj M} (a : A' ⟶ A) (_ha : IsIso a.base)
    (φ : A' ⟶ B) : Gp (M.bmon.val A.base) :=
  haveI := _ha
  gpMap _ (M.bmon.map (inv a.base))
    (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u)

theorem sliceUGpOf_eq {A B A' : Obj M} (a : A' ⟶ A) (ha : IsIso a.base) (φ : A' ⟶ B) :
    sliceUGpOf a ha φ = haveI := ha
      gpMap _ (M.bmon.map (inv a.base))
        (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u) := rfl

theorem sliceUGpOf_congr {A B A' : Obj M} {a a' : A' ⟶ A} (h : a = a') (ha : IsIso a.base)
    (ha' : IsIso a'.base) (φ : A' ⟶ B) :
    sliceUGpOf a ha φ = sliceUGpOf a' ha' φ := by
  subst h; rfl

variable (M) in
theorem gpMap_bmon_comp {A B E : D} (β : A ⟶ B) (α : B ⟶ E) (x : Gp (M.bmon.val E)) :
    gpMap _ (M.bmon.map β) (gpMap _ (M.bmon.map α) x) = gpMap _ (M.bmon.map (β ≫ α)) x := by
  have h : (M.bmon.map β).comp (M.bmon.map α) = M.bmon.map (β ≫ α) := by
    ext y
    exact (M.bmon.map_comp α β y).symm
  rw [← h, gpMap_comp]
  rfl

theorem gpMap_bmon_id {X : D} (x : Gp (M.bmon.val X)) :
    gpMap _ (M.bmon.map (𝟙 X)) x = x := by
  have hid : M.bmon.map (𝟙 X) = AddMonoidHom.id _ := by ext z; exact M.bmon.map_id _ z
  rw [hid]; simp [gpMap_id]

theorem gpMap_bmon_inv_left {X Y : D} (g : X ⟶ Y) [IsIso g] (x : Gp (M.bmon.val X)) :
    gpMap _ (M.bmon.map g) (gpMap _ (M.bmon.map (inv g)) x) = x := by
  rw [gpMap_bmon_comp, IsIso.hom_inv_id]; exact gpMap_bmon_id x

/-- ★★★**不変量である** —— 遷移 `c` で代表元を取り替えても値が変わらない。 -/
theorem sliceUGpOf_precomp {A B A' A'' : Obj M} (a : A' ⟶ A) (ha : IsIso a.base)
    (has : a.deg = 1) (c : A'' ⟶ A') (_hc : IsIso c.base)
    (hca : IsIso (c ≫ a).base) (φ : A' ⟶ B) (hcs : c.deg = 1) :
    sliceUGpOf (c ≫ a) hca (c ≫ φ) = sliceUGpOf a ha φ := by
  haveI := ha
  haveI := hca
  have hnum : toGp _ (c ≫ φ).u - (((c ≫ φ).deg : ℕ+) : ℕ) • toGp _ (c ≫ a).u
      = gpMap _ (M.bmon.map c.base)
        (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u) := by
    show toGp _ (M.bmon.map c.base φ.u + ((φ.deg : ℕ+) : ℕ) • c.u)
        - (((φ.deg * c.deg : ℕ+)) : ℕ)
          • toGp _ (M.bmon.map c.base a.u + ((a.deg : ℕ+) : ℕ) • c.u)
      = _
    rw [hcs, has, mul_one, PNat.one_coe, one_smul,
      toGp_add, toGp_add, toGp_nsmul, map_sub, gpMap_toGp, map_nsmul, gpMap_toGp]
    simp only [smul_add]
    abel
  have hkey : inv ((c ≫ a).base) ≫ c.base = inv a.base := by
    show inv (c.base ≫ a.base) ≫ c.base = inv a.base
    rw [IsIso.inv_comp_eq, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  rw [sliceUGpOf_eq, sliceUGpOf_eq, hnum, gpMap_bmon_comp]
  exact congrArg (fun t => gpMap _ (M.bmon.map t)
    (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u)) hkey

/-- ★★**代表元によらない**。 -/
theorem sliceUGpOf_eq_of_mk_eq (h : Hyp M) {A B : Obj M}
    {Z W : IdxBirat (modelPre h) (model_frobenioid h) A}
    {φ : Z.unop.left.obj ⟶ B} {ψ : W.unop.left.obj ⟶ B}
    (heq : HomBirat.mk Z φ = HomBirat.mk W ψ) :
    sliceUGpOf Z.unop.hom.hom Z.unop.hom.property.2.2 φ
      = sliceUGpOf W.unop.hom.hom W.unop.hom.property.2.2 ψ := by
  obtain ⟨V, u, v, huv⟩ := HomBirat.eq_iff.mp heq
  have hZ : u.unop.left.hom ≫ Z.unop.hom.hom = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj (modelPre h) (model_frobenioid h) A => t.hom)
      (Over.w u.unop)
  have hW : v.unop.left.hom ≫ W.unop.hom.hom = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj (modelPre h) (model_frobenioid h) A => t.hom)
      (Over.w v.unop)
  have hcaZ : IsIso ((u.unop.left.hom ≫ Z.unop.hom.hom).base) := by
    rw [hZ]; exact V.unop.hom.property.2.2
  have hcaW : IsIso ((v.unop.left.hom ≫ W.unop.hom.hom).base) := by
    rw [hW]; exact V.unop.hom.property.2.2
  have h1 : sliceUGpOf V.unop.hom.hom V.unop.hom.property.2.2 (u.unop.left.hom ≫ φ)
      = sliceUGpOf Z.unop.hom.hom Z.unop.hom.property.2.2 φ :=
    (sliceUGpOf_congr hZ hcaZ V.unop.hom.property.2.2 (u.unop.left.hom ≫ φ)).symm.trans
      (sliceUGpOf_precomp Z.unop.hom.hom Z.unop.hom.property.2.2
        Z.unop.hom.property.2.1 u.unop.left.hom u.unop.left.property.2.2 hcaZ φ
        u.unop.left.property.2.1)
  have h2 : sliceUGpOf V.unop.hom.hom V.unop.hom.property.2.2 (v.unop.left.hom ≫ ψ)
      = sliceUGpOf W.unop.hom.hom W.unop.hom.property.2.2 ψ :=
    (sliceUGpOf_congr hW hcaW V.unop.hom.property.2.2 (v.unop.left.hom ≫ ψ)).symm.trans
      (sliceUGpOf_precomp W.unop.hom.hom W.unop.hom.property.2.2
        W.unop.hom.property.2.1 v.unop.left.hom v.unop.left.property.2.2 hcaW ψ
        v.unop.left.property.2.1)
  rw [← h1, ← h2, huv]

/-- ★★★**`Hom^birat` の `u` 成分**。 -/
noncomputable def biratU (h : Hyp M) {A B : Obj M}
    (f : HomBirat (modelPre h) (model_frobenioid h) A B) : Gp (M.bmon.val A.base) :=
  sliceUGpOf (HomBirat.exists_rep f).choose.unop.hom.hom
    (HomBirat.exists_rep f).choose.unop.hom.property.2.2
    (HomBirat.exists_rep f).choose_spec.choose

@[simp] theorem biratU_mk (h : Hyp M) {A B : Obj M}
    (Z : IdxBirat (modelPre h) (model_frobenioid h) A) (φ : Z.unop.left.obj ⟶ B) :
    biratU h (HomBirat.mk Z φ)
      = sliceUGpOf Z.unop.hom.hom Z.unop.hom.property.2.2 φ :=
  sliceUGpOf_eq_of_mk_eq h (HomBirat.exists_rep (HomBirat.mk Z φ)).choose_spec.choose_spec

/-! ## ★3. 4 不変量による忠実性 -/

/-- ★★★★**model の `𝒞^birat` は 4 不変量 `(Base, deg, Div^gp, u)` で忠実**。

★`𝒞^model` の射は 4 成分そのものなので、`𝒞` の側の忠実性は**定義そのもの**
(`ModelData.Hom.ext`)である。 -/
theorem model_birat_faithful (h : Hyp M) {A B : Obj M}
    (f g : HomBirat (modelPre h) (model_frobenioid h) A B)
    (hb : biratBase f = biratBase g) (hd : biratDeg f = biratDeg g)
    (hz : biratDivGp f = biratDivGp g) (hu : biratU h f = biratU h g) : f = g := by
  obtain ⟨W, φ, ψ, hφ, hψ⟩ := HomBirat.exists_rep_pair f g
  subst hφ; subst hψ
  haveI hba : IsIso ((modelPre h).Base W.unop.hom.hom) := W.unop.hom.property.2.2
  haveI hba2 : IsIso (W.unop.hom.hom.base) := W.unop.hom.property.2.2
  rw [biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq] at hb
  rw [biratDeg_mk, biratDeg_mk] at hd
  have hd' : φ.deg = ψ.deg := hd
  rw [biratDivGp_mk, biratDivGp_mk, sliceDivGpOf_eq, sliceDivGpOf_eq, hd] at hz
  rw [biratU_mk, biratU_mk, sliceUGpOf_eq, sliceUGpOf_eq, hd'] at hu
  have hbase : φ.base = ψ.base := (cancel_epi (inv ((modelPre h).Base W.unop.hom.hom))).mp hb
  have hzz := congrArg (gpMap _ (M.phi.map ((modelPre h).Base W.unop.hom.hom))) hz
  rw [gpMap_phi_inv_left, gpMap_phi_inv_left] at hzz
  have hdiv : φ.div = ψ.div := ((modelPre h).divisorial _).1.1 (sub_left_injective hzz)
  have huu := congrArg (gpMap _ (M.bmon.map (W.unop.hom.hom.base))) hu
  rw [gpMap_bmon_inv_left, gpMap_bmon_inv_left] at huu
  have hu' : φ.u = ψ.u := by
    obtain ⟨c, hc⟩ := toGp_eq_iff.mp (sub_left_injective huu)
    exact groupLike_add_right_cancel (h.bmonGroupLike _) hc
  exact congrArg (HomBirat.mk W) (Hom.ext hbase hdiv hd hu')

/-! ## ★4. `u` の合成則 -/

/-- ★`biratBase` を `Obj M` の底の型で書いたもの。 -/
noncomputable def biratBaseC (h : Hyp M) {A B : Obj M}
    (f : HomBirat (modelPre h) (model_frobenioid h) A B) : A.base ⟶ B.base := biratBase f

theorem biratBaseC_mk (h : Hyp M) {A B : Obj M}
    (Z : IdxBirat (modelPre h) (model_frobenioid h) A) (φ : Z.unop.left.obj ⟶ B)
    (hz : IsIso (Z.unop.hom.hom.base)) :
    biratBaseC h (HomBirat.mk Z φ) = inv (Z.unop.hom.hom.base) ≫ φ.base := by
  show biratBase (HomBirat.mk Z φ) = _
  rw [biratBase_mk, sliceBaseOf_eq]
  rfl

/-- ★★★★**代表元レベルの合成則**(`u` 版)。 -/
theorem sliceUGpOf_comp {A B E A' B' Dd : Obj M}
    (a : A' ⟶ A) (ha : IsIso a.base) (has : a.deg = 1)
    (b : B' ⟶ B) (hb : IsIso b.base) (hbs : b.deg = 1)
    (φ : A' ⟶ B) (ψ : B' ⟶ E)
    (γ : Dd ⟶ A') (hγ : IsIso γ.base) (hγs : γ.deg = 1)
    (α : Dd ⟶ B') (hsq : γ ≫ φ = α ≫ b)
    (hca : IsIso (γ ≫ a).base) :
    sliceUGpOf (γ ≫ a) hca (α ≫ ψ)
      = gpMap _ (M.bmon.map (inv a.base ≫ φ.base)) (sliceUGpOf b hb ψ)
        + ((ψ.deg : ℕ+) : ℕ) • sliceUGpOf a ha φ := by
  haveI := ha; haveI := hb; haveI := hγ; haveI := hca
  have hdα : α.deg = φ.deg := by
    have h : φ.deg * γ.deg = b.deg * α.deg := congrArg Hom.deg hsq
    rw [hγs, hbs, mul_one, one_mul] at h
    exact h.symm
  have hbsq : γ.base ≫ φ.base = α.base ≫ b.base := congrArg Hom.base hsq
  have hb1 : inv γ.base ≫ α.base = φ.base ≫ inv b.base := by
    rw [IsIso.inv_comp_eq, ← Category.assoc, hbsq, Category.assoc, IsIso.hom_inv_id,
      Category.comp_id]
  have hinvca : inv ((γ ≫ a).base) = inv a.base ≫ inv γ.base := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    show (γ.base ≫ a.base) ≫ (inv a.base ≫ inv γ.base) = 𝟙 _
    rw [Category.assoc, ← Category.assoc a.base, IsIso.hom_inv_id,
      Category.id_comp, IsIso.hom_inv_id]
  have hu : M.bmon.map γ.base φ.u + ((φ.deg : ℕ+) : ℕ) • γ.u
      = M.bmon.map α.base b.u + α.u := by
    have h : M.bmon.map γ.base φ.u + ((φ.deg : ℕ+) : ℕ) • γ.u
        = M.bmon.map α.base b.u + ((b.deg : ℕ+) : ℕ) • α.u := congrArg Hom.u hsq
    rwa [hbs, PNat.one_coe, one_smul] at h
  have hB : toGp _ α.u
      = toGp _ (M.bmon.map γ.base φ.u) + ((φ.deg : ℕ+) : ℕ) • toGp _ γ.u
        - toGp _ (M.bmon.map α.base b.u) := by
    have hgp : toGp (M.bmon.val Dd.base) (M.bmon.map γ.base φ.u)
          + ((φ.deg : ℕ+) : ℕ) • toGp _ γ.u
        = toGp _ (M.bmon.map α.base b.u) + toGp _ α.u := by
      rw [← toGp_nsmul, ← toGp_add, ← toGp_add, hu]
    exact eq_sub_of_add_eq' hgp.symm
  have hkey : toGp _ (α ≫ ψ).u - (((α ≫ ψ).deg : ℕ+) : ℕ) • toGp _ (γ ≫ a).u
      = gpMap _ (M.bmon.map α.base) (toGp _ ψ.u - ((ψ.deg : ℕ+) : ℕ) • toGp _ b.u)
        + ((ψ.deg : ℕ+) : ℕ) • gpMap _ (M.bmon.map γ.base)
          (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u) := by
    show toGp _ (M.bmon.map α.base ψ.u + ((ψ.deg : ℕ+) : ℕ) • α.u)
        - ((ψ.deg * α.deg : ℕ+) : ℕ)
          • toGp _ (M.bmon.map γ.base a.u + ((a.deg : ℕ+) : ℕ) • γ.u) = _
    rw [hdα, has, PNat.one_coe, one_smul,
      toGp_add, toGp_add, toGp_nsmul, map_sub, map_sub, gpMap_toGp, gpMap_toGp,
      map_nsmul, map_nsmul, gpMap_toGp, gpMap_toGp, hB, PNat.mul_coe]
    simp only [smul_sub, smul_add, smul_smul]
    abel
  have e1 : (inv a.base ≫ inv γ.base) ≫ α.base
      = (inv a.base ≫ φ.base) ≫ inv b.base := by
    rw [Category.assoc, hb1, ← Category.assoc]
  have e2 : (inv a.base ≫ inv γ.base) ≫ γ.base = inv a.base := by
    rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  rw [sliceUGpOf_eq, sliceUGpOf_eq, sliceUGpOf_eq, hkey, hinvca,
    map_add, map_nsmul, gpMap_bmon_comp, gpMap_bmon_comp, gpMap_bmon_comp, e1, e2]

/-- ★★★★**`Hom^birat` の `u` 成分は合成を保つ**。 -/
theorem biratU_comp (h : Hyp M) {A B E : Obj M}
    (f : HomBirat (modelPre h) (model_frobenioid h) A B)
    (g : HomBirat (modelPre h) (model_frobenioid h) B E) :
    biratU h (compBirat (modelPre h) (model_frobenioid h) (model_frobenioid h).core f g)
      = gpMap _ (M.bmon.map (biratBaseC h f)) (biratU h g)
        + ((biratDeg g : ℕ+) : ℕ) • biratU h f := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  haveI hz : IsIso (Z.unop.hom.hom.base) := Z.unop.hom.property.2.2
  rw [compBirat_mk, biratU_mk, biratU_mk, biratU_mk, biratBaseC_mk h Z φ hz, biratDeg_mk]
  exact sliceUGpOf_comp Z.unop.hom.hom Z.unop.hom.property.2.2 Z.unop.hom.property.2.1
    W.unop.hom.hom W.unop.hom.property.2.2 W.unop.hom.property.2.1 φ ψ
    (biratPullGamma (model_frobenioid h).core Z φ W)
    (biratPullGamma_preStep (model_frobenioid h).core Z φ W).2
    (biratPullGamma_preStep (model_frobenioid h).core Z φ W).1
    (biratPullAlpha (model_frobenioid h).core Z φ W)
    (biratPull_sq (model_frobenioid h).core Z φ W) _

theorem biratU_comp' (h : Hyp M) {A B E : BiratCat (modelPre h) (model_frobenioid h)}
    (f : A ⟶ B) (g : B ⟶ E) :
    biratU h (f ≫ g) = gpMap _ (M.bmon.map (biratBaseC h f)) (biratU h g)
      + ((biratDeg g : ℕ+) : ℕ) • biratU h f :=
  biratU_comp h f g

theorem biratU_toHomBirat (h : Hyp M) {A B : Obj M} (φ : A ⟶ B) :
    biratU h (toHomBirat (P := modelPre h) (G := model_frobenioid h) φ)
      = toGp _ φ.u := by
  rw [toHomBirat, biratU_mk]
  haveI hI : IsIso (Hom.base (𝟙 A)) := by
    rw [id_base]; infer_instance
  show sliceUGpOf (𝟙 A) hI φ = _
  rw [sliceUGpOf_eq]
  have hb : inv (Hom.base (𝟙 A)) = 𝟙 A.base := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [id_base, Category.id_comp]
  rw [show Hom.u (𝟙 A) = 0 from id_u A, hb, toGp_zero, smul_zero, sub_zero, gpMap_bmon_id]

theorem biratU_id (h : Hyp M) (A : BiratCat (modelPre h) (model_frobenioid h)) :
    biratU h (𝟙 A) = 0 := by
  show biratU h (toHomBirat (P := modelPre h) (G := model_frobenioid h)
    (𝟙 (biratDown (modelPre h) (model_frobenioid h) A))) = 0
  rw [biratU_toHomBirat,
    show Hom.u (𝟙 (biratDown (modelPre h) (model_frobenioid h) A)) = 0 from
      id_u (biratDown (modelPre h) (model_frobenioid h) A)]
  exact toGp_zero _

theorem gpMap_bmon_biratBaseC_of_baseIdentity (h : Hyp M)
    {A : BiratCat (modelPre h) (model_frobenioid h)} {δ : A ⟶ A}
    (hδ : IsBaseIdentity (biratPre (modelPre h) (model_frobenioid h)) δ)
    (x : Gp (M.bmon.val (biratDown (modelPre h) (model_frobenioid h) A).base)) :
    gpMap _ (M.bmon.map (biratBaseC h δ)) x = x := by
  have hb : biratBaseC h δ = 𝟙 _ := by
    show biratBase δ = 𝟙 _
    refine hδ.trans ?_
    show (biratPre (modelPre h) (model_frobenioid h)).Base (𝟙 A) = 𝟙 _
    rw [(biratPre (modelPre h) (model_frobenioid h)).Base_id]
  rw [hb, gpMap_bmon_id]

theorem biratU_pow_otri (h : Hyp M)
    {X : BiratCat (modelPre h) (model_frobenioid h)} {α : End X}
    (hα : α ∈ OTri (biratPre (modelPre h) (model_frobenioid h)) X) (n : ℕ) :
    biratU h (((α ^ n : End X) : X ⟶ X)) = n • biratU h ((α : X ⟶ X)) := by
  induction n with
  | zero => simpa using biratU_id h X
  | succ k ih =>
    have hmem : (α ^ k : End X) ∈ OTri (biratPre (modelPre h) (model_frobenioid h)) X :=
      Submonoid.pow_mem _ hα k
    have hcomp : ((α ^ (k + 1) : End X) : X ⟶ X)
        = ((α : X ⟶ X)) ≫ ((α ^ k : End X) : X ⟶ X) := by
      rw [pow_succ]; rfl
    rw [hcomp, biratU_comp' h, gpMap_bmon_biratBaseC_of_baseIdentity h hα.1,
      show ((biratDeg ((α ^ k : End X) : X ⟶ X) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg ((α ^ k : End X) : X ⟶ X) = 1 from hmem.2]; rfl,
      one_smul, ih, succ_nsmul]

/-! ## ★5. `Theorem 5.2, (ii)` の残り 1 —— birationally Frobenius-normalized 型 -/

/-- ★★★★★**model Frobenioid の `𝒞^birat` は Frobenius-normalized**。

★`birat_frobNormalized_of_unitTrivial`(`Thm51Sec.lean`)と**同じ 4 行の計算**である
——忠実性が `u` 込みの 4 不変量になっただけ。 -/
theorem model_birat_frobNormalized (h : Hyp M)
    (X : BiratCat (modelPre h) (model_frobenioid h)) :
    IsFrobeniusNormalized (biratPre (modelPre h) (model_frobenioid h)) X := by
  intro φ hφb α hα
  set Pb := biratPre (modelPre h) (model_frobenioid h) with hPb
  set d : ℕ := ((Pb.degFr ((φ : End X) : X ⟶ X) : ℕ+) : ℕ) with hddef
  have hβmem : (α ^ d : End X) ∈ OTri Pb X := Submonoid.pow_mem _ hα d
  have hbφ : Pb.Base ((φ : End X) : X ⟶ X) = 𝟙 _ := hφb.trans (Pb.Base_id X)
  have hbα : Pb.Base ((α : End X) : X ⟶ X) = 𝟙 _ := hα.1.trans (Pb.Base_id X)
  have hbβ : Pb.Base (((α ^ d : End X)) : X ⟶ X) = 𝟙 _ := hβmem.1.trans (Pb.Base_id X)
  refine model_birat_faithful h _ _ ?_ ?_ ?_ ?_
  · show Pb.Base _ = Pb.Base _
    rw [Pb.Base_comp, Pb.Base_comp, hbφ, hbα, hbβ, Category.id_comp]
  · show Pb.degFr _ = Pb.degFr _
    rw [Pb.degFr_comp, Pb.degFr_comp,
      show Pb.degFr ((α : End X) : X ⟶ X) = 1 from hα.2,
      show Pb.degFr (((α ^ d : End X)) : X ⟶ X) = 1 from hβmem.2,
      one_mul, mul_one]
  · rw [biratDivGp_comp', biratDivGp_comp', gpMap_biratBase_of_baseIdentity hφb,
      gpMap_biratBase_of_baseIdentity hα.1,
      show ((biratDeg (((α ^ d : End X)) : X ⟶ X) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg (((α ^ d : End X)) : X ⟶ X) = 1 from hβmem.2]; rfl,
      one_smul, biratDivGp_pow_otri (model_frobenioid h) hα d,
      show ((biratDeg ((φ : End X) : X ⟶ X) : ℕ+) : ℕ) = d from rfl]
    abel
  · rw [biratU_comp' h, biratU_comp' h, gpMap_bmon_biratBaseC_of_baseIdentity h hφb,
      gpMap_bmon_biratBaseC_of_baseIdentity h hα.1,
      show ((biratDeg (((α ^ d : End X)) : X ⟶ X) : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg (((α ^ d : End X)) : X ⟶ X) = 1 from hβmem.2]; rfl,
      one_smul, biratU_pow_otri h hα d,
      show ((biratDeg ((φ : End X) : X ⟶ X) : ℕ+) : ℕ) = d from rfl]
    abel

/-- ★★★★★**[FrdI] Theorem 5.2, (ii)** —— model Frobenioid は
**birationally Frobenius-normalized 型**。 -/
theorem model_isOfBiratFrobNormalizedType (h : Hyp M) :
    IsOfBirationallyFrobeniusNormalizedType (Obj M) (modelPre h) (model_frobenioid h) :=
  fun A => model_birat_frobNormalized h ((toBiratCat (modelPre h) (model_frobenioid h)).obj A)

/-- ★★★★★**[FrdI] Theorem 5.2, (ii)** —— model Frobenioid は **model 型**。 -/
theorem model_isOfModelType (h : Hyp M) (hskel : Skeletal D) :
    haveI := h.connectedD
    IsOfModelType (Obj M) (modelPre h) (model_frobenioid h) :=
  haveI := h.connectedD
  ⟨model_isPreModelType h hskel, model_isOfBiratFrobNormalizedType h⟩

/-! ## ★6. `Theorem 5.2, (ii)` の残り 2 —— `𝒪^×(−^birat) ≅ B` -/

/-- ★`Div_B` を群化の上へ延ばしたもの。 -/
noncomputable def divBGp (h : Hyp M) (A : D) :
    Gp (M.bmon.val A) →+ Gp (M.phi.val A) :=
  (M.divB A).comp (gpEquivOfGroupLike (h.bmonGroupLike A)).toAddMonoidHom

theorem divBGp_toGp (h : Hyp M) (A : D) (x : M.bmon.val A) :
    divBGp h A (toGp _ x) = M.divB A x := by
  show M.divB A (gpEquivOfGroupLike (h.bmonGroupLike A) (toGp _ x)) = _
  rw [gpEquivOfGroupLike_toGp]

theorem divBGp_nat (h : Hyp M) {A B : D} (θ : A ⟶ B) (z : Gp (M.bmon.val B)) :
    divBGp h A (gpMap _ (M.bmon.map θ) z) = gpMap _ (M.phi.map θ) (divBGp h B z) := by
  obtain ⟨y, rfl⟩ := toGp_surjective_of_groupLike (h.bmonGroupLike B) z
  rw [gpMap_toGp, divBGp_toGp, divBGp_toGp]
  exact M.divB_nat θ y

/-- ★★同じ底・次数 1 の 2 射の `Div` の差は `Div_B` の差。 -/
theorem div_sub_eq_divB_sub {X E : Obj M} (a : E ⟶ X) (had : a.deg = 1)
    (φ : E ⟶ X) (hfd : φ.deg = 1) (hfb : φ.base = a.base) :
    toGp _ φ.div - toGp _ a.div = M.divB _ φ.u - M.divB _ a.u := by
  have h1 := φ.cond
  have h2 := a.cond
  rw [hfd, hfb, PNat.one_coe, one_smul] at h1
  rw [had, PNat.one_coe, one_smul] at h2
  have h3 : (E.cls + toGpHom (M.phi.val E.base) φ.div)
      - (E.cls + toGpHom (M.phi.val E.base) a.div)
      = (M.phi.gpMapOn a.base X.cls + M.divB _ φ.u)
        - (M.phi.gpMapOn a.base X.cls + M.divB _ a.u) := by
    rw [h1, h2]
  rwa [add_sub_add_left_eq_sub, add_sub_add_left_eq_sub] at h3

/-- ★★★**`Div_B ∘ κ = Div^gp`(代表元レベル)**。 -/
theorem divBGp_sliceUGpOf (h : Hyp M) {X E : Obj M} (a : E ⟶ X) (ha : IsIso a.base)
    (had : a.deg = 1) (φ : E ⟶ X) (hfd : φ.deg = 1) (hfb : φ.base = a.base) :
    divBGp h X.base (sliceUGpOf a ha φ) = sliceDivGpOf (P := modelPre h) a ha φ := by
  haveI := ha
  rw [sliceUGpOf_eq, divBGp_nat]
  show gpMap _ (M.phi.map (inv a.base))
      (divBGp h E.base (toGp _ φ.u - ((φ.deg : ℕ+) : ℕ) • toGp _ a.u)) = _
  rw [hfd, PNat.one_coe, one_smul, map_sub, divBGp_toGp, divBGp_toGp,
    ← div_sub_eq_divB_sub a had φ hfd hfb]
  show _ = gpMap _ (M.phi.map (inv a.base))
      (toGp _ φ.div - ((φ.deg : ℕ+) : ℕ) • toGp _ a.div)
  rw [hfd, PNat.one_coe, one_smul]

/-- ★★★★**`Div_B ∘ κ = Div^gp`**(`𝒪^▷(X^birat)` の上で)。 -/
theorem divBGp_biratU_otri (h : Hyp M) (X : BiratCat (modelPre h) (model_frobenioid h))
    (δ : OTri (biratPre (modelPre h) (model_frobenioid h)) X) :
    divBGp h _ (biratU h ((δ : End X) : X ⟶ X)) = biratDivGp ((δ : End X) : X ⟶ X) := by
  obtain ⟨W, φ, hW⟩ := HomBirat.exists_rep ((δ : End X) : X ⟶ X)
  obtain ⟨hfd, hfb⟩ := otri_rep_base_deg X δ W φ hW
  have hfd' : φ.deg = 1 := hfd
  have hfb' : φ.base = W.unop.hom.hom.base := hfb
  have had : W.unop.hom.hom.deg = 1 := W.unop.hom.property.2.1
  haveI hia : IsIso (W.unop.hom.hom.base) := W.unop.hom.property.2.2
  rw [← hW, biratU_mk, biratDivGp_mk]
  exact divBGp_sliceUGpOf h W.unop.hom.hom hia had φ hfd' hfb'

/-- ★★`𝒪^×(X^birat) → B(Base X)`。 -/
noncomputable def modelKappaHom (h : Hyp M)
    (X : BiratCat (modelPre h) (model_frobenioid h)) :
    Additive ↥(OTimes (biratPre (modelPre h) (model_frobenioid h)) X)
      →+ M.bmon.val (biratDown (modelPre h) (model_frobenioid h) X).base where
  toFun δ := gpEquivOfGroupLike (h.bmonGroupLike _)
    (biratU h (((Additive.toMul δ :
      OTimes (biratPre (modelPre h) (model_frobenioid h)) X) : End X) : X ⟶ X))
  map_zero' := by
    show gpEquivOfGroupLike (h.bmonGroupLike _) (biratU h (𝟙 X)) = 0
    rw [biratU_id]
    exact map_zero _
  map_add' x y := by
    show gpEquivOfGroupLike (h.bmonGroupLike _)
        (biratU h ((((Additive.toMul y :
              OTimes (biratPre (modelPre h) (model_frobenioid h)) X) : End X) : X ⟶ X)
          ≫ (((Additive.toMul x :
              OTimes (biratPre (modelPre h) (model_frobenioid h)) X) : End X) : X ⟶ X))) = _
    rw [biratU_comp' h,
      gpMap_bmon_biratBaseC_of_baseIdentity h
        (Additive.toMul y : OTimes (biratPre (modelPre h) (model_frobenioid h)) X).2.1.1,
      show ((biratDeg ((((Additive.toMul x :
            OTimes (biratPre (modelPre h) (model_frobenioid h)) X) : End X) : X ⟶ X))
          : ℕ+) : ℕ) = 1 from by
        rw [show biratDeg ((((Additive.toMul x :
            OTimes (biratPre (modelPre h) (model_frobenioid h)) X) : End X) : X ⟶ X)) = 1 from
          (Additive.toMul x :
            OTimes (biratPre (modelPre h) (model_frobenioid h)) X).2.1.2]; rfl,
      one_smul, map_add]

theorem modelKappaHom_apply (h : Hyp M) (X : BiratCat (modelPre h) (model_frobenioid h))
    (δ : OTimes (biratPre (modelPre h) (model_frobenioid h)) X) :
    modelKappaHom h X (Additive.ofMul δ)
      = gpEquivOfGroupLike (h.bmonGroupLike _) (biratU h (((δ : End X) : X ⟶ X))) := rfl

/-- ★`[a]⁻¹ ≫ [f]`(`a`・`f` とも co-angular pre-step)は `𝒞^birat` の同型。 -/
theorem birat_mk_isIso_model (h : Hyp M) {X E : Obj M} (a : E ⟶ X) (f : E ⟶ X)
    (hac : IsCoAngular (modelPre h) a) (has : IsPreStep (modelPre h) a)
    (hfc : IsCoAngular (modelPre h) f) (hfs : IsPreStep (modelPre h) f) :
    @IsIso (BiratCat (modelPre h) (model_frobenioid h)) _ X X
      (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h) a hac has) f) := by
  haveI h1 : IsIso ((toBiratCat (modelPre h) (model_frobenioid h)).map a) :=
    birat_isIso_of_coaPre a hac has
  haveI h2 : IsIso ((toBiratCat (modelPre h) (model_frobenioid h)).map f) :=
    birat_isIso_of_coaPre f hfc hfs
  have hcomp : (toBiratCat (modelPre h) (model_frobenioid h)).map a
      ≫ (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h) a hac has) f)
      = (toBiratCat (modelPre h) (model_frobenioid h)).map f := birat_toHom_comp_mk a hac has f
  have h3 : IsIso ((toBiratCat (modelPre h) (model_frobenioid h)).map a
      ≫ (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h) a hac has) f)) := by
    rw [hcomp]; exact h2
  exact @IsIso.of_isIso_comp_left (BiratCat (modelPre h) (model_frobenioid h)) _ E X X
    ((toBiratCat (modelPre h) (model_frobenioid h)).map a)
    (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h) a hac has) f) h1 h3

/-- ★与えられた `u` を実現する射(`Div_B u = toGp d₁ - toGp d₂` の分解から)。 -/
def uHom (X : Obj M) (d₁ d₂ : M.phi.val X.base) (y : M.bmon.val X.base)
    (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂) : objDown X d₂ ⟶ X :=
  ⟨𝟙 X.base, d₁, 1, y, by
    show ((1 : ℕ+) : ℕ) • (X.cls - toGp (M.phi.val X.base) d₂) + toGp (M.phi.val X.base) d₁
      = M.phi.gpMapOn (𝟙 X.base) X.cls + M.divB _ y
    rw [MonoidOn.gpMapOn_id, hy, PNat.one_coe, one_smul]
    abel⟩

@[simp] theorem uHom_base (X : Obj M) (d₁ d₂ : M.phi.val X.base)
    (y : M.bmon.val X.base) (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂) :
    (uHom X d₁ d₂ y hy).base = 𝟙 X.base := rfl
@[simp] theorem uHom_deg (X : Obj M) (d₁ d₂ : M.phi.val X.base)
    (y : M.bmon.val X.base) (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂) :
    (uHom X d₁ d₂ y hy).deg = 1 := rfl
@[simp] theorem uHom_u (X : Obj M) (d₁ d₂ : M.phi.val X.base)
    (y : M.bmon.val X.base) (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂) :
    (uHom X d₁ d₂ y hy).u = y := rfl

theorem uHom_isPreStep (h : Hyp M) (X : Obj M) (d₁ d₂ : M.phi.val X.base)
    (y : M.bmon.val X.base) (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂) :
    IsPreStep (modelPre h) (uHom X d₁ d₂ y hy) := by
  refine ⟨rfl, ?_⟩
  show IsIso (𝟙 X.base)
  infer_instance

theorem sliceUGpOf_preStepDown (X0 : Obj M) (d₁ d₂ : M.phi.val X0.base)
    (y : M.bmon.val X0.base) (hy : M.divB _ y = toGp _ d₁ - toGp _ d₂)
    (hab : IsIso ((preStepDown X0 d₂).base)) :
    sliceUGpOf (preStepDown X0 d₂) hab (uHom X0 d₁ d₂ y hy) = toGp _ y := by
  haveI := hab
  have hinv : inv ((preStepDown X0 d₂).base) = 𝟙 X0.base := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    show (𝟙 X0.base) ≫ 𝟙 X0.base = 𝟙 X0.base
    rw [Category.id_comp]
  rw [sliceUGpOf_eq]
  show gpMap _ (M.bmon.map (inv ((preStepDown X0 d₂).base)))
      (toGp (M.bmon.val (objDown X0 d₂).base) y
        - (((1 : ℕ+)) : ℕ) • toGp (M.bmon.val (objDown X0 d₂).base) 0) = toGp _ y
  rw [PNat.one_coe, one_smul, toGp_zero, sub_zero, hinv]
  exact gpMap_bmon_id _

theorem modelKappaHom_surjective (h : Hyp M) (X0 : Obj M) :
    Function.Surjective
      (modelKappaHom h ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0)) := by
  intro y
  obtain ⟨d₁, d₂, hy⟩ := exists_toGp_sub (M.divB X0.base y)
  have hac : IsCoAngular (modelPre h) (preStepDown X0 d₂) := model_coAngular h _
  have has : IsPreStep (modelPre h) (preStepDown X0 d₂) := preStepDown_isPreStep h X0 d₂
  have hfc : IsCoAngular (modelPre h) (uHom X0 d₁ d₂ y hy) := model_coAngular h _
  have hfs : IsPreStep (modelPre h) (uHom X0 d₁ d₂ y hy) := uHom_isPreStep h X0 d₁ d₂ y hy
  have hZh : (idxBiratMk (modelPre h) (model_frobenioid h)
      (preStepDown X0 d₂) hac has).unop.hom.hom = preStepDown X0 d₂ :=
    idxBiratMk_hom (modelPre h) (model_frobenioid h) (preStepDown X0 d₂) hac has
  haveI hab : IsIso ((preStepDown X0 d₂).base) := by
    show IsIso (𝟙 X0.base); infer_instance
  haveI hZb : IsIso ((idxBiratMk (modelPre h) (model_frobenioid h)
      (preStepDown X0 d₂) hac has).unop.hom.hom.base) := by rw [hZh]; exact hab
  have hinv : inv ((preStepDown X0 d₂).base) = 𝟙 X0.base := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    show (𝟙 X0.base) ≫ 𝟙 X0.base = 𝟙 X0.base
    rw [Category.id_comp]
  have hdeg : biratDeg (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h)
      (preStepDown X0 d₂) hac has) (uHom X0 d₁ d₂ y hy)) = 1 := by
    rw [biratDeg_mk]; rfl
  have hbase : biratBase (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h)
      (preStepDown X0 d₂) hac has) (uHom X0 d₁ d₂ y hy)) = 𝟙 _ := by
    rw [biratBase_mk, sliceBaseOf_eq]
    show inv ((preStepDown X0 d₂).base) ≫ (uHom X0 d₁ d₂ y hy).base = 𝟙 _
    rw [hinv]
    exact Category.id_comp _
  have hmem : (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h)
        (preStepDown X0 d₂) hac has) (uHom X0 d₁ d₂ y hy))
      ∈ OTimes (biratPre (modelPre h) (model_frobenioid h))
        ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0) := by
    refine ⟨⟨?_, hdeg⟩, ?_⟩
    · show (biratPre (modelPre h) (model_frobenioid h)).Base _
        = (biratPre (modelPre h) (model_frobenioid h)).Base
          (𝟙 ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0))
      rw [(biratPre (modelPre h) (model_frobenioid h)).Base_id]
      exact hbase
    · exact (isUnit_iff_isIso _).mpr (birat_mk_isIso_model h (preStepDown X0 d₂)
        (uHom X0 d₁ d₂ y hy) hac has hfc hfs)
  refine ⟨Additive.ofMul ⟨_, hmem⟩, ?_⟩
  rw [modelKappaHom_apply]
  show gpEquivOfGroupLike (h.bmonGroupLike _)
    (biratU h (HomBirat.mk (idxBiratMk (modelPre h) (model_frobenioid h)
      (preStepDown X0 d₂) hac has) (uHom X0 d₁ d₂ y hy))) = y
  rw [biratU_mk, sliceUGpOf_congr hZh hZb hab (uHom X0 d₁ d₂ y hy)]
  exact (congrArg (gpEquivOfGroupLike (h.bmonGroupLike X0.base))
    (sliceUGpOf_preStepDown X0 d₁ d₂ y hy hab)).trans
    (gpEquivOfGroupLike_toGp (h.bmonGroupLike X0.base) y)

theorem modelKappaHom_injective (h : Hyp M)
    (X : BiratCat (modelPre h) (model_frobenioid h)) :
    Function.Injective (modelKappaHom h X) := by
  intro x z hxz
  set dx : OTimes (biratPre (modelPre h) (model_frobenioid h)) X := Additive.toMul x with hdx
  set dz : OTimes (biratPre (modelPre h) (model_frobenioid h)) X := Additive.toMul z with hdz
  have hu : biratU h ((dx : End X) : X ⟶ X) = biratU h ((dz : End X) : X ⟶ X) :=
    (gpEquivOfGroupLike (h.bmonGroupLike _)).injective hxz
  have hzz : biratDivGp ((dx : End X) : X ⟶ X) = biratDivGp ((dz : End X) : X ⟶ X) := by
    rw [← divBGp_biratU_otri h X ⟨(dx : End X), OTimes_le_OTri _ _ dx.2⟩,
      ← divBGp_biratU_otri h X ⟨(dz : End X), OTimes_le_OTri _ _ dz.2⟩, hu]
  have hb : biratBase ((dx : End X) : X ⟶ X) = biratBase ((dz : End X) : X ⟶ X) :=
    (dx.2.1.1).trans (dz.2.1.1).symm
  have hd : biratDeg ((dx : End X) : X ⟶ X) = biratDeg ((dz : End X) : X ⟶ X) :=
    (dx.2.1.2).trans (dz.2.1.2).symm
  have heq := model_birat_faithful h _ _ hb hd hzz hu
  exact Additive.toMul.injective (Subtype.ext heq)

/-- ★★★★★**`𝒪^×(X^birat) ≅ B(Base X)`**(`Theorem 5.2, (ii)` の「Moreover」)。 -/
noncomputable def modelKappa (h : Hyp M) (X0 : Obj M) :
    Additive ↥(OTimes (biratPre (modelPre h) (model_frobenioid h))
        ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0))
      ≃+ M.bmon.val ((modelPre h).toElem.obj X0).base :=
  AddEquiv.ofBijective (modelKappaHom h _)
    ⟨modelKappaHom_injective h _, modelKappaHom_surjective h X0⟩

/-- ★★`Div_B ∘ κ = Div^gp`。 -/
theorem modelKappa_divB (h : Hyp M) (X0 : Obj M)
    (δ : OTimes (biratPre (modelPre h) (model_frobenioid h))
      ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0)) :
    M.divB _ (modelKappa h X0 (Additive.ofMul δ))
      = biratDivGp ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0))
        : _ ⟶ _) := by
  show divBGp h _ (biratU h ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0))
    : _ ⟶ _)) = _
  exact divBGp_biratU_otri h ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0)
    ⟨(δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj X0)),
      OTimes_le_OTri (biratPre (modelPre h) (model_frobenioid h)) _ δ.2⟩

/-- ★★★★**`κ` の pull-back 射に沿った自然性**。 -/
theorem modelKappa_pull (h : Hyp M) {A B : Obj M}
    (φ : (toBiratCat (modelPre h) (model_frobenioid h)).obj A
      ⟶ (toBiratCat (modelPre h) (model_frobenioid h)).obj B)
    (θ : ((modelPre h).toElem.obj A).base ⟶ ((modelPre h).toElem.obj B).base)
    (δ : OTimes (biratPre (modelPre h) (model_frobenioid h))
      ((toBiratCat (modelPre h) (model_frobenioid h)).obj B))
    (ε : OTimes (biratPre (modelPre h) (model_frobenioid h))
      ((toBiratCat (modelPre h) (model_frobenioid h)).obj A))
    (hpb : IsPullBack (biratPre (modelPre h) (model_frobenioid h)) φ)
    (hθ : biratBase φ = θ)
    (hsq : φ ≫ ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj B)) : _ ⟶ _)
      = ((ε : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj A)) : _ ⟶ _) ≫ φ) :
    modelKappa h A (Additive.ofMul ε)
      = M.bmon.map θ (modelKappa h B (Additive.ofMul δ)) := by
  have hθ' : biratBaseC h φ = θ := hθ
  have hdφ : biratDeg φ = 1 := ((model_birat_frobenioidCore h).pullBackLB φ hpb).2
  have hkey := congrArg (biratU h) hsq
  rw [biratU_comp' h, biratU_comp' h,
    gpMap_bmon_biratBaseC_of_baseIdentity h ε.2.1.1,
    show ((biratDeg ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj B))
        : _ ⟶ _) : ℕ+) : ℕ) = 1 from by rw [show biratDeg _ = 1 from δ.2.1.2]; rfl,
    show ((biratDeg φ : ℕ+) : ℕ) = 1 from by rw [hdφ]; rfl,
    one_smul, one_smul] at hkey
  have hε : biratU h ((ε : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj A)) : _ ⟶ _)
      = gpMap _ (M.bmon.map (biratBaseC h φ))
        (biratU h ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj B))
          : _ ⟶ _)) := by
    have h2 : gpMap _ (M.bmon.map (biratBaseC h φ))
          (biratU h ((δ : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj B))
            : _ ⟶ _))
        + biratU h φ
        = biratU h ((ε : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj A)) : _ ⟶ _)
          + biratU h φ := by
      rw [hkey]; abel
    exact (add_right_cancel h2).symm
  show gpEquivOfGroupLike (h.bmonGroupLike _)
      (biratU h ((ε : End ((toBiratCat (modelPre h) (model_frobenioid h)).obj A)) : _ ⟶ _)) = _
  rw [hε, gpEquivOfGroupLike_map (h.bmonGroupLike _) (h.bmonGroupLike _), hθ']
  rfl

/-- ★★★★★★**[FrdI] Theorem 5.2, (ii) の「Moreover」** ——
model Frobenioid では `RatFnData`(`𝒪^×(−^birat) ≅ B`、`Div_B` と両立、
pull-back 射に沿った自然性)が**構成できる**。

★これで `Theorem 5.2, (iv)` の interface が空でないことも同時に分かる。 -/
noncomputable def modelRatFnData (h : Hyp M) :
    RatFnData (modelPre h) (model_frobenioid h) where
  bmon := M.bmon
  divB := M.divB
  divB_nat := M.divB_nat
  kappa := modelKappa h
  kappa_divB := modelKappa_divB h
  kappa_pull φ θ δ ε hpb hθ hsq := modelKappa_pull h φ θ δ ε hpb hθ hsq
  bneg A x := negOfGroupLike (h.bmonGroupLike A) x
  bneg_add A x := negOfGroupLike_add (h.bmonGroupLike A) x

/-- ★★★★★**[FrdI] Theorem 5.2, (iii) の後半** —— **rationally standard 型の判定**。

原文 (FrdI p.101):
> (iii) C is of standard type if and only if the following conditions are satisﬁed: (a) if

★★原文の (a) のうち **birationally Frobenius-normalized 型**は model Frobenioid では
**自動である**(`model_isOfBiratFrobNormalizedType`)。それがこの同値の中身である。 -/
theorem model_ratStandardType_iff (h : Hyp M)
    (ι : ∀ Y : D, Prime (M.phi.val Y) → Pf (M.phi.val Y) → NNReal) :
    IsOfRationallyStandardType (modelPre h) (model_frobenioid h) ι ↔
      (IsOfRationalType (Obj M) (modelPre h) (model_frobenioid h) ι
        ∧ IsOfStandardType D (Obj M) (modelPre h) (model_frobenioid h).core
        ∧ ∃ X : BiratCat (unTrPre (modelPre h) (model_frobenioid h).core)
            (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)),
          IsFrobeniusCompact
            (unTrBiratPre (modelPre h) (model_frobenioid h).core (model_frobenioid h)) X) := by
  constructor
  · intro hs
    exact ⟨hs.rational, hs.standard, hs.unTrBiratCompact⟩
  · rintro ⟨h1, h2, h3⟩
    exact
      { biratFrobNormalized := model_isOfBiratFrobNormalizedType h
        rational := h1
        standard := h2
        unTrBiratCompact := h3 }

end ModelData

/-! ### ★出典の紐付け(`.src`) -/

/-- ★★★locator —— `Theorem 5.2, (ii)` の残り 2 点(model 型 ＋ `𝒪^×(−^birat) ≅ B`)。 -/
def ModelData.model_isOfModelType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (ii) — model 型 と 𝒪^×(−^birat) ≅ B の自然同型",
    sectionId := "frdi-thm-5-2" }

/-- ★★★★★★★**[FrdI] Theorem 5.2** —— 4 条がすべて実装された。

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | 圏の構成と `F_Phi` への関手 | `ModelData.Obj` / `ModelData.toElem` |
| (ii) | Frobenioid・isotropic 型 | `ModelData.model_frobenioid` |
| (ii) | model 型 | `ModelData.model_isOfModelType` |
| (ii) | `O^x(-^birat) = B`(`Div_B` と両立) | `ModelData.modelRatFnData` |
| (iii) | standard 型の判定 | `ModelData.model_standardType_iff` |
| (iii) | rationally standard 型の判定 | `ModelData.model_ratStandardType_iff` |
| (iv) | model 型の Frobenioid は model Frobenioid と圏同値 | `thm_5_2_iv` | -/
def thm_5_2.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 100, item := "Theorem 5.2", sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
