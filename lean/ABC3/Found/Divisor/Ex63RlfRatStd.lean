/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Rlf
import ABC3.Found.Divisor.Ex63RatStd

/-!
# [FrdI] Theorem 6.4, (i) —— 算術の `𝒞^rlf` は **rationally standard 型**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> Then the Frobenioids C, Cpf, Crlf, Cun-tr, (Cpf)un-tr are of isotropic and
> rationally standard type, but not of group-like type

## ★★★★★★★これで 5 圏がそろった

`Ex63RatStd.lean` は `𝒞`・`𝒞^un-tr`・`𝒞^pf`・`(𝒞^pf)^un-tr` の 4 圏を出した。
★本ファイルが**残る `𝒞^rlf`** を出す。

★★中身は `Ex63RatStd.lean` の証明を `Γ_L = arithDivGroup L` から
**`Γ_L = ⊤`**(`Ex63Rlf.lean` の `arithDatumRlf`)に取り替えただけである ——
`gen` / `locMono` / `coord` が自明になるぶん、むしろ短い。

★★★**テンソル `ScT ℝ≥0` も錐 `rlfCone` も、`ℝ≥0` の ℕ 上の平坦性も使わない。**
逸脱(`Crlf(𝒞) ≅ 𝒞(Δ^rlf)` を証明していないこと)の記録は `Ex63Rlf.lean` にある。
-/

namespace ABC3.Found.Divisor

open ABC3.Found.FrdI CategoryTheory _root_.NumberField
open ABC3.Meta

/-! ## ★0. `⊤` についての自明な条件 -/

/-- ★`⊤` は各素点に正の元をもつ(自明)。 -/
theorem isGenSubgroupR_top {S : Type} : IsGenSubgroupR (⊤ : AddSubgroup (S →₀ ℝ)) :=
  fun _ => ⟨1, one_pos, AddSubgroup.mem_top _⟩

/-- ★`⊤` は座標ごとに閉じている(自明)。 -/
theorem isCoordwiseR_top {S : Type} : IsCoordwiseR (⊤ : AddSubgroup (S →₀ ℝ)) :=
  fun _ _ _ => AddSubgroup.mem_top _

section

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-! ## ★1. non-dilating -/

variable (F Kbar) in
/-- ★★★★★★**[FrdI] Theorem 6.4, (i)** —— **`Φ^rlf` は non-dilating**。

★証明は `arithDatumGalois_isNonDilatingOn` と同じ(`Δ.coord` と `pullOf_apply` しか
使っていないので、`grp` を替えても通る)。 -/
theorem arithDatumRlf_isNonDilatingOn :
    MonoidOn.IsNonDilatingOn ((arithDatumRlf F Kbar).phi finSubOp_isOfFSMType) := by
  intro A e
  set Δ := arithDatumRlf F Kbar with hΔ
  have hmap : (Δ.phi finSubOp_isOfFSMType).map e = Δ.mapHom e := rfl
  rw [hmap]
  refine isNonDilating_of_primary_sharp (isSharp_effR _) _
    (closure_primary_effR_eq_top (Δ.coord A)) ?_
  intro a ha hprec
  obtain ⟨v, hv⟩ := effR_single_supported_of_primary (Δ.coord A) ha
  set R : Δ.primes A → Δ.primes A :=
    fun V => @resPlace _ _ _ _ _ _ (algOfHom e.unop) V with hR
  have happ : ∀ (x : effR (Δ.grp A)) (V : Δ.primes A),
      ((Δ.mapHom e x : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) V
        = ((x : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) (R V) :=
    fun x V => pullOf_apply e.unop _ V
  obtain ⟨n, hn, c, hc⟩ := hprec
  have hcoe : ((Δ.mapHom e a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ)
      + ((c : effR (Δ.grp A)) : Δ.primes A →₀ ℝ)
      = n • ((a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) :=
    congrArg Subtype.val hc
  have hzero : ∀ V : Δ.primes A, V ≠ v →
      ((a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) (R V) = 0 := by
    intro V hV
    have h1 := congrArg (fun f : Δ.primes A →₀ ℝ => f V) hcoe
    simp only [Finsupp.add_apply, Finsupp.smul_apply] at h1
    rw [hv, Finsupp.single_eq_of_ne hV, smul_zero] at h1
    rw [← happ a V]
    have h2 := (mem_effR.mp (Δ.mapHom e a).2).2 V
    have h3 := (mem_effR.mp c.2).2 V
    linarith
  have hane : (Δ.mapHom e a) ≠ 0 := fun h =>
    ha.1 (Δ.mapHom_injective e (by rw [h, map_zero]))
  have hex : ∃ V : Δ.primes A, ((Δ.mapHom e a : effR (Δ.grp A)) : Δ.primes A →₀ ℝ) V ≠ 0 := by
    by_contra hcon
    simp only [not_exists, Classical.not_not] at hcon
    exact hane (Subtype.ext (Finsupp.ext hcon))
  obtain ⟨V₀, hV₀⟩ := hex
  have hV₀v : V₀ = v := by
    by_contra hne
    exact hV₀ ((happ a V₀).trans (hzero V₀ hne))
  rw [hV₀v] at hV₀
  have hRv : R v = v := by
    by_contra hne
    refine hV₀ ((happ a v).trans ?_)
    rw [hv]
    exact Finsupp.single_eq_of_ne hne
  refine Subtype.ext (Finsupp.ext fun V => ?_)
  rcases eq_or_ne V v with rfl | hV
  · rw [happ a V, hRv]
  · rw [happ a V, hzero V hV, hv]
    exact (Finsupp.single_eq_of_ne hV).symm

/-! ## ★2. 有理関数の除数まわり -/

/-- ★★**有効因子の差が `arithDiv y` なら、その差は `Div_B(y)` である**。 -/
theorem toGp_sub_eq_divBRlf (A : (FinSub F Kbar)ᵒᵖ) (y : (bmonGalois F Kbar).val A)
    (a b : effR (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ)))
    (hab : ((a : ArithPlace A.unop.toIF →₀ ℝ)) - ((b : ArithPlace A.unop.toIF →₀ ℝ))
      = arithDiv (Additive.toMul y)) :
    toGp _ a - toGp _ b = divBRlf A y := by
  refine effRGpHom_injective (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ)) ?_
  have h1 : effRGpHom (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ)) (toGp _ a - toGp _ b)
      = ⟨(a : ArithPlace A.unop.toIF →₀ ℝ) - (b : ArithPlace A.unop.toIF →₀ ℝ),
          AddSubgroup.mem_top _⟩ := by
    rw [effRGpHom, map_sub, gpLift_toGp, gpLift_toGp]
    rfl
  rw [h1]
  show _ = (effRGpEquiv (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ))
      ((arithDatumRlf F Kbar).gen A))
      ((effRGpEquiv (⊤ : AddSubgroup (ArithPlace A.unop.toIF →₀ ℝ))
        ((arithDatumRlf F Kbar).gen A)).symm
        (arithDivGroupHomTop A.unop.toIF y))
  rw [AddEquiv.apply_symm_apply]
  exact Subtype.ext hab

/-- ★★**`Div_B(y)` は `arithDiv y ≠ 0` なら無限位数**。 -/
theorem divBRlf_nsmul_ne_zero (A : (FinSub F Kbar)ᵒᵖ)
    (y : ((A.unop.toIF))ˣ) (v : ArithPlace A.unop.toIF) (hy : arithDiv y v ≠ 0)
    (n : ℕ) (hn : 0 < n) :
    n • (divBRlf A (Additive.ofMul y)) ≠ 0 := by
  intro hc
  have h1 := congrArg (phiGpHom (arithDatumRlf F Kbar) finSubOp_isOfFSMType) hc
  rw [map_nsmul, map_zero] at h1
  have h2 : Subtype.val (phiGpHom (arithDatumRlf F Kbar) finSubOp_isOfFSMType
      (divBRlf A (Additive.ofMul y))) = arithDiv y :=
    phiGpHom_divBRlf A (Additive.ofMul y)
  have h4 := congrArg (fun w : ((arithDatumRlf F Kbar).grp A) => (Subtype.val w) v) h1
  simp only [AddSubgroup.coe_nsmul, Finsupp.smul_apply, ZeroMemClass.coe_zero,
    Finsupp.coe_zero, nsmul_eq_mul, h2] at h4
  rcases mul_eq_zero.mp h4 with hz | hz
  · exact absurd (Nat.cast_eq_zero.mp hz) hn.ne'
  · exact hy hz

end

/-! ## ★3. `rationally standard` の 3 つの入力 -/

section Std

variable (F Kbar : Type) [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-- ★★★**`Φ^rlf` は `⊥` の上では恒等しか誘導しない**。 -/
theorem rlf_phi_map_bot_eq_id
    (e : (Opposite.op (botSub F Kbar)) ⟶ (Opposite.op (botSub F Kbar))) :
    (rlfModelData F Kbar).phi.map e
      = AddMonoidHom.id ((rlfModelData F Kbar).phi.val (Opposite.op (botSub F Kbar))) := by
  have he : e = 𝟙 _ := Quiver.Hom.unop_inj (botSub_end_eq_id e.unop)
  rw [he]
  ext x
  exact MonoidOn.map_id _ _ x

/-- ★★**因子分解写像 `ι`**(`Definition 2.4, (i)` の (c)(d))。 -/
noncomputable def rlfIota :
    ∀ Y : (FinSub F Kbar)ᵒᵖ, Prime ((rlfModelData F Kbar).phi.val Y)
      → Pf ((rlfModelData F Kbar).phi.val Y) → NNReal :=
  fun Y => iotaR (Γ := (⊤ : AddSubgroup (ArithPlace Y.unop.toIF →₀ ℝ))) isGenSubgroupR_top

/-- ★★★★**[FrdI] Definition 4.5, (ii) の入力**(`𝒞^rlf` が strictly rational)。 -/
theorem rlf_hsp :
    ∀ (A : ModelData.Obj (rlfModelData F Kbar))
      (p : Prime ((rlfModelData F Kbar).phi.val
        ((ModelData.modelPre (rlfHyp F Kbar)).toElem.obj A).base)),
      ∃ (a b : (rlfModelData F Kbar).phi.val
          ((ModelData.modelPre (rlfHyp F Kbar)).toElem.obj A).base)
        (y : (ModelData.modelRatFnData (rlfHyp F Kbar)).bmon.val
          ((ModelData.modelPre (rlfHyp F Kbar)).toElem.obj A).base),
        (toGp _ a - toGp _ b = (ModelData.modelRatFnData (rlfHyp F Kbar)).divB _ y ∨
          toGp _ a - toGp _ b = -((ModelData.modelRatFnData (rlfHyp F Kbar)).divB _ y)) ∧
        p ∈ SuppElt (rlfIota F Kbar _) a ∧ p ∉ SuppElt (rlfIota F Kbar _) b := by
  classical
  intro A p
  obtain ⟨y, hy⟩ := exists_arithDiv_ne_zero_at A.base.unop.toIF
    (effRPrimeEquiv (Γ := (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ)))
      isGenSubgroupR_top p)
  obtain ⟨a, b, z, hzmem, hzeq, hab, hap, hbp⟩ :=
    exists_split_suppElt_R (Γ := (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ)))
      isCoordwiseR_top isGenSubgroupR_top (AddSubgroup.mem_top (arithDiv y)) p hy
  refine ⟨a, b, Additive.ofMul y, ?_, hap, hbp⟩
  have hsub : ((a : effR (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ)))
        : ArithPlace _ →₀ ℝ)
      - ((b : effR (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ)))
        : ArithPlace _ →₀ ℝ) = z := by
    rw [hab]; abel
  rcases hzeq with rfl | rfl
  · exact Or.inl (toGp_sub_eq_divBRlf A.base (Additive.ofMul y) a b hsub)
  · have h2 := toGp_sub_eq_divBRlf A.base (Additive.ofMul y) b a
      (by rw [← neg_sub, hsub]; abel)
    have h3 : toGp (effR (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ))) a
        - toGp (effR (⊤ : AddSubgroup (ArithPlace A.base.unop.toIF →₀ ℝ))) b
        = -(divBRlf A.base (Additive.ofMul y)) := by
      rw [← h2]; abel
    exact Or.inr h3

/-- ★★★**[FrdI] Theorem 6.4, (i)** —— `𝒞^rlf` は **group-like 型でない**。

★中身は「`Φ^rlf(L) ≠ 0` かつ sharp」だけ(無限素点に台をもつ有効因子がある)。 -/
theorem rlf_not_isOfGroupLikeType :
    ¬ IsOfGroupLikeType (ModelData.modelPre (rlfHyp F Kbar)) := by
  classical
  intro hgl
  have hcc := hgl ⟨Opposite.op (botSub F Kbar), 0⟩
  obtain ⟨w⟩ := (inferInstance : Nonempty (InfinitePlace ((botSub F Kbar).toIF)))
  set v : ArithPlace ((botSub F Kbar).toIF) := Sum.inr w with hv
  set xx : effR (⊤ : AddSubgroup (ArithPlace ((botSub F Kbar).toIF) →₀ ℝ)) :=
    ⟨Finsupp.single v 1, mem_effR.mpr ⟨AddSubgroup.mem_top _, fun s => by
      rcases eq_or_ne v s with rfl | hvs
      · simp
      · simp [hvs]⟩⟩ with hxx
  have h0 : xx = 0 := eq_zero_of_isGroupLike_of_isSharp hcc (isSharp_effR _) xx
  have h1 : (Finsupp.single v (1:ℝ)) = 0 := congrArg Subtype.val h0
  have h2 := congrArg (fun f : ArithPlace ((botSub F Kbar).toIF) →₀ ℝ => f v) h1
  simp at h2

/-! ## ★4. ★★★★★★★主定理 -/

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)** —— 算術の **`𝒞^rlf` は rationally standard 型**。

原文 (FrdI p.114):
> Then the Frobenioids C, Cpf, Crlf, Cun-tr, (Cpf)un-tr are of isotropic and
> rationally standard type, but not of group-like type

★★★これで原文の **5 圏がそろった**(残る 4 圏は `Ex63RatStd.lean`)。

★★証明は `ex63_ratStd` と同じ形で、`Γ_L` を `⊤` に替えただけである。 -/
theorem rlf_ratStd :
    IsOfRationallyStandardType (ModelData.modelPre (rlfHyp F Kbar))
      (ModelData.model_frobenioid (rlfHyp F Kbar)) (rlfIota F Kbar) := by
  classical
  set L := (botSub F Kbar).toIF with hL
  set Y : (FinSub F Kbar)ᵒᵖ := Opposite.op (botSub F Kbar) with hY
  set M := rlfModelData F Kbar with hM
  set h := rlfHyp F Kbar with hh
  set P := ModelData.modelPre h with hP
  set Gm := ModelData.model_frobenioid h with hGm
  set Fc := Gm.core with hFc
  set v0 : ArithPlace L := Sum.inr (Classical.arbitrary (InfinitePlace L)) with hv0
  obtain ⟨y, hy0⟩ := exists_arithDiv_ne_zero_at L v0
  obtain ⟨a0, b0, ha0, hb0, hab0⟩ := exists_effR_sub
    (⊤ : AddSubgroup (ArithPlace L →₀ ℝ)) isGenSubgroupR_top (arithDiv y)
    (AddSubgroup.mem_top _)
  set a : M.phi.val Y := ⟨a0, ha0⟩ with hadef
  set b : M.phi.val Y := ⟨b0, hb0⟩ with hbdef
  have hkey : toGp (M.phi.val Y) a - toGp (M.phi.val Y) b = M.divB Y (Additive.ofMul y) :=
    toGp_sub_eq_divBRlf Y (Additive.ofMul y) a b (by
      show a0 - b0 = arithDiv y
      rw [hab0])
  set X : ModelData.Obj M := ⟨Y, 0⟩ with hXdef
  set Am : ModelData.Obj M := ⟨Y, toGp (M.phi.val Y) b⟩ with hAmdef
  set gm : X ⟶ Am :=
    { base := 𝟙 Y, div := b, deg := 1, u := 0,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) b
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y 0
        rw [M.phi.gpMapOn_id, map_zero, smul_zero, zero_add, add_zero] } with hgmdef
  set fm : X ⟶ Am :=
    { base := 𝟙 Y, div := a, deg := 1, u := Additive.ofMul y,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) a
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y (Additive.ofMul y)
        rw [M.phi.gpMapOn_id, smul_zero, zero_add, ← hkey]
        abel } with hfmdef
  set Xi : Istr P := ⟨X, ModelData.model_isotropicType h X⟩ with hXi
  set Ami : Istr P := ⟨Am, ModelData.model_isotropicType h Am⟩ with hAmi
  set Q := unTrPre P Fc with hQ
  set GQ := unTr_frobenioid P Fc Gm with hGQ
  set fu : (show UnTr P from Xi) ⟶ (show UnTr P from Ami) :=
    (istrToUnTr P).map (ObjectProperty.homMk fm) with hfu
  set gu : (show UnTr P from Xi) ⟶ (show UnTr P from Ami) :=
    (istrToUnTr P).map (ObjectProperty.homMk gm) with hgu
  haveI : IsIso (Q.Base fu) := by
    show IsIso (𝟙 Y)
    infer_instance
  have hmem := toGp_div_sub_mem_phiBiratAt_of_preStep Q GQ
    (fun Z => (unTr_isOfModelType Fc Gm).2 Z)
    (X := show UnTr P from Xi) (A := show UnTr P from Ami)
    (unTr_isotropic P Fc _) fu gu rfl rfl rfl
  refine ModelData.model_isOfRationallyStandardType_of_baseId h (rlfIota F Kbar)
    (ModelData.modelRatFnData h) (rlf_hsp F Kbar)
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    (arithDatumRlf_isNonDilatingOn F Kbar)
    (fun hgl => absurd hgl (rlf_not_isOfGroupLikeType F Kbar))
    (show BiratCat Q GQ from (show UnTr P from Xi))
    (fun θ => rlf_phi_map_bot_eq_id F Kbar (biratBase θ.inv))
    (M.divB Y (Additive.ofMul y)) ?_ ?_
  · rw [← hkey]
    exact hmem
  · intro n hn
    exact divBRlf_nsmul_ne_zero Y y v0 hy0 n hn

/-- ★★★★**[FrdI] Theorem 6.4, (i)** —— 算術の `𝒞^rlf` は **standard 型**。 -/
theorem rlf_standardType :
    IsOfStandardType ((FinSub F Kbar)ᵒᵖ) _
      (ModelData.modelPre (rlfHyp F Kbar))
      (ModelData.model_frobenioid (rlfHyp F Kbar)).core :=
  ModelData.model_isOfStandardType (rlfHyp F Kbar)
    (ModelData.model_frobenioid (rlfHyp F Kbar)).core
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    (arithDatumRlf_isNonDilatingOn F Kbar)
    (fun hgl => absurd hgl (rlf_not_isOfGroupLikeType F Kbar))

/-- ★★★★**[FrdI] Theorem 6.4, (i)** —— `𝒞^rlf` とその派生圏は **group-like 型でない**。 -/
theorem rlf_not_groupLike_family :
    ¬ IsOfGroupLikeType (ModelData.modelPre (rlfHyp F Kbar))
      ∧ ¬ IsOfGroupLikeType (pfRootPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (pfRootPre (ModelData.modelPre (rlfHyp F Kbar))
            (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
          (pfRootCore (F := (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
            (isOfFrobeniusIsotropicType_of_isotropic
              (ModelData.model_isotropicType (rlfHyp F Kbar))))) := by
  haveI := (rlfHyp F Kbar).connectedD
  have h0 := rlf_not_isOfGroupLikeType F Kbar
  exact ⟨h0, pfRoot_not_isOfGroupLikeType h0, unTr_not_isOfGroupLikeType _ h0,
    unTr_not_isOfGroupLikeType _ (pfRoot_not_isOfGroupLikeType h0)⟩

/-- ★★★★★**[FrdI] Theorem 6.4, (i)** —— `𝒞^rlf` とその派生圏の isotropic 型。 -/
theorem rlf_isotropic_family :
    IsOfIsotropicType (ModelData.modelPre (rlfHyp F Kbar))
      ∧ IsOfFrobeniusIsotropicType (ModelData.modelPre (rlfHyp F Kbar))
      ∧ IsOfIsotropicType (pfRootPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ IsOfIsotropicType (istrPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ IsOfFrobeniusIsotropicType (unTrPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ IsOfIsotropicType (unTrPre (ModelData.modelPre (rlfHyp F Kbar))
          (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
      ∧ IsOfIsotropicType
          (unTrPre (pfRootPre (ModelData.modelPre (rlfHyp F Kbar))
              (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
            (pfRootCore (F := (ModelData.model_frobenioid (rlfHyp F Kbar)).core)
              (isOfFrobeniusIsotropicType_of_isotropic
                (ModelData.model_isotropicType (rlfHyp F Kbar))))) :=
  ModelIso.model_isotropic_family (rlfHyp F Kbar)

end Std

/-! ### ★出典の紐付け -/

def rlf_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf は rationally standard 型(5 圏の最後の 1 つ)",
    sectionId := "frdi-thm-6-4" }

def rlf_ratStd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ModelData.model_isOfRationallyStandardType_of_baseId(一般形)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.model_isOfRationallyStandardType_of_baseId") 114,
    .citation "[ABC3]" "rlf_hsp(strictly rational の中身)"
      (.inProject "ABC3" "ABC3.Found.Divisor.rlf_hsp") 115,
    .citation "[ABC3]" "arithDatumRlf_isNonDilatingOn(non-dilating)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDatumRlf_isNonDilatingOn") 114,
    .implicitStep
      ("★逸脱: 𝒞^rlf を「𝒞 の実現化」ではなく grp = ⊤ の ArithDatum として建てている" ++
       "(Ex63Rlf.lean に記録)。ScT も rlfCone も ℝ≥0 の平坦性も使わない") 114 ]

def arithDatumRlf_isNonDilatingOn.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^rlf は non-dilating",
    sectionId := "frdi-thm-6-4" }

def rlf_hsp.src : Source :=
  { paper := "FrdI", pdfPage := 115,
    item := "Theorem 6.4, (i) — 𝒞^rlf は strictly rational 型",
    sectionId := "frdi-thm-6-4" }

def rlf_not_isOfGroupLikeType.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf は group-like 型でない",
    sectionId := "frdi-thm-6-4" }

def rlf_standardType.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf は standard 型",
    sectionId := "frdi-thm-6-4" }

def rlf_isotropic_family.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf とその派生圏の isotropic 型",
    sectionId := "frdi-thm-6-4" }

def rlf_not_groupLike_family.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 𝒞^rlf とその派生圏は group-like 型でない",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.Divisor
