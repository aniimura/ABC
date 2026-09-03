/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex61Units
import ABC3.Found.FrdI.Thm64RatStd

/-!
# Ex61RatStd —— `[FrdI] Theorem 6.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta
universe u
variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]
  [IsGalois V.functionField Kbar]
  (hkq : IsKQCartier V DK
    (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-! ## ★1. `Div_B` と有効因子の差の橋 -/

/-- ★★**有効因子の差が `div(y)` なら、その差は `Div_B(y)` である**。

★`Φ(L)^gp ≃ Γ`(`effSubGpEquiv`)の単射性 1 本。 -/
theorem toGp_sub_eq_divBGeom (A : (FinSub V.functionField Kbar)ᵒᵖ)
    (y : (bmonGeom V DK).val A)
    (a b : effSub ((cartierDatumGeom V DK hkq).grp A))
    (hab : ((a : ((cartierDatumGeom V DK hkq).primes A) →₀ ℤ))
        - ((b : ((cartierDatumGeom V DK hkq).primes A) →₀ ℤ))
      = (((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
          : (DLSet V DK A.unop) →₀ ℤ)
          : ((cartierDatumGeom V DK hkq).primes A) →₀ ℤ)) :
    toGp _ a - toGp _ b = divBGeom V DK hkq A y := by
  refine effSubGpHom_injective ((cartierDatumGeom V DK hkq).grp A) ?_
  have h1 : effSubGpHom ((cartierDatumGeom V DK hkq).grp A) (toGp _ a - toGp _ b)
      = divBHom V DK A.unop y := by
    rw [effSubGpHom, map_sub, gpLift_toGp, gpLift_toGp]
    exact Subtype.ext hab
  rw [h1]
  show _ = (effSubGpEquiv ((cartierDatumGeom V DK hkq).grp A)
      ((cartierDatumGeom V DK hkq).qc A))
      ((effSubGpEquiv ((cartierDatumGeom V DK hkq).grp A)
        ((cartierDatumGeom V DK hkq).qc A)).symm (divBHom V DK A.unop y))
  rw [AddEquiv.apply_symm_apply]

/-- ★★**`Div_B(y)` は `div(y)` がある素因子で非零なら無限位数**。 -/
theorem divBGeom_nsmul_ne_zero (A : (FinSub V.functionField Kbar)ᵒᵖ)
    (y : (bmonGeom V DK).val A) (s : (DLSet V DK A.unop : Type u))
    (hy : ((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
        : (DLSet V DK A.unop) →₀ ℤ) s ≠ 0)
    (n : ℕ) (hn : 0 < n) :
    n • (divBGeom V DK hkq A y) ≠ 0 := by
  intro hc
  have h1 := congrArg (phiGpHomC (cartierDatumGeom V DK hkq) finSubOp_isOfFSMType) hc
  rw [map_nsmul, map_zero, phiGpHomC_divBGeom] at h1
  have h4 := congrArg
    (fun w : ((cartierDatumGeom V DK hkq).grp A) => (Subtype.val w) s) h1
  simp only [ZeroMemClass.coe_zero, Finsupp.coe_zero] at h4
  rcases mul_eq_zero.mp h4 with hz | hz
  · exact absurd (Nat.cast_eq_zero.mp hz) hn.ne'
  · exact hy hz

/-! ## ★2. group-like でないこと -/

/-- ★★★**[FrdI] Theorem 6.2, (iii)** —— 幾何 Frobenioid は **group-like 型でない**。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized

★原文の証明は「`D_K ≠ ∅`(と `Φ` の定義)から immediate」。
★★**逸脱の記録**: 原文の `D_K ≠ ∅` から `D_L ≠ ∅` を出す段
(有限射 `V[L] → V` の余次元 1 点での持ち上げ)は畳んだ。
本形式化は**その 1 点(`D_L ≠ ∅` なる `L` が 1 つある)を仮定として受ける**。
★`IsOfGroupLikeType` は全対象への `∀` なので、否定には対象 1 つで足りる。 -/
theorem ex61_not_isOfGroupLikeType (A₀ : (FinSub V.functionField Kbar)ᵒᵖ)
    (hne : Nonempty (DLSet V DK A₀.unop)) :
    ¬ IsOfGroupLikeType (ModelData.modelPre (ex61Hyp V DK hkq)) := by
  classical
  intro hgl
  obtain ⟨s⟩ := hne
  have hcc := hgl ⟨A₀, 0⟩
  obtain ⟨a, ha⟩ := exists_effSub_support_eq
    ((cartierDatumGeom V DK hkq).qc A₀) ({s} : Finset (DLSet V DK A₀.unop))
  have h0 : a = 0 := eq_zero_of_isGroupLike_of_isSharp hcc (isSharp_effSub _) a
  rw [h0] at ha
  simp only [ZeroMemClass.coe_zero, Finsupp.support_zero] at ha
  exact (Finset.singleton_ne_empty s) ha.symm

/-- ★★**`Φ` は `⊥` の上では恒等しか誘導しない** —— これが `σ = id` の中身。 -/
theorem ex61_phi_map_bot_eq_id
    (e : (Opposite.op (botSub V.functionField Kbar))
      ⟶ (Opposite.op (botSub V.functionField Kbar))) :
    (ex61ModelData V DK hkq).phi.map e
      = AddMonoidHom.id ((ex61ModelData V DK hkq).phi.val
          (Opposite.op (botSub V.functionField Kbar))) := by
  have he : e = 𝟙 _ := Quiver.Hom.unop_inj (botSub_end_eq_id e.unop)
  rw [he]
  ext x
  exact MonoidOn.map_id _ _ x

/-! ## ★3. `Definition 4.5, (ii)` の入力 -/

/-- ★★**因子分解写像 `ι`** —— `Definition 2.4, (i)` の (c)(d)。`iotaAt` そのもの。 -/
noncomputable def ex61Iota :
    ∀ Y : (FinSub V.functionField Kbar)ᵒᵖ,
      Prime ((ex61ModelData V DK hkq).phi.val Y)
      → Pf ((ex61ModelData V DK hkq).phi.val Y) → NNReal :=
  fun Y => iotaAt ((cartierDatumGeom V DK hkq).qc Y)

/-- ★★★★**[FrdI] Definition 4.5, (ii) の入力**(`strictly rational` の中身)。

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type. -/
theorem ex61_hsp
    (hsupp : ∀ (A : (FinSub V.functionField Kbar)ᵒᵖ) (s : (DLSet V DK A.unop : Type u)),
      ∃ y : (bmonGeom V DK).val A,
        ((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
          : (DLSet V DK A.unop) →₀ ℤ) s ≠ 0) :
    ∀ (A : ModelData.Obj (ex61ModelData V DK hkq))
      (p : Prime ((ex61ModelData V DK hkq).phi.val
        ((ModelData.modelPre (ex61Hyp V DK hkq)).toElem.obj A).base)),
      ∃ (a b : (ex61ModelData V DK hkq).phi.val
          ((ModelData.modelPre (ex61Hyp V DK hkq)).toElem.obj A).base)
        (y : (ModelData.modelRatFnData (ex61Hyp V DK hkq)).bmon.val
          ((ModelData.modelPre (ex61Hyp V DK hkq)).toElem.obj A).base),
        (toGp _ a - toGp _ b
            = (ModelData.modelRatFnData (ex61Hyp V DK hkq)).divB _ y ∨
          toGp _ a - toGp _ b
            = -((ModelData.modelRatFnData (ex61Hyp V DK hkq)).divB _ y)) ∧
        p ∈ SuppElt (ex61Iota V DK hkq _) a ∧ p ∉ SuppElt (ex61Iota V DK hkq _) b := by
  classical
  intro A p
  obtain ⟨y, hy⟩ := hsupp A.base
    (effSubPrimeEquiv ((cartierDatumGeom V DK hkq).qc A.base) p)
  obtain ⟨a, b, z, hzmem, hzeq, hab, hap, hbp⟩ :=
    exists_split_suppElt_of_qc ((cartierDatumGeom V DK hkq).qc A.base)
      (divBHom V DK A.base.unop y).2 p hy
  refine ⟨a, b, y, ?_, hap, hbp⟩
  have hsub : ((a : effSub ((cartierDatumGeom V DK hkq).grp A.base))
        : ((cartierDatumGeom V DK hkq).primes A.base) →₀ ℤ)
      - ((b : effSub ((cartierDatumGeom V DK hkq).grp A.base))
        : ((cartierDatumGeom V DK hkq).primes A.base) →₀ ℤ) = z := by
    rw [hab]; abel
  rcases hzeq with rfl | rfl
  · exact Or.inl (toGp_sub_eq_divBGeom V DK hkq A.base y a b hsub)
  · have h2 := toGp_sub_eq_divBGeom V DK hkq A.base y b a (by rw [← neg_sub, hsub]; abel)
    have h3 : toGp (effSub ((cartierDatumGeom V DK hkq).grp A.base)) a
        - toGp (effSub ((cartierDatumGeom V DK hkq).grp A.base)) b
        = -(divBGeom V DK hkq A.base y) := by rw [← h2]; abel
    exact Or.inr h3

/-! ## ★4. `Theorem 6.2, (iii)` の 5 条 -/

/-- ★★★★**[FrdI] Theorem 6.2, (iii)** —— 幾何 Frobenioid は **standard 型**。 -/
theorem ex61_standardType
    (A₀ : (FinSub V.functionField Kbar)ᵒᵖ) (hne : Nonempty (DLSet V DK A₀.unop)) :
    IsOfStandardType ((FinSub V.functionField Kbar)ᵒᵖ) _
      (ModelData.modelPre (ex61Hyp V DK hkq))
      (ModelData.model_frobenioid (ex61Hyp V DK hkq)).core :=
  ModelData.model_isOfStandardType (ex61Hyp V DK hkq)
    (ModelData.model_frobenioid (ex61Hyp V DK hkq)).core
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    ((cartierDatumGeom V DK hkq).isNonDilatingOn finSubOp_isOfFSMType
      (fun _ e => finSubOp_isIso_of_endo e))
    (fun hgl => absurd hgl (ex61_not_isOfGroupLikeType V DK hkq A₀ hne))

/-- ★★★★★★**[FrdI] Theorem 6.2, (iii)** ——
**幾何 Frobenioid `C_{V,K̄,D_K}` は rationally standard 型**。

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type.

★★手順は `Example 6.3`(算術)と同じ:
底 `⊥`(= `K` 自身)の上の対象を選んで `σ = id` を無料にし、
`Φ^birat` の無限位数の元は `Div_B(y)` を **2 本の pre-step の差**として作る。 -/
theorem ex61_ratStd
    (hsupp : ∀ (A : (FinSub V.functionField Kbar)ᵒᵖ) (s : (DLSet V DK A.unop : Type u)),
      ∃ y : (bmonGeom V DK).val A,
        ((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
          : (DLSet V DK A.unop) →₀ ℤ) s ≠ 0)
    (s0 : (DLSet V DK (botSub V.functionField Kbar) : Type u)) :
    IsOfRationallyStandardType (ModelData.modelPre (ex61Hyp V DK hkq))
      (ModelData.model_frobenioid (ex61Hyp V DK hkq)) (ex61Iota V DK hkq) := by
  classical
  set Y : (FinSub V.functionField Kbar)ᵒᵖ := Opposite.op (botSub V.functionField Kbar) with hY
  set M := ex61ModelData V DK hkq with hM
  set h := ex61Hyp V DK hkq with hh
  set P := ModelData.modelPre h with hP
  set G' := ModelData.model_frobenioid h with hG'
  set Fc := G'.core with hFc
  obtain ⟨y, hy⟩ := hsupp Y s0
  obtain ⟨a0, b0, ha0, hb0, hab0⟩ := exists_effSub_sub
    ((cartierDatumGeom V DK hkq).grp Y) ((cartierDatumGeom V DK hkq).qc Y)
    (((divBHom V DK Y.unop y : cartierOnDL V DK Y.unop _)
      : (DLSet V DK Y.unop) →₀ ℤ) : ((cartierDatumGeom V DK hkq).primes Y) →₀ ℤ)
    (divBHom V DK Y.unop y).2
  set a : M.phi.val Y := ⟨a0, ha0⟩ with hadef
  set b : M.phi.val Y := ⟨b0, hb0⟩ with hbdef
  have hkey : toGp (M.phi.val Y) a - toGp (M.phi.val Y) b = M.divB Y y :=
    toGp_sub_eq_divBGeom V DK hkq Y y a b (by rw [hab0])
  set X : ModelData.Obj M := ⟨Y, 0⟩ with hXdef
  set A' : ModelData.Obj M := ⟨Y, toGp (M.phi.val Y) b⟩ with hA'def
  set gm : X ⟶ A' :=
    { base := 𝟙 Y, div := b, deg := 1, u := 0,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) b
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y 0
        rw [M.phi.gpMapOn_id, map_zero, smul_zero, zero_add, add_zero] } with hgmdef
  set fm : X ⟶ A' :=
    { base := 𝟙 Y, div := a, deg := 1, u := y,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) a
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y y
        rw [M.phi.gpMapOn_id, smul_zero, zero_add, ← hkey]
        abel } with hfmdef
  set Xi : Istr P := ⟨X, ModelData.model_isotropicType h X⟩ with hXi
  set A'i : Istr P := ⟨A', ModelData.model_isotropicType h A'⟩ with hA'i
  set Q := unTrPre P Fc with hQ
  set GQ := unTr_frobenioid P Fc G' with hGQ
  set fu : (show UnTr P from Xi) ⟶ (show UnTr P from A'i) :=
    (istrToUnTr P).map (ObjectProperty.homMk fm) with hfu
  set gu : (show UnTr P from Xi) ⟶ (show UnTr P from A'i) :=
    (istrToUnTr P).map (ObjectProperty.homMk gm) with hgu
  haveI : IsIso (Q.Base fu) := by
    show IsIso (𝟙 Y)
    infer_instance
  have hmem := toGp_div_sub_mem_phiBiratAt_of_preStep Q GQ
    (fun Z => (unTr_isOfModelType Fc G').2 Z)
    (X := show UnTr P from Xi) (A := show UnTr P from A'i)
    (unTr_isotropic P Fc _) fu gu rfl rfl rfl
  refine ModelData.model_isOfRationallyStandardType_of_baseId h (ex61Iota V DK hkq)
    (ModelData.modelRatFnData h) (ex61_hsp V DK hkq hsupp)
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    ((cartierDatumGeom V DK hkq).isNonDilatingOn finSubOp_isOfFSMType
      (fun _ e => finSubOp_isIso_of_endo e))
    (fun hgl => absurd hgl (ex61_not_isOfGroupLikeType V DK hkq Y ⟨s0⟩))
    (show BiratCat Q GQ from (show UnTr P from Xi))
    (fun θ => ex61_phi_map_bot_eq_id V DK hkq (biratBase θ.inv))
    (M.divB Y y) ?_ ?_
  · rw [← hkey]
    exact hmem
  · intro n hn
    exact divBGeom_nsmul_ne_zero V DK hkq Y y s0 hy n hn

/-- ★★★★★★★**[FrdI] Theorem 6.2, (iii)** —— 5 条まとめ。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized -/
theorem ex61_thm62_iii_types
    (hsupp : ∀ (A : (FinSub V.functionField Kbar)ᵒᵖ) (s : (DLSet V DK A.unop : Type u)),
      ∃ y : (bmonGeom V DK).val A,
        ((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
          : (DLSet V DK A.unop) →₀ ℤ) s ≠ 0)
    (s0 : (DLSet V DK (botSub V.functionField Kbar) : Type u)) :
    IsOfIsotropicType (ModelData.modelPre (ex61Hyp V DK hkq))
      ∧ IsOfStandardType ((FinSub V.functionField Kbar)ᵒᵖ) _
          (ModelData.modelPre (ex61Hyp V DK hkq))
          (ModelData.model_frobenioid (ex61Hyp V DK hkq)).core
      ∧ IsOfBirationallyFrobeniusNormalizedType _
          (ModelData.modelPre (ex61Hyp V DK hkq))
          (ModelData.model_frobenioid (ex61Hyp V DK hkq))
      ∧ ¬ IsOfGroupLikeType (ModelData.modelPre (ex61Hyp V DK hkq))
      ∧ IsOfRationallyStandardType (ModelData.modelPre (ex61Hyp V DK hkq))
          (ModelData.model_frobenioid (ex61Hyp V DK hkq)) (ex61Iota V DK hkq) :=
  ⟨ModelData.model_isotropicType (ex61Hyp V DK hkq),
   ex61_standardType V DK hkq (Opposite.op (botSub V.functionField Kbar)) ⟨s0⟩,
   ModelData.model_isOfBiratFrobNormalizedType (ex61Hyp V DK hkq),
   ex61_not_isOfGroupLikeType V DK hkq (Opposite.op (botSub V.functionField Kbar)) ⟨s0⟩,
   ex61_ratStd V DK hkq hsupp s0⟩

def ex61_not_isOfGroupLikeType.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 幾何 Frobenioid は group-like 型でない",
    sectionId := "frdi-thm-6-2" }

def ex61_not_isOfGroupLikeType.needs : List ProofObligation :=
  [ .implicitStep
      ("★原文は「D_K ≠ ∅ から immediate」と書く。★D_K ≠ ∅ から D_L ≠ ∅ を出す段" ++
       "(有限射 V[L] → V の余次元 1 点での持ち上げ)は畳まれている。" ++
       "本形式化はその 1 点(D_L ≠ ∅ なる L が 1 つある)を仮定として受ける") 111 ]

def ex61_standardType.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 幾何 Frobenioid は standard 型",
    sectionId := "frdi-thm-6-2" }

def ex61_hsp.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 幾何 Frobenioid は strictly rational 型",
    sectionId := "frdi-thm-6-2" }

def ex61_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 幾何 Frobenioid は rationally standard 型",
    sectionId := "frdi-thm-6-2" }

def ex61_ratStd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "model_isOfRationallyStandardType_of_baseId"
      (.inProject "ABC3"
        "ABC3.Found.FrdI.ModelData.model_isOfRationallyStandardType_of_baseId") 111,
    .citation "[ABC3]" "toGp_div_sub_mem_phiBiratAt_of_preStep(Φ^birat の元を作る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.toGp_div_sub_mem_phiBiratAt_of_preStep") 111,
    .implicitStep
      ("★★訂正(2026-08-25): 在庫の birat_isFrobeniusCompact_of_primary は" ++
       "「Φ の非零元が Φ^birat に入る」を要求するが、V が proper なので" ++
       "有効な主因子は 0 だけ(forall_ordPt_eq_zero_of_forall_nonneg)であり、" ++
       "幾何でも満たしようがない。底 ⊥ の対象を選ぶ版に差し替えた") 111 ]

def ex61_thm62_iii_types.src : Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 5 条(isotropic / standard / birat-Frob-normalized / not group-like / rationally standard)",
    sectionId := "frdi-thm-6-2" }

/-- ★★★★★★★★**[FrdI] Theorem 6.2** —— 4 条がすべて実装された。

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | `L₂ = L₁·K₂` が次数を保つ | `SepClosedFunctor.lean` |
| (i) | `ψ` が `𝒟₁ → 𝒟₂` を定める | `SepClosedFunctor.lean` |
| (i) | `Φ₁ → Φ₂\|𝒟₁`、`B₁ → B₂\|𝒟₁` | `Thm62Pull.lean` / `Thm62Hom.lean` |
| (i) | `Ψ : 𝒞₁ ⥤ 𝒞₂`(Frobenius 次数と両立) | `GeomPullDatum.psiFunctor` / `ModelDataHom.hom_deg` |
| (ii) | 標数 `p` の Frobenius が定める `Ψ` は naive Frobenius 関手と同型 | `Thm62Frob.lean` |
| (iii) | isotropic / standard / birat-Frob-normalized / not group-like / rationally standard | `ex61_thm62_iii_types` |
| (iv) | `𝒟` は Frobenius-slim | `isFrobeniusSlim_of_mulEquiv_subgroup` |
| (iv) | `𝒟` が slim ⟺ `Z = {1}` | `isSlimCat_iff_of_mulEquiv` |
| (iv) | `𝒟` が Div-slim ⟺ `z` が `Φ(L)` に非自明に作用 | `isDivSlim_iff_of_mulEquiv` |

★(iv) は原文と同じく [Mzk7] `Corollary 1.1.6`(`Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)`)を
**仮定として受ける**(原文の証明もそこを引いている)。 -/
def thm_6_2.src : Source :=
  { paper := "FrdI", pdfPage := 110, item := "Theorem 6.2",
    sectionId := "frdi-thm-6-2" }

def thm_6_2.needs : List ProofObligation :=
  [ .otherPaper "[Mzk7]" "Corollary 1.1.6 — Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)。★(iv) が仮定として受ける" 112,
    .citation "[ABC3]" "(iii) の 5 条"
      (.inProject "ABC3" "ABC3.Found.Divisor.ex61_thm62_iii_types") 111,
    .citation "[ABC3]" "(i) の Ψ"
      (.inProject "ABC3" "ABC3.Found.Divisor.GeomPullDatum.psiFunctor") 110,
    .citation "[ABC3]" "(ii) の naive Frobenius 関手との同型"
      (.inProject "ABC3" "ABC3.Found.FrdI.naiveFrobIso_of_unit") 111 ]

end ABC3.Found.Divisor
