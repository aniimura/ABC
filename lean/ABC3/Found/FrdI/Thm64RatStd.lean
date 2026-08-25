/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm62RatStd

/-!
# [FrdI] Theorem 6.4, (i) —— `(𝒞^un-tr)^birat` の Frobenius-compact 対象(算術版)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.115):
> C is of [strictly] rational type, and that every object of (Cun-tr)birat is Frobenius-

## ★★★★測定(2026-08-25)—— 幾何版の入口は算術では**閉じている**

在庫の `birat_isFrobeniusCompact_of_primary`(`Thm62RatStd.lean`)は 2 つを要求する:

| 仮引数 | 中身 |
|---|---|
| `hprim` | **準素元**が `Φ^birat` に入る |
| `x₀` / `hx₀mem` / `hx₀` | `Φ` の**非零元**が `Φ^birat` に入る |

★★★どちらも**算術 (`Example 6.3`) では成り立たない**。理由は 1 つ:

* `RatFnData` の `κ` が全単射なので `Φ^birat(L) = Div_B(B(L))`(主因子)ちょうどである。
* 算術の主因子は**積公式で次数 0**、`Φ(L)` は**有効**(全成分 ≥ 0)。
* `arithDegree` は全成分の**正係数和**なので、`Φ(L) ∩ Φ^birat(L) = {0}`。

★つまり `x₀ ∈ Φ(L)`、`x₀ ≠ 0`、`toGp x₀ ∈ Φ^birat(L)` を同時に満たす元は**無い**。

## ★★★★★訂正(2026-08-25、**この節の初版は誤っていた**)

初版はここに「幾何 (`Example 6.1`) では `V[L]` 上の正則関数の因子が有効な主因子を
与えるので通る」と書いたが、**誤りである**。
★`V` は **proper** なので `Γ(V[L], 𝒪)^× = k_L^×` であり、
**有効な主因子は 0 だけ**である —— 在庫の

    forall_ordPt_eq_zero_of_forall_nonneg   (`Ex61Units.lean`)
    ex61_otri_le_otimes(`𝒪^▷(A) = 𝒪^×(A)`)

がまさにそれを言っている。
★★したがって `Φ(L) ∩ Φ^birat(L) = {0}` は **`Example 6.1` でも成り立つ**。
★★★**本ファイルの逃げ道は幾何・算術の両方で必要**である
(`Ex61RatStd.lean` / `Ex63RatStd.lean` の両方が使う)。

## ★★本ファイルの逃げ道 —— 「底の自己射が恒等しかない対象」を選ぶ

`birat_cd_eq_of_primary` が `hprim` を使うのは **`σ := Φ(Base θ⁻¹) = id`** を出す 1 箇所だけ
(`addMonoidHom_eq_id_of_primary_mprec`)。★`Definition 4.5, (iii), (b)` は
**「Frobenius-compact な対象が 1 つある」**としか言わないので、

    Base(A) の自己射が恒等だけ

であるような `A` を選べば `σ = id` は**無料**である。
★★`Example 6.3` の `𝒟 = B(G)⁰` では **`⊥`(= `F` 自身)** がそれである
(`F`-代数の自己射は恒等しかない)。

★もう 1 つ、`x₀` を `Φ` の元から **`Φ^gp` の元 `z`** へ緩める。
`birat_cd_eq_of_primary` は `x₀` を `toGp x₀` の形でしか使っていないので、
これは**証明を 1 行も損なわない**。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `toGp_div_sub_mem_phiBiratAt_of_preStep` | ★2 本の pre-step の `Div` の差は `Φ^birat` に入る |
| `birat_frobeniusCompact_cond2_of_mem_gp` | ★第 2 条(`Φ^gp` の元でよい) |
| `birat_cd_eq_of_baseId` | ★★`σ = id` なら `c = d`(準素元の仮定なし) |
| `birat_isFrobeniusCompact_of_baseId` | ★★★`𝒞^birat` の Frobenius-compact 対象 |
| `unTrBiratCompact_of_baseId` | ★★★★`(𝒞^un-tr)^birat` の Frobenius-compact 対象 |
| `ModelData.model_isOfRationallyStandardType_of_baseId` | ★★★★★★**model は rationally standard** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)

/-! ## ★1. pre-step の対から `Φ^birat` の元を作る -/

section PreStep

variable (G : Frobenioid P)

/-- ★★★★**同じ底を持つ 2 本の pre-step の `Div` の差は `Φ^birat` に入る**。

★`toGp_div_sub_mem_phiBiratAt`(在庫)は `𝒪^▷(X)`(始域 = 終域)の場合だが、
**始域と終域が違ってよい**のが本補題である。
★★これが `Φ^birat = Div_B(B)` を「有効因子の差」として掴む唯一の入口になる。 -/
theorem toGp_div_sub_mem_phiBiratAt_of_preStep
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {X A : C} (hX : IsIsotropic P X) (f g : X ⟶ A)
    (hfl : P.degFr f = 1) (hgl : P.degFr g = 1)
    (hb : P.Base f = P.Base g) [IsIso (P.Base f)] :
    toGp _ (P.Div f) - toGp _ (P.Div g) ∈ phiBiratAt P G (show BiratCat P G from X) := by
  haveI hXb : IsIsotropic (biratPre P G) (show BiratCat P G from X) :=
    birat_isotropic_up P G hX
  haveI hbg : IsIso (P.Base g) := by rw [← hb]; infer_instance
  have hpf : IsPreStep (biratPre P G) ((toBiratCat P G).map f) := by
    refine ⟨?_, ?_⟩
    · show biratDeg (toHomBirat (P := P) (G := G) f) = 1
      rw [biratDeg_toHomBirat, hfl]
    · show IsIso (biratBase (toHomBirat (P := P) (G := G) f))
      rw [biratBase_toHomBirat]
      infer_instance
  have hpg : IsPreStep (biratPre P G) ((toBiratCat P G).map g) := by
    refine ⟨?_, ?_⟩
    · show biratDeg (toHomBirat (P := P) (G := G) g) = 1
      rw [biratDeg_toHomBirat, hgl]
    · show IsIso (biratBase (toHomBirat (P := P) (G := G) g))
      rw [biratBase_toHomBirat]
      infer_instance
  haveI : IsIso ((toBiratCat P G).map f) :=
    birat_isIso_of_preStep_of_isotropic P G hfn hXb hpf
  haveI : IsIso ((toBiratCat P G).map g) :=
    birat_isIso_of_preStep_of_isotropic P G hfn hXb hpg
  have hbb : biratBase ((toBiratCat P G).map f) = biratBase ((toBiratCat P G).map g) := by
    show biratBase (toHomBirat (P := P) (G := G) f) = biratBase (toHomBirat (P := P) (G := G) g)
    rw [biratBase_toHomBirat, biratBase_toHomBirat, hb]
  have h := birat_divGp_sub_mem P G
    ((toBiratCat P G).map f) ((toBiratCat P G).map g) hbb
  rwa [show biratDivGp ((toBiratCat P G).map f) = toGp _ (P.Div f) from
      biratDivGp_toHomBirat _,
    show biratDivGp ((toBiratCat P G).map g) = toGp _ (P.Div g) from
      biratDivGp_toHomBirat _] at h

/-- ★★★locator —— `Proposition 4.4, (iii)` の pre-step 対版。 -/
def toGp_div_sub_mem_phiBiratAt_of_preStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 同じ底の 2 本の pre-step の Div の差は Φ^birat に入る",
    sectionId := "frdi-prop-4-4" }

end PreStep

/-! ## ★2. `Φ^gp` の元で足りる 2 条 -/

section BaseId

variable (G : Frobenioid P)

/-- ★★**第 2 条(無限位数の単元)** —— `Φ^birat` の**無限位数の元**があればよい。

★在庫の `birat_frobeniusCompact_cond2_of_mem` は `x₀ ∈ Φ` の像として書かれているが、
証明は `toGp x₀` の形しか使っていない。★算術では `Φ ∩ Φ^birat = {0}` なので
**`Φ^gp` の元へ緩めることが不可欠**である。 -/
theorem birat_frobeniusCompact_cond2_of_mem_gp (A : BiratCat P G)
    (z : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base))
    (hzmem : z ∈ phiBiratAt P G A)
    (hz : ∀ n : ℕ, 0 < n → n • z ≠ 0) :
    ∃ u : End A, u ∈ OTimes (biratPre P G) A ∧
      ∀ k : ℕ+, (u ^ ((k : ℕ+) : ℕ) : End A) ≠ 1 := by
  obtain ⟨u, hu, huz⟩ := hzmem
  refine ⟨u, hu, fun k hk => ?_⟩
  have h := congrArg (fun t : End A => biratDivGp ((t : A ⟶ A))) hk
  rw [biratDivGp_pow hu, huz] at h
  refine hz ((k : ℕ+) : ℕ) k.2 ?_
  rw [h]
  exact biratDivGp_id A

/-- ★★★★**`c = d`(底の自己射が恒等の版)**。

★在庫の `birat_cd_eq_of_primary` が `hprim`(準素元が `Φ^birat` に入る)を使うのは
**`σ := Φ(Base θ⁻¹) = id` を出す 1 箇所だけ**である。
★★底の自己射が恒等しかない対象を選べば `σ = id` は仮定として直接置ける。 -/
theorem birat_cd_eq_of_baseId
    (A : BiratCat P G) (θ : A ≅ A) (c d : ℕ+)
    (hσ : Φ.map (biratBase θ.inv)
      = AddMonoidHom.id (Φ.val (P.toElem.obj (biratDown P G A)).base))
    (z : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base))
    (hzmem : z ∈ phiBiratAt P G A)
    (hz : ∀ n : ℕ, 0 < n → n • z ≠ 0)
    (hyp : ∀ u : End A, u ∈ OTimes (biratPre P G) A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) :
    c = d := by
  obtain ⟨u, hu, huz⟩ := hzmem
  obtain ⟨k, hk⟩ := hyp u hu
  have hcj : (endConj θ u) ∈ OTimes (biratPre P G) A :=
    endConj_mem_otimes (biratPre P G) θ hu
  have h := congrArg (fun t : End A => biratDivGp ((t : A ⟶ A))) hk
  rw [biratDivGp_pow hcj, biratDivGp_pow hu, biratDivGp_endConj θ hu, huz, hσ, gpMap_id] at h
  have hk' : (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • z
      = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • z := h
  by_contra hne
  have hkpos : 0 < ((k : ℕ+) : ℕ) := k.2
  rcases Nat.lt_or_ge (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
      (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) with hlt | hge
  · exact hz _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hlt) hk')
  · have hne' : ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) ≠ ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by
      intro hcc
      exact hne (PNat.coe_injective (Nat.eq_of_mul_eq_mul_right hkpos hcc)).symm
    have hgt : ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) < ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by omega
    exact hz _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hgt) hk'.symm)

/-- ★★★★★**`𝒞^birat` の Frobenius-compact 対象(底の自己射が恒等の版)**。 -/
theorem birat_isFrobeniusCompact_of_baseId
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (A : BiratCat P G)
    (hbase : ∀ θ : A ≅ A, Φ.map (biratBase θ.inv)
      = AddMonoidHom.id (Φ.val (P.toElem.obj (biratDown P G A)).base))
    (z : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base))
    (hzmem : z ∈ phiBiratAt P G A)
    (hz : ∀ n : ℕ, 0 < n → n • z ≠ 0) :
    IsFrobeniusCompact (biratPre P G) A :=
  ⟨birat_frobeniusCompact_cond1 hfn A,
   birat_frobeniusCompact_cond2_of_mem_gp P G A z hzmem hz,
   fun θ c d hyp =>
     frobeniusCompact_cond3_of_eq (biratPre P G) θ c d
       (birat_cd_eq_of_baseId P G A θ c d (hbase θ) z hzmem hz hyp) hyp⟩

/-- ★★★locator —— `Definition 1.2, (iv)` の Frobenius-compact 対象(底の自己射が恒等)。 -/
def birat_isFrobeniusCompact_of_baseId.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — (C^un-tr)^birat の Frobenius-compact 対象",
    sectionId := "frdi-thm-6-4" }

/-- ★★★**`(𝒞^un-tr)^birat` の Frobenius-compact 対象(底の自己射が恒等の版)**。 -/
theorem unTrBiratCompact_of_baseId (G' : Frobenioid P)
    (hfn : ∀ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'),
      IsFrobeniusNormalized (unTrBiratPre P Fc G') X)
    (A : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'))
    (hbase : ∀ θ : A ≅ A, Φ.map (biratBase θ.inv)
      = AddMonoidHom.id (Φ.val ((unTrPre P Fc).toElem.obj
          (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G') A)).base))
    (z : Gp (Φ.val ((unTrPre P Fc).toElem.obj
      (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G') A)).base))
    (hzmem : z ∈ phiBiratAt (unTrPre P Fc) (unTr_frobenioid P Fc G') A)
    (hz : ∀ n : ℕ, 0 < n → n • z ≠ 0) :
    ∃ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'),
      IsFrobeniusCompact (unTrBiratPre P Fc G') X :=
  ⟨A, birat_isFrobeniusCompact_of_baseId (unTrPre P Fc) (unTr_frobenioid P Fc G')
    hfn A hbase z hzmem hz⟩

end BaseId

/-! ## ★3. model Frobenioid が rationally standard 型(算術版) -/

namespace ModelData

variable {M : ModelData.{v, u, w} D}

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)** —— model Frobenioid が
**rationally standard 型**(算術で通る版)。

★在庫の `model_isOfRationallyStandardType_of_primary` との違いは
`hprim` / `hx₀mem` を **`hbase`(底の自己射が恒等)＋ `Φ^gp` の無限位数の元** に
差し替えたことだけである。

原文 (FrdI p.115):
> C is of [strictly] rational type, and that every object of (Cun-tr)birat is Frobenius- -/
theorem model_isOfRationallyStandardType_of_baseId (h : Hyp M)
    (ι : ∀ Y : D, Prime (M.phi.val Y) → Pf (M.phi.val Y) → NNReal)
    (R : RatFnData (modelPre h) (model_frobenioid h))
    (hsp : ∀ (A : Obj M) (p : Prime (M.phi.val ((modelPre h).toElem.obj A).base)),
      ∃ (a b : M.phi.val ((modelPre h).toElem.obj A).base)
        (y : R.bmon.val ((modelPre h).toElem.obj A).base),
        (toGp _ a - toGp _ b = R.divB _ y ∨ toGp _ a - toGp _ b = -(R.divB _ y)) ∧
        p ∈ SuppElt (ι _) a ∧ p ∉ SuppElt (ι _) b)
    (hfsmff : IsOfFSMFFType D) (hnd : M.phi.IsNonDilatingOn)
    (hgl : IsOfGroupLikeType (modelPre h) →
      ∃ A, IsFrobeniusCompact (istrPre (modelPre h) (model_frobenioid h).core) A)
    (A : BiratCat (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)))
    (hbase : ∀ θ : A ≅ A, M.phi.map (biratBase θ.inv)
      = AddMonoidHom.id (M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
          (biratDown _ _ A)).base))
    (z : Gp (M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
      (biratDown _ _ A)).base))
    (hzmem : z ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)) A)
    (hz : ∀ n : ℕ, 0 < n → n • z ≠ 0) :
    haveI := h.connectedD
    IsOfRationallyStandardType (modelPre h) (model_frobenioid h) ι :=
  haveI := h.connectedD
  (model_ratStandardType_iff h ι).mpr
    ⟨isOfRationalType_of_divB R ι hsp,
     model_isOfStandardType h (model_frobenioid h).core hfsmff hnd hgl,
     unTrBiratCompact_of_baseId (modelPre h) (model_frobenioid h).core (model_frobenioid h)
       (fun X => (unTr_isOfModelType (model_frobenioid h).core (model_frobenioid h)).2 X)
       A hbase z hzmem hz⟩

/-- ★★★★★★locator —— `Theorem 6.4, (i)` の rationally standard 型。 -/
def model_isOfRationallyStandardType_of_baseId.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 model Frobenioid は rationally standard 型",
    sectionId := "frdi-thm-6-4" }

end ModelData

/-! ## ★4. `𝒞^un-tr` の側 —— `Proposition 5.5, (iii)` の rationally standard の入口

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★`𝒞^un-tr` の `Φ^birat` は **`Div_B(B)` を丸ごと含む**。
中身は「同じ底の 2 本の pre-step の差」1 本で、model の `Φ^gp = Φ − Φ` を使うだけである
(`gp_eq_sub` に相当する `AddLocalization.induction_on`)。 -/

section UnTrRational

variable {M : ModelData.{v, u, w} D}

open ModelData in
/-- ★★★★★**`Div_B` の像は `𝒞^un-tr` の `Φ^birat` に必ず入る**。

★`X` を始域とする 2 本の pre-step

    f : X ⟶ A'   (Base = 𝟙、degFr = 1、Div = a、u = y)
    g : X ⟶ A'   (Base = 𝟙、degFr = 1、Div = b、u = 0)

を `A' := (Base X, cls X + toGp b)` の上に作る。`Div_B(y) = toGp a − toGp b` は
`Φ^gp = Φ − Φ` から取れるので、**追加の仮定はいらない**。 -/
theorem divB_mem_phiBiratAt_unTr (h : M.Hyp)
    (X : UnTr (modelPre h)) (y : M.bmon.val X.obj.base) :
    M.divB X.obj.base y ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h))
      (show BiratCat _ _ from X) := by
  haveI := h.connectedD
  obtain ⟨a, b, hab⟩ : ∃ a b : M.phi.val X.obj.base,
      M.divB X.obj.base y = toGp _ a - toGp _ b := by
    refine AddLocalization.induction_on (M.divB X.obj.base y) ?_
    rintro ⟨u, v⟩
    exact ⟨u, (v : M.phi.val X.obj.base), by
      rw [eq_sub_iff_add_eq]; exact mk_add_toGp _ u v⟩
  set A' : ModelData.Obj M :=
    ⟨X.obj.base, X.obj.cls + toGp (M.phi.val X.obj.base) b⟩ with hA'
  set gm : X.obj ⟶ A' :=
    { base := 𝟙 X.obj.base, div := b, deg := 1, u := 0,
      cond := by
        show ((1 : ℕ+) : ℕ) • X.obj.cls + toGp (M.phi.val X.obj.base) b
          = M.phi.gpMapOn (𝟙 X.obj.base)
              (X.obj.cls + toGp (M.phi.val X.obj.base) b) + M.divB X.obj.base 0
        rw [M.phi.gpMapOn_id, map_zero, add_zero]
        show (1 : ℕ) • X.obj.cls + _ = _
        rw [one_nsmul] } with hgm
  set fm : X.obj ⟶ A' :=
    { base := 𝟙 X.obj.base, div := a, deg := 1, u := y,
      cond := by
        show ((1 : ℕ+) : ℕ) • X.obj.cls + toGp (M.phi.val X.obj.base) a
          = M.phi.gpMapOn (𝟙 X.obj.base)
              (X.obj.cls + toGp (M.phi.val X.obj.base) b) + M.divB X.obj.base y
        rw [M.phi.gpMapOn_id, hab]
        show (1 : ℕ) • X.obj.cls + _ = _
        rw [one_nsmul]
        abel } with hfm
  set A'i : Istr (modelPre h) := ⟨A', ModelData.model_isotropicType h A'⟩ with hA'i
  set fu : X ⟶ (show UnTr (modelPre h) from A'i) :=
    (istrToUnTr (modelPre h)).map (ObjectProperty.homMk fm) with hfu
  set gu : X ⟶ (show UnTr (modelPre h) from A'i) :=
    (istrToUnTr (modelPre h)).map (ObjectProperty.homMk gm) with hgu
  haveI : IsIso ((unTrPre (modelPre h) (model_frobenioid h).core).Base fu) := by
    show IsIso (𝟙 X.obj.base)
    infer_instance
  have hmem := toGp_div_sub_mem_phiBiratAt_of_preStep
    (unTrPre (modelPre h) (model_frobenioid h).core)
    (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h))
    (fun Z => (unTr_isOfModelType (model_frobenioid h).core (model_frobenioid h)).2 Z)
    (X := X) (A := show UnTr (modelPre h) from A'i)
    (unTr_isotropic (modelPre h) (model_frobenioid h).core _) fu gu rfl rfl rfl
  rw [hab]
  exact hmem

open ModelData in
/-- ★★★★★★**`𝒞^un-tr` は strictly rational 型** —— model の `hsp` がそのまま効く。

★`Φ` は `𝒞` と `𝒞^un-tr` で**同じ**なので `ι` もそのまま使える。
★★`Proposition 5.5, (iii)` の rationally standard の 4 条のうち、
**中身のある 1 条**がこれである(他の 3 条は在庫と移送)。 -/
theorem unTr_isOfStrictlyRationalType_of_hsp (h : M.Hyp)
    (ι : ∀ Y : D, Prime (M.phi.val Y) → Pf (M.phi.val Y) → NNReal)
    (hsp : ∀ (A : ModelData.Obj M)
        (p : Prime (M.phi.val ((modelPre h).toElem.obj A).base)),
      ∃ (a b : M.phi.val ((modelPre h).toElem.obj A).base)
        (y : M.bmon.val ((modelPre h).toElem.obj A).base),
        (toGp _ a - toGp _ b = M.divB _ y ∨ toGp _ a - toGp _ b = -(M.divB _ y)) ∧
        p ∈ SuppElt (ι _) a ∧ p ∉ SuppElt (ι _) b) :
    IsOfStrictlyRationalType (UnTr (modelPre h))
      (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)) ι := by
  haveI := h.connectedD
  intro X p
  obtain ⟨a, b, y, hor, hap, hbp⟩ := hsp X.obj p
  refine ⟨a, b, ?_, hap, hbp⟩
  have hy := divB_mem_phiBiratAt_unTr h X y
  rcases hor with heq | heq
  · have h1 : toGp (M.phi.val X.obj.base) a - toGp (M.phi.val X.obj.base) b
        = M.divB X.obj.base y := heq
    have h2 : toGp (M.phi.val X.obj.base) a - toGp (M.phi.val X.obj.base) b
        ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
          (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h))
          (show BiratCat _ _ from X) := by rw [h1]; exact hy
    exact h2
  · have h1 : toGp (M.phi.val X.obj.base) a - toGp (M.phi.val X.obj.base) b
        = -(M.divB X.obj.base y) := heq
    have h2 : toGp (M.phi.val X.obj.base) a - toGp (M.phi.val X.obj.base) b
        ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
          (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h))
          (show BiratCat _ _ from X) := by
      rw [h1]; exact AddSubgroup.neg_mem _ hy
    exact h2

/-- ★★★★locator —— `Proposition 5.5, (iii)` の rationally standard の入口。 -/
def divB_mem_phiBiratAt_unTr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — Div_B の像は 𝒞^un-tr の Φ^birat に入る",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の `𝒞^un-tr` が strictly rational 型。 -/
def unTr_isOfStrictlyRationalType_of_hsp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr は strictly rational 型",
    sectionId := "frdi-prop-5-5" }

end UnTrRational

end ABC3.Found.FrdI
