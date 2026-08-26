/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32
import ABC3.Found.FrdI.Prop25
import ABC3.Found.FrdI.Prop55UntrIdx

/-!
# [FrdI] Proposition 5.5, (ii) —— 左辺を `(𝒞^pf)^un-tr` そのものに直す

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C

## ★★何が足りなかったか

`Prop55Untr.lean` / `Prop55UntrIdx.lean` は
「`Hom^pf(A,B)` を**ある段で `𝔽_Φ` の像が一致する**という関係で割ったもの」を左辺に取った。
★しかし原文の `(𝒞^pf)^un-tr` は、**`𝒞^pf` 自身の pre-Frobenioid 構造
(`pfPre`、`Proposition 3.2, (i)`)で unit-trivial 化したもの**である。
★★**その 2 つの同値関係が一致すること**が要る。

## ★★★一致の中身 —— divisorial から出る `n` 倍の単射性

段 `V` の脚を `a`, `b`、根を `r = deg_Fr(a)` とすると

| `𝒞^pf` の量 | 代表元での式 |
|---|---|
| `pfDeg` | `deg_Fr(φ)` |
| `pfBase` | `Base(a) ≫ Base(φ) ≫ Base(b)⁻¹` |
| `pfDiv` | `Base(a)^*(Div(φ)) / r`(`Φ^pf` の元) |

★`pfDeg`・`pfBase` からは `deg_Fr(φ)`・`Base(φ)` が**そのまま**戻る
(`Base(a)`, `Base(b)` は同型)。
★★**`pfDiv` から `Div(φ)` を戻すところだけが非自明**である ——
`Φ^pf` で `x/r = y/r` は「ある `k` で `(kr)·x = (kr)·y`」しか言わない。
★★★**そこで `n` 倍の単射性(`nsmul_injective_of_divisorial`)が効く** ——
integral・saturated・sharp の 3 つ、すなわち **`Φ` が divisorial であること**である。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `Pf.mk_eq_mk_same_iff` | `x/r = y/r ⟺ x = y`(divisorial のとき) |
| `exists_common_rep` | 2 元を同じ段の代表元で表す |
| `toQuot_eq_iff_pfToElem` | ★★**2 つの同値関係の一致** |
| `unTrPf_equiv_quot` | `Hom_{(𝒞^pf)^un-tr}(A,B) ≃ (前ファイルの商)` |
| `prop_5_5_ii_untr_cat` | ★★★★**`Hom_{(𝒞^pf)^un-tr} ≃ Hom_{(𝒞^un-tr)^pf}`** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

/-! ## ★1. `Φ^pf` の中の同分母の等式 -/

/-- ★★**divisorial なモノイドでは `x/r = y/r ⟺ x = y`**。

★`⟹` は `n` 倍の単射性(`nsmul_injective_of_divisorial`)そのもの。 -/
theorem Pf.mk_eq_mk_same_iff {M : Type w} [AddCommMonoid M] (h : IsDivisorial M)
    {x y : M} {r : ℕ+} : Pf.mk x r = Pf.mk y r ↔ x = y := by
  constructor
  · intro hxy
    obtain ⟨k, hk⟩ := Quotient.exact hxy
    exact nsmul_injective_of_divisorial h.1.1 h.1.2.1 h.2
      (n := (k : ℕ) * (r : ℕ)) (Nat.mul_pos k.pos r.pos) hk
  · rintro rfl
    rfl

/-! ## ★2. 2 元を同じ段に揃える -/

section CommonRep

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P)

/-- ★**2 元は同じ段の代表元で書ける** —— 添字圏が filtered だから。 -/
theorem exists_common_rep {A B : C} (f g : HomPf P F A B) :
    ∃ (V : IdxPf P F A B) (φ ψ : V.right.obj.1 ⟶ V.right.obj.2),
      f = HomPf.mk V φ ∧ g = HomPf.mk V ψ := by
  obtain ⟨Z, φ₀, hZ⟩ := HomColim.exists_rep (homFunctorPf P F A B) f
  obtain ⟨W, ψ₀, hW⟩ := HomColim.exists_rep (homFunctorPf P F A B) g
  refine ⟨IsFiltered.max Z W,
    idxTransport P F (IsFiltered.leftToMax Z W) φ₀.down,
    idxTransport P F (IsFiltered.rightToMax Z W) ψ₀.down, ?_, ?_⟩
  · rw [← hZ]
    exact (HomColim.mk_map (homFunctorPf P F A B) (IsFiltered.leftToMax Z W) φ₀).symm
  · rw [← hW]
    exact (HomColim.mk_map (homFunctorPf P F A B) (IsFiltered.rightToMax Z W) ψ₀).symm

/-! ## ★3. 2 つの同値関係の一致 -/

variable {P F}

/-- ★`𝔽_Φ` で一致すれば `𝒞^pf` でも一致する(易しい向き)。 -/
theorem pfToElem_mk_congr {A B : C} (V : IdxPf P F A B)
    {φ ψ : V.right.obj.1 ⟶ V.right.obj.2} (h : P.toElem.map φ = P.toElem.map ψ) :
    (pfToElem P F).map (HomPf.mk V φ) = (pfToElem P F).map (HomPf.mk V ψ) := by
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · show pfBase (HomPf.mk V φ) = pfBase (HomPf.mk V ψ)
    rw [pfBase_mk, pfBase_mk]
    simp only [repBase, repBaseOf]
    rw [show P.Base φ = P.Base ψ from congrArg ElemFrobCat.Hom.base h]
  · show pfDiv (HomPf.mk V φ) = pfDiv (HomPf.mk V ψ)
    rw [pfDiv_mk, pfDiv_mk]
    simp only [repDiv]
    rw [show P.Div φ = P.Div ψ from congrArg ElemFrobCat.Hom.div h]
  · show pfDeg (HomPf.mk V φ) = pfDeg (HomPf.mk V ψ)
    rw [pfDeg_mk, pfDeg_mk]
    exact congrArg ElemFrobCat.Hom.deg h

/-- ★★★**`𝒞^pf` で一致すれば `𝔽_Φ` でも一致する**(中身のある向き)。

★`pfDiv` から `Div` を戻すところで **`Φ` の divisorial 性**が効く。 -/
theorem toElem_eq_of_pfToElem_mk {A B : C} (V : IdxPf P F A B)
    {φ ψ : V.right.obj.1 ⟶ V.right.obj.2}
    (h : (pfToElem P F).map (HomPf.mk V φ) = (pfToElem P F).map (HomPf.mk V ψ)) :
    P.toElem.map φ = P.toElem.map ψ := by
  haveI ha : IsIso (P.Base V.hom.hom.1) := V.hom.property.1.2
  haveI hb : IsIso (P.Base V.hom.hom.2) := V.hom.property.2.1.2
  have hdeg : P.degFr φ = P.degFr ψ := by
    have h3 : pfDeg (HomPf.mk V φ) = pfDeg (HomPf.mk V ψ) := congrArg ElemFrobCat.Hom.deg h
    rw [pfDeg_mk, pfDeg_mk] at h3
    exact h3
  have hbase : P.Base φ = P.Base ψ := by
    have h1 : pfBase (HomPf.mk V φ) = pfBase (HomPf.mk V ψ) := congrArg ElemFrobCat.Hom.base h
    rw [pfBase_mk, pfBase_mk] at h1
    simp only [repBase, repBaseOf] at h1
    have h2 := (cancel_epi (P.Base V.hom.hom.1)).mp h1
    exact (cancel_mono (inv (P.Base V.hom.hom.2))).mp h2
  have hdiv : P.Div φ = P.Div ψ := by
    have h2 : pfDiv (HomPf.mk V φ) = pfDiv (HomPf.mk V ψ) := congrArg ElemFrobCat.Hom.div h
    rw [pfDiv_mk, pfDiv_mk] at h2
    simp only [repDiv] at h2
    exact Φ.map_injective _
      ((Pf.mk_eq_mk_same_iff (P.divisorial _)).mp h2)
  exact ElemFrobCat.Hom.ext hbase hdiv hdeg

/-- ★★★★**2 つの同値関係は一致する** ——
「`Hom^pf` の 2 元がある段で `𝔽_Φ` の像を共有する」ことと
「`𝒞^pf` の pre-Frobenioid 構造で同じ像を持つ」ことは同値。 -/
theorem toQuot_eq_iff_pfToElem {A B : C} (f g : HomPf P F A B) :
    HomColim.toQuot (unTrQuotPf P F A B) f = HomColim.toQuot (unTrQuotPf P F A B) g
      ↔ (pfToElem P F).map f = (pfToElem P F).map g := by
  obtain ⟨V, φ, ψ, hf, hg⟩ := exists_common_rep P F f g
  subst hf
  subst hg
  constructor
  · intro h
    obtain ⟨W, u, hu⟩ :=
      (HomColim.toQuot_eq_iff_same (F := homFunctorPf P F A B) (unTrQuotPf P F A B)
        (i := V) (x := ULift.up φ) (y := ULift.up ψ)).mp h
    have e₁ : HomPf.mk W (idxTransport P F u φ) = HomPf.mk V φ :=
      HomColim.mk_map (homFunctorPf P F A B) u (ULift.up φ)
    have e₂ : HomPf.mk W (idxTransport P F u ψ) = HomPf.mk V ψ :=
      HomColim.mk_map (homFunctorPf P F A B) u (ULift.up ψ)
    rw [← e₁, ← e₂]
    exact pfToElem_mk_congr W hu
  · intro h
    refine (HomColim.toQuot_eq_iff_same (F := homFunctorPf P F A B) (unTrQuotPf P F A B)
      (i := V) (x := ULift.up φ) (y := ULift.up ψ)).mpr ⟨V, 𝟙 V, ?_⟩
    show P.toElem.map (idxTransport P F (𝟙 V) φ) = P.toElem.map (idxTransport P F (𝟙 V) ψ)
    rw [idxTransport_id, idxTransport_id]
    exact toElem_eq_of_pfToElem_mk V h

end CommonRep

/-! ## ★4. 主定理 -/

section Main

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)
  (F₁ : FrobenioidCore (istrPre P Fc)) (F₂ : FrobenioidCore (unTrPre P Fc))

/-- ★★**`(𝒞^pf)^un-tr` の射の集合は、前ファイルの商そのもの**。 -/
noncomputable def unTrPf_equiv_quot (A B : Istr P) :
    HomUnTr (pfPre (istrPre P Fc) F₁) A B
      ≃ Quotient (HomColim.quotKer (unTrQuotPf (istrPre P Fc) F₁ A B)) :=
  Quotient.congr (Equiv.refl _) (fun a b =>
    (toQuot_eq_iff_pfToElem (P := istrPre P Fc) (F := F₁) a b).symm)

/-- ★★★★★★**[FrdI] Proposition 5.5, (ii)**(un-tr の側、原文の形)——

  `Hom_{(𝒞^pf)^un-tr}(A, B) ≃ Hom_{(𝒞^un-tr)^pf}(A, B)`

原文 (FrdI p.105):
> tween the respective sets of morphisms between the images of two given objects of C -/
noncomputable def prop_5_5_ii_untr_cat (A B : Istr P) :
    HomUnTr (pfPre (istrPre P Fc) F₁) A B ≃ HomPf (unTrPre P Fc) F₂ A B :=
  (unTrPf_equiv_quot P Fc F₁ A B).trans (prop_5_5_ii_untr P Fc F₁ F₂ A B)

end Main

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (ii)`(un-tr の側、原文の形)。 -/
def prop_5_5_ii_untr_cat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (ii) — Hom_{(𝒞^pf)^un-tr} ≃ Hom_{(𝒞^un-tr)^pf}",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
