/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.TypeTransport
import ABC3.Found.FrdI.Cor54Rigid
import ABC3.Found.FrdI.Cor57Model

/-!
# [FrdI] Proposition 5.5, (iii) —— standard 型が `𝒞^un-tr` へ移る

原文 (FrdI p.105):
> the pull-back morphisms of Cun-tr, Crlf are precisely the linear isometries [cf. Proposition 1.4,

★★原文 (iii) の最後の文は
「`𝒞` が not group-like で standard(resp. rationally standard)型なら
`𝒞^un-tr`・`𝒞^rlf` もそう」である。

`IsOfStandardType` は 6 条:

| 条 | `𝒞^un-tr` での根拠 |
|---|---|
| quasi-isotropic | `unTr_isotropic`(すべての対象が isotropic) |
| Frobenius-isotropic | 恒等射が Frobenius 型 ＋ isotropic |
| **Frobenius-normalized** | ★本ファイル(`unTr_frobNormalizedType`) |
| `𝒟` が FSMFF | `𝒞` と同じ `𝒟` |
| `Φ` が non-dilating | `𝒞` と同じ `Φ` |
| group-like ⟹ compact 対象 | ★未(`𝒞` が not group-like なら前件が偽になることを要する) |

## ★★Frobenius-normalized の筋

★`𝒞^un-tr` は **model 型**(`unTr_isOfModelType`)なので、
`Theorem 5.2, (iv)` の圏同値 `unTr_modelFrobenioid` で model Frobenioid と結べる。
★model Frobenioid では Frobenius-normalized は**無条件**である
(`ModelData.model_frobNormalizedType` —— 射が 4 成分の明示的な組だから)。

★★★そこで**逆向きに引き戻す**。要る道具は `TypeTransport.lean` の

  `isFrobeniusNormalized_map_of_toElem_iso`

で、仮定は「充満忠実 ＋ `𝔽_Φ` への関手と 1-可換」だけである。
1-可換性は在庫の `modelTypeEquiv_comp_toElem`(`Cor54Rigid.lean`)がそのまま与える。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u4 v4

/-! ## ★1. `RatFnData` からの `Hyp` と、1-可換性の逆向き -/

section Generic

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-- ★★`RatFnData` からは `ModelData.Hyp` が**無料で出る**。

★`B` が group-like であることは `RatFnData` の `bneg` がそのまま与え、
残る 3 条(`Φ` が divisorial・`𝒟` が totally epimorphic・`𝒟` が connected)は
`PreFrobenioid` のフィールドそのものである。 -/
theorem RatFnData.hyp {G : Frobenioid P} (R : RatFnData P G) : ModelData.Hyp R.model where
  divisorial := P.divisorial
  bmonGroupLike A := (isGroupLike_iff _).mpr fun a =>
    ⟨⟨a, R.bneg A a, (add_comm a (R.bneg A a)).trans (R.bneg_add A a), R.bneg_add A a⟩, rfl⟩
  totEpiD := P.totEpiD
  connectedD := P.connectedD

end Generic

section InverseIso

variable {D₀ : Type u} [Category.{v} D₀] {Φ₀ : MonoidOn.{v, u, w} D₀}
  {C₁ : Type u2} [Category.{v2} C₁] {C₂ : Type u4} [Category.{v4} C₂]
  {P₁ : PreFrobenioid C₁ Φ₀} {P₂ : PreFrobenioid C₂ Φ₀}

/-- ★★**`𝔽_Φ` への関手との 1-可換性は圏同値の逆向きへ移る**。

★`η` を `e.inverse` で左から whisker し、counit で `𝟭` に潰すだけ。 -/
noncomputable def inverse_comp_toElem_iso (e : C₁ ≌ C₂)
    (η : e.functor ⋙ P₂.toElem ≅ P₁.toElem) :
    e.inverse ⋙ P₁.toElem ≅ P₂.toElem :=
  Functor.isoWhiskerLeft e.inverse η.symm ≪≫ (Functor.associator _ _ _).symm
    ≪≫ Functor.isoWhiskerRight e.counitIso P₂.toElem ≪≫ Functor.leftUnitor _

end InverseIso

/-! ## ★2. `𝒞^un-tr` は Frobenius-normalized 型 -/

section UnTr

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D]

set_option maxHeartbeats 800000 in
/-- ★★★★**[FrdI] Proposition 5.5, (iii) の一条** —— `𝒞^un-tr` は **Frobenius-normalized 型**。

★★筋は 3 段:
1. `𝒞^un-tr` は model 型(`unTr_isOfModelType`)なので model Frobenioid と圏同値。
2. model Frobenioid では Frobenius-normalized は**無条件**
   (`ModelData.model_frobNormalizedType`)。
3. 1-可換性(`modelTypeEquiv_comp_toElem`)を**逆向き**にして引き戻す
   (`isFrobeniusNormalized_map_of_toElem_iso`)。 -/
theorem unTr_frobNormalizedType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    IsOfFrobeniusNormalizedType (unTrPre P Fc) := by
  intro A
  have h : ModelData.Hyp (unTr_ratFnData Fc G hint hfsmD).model :=
    RatFnData.hyp (unTr_ratFnData Fc G hint hfsmD)
  have η : (unTr_modelFrobenioid Fc G hint hfsmD).functor
      ⋙ (ModelData.modelPre h).toElem ≅ (unTrPre P Fc).toElem :=
    modelTypeEquiv_comp_toElem (unTr_ratFnData Fc G hint hfsmD)
      (unTr_isotropic P Fc) (unTr_isOfModelType Fc G)
  have η' := inverse_comp_toElem_iso (unTr_modelFrobenioid Fc G hint hfsmD) η
  have hback : IsFrobeniusNormalized (unTrPre P Fc)
      ((unTr_modelFrobenioid Fc G hint hfsmD).inverse.obj
        ((unTr_modelFrobenioid Fc G hint hfsmD).functor.obj A)) :=
    isFrobeniusNormalized_map_of_toElem_iso _ η' _
      (ModelData.model_frobNormalizedType h _)
  exact isFrobeniusNormalized_of_iso
    ((unTr_modelFrobenioid Fc G hint hfsmD).unitIso.app A).symm hback

/-- ★★`𝒞^un-tr` は **Frobenius-isotropic 型**(すべての対象が isotropic だから恒等射でよい)。 -/
theorem unTr_frobIsotropicType (Fc : FrobenioidCore P) :
    IsOfFrobeniusIsotropicType (unTrPre P Fc) :=
  fun A => ⟨A, 𝟙 A, isFrobeniusType_of_isIso (unTrPre P Fc) (𝟙 A), unTr_isotropic P Fc A⟩

/-! ## ★3. standard 型が `𝒞^un-tr` へ移る

★★`IsOfStandardType` の 6 条のうち 5 条は上と在庫で埋まり、残る
`groupLikeCompact` は**前件が偽**になる ——
原文の「suppose that `𝒞`, hence also `𝒞^un-tr`, `𝒞^rlf`, are not of group-like type」が
それである。

★`𝒞^un-tr` の対象は `𝒞^istr` の対象そのもの(`UnTr P := Istr P`)で、
`𝔽_Φ` への関手も `A ↦ P.toElem.obj A.obj` なので、
**group-like 性は `𝒞` のそれと同じ条件**である(`IsGroupLikeObj` は `Φ` と底だけで決まる)。
★したがって `𝒞` が group-like 型なら `𝒞^un-tr` もそう(`unTr_isOfGroupLikeType`)。
★★逆(原文の「hence also」)は `𝒞` の各対象を isotropic 包へ送る段が要るので、
ここでは **`¬ IsOfGroupLikeType (unTrPre P Fc)` を仮引数で受ける**。 -/

omit [IsConnected D] in
/-- ★**`𝒞` が group-like 型なら `𝒞^un-tr` も** —— 対象も底も同じだから。 -/
theorem unTr_isOfGroupLikeType (Fc : FrobenioidCore P) (h : IsOfGroupLikeType P) :
    IsOfGroupLikeType (unTrPre P Fc) :=
  fun A => h (show Istr P from A).obj

set_option maxHeartbeats 800000 in
/-- ★★★★★**[FrdI] Proposition 5.5, (iii)** —— **standard 型は `𝒞^un-tr` へ移る**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★6 条の内訳:
| 条 | 根拠 |
|---|---|
| quasi-isotropic | `unTr_isotropic`(isotropic 型 ⟹ quasi-isotropic) |
| Frobenius-isotropic | `unTr_frobIsotropicType` |
| group-like ⟹ compact 対象 | ★**前件が偽**(`hngl`) |
| Frobenius-normalized | `unTr_frobNormalizedType` |
| `𝒟` が FSMFF | `𝒞` と同じ `𝒟` |
| `Φ` が non-dilating | `𝒞` と同じ `Φ` | -/
theorem unTr_standardType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D)
    (F' : FrobenioidCore (unTrPre P Fc))
    (hngl : ¬ IsOfGroupLikeType (unTrPre P Fc))
    (hstd : IsOfStandardType D C P Fc) :
    IsOfStandardType D (UnTr P) (unTrPre P Fc) F' where
  quasiIsotropic :=
    isOfQuasiIsotropicType_of_isOfIsotropicType (unTrPre P Fc) F' (unTr_isotropic P Fc)
  frobIsotropic := unTr_frobIsotropicType Fc
  groupLikeCompact hg := absurd hg hngl
  frobNormalized := unTr_frobNormalizedType Fc G hint hfsmD
  baseFSMFF := hstd.baseFSMFF
  phiNonDilating := hstd.phiNonDilating

end UnTr

/-! ## ★4. standard 型が `𝒞^rlf` へ移る

★★`𝒞^rlf` は **model Frobenioid そのもの**(`ScModelObj S G … = ModelData.Obj (scModel …)`)
なので、`ModelData.model_isOfStandardType` が**直接当たる** —— 移送は要らない。
★残る入力は 3 つ:
* `hfsmff` —— `𝒟` が FSMFF(`𝒞` と同じ `𝒟` なので `𝒞` の standard 性から)
* `hnd` —— **係数拡大した `Φ` が non-dilating**(★実化が non-dilating を保つこと。未)
* `hngl` —— `𝒞^rlf` が not group-like(原文の「hence also」。未)
-/

section Rlf

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {S : Type} [CommSemiring S]

set_option maxHeartbeats 800000 in
/-- ★★★★★**[FrdI] Proposition 5.5, (iii)** —— **standard 型は `𝒞^rlf` へ移る**。

★`𝒞^rlf` は model Frobenioid そのものなので、`model_isOfStandardType` の直接適用である
(quasi-isotropic・Frobenius-isotropic・Frobenius-normalized の 3 条は
model Frobenioid では**無条件**)。 -/
theorem scModel_standardType (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdiv : (scModel S G hiso hfn hcharInj hint hfsmD).phi.IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D)
    (F' : FrobenioidCore (ModelData.modelPre
      (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)))
    (hfsmff : IsOfFSMFFType D)
    (hnd : (scModel S G hiso hfn hcharInj hint hfsmD).phi.IsNonDilatingOn)
    (hngl : ¬ IsOfGroupLikeType (ModelData.modelPre
      (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn))) :
    IsOfStandardType D (ScModelObj S G hiso hfn hcharInj hint hfsmD)
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)) F' :=
  ModelData.model_isOfStandardType _ F' hfsmff hnd (fun hg => absurd hg hngl)

end Rlf

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Proposition 5.5, (iii)` の「`𝒞^un-tr` は Frobenius-normalized 型」の条。 -/
def unTr_frobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr は Frobenius-normalized 型",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★locator —— `Proposition 5.5, (iii)` の「standard 型が `𝒞^un-tr` へ移る」の条
(★**条つき**: `𝒞^un-tr` が not group-like であることを仮引数で受けている)。 -/
def unTr_standardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — standard 型は 𝒞^un-tr へ移る",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
