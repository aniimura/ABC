/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54
import ABC3.Found.FrdI.Thm52NatIso

/-!
# [FrdI] Proposition 5.5, (iii)(iv) の実化の側

原文 (FrdI p.104):
> then so is Cpf. Moreover, Cun-tr, Crlf are always of model type. Finally, suppose

★`(iii)` の「`𝒞^un-tr`, `𝒞^rlf` はつねに model 型」のうち、
**`𝒞^un-tr` は済**(`unTr_isOfModelType`)。本ファイルは **`𝒞^rlf` の側**を埋める。

★`𝒞^rlf` は定義からして model Frobenioid(`Proposition 5.3` の第 1 文)なので、
`Theorem 5.2, (ii)`(`model_isOfModelType`、実装済み)を当てるだけである。
★★**そこで実際に効くのは `Hyp` の 4 条**で、そのうち
**`B` が group-like** は無料で出る —— `S·Φ^birat` は**部分群**だからである。

## ★★残る仮定(記録)

`Hyp` の残り 3 条(`Φ^rlf` が divisorial・`𝒟` が totally epimorphic・connected)と
`Skeletal 𝒟` は仮定として受け取る。
★このうち **`Φ^rlf` が divisorial** は `rlf-flat` と同じ性質の宿題である
(`ℝ≥0 ⊗_ℕ Φ` が integral・saturated・characteristic 型であること)。
★`Skeletal 𝒟` は `model_isPreModelType` が要求するもので、
`index.html` に逸脱として開示済み。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {S : Type} [CommSemiring S]

/-! ## ★1. `S·Φ^birat` は group-like -/

/-- ★部分群は group-like。 -/
theorem isGroupLike_addSubgroup {N : Type w} [AddCommGroup N] (T : AddSubgroup N) :
    IsGroupLike ↥T :=
  (isGroupLike_iff _).mpr fun a => ⟨⟨a, -a, add_neg_cancel a, neg_add_cancel a⟩, rfl⟩

/-- ★★`S·Φ^birat` は group-like —— 原文が有理関数の単系に課す条件は**無料で出る**。 -/
theorem sPhiBiratOn_isGroupLike (G : Frobenioid P) (d : D) :
    IsGroupLike ↥(sPhiBiratOn S G d) :=
  isGroupLike_addSubgroup _

/-! ## ★2. `Hyp` と model 型 -/

/-- ★★`(Φ ⊗_S, S·Φ^birat)` の `Hyp` —— **group-like の条だけは無料**。 -/
theorem scModel_hyp (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdiv : (phiScOn S Φ hcharInj).IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D) :
    ModelData.Hyp (scModel S G hiso hfn hcharInj hint hfsmD) where
  divisorial := hdiv
  bmonGroupLike A := sPhiBiratOn_isGroupLike G A
  totEpiD := htot
  connectedD := hconn

/-- ★★★★**[FrdI] Proposition 5.5, (iii) の一節** —— `𝒞^rlf` は**つねに model 型**。

★`S = ℝ≥0` が実化、`S = ℚ≥0` が完全化の側。 -/
theorem scModel_isOfModelType (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdiv : (phiScOn S Φ hcharInj).IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D) (hskel : Skeletal D) :
    haveI := hconn
    IsOfModelType (ModelData.Obj (scModel S G hiso hfn hcharInj hint hfsmD))
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn))
      (ModelData.model_frobenioid
        (scModel_hyp G hiso hfn hcharInj hint hfsmD hdiv htot hconn)) :=
  ModelData.model_isOfModelType _ hskel

/-- ★★**[FrdI] Proposition 5.5, (iv) の実化の場合** ——
`𝒞^rlf` は `(Φ^rlf, ℝ·Φ^birat)` の model Frobenioid **そのもの**。

★★原文は「there is a natural equivalence of categories between `𝒞^rlf` and the model
Frobenioid associated to `Φ^rlf`, `ℝ·Φ^birat`」と述べるが、
**我々は `Proposition 5.3` の第 1 文(定義)をそのまま採ったので、これは等式である**。
★これは逃げではない —— 原文自身が第 1 文で `𝒞^rlf` を model Frobenioid と**定義**しており、
(iv) の実化の場合はその定義の再確認だからである。
★★内容が残るのは **`𝒞^pf` の場合**(鎖 `prop55` の `p55iv-pf`)である。 -/
theorem crlf_eq_model (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    (hfsmD : IsOfFSMType D) :
    Crlf G hiso hfn hcharInj hint hfsmD
      = ModelData.Obj (scModel NNReal G hiso hfn hcharInj hint hfsmD) := rfl

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.5, (iii)` の「`𝒞^rlf` はつねに model 型」の条
(★**条つき**: `Φ^rlf` の divisorial 性と `Skeletal 𝒟` を仮定に置いている)。 -/
def scModel_isOfModelType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (iii) — 𝒞^rlf はつねに model 型",
    sectionId := "frdi-prop-5-5" }

/-- ★locator —— `Proposition 5.5, (iv)` の実化の場合。 -/
def crlf_eq_model.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (iv) — 𝒞^rlf は (Φ^rlf, ℝ·Φ^birat) の model Frobenioid",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
