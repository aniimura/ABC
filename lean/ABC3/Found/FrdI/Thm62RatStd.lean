/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Birat
import ABC3.Found.FrdI.Prop48Cond3
import ABC3.Found.FrdI.Prop48Cpt
import ABC3.Found.FrdI.Thm52NatIso
import ABC3.Found.FrdI.Thm52Frob

/-!
# [FrdI] Theorem 6.2, (iii) —— rationally standard 型の組み立て

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type.

## ★★何が残っていたか

`Definition 4.5, (iii)` の `IsOfRationallyStandardType` は 4 条:

| 条 | 状況 |
|---|---|
| (a) birationally Frobenius-normalized 型 | ★model では**自動**(`model_isOfBiratFrobNormalizedType`) |
| (a) rational 型 | `isOfRationalType_of_divB`(`Prop53Birat.lean`、済) |
| (a) standard 型 | `model_isOfStandardType`(`Thm52Frob.lean`、6 条中 3 条は自動) |
| (b) `(𝒞^un-tr)^birat` の Frobenius-compact 対象 | ★**本ファイル** |

★(b) の中身は `birat_isFrobeniusCompact_of_ne_zero`(`Prop48Cond3.lean`、在庫)であり、
入力は「`Div : 𝒪^▷ ↠ Φ`」「`(𝒞^un-tr)^birat` が Frobenius-normalized」
「`Φ` が non-dilating」「`Φ` に 0 でない元が 1 つある」の 4 つである。

★★**`Div : 𝒪^▷ ↠ Φ` は `𝒞^istr` の側で持っていれば足りる**
(`unTr_divSurj_of_istr`)—— `𝒞^un-tr` は射を `𝔽_Φ` の像で割った商なので、
`Base` も `degFr` も `Div` も**代表元のものがそのまま降りる**(3 つとも `rfl` で済んだ)。

★★「`(𝒞^un-tr)^birat` が Frobenius-normalized」は
`unTr_isOfModelType` の第 2 成分がそのまま与える(在庫)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `unTr_divSurj_of_istr` | ★`Div : 𝒪^▷ ↠ Φ` が `𝒞^istr → 𝒞^un-tr` で降りる |
| `unTrBiratCompact_of_ne_zero` | ★★`(𝒞^un-tr)^birat` の Frobenius-compact 対象 |
| `ModelData.model_isOfRationallyStandardType` | ★★★★★★**model Frobenioid は rationally standard 型** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)

/-! ## ★1. `Div : 𝒪^▷ ↠ Φ` は `𝒞^un-tr` へ降りる -/

/-- ★★**`Div : 𝒪^▷ ↠ Φ` は `𝒞^istr` から `𝒞^un-tr` へ移る**。

★`𝒞^un-tr` は射を `𝔽_Φ` の像で割った商なので、
`Base`・`degFr`・`Div` はどれも**代表元のものがそのまま降りる**。 -/
theorem unTr_divSurj_of_istr
    (hdivS : ∀ (Y : Istr P) (a : Φ.val ((istrPre P Fc).toElem.obj Y).base),
      ∃ u : OTri (istrPre P Fc) Y, (istrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a)
    (Y : UnTr P) (a : Φ.val ((unTrPre P Fc).toElem.obj Y).base) :
    ∃ u : OTri (unTrPre P Fc) Y, (unTrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a := by
  obtain ⟨u, hu⟩ := hdivS (show Istr P from Y) a
  refine ⟨⟨(istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _), ?_, ?_⟩, hu⟩
  · show (unTrPre P Fc).Base ((istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _))
      = (unTrPre P Fc).Base (𝟙 Y)
    exact u.2.1
  · show (unTrPre P Fc).degFr ((istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _)) = 1
    exact u.2.2

/-! ## ★2. `(𝒞^un-tr)^birat` の Frobenius-compact 対象 -/

/-- ★★★**`Definition 4.5, (iii), (b)`** —— `(𝒞^un-tr)^birat` の Frobenius-compact 対象。

★中身は `birat_isFrobeniusCompact_of_ne_zero`(在庫)を `𝒞^un-tr` に当てるだけ。 -/
theorem unTrBiratCompact_of_ne_zero (G : Frobenioid P)
    (hdivS : ∀ (Y : Istr P) (a : Φ.val ((istrPre P Fc).toElem.obj Y).base),
      ∃ u : OTri (istrPre P Fc) Y, (istrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a)
    (hfn : ∀ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G),
      IsFrobeniusNormalized (unTrBiratPre P Fc G) X)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (A : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G))
    (x₀ : Φ.val ((unTrPre P Fc).toElem.obj
      (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G) A)).base) (hx₀ : x₀ ≠ 0) :
    ∃ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G),
      IsFrobeniusCompact (unTrBiratPre P Fc G) X :=
  ⟨A, birat_isFrobeniusCompact_of_ne_zero (unTr_divSurj_of_istr P Fc hdivS) hfn hndOn A x₀ hx₀⟩

/-! ## ★3. model Frobenioid は rationally standard 型 -/

namespace ModelData

variable {M : ModelData.{v, u, w} D}

/-- ★★★★★★**[FrdI] Theorem 6.2, (iii)** —— model Frobenioid が **rationally standard 型**。

★4 条のうち **(a) の第 1 条は model では自動**である
(`model_isOfBiratFrobNormalizedType`、`model_ratStandardType_iff` が吸収済み)。

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type. -/
theorem model_isOfRationallyStandardType (h : Hyp M)
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
    (hdivS : ∀ (Y : Istr (modelPre h))
        (a : M.phi.val ((istrPre (modelPre h) (model_frobenioid h).core).toElem.obj Y).base),
      ∃ u : OTri (istrPre (modelPre h) (model_frobenioid h).core) Y,
        (istrPre (modelPre h) (model_frobenioid h).core).Div ((u : End Y) : Y ⟶ Y) = a)
    (A : BiratCat (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)))
    (x₀ : M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
      (biratDown _ _ A)).base) (hx₀ : x₀ ≠ 0) :
    haveI := h.connectedD
    IsOfRationallyStandardType (modelPre h) (model_frobenioid h) ι :=
  haveI := h.connectedD
  (model_ratStandardType_iff h ι).mpr
    ⟨isOfRationalType_of_divB R ι hsp,
     model_isOfStandardType h (model_frobenioid h).core hfsmff hnd hgl,
     unTrBiratCompact_of_ne_zero (modelPre h) (model_frobenioid h).core (model_frobenioid h)
       hdivS
       (fun X => (unTr_isOfModelType (model_frobenioid h).core (model_frobenioid h)).2 X)
       hnd A x₀ hx₀⟩

end ModelData

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Theorem 6.2, (iii)` の rationally standard 型。 -/
def ModelData.model_isOfRationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — model Frobenioid は rationally standard 型",
    sectionId := "frdi-thm-6-2" }

/-- ★★★locator —— `Definition 4.5, (iii), (b)` の Frobenius-compact 対象。 -/
def unTrBiratCompact_of_ne_zero.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — (𝒞^un-tr)^birat の Frobenius-compact 対象",
    sectionId := "frdi-thm-6-2" }

end ABC3.Found.FrdI
