/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Model
import ABC3.Found.FrdI.Thm62RatStd
import ABC3.Found.FrdI.Prop55Std

/-!
# [FrdI] Theorem 6.4, (i) —— 5 つの圏の型(`isotropic` の側)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> and rationally standard type, but not of group-like type; D is Frobenius-slim

## ★★原文の主張を 3 つに割る

原文 `Theorem 6.4, (i)` は `C, C^pf, C^rlf, C^un-tr, (C^pf)^un-tr` について

1. **isotropic 型**
2. **rationally standard 型**
3. **group-like 型でない**

の 3 つを言う。★本ファイルは **1 を model Frobenioid 一般で閉じ**、
`Example 6.3`(算術 Frobenioid)へ当てる。

## ★★★測定(2026-08-25)—— 在庫の位置

| 要るもの | 在庫 |
|---|---|
| `C` が isotropic | `ModelData.model_isotropicType`(在庫) |
| `C^pf` が isotropic | `pfRoot_isOfIsotropicType`(在庫) |
| `C^istr` が isotropic | `istr_isOfIsotropicType`(在庫) |
| `C^un-tr` が **Frobenius**-isotropic | `unTr_frobIsotropicType`(在庫) |
| `C^un-tr` が isotropic | ★**無い**(2026-08-25 実測、`exact?` 空振り) |
| rationally standard(一般形) | `ModelData.model_isOfRationallyStandardType_of_primary`(在庫、**逸脱なし**) |

★★`Example 6.3` の側は **`ex63Frobenioid` まで組み上がっている**
(`Found/Divisor/ArithFrobenioid.lean` ほか 17 ファイル、4,291 行、`sorry` 0)。
★★★`Skeleton/Divisor/ArithDivisor.lean` の 11 個の `sorry` は
**実装済みの `Found` を複製した古い写し**であり、この節点を止めているものではない。

## ★★★★残っているもの(`thm64-i-types`)

`model_isOfRationallyStandardType_of_primary` の仮引数のうち、
**算術の側で作らねばならない入力**が残る:

* `hsp` —— どの素点 `p` にも `a, b ∈ Φ(L)` と `y ∈ B(L)` があって
  `toGp a − toGp b = ±Div_B(y)`、`p ∈ Supp a`、`p ∉ Supp b`。
  ★有限素点は `exists_ordFin_eq_one`(在庫)で一様化元を取ればよい。
  ★★無限素点は `y = 2` のように `|y|_w ≠ 1` の元を取る段が要る。
* `hprim` / `hx₀mem` —— 準素元が `Φ^birat` に入ること。
* `hgl` —— group-like のときのコンパクト対象(`Φ ≠ 0` なので前件が偽の見込み)。

★`hnd`(non-dilating)は **`arithDatumGalois_isNonDilatingOn` で済**、
`hfsmff` は `finSubOp_isOfFSMType` から出る。
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI ABC3.Meta

/-! ## ★1. model Frobenioid 一般 —— isotropic の 4 条 -/

namespace ModelIso

open ABC3.Found.FrdI.ModelData

universe v u w

variable {D : Type u} [Category.{v} D] {M : ModelData.{v, u, w} D}

/-- ★★★★★**model Frobenioid とその派生圏の isotropic 性**。

★原文 `Theorem 6.4, (i)`(および `Theorem 6.2, (iii)`)が
`C, C^pf, C^un-tr, (C^pf)^un-tr, …` について主張する isotropic 型。

## ★★★★見立て違いの訂正(2026-08-25、**9 回目**)

前版はここに「`C^un-tr` は**Frobenius**-isotropic までしか在庫が無い」と書いていたが、
**誤り**だった —— `Prop33UnTr.lean` の

    unTr_isotropic P Fc B : IsIsotropic (unTrPre P Fc) B

が**すべての対象**について成り立つので、`IsOfIsotropicType (unTrPre P Fc)` は
`fun B => unTr_isotropic P Fc B` の 1 行である。
★実際 `unTr_frobIsotropicType` 自身がこの補題を使って書かれていた。
★★`exact?` が空振りしたのは「`IsOfIsotropicType` を展開しないと当たらない」形だから
であって、在庫が無いからではない。

★教訓は毎回同じ: **「無い」と書く前に、使っている補題の中身を読む。** -/
theorem model_isotropic_family (h : M.Hyp) :
    IsOfIsotropicType (modelPre h)
      ∧ IsOfFrobeniusIsotropicType (modelPre h)
      ∧ IsOfIsotropicType (pfRootPre (modelPre h) (model_frobenioid h).core)
      ∧ IsOfIsotropicType (istrPre (modelPre h) (model_frobenioid h).core)
      ∧ IsOfFrobeniusIsotropicType (unTrPre (modelPre h) (model_frobenioid h).core)
      ∧ IsOfIsotropicType (unTrPre (modelPre h) (model_frobenioid h).core)
      ∧ IsOfIsotropicType (unTrPre (pfRootPre (modelPre h) (model_frobenioid h).core)
          (pfRootCore (F := (model_frobenioid h).core)
            (isOfFrobeniusIsotropicType_of_isotropic (model_isotropicType h)))) := by
  haveI := h.connectedD
  have hiso : IsOfIsotropicType (modelPre h) := model_isotropicType h
  have hfiso : IsOfFrobeniusIsotropicType (modelPre h) :=
    isOfFrobeniusIsotropicType_of_isotropic hiso
  exact ⟨hiso, hfiso, pfRoot_isOfIsotropicType hfiso, istr_isOfIsotropicType,
    unTr_frobIsotropicType _, fun A => unTr_isotropic _ _ A, fun A => unTr_isotropic _ _ A⟩

end ModelIso

/-! ## ★2. `Example 6.3` への当て込み -/

section Ex63

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

variable (F Kbar) in
/-- ★`Example 6.3` の model Frobenioid の仮説束。 -/
theorem ex63Hyp : (ex63ModelData (F := F) (Kbar := Kbar)).Hyp :=
  (arithDatumGalois F Kbar).modelHyp finSubOp_isOfFSMType (bmonGalois F Kbar)
    (fun A => divBGalois A) (fun f x => divBGalois_nat f x)
    (bmonGalois_isGroupLike F Kbar) finSubOp_totallyEpimorphic finSubOp_isConnected

variable (F Kbar) in
/-- ★★★★★**[FrdI] Theorem 6.4, (i)** —— 算術 Frobenioid とその派生圏の isotropic 性。

原文 (FrdI p.114):
> these indices]. Then the Frobenioids C, Cpf, Crlf, Cun-tr, (Cpf)un-tr are of isotropic

★★rationally standard 型の側は `model_isOfRationallyStandardType_of_primary`
(在庫、`sorry` 無し)に算術の入力(`hsp` / `hprim` / `hgl`)を渡す段が残る。 -/
theorem ex63_isotropic_family :
    IsOfIsotropicType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ IsOfFrobeniusIsotropicType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ IsOfIsotropicType (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType (istrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfFrobeniusIsotropicType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType
          (unTrPre (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
              (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
            (pfRootCore (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
              (isOfFrobeniusIsotropicType_of_isotropic
                (ModelData.model_isotropicType (ex63Hyp F Kbar))))) :=
  ModelIso.model_isotropic_family (ex63Hyp F Kbar)

end Ex63

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def ModelIso.model_isotropic_family.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — C, C^pf, C^un-tr … の isotropic 型(model 一般)",
    sectionId := "frdi-thm-6-4" }

def ModelIso.model_isotropic_family.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ModelData.model_isotropicType"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.model_isotropicType") 114,
    .citation "[ABC3]" "pfRoot_isOfIsotropicType"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_isOfIsotropicType") 114,
    .citation "[ABC3]" "istr_isOfIsotropicType"
      (.inProject "ABC3" "ABC3.Found.FrdI.istr_isOfIsotropicType") 114,
    .citation "[ABC3]" "unTr_frobIsotropicType"
      (.inProject "ABC3" "ABC3.Found.FrdI.unTr_frobIsotropicType") 114,
    .citation "[ABC3]" "C^un-tr が(Frobenius でない)isotropic であること"
      (.absent "2026-08-25 実測。unTrPre に対する IsOfIsotropicType は在庫に無い(exact? 空振り)") 114 ]

def ex63_isotropic_family.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 Frobenioid の isotropic 型",
    sectionId := "frdi-thm-6-4" }

def ex63_isotropic_family.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ex63ModelData / ex63Frobenioid(Example 6.3、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ex63Frobenioid") 114,
    .citation "[ABC3]" "model_isotropic_family"
      (.inProject "ABC3" "ABC3.Found.Divisor.ModelIso.model_isotropic_family") 114,
    .citation "[ABC3]" "rationally standard 型の一般形(算術の入力待ち)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.model_isOfRationallyStandardType_of_primary") 114,
    .implicitStep
      "★原文は 5 つの圏の型を 1 文で並べる。rationally standard の側は hsp / hprim / hgl の算術入力が残る" 114 ]

def ex63Hyp.src : Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — model Frobenioid の仮説束",
    sectionId := "frdi-example-6-3" }

def ex63Hyp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "arithDatumGalois(算術の ArithDatum)"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDatumGalois") 113 ]

end ABC3.Found.Divisor
