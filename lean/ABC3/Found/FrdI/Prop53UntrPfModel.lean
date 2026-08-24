/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfArbFull
import ABC3.Found.FrdI.Prop53UntrBirat

/-!
# [FrdI] Proposition 5.3 の第 2 文 —— `(𝒞^un-tr)^pf` は model 型

原文 (FrdI p.103):
> tively, Φpf) and the rational function monoid Φbirat (respectively, Q · Φbirat =

★★`Proposition 5.3` の第 2 文は 2 つのことを言う:

1. **`𝒞^un-tr`** は model 型で、`(Φ, Φ^birat)` の model Frobenioid である
   —— 済(`unTr_isOfModelType` / `unTr_modelFrobenioid`)。
2. **`(𝒞^un-tr)^pf`** は model 型で、`(Φ^pf, ℚ·Φ^birat)` の model Frobenioid である
   —— **本ファイル**。

## ★なぜ今まで書けなかったか(記録)

一般形 `pfRoot_isOfModelType` は仮定
**`hftr : ∀ X, IsFrobeniusTrivial P X`**(`𝒞` が Frobenius-trivial 型)を
引きずっており、これは `𝒞^un-tr` では成り立たない
(`𝔽_Φ` で見れば `n·α = α`、すなわち全対象の零因子類が `0`)。

★`Prop55PfArb.lean` が原文の 1 文

> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

を埋めて `hftr` を落としたので、**`𝒞^un-tr` へそのまま当たる**ようになった。

★★仮引数として残るのは `hint`(`Φ` が整域単系)と `hfsmD`(`𝒟` が FSM 型)だけで、
どちらも `Proposition 5.3` の設定の常備仮定である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UnTrPf

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

variable (P) in
/-- ★`(𝒞^un-tr)^pf` の Frobenioid 構造(`𝒞^un-tr` は isotropic なので
`Proposition 3.2` がそのまま当たる)。 -/
theorem untrPf_frobenioid (Fc : FrobenioidCore P) (G : Frobenioid P) :
    Frobenioid (pfRootPre (unTrPre P Fc) (unTr_frobenioidCore P Fc)) :=
  pfRoot_frobenioid (F := unTr_frobenioidCore P Fc) (unTr_frobIsotropicType Fc)
    (unTr_frobenioid P Fc G)

variable (P) in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の第 2 文(`(𝒞^un-tr)^pf` の場合)** ——
**`(𝒞^un-tr)^pf` は model 型**。

原文 (FrdI p.103):
> tively, Φpf) and the rational function monoid Φbirat (respectively, Q · Φbirat =

★中身は `pfRoot_isOfModelType_of_arb`(`hftr` なしの版)を
`𝒞 := 𝒞^un-tr` に当てるだけである。4 つの仮定はすべて在庫:

| 仮定 | 在庫 |
|---|---|
| Frobenius-isotropic 型 | `unTr_frobIsotropicType` |
| isotropic | `unTr_isotropic` |
| Frobenius-normalized 型 | `unTr_frobNormalizedType` |
| unit-trivial 型 | `unTr_unitTrivial` | -/
theorem untrPf_isOfModelType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    IsOfModelType (PfRootObj (unTrPre P Fc) (unTr_frobenioidCore P Fc))
      (pfRootPre (unTrPre P Fc) (unTr_frobenioidCore P Fc))
      (untrPf_frobenioid P Fc G) :=
  pfRoot_isOfModelType_of_arb (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
    (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc) _

variable (P) in
/-- ★★★★★★**`(𝒞^un-tr)^pf` の有理関数の単系の interface**。 -/
noncomputable def untrPf_ratFnData (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    RatFnData (pfRootPre (unTrPre P Fc) (unTr_frobenioidCore P Fc)) (untrPf_frobenioid P Fc G) :=
  pfRoot_ratFnData_of_arb (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
    (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc) _ hfsmD

variable (P) in
set_option maxHeartbeats 800000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の第 2 文(`(𝒞^un-tr)^pf` の場合)** ——
**`(𝒞^un-tr)^pf` は model Frobenioid と圏同値**。 -/
noncomputable def untrPf_modelFrobenioid (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    PfRootObj (unTrPre P Fc) (unTr_frobenioidCore P Fc)
      ≌ ModelData.Obj (untrPf_ratFnData P Fc G hint hfsmD).model :=
  pfRoot_modelFrobenioid_of_arb (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
    (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc) _ hfsmD

variable (P) in
set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の第 2 文(`(𝒞^un-tr)^pf` の場合)の単系の側** ——
`(𝒞^un-tr)^pf` の有理関数の単系は各 `d ∈ Ob(𝒟)` で
**`ℚ·(Φ^birat of 𝒞^un-tr)(d)`** である。

原文 (FrdI p.103):
> tively, Φpf) and the rational function monoid Φbirat (respectively, Q · Φbirat =

★★`𝒞^un-tr` は unit-trivial なので `pfRoot_ratFnData_bmon_val_full`(条なし)が
そのまま当たる。★`Φ^birat` を `𝒞` のものへ書き換えた版は
`untrPf_ratFnData_bmon_val_of_C`(下)である。 -/
theorem untrPf_ratFnData_bmon_val (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) (d : D) :
    (untrPf_ratFnData P Fc G hint hfsmD).bmon.val d
      = ↥(qPhiBiratOn (unTrPre P Fc) (unTr_frobenioid P Fc G) d) :=
  pfRoot_ratFnData_bmon_val_full (unTr_frobIsotropicType Fc) (unTr_isotropic P Fc)
    (unTr_frobNormalizedType Fc G hint hfsmD) (unTr_unitTrivial P Fc)
    (unTr_frobenioid P Fc G) (untrPf_frobenioid P Fc G) hfsmD d

variable (P) in
/-- ★`ℚ·Φ^birat` も `𝒞` と `𝒞^un-tr` で一致する(`phiBiratOn_unTr_eq` の飽和版)。 -/
theorem qPhiBiratOn_unTr_eq (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (d : D) :
    qPhiBiratOn (unTrPre P Fc) (unTr_frobenioid P Fc G) d = qPhiBiratOn P G d := by
  have himg : phiBiratPfImage (unTrPre P Fc) (unTr_frobenioid P Fc G) d
      = phiBiratPfImage P G d :=
    congrArg (AddSubgroup.map (gpMap _ (Pf.of (M := Φ.val d))))
      (phiBiratOn_unTr_eq Fc G hiso d)
  ext x
  show (∃ k : ℕ+, ((k : ℕ+) : ℕ) • x
      ∈ phiBiratPfImage (unTrPre P Fc) (unTr_frobenioid P Fc G) d)
    ↔ (∃ k : ℕ+, ((k : ℕ+) : ℕ) • x ∈ phiBiratPfImage P G d)
  rw [himg]

variable (P) in
set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.3 の第 2 文(`(𝒞^un-tr)^pf` の場合)** ——
`(𝒞^un-tr)^pf` の有理関数の単系は各 `d ∈ Ob(𝒟)` で
**`ℚ·Φ^birat(d)`(`𝒞` の `Φ^birat`)** である。

原文 (FrdI p.103):
> tively, Φpf) and the rational function monoid Φbirat (respectively, Q · Φbirat =

★`untrPf_ratFnData_bmon_val` が `𝒞^un-tr` の `Φ^birat` で述べたものを、
`qPhiBiratOn_unTr_eq` で **`𝒞` の `Φ^birat`** へ書き換えたものである。 -/
theorem untrPf_ratFnData_bmon_val_of_C (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) (d : D) :
    (untrPf_ratFnData P Fc G hint hfsmD).bmon.val d = ↥(qPhiBiratOn P G d) := by
  rw [untrPf_ratFnData_bmon_val P Fc G hint hfsmD d, qPhiBiratOn_unTr_eq P Fc G hiso d]

variable (P) in
/-- ★★★★★**[FrdI] Proposition 5.3 の第 2 文(`𝒞^un-tr` の場合)の単系の側** ——
`𝒞^un-tr` の有理関数の単系は各 `d ∈ Ob(𝒟)` で **`Φ^birat(d)`(`𝒞` のもの)**。 -/
theorem unTr_ratFnData_bmon_val_of_C (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) (d : D) :
    (unTr_ratFnData Fc G hint hfsmD).bmon.val d = ↥(phiBiratOn G d) := by
  show ↥(phiBiratOn (unTr_frobenioid P Fc G) d) = ↥(phiBiratOn G d)
  rw [phiBiratOn_unTr_eq Fc G hiso d]

end UnTrPf

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.3` の第 2 文の後半
「`(𝒞^un-tr)^pf` は `(Φ^pf, ℚ·Φ^birat)` の model Frobenioid」の**圏の側**。 -/
def untrPf_modelFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — (𝒞^un-tr)^pf は (Φ^pf, ℚ·Φ^birat) の model Frobenioid",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
