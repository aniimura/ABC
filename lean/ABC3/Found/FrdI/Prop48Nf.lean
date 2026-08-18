/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop21
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 4.8, (ii) —— naive Frobenius 関手を `𝒞^birat` へ降ろす準備

原文 (FrdI p.85):
> naive Frobenius functor

★★`𝒞^birat` の射は `HomBirat A B = HomColim (homFunctorBirat P G A B)` で、
添字は `IdxBirat P G A`(`A` への **co-angular pre-step** の反対圏)、
遷移は**前合成**である。

★★★したがって `naiveFrob` を降ろすのに要るのは
**「`naiveFrob` が co-angular pre-step を保つ」**の 1 本だけ。

## ★★★在庫の見落とし 10 件目(2026-08-19)

最初、`nfMap` が linear / base-isomorphism / pre-step を保つことを
**手で書き下してしまった**。★実際には `Prop21.lean` に
`nfMap_preStep` / `nfMap_frobType` / `nfMap_pullBack` / `prop_2_1_ii_degFr` があり、
さらに `Prop110.lean` に **`prop_1_10_i_coAngular_of`** まであった。

★★**検索語を `naiveFrob` / `prop_2_1` に限ったのが原因**である ——
`nfMap_*` という実際の名前で引いていれば 1 回で当たった。
★対策: **「性質 P が射のクラス Q で保たれる」を書く前に、
`prop_<原典番号>_.*_of` と `<構成の接頭辞>_<性質>` の両方で引く**こと。

★下は在庫を組むだけの 2 本になった。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (F : FrobenioidCore P) (d : ℕ+)

/-- ★★★★**`naiveFrob` は co-angular 射を保つ** —— `prop_1_10_i_coAngular_of` に
`nfMap` の四角形を当てるだけ。 -/
theorem nfMap_coAngular {A B : C} (φ : A ⟶ B) (h : IsCoAngular P φ) :
    IsCoAngular P (nfMap P F d φ) :=
  prop_1_10_i_coAngular_of P F (nfHom_frobType P F d A) (nfHom_frobType P F d B)
    (nfMap_sq P F d φ) h

/-- ★★★★★**`naiveFrob` は co-angular pre-step を保つ**。

★★これが `𝒞^birat` の添字圏 `IdxBirat`(co-angular pre-step のなす圏)を
`naiveFrob` で送るために要る唯一の性質である。 -/
theorem nfMap_coaPreStep {A B : C} (φ : A ⟶ B) (hca : IsCoAngular P φ) (hps : IsPreStep P φ) :
    IsCoAngular P (nfMap P F d φ) ∧ IsPreStep P (nfMap P F d φ) :=
  ⟨nfMap_coAngular P F d φ hca, nfMap_preStep P F d φ hps⟩

/-- ★**`𝒞^coa-pre` の射のクラスとして書いた形** —— `coaPreProp` を保つ。 -/
theorem nfMap_coaPreProp {A B : C} (φ : A ⟶ B) (h : coaPreProp P φ) :
    coaPreProp P (nfMap P F d φ) :=
  nfMap_coaPreStep P F d φ h.1 h.2

def nfMap_coaPreStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85,
    item := "Proposition 4.8, (ii) — naive Frobenius 関手は co-angular pre-step を保つ",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
