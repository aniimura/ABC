/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop21
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] Proposition 4.8, (ii) —— naive Frobenius 関手を `𝒞^birat` へ降ろす準備

原文 (FrdI p.88):
> (ii), observe that the naive Frobenius functor [cf. Proposition 2.1] determines a

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
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — naive Frobenius 関手は co-angular pre-step を保つ",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★段 2 —— 添字圏の写像

★`IdxBirat P G A = (Over (coaPreObj P G A))ᵒᵖ` なので、
`coaPreProp` を保つ関手から `Over.post` ＋ `.op` で降ろすだけ。
★★これは `iii-psipf` の `idxPfMap`(`Under.post`)と**同じ形**で、
向きが反対圏になるだけである。 -/

section Birat

variable (G : Frobenioid P)

/-- ★★★★**`naiveFrob` は `𝒞^{coa-pre}` の間の関手を誘導する**。 -/
noncomputable def coaPreNfMap : CoaPre P G ⥤ CoaPre P G where
  obj X := ⟨nfObj P G.core d X.obj⟩
  map f := ⟨nfMap P G.core d f.hom, nfMap_coaPreProp P G.core d f.hom f.property⟩
  map_id X := WideSubcategory.hom_ext _ ((naiveFrob P G.core d).map_id X.obj)
  map_comp f g := WideSubcategory.hom_ext _ ((naiveFrob P G.core d).map_comp f.hom g.hom)

/-- ★★**添字圏の写像** —— スライスへ降ろして反対圏を取る。 -/
noncomputable def idxBiratNfMap (A : C) :
    IdxBirat P G A ⥤ IdxBirat P G (nfObj P G.core d A) :=
  (Over.post (X := coaPreObj P G A) (coaPreNfMap P d G)).op

/-! ## ★★★段 3 —— 余錐と射の写像

★遷移は**前合成**なので、naturality は `nfMap` の関手性(`nfMap_comp`)そのもの。
★★`Ψ^pf` のときの `idxTransport_map` に当たる補題が**要らない**。 -/

/-- ★★★★**余錐** —— 各添字で `nfMap` を当てるだけ。 -/
noncomputable def biratNfCocone (A B : C) :
    Limits.Cocone (homFunctorBirat P G A B) where
  pt := HomBirat P G (nfObj P G.core d A) (nfObj P G.core d B)
  ι :=
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomBirat.mk ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ.down)
      naturality := by
        intro Z W u
        ext φ
        show HomBirat.mk ((idxBiratNfMap P d G A).obj W)
            (nfMap P G.core d (u.unop.left.hom ≫ φ.down))
          = HomBirat.mk ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ.down)
        rw [nfMap_comp]
        exact HomBirat.mk_map ((idxBiratNfMap P d G A).map u) (nfMap P G.core d φ.down) }

/-- ★★★★★**射の写像** —— 余極限の普遍性で降ろす。 -/
noncomputable def biratNfMap (A B : C) :
    HomBirat P G A B → HomBirat P G (nfObj P G.core d A) (nfObj P G.core d B) :=
  (Limits.colimit.desc (homFunctorBirat P G A B) (biratNfCocone P d G A B)).hom

/-- ★**代表元での計算則**。 -/
@[simp] theorem biratNfMap_mk {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    biratNfMap P d G A B (HomBirat.mk Z φ)
      = HomBirat.mk ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) :=
  congrFun (congrArg (fun t => t.hom)
    (Limits.colimit.ι_desc (biratNfCocone P d G A B) Z)) (ULift.up φ)

end Birat

def biratNfMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — naive Frobenius を 𝒞^birat の射へ降ろす",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
