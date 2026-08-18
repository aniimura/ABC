/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Rmk451
import ABC3.Found.FrdI.Prop44

/-!
# [FrdI] `Remark 4.5.1` —— birationalization の比較

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.86。

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

## ★何をするファイルか

`Rmk451.lean` で **`𝒞^un-tr ≌ (𝒞^istr)^un-tr`**(`unTr2Equiv`)まで作った。
★★`Definition 4.5, (iii), (b)` は **`(𝒞^un-tr)^birat` の Frobenius-compact 対象**を
要求するので、その birationalization どうしを繋ぐ必要がある。

## ★段取りと現在地

| 段 | 内容 | 宣言 | 状態 |
|---|---|---|---|
| 1 | `𝒞^{coa-pre}` の対応 | `coaPre2` | ★済 |
| 2 | スライス(＝添字圏)の対応 | `slice2` / `idx2` | ★済 |
| 3 | `Hom^birat` の写像 | `biratMapCocone` / `biratMap` / `biratMap_mk` | ★済 |
| 4 | 関手性・全単射・`Base`/`degFr` の保存 | —— | 未 |

★★**段1 が効くのは `coaPreProp` が `unTr2` で両向きに移るから**である
(`unTr2_coaPreProp_iff`)。★その根拠は
**`𝒞^un-tr` も `(𝒞^istr)^un-tr` も isotropic 型**(`unTr_isotropic`)であり、
co-angular がすべての射について自動になること。

★★段3 は既存の `HomColim` に「添字圏をまたぐ写像」の API が無いので、
`bimapCocone` と同じ型で**余錐を組んで `colimit.desc` で降ろす**形にした。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w v2 u2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★1. `(𝒞^istr)^un-tr` 側の Frobenioid と `𝒞^{coa-pre}` の対応 -/

/-- ★`(𝒞^istr)^un-tr` の Frobenioid 構造。 -/
noncomputable abbrev unTr2G (G : Frobenioid P) :
    Frobenioid (unTrPre (istrPre P F) (istr_frobenioidCore P F)) :=
  unTr_frobenioid (istrPre P F) (istr_frobenioidCore P F) (istr_frobenioid P F G)

/-- ★★**段1** —— `𝒞^{coa-pre}` の対応。

★`coaPreProp` が `unTr2` で両向きに移る(`unTr2_coaPreProp_iff`)ので、
ワイド部分圏がそのまま対応する。 -/
noncomputable def coaPre2 (G : Frobenioid P) :
    CoaPre (unTrPre P F) (unTr_frobenioid P F G) ⥤
      CoaPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) where
  obj X := ⟨(unTr2 (P := P) (F := F)).obj X.obj⟩
  map {X Y} f := ⟨(unTr2 (P := P) (F := F)).map f.hom,
    (unTr2_coaPreProp_iff (F := F) G f.hom).mpr f.property⟩
  map_id X := by
    refine WideSubcategory.hom_ext _ ?_
    show (unTr2 (P := P) (F := F)).map (𝟙 X.obj) = 𝟙 _
    exact CategoryTheory.Functor.map_id _ _
  map_comp {X Y Z} f g := by
    refine WideSubcategory.hom_ext _ ?_
    show (unTr2 (P := P) (F := F)).map (f.hom ≫ g.hom)
      = (unTr2 (P := P) (F := F)).map f.hom ≫ (unTr2 (P := P) (F := F)).map g.hom
    exact CategoryTheory.Functor.map_comp _ _ _

/-! ## ★2. スライスと添字圏の対応 -/

/-- ★★**段2** —— スライス `(𝒞^{coa-pre})_A` の対応。 -/
noncomputable abbrev slice2 (G : Frobenioid P) (A : UnTr P) :
    SliceA (unTrPre P F) (unTr_frobenioid P F G) A ⥤
      SliceA (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
        ((unTr2 (P := P) (F := F)).obj A) :=
  Over.post (coaPre2 (F := F) G)

/-- ★★**段2'** —— `Hom^birat` の添字圏(スライスの反対圏)の対応。 -/
noncomputable abbrev idx2 (G : Frobenioid P) (A : UnTr P) :
    IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A ⥤
      IdxBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
        ((unTr2 (P := P) (F := F)).obj A) :=
  (slice2 (F := F) G A).op

/-- ★添字対象の底は `unTr2` で移ったもの(`rfl`)。 -/
theorem idx2_obj_left (G : Frobenioid P) (A : UnTr P)
    (Z : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A) :
    ((idx2 (F := F) G A).obj Z).unop.left.obj
      = (unTr2 (P := P) (F := F)).obj Z.unop.left.obj := rfl

/-! ## ★3. `Hom^birat` の写像

★`HomColim` には「添字圏をまたぐ写像」の API が無いので、
`bimapCocone` と同じ型で**余錐を組んで `colimit.desc` で降ろす**。 -/

/-- ★★**段3** —— 降ろすための余錐。自然性は `HomBirat.mk_map` 1 本。 -/
noncomputable def biratMapCocone (G : Frobenioid P) (A B : UnTr P) :
    Cocone (homFunctorBirat (unTrPre P F) (unTr_frobenioid P F G) A B) :=
  Cocone.mk (HomBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
      ((unTr2 (P := P) (F := F)).obj A) ((unTr2 (P := P) (F := F)).obj B))
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomBirat.mk ((idx2 (F := F) G A).obj Z) ((unTr2 (P := P) (F := F)).map φ.down)
      naturality := fun Z W u => by
        ext φ
        show HomBirat.mk ((idx2 (F := F) G A).obj W)
              ((unTr2 (P := P) (F := F)).map (u.unop.left.hom ≫ φ.down))
          = HomBirat.mk ((idx2 (F := F) G A).obj Z)
              ((unTr2 (P := P) (F := F)).map φ.down)
        rw [CategoryTheory.Functor.map_comp]
        exact HomBirat.mk_map ((idx2 (F := F) G A).map u)
          ((unTr2 (P := P) (F := F)).map φ.down) }

/-- ★★★**段3** —— `Hom^birat` の写像。 -/
noncomputable def biratMap (G : Frobenioid P) (A B : UnTr P) :
    HomBirat (unTrPre P F) (unTr_frobenioid P F G) A B →
      HomBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
        ((unTr2 (P := P) (F := F)).obj A) ((unTr2 (P := P) (F := F)).obj B) :=
  fun z => colimit.desc _ (biratMapCocone (F := F) G A B) z

/-- ★代表元での値 —— これが以後の計算の入口になる。 -/
theorem biratMap_mk (G : Frobenioid P) (A B : UnTr P)
    (Z : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B) :
    biratMap (F := F) G A B (HomBirat.mk Z φ)
      = HomBirat.mk ((idx2 (F := F) G A).obj Z) ((unTr2 (P := P) (F := F)).map φ) := by
  show colimit.desc _ (biratMapCocone (F := F) G A B)
      (colimit.ι (homFunctorBirat (unTrPre P F) (unTr_frobenioid P F G) A B) Z (ULift.up φ)) = _
  rw [← types_comp_apply
    (colimit.ι (homFunctorBirat (unTrPre P F) (unTr_frobenioid P F G) A B) Z)
    (colimit.desc _ (biratMapCocone (F := F) G A B)), colimit.ι_desc]
  rfl

/-! ## ★4. 残り —— 関手性・全単射・`Base`/`degFr` の保存

★`biratMap` が関手をなすこと(恒等と合成)、全単射であること
(`unTr2Inv` から逆向きを作り、往復が恒等であることを `HomBirat.eq_iff` で見る)、
`Base`・`degFr` を保つことを示せば、
`isFrobeniusCompact_transport`(`Rmk451.lean`)がそのまま当たり、
`Definition 4.5, (iii), (b)` が `𝒞^istr` へ移る。
-/

end ABC3.Found.FrdI
