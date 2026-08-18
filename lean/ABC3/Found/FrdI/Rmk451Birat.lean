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


/-! ## ★4. ★★★★比較関手 —— 恒等と合成の保存

★★合成の保存が要点である。`compBirat` は `Proposition 1.11, (vii)` の**引き戻しの選択**を
使うので、両側の選択は**独立**であり、定義的には一致しない。
★★★しかし添字圏は**細い**(`idxBirat_hom_ext`)ので、
**共通の上界を取れば図式は自動的に可換**になり、選択の違いが消える。
★これが `compBirat_natural_right` と同じ骨である。 -/

/-- ★★`biratMap` は `𝒞 → 𝒞^birat` の像と可換(`toHomBirat` を保つ)。 -/
theorem biratMap_toHomBirat (G : Frobenioid P) {A B : UnTr P} (φ : A ⟶ B) :
    biratMap (F := F) G A B (toHomBirat φ)
      = toHomBirat ((unTr2 (P := P) (F := F)).map φ) := by
  show biratMap (F := F) G A B
      (HomBirat.mk (idxBiratOne (unTrPre P F) (unTr_frobenioid P F G) A) φ) = _
  rw [biratMap_mk]
  rfl

/-- ★`biratMap` は恒等を保つ。 -/
theorem biratMap_id (G : Frobenioid P) (A : UnTr P) :
    biratMap (F := F) G A A (𝟙 (show BiratCat (unTrPre P F) (unTr_frobenioid P F G) from A))
      = 𝟙 (show BiratCat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) from
          (unTr2 (P := P) (F := F)).obj A) := by
  show biratMap (F := F) G A A (toHomBirat (𝟙 A)) = _
  rw [biratMap_toHomBirat]
  rfl

set_option maxHeartbeats 1000000 in
/-- ★★★**引き戻しの `α` は共通の上界の上で移送される**。

★`W₂` の構造射が **mono**(pre-step)なので消約でき、
`biratPull_sq` を両側で使うと `γ` の等式(`hkey`、細さから)に帰着する。 -/
theorem biratPull_alpha_transport (G : Frobenioid P) {A B : UnTr P}
    (Z : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B)
    (W : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) B)
    {V : IdxBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
      ((unTr2 (P := P) (F := F)).obj A)}
    (c : (idx2 (F := F) G A).obj (biratPullIdx (unTr_frobenioid P F G).core Z φ W) ⟶ V)
    (c' : biratPullIdx (unTr2G G).core ((idx2 (F := F) G A).obj Z)
        ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W) ⟶ V)
    (hkey : c.unop.left.hom ≫ (unTr2 (P := P) (F := F)).map
          (biratPullGamma (unTr_frobenioid P F G).core Z φ W)
        = c'.unop.left.hom ≫ biratPullGamma (unTr2G G).core ((idx2 (F := F) G A).obj Z)
          ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W)) :
    c.unop.left.hom ≫ (unTr2 (P := P) (F := F)).map
        (biratPullAlpha (unTr_frobenioid P F G).core Z φ W)
      = c'.unop.left.hom ≫ biratPullAlpha (unTr2G G).core ((idx2 (F := F) G A).obj Z)
          ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W) := by
  haveI hb : Mono ((idx2 (F := F) G B).obj W).unop.hom.hom :=
    (unTr2G G).core.preStepMono _ ((idx2 (F := F) G B).obj W).unop.hom.property.2
  refine (cancel_mono ((idx2 (F := F) G B).obj W).unop.hom.hom).mp ?_
  have h1 := biratPull_sq (unTr_frobenioid P F G).core Z φ W
  have h2 := biratPull_sq (unTr2G G).core ((idx2 (F := F) G A).obj Z)
    ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W)
  have e1 : (unTr2 (P := P) (F := F)).map
        (biratPullAlpha (unTr_frobenioid P F G).core Z φ W)
      ≫ ((idx2 (F := F) G B).obj W).unop.hom.hom
      = (unTr2 (P := P) (F := F)).map
          (biratPullGamma (unTr_frobenioid P F G).core Z φ W)
        ≫ (unTr2 (P := P) (F := F)).map φ := by
    have h := congrArg (unTr2 (P := P) (F := F)).map h1
    rw [CategoryTheory.Functor.map_comp] at h
    exact (h.trans (CategoryTheory.Functor.map_comp _ _ _)).symm
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c.unop.left.hom ≫ t) e1) ?_
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  refine Eq.trans (congrArg (fun t => t ≫ (unTr2 (P := P) (F := F)).map φ) hkey) ?_
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c'.unop.left.hom ≫ t) h2) ?_
  exact (Category.assoc _ _ _).symm

set_option maxHeartbeats 1000000 in
/-- ★★★代表元での合成の保存。★細さ(`idxBirat_hom_ext`)がすべてを潰す。 -/
theorem biratMap_comp_mk (G : Frobenioid P) {A B E : UnTr P}
    (Z : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B)
    (W : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) B) (ψ : W.unop.left.obj ⟶ E) :
    HomBirat.mk ((idx2 (F := F) G A).obj
        (biratPullIdx (unTr_frobenioid P F G).core Z φ W))
      ((unTr2 (P := P) (F := F)).map
        (biratPullAlpha (unTr_frobenioid P F G).core Z φ W ≫ ψ))
    = HomBirat.mk
        (biratPullIdx (unTr2G G).core ((idx2 (F := F) G A).obj Z)
          ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W))
        (biratPullAlpha (unTr2G G).core ((idx2 (F := F) G A).obj Z)
          ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W)
          ≫ (unTr2 (P := P) (F := F)).map ψ) := by
  refine HomBirat.sound _
    (IsFiltered.leftToMax ((idx2 (F := F) G A).obj
      (biratPullIdx (unTr_frobenioid P F G).core Z φ W))
      (biratPullIdx (unTr2G G).core ((idx2 (F := F) G A).obj Z)
        ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W)))
    (IsFiltered.rightToMax ((idx2 (F := F) G A).obj
      (biratPullIdx (unTr_frobenioid P F G).core Z φ W))
      (biratPullIdx (unTr2G G).core ((idx2 (F := F) G A).obj Z)
        ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W))) ?_
  have hkey := congrArg (fun t => t.unop.left.hom)
    (idxBirat_hom_ext
      ((idx2 (F := F) G A).map (biratPullHom (unTr_frobenioid P F G).core Z φ W)
        ≫ IsFiltered.leftToMax _ _)
      (biratPullHom (unTr2G G).core ((idx2 (F := F) G A).obj Z)
        ((unTr2 (P := P) (F := F)).map φ) ((idx2 (F := F) G B).obj W)
        ≫ IsFiltered.rightToMax _ _))
  have hmid := biratPull_alpha_transport (F := F) G Z φ W _ _ hkey
  simp only [CategoryTheory.Functor.map_comp, ← Category.assoc]
  exact Eq.trans (Category.assoc _ _ _).symm
    (congrArg (fun t => t ≫ (unTr2 (P := P) (F := F)).map ψ) hmid)

set_option maxHeartbeats 1000000 in
/-- ★★★`biratMap` は合成を保つ。 -/
theorem biratMap_comp (G : Frobenioid P) {A B E : UnTr P}
    (f : HomBirat (unTrPre P F) (unTr_frobenioid P F G) A B)
    (g : HomBirat (unTrPre P F) (unTr_frobenioid P F G) B E) :
    biratMap (F := F) G A E
        (compBirat (unTrPre P F) (unTr_frobenioid P F G) (unTr_frobenioid P F G).core f g)
      = compBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)
          (unTr2G G).core (biratMap (F := F) G A B f) (biratMap (F := F) G B E g) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratMap_mk, biratMap_mk, biratMap_mk, compBirat_mk]
  exact biratMap_comp_mk (F := F) G Z φ W ψ

/-- ★★★★**段4** —— **birationalization の比較関手**。

★★これで `(𝒞^un-tr)^birat` から `((𝒞^istr)^un-tr)^birat` へ渡れる。 -/
noncomputable def biratFunctor (G : Frobenioid P) :
    BiratCat (unTrPre P F) (unTr_frobenioid P F G) ⥤
      BiratCat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) where
  obj A := (unTr2 (P := P) (F := F)).obj A
  map {A B} f := biratMap (F := F) G A B f
  map_id A := biratMap_id (F := F) G A
  map_comp f g := biratMap_comp (F := F) G f g


/-! ## ★5. 逆向きの写像

★★`unTr2` は**圏の同型**である —— 対象も射も `rfl` で往復する
(`unTr2Inv.obj (unTr2.obj A) = A` および `unTr2Inv.map (unTr2.map f) = f`)。
★そこで逆向きの比較も同じ手順で組める。 -/

/-- ★逆向きの `coaPreProp` の対応。 -/
theorem unTr2Inv_coaPreProp (G : Frobenioid P) {X Y : UnTr (istrPre P F)}
    (f : X ⟶ Y) (hf : coaPreProp (unTrPre (istrPre P F) (istr_frobenioidCore P F)) f) :
    coaPreProp (unTrPre P F) ((unTr2Inv (P := P) (F := F)).map f) := by
  refine (unTr2_coaPreProp_iff (F := F) G ((unTr2Inv (P := P) (F := F)).map f)).mp ?_
  have h : (unTr2 (P := P) (F := F)).map ((unTr2Inv (P := P) (F := F)).map f) = f := by
    refine Quotient.inductionOn f (fun α => ?_); rfl
  rw [h]
  exact hf

/-- ★★逆向きの `𝒞^{coa-pre}` の対応。 -/
noncomputable def coaPre2Inv (G : Frobenioid P) :
    CoaPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) ⥤
      CoaPre (unTrPre P F) (unTr_frobenioid P F G) where
  obj X := ⟨(unTr2Inv (P := P) (F := F)).obj X.obj⟩
  map {X Y} f := ⟨(unTr2Inv (P := P) (F := F)).map f.hom,
    unTr2Inv_coaPreProp (F := F) G f.hom f.property⟩
  map_id X := by
    refine WideSubcategory.hom_ext _ ?_
    show (unTr2Inv (P := P) (F := F)).map (𝟙 X.obj) = 𝟙 _
    exact CategoryTheory.Functor.map_id _ _
  map_comp {X Y Z} f g := by
    refine WideSubcategory.hom_ext _ ?_
    show (unTr2Inv (P := P) (F := F)).map (f.hom ≫ g.hom)
      = (unTr2Inv (P := P) (F := F)).map f.hom ≫ (unTr2Inv (P := P) (F := F)).map g.hom
    exact CategoryTheory.Functor.map_comp _ _ _

/-- ★★逆向きの添字圏の対応。 -/
noncomputable abbrev idx2Inv (G : Frobenioid P) (X : UnTr (istrPre P F)) :
    IdxBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) X ⥤
      IdxBirat (unTrPre P F) (unTr_frobenioid P F G)
        ((unTr2Inv (P := P) (F := F)).obj X) :=
  (Over.post (coaPre2Inv (F := F) G)).op

/-- ★逆向きの余錐。 -/
noncomputable def biratMapInvCocone (G : Frobenioid P) (X Y : UnTr (istrPre P F)) :
    Cocone (homFunctorBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F))
      (unTr2G G) X Y) :=
  Cocone.mk (HomBirat (unTrPre P F) (unTr_frobenioid P F G)
      ((unTr2Inv (P := P) (F := F)).obj X) ((unTr2Inv (P := P) (F := F)).obj Y))
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomBirat.mk ((idx2Inv (F := F) G X).obj Z) ((unTr2Inv (P := P) (F := F)).map φ.down)
      naturality := fun Z W u => by
        ext φ
        show HomBirat.mk ((idx2Inv (F := F) G X).obj W)
              ((unTr2Inv (P := P) (F := F)).map (u.unop.left.hom ≫ φ.down))
          = HomBirat.mk ((idx2Inv (F := F) G X).obj Z)
              ((unTr2Inv (P := P) (F := F)).map φ.down)
        rw [CategoryTheory.Functor.map_comp]
        exact HomBirat.mk_map ((idx2Inv (F := F) G X).map u)
          ((unTr2Inv (P := P) (F := F)).map φ.down) }

/-- ★★逆向きの `Hom^birat` の写像。 -/
noncomputable def biratMapInv (G : Frobenioid P) (X Y : UnTr (istrPre P F)) :
    HomBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) X Y →
      HomBirat (unTrPre P F) (unTr_frobenioid P F G)
        ((unTr2Inv (P := P) (F := F)).obj X) ((unTr2Inv (P := P) (F := F)).obj Y) :=
  fun z => colimit.desc _ (biratMapInvCocone (F := F) G X Y) z

/-- ★代表元での値。 -/
theorem biratMapInv_mk (G : Frobenioid P) (X Y : UnTr (istrPre P F))
    (Z : IdxBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) X)
    (φ : Z.unop.left.obj ⟶ Y) :
    biratMapInv (F := F) G X Y (HomBirat.mk Z φ)
      = HomBirat.mk ((idx2Inv (F := F) G X).obj Z)
          ((unTr2Inv (P := P) (F := F)).map φ) := by
  show colimit.desc _ (biratMapInvCocone (F := F) G X Y)
      (colimit.ι (homFunctorBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F))
        (unTr2G G) X Y) Z (ULift.up φ)) = _
  rw [← types_comp_apply
    (colimit.ι (homFunctorBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F))
      (unTr2G G) X Y) Z)
    (colimit.desc _ (biratMapInvCocone (F := F) G X Y)), colimit.ι_desc]
  rfl

/-! ## ★6. 残り —— 全単射と `Base`/`degFr` の保存

★★**数学は済んでいる。残りは構造の畳み込み(eta)である。**

★往復を代表元で計算すると
`biratMapInv (biratMap (mk Z φ)) = mk ((idx2Inv).obj ((idx2).obj Z)) (unTr2Inv.map (unTr2.map φ))`
になる。射の側は `rfl`(`unTr2` が圏の同型だから)。★**添字対象の側が `rfl` にならない**:

| 部分 | 状態 |
|---|---|
| `((idx2Inv).obj ((idx2).obj Z)).unop.left = Z.unop.left` | ★`rfl` |
| `((idx2Inv).obj ((idx2).obj Z)).unop.hom.hom = Z.unop.hom.hom` | ★**`rfl` にならない** |

★`Over.post` の合成が `Over.mk` を作り直すため、`Over`(= `Comma`)の構造 eta を
Lean が畳まない。★★**次の一手**: 恒等射を左成分に持つ射
`Z ⟶ (idx2Inv).obj ((idx2).obj Z)` を `Over.homMk (𝟙 _)` で作り、
`HomBirat.mk_map` と `Category.id_comp` で潰す(添字圏は細いので射の一意性は自動)。
★これが済めば `biratMap_bijective` が出て、`Base`/`degFr` の保存(代表元で `rfl` の見込み)と
併せて `isFrobeniusCompact_transport` が当たり、`Remark 4.5.1` が閉じる。
-/

end ABC3.Found.FrdI
