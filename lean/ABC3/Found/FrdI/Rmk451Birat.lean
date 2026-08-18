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


/-! ## ★6. ★★★★全単射

★★**罠の記録**: `unTr2Inv.map (unTr2.map f) = f` は **`rfl` にならない**。
`Quotient.liftOn` は**変数のままでは簡約しない**ので、代表元での場合分けが要る
(`unTr2_round`)。★対象の側(`unTr2Inv.obj (unTr2.obj A) = A`)は `rfl` である。

★そのため添字対象の往復も `rfl` にならず、
**恒等射を左成分に持つ射**(`idxRoundHom`)を作って `HomBirat.mk_map` で潰す。
★添字圏は細いので、この射の取り方に任意性は無い。 -/

/-- ★往復(合成関手を経由しない形。`rw` / `congrArg` で使うにはこちらが要る)。 -/
theorem unTr2_round {A B : UnTr P} (f : A ⟶ B) :
    (unTr2Inv (P := P) (F := F)).map ((unTr2 (P := P) (F := F)).map f) = f := by
  refine Quotient.inductionOn f (fun α => ?_); rfl

/-- ★逆向きの往復。 -/
theorem unTr2Inv_round {X Y : UnTr (istrPre P F)} (g : X ⟶ Y) :
    (unTr2 (P := P) (F := F)).map ((unTr2Inv (P := P) (F := F)).map g) = g := by
  refine Quotient.inductionOn g (fun α => ?_); rfl

/-- ★★往復した添字対象への射(左成分は恒等)。 -/
noncomputable def idxRoundHom (G : Frobenioid P) (A : UnTr P)
    (Z : IdxBirat (unTrPre P F) (unTr_frobenioid P F G) A) :
    Z ⟶ (idx2Inv (F := F) G ((unTr2 (P := P) (F := F)).obj A)).obj
        ((idx2 (F := F) G A).obj Z) :=
  Quiver.Hom.op (Over.homMk
    (𝟙 (show CoaPre (unTrPre P F) (unTr_frobenioid P F G) from Z.unop.left))
    ((Category.id_comp Z.unop.hom).trans
      (WideSubcategory.hom_ext (coaPrePropOf (unTrPre P F) (unTr_frobenioid P F G))
        (unTr2_round (F := F) Z.unop.hom.hom).symm)))

/-- ★★逆向きの往復射。 -/
noncomputable def idxRoundHom' (G : Frobenioid P) (X : UnTr (istrPre P F))
    (Z : IdxBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) X) :
    Z ⟶ (idx2 (F := F) G ((unTr2Inv (P := P) (F := F)).obj X)).obj
        ((idx2Inv (F := F) G X).obj Z) :=
  Quiver.Hom.op (Over.homMk
    (𝟙 (show CoaPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) from
      Z.unop.left))
    ((Category.id_comp Z.unop.hom).trans
      (WideSubcategory.hom_ext
        (coaPrePropOf (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G))
        (unTr2Inv_round (F := F) Z.unop.hom.hom).symm)))

/-- ★★★往復は恒等(`𝒞^un-tr` 側)。 -/
theorem biratMapInv_biratMap (G : Frobenioid P) (A B : UnTr P)
    (z : HomBirat (unTrPre P F) (unTr_frobenioid P F G) A B) :
    biratMapInv (F := F) G _ _ (biratMap (F := F) G A B z) = z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratMap_mk, biratMapInv_mk]
  have h2 := HomBirat.mk_map (idxRoundHom (F := F) G A Z) φ
  rw [show (idxRoundHom (F := F) G A Z).unop.left.hom = 𝟙 Z.unop.left.obj from rfl] at h2
  have h3 := (congrArg (HomBirat.mk
      ((idx2Inv (F := F) G ((unTr2 (P := P) (F := F)).obj A)).obj
        ((idx2 (F := F) G A).obj Z))) (Category.id_comp φ)).symm.trans h2
  exact (congrArg (HomBirat.mk
      ((idx2Inv (F := F) G ((unTr2 (P := P) (F := F)).obj A)).obj
        ((idx2 (F := F) G A).obj Z))) (unTr2_round (F := F) φ)).trans h3

/-- ★★★往復は恒等(`(𝒞^istr)^un-tr` 側)。 -/
theorem biratMap_biratMapInv (G : Frobenioid P) (X Y : UnTr (istrPre P F))
    (z : HomBirat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G) X Y) :
    biratMap (F := F) G _ _ (biratMapInv (F := F) G X Y z) = z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratMapInv_mk, biratMap_mk]
  have h2 := HomBirat.mk_map (idxRoundHom' (F := F) G X Z) φ
  rw [show (idxRoundHom' (F := F) G X Z).unop.left.hom = 𝟙 Z.unop.left.obj from rfl] at h2
  have h3 := (congrArg (HomBirat.mk
      ((idx2 (F := F) G ((unTr2Inv (P := P) (F := F)).obj X)).obj
        ((idx2Inv (F := F) G X).obj Z))) (Category.id_comp φ)).symm.trans h2
  exact (congrArg (HomBirat.mk
      ((idx2 (F := F) G ((unTr2Inv (P := P) (F := F)).obj X)).obj
        ((idx2Inv (F := F) G X).obj Z))) (unTr2Inv_round (F := F) φ)).trans h3

/-- ★★★★**`biratMap` は全単射** —— これで比較関手は充満忠実になる。 -/
theorem biratMap_bijective (G : Frobenioid P) (A B : UnTr P) :
    Function.Bijective (biratMap (F := F) G A B) := by
  constructor
  · intro z z' h
    have h1 := congrArg (biratMapInv (F := F) G _ _) h
    rwa [biratMapInv_biratMap, biratMapInv_biratMap] at h1
  · intro z
    exact ⟨biratMapInv (F := F) G _ _ z, biratMap_biratMapInv (F := F) G _ _ z⟩


/-! ## ★7. ★★★★`Base` と `degFr` の保存

★これで比較関手は **pre-Frobenioid の構造まで保つ圏同値**になる。 -/

/-- ★★`biratMap` は Frobenius 次数を保つ。 -/
theorem biratMap_degFr (G : Frobenioid P) (A B : UnTr P)
    (z : HomBirat (unTrPre P F) (unTr_frobenioid P F G) A B) :
    biratDeg (biratMap (F := F) G A B z) = biratDeg z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratMap_mk, biratDeg_mk, biratDeg_mk]
  exact unTr2_degFr (F := F) φ

/-- ★★`biratMap` は底を保つ。

★`sliceBaseOf = inv (Base a) ≫ Base φ` に開いてから `congr` で 2 成分に分け、
どちらも `unTr2_Base` で閉じる(`inv` の側はさらに 1 段 `congr` する)。 -/
theorem biratMap_base (G : Frobenioid P) (A B : UnTr P)
    (z : HomBirat (unTrPre P F) (unTr_frobenioid P F G) A B) :
    biratBase (biratMap (F := F) G A B z) = biratBase z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratMap_mk, biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq]
  congr 1
  · congr 1
    exact unTr2_Base (F := F) Z.unop.hom.hom
  · exact unTr2_Base (F := F) φ

/-- ★★`biratFunctor` は `degFr` を保つ(pre-Frobenioid の言葉で)。 -/
theorem biratFunctor_degFr (G : Frobenioid P) {A B : BiratCat (unTrPre P F)
    (unTr_frobenioid P F G)} (f : A ⟶ B) :
    (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)).degFr
        ((biratFunctor (F := F) G).map f)
      = (biratPre (unTrPre P F) (unTr_frobenioid P F G)).degFr f :=
  biratMap_degFr (F := F) G A B f

/-- ★★`biratFunctor` は `Base` を保つ。 -/
theorem biratFunctor_base (G : Frobenioid P) {A B : BiratCat (unTrPre P F)
    (unTr_frobenioid P F G)} (f : A ⟶ B) :
    (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)).Base
        ((biratFunctor (F := F) G).map f)
      = (biratPre (unTrPre P F) (unTr_frobenioid P F G)).Base f :=
  biratMap_base (F := F) G A B f


/-! ## ★8. ★★★★★`Frobenius-compact` 対象の移送 —— `Definition 4.5, (iii), (b)`

★`isFrobeniusCompact_transport`(`Rmk451.lean`)の 4 材料を、
充満忠実性と構造保存から組む。 -/

instance biratFunctor_full (G : Frobenioid P) : (biratFunctor (F := F) G).Full where
  map_surjective {A B} f := (biratMap_bijective (F := F) G A B).2 f

instance biratFunctor_faithful (G : Frobenioid P) : (biratFunctor (F := F) G).Faithful where
  map_injective {A B} {_ _} h := (biratMap_bijective (F := F) G A B).1 h

/-- ★★材料 1 —— `End` の単系同型。 -/
noncomputable def biratEndEquiv (G : Frobenioid P)
    (X : BiratCat (unTrPre P F) (unTr_frobenioid P F G)) :
    End X ≃* End ((biratFunctor (F := F) G).obj X) :=
  MulEquiv.ofBijective (CategoryTheory.Functor.mapEnd X (biratFunctor (F := F) G))
    (biratMap_bijective (F := F) G X X)

/-- ★★材料 2 —— 自己同型の対応。 -/
noncomputable def biratIsoEquiv (G : Frobenioid P)
    (X : BiratCat (unTrPre P F) (unTr_frobenioid P F G)) :
    (X ≅ X) ≃ ((biratFunctor (F := F) G).obj X ≅ (biratFunctor (F := F) G).obj X) where
  toFun θ := (biratFunctor (F := F) G).mapIso θ
  invFun θ' := (Functor.FullyFaithful.ofFullyFaithful (biratFunctor (F := F) G)).preimageIso θ'
  left_inv θ := by
    refine Iso.ext ?_
    exact (biratFunctor (F := F) G).map_injective
      (by simp [Functor.FullyFaithful.preimageIso])
  right_inv θ' := by
    refine Iso.ext ?_
    simp [Functor.FullyFaithful.preimageIso]

theorem bf_base (G : Frobenioid P) {A B : BiratCat (unTrPre P F) (unTr_frobenioid P F G)}
    (f : A ⟶ B) :
    (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)).Base
        ((biratFunctor (F := F) G).map f)
      = (biratPre (unTrPre P F) (unTr_frobenioid P F G)).Base f :=
  biratMap_base (F := F) G A B f

theorem bf_degFr (G : Frobenioid P) {A B : BiratCat (unTrPre P F) (unTr_frobenioid P F G)}
    (f : A ⟶ B) :
    (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)).degFr
        ((biratFunctor (F := F) G).map f)
      = (biratPre (unTrPre P F) (unTr_frobenioid P F G)).degFr f :=
  biratMap_degFr (F := F) G A B f

/-- ★★材料 3 —— `OTimes` の対応。★`IsBaseIdentity`・`IsLinear` は構造保存から。 -/
theorem birat_otimes_iff (G : Frobenioid P)
    (X : BiratCat (unTrPre P F) (unTr_frobenioid P F G)) (u : End X) :
    u ∈ OTimes (biratPre (unTrPre P F) (unTr_frobenioid P F G)) X ↔
      biratEndEquiv (F := F) G X u
        ∈ OTimes (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G))
          ((biratFunctor (F := F) G).obj X) := by
  have hid : (biratFunctor (F := F) G).map (𝟙 X) = 𝟙 ((biratFunctor (F := F) G).obj X) :=
    CategoryTheory.Functor.map_id _ _
  constructor
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, hu.map (biratEndEquiv (F := F) G X).toMonoidHom⟩
    · show (biratPre _ (unTr2G G)).Base ((biratFunctor (F := F) G).map (u : X ⟶ X))
        = (biratPre _ (unTr2G G)).Base (𝟙 ((biratFunctor (F := F) G).obj X))
      rw [bf_base, ← hid, bf_base]
      exact hb
    · show (biratPre _ (unTr2G G)).degFr ((biratFunctor (F := F) G).map (u : X ⟶ X)) = 1
      rw [bf_degFr]
      exact hl
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, ?_⟩
    · show (biratPre (unTrPre P F) (unTr_frobenioid P F G)).Base (u : X ⟶ X)
        = (biratPre (unTrPre P F) (unTr_frobenioid P F G)).Base (𝟙 X)
      rw [← bf_base (F := F) G (u : X ⟶ X), ← bf_base (F := F) G (𝟙 X), hid]
      exact hb
    · show (biratPre (unTrPre P F) (unTr_frobenioid P F G)).degFr (u : X ⟶ X) = 1
      rw [← bf_degFr (F := F) G (u : X ⟶ X)]
      exact hl
    · have h2 := hu.map (biratEndEquiv (F := F) G X).symm.toMonoidHom
      simpa using h2

/-- ★★材料 4 —— 共役との両立(関手性から)。 -/
theorem birat_endConj (G : Frobenioid P)
    (X : BiratCat (unTrPre P F) (unTr_frobenioid P F G)) (θ : X ≅ X) (u : End X) :
    biratEndEquiv (F := F) G X (endConj θ u)
      = endConj (biratIsoEquiv (F := F) G X θ) (biratEndEquiv (F := F) G X u) := by
  show (biratFunctor (F := F) G).map (θ.inv ≫ (u : X ⟶ X) ≫ θ.hom)
    = ((biratFunctor (F := F) G).mapIso θ).inv
      ≫ (biratFunctor (F := F) G).map (u : X ⟶ X)
      ≫ ((biratFunctor (F := F) G).mapIso θ).hom
  rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp]
  rfl

/-- ★★★★★**[FrdI] Definition 4.5, (iii), (b) が `𝒞^istr` へ移る** ——
`(𝒞^un-tr)^birat` の Frobenius-compact 対象から
`((𝒞^istr)^un-tr)^birat` の Frobenius-compact 対象が得られる。

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

★★これが `Remark 4.5.1` の 4 条のうち (b) である。 -/
theorem istr_unTrBiratCompact (G : Frobenioid P)
    (h : ∃ X : BiratCat (unTrPre P F) (unTr_frobenioid P F G),
      IsFrobeniusCompact (biratPre (unTrPre P F) (unTr_frobenioid P F G)) X) :
    ∃ X : BiratCat (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G),
      IsFrobeniusCompact
        (biratPre (unTrPre (istrPre P F) (istr_frobenioidCore P F)) (unTr2G G)) X := by
  obtain ⟨X, hX⟩ := h
  refine ⟨(biratFunctor (F := F) G).obj X, ?_⟩
  exact isFrobeniusCompact_transport _ _ (biratEndEquiv (F := F) G X)
    (biratIsoEquiv (F := F) G X) (birat_otimes_iff (F := F) G X)
    (birat_endConj (F := F) G X) hX

/-! ## ★9. 残り —— `Remark 4.5.1` の (a) の 2 条

★★**測定の訂正(2026-08-18)**: `Remark 4.5.1` の 4 条のうち
**(b) と standard 型は済んだ**が、(a) の残り 2 条
(**birationally Frobenius-normalized 型**・**rational 型**)は
**別の比較**を要する。

★これらは `𝒞^birat` と `(𝒞^istr)^birat` の比較であって、本ファイルで作った
`𝒞^un-tr` と `(𝒞^istr)^un-tr` の比較ではない。
★★★**そして `𝒞^istr` は `𝒞` の真の充満部分圏なので、その比較は圏同値ではない。**
`A` が isotropic でも `A′ → A` が pre-step のとき `A′` が isotropic とは限らないため、
添字圏 `SliceA` が一致しない。
★したがって (a) の 2 条は本ファイルの手法をそのまま流用できず、別立ての作業になる。
-/

end ABC3.Found.FrdI
