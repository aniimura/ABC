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

/-! ## ★★段 4 —— 関手則の `map_id` -/

/-- ★★**添字圏の頂点は頂点に写る** —— `coaPreNfMap` の `map_id` だけ。 -/
theorem idxBiratNfMap_one (A : C) :
    (idxBiratNfMap P d G A).obj (idxBiratOne P G A)
      = idxBiratOne P G (nfObj P G.core d A) :=
  congrArg Opposite.op (congrArg Over.mk ((coaPreNfMap P d G).map_id (coaPreObj P G A)))

/-- ★★★★**`𝒞 → 𝒞^birat` と可換** —— これが `map_id` を与える。

★対象の等式で `rw` すると motive が壊れるので、
`iii-psipf` と同じく**添字圏の射を 1 本作って `HomBirat.mk_map` で移す**。 -/
theorem biratNfMap_toHomBirat {A B : C} (φ : A ⟶ B) :
    biratNfMap P d G A B (toHomBirat (P := P) (G := G) φ)
      = toHomBirat (P := P) (G := G) (nfMap P G.core d φ) := by
  let u : idxBiratOne P G (nfObj P G.core d A)
      ⟶ (idxBiratNfMap P d G A).obj (idxBiratOne P G A) :=
    Quiver.Hom.op (Over.homMk (𝟙 _) (by
      refine (Category.comp_id _).trans ?_
      exact ((coaPreNfMap P d G).map_id (coaPreObj P G A)).symm))
  have h1 : biratNfMap P d G A B (HomBirat.mk (idxBiratOne P G A) φ)
      = HomBirat.mk ((idxBiratNfMap P d G A).obj (idxBiratOne P G A))
          (nfMap P G.core d φ) :=
    biratNfMap_mk P d G (idxBiratOne P G A) φ
  refine h1.trans (Eq.trans ?_ (HomBirat.mk_map u (nfMap P G.core d φ)))
  exact congrArg (HomBirat.mk ((idxBiratNfMap P d G A).obj (idxBiratOne P G A)))
    (Category.id_comp (nfMap P G.core d φ)).symm

def biratNfMap_toHomBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — Ψ^birat は 𝒞 → 𝒞^birat と可換",
    sectionId := "frdi-prop-4-8" }

/-- ★★★★★**[FrdI] Proposition 4.8, (ii)** の `Ψ^birat : 𝒞^birat ⥤ 𝒞^birat`。

★★`section Birat` の中(`G` が section 変数)で宣言する ——
先行例 `Rmk451Birat.lean` の `biratFunctor` と同じ形。
★`map_id` は `biratNfMap_toHomBirat`(実装済)、`map_comp` は仮定として受ける。 -/
noncomputable def psiBiratNf
    (hcomp : ∀ {A B E : C} (f : HomBirat P G A B) (g : HomBirat P G B E),
      biratNfMap P d G A E (compBirat P G G.core f g)
        = compBirat P G G.core (biratNfMap P d G A B f) (biratNfMap P d G B E g)) :
    BiratCat P G ⥤ BiratCat P G where
  obj A := (nfObj P G.core d (biratDown P G A) : C)
  map {_ _} f := biratNfMap P d G _ _ f
  map_id A := by
    have h := biratNfMap_toHomBirat P d G (𝟙 (biratDown P G A))
    have hid : nfMap P G.core d (𝟙 (biratDown P G A)) = 𝟙 _ :=
      (naiveFrob P G.core d).map_id (biratDown P G A)
    rw [hid] at h
    exact h
  map_comp f g := hcomp f g


/-- ★★★★**1-可換図式** —— `𝒞 → 𝒞^birat` の四角形。

★対象は定義的に一致するので**成分は恒等**で済み、
自然性がちょうど `biratNfMap_toHomBirat` である。 -/
noncomputable def psiBiratNfSquare
    (hcomp : ∀ {A B E : C} (f : HomBirat P G A B) (g : HomBirat P G B E),
      biratNfMap P d G A E (compBirat P G G.core f g)
        = compBirat P G G.core (biratNfMap P d G A B f) (biratNfMap P d G B E g)) :
    toBiratCat P G ⋙ psiBiratNf P d G hcomp
      ≅ naiveFrob P G.core d ⋙ toBiratCat P G :=
  NatIso.ofComponents (fun _ => Iso.refl _) (fun {A B} φ => by
    show compBirat P G G.core (biratNfMap P d G A B (toHomBirat (P := P) (G := G) φ))
        (toHomBirat (P := P) (G := G) (𝟙 (nfObj P G.core d B)))
      = compBirat P G G.core (toHomBirat (P := P) (G := G) (𝟙 (nfObj P G.core d A)))
          (toHomBirat (P := P) (G := G) (nfMap P G.core d φ))
    rw [biratNfMap_toHomBirat]
    rw [show compBirat P G G.core (toHomBirat (P := P) (G := G) (nfMap P G.core d φ))
        (toHomBirat (P := P) (G := G) (𝟙 (nfObj P G.core d B)))
        = toHomBirat (P := P) (G := G) (nfMap P G.core d φ) from
      compBirat_id_right G.core _]
    rw [show compBirat P G G.core (toHomBirat (P := P) (G := G) (𝟙 (nfObj P G.core d A)))
        (toHomBirat (P := P) (G := G) (nfMap P G.core d φ))
        = toHomBirat (P := P) (G := G) (nfMap P G.core d φ) from
      compBirat_id_left G.core _])

def psiBiratNf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — Ψ^birat の構成と 1-可換図式",
    sectionId := "frdi-prop-4-8" }


/-! ## ★★★★★段 4 の `map_comp`

原文 (FrdI p.88):
> (ii), observe that the naive Frobenius functor [cf. Proposition 2.1] determines a

★★★**梃子は `Definition 1.3, (v), (a)`(`preStepMono`)—— pre-step は mono**。
共通の上界へ送ったあと **2 回 mono を消す**だけで等式が出る。
★`calc` は型合わせで詰まるので、**各段を独立の `have` にして `.trans` で繋ぐ**。 -/

/-- ★添字圏の射が満たす三角形(`Over.w` を `𝒞` の射の等式に落としたもの)。 -/
theorem idxBirat_w {A : C} {Z V : IdxBirat P G A} (u : Z ⟶ V) :
    u.unop.left.hom ≫ Z.unop.hom.hom = V.unop.hom.hom :=
  congrArg (fun t : V.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w u.unop)

/-- ★★★★★**`Ψ^birat` は合成を保つ**。 -/
theorem biratNfMap_compBirat {A B E : C} (f : HomBirat P G A B) (g : HomBirat P G B E) :
    biratNfMap P d G A E (compBirat P G G.core f g)
      = compBirat P G G.core (biratNfMap P d G A B f) (biratNfMap P d G B E g) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratNfMap_mk, biratNfMap_mk, biratNfMap_mk, compBirat_mk, nfMap_comp]
  refine HomBirat.sound
    (IsFiltered.max ((idxBiratNfMap P d G A).obj (biratPullIdx G.core Z φ W))
      (biratPullIdx G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)))
    (IsFiltered.leftToMax _ _) (IsFiltered.rightToMax _ _) ?_
  set u := IsFiltered.leftToMax ((idxBiratNfMap P d G A).obj (biratPullIdx G.core Z φ W))
    (biratPullIdx G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) with hu
  set v := IsFiltered.rightToMax ((idxBiratNfMap P d G A).obj (biratPullIdx G.core Z φ W))
    (biratPullIdx G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) with hv
  haveI hmZ : Mono (nfMap P G.core d Z.unop.hom.hom) :=
    G.core.preStepMono _ (nfMap_preStep P G.core d _ Z.unop.hom.property.2)
  haveI hmW : Mono (nfMap P G.core d W.unop.hom.hom) :=
    G.core.preStepMono _ (nfMap_preStep P G.core d _ W.unop.hom.property.2)
  have hwu := idxBirat_w P G u
  have hwv := idxBirat_w P G v
  -- ★段 1: `γ` の側。
  have hgamma : u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W) = v.unop.left.hom ≫ biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) := by
    refine (cancel_mono (nfMap P G.core d Z.unop.hom.hom)).mp ?_
    have e1 : (u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W)) ≫ nfMap P G.core d Z.unop.hom.hom
        = u.unop.left.hom ≫ (nfMap P G.core d (biratPullGamma G.core Z φ W) ≫ nfMap P G.core d Z.unop.hom.hom) := Category.assoc _ _ _
    have e2 : u.unop.left.hom ≫ (nfMap P G.core d (biratPullGamma G.core Z φ W) ≫ nfMap P G.core d Z.unop.hom.hom)
        = u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W ≫ Z.unop.hom.hom) :=
      congrArg (fun t => u.unop.left.hom ≫ t) (nfMap_comp P G.core d _ _).symm
    have e3 : u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W ≫ Z.unop.hom.hom)
        = v.unop.left.hom ≫ (biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d Z.unop.hom.hom) := hwu.trans hwv.symm
    have e4 : v.unop.left.hom ≫ (biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d Z.unop.hom.hom)
        = (v.unop.left.hom ≫ biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d Z.unop.hom.hom := (Category.assoc _ _ _).symm
    exact e1.trans (e2.trans (e3.trans e4))
  -- ★段 2: `α` の側。
  have halpha : u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W) = v.unop.left.hom ≫ biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) := by
    refine (cancel_mono (nfMap P G.core d W.unop.hom.hom)).mp ?_
    have f1 : (u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W)) ≫ nfMap P G.core d W.unop.hom.hom
        = u.unop.left.hom ≫ (nfMap P G.core d (biratPullAlpha G.core Z φ W) ≫ nfMap P G.core d W.unop.hom.hom) := Category.assoc _ _ _
    have f2 : u.unop.left.hom ≫ (nfMap P G.core d (biratPullAlpha G.core Z φ W) ≫ nfMap P G.core d W.unop.hom.hom)
        = u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W ≫ W.unop.hom.hom) :=
      congrArg (fun t => u.unop.left.hom ≫ t) (nfMap_comp P G.core d _ _).symm
    have f3 : u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W ≫ W.unop.hom.hom)
        = u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W ≫ φ) :=
      congrArg (fun t => u.unop.left.hom ≫ nfMap P G.core d t) (biratPull_sq G.core Z φ W).symm
    have f4 : u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W ≫ φ)
        = u.unop.left.hom ≫ (nfMap P G.core d (biratPullGamma G.core Z φ W) ≫ nfMap P G.core d φ) :=
      congrArg (fun t => u.unop.left.hom ≫ t) (nfMap_comp P G.core d _ _)
    have f5 : u.unop.left.hom ≫ (nfMap P G.core d (biratPullGamma G.core Z φ W) ≫ nfMap P G.core d φ)
        = (u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W)) ≫ nfMap P G.core d φ := (Category.assoc _ _ _).symm
    have f6 : (u.unop.left.hom ≫ nfMap P G.core d (biratPullGamma G.core Z φ W)) ≫ nfMap P G.core d φ = (v.unop.left.hom ≫ biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d φ :=
      congrArg (fun t => t ≫ nfMap P G.core d φ) hgamma
    have f7 : (v.unop.left.hom ≫ biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d φ = v.unop.left.hom ≫ (biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d φ) := Category.assoc _ _ _
    have f8 : v.unop.left.hom ≫ (biratPullGamma G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d φ)
        = v.unop.left.hom ≫ (biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d W.unop.hom.hom) :=
      congrArg (fun t => v.unop.left.hom ≫ t)
        (biratPull_sq G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W))
    have f9 : v.unop.left.hom ≫ (biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d W.unop.hom.hom)
        = (v.unop.left.hom ≫ biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d W.unop.hom.hom := (Category.assoc _ _ _).symm
    exact f1.trans (f2.trans (f3.trans (f4.trans (f5.trans
      (f6.trans (f7.trans (f8.trans f9)))))))
  have g1 : u.unop.left.hom ≫ (nfMap P G.core d (biratPullAlpha G.core Z φ W) ≫ nfMap P G.core d ψ)
      = (u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W)) ≫ nfMap P G.core d ψ := (Category.assoc _ _ _).symm
  have g2 : (u.unop.left.hom ≫ nfMap P G.core d (biratPullAlpha G.core Z φ W)) ≫ nfMap P G.core d ψ = (v.unop.left.hom ≫ biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d ψ :=
    congrArg (fun t => t ≫ nfMap P G.core d ψ) halpha
  have g3 : (v.unop.left.hom ≫ biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W)) ≫ nfMap P G.core d ψ = v.unop.left.hom ≫ (biratPullAlpha G.core ((idxBiratNfMap P d G A).obj Z) (nfMap P G.core d φ) ((idxBiratNfMap P d G B).obj W) ≫ nfMap P G.core d ψ) := Category.assoc _ _ _
  exact g1.trans (g2.trans g3)

def biratNfMap_compBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — Ψ^birat は合成を保つ",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★★★`naive-frob-birat` 完了 —— 仮定なしの `Ψ^birat` -/

/-- ★★★★★**[FrdI] Proposition 4.8, (ii)** の `Ψ^birat : 𝒞^birat ⥤ 𝒞^birat`(**仮定なし**)。

★`map_comp` が `biratNfMap_compBirat` で証明できたので、仮定を外した版。 -/
noncomputable def psiBirat : BiratCat P G ⥤ BiratCat P G :=
  psiBiratNf P d G (fun f g => biratNfMap_compBirat P d G f g)

/-- ★★★★**1-可換図式**(仮定なし版)。 -/
noncomputable def psiBiratSquare :
    toBiratCat P G ⋙ psiBirat P d G ≅ naiveFrob P G.core d ⋙ toBiratCat P G :=
  psiBiratNfSquare P d G (fun f g => biratNfMap_compBirat P d G f g)

def psiBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — Ψ^birat(仮定なし)と 1-可換図式",
    sectionId := "frdi-prop-4-8" }


/-! ## ★★★`psi-birat-equiv` —— `Ψ^birat` が圏同値であること

原文 (FrdI p.88):
> (ii), observe that the naive Frobenius functor [cf. Proposition 2.1] determines a

★★在庫の `prop_2_1_iii_mp`(perfect 型なら `naiveFrob` は圏同値)を
`𝒞^birat` へ移す。★対象は `𝒞^birat` が `𝒞` と同じ型なので、
**本質的全射性はそのまま移る**。 -/

/-- ★★★★**`Ψ^birat` は本質的全射** —— 対象が同じ型なので `naiveFrob` の
本質的全射性を `toBiratCat` で送るだけ。 -/
theorem psiBirat_essSurj (hpt : IsOfPerfectType P) : (psiBirat P d G).EssSurj := by
  haveI := naiveFrob_essSurj P G.core d hpt
  constructor
  intro B
  obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := naiveFrob P G.core d)
    (biratDown P G B)
  exact ⟨A, ⟨(toBiratCat P G).mapIso ε⟩⟩

def psiBirat_essSurj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — Ψ^birat は本質的全射",
    sectionId := "frdi-prop-4-8" }


end Birat

def biratNfMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — naive Frobenius を 𝒞^birat の射へ降ろす",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★`Ψ^birat` の組み上げ —— universe の測定(2026-08-19)

原文 (FrdI p.88):
> (ii), observe that the naive Frobenius functor [cf. Proposition 2.1] determines a

★★`BiratCat P G ⥤ BiratCat P G` を書こうとすると
**`Category.{v2, u2} (BiratCat P G)` の合成に失敗する** ——
`biratCategory` が与えるのは `Category.{max v u2 v2}` だからである。

★`CategoryTheory.Functor.{max v u2 v2, max v u2 v2, u2, u2}` と明示しても
`obj` フィールドの位置で同じ要求が出る(2 回試して同じ)。

★★★**先行例 `Rmk451Birat.lean` の `biratFunctor` は通っている**ので、
そちらの宣言の形(`variable` の並び・`section` の切り方)を写すのが次の手である。
★★★**3 回目の試行で本当の原因が見えた(2026-08-19)**:
`biratNfObj (G : Frobenioid P) (A : BiratCat P G) : BiratCat P G` と
**対象の写像に名前を付けた瞬間**、同じ `Category.{v2, u2}` 要求が出た。
★同時に「`biratNfObj P d` の `d` が `Frobenioid P`(sort `Prop`)を期待される」
というエラーも出た —— ★★**`d` は statement に現れないので `variable` の
自動包含に入らない**(`include hPB hPB' in` と同じ罠)。

★★したがって手は 2 つ:
1. `d` を `include d in` で明示的に包含する
2. `BiratCat` の universe を `Rmk451Birat.lean` と同じ形に揃える
   —— そちらは `universe v u w v2 u2`(**`v2` が `u2` の前**)で宣言している。
   本ファイルは `universe v u w u2 v2` である。★この順序差を先に潰すこと。

★★本ファイルは**ここまでを緑のまま**にしておく ——
並行セッション(ABC3b)の `git add -A` が赤い状態を巻き込むため。 -/

end ABC3.Found.FrdI
