/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop21
import ABC3.Found.FrdI.Prop44
import ABC3.Found.FrdI.Thm34Quasi

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

/-! ## ★★★`naiveFrob` の反射 —— `full` / `faithful` の下ごしらえ

原文 (FrdI p.88):
> (ii), observe that the naive Frobenius functor [cf. Proposition 2.1] determines a

★★添字圏 `IdxBirat` は co-angular pre-step のなす圏なので、
`Ψ^birat` の `full` / `faithful` には
**`naiveFrob` が `coaPreProp` を反射する**ことが要る。

★★★材料はすべて在庫だった:
- pre-step の反射 —— `nfMap_preStep_of`(`Prop21.lean`)
- co-angular の反射 —— `isCoAngular_of_map`(`Thm34Quasi.lean`、本セッション)
  が要求する 4 つの保存も、`prop_1_10_i_*_of` を `nfMap` の四角形に当てるだけ。 -/

/-- ★**`naiveFrob` は linear を保つ**。 -/
theorem nfMap_linear {A B : C} (φ : A ⟶ B) (h : IsLinear P φ) :
    IsLinear P (nfMap P F d φ) :=
  prop_1_10_i_linear_of P (by rw [nfHom_degFr, nfHom_degFr]) (nfMap_sq P F d φ) h

/-- ★**`naiveFrob` は isometric を保つ** —— Frobenius 型射は isometric だから。 -/
theorem nfMap_isometric {A B : C} (φ : A ⟶ B) (h : IsIsometric P φ) :
    IsIsometric P (nfMap P F d φ) :=
  prop_1_10_i_isometric_of P (nfHom_frobType P F d A).1.2 (nfHom_frobType P F d B).1.2
    (nfMap_sq P F d φ) h

/-- ★**`naiveFrob` は base-isomorphism を保つ**。 -/
theorem nfMap_baseIso {A B : C} (φ : A ⟶ B) (h : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism P (nfMap P F d φ) :=
  prop_1_10_i_baseIso_of P (nfHom_frobType P F d A).2 (nfHom_frobType P F d B).2
    (nfMap_sq P F d φ) h

/-- ★★★★**`naiveFrob` は co-angular を反射する**(perfect 型)。

★`isCoAngular_of_map` に 4 つの保存を渡すだけ。 -/
theorem nfMap_coAngular_of (hpt : IsOfPerfectType P) {A B : C} (φ : A ⟶ B)
    (h : IsCoAngular P (nfMap P F d φ)) : IsCoAngular P φ := by
  haveI := prop_2_1_iii_mp P F d hpt
  exact isCoAngular_of_map (P₁ := P) (P₂ := P) (naiveFrob P F d)
    (fun f hf => nfMap_linear P F d f hf) (fun f hf => nfMap_isometric P F d f hf)
    (fun f hf => nfMap_preStep P F d f hf) (fun f hf => nfMap_baseIso P F d f hf) φ h

/-- ★★★★★**`naiveFrob` は `coaPreProp` を反射する**(perfect 型)。

★★これで `𝒞^{coa-pre}` の間の関手が**充満忠実**になる下ごしらえが揃う。 -/
theorem nfMap_coaPreProp_of (hpt : IsOfPerfectType P) {A B : C} (φ : A ⟶ B)
    (h : coaPreProp P (nfMap P F d φ)) : coaPreProp P φ :=
  ⟨nfMap_coAngular_of P F d hpt φ h.1, nfMap_preStep_of P F d φ h.2⟩

def nfMap_coaPreProp_of.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — naive Frobenius は co-angular pre-step を反射する",
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

/-! ## ★★★★`𝒞^{coa-pre}` の圏同値

★★保存(`nfMap_coaPreProp`)と反射(`nfMap_coaPreProp_of`)が揃ったので、
`coaPreNfMap` は**充満忠実**になる。★本質的全射性は同型を包むだけ。 -/

/-- ★**同型は `𝒞^{coa-pre}` の同型に包める** —— 同型は co-angular pre-step だから。 -/
def coaPreIsoOfIso {X Y : C} (θ : X ≅ Y) :
    (⟨X⟩ : CoaPre P G) ≅ (⟨Y⟩ : CoaPre P G) where
  hom := ⟨θ.hom, isCoAngular_of_isIso P θ.hom, isPreStep_of_isIso P θ.hom⟩
  inv := ⟨θ.inv, isCoAngular_of_isIso P θ.inv, isPreStep_of_isIso P θ.inv⟩
  hom_inv_id := WideSubcategory.hom_ext _ θ.hom_inv_id
  inv_hom_id := WideSubcategory.hom_ext _ θ.inv_hom_id

/-- ★★**忠実** —— `naiveFrob` の忠実性そのもの。 -/
theorem coaPreNfMap_faithful (hpt : IsOfPerfectType P) :
    (coaPreNfMap P d G).Faithful where
  map_injective {_ _ f g} h := by
    haveI := naiveFrob_faithful P G.core d hpt
    refine WideSubcategory.hom_ext _ ((naiveFrob P G.core d).map_injective ?_)
    exact congrArg (fun t : (coaPreNfMap P d G).obj _ ⟶ _ => t.hom) h

/-- ★★★**充満** —— `naiveFrob` の充満性 ＋ `coaPreProp` の**反射**。 -/
theorem coaPreNfMap_full (hpt : IsOfPerfectType P) : (coaPreNfMap P d G).Full where
  map_surjective {X Y} g := by
    haveI := naiveFrob_full P G.core d hpt
    obtain ⟨f₀, hf₀⟩ := (naiveFrob P G.core d).map_surjective g.hom
    refine ⟨⟨f₀, nfMap_coaPreProp_of P G.core d hpt f₀ ?_⟩,
      WideSubcategory.hom_ext _ hf₀⟩
    show coaPreProp P (nfMap P G.core d f₀)
    rw [show nfMap P G.core d f₀ = g.hom from hf₀]
    exact g.property

/-- ★★**本質的全射** —— `naiveFrob` の本質的全射性を同型ごと包む。 -/
theorem coaPreNfMap_essSurj (hpt : IsOfPerfectType P) :
    (coaPreNfMap P d G).EssSurj where
  mem_essImage Y := by
    haveI := naiveFrob_essSurj P G.core d hpt
    obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := naiveFrob P G.core d) Y.obj
    exact ⟨(⟨A⟩ : CoaPre P G), ⟨coaPreIsoOfIso P G ε⟩⟩

/-- ★★★★★**`𝒞^{coa-pre}` の間の関手は圏同値**(perfect 型)。 -/
theorem coaPreNfMap_isEquivalence (hpt : IsOfPerfectType P) :
    (coaPreNfMap P d G).IsEquivalence where
  faithful := coaPreNfMap_faithful P d G hpt
  full := coaPreNfMap_full P d G hpt
  essSurj := coaPreNfMap_essSurj P d G hpt

/-- ★★★★**添字圏の写像も圏同値** —— `Over.post` と `.op` が圏同値を保つ。 -/
theorem idxBiratNfMap_isEquivalence (hpt : IsOfPerfectType P) (A : C) :
    (idxBiratNfMap P d G A).IsEquivalence := by
  haveI := coaPreNfMap_isEquivalence P d G hpt
  exact inferInstanceAs
    (Functor.IsEquivalence (Over.post (X := coaPreObj P G A) (coaPreNfMap P d G)).op)

def coaPreNfMap_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (ii) — 𝒞^{coa-pre} の間の関手は圏同値",
    sectionId := "frdi-prop-4-8" }

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

/-- ★★★★★**`Ψ^birat` は忠実**。

★★添字圏が圏同値なので、共通の上界を**像の中に取り直せる** ——
そのあとは `naiveFrob` の忠実性 1 本。 -/
theorem psiBirat_faithful (hpt : IsOfPerfectType P) : (psiBirat P d G).Faithful where
  map_injective {X Y} {x y} h := by
    haveI := idxBiratNfMap_isEquivalence P d G hpt (biratDown P G X)
    haveI := naiveFrob_faithful P G.core d hpt
    obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep x
    obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep y
    have h' : HomBirat.mk ((idxBiratNfMap P d G (biratDown P G X)).obj Z) (nfMap P G.core d φ)
        = HomBirat.mk ((idxBiratNfMap P d G (biratDown P G X)).obj W) (nfMap P G.core d ψ) :=
      (biratNfMap_mk P d G Z φ).symm.trans (h.trans (biratNfMap_mk P d G W ψ))
    obtain ⟨V, u, v, hV⟩ := HomBirat.eq_iff.mp h'
    obtain ⟨V₀, ⟨θ⟩⟩ := Functor.EssSurj.mem_essImage (F := (idxBiratNfMap P d G (biratDown P G X))) V
    obtain ⟨u₀, hu₀⟩ := (idxBiratNfMap P d G (biratDown P G X)).map_surjective (u ≫ θ.inv)
    obtain ⟨v₀, hv₀⟩ := (idxBiratNfMap P d G (biratDown P G X)).map_surjective (v ≫ θ.inv)
    have hu : nfMap P G.core d u₀.unop.left.hom = θ.inv.unop.left.hom ≫ u.unop.left.hom :=
      congrArg (fun t : (idxBiratNfMap P d G (biratDown P G X)).obj Z ⟶ (idxBiratNfMap P d G (biratDown P G X)).obj V₀ => t.unop.left.hom) hu₀
    have hv : nfMap P G.core d v₀.unop.left.hom = θ.inv.unop.left.hom ≫ v.unop.left.hom :=
      congrArg (fun t : (idxBiratNfMap P d G (biratDown P G X)).obj W ⟶ (idxBiratNfMap P d G (biratDown P G X)).obj V₀ => t.unop.left.hom) hv₀
    refine HomBirat.eq_iff.mpr ⟨V₀, u₀, v₀, (naiveFrob P G.core d).map_injective ?_⟩
    show nfMap P G.core d (u₀.unop.left.hom ≫ φ) = nfMap P G.core d (v₀.unop.left.hom ≫ ψ)
    have e3 : nfMap P G.core d (u₀.unop.left.hom ≫ φ)
        = nfMap P G.core d u₀.unop.left.hom ≫ nfMap P G.core d φ := nfMap_comp P G.core d _ _
    have e4 : nfMap P G.core d (v₀.unop.left.hom ≫ ψ)
        = nfMap P G.core d v₀.unop.left.hom ≫ nfMap P G.core d ψ := nfMap_comp P G.core d _ _
    have g1 : nfMap P G.core d u₀.unop.left.hom ≫ nfMap P G.core d φ
        = (θ.inv.unop.left.hom ≫ u.unop.left.hom) ≫ nfMap P G.core d φ :=
      congrArg (fun t => t ≫ nfMap P G.core d φ) hu
    have g2 : (θ.inv.unop.left.hom ≫ u.unop.left.hom) ≫ nfMap P G.core d φ
        = θ.inv.unop.left.hom ≫ (u.unop.left.hom ≫ nfMap P G.core d φ) := Category.assoc _ _ _
    have g3 : θ.inv.unop.left.hom ≫ (u.unop.left.hom ≫ nfMap P G.core d φ)
        = θ.inv.unop.left.hom ≫ (v.unop.left.hom ≫ nfMap P G.core d ψ) :=
      congrArg (fun t => θ.inv.unop.left.hom ≫ t) hV
    have g4 : θ.inv.unop.left.hom ≫ (v.unop.left.hom ≫ nfMap P G.core d ψ)
        = (θ.inv.unop.left.hom ≫ v.unop.left.hom) ≫ nfMap P G.core d ψ := (Category.assoc _ _ _).symm
    have g5 : (θ.inv.unop.left.hom ≫ v.unop.left.hom) ≫ nfMap P G.core d ψ
        = nfMap P G.core d v₀.unop.left.hom ≫ nfMap P G.core d ψ :=
      congrArg (fun t => t ≫ nfMap P G.core d ψ) hv.symm
    exact e3.trans (g1.trans (g2.trans (g3.trans (g4.trans (g5.trans e4.symm)))))

/-- ★★★★★**`Ψ^birat` は充満**。

★添字を像の中へ取り直し、`naiveFrob` の充満性で射を持ち上げる。 -/
theorem psiBirat_full (hpt : IsOfPerfectType P) : (psiBirat P d G).Full where
  map_surjective {X Y} h := by
    haveI := idxBiratNfMap_isEquivalence P d G hpt (biratDown P G X)
    haveI := naiveFrob_full P G.core d hpt
    obtain ⟨V, χ, rfl⟩ := HomBirat.exists_rep h
    obtain ⟨V₀, ⟨θ⟩⟩ := Functor.EssSurj.mem_essImage (F := (idxBiratNfMap P d G (biratDown P G X))) V
    obtain ⟨χ₀, hχ₀⟩ := (naiveFrob P G.core d).map_surjective (θ.inv.unop.left.hom ≫ χ)
    refine ⟨HomBirat.mk V₀ χ₀, ?_⟩
    refine (biratNfMap_mk P d G V₀ χ₀).trans ?_
    have e2 : nfMap P G.core d χ₀ = θ.inv.unop.left.hom ≫ χ := hχ₀
    exact (congrArg (HomBirat.mk ((idxBiratNfMap P d G (biratDown P G X)).obj V₀)) e2).trans (HomBirat.mk_map θ.inv χ)

/-- ★★★★★**[FrdI] Proposition 4.8, (ii)** —— `Ψ^birat` は**圏同値**。 -/
theorem psiBirat_isEquivalence (hpt : IsOfPerfectType P) :
    (psiBirat P d G).IsEquivalence where
  faithful := psiBirat_faithful P d G hpt
  full := psiBirat_full P d G hpt
  essSurj := psiBirat_essSurj P d G hpt

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

/-! ## ★★★★★`Proposition 4.8, (ii)` —— 条として揃った

原文 (FrdI p.88):
> (ii) If C is of perfect and isotropic type, then so is Cbirat.

★実装した部品:

| 主張 | 宣言 |
|---|---|
| naive Frobenius が co-angular pre-step を保つ／反射する | `nfMap_coaPreProp` / `nfMap_coaPreProp_of` |
| `𝒞^{coa-pre}` と添字圏の圏同値 | `coaPreNfMap_isEquivalence` / `idxBiratNfMap_isEquivalence` |
| `Ψ^birat` の構成と 1-可換図式 | `psiBirat` / `psiBiratSquare` |
| `Ψ^birat` が圏同値 | `psiBirat_isEquivalence` |

★★★原文が「the naive Frobenius functor determines a natural 1-commutative diagram」
と一行で書く部分が、余極限の移送 4 段＋圏同値 3 性質になった。 -/
def prop_4_8_ii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (ii)",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
