/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34EndBs
import ABC3.Found.FrdI.Def31Pf

/-!
# [FrdI] Theorem 3.4, (iii) —— `Ψ^pf` の構成

原文 (FrdI p.64):
> Definition 3.1, (iii)] that we obtain a 1-unique 1-commutative diagram as in the

## ★★測ったこと(2026-08-19)

`𝒞^pf` は `Def31Pf.lean` で**余極限**として構成されている:

- 対象は `𝒞` の対象そのもの(`PfCat _P _F := C`)
- 射は添字圏 `IdxPf P F A B = Under (biFrObj P F A B)` 上の余極限
- 添字圏の底は `BiFr = WideSubcategory (biFrPropOf P F)`、すなわち
  「同次数の Frobenius 型射の対」がなす `𝒞 × 𝒞` の広い部分圏

★★★**`iv-endbsiso` の `plBkMap` と同じ形**である ——
射のクラスを保つ関手が広い部分圏の間の関手を誘導し、
`Under.post` / `Over.post` でスライスへ降りる。
★本ファイルはその第 1 段(`biFrMap`)を取る。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (F₁ : FrobenioidCore P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (F₂ : FrobenioidCore P₂)

/-! ## ★段 1 —— `𝒞^{bi-Fr}` の間の関手 -/

/-- ★★★★**Frobenius 型と次数の一致を保つ関手は `𝒞^{bi-Fr}` の間の関手を誘導する**。

★`biFrProp` は「両成分が Frobenius 型」かつ「次数が一致」なので、
`Ψ` について要るのはちょうどその 2 つ。
★★これは `iv-endbsiso` の `plBkMap` と**同じ形**である。 -/
def biFrMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) :
    BiFr P₁ F₁ ⥤ BiFr P₂ F₂ where
  obj X := ⟨(Ψ.obj X.obj.1, Ψ.obj X.obj.2)⟩
  map f := ⟨(Ψ.map f.hom.1, Ψ.map f.hom.2),
    hFT _ f.property.1, hFT _ f.property.2.1, hdegEq _ _ f.property.2.2⟩
  map_id X := WideSubcategory.hom_ext _ (Prod.ext (Ψ.map_id X.obj.1) (Ψ.map_id X.obj.2))
  map_comp f g := WideSubcategory.hom_ext _
    (Prod.ext (Ψ.map_comp f.hom.1 g.hom.1) (Ψ.map_comp f.hom.2 g.hom.2))

@[simp]
theorem biFrMap_obj (Ψ : C₁ ⥤ C₂) (hFT) (hdegEq) (X : BiFr P₁ F₁) :
    (biFrMap F₁ F₂ Ψ hFT hdegEq).obj X = ⟨(Ψ.obj X.obj.1, Ψ.obj X.obj.2)⟩ := rfl

/-- ★★**添字圏の底が対応する** —— `biFrObj` は `biFrMap` で `biFrObj` に写る。 -/
theorem biFrMap_biFrObj (Ψ : C₁ ⥤ C₂) (hFT) (hdegEq) (A B : C₁) :
    (biFrMap F₁ F₂ Ψ hFT hdegEq).obj (biFrObj P₁ F₁ A B)
      = biFrObj P₂ F₂ (Ψ.obj A) (Ψ.obj B) := rfl

/-- ★★★★**添字圏の間の関手** —— `Under.post` で降ろすだけ。

★`IdxPf P F A B = Under (biFrObj P F A B)` なので、
`biFrMap` を `Under.post` に通せばそのまま添字圏の写像になる。 -/
def idxPfMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) (A B : C₁) :
    IdxPf P₁ F₁ A B ⥤ IdxPf P₂ F₂ (Ψ.obj A) (Ψ.obj B) :=
  Under.post (X := biFrObj P₁ F₁ A B) (biFrMap F₁ F₂ Ψ hFT hdegEq)

def biFrMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf の添字圏の写像",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★段 2 の要 —— `Ψ` は遷移写像と可換

原文 (FrdI p.64):
> Definition 3.1, (iii)] that we obtain a 1-unique 1-commutative diagram as in the

★★`idxTransport` は「四角形を可換にする**唯一の**射」として定まる
(`frobTransport_eq`)。`Ψ` を四角形に当てれば同じ形の四角形になるので、
**一意性 1 本**で可換性が出る。
★★★これが段 2(余錐)の naturality そのものである。 -/

/-- ★★★★★**`Ψ` は遷移写像と可換**。

★`Ψ` を四角形 `φ ≫ β = α ≫ φ′` に当てると
`Ψφ ≫ Ψβ = Ψα ≫ Ψφ′` になり、`frobTransport_eq` の一意性が効く。 -/
theorem idxTransport_map (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {A B : C₁} {Z W : IdxPf P₁ F₁ A B} (u : Z ⟶ W)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    idxTransport P₂ F₂ ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).map u) (Ψ.map φ)
      = Ψ.map (idxTransport P₁ F₁ u φ) := by
  refine frobTransport_eq _ _ _ _ _ _ _ ?_
  have hsq : φ ≫ u.right.hom.2 = u.right.hom.1 ≫ idxTransport P₁ F₁ u φ :=
    idxTransport_spec u φ
  have hmap := congrArg (fun f : Z.right.obj.1 ⟶ W.right.obj.2 => Ψ.map f) hsq
  simp only [Ψ.map_comp] at hmap
  exact hmap

def idxTransport_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ は遷移写像と可換",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★段 2 —— 余錐と、射の写像 -/

/-- ★★★★**余錐** —— 各添字で `Ψ` を当て、添字は `idxPfMap` で送る。

★naturality はちょうど `idxTransport_map` ＋ `HomPf.mk_map` の 2 本。 -/
noncomputable def homPfCocone (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    (A B : C₁) : Limits.Cocone (homFunctorPf P₁ F₁ A B) where
  pt := HomPf P₂ F₂ (Ψ.obj A) (Ψ.obj B)
  ι :=
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj Z) (Ψ.map φ.down)
      naturality := by
        intro Z W u
        ext φ
        show HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj W)
            (Ψ.map (idxTransport P₁ F₁ u φ.down))
          = HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj Z) (Ψ.map φ.down)
        rw [← idxTransport_map F₁ F₂ Ψ hFT hdegEq u φ.down]
        exact HomPf.mk_map ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).map u) (Ψ.map φ.down) }

/-- ★★★★★**射の写像** —— 余極限の普遍性で降ろす。 -/
noncomputable def homPfMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    (A B : C₁) : HomPf P₁ F₁ A B → HomPf P₂ F₂ (Ψ.obj A) (Ψ.obj B) :=
  (Limits.colimit.desc (homFunctorPf P₁ F₁ A B) (homPfCocone F₁ F₂ Ψ hFT hdegEq A B)).hom

/-- ★**代表元での計算則** —— `mk` は `mk` に写る。 -/
@[simp] theorem homPfMap_mk (Ψ : C₁ ⥤ C₂) (hFT) (hdegEq) {A B : C₁}
    (Z : IdxPf P₁ F₁ A B) (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    homPfMap F₁ F₂ Ψ hFT hdegEq A B (HomPf.mk Z φ)
      = HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj Z) (Ψ.map φ) :=
  congrFun (congrArg (fun t => t.hom)
    (Limits.colimit.ι_desc (homPfCocone F₁ F₂ Ψ hFT hdegEq A B) Z)) (ULift.up φ)

def homPfMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf の射の写像",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★段 3 —— 関手則 -/

/-- ★★**添字圏の始対象は始対象に写る** —— `Ψ.map 𝟙 = 𝟙` の分だけの話。 -/
theorem idxPfMap_idxOne (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    (A B : C₁) :
    (idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj (idxOne P₁ F₁ A B)
      = idxOne P₂ F₂ (Ψ.obj A) (Ψ.obj B) :=
  congrArg Under.mk (WideSubcategory.hom_ext _ (Prod.ext (Ψ.map_id A) (Ψ.map_id B)))

/-- ★★★★**`𝒞 → 𝒞^pf` と可換** —— `toHomPf` は `toHomPf` に写る。

★これがそのまま `map_id`(`𝟙 = toHomPf 𝟙`)を与える。 -/
theorem homPfMap_toHomPf (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {A B : C₁} (φ : A ⟶ B) :
    homPfMap F₁ F₂ Ψ hFT hdegEq A B (toHomPf (F := F₁) φ)
      = toHomPf (F := F₂) (Ψ.map φ) := by
  -- ★対象の等式で `rw` すると motive が壊れるので、
  --   **添字圏の射を 1 本作って `HomPf.mk_map` で移す**。
  let hv : idxOne P₂ F₂ (Ψ.obj A) (Ψ.obj B)
      ⟶ (idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj (idxOne P₁ F₁ A B) :=
    Under.homMk (𝟙 _) (by
      refine (Category.comp_id _).trans ?_
      exact (WideSubcategory.hom_ext _ (Prod.ext (Ψ.map_id A) (Ψ.map_id B))).symm)
  have ht : idxTransport P₂ F₂ hv (Ψ.map φ) = Ψ.map φ :=
    frobTransport_eq _ _ _ _ _ _ _ (by
      show Ψ.map φ ≫ 𝟙 _ = 𝟙 _ ≫ Ψ.map φ
      rw [Category.comp_id, Category.id_comp])
  calc homPfMap F₁ F₂ Ψ hFT hdegEq A B (HomPf.mk (idxOne P₁ F₁ A B) φ)
      = HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj (idxOne P₁ F₁ A B)) (Ψ.map φ) :=
        homPfMap_mk F₁ F₂ Ψ hFT hdegEq (idxOne P₁ F₁ A B) φ
    _ = HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj (idxOne P₁ F₁ A B))
          (idxTransport P₂ F₂ hv (Ψ.map φ)) :=
        congrArg (HomPf.mk ((idxPfMap F₁ F₂ Ψ hFT hdegEq A B).obj (idxOne P₁ F₁ A B))) ht.symm
    _ = HomPf.mk (idxOne P₂ F₂ (Ψ.obj A) (Ψ.obj B)) (Ψ.map φ) := HomPf.mk_map hv (Ψ.map φ)

def homPfMap_toHomPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf は 𝒞 → 𝒞^pf と可換",
    sectionId := "frdi-thm-3-4" }

/-! ### ★★3 脚版 —— `map_comp` のために -/

/-- ★★★**`𝒞^{tri-Fr}` の間の関手** —— `biFrMap` の 3 脚版。 -/
def triFrMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) :
    TriFr P₁ F₁ ⥤ TriFr P₂ F₂ where
  obj X := ⟨(Ψ.obj X.obj.1, Ψ.obj X.obj.2.1, Ψ.obj X.obj.2.2)⟩
  map f := ⟨(Ψ.map f.hom.1, Ψ.map f.hom.2.1, Ψ.map f.hom.2.2),
    hFT _ f.property.1, hFT _ f.property.2.1, hFT _ f.property.2.2.1,
    hdegEq _ _ f.property.2.2.2.1, hdegEq _ _ f.property.2.2.2.2⟩
  map_id X := WideSubcategory.hom_ext _
    (Prod.ext (Ψ.map_id X.obj.1) (Prod.ext (Ψ.map_id X.obj.2.1) (Ψ.map_id X.obj.2.2)))
  map_comp f g := WideSubcategory.hom_ext _
    (Prod.ext (Ψ.map_comp f.hom.1 g.hom.1)
      (Prod.ext (Ψ.map_comp f.hom.2.1 g.hom.2.1) (Ψ.map_comp f.hom.2.2 g.hom.2.2)))

/-- ★★★**3 つ組の添字圏の間の関手**。 -/
def idxPf3Map (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    (A B E : C₁) :
    IdxPf3 P₁ F₁ A B E ⥤ IdxPf3 P₂ F₂ (Ψ.obj A) (Ψ.obj B) (Ψ.obj E) :=
  Under.post (X := triFrObj P₁ F₁ A B E) (triFrMap F₁ F₂ Ψ hFT hdegEq)

/-- ★★★★★**`map_comp`** —— 3 脚の共通代表元を取って `compPf_mk` を 2 回当てる。

★★`idx12`/`idx23`/`idx13` はすべて `Under.post` なので、
`biFrMap` との可換性は**成分ごとに `rfl`** である。 -/
theorem homPfMap_compPf (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {A B E : C₁} (f : HomPf P₁ F₁ A B) (g : HomPf P₁ F₁ B E) :
    homPfMap F₁ F₂ Ψ hFT hdegEq A E (compPf P₁ F₁ f g)
      = compPf P₂ F₂ (homPfMap F₁ F₂ Ψ hFT hdegEq A B f)
          (homPfMap F₁ F₂ Ψ hFT hdegEq B E g) := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 f g
  rw [compPf_mk V φ ψ, homPfMap_mk, homPfMap_mk, homPfMap_mk, Ψ.map_comp]
  exact (compPf_mk ((idxPf3Map F₁ F₂ Ψ hFT hdegEq A B E).obj V) (Ψ.map φ) (Ψ.map ψ)).symm

def homPfMap_compPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf は合成を保つ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★`Ψ^pf` —— 組み上げ

原文 (FrdI p.64):
> Definition 3.1, (iii)] that we obtain a 1-unique 1-commutative diagram as in the
-/

/-- ★★★★★**[FrdI] Theorem 3.4, (iii)** —— `Ψ^pf : 𝒞₁^pf ⥤ 𝒞₂^pf`。

★対象は `𝒞^pf` が `𝒞` と同じ型なのでそのまま、射は `homPfMap`。
★関手則は `homPfMap_toHomPf`(単位)と `homPfMap_compPf`(合成)。 -/
noncomputable def psiPf (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) :
    PfCat P₁ F₁ ⥤ PfCat P₂ F₂ where
  obj A := Ψ.obj (pfDown P₁ F₁ A)
  map {A B} f := homPfMap F₁ F₂ Ψ hFT hdegEq (pfDown P₁ F₁ A) (pfDown P₁ F₁ B) f
  map_id A := by
    have h := homPfMap_toHomPf F₁ F₂ Ψ hFT hdegEq (𝟙 (pfDown P₁ F₁ A))
    rw [Ψ.map_id] at h
    exact h
  map_comp f g := homPfMap_compPf F₁ F₂ Ψ hFT hdegEq f g

/-- ★★★★**1-可換図式** —— `𝒞 → 𝒞^pf` の四角形。

★★対象は定義的に一致するので**成分は恒等**で済み、
自然性がちょうど `homPfMap_toHomPf` である。
★原文の「1-commutative」は自然同型のことなので、これがその形。 -/
noncomputable def psiPfSquare (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) :
    toPfCat P₁ F₁ ⋙ psiPf F₁ F₂ Ψ hFT hdegEq ≅ Ψ ⋙ toPfCat P₂ F₂ :=
  NatIso.ofComponents (fun _ => Iso.refl _) (fun {A B} φ => by
    show compPf P₂ F₂ (homPfMap F₁ F₂ Ψ hFT hdegEq A B (toHomPf (F := F₁) φ))
        (toHomPf (F := F₂) (𝟙 (Ψ.obj B)))
      = compPf P₂ F₂ (toHomPf (F := F₂) (𝟙 (Ψ.obj A))) (toHomPf (F := F₂) (Ψ.map φ))
    rw [homPfMap_toHomPf]
    rw [show compPf P₂ F₂ (toHomPf (F := F₂) (Ψ.map φ)) (toHomPf (F := F₂) (𝟙 (Ψ.obj B)))
        = toHomPf (F := F₂) (Ψ.map φ) from compPf_id_right _]
    rw [show compPf P₂ F₂ (toHomPf (F := F₂) (𝟙 (Ψ.obj A))) (toHomPf (F := F₂) (Ψ.map φ))
        = toHomPf (F := F₂) (Ψ.map φ) from compPf_id_left _])

def psiPf.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf の構成と 1-可換図式",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★`Theorem 3.4, (iii)` —— 条として揃った

原文 (FrdI p.62):
> (iii) Suppose that: (a) C1, C2 are of standard type; (b) if C1, C2 are of group-

★実装した部品(場所は `Thm34Quasi.lean` / `Thm34Pre.lean` / `Thm34EndBs.lean` / 本ファイル):

| 主張 | 宣言 |
|---|---|
| `Ψ_{N≥1}` の存在((F1)(F2)(F3)＋連結性) | `exists_admissible_F1` / `admissible_of_hom` / `degFr_order_preserve` / `admissible_all_of_connected` |
| 8 つの射のクラスの保存 | `isPrimeFrobComposite_map` / `isFrobeniusType_map_of_primeFrob` / `isLinear_map_iff_of_degMap` / `isBaseIsomorphism_map_of_classes` / `isCoAngular_map_of_equiv` / `isPullBack_map_of_equiv` / `isIsometric_map_of_classes` / `isLBInvertible_map_of_classes` |
| 非 group-like なら `Ψ_{N≥1}` は恒等 | `psiN_eq_id_of_orderPreserve` |
| `Ψ^pf` の構成 | `psiPf`(本ファイル) |
| 1-可換図式 | `psiPfSquare`(本ファイル) |
| 図式の合成関手は rigid | `isRigidFunctor_of_equivalence_comp` / `thm_3_4_iii_rigid` |

★★★**`Ψ^pf` の構成が最後の 1 本**だった。余極限の移送は
「`Ψ` が遷移写像と可換」(`idxTransport_map`)1 本に帰着し、
それは `frobTransport_eq` の**一意性**から出た。 -/
def thm_3_4_iii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62, item := "Theorem 3.4, (iii)",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
