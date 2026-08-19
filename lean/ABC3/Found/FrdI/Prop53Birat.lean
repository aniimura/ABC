/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm52ModelType
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] `Φ^birat` を `𝒟` 上の単系として立てる

原文 (FrdI p.83):
> (iii) There exists a unique subfunctor of groups Φbirat ⊆Φgp such that

★`phiBiratAt P G A`(`Prop44Gp.lean`)は **`𝒞` の対象 `A`** で添字づけられている。
原文の `Φ^birat` は **`𝒟` 上の**部分関手なので、`𝒟` の対象へ移す必要がある。

## ★段取り

1. **底の同型に沿う輸送**(`mem_phiBiratAt_baseIso`)——
   `Definition 1.3, (i)(b)` の pre-step の span で挟み、`mem_phiBiratAt_transport` を 2 回使う。
2. **pull-back 射に沿う輸送**(`mem_phiBiratAt_pullBack`)——
   `𝒪^×` の pull-back(`otimesPull`)を `𝒞^birat` で行い、四角形から `Div^gp` の関係を読む
   (`biratDivGp_of_square`)。
3. `Definition 1.3, (i)(a)` の `baseSurj` で各 `d ∈ Ob(𝒟)` の上の対象を選び、
   1. で**選び方によらない**ことを言う(`phiBiratOn` / `mem_phiBiratOn_iff`)。
4. 2. と `plBk_baseChange` で**部分関手性**(`phiBiratOn_map`)。
5. `MonoidOn` の 3 フィールドを埋める(`phiBiratMonoidOn`)。

## ★★逸脱(記録)

`MonoidOn`(= `Definition 1.1` の「monoid on `𝒟`」)は
**(a) characteristic injectivity** と **(b) FSM-morphism なら同型** を要求する。
★原文は `Proposition 4.4, (iii)` で `Φ^birat` を「subfunctor of groups」としか呼ばず、
(a)(b) を主張していない。しかし `Theorem 5.2` の有理関数の単系 `B` は
「group-like monoid on `𝒟`」でなければならないので、`Proposition 5.3` の
「`𝒞^un-tr` は `(Φ, Φ^birat)` の model Frobenioid」を言うには (a)(b) が要る。

* **(a)** は `Φ^gp` から**無料で**出る(群なので `charMap` の単射性は自動)。
* **(b)** は `𝒟` が **of FSM-type**(§0)であることを**追加の仮定** `hfsmD` に置いた。
  ★★これは原文にない仮定である。★幾何的な `𝒟`(連結有限次被覆の圏など)では
  FSM-morphism は同型なので成り立つが、一般の連結・全射的圏では
  「pull-back で自明化される因子はもとから自明」という主張(= `Pic` の単射性)に
  相当し、一般には期待できない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. 底の同型に沿う輸送 -/

/-- ★★`Φ^birat` は**底の同型に沿って移る**(pre-step の span で挟む)。

★`Definition 1.3, (i)(b)`(`preStepSpan`)が `α = Base(φ)⁻¹ ∘ Base(ψ)` と書き、
`Proposition 4.4, (iii)` の輸送(`mem_phiBiratAt_transport`)を `φ`, `ψ` の両方で使う。 -/
theorem mem_phiBiratAt_baseIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A B : C} (α : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) (hα : IsIso α)
    (y : Gp (Φ.val (P.toElem.obj B).base)) :
    gpMap _ (Φ.map α) y ∈ phiBiratAt P G A ↔ y ∈ phiBiratAt P G B := by
  obtain ⟨X, φ, ψ, hφ, hψ, hspan⟩ := G.core.preStepSpan A B α hα
  have hcφ : IsCoAngular P φ := prop_1_4_i P φ (fun Y _ => hiso Y)
  have hcψ : IsCoAngular P ψ := prop_1_4_i P ψ (fun Y _ => hiso Y)
  haveI := hφ.2
  have hbe : P.Base φ ≫ α = P.Base ψ := by
    rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  rw [← mem_phiBiratAt_transport G φ hcφ hφ (gpMap _ (Φ.map α) y),
    ← mem_phiBiratAt_transport G ψ hcψ hψ y, gpMap_phi_comp, hbe]

/-! ## ★2. pull-back 射に沿う輸送 -/

/-- ★★★**四角形 `φ ≫ δ = ε ≫ φ` から `Div^gp` の関係**(`φ` の次数が 1 のとき)。

★`δ`, `ε` は base-identity かつ linear(`OTri`)なので、
`biratDivGp` の合成則の両辺が `Div^gp(φ)` を 1 つずつ含み、消去できる。 -/
theorem biratDivGp_of_square (G : Frobenioid P) {A B : BiratCat P G} (φ : A ⟶ B)
    (hd : biratDeg φ = 1)
    (δ : OTri (biratPre P G) B) (ε : OTri (biratPre P G) A)
    (hsq : φ ≫ ((δ : End B) : B ⟶ B) = ((ε : End A) : A ⟶ A) ≫ φ) :
    biratDivGp ((ε : End A) : A ⟶ A)
      = gpMap _ (Φ.map (biratBase φ)) (biratDivGp ((δ : End B) : B ⟶ B)) := by
  have hkey := congrArg biratDivGp hsq
  rw [biratDivGp_comp', biratDivGp_comp',
    gpMap_biratBase_of_baseIdentity ε.2.1,
    show ((biratDeg (((δ : End B) : B ⟶ B)) : ℕ+) : ℕ) = 1 from by
      rw [show biratDeg (((δ : End B) : B ⟶ B)) = 1 from δ.2.2]; rfl,
    show ((biratDeg φ : ℕ+) : ℕ) = 1 from by rw [hd]; rfl,
    one_smul, one_smul] at hkey
  have h2 : gpMap _ (Φ.map (biratBase φ)) (biratDivGp ((δ : End B) : B ⟶ B))
      + biratDivGp φ = biratDivGp ((ε : End A) : A ⟶ A) + biratDivGp φ := by
    rw [hkey]; abel
  exact (add_right_cancel h2).symm

/-- ★★★★**`Φ^birat` は pull-back 射に沿って引き戻る**。

★`Definition 1.3, (iv)` の `𝒪^×` の pull-back(`otimesPull`)を
`𝒞^birat`(`birat_isPullBack_iff` で pull-back に移る)で行う。 -/
theorem mem_phiBiratAt_pullBack (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {A B : C} (ν : A ⟶ B) (hν : IsPullBack P ν)
    {z : Gp (Φ.val (P.toElem.obj B).base)} (hz : z ∈ phiBiratAt P G B) :
    gpMap _ (Φ.map (P.Base ν)) z ∈ phiBiratAt P G A := by
  classical
  set Fc : FrobenioidCore (biratPre P G) := birat_frobenioidCore_of_frobNormalized P G hfn
  have hisob : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X :=
    fun X => birat_isOfIsotropicType hiso X
  obtain ⟨hlb, hlin⟩ := G.core.pullBackLB ν hν
  have hpbb : IsPullBack (biratPre P G) ((toBiratCat P G).map ν) :=
    (birat_isPullBack_iff P G Fc ν).mpr ⟨hlb.1, hlin⟩
  obtain ⟨δ, hδ, rfl⟩ := hz
  have hεmem := otimesPull_mem (biratPre P G) Fc hisob _ hpbb ⟨δ, hδ⟩
  have hsq := otimesPull_spec (biratPre P G) Fc hisob _ hpbb ⟨δ, hδ⟩
  have hd : biratDeg ((toBiratCat P G).map ν) = 1 := by rw [biratDeg_toBiratMap]; exact hlin
  have key := biratDivGp_of_square G ((toBiratCat P G).map ν) hd ⟨δ, hδ.1⟩
    (otimesPull (biratPre P G) Fc hisob _ hpbb ⟨δ, hδ⟩) hsq
  rw [biratBase_toBiratMap] at key
  exact ⟨_, hεmem, key⟩

/-! ## ★3. `𝒟` の対象で添字づける -/

/-- ★`Definition 1.3, (i)(a)` が与える「`d` の上の対象」。 -/
noncomputable def biratBaseObj (G : Frobenioid P) (d : D) : C := (G.core.baseSurj d).choose

/-- ★その同型 `Base (biratBaseObj d) ≅ d`。 -/
noncomputable def biratBaseIso (G : Frobenioid P) (d : D) :
    (P.toElem.obj (biratBaseObj G d)).base ≅ d :=
  Classical.choice (G.core.baseSurj d).choose_spec.2

/-- ★★★**`Φ^birat` を `𝒟` の対象で添字づけたもの**。 -/
noncomputable def phiBiratOn (G : Frobenioid P) (d : D) : AddSubgroup (Gp (Φ.val d)) :=
  AddSubgroup.comap (gpMap _ (Φ.map (biratBaseIso G d).hom))
    (phiBiratAt P G ((biratBaseObj G d : C) : BiratCat P G))

/-- ★★★★**選び方によらない** —— `Base A ≅ d` なる**どの** `A` を使ってもよい。 -/
theorem mem_phiBiratOn_iff (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {d : D} {A : C} (e : (P.toElem.obj A).base ≅ d) (y : Gp (Φ.val d)) :
    y ∈ phiBiratOn G d ↔ gpMap _ (Φ.map e.hom) y ∈ phiBiratAt P G A := by
  have hα : IsIso ((biratBaseIso G d).hom ≫ e.inv) := inferInstance
  have h := mem_phiBiratAt_baseIso G hiso ((biratBaseIso G d).hom ≫ e.inv) hα
    (gpMap _ (Φ.map e.hom) y)
  rw [gpMap_phi_comp, Category.assoc, e.inv_hom_id, Category.comp_id] at h
  exact h

/-- ★★**`𝒞` の対象の上では `Φ^birat` はもとの `phiBiratAt` そのもの**。 -/
theorem phiBiratOn_base (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) (A : C) :
    phiBiratOn G (P.toElem.obj A).base = phiBiratAt P G A := by
  ext y
  rw [mem_phiBiratOn_iff G hiso (Iso.refl (P.toElem.obj A).base) y]
  exact iff_of_eq (congrArg (fun z => z ∈ phiBiratAt P G A) (gpMap_phi_id y))

/-! ## ★4. 部分関手性 -/

/-- ★★★★**`Φ^birat` は `𝒟` 上の部分関手** —— `𝒟` の任意の射で引き戻る。

★`plBk_baseChange`(`Definition 1.3, (i)(c)` の本質的全射性)で
`𝒟` の射 `f` を pull-back 射 `α̃` の底に持ち上げ、`mem_phiBiratAt_pullBack` を当てる。 -/
theorem phiBiratOn_map (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') {y : Gp (Φ.val d')} (hy : y ∈ phiBiratOn G d') :
    gpMap _ (Φ.map f) y ∈ phiBiratOn G d := by
  set A' := biratBaseObj G d' with hA'
  set e' := biratBaseIso G d' with he'
  obtain ⟨Yt, αt, k, hαt, hbase⟩ := plBk_baseChange P G.core A' (f ≫ e'.inv)
  rw [mem_phiBiratOn_iff G hiso k, gpMap_phi_comp]
  have hcomp : k.hom ≫ f = P.Base αt ≫ e'.hom := by
    rw [hbase, Category.assoc, Category.assoc, e'.inv_hom_id, Category.comp_id]
  rw [hcomp, ← gpMap_phi_comp]
  exact mem_phiBiratAt_pullBack G hiso hfn αt hαt
    ((mem_phiBiratOn_iff G hiso e' y).mp hy)

/-- ★同型に沿っては**両向き**に移る。 -/
theorem mem_phiBiratOn_iso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') [IsIso f] (y : Gp (Φ.val d')) :
    gpMap _ (Φ.map f) y ∈ phiBiratOn G d ↔ y ∈ phiBiratOn G d' := by
  refine ⟨fun h => ?_, fun h => phiBiratOn_map G hiso hfn f h⟩
  have := phiBiratOn_map G hiso hfn (inv f) h
  rwa [gpMap_phi_inv_right] at this

/-- ★`Φ^birat` の誘導射(`Φ^gp` の誘導射の制限)。 -/
noncomputable def phiBiratMap (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') :
    ↥(phiBiratOn G d') →+ ↥(phiBiratOn G d) :=
  AddMonoidHom.codRestrict ((gpMap _ (Φ.map f)).comp (phiBiratOn G d').subtype) _
    (fun x => phiBiratOn_map G hiso hfn f x.2)

@[simp] theorem phiBiratMap_coe (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    {d d' : D} (f : d ⟶ d') (x : ↥(phiBiratOn G d')) :
    ((phiBiratMap G hiso hfn f x : Gp (Φ.val d))) = gpMap _ (Φ.map f) (x : Gp (Φ.val d')) := rfl

/-- ★★**`Φ^birat` の反変関手** `𝒟ᵒᵖ ⥤ AddCommMonCat`。 -/
noncomputable def phiBiratFunctor (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of ↥(phiBiratOn G X.unop)
  map f := AddCommMonCat.ofHom (phiBiratMap G hiso hfn f.unop)
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    exact Subtype.ext (gpMap_phi_id (x : Gp (Φ.val X.unop)))
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    exact Subtype.ext (gpMap_phi_comp Φ g.unop f.unop (x : Gp (Φ.val X.unop))).symm

/-! ## ★5. `monoid on 𝒟` として立てる -/

/-- ★★★★★**`Φ^birat` は `𝒟` 上の単系(monoid on `𝒟`)**。

★`hint` は `Φ` が整域的であること(divisorial なら従う)。
★★`hfsmD`(`𝒟` が of FSM-type)は**原文にない追加仮定**である
—— ファイル冒頭の「逸脱」を見よ。 -/
noncomputable def phiBiratMonoidOn (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hfsmD : IsOfFSMType D) : MonoidOn.{v, u, w} D where
  functor := phiBiratFunctor G hiso hfn
  charInj {A B} α := by
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val A) (hint A)
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val B) (hint B)
    show IsCharacteristicallyInjective (phiBiratMap G hiso hfn α)
    refine ⟨fun x y h => Subtype.ext (gpMap_injective _ (Φ.map_injective α) ?_),
      charMap_injective_of_addGroup _⟩
    exact congrArg Subtype.val h
  fsmIso {A B} α hα := by
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val A) (hint A)
    letI := isCancelAdd_of_isIntegralMonoid (Φ.val B) (hint B)
    haveI : IsIso α := hfsmD _ _ α hα
    show Function.Bijective (phiBiratMap G hiso hfn α)
    refine ⟨fun x y h => Subtype.ext (gpMap_injective _ (Φ.map_injective α)
      (congrArg Subtype.val h)), fun y => ?_⟩
    refine ⟨⟨gpMap _ (Φ.map (inv α)) (y : Gp (Φ.val B)),
      phiBiratOn_map G hiso hfn (inv α) y.2⟩, Subtype.ext ?_⟩
    exact gpMap_phi_inv_left α (y : Gp (Φ.val B))

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 4.4, (iii)` の「`𝒟` 上の部分関手」の条。 -/
def phiBiratOn_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — Φ^birat は 𝒟 上の部分関手",
    sectionId := "frdi-prop-4-4" }

/-! ## ★6. `Proposition 5.3` の第 2 文 —— `𝒞^un-tr` は `(Φ, Φ^birat)` の model Frobenioid

原文 (FrdI p.103):
> monoid Φbirat (respectively, Q · Φbirat = Φbirat ⊗Z Q = (Φbirat)pf). In particular, if
-/

/-- ★★★★★**unit-trivial 型の Frobenioid の有理関数の単系は `Φ^birat`**。

★`κ` は `Proposition 4.4, (iii)` の全射(unit-trivial 型なので**全単射**
—— `otimesBiratEquivPhiBirat`)であり、`Div_B` は `Φ^birat ⊆ Φ^gp` の包含そのもの。
★`κ` の pull-back 自然性は `biratDivGp_of_square` そのものである。 -/
noncomputable def biratRatFnData (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hut : IsOfUnitTrivialType P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hfsmD : IsOfFSMType D) : RatFnData P G where
  bmon := phiBiratMonoidOn G hiso hfn hint hfsmD
  divB := fun A => (phiBiratOn G A).subtype
  divB_nat := fun f x => rfl
  kappa := fun A => (otimesBiratEquivPhiBirat G hiso hut ((toBiratCat P G).obj A)).trans
    (AddEquiv.addSubgroupCongr (phiBiratOn_base G hiso A)).symm
  kappa_divB := fun A d => rfl
  kappa_pull := by
    intro A B f th d e hpb hth hsq
    refine Subtype.ext ?_
    have hd : biratDeg f = 1 := (birat_pullBackLB P G f hpb).2
    have key := biratDivGp_of_square G f hd ⟨d, d.2.1⟩ ⟨e, e.2.1⟩ hsq
    rw [hth] at key
    exact key
  bneg := fun A (x : ↥(phiBiratOn G A)) => (-x : ↥(phiBiratOn G A))
  bneg_add := fun A (x : ↥(phiBiratOn G A)) => neg_add_cancel x

variable [IsConnected D]

/-- ★`𝒞^un-tr` の有理関数の単系 —— `Φ^birat`。 -/
noncomputable def unTr_ratFnData (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    RatFnData (unTrPre P Fc) (unTr_frobenioid P Fc G) :=
  biratRatFnData (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
    (fun Z => (unTr_isOfModelType Fc G).2 Z) (unTr_unitTrivial P Fc) hint hfsmD

/-- ★★★★★★**[FrdI] Proposition 5.3 の第 2 文(`𝒞^un-tr` の場合)** ——
`𝒞^un-tr` は**因子の単系 `Φ`、有理関数の単系 `Φ^birat`** の model Frobenioid
と圏同値。

★`𝒞^un-tr` が model 型であること(`unTr_isOfModelType`)と
`Theorem 5.2, (iv)`(`modelType_equiv`)の合成である。 -/
noncomputable def unTr_modelFrobenioid (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) (hfsmD : IsOfFSMType D) :
    UnTr P ≌ (unTr_ratFnData Fc G hint hfsmD).model.Obj :=
  modelType_equiv _ (unTr_isotropic P Fc) (unTr_isOfModelType Fc G)

/-- ★locator —— `Proposition 5.3` の第 2 文(`𝒞^un-tr` の場合)。 -/
def unTr_modelFrobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 𝒞^un-tr は (Φ, Φ^birat) の model Frobenioid",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
