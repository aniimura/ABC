import ABC3.Found.FrdI.Prop33

/-!
# [FrdI] Definition 3.1, (iv) の `𝒞^un-tr` と Proposition 3.3, (iii)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.57–p.60。

原文 (FrdI p.57):
> for the set of unit-equivalence classes of morphisms A

原文 (FrdI p.60):
> which is full and essentially surjective; moreover, this functor is an equivalence

## ★`Proposition 3.3, (ii)` が効くところ

原文は `Definition 3.1, (iv)` で
「`Proposition 3.3, (ii)` により `≈` は同値関係であり、しかも合成で閉じている」
と**前方参照**する。★**その (ii) を取ったので、ここで商が作れる。**

★★**(ii) は `α₁ ≈ α₂ ⟺ 𝔽_Φ へ同じ射に写る` と言っている**ので、
商は「`𝔽_Φ` への像で同一視する」ことに他ならない。
★**同値関係性も合成で閉じることも、`toElem` が関手であることから直ちに従う。**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★`Hom^un-tr` -/

/-- ★**`Proposition 3.3, (ii)` が与える同値関係** —— `𝔽_Φ` へ同じ射に写ること。 -/
def unTrSetoid (A B : C) : Setoid (A ⟶ B) where
  r α₁ α₂ := P.toElem.map α₁ = P.toElem.map α₂
  iseqv := ⟨fun _ => rfl, Eq.symm, Eq.trans⟩

/-- ★★★**[FrdI] Definition 3.1, (iv)** —— `Hom^un-tr_{𝒞^istr}(A, B)`。 -/
def HomUnTr (A B : C) : Type v2 := Quotient (unTrSetoid P A B)

/-- ★`Hom` から `Hom^un-tr` への自然な全射。 -/
def toHomUnTr {A B : C} (α : A ⟶ B) : HomUnTr P A B := Quotient.mk (unTrSetoid P A B) α

theorem toHomUnTr_eq_iff {A B : C} (α₁ α₂ : A ⟶ B) :
    toHomUnTr P α₁ = toHomUnTr P α₂ ↔ P.toElem.map α₁ = P.toElem.map α₂ :=
  Quotient.eq (r := unTrSetoid P A B)

/-- ★★**`unit-equivalent` との一致**(`Proposition 3.3, (ii)`)。 -/
theorem toHomUnTr_eq_iff_unitEquivalent (Fc : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (α₁ α₂ : A ⟶ B) :
    toHomUnTr P α₁ = toHomUnTr P α₂ ↔ IsUnitEquivalent P α₁ α₂ := by
  rw [toHomUnTr_eq_iff]
  constructor
  · intro h
    refine (prop_3_3_ii P Fc hiso α₁ α₂).mpr ⟨?_, ?_, ?_⟩
    · exact congrArg ElemFrobCat.Hom.deg h
    · exact congrArg ElemFrobCat.Hom.div h
    · exact congrArg ElemFrobCat.Hom.base h
  · exact prop_3_3_ii_toElem P

/-! ## ★`𝒞^un-tr`

★対象は `𝒞^istr` の対象、射は `Hom^un-tr`。
★**合成が well-defined なのは `toElem` が関手だから。**
-/

/-- ★★★**[FrdI] Definition 3.1, (iv)** —— `𝒞^un-tr`。 -/
def UnTr : Type u2 := Istr P

instance : Category.{v2} (UnTr P) where
  Hom A B := HomUnTr P (A : Istr P).obj (B : Istr P).obj
  id A := toHomUnTr P (𝟙 (A : Istr P).obj)
  comp {A B E} f g :=
    Quotient.liftOn₂ f g (fun α β => toHomUnTr P (α ≫ β))
      (fun α₁ β₁ α₂ β₂ ha hb => by
        refine (toHomUnTr_eq_iff P _ _).mpr ?_
        rw [P.toElem.map_comp, P.toElem.map_comp,
          show P.toElem.map α₁ = P.toElem.map α₂ from ha,
          show P.toElem.map β₁ = P.toElem.map β₂ from hb])
  id_comp {A B} f := by
    refine Quotient.inductionOn f (fun α => ?_)
    show toHomUnTr P (𝟙 _ ≫ α) = toHomUnTr P α
    rw [Category.id_comp]
  comp_id {A B} f := by
    refine Quotient.inductionOn f (fun α => ?_)
    show toHomUnTr P (α ≫ 𝟙 _) = toHomUnTr P α
    rw [Category.comp_id]
  assoc {A B E G} f g h := by
    refine Quotient.inductionOn₃ f g h (fun α β γ => ?_)
    show toHomUnTr P ((α ≫ β) ≫ γ) = toHomUnTr P (α ≫ β ≫ γ)
    rw [Category.assoc]

/-- ★★**`𝒞^istr → 𝒞^un-tr`** —— 対象は同じ、射は商へ落とす。 -/
def istrToUnTr : Istr P ⥤ UnTr P where
  obj A := A
  map {_ _} f := toHomUnTr P f.hom
  map_id _ := rfl
  map_comp _ _ := rfl

/-! ## ★`Proposition 3.3, (iii)`

原文 (FrdI p.60):
> of categories if and only if Cistr is of unit-trivial type.
-/

/-- ★★**`𝒞^istr → 𝒞^un-tr` は full** —— 商への全射なので構成から。 -/
instance : (istrToUnTr P).Full where
  map_surjective {A B} f := Quotient.inductionOn f (fun α => ⟨ObjectProperty.homMk α, rfl⟩)

/-- ★★**`𝒞^istr → 𝒞^un-tr` は本質的全射** —— 対象を変えないので構成から。 -/
instance : (istrToUnTr P).EssSurj where
  mem_essImage A := ⟨show Istr P from A, ⟨Iso.refl _⟩⟩

/-- ★★★**[FrdI] Proposition 3.3, (iii)** —— `𝒞^istr → 𝒞^un-tr` が**忠実**
(したがって圏同値)であるのは、ちょうど `𝒞^istr` が **unit-trivial 型**のとき。

★`⟸` は `Proposition 3.3, (ii)` から: `𝔽_Φ` で一致すれば unit-equivalent、
`𝒪^× = ⊥` なら `δ = 1` なので 2 射は等しい。
★`⟹` は `𝟙` と単元 `δ` がつねに unit-equivalent なことから。 -/
theorem prop_3_3_iii (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X) :
    (istrToUnTr P).Faithful ↔ ∀ A : Istr P, IsUnitTrivial P A.obj := by
  constructor
  · intro hf A
    refine le_antisymm (fun δ hδ => ?_) bot_le
    haveI := hf
    -- `𝟙` と `δ` は unit-equivalent
    have hue : IsUnitEquivalent P (𝟙 A.obj) ((δ : A.obj ⟶ A.obj)) :=
      ⟨A.obj, 𝟙 _, 𝟙 _, δ, hδ, by simp, by simp⟩
    have h : toHomUnTr P (𝟙 A.obj) = toHomUnTr P ((δ : A.obj ⟶ A.obj)) :=
      (toHomUnTr_eq_iff_unitEquivalent P Fc hiso _ _).mpr hue
    have h2 : (ObjectProperty.homMk (𝟙 A.obj) : A ⟶ A)
        = ObjectProperty.homMk ((δ : A.obj ⟶ A.obj)) := (istrToUnTr P).map_injective h
    exact Submonoid.mem_bot.mpr (congrArg (fun z : A ⟶ A => z.hom) h2).symm
  · intro hut
    refine ⟨fun {A B} {f g} h => ?_⟩
    have h' : toHomUnTr P f.hom = toHomUnTr P g.hom := h
    obtain ⟨Cc, γ, β, δ, hδ, h₁, h₂⟩ :=
      (toHomUnTr_eq_iff_unitEquivalent P Fc hiso f.hom g.hom).mp h'
    have hd1 : (δ : Cc ⟶ Cc) = 𝟙 Cc := by
      have hb : δ ∈ (⊥ : Submonoid (End Cc)) := by rw [← hut ⟨Cc, hiso Cc⟩]; exact hδ
      simpa using hb
    refine InducedCategory.hom_ext ?_
    rw [h₁, h₂, hd1]
    simp

/-! ## ★`Proposition 3.3, (iv)` の前半 —— `𝒞^un-tr → 𝔽_Φ` は忠実

原文 (FrdI p.60):
> which is faithful and essentially surjective; moreover, this functor determines

★**忠実性は構成そのもの** —— `𝒞^un-tr` は `𝔽_Φ` への像で同一視して作ったから。
★**`𝒞^un-tr` に Frobenioid の構造を入れる**部分と、射の類型が
`𝒞^istr` から来ること(原文の後半)は別途取る。
-/

/-- ★★**`𝒞^un-tr → 𝔽_Φ`** —— `𝒞^istr → 𝔽_Φ` が商を経由する。 -/
def unTrToElem : UnTr P ⥤ ElemFrobCat Φ where
  obj A := P.toElem.obj (show Istr P from A).obj
  map {_ _} f := Quotient.liftOn f (fun α => P.toElem.map α) (fun _ _ h => h)
  map_id A := P.toElem.map_id (show Istr P from A).obj
  map_comp {_ _ _} f g := by
    refine Quotient.inductionOn₂ f g (fun α β => ?_)
    exact P.toElem.map_comp α β

/-- ★★★**[FrdI] Proposition 3.3, (iv) の前半** —— `𝒞^un-tr → 𝔽_Φ` は**忠実**。

★**構成そのもの**である —— `𝒞^un-tr` の射は `𝔽_Φ` への像で同一視した類だから。 -/
instance : (unTrToElem P).Faithful where
  map_injective {A B} {f g} := by
    refine Quotient.inductionOn₂ f g (fun α β h => ?_)
    exact (toHomUnTr_eq_iff P α β).mpr h

/-- ★**`𝒞^istr → 𝒞^un-tr → 𝔽_Φ` は `𝒞^istr → 𝔽_Φ` に等しい**
(原文の「The functor `𝒞^istr → 𝔽_Φ` factors naturally through `𝒞^un-tr`」)。 -/
theorem istrToUnTr_comp_unTrToElem :
    istrToUnTr P ⋙ unTrToElem P = (isotropicProp P).ι ⋙ P.toElem := rfl

include P in
/-- ★★★**[FrdI] Proposition 3.3, (iv)** —— `𝒞^un-tr → 𝔽_Φ` は**本質的全射**。

★`Definition 1.3, (i), (a)`(`baseSurj`)で底を実現し、
`Definition 1.3, (vii), (a)`(isotropic hull)で `𝒞^istr` へ移す
——`Proposition 2.2, (i)` と同じ 2 手である。 -/
theorem unTrToElem_essSurj (Fc : FrobenioidCore P) : (unTrToElem P).EssSurj := by
  refine ⟨fun X => ?_⟩
  obtain ⟨A₀, -, ⟨e⟩⟩ := Fc.baseSurj X.base
  obtain ⟨B, φ, hφ⟩ := Fc.isotropicHullExists A₀
  haveI hbi : IsIso (P.Base φ) := hφ.2.1.2
  have e' : (P.toElem.obj B).base ≅ X.base := (@asIso _ _ _ _ (P.Base φ) hbi).symm ≪≫ e
  obtain ⟨k, hk⟩ : ∃ k : P.toElem.obj B ⟶ X, k = ⟨e'.hom, 0, 1⟩ := ⟨_, rfl⟩
  have hkiso : IsIso k := by
    refine (ElemFrobCat.isIso_iff k).mpr ⟨?_, ?_, ?_⟩
    · rw [hk]; exact inferInstanceAs (IsIso e'.hom)
    · rw [hk]; exact isAddUnit_zero
    · rw [hk]
  exact ⟨(show UnTr P from (⟨B, hφ.2.2.1⟩ : Istr P)), ⟨@asIso _ _ _ _ k hkiso⟩⟩

/-! ## ★`Proposition 3.3, (iv)` の中盤 —— `𝒞^un-tr` の **pre-Frobenioid 構造**

原文 (FrdI p.60):
> which is faithful and essentially surjective; moreover, this functor determines

★★原文は「determines a natural structure of Frobenioid」と 1 行で言うが、
**まず `PreFrobenioid` の 6 フィールドを埋める**必要がある。★4 つは `P` のものが
そのまま通り、要るのは **`totEpiC` と `connectedC` の 2 つ**である。

★★**`totEpiC` は在庫で済んだ** —— `isTotallyEpimorphic_elemFrobCat`
(`𝒟` が totally epimorphic ＋ `Φ(A)` が integral ⟹ `𝔽_Φ` が totally epimorphic)。
★★★**そこで効いているのは `Definition 1.1, (ii), (a)` の
「`Φ(A) → Φ(B)` は characteristically injective」**である。
★私は一度「`Φ.map` の単射性は無いので `𝒞^un-tr` の全射性は原文の隠れ段だ」と
測定しかけたが、**単射性は `MonoidOn` の構造フィールドとして最初から入っていた**
(`MonoidOn.charInj` → `MonoidOn.map_injective`)。
★★**「原文の隠れ段」と書く前に在庫を引く**、の 3 例目である。
-/

include P in
/-- ★★**`𝒞^un-tr` は totally epimorphic**。

★`𝒞^un-tr → 𝔽_Φ` は**忠実**なので、`𝔽_Φ` 側の epi 性がそのまま降りてくる。 -/
theorem unTr_totEpi : IsTotallyEpimorphic (UnTr P) := by
  intro A B f
  refine ⟨fun {Z} g h hgh => ?_⟩
  haveI : Epi ((unTrToElem P).map f) :=
    isTotallyEpimorphic_elemFrobCat P.totEpiD (fun X => (P.divisorial X).1.1) _ _ _
  refine (unTrToElem P).map_injective ?_
  refine (cancel_epi ((unTrToElem P).map f)).mp ?_
  rw [← (unTrToElem P).map_comp, ← (unTrToElem P).map_comp, hgh]

include P in
/-- ★★**`𝒞^un-tr` は connected** —— `istrToUnTr` は**対象について恒等**なので、
`𝒞^istr` の zigzag をそのまま送るだけ。 -/
theorem isConnected_unTr (Fc : FrobenioidCore P) : IsConnected (UnTr P) := by
  haveI := isConnected_istr P Fc
  haveI : Nonempty (UnTr P) := (inferInstance : Nonempty (Istr P))
  refine zigzag_isConnected ?_
  intro X Y
  exact zigzag_obj_of_zigzag (istrToUnTr P)
    (@isPreconnected_zigzag (Istr P) _ _ (X : Istr P) (Y : Istr P))

/-- ★★★**[FrdI] Proposition 3.3, (iv)** —— `𝒞^un-tr` の **pre-Frobenioid 構造**。

★6 フィールドのうち **4 つは `P` のものがそのまま通る**
(`divisorial` は `Φ` のものだから、`totEpiD`・`connectedD` は `𝒟` が変わらないから)。 -/
noncomputable def unTrPre (Fc : FrobenioidCore P) : PreFrobenioid (UnTr P) Φ where
  toElem := unTrToElem P
  divisorial := P.divisorial
  totEpiC := unTr_totEpi P
  totEpiD := P.totEpiD
  connectedC := isConnected_unTr P Fc
  connectedD := P.connectedD

/-! ### ★射の 3 つの不変量は `𝒞^istr` のものと一致する

原文 (FrdI p.60):
> phism; morphism of a given Frobenius degree) if and only if it arises from such

★★これは **`rfl`** である —— `unTrToElem` は代表元に `P.toElem` を当てるだけだから。
★原文の「if and only if it arises from such an arrow of `𝒞^istr`」の
**「only if」の側が構成から自明になる**のはこのためである。 -/

theorem unTrPre_degFr (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    (unTrPre P Fc).degFr (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      = P.degFr α := rfl

theorem unTrPre_Base (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    (unTrPre P Fc).Base (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      = P.Base α := rfl

theorem unTrPre_Div (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    (unTrPre P Fc).Div (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      = P.Div α := rfl

/-- ★★**`𝒞^un-tr` は unit-trivial 型** —— `𝒪^×` は `𝔽_Φ` で恒等に写る元の集合であり、
`𝒞^un-tr` ではそれらがすべて `𝟙` に潰れている。

★これが原文の「with respect to which `𝒞^un-tr` is of ... unit-trivial type」の中身。 -/
theorem unTr_unitTrivial (Fc : FrobenioidCore P) (A : UnTr P) :
    IsUnitTrivial (unTrPre P Fc) A := by
  refine le_antisymm (fun δ hδ => ?_) bot_le
  refine Submonoid.mem_bot.mpr ?_
  obtain ⟨⟨hbase, hdeg⟩, hunit⟩ := hδ
  haveI : IsIso ((δ : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso _).mp hunit
  haveI : IsIso ((unTrToElem P).map (δ : A ⟶ A)) := inferInstance
  refine (unTrToElem P).map_injective ?_
  have hid : (unTrToElem P).map ((1 : End A) : A ⟶ A) = 𝟙 ((unTrToElem P).obj A) :=
    (unTrToElem P).map_id A
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · exact hbase
  · rw [show ((unTrToElem P).map ((1 : End A) : A ⟶ A)).div = 0 from
      congrArg ElemFrobCat.Hom.div hid]
    exact ElemFrobCat.div_eq_zero_of_isIso (fun X => (P.divisorial X).2) _
  · rw [show ((unTrToElem P).map ((1 : End A) : A ⟶ A)).deg = 1 from
      congrArg ElemFrobCat.Hom.deg hid]
    exact hdeg

/-! ### ★`Proposition 3.3, (iv)` の最終文 —— 射の類型は `𝒞^istr` から来る

原文 (FrdI p.60):
> phism; morphism of a given Frobenius degree) if and only if it arises from such

★原文は **9 クラス**を挙げる: Frobenius 型 / pre-step / base-isomorphism / 同型 /
pull-back / isometry / co-angular / LB-invertible / 与えられた Frobenius 次数。

★★**仕分け(測定)**: このうち **4 クラスは `Base`・`Div`・`degFr` だけで決まる**ので
`Iff.rfl` である(`unTrToElem` は代表元に `P.toElem` を当てるだけだから)。
★**同型は `Iff.rfl` ではないが取れる** —— `𝒞^istr` が isotropic 型であることが効く。
★残る 4 クラス(pull-back / co-angular / LB-invertible / Frobenius 型)は
**圏論的な条件**なので、商を跨ぐには「unit-equivalent な 2 射で条件が一致すること」が要る。
-/

theorem unTr_degFr_iff (Fc : FrobenioidCore P) (n : ℕ+) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    (unTrPre P Fc).degFr (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B)) = n
      ↔ P.degFr α = n := Iff.rfl

theorem unTr_isBaseIsomorphism_iff (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    IsBaseIsomorphism (unTrPre P Fc)
        (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsBaseIsomorphism P α := Iff.rfl

theorem unTr_isIsometric_iff (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    IsIsometric (unTrPre P Fc)
        (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsIsometric P α := Iff.rfl

theorem unTr_isPreStep_iff (Fc : FrobenioidCore P) {A B : Istr P} (α : A.obj ⟶ B.obj) :
    IsPreStep (unTrPre P Fc) (toHomUnTr P α : (show UnTr P from A) ⟶ (show UnTr P from B))
      ↔ IsPreStep P α := Iff.rfl

/-- ★★★**同型も `𝒞^istr` から来る**。

★`⟸` は関手性。★★**`⟹` が中身のある向き**である ——
`𝒞^un-tr` で同型なら `𝔽_Φ` で同型、したがって `Div = 0`・`degFr = 1`・`Base` が同型、
すなわち `α` は **isometric な pre-step**。★`𝒞^istr` は isotropic 型なので
`Definition 1.2, (iv)` により `α` は同型である。

★★**商が同型を「新たに作らない」ことの証明**であり、`Proposition 3.3, (iii)` の
「圏同値 ⟺ unit-trivial 型」とは別の主張である(あちらは忠実性の話)。 -/
theorem unTr_isIso_iff (Fc : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : Istr P} (α : A.obj ⟶ B.obj) :
    @IsIso (UnTr P) _ A B (toHomUnTr P α) ↔ IsIso α := by
  constructor
  · intro h
    haveI := h
    have hE : @IsIso (ElemFrobCat Φ) _ _ _
        ((unTrToElem P).map (X := A) (Y := B) (toHomUnTr P α)) := inferInstance
    have hE' : IsIso (P.toElem.map α) := hE
    obtain ⟨hb, hu, hd⟩ := (ElemFrobCat.isIso_iff (P.toElem.map α)).mp hE'
    have hdiv : IsIsometric P α := by
      have hau : IsAddUnit (P.Div α) := hu
      exact (P.divisorial (P.toElem.obj A.obj).base).2 _ hau
    exact hiso A.obj B.obj α hdiv ⟨hd, hb⟩
  · intro h
    haveI := h
    haveI hmap : IsIso ((isotropicProp P).ι.map (ObjectProperty.homMk α : A ⟶ B)) := h
    haveI : IsIso (ObjectProperty.homMk α : A ⟶ B) :=
      (ObjectProperty.fullyFaithfulι (isotropicProp P)).isIso_of_isIso_map _
    let e := (istrToUnTr P).mapIso (asIso (ObjectProperty.homMk α : A ⟶ B))
    exact ⟨⟨e.inv, e.hom_inv_id, e.inv_hom_id⟩⟩

end ABC3.Found.FrdI
