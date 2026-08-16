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

end ABC3.Found.FrdI
