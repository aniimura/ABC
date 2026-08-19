/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Sec

/-!
# [FrdI] Theorem 5.1, (iii) —— `𝒞^Fr-tr` は Frobenioid

原文 (FrdI p.97):
> determined by the Frobenius-trivial objects and isometric morphisms is a

原文 (FrdI p.97):
> Frobenioid of isotropic, group-like, base-trivial, and Aut-ample type. In particular, the isomorphism class of a Frobenius-trivial object of C is completely determined by the isomorphism class of its projection to D; all Frobenius-trivial objects of C are Aut-ample.

★`𝒞^Fr-tr` の射はすべて等長なので、単系は `0_𝒟` になる。すると:

| `𝒞^Fr-tr` で | 中身 |
|---|---|
| isometric | **すべての射** |
| pre-step | 次数 1 の base-isomorphism ＝ **同型**(isotropic 型だから) |
| co-angular | **すべての射**(上の帰結) |
| LB-invertible | **すべての射** |
| Frobenius 型 | base-isomorphism |
| pull-back | `𝒞` の pull-back 射(`frTr_isPullBack_iff`) |

★★これで `Definition 1.3` の 21 + 2 条件は、ほとんどが `𝒞` 側から降りてくるか、
「pre-step は同型」から自明になる。★原文の「it follows immediately」の中身である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-- ★sharp な単系では `n • x = 0`(`n ≠ 0`)から `x = 0`。 -/
theorem eq_zero_of_nsmul_eq_zero {M : Type w} [AddCommMonoid M] (hs : IsSharp M) {n : ℕ}
    (hn : n ≠ 0) {x : M} (h : n • x = 0) : x = 0 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_succ_of_ne_zero hn
  rw [succ_nsmul] at h
  exact hs x (IsAddUnit.of_add_eq_zero (k • x) (by rw [add_comm]; exact h))

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★1. 圏 `𝒞^Fr-tr` -/

/-- ★★**`𝒞^Fr-tr`** —— Frobenius-trivial 対象と isometric 射のなす部分圏。 -/
def FrTrCat (P : PreFrobenioid C Φ) : Type u2 := {A : C // IsFrobeniusTrivial P A}

variable {P : PreFrobenioid C Φ}

instance : Category.{v2} (FrTrCat P) where
  Hom A B := {f : A.1 ⟶ B.1 // IsIsometric P f}
  id A := ⟨𝟙 A.1, isIsometric_of_isIso P _⟩
  comp {A B E} f g := ⟨f.1 ≫ g.1, by
    show P.Div (f.1 ≫ g.1) = 0
    rw [P.Div_comp, show P.Div g.1 = 0 from g.2, show P.Div f.1 = 0 from f.2,
      map_zero, smul_zero, add_zero]⟩
  id_comp _ := Subtype.ext (Category.id_comp _)
  comp_id _ := Subtype.ext (Category.comp_id _)
  assoc _ _ _ := Subtype.ext (Category.assoc _ _ _)

@[simp] theorem frTr_id_val (A : FrTrCat P) : (𝟙 A : A ⟶ A).1 = 𝟙 A.1 := rfl

@[simp] theorem frTr_comp_val {A B E : FrTrCat P} (f : A ⟶ B) (g : B ⟶ E) :
    (f ≫ g).1 = f.1 ≫ g.1 := rfl

/-- `𝒞^Fr-tr → 𝔽_{0_𝒟}`。 -/
def frTrToElem (P : PreFrobenioid C Φ) : FrTrCat P ⥤ ElemFrobCat (trivialOn.{v, u, w} D) where
  obj A := ⟨(P.toElem.obj A.1).base⟩
  map {A B} f := ⟨P.Base f.1, PUnit.unit, P.degFr f.1⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ rfl ?_
    · show P.Base (𝟙 A.1) = _; rw [P.Base_id]; rfl
    · show P.degFr (𝟙 A.1) = 1; rw [P.degFr_id]
  map_comp {A B E} f g := by
    refine ElemFrobCat.Hom.ext ?_ rfl ?_
    · show P.Base (f.1 ≫ g.1) = _; rw [P.Base_comp]; rfl
    · show P.degFr (f.1 ≫ g.1) = _; rw [P.degFr_comp]; rfl

/-! ## ★2. pre-Frobenioid 構造 -/

theorem frTr_totEpi : IsTotallyEpimorphic (FrTrCat P) := by
  intro A B f
  refine ⟨fun {Z} g h heq => ?_⟩
  haveI : Epi f.1 := P.totEpiC _ _ f.1
  exact Subtype.ext ((cancel_epi f.1).mp (congrArg Subtype.val heq))

/-- ★底の対象から Frobenius-trivial な持ち上げを選ぶ。 -/
noncomputable def ftLift (G : Frobenioid P) (Y : D) : FrTrCat P :=
  ⟨(G.core.baseSurj Y).choose, (G.core.baseSurj Y).choose_spec.1⟩

noncomputable def ftLiftIso (G : Frobenioid P) (Y : D) :
    (P.toElem.obj (ftLift G Y).1).base ≅ Y :=
  ((G.core.baseSurj Y).choose_spec.2).some

/-- ★`𝒞^Fr-tr` では、底の射が射に持ち上がる。 -/
theorem frTr_hom_of_baseHom (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A B : FrTrCat P} (ψ : (P.toElem.obj A.1).base ⟶ (P.toElem.obj B.1).base) :
    Nonempty (A ⟶ B) := by
  obtain ⟨f, hpb, hb⟩ := exists_pullBack_frobTrivial G hiso A.2 B.2 ψ
  exact ⟨⟨f, (G.core.pullBackLB f hpb).1.2⟩⟩

theorem frTr_isConnected (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) :
    IsConnected (FrTrCat P) := by
  haveI := P.connectedD
  haveI : Nonempty (FrTrCat P) := ⟨ftLift G (Classical.arbitrary D)⟩
  refine zigzag_isConnected ?_
  intro A B
  have hstep : ∀ Z : FrTrCat P, Zigzag Z (ftLift G ((P.toElem.obj Z.1).base)) := by
    intro Z
    obtain ⟨f⟩ := frTr_hom_of_baseHom G hiso (A := Z)
      (B := ftLift G ((P.toElem.obj Z.1).base)) (ftLiftIso G _).inv
    exact Relation.ReflTransGen.single (Or.inl ⟨f⟩)
  have hlift : ∀ {X Y : D}, (X ⟶ Y) → Nonempty (ftLift G X ⟶ ftLift G Y) := by
    intro X Y f
    exact frTr_hom_of_baseHom G hiso ((ftLiftIso G X).hom ≫ f ≫ (ftLiftIso G Y).inv)
  have hzz : ∀ X Y : D, Zigzag X Y → Zigzag (ftLift G X) (ftLift G Y) := by
    intro X Y h
    induction h with
    | refl => exact Relation.ReflTransGen.refl
    | tail _ hstep' ih =>
      refine ih.trans (Relation.ReflTransGen.single ?_)
      cases hstep' with
      | inl h1 => obtain ⟨f⟩ := h1; exact Or.inl (hlift f)
      | inr h1 => obtain ⟨f⟩ := h1; exact Or.inr (hlift f)
  refine (hstep A).trans ((hzz _ _ (isPreconnected_zigzag _ _)).trans (hstep B).symm)

/-- ★★★**`𝒞^Fr-tr → 𝔽_{0_𝒟}`** は pre-Frobenioid 構造。 -/
noncomputable def frTrPre (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) :
    PreFrobenioid (FrTrCat P) (trivialOn.{v, u, w} D) where
  toElem := frTrToElem P
  divisorial := trivialOn_isDivisorialOn D
  totEpiC := frTr_totEpi
  totEpiD := P.totEpiD
  connectedC := frTr_isConnected G hiso
  connectedD := P.connectedD

variable {G : Frobenioid P} {hiso : ∀ Y : C, IsIsotropic P Y}

@[simp] theorem frTrPre_Base {A B : FrTrCat P} (f : A ⟶ B) :
    (frTrPre G hiso).Base f = P.Base f.1 := rfl

@[simp] theorem frTrPre_degFr {A B : FrTrCat P} (f : A ⟶ B) :
    (frTrPre G hiso).degFr f = P.degFr f.1 := rfl

theorem frTrPre_Div {A B : FrTrCat P} (f : A ⟶ B) : (frTrPre G hiso).Div f = 0 :=
  Subsingleton.elim (α := PUnit.{w + 1}) _ _

/-! ## ★3. 辞書 -/

theorem frTr_isIsometric {A B : FrTrCat P} (f : A ⟶ B) :
    IsIsometric (frTrPre G hiso) f := frTrPre_Div f

/-- ★★★**`𝒞^Fr-tr` では pre-step は同型**。 -/
theorem frTr_isIso_of_preStep {A B : FrTrCat P} (f : A ⟶ B)
    (hf : IsPreStep (frTrPre G hiso) f) : IsIso f := by
  haveI : IsIso f.1 := hiso A.1 B.1 f.1 f.2 ⟨hf.1, hf.2⟩
  refine ⟨⟨inv f.1, isIsometric_of_isIso P _⟩, ?_, ?_⟩
  · exact Subtype.ext (IsIso.hom_inv_id f.1)
  · exact Subtype.ext (IsIso.inv_hom_id f.1)

/-- ★★★**`𝒞^Fr-tr` ではすべての射が co-angular**。 -/
theorem frTr_isCoAngular {A B : FrTrCat P} (f : A ⟶ B) :
    IsCoAngular (frTrPre G hiso) f := by
  intro X Y γ β α _ _ _ hβ _
  exact frTr_isIso_of_preStep β hβ

/-- ★★`𝒞` の pull-back 射は `𝒞^Fr-tr` でも pull-back 射。 -/
theorem frTr_isPullBack_of {A B : FrTrCat P} (f : A ⟶ B) (h : IsPullBack P f.1) :
    IsPullBack (frTrPre G hiso) f := by
  intro Z
  constructor
  · intro g g' heq
    have h1 := congrArg (fun t => t.1) heq
    refine Subtype.ext (pullBack_uniq h ?_ ?_)
    · exact congrArg (fun t => (t.1 : Z.1 ⟶ B.1)) (congrArg Prod.fst h1)
    · exact congrArg Prod.snd h1
  · rintro ⟨⟨hh, k⟩, hcomm⟩
    obtain ⟨g₀, hg₀1, hg₀2⟩ := pullBack_lift h hh.1 k hcomm
    have hdiv : P.Div g₀ = 0 := by
      have hz : P.Div hh.1 = 0 := hh.2
      rw [← hg₀1, P.Div_comp, show P.Div f.1 = 0 from f.2, map_zero, zero_add] at hz
      exact eq_zero_of_nsmul_eq_zero (P.divisorial _).2 (by positivity) hz
    exact ⟨⟨g₀, hdiv⟩, Subtype.ext (Prod.ext (Subtype.ext hg₀1) hg₀2)⟩

/-- ★★★**`𝒞^Fr-tr` の pull-back 射は `𝒞` の pull-back 射**(逆も上で示した)。 -/
theorem frTr_isPullBack_iff {A B : FrTrCat P} (f : A ⟶ B) :
    IsPullBack (frTrPre G hiso) f ↔ IsPullBack P f.1 := by
  refine ⟨fun hf => ?_, frTr_isPullBack_of f⟩
  obtain ⟨π, hπpb, hπb⟩ := exists_pullBack_frobTrivial G hiso A.2 B.2 (P.Base f.1)
  set πt : A ⟶ B := ⟨π, (G.core.pullBackLB π hπpb).1.2⟩ with hπt
  have hπt' : IsPullBack (frTrPre G hiso) πt := frTr_isPullBack_of πt hπpb
  obtain ⟨w, hw1, hw2⟩ := pullBack_lift hπt' f (𝟙 _) (by
    show P.Base f.1 = 𝟙 _ ≫ P.Base π
    rw [hπb, Category.id_comp])
  obtain ⟨w', hw'1, hw'2⟩ := pullBack_lift hf πt (𝟙 _) (by
    show P.Base π = 𝟙 _ ≫ P.Base f.1
    rw [hπb, Category.id_comp])
  have hw2' : P.Base w.1 = 𝟙 ((P.toElem.obj A.1).base) := hw2
  have hw'2' : P.Base w'.1 = 𝟙 ((P.toElem.obj A.1).base) := hw'2
  have hww' : w ≫ w' = 𝟙 A := by
    refine pullBack_uniq hf (X := A) ?_ ?_
    · rw [Category.assoc, hw'1, hw1, Category.id_comp]
    · show P.Base (w.1 ≫ w'.1) = P.Base (𝟙 A.1)
      rw [P.Base_comp, hw2', hw'2', Category.id_comp, P.Base_id]
  have hw'w : w' ≫ w = 𝟙 A := by
    refine pullBack_uniq hπt' (X := A) ?_ ?_
    · rw [Category.assoc, hw1, hw'1, Category.id_comp]
    · show P.Base (w'.1 ≫ w.1) = P.Base (𝟙 A.1)
      rw [P.Base_comp, hw2', hw'2', Category.id_comp, P.Base_id]
  haveI : IsIso w.1 := ⟨w'.1, congrArg Subtype.val hww', congrArg Subtype.val hw'w⟩
  have hfeq : f.1 = w.1 ≫ π := congrArg Subtype.val hw1.symm
  rw [hfeq]
  exact IsPullBack.comp P (isPullBack_of_isIso P w.1) hπpb

/-- ★`𝒞^Fr-tr` の対象は Frobenius-trivial。 -/
theorem frTr_isFrobeniusTrivial (A : FrTrCat P) :
    IsFrobeniusTrivial (frTrPre G hiso) A := by
  obtain ⟨ζ, hdeg, hprop⟩ := A.2
  refine ⟨{ toFun := fun n => (⟨((ζ n : End A.1) : A.1 ⟶ A.1), (hprop n).2.1.2⟩ : End A),
            map_one' := ?_, map_mul' := ?_ }, ?_, ?_⟩
  · exact Subtype.ext (congrArg (fun t : End A.1 => (t : A.1 ⟶ A.1)) (map_one ζ))
  · intro m n
    exact Subtype.ext (congrArg (fun t : End A.1 => (t : A.1 ⟶ A.1)) (map_mul ζ m n))
  · intro n; exact hdeg n
  · intro n
    exact ⟨(hprop n).1, ⟨⟨frTr_isCoAngular _, frTrPre_Div _⟩, (hprop n).2.2⟩⟩

/-- ★`𝒞^Fr-tr` はすべての対象が isotropic。 -/
theorem frTr_isIsotropic (A : FrTrCat P) : IsIsotropic (frTrPre G hiso) A := by
  intro Dd φ _ hs
  exact frTr_isIso_of_preStep φ hs

/-- ★`𝒞^Fr-tr` はすべての対象が group-like。 -/
theorem frTr_isGroupLikeObj (A : FrTrCat P) : IsGroupLikeObj (frTrPre G hiso) A := by
  refine (isGroupLike_iff _).mpr (fun a => ?_)
  have ha : a = 0 := Subsingleton.elim (α := PUnit.{w + 1}) _ _
  rw [ha]
  exact isAddUnit_zero

/-- ★★`𝒞^Fr-tr` はすべての対象が base-trivial(`Theorem 5.1, (iii)` の中心)。 -/
theorem frTr_isBaseTrivial (A : FrTrCat P) : IsBaseTrivial (frTrPre G hiso) A := by
  intro Dd hbi
  obtain ⟨e⟩ := hbi
  let e' : (P.toElem.obj A.1).base ≅ (P.toElem.obj Dd.1).base := e
  haveI : IsIso e'.inv := inferInstance
  obtain ⟨θ, hθ⟩ := frobTrivial_iso_of_baseIso G hiso A.2 Dd.2 e'.inv
  refine ⟨⟨⟨θ.hom, isIsometric_of_isIso P _⟩, ⟨θ.inv, isIsometric_of_isIso P _⟩, ?_, ?_⟩⟩
  · exact Subtype.ext θ.hom_inv_id
  · exact Subtype.ext θ.inv_hom_id

/-- ★★`𝒞^Fr-tr` はすべての対象が Aut-ample(`Theorem 5.1, (iii)` の中心)。 -/
theorem frTr_isAutAmple (A : FrTrCat P) : IsAutAmple (frTrPre G hiso) A := by
  intro g hg
  let g' : (P.toElem.obj A.1).base ⟶ (P.toElem.obj A.1).base := g
  haveI : IsIso g' := hg
  obtain ⟨θ, hθ⟩ := frobTrivial_iso_of_baseIso G hiso A.2 A.2 g'
  refine ⟨⟨θ.hom, isIsometric_of_isIso P _⟩, ⟨⟨θ.inv, isIsometric_of_isIso P _⟩,
    Subtype.ext θ.hom_inv_id, Subtype.ext θ.inv_hom_id⟩, hθ⟩

/-! ## ★4. `Definition 1.3` の 21 条件 -/

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** —— `𝒞^Fr-tr` は Frobenioid(の核部分)。 -/
theorem frTr_frobenioidCore : FrobenioidCore (frTrPre G hiso) where
  baseSurj Y := ⟨ftLift G Y, frTr_isFrobeniusTrivial _, ⟨ftLiftIso G Y⟩⟩
  preStepSpan A B α hα := by
    have hb : (frTrPre G hiso).Base (𝟙 A) = 𝟙 _ := (frTrPre G hiso).Base_id A
    have hid : IsPreStep (frTrPre G hiso) (𝟙 A) := by
      refine ⟨(frTrPre G hiso).degFr_id A, ?_⟩
      show IsIso ((frTrPre G hiso).Base (𝟙 A))
      rw [hb]
      infer_instance
    haveI hbi : IsIso ((frTrPre G hiso).Base (𝟙 A)) := hid.2
    have hinv : inv ((frTrPre G hiso).Base (𝟙 A)) = 𝟙 _ := by
      refine IsIso.inv_eq_of_hom_inv_id ?_
      rw [hb, Category.id_comp]
    obtain ⟨θ, hθ⟩ := @frobTrivial_iso_of_baseIso _ _ _ _ _ _ G hiso B.1 A.1 B.2 A.2 α hα
    refine ⟨A, 𝟙 A, ⟨θ.hom, isIsometric_of_isIso P _⟩, hid, ⟨?_, ?_⟩, ?_⟩
    · show P.degFr θ.hom = 1; exact degFr_of_isIso P θ.hom
    · show IsIso (P.Base θ.hom); rw [hθ]; exact hα
    · rw [hinv, Category.id_comp]
      exact hθ.symm
  plBkEquiv A := by
    haveI hfa : (plBkOverFunctor (frTrPre G hiso) A).Faithful := by
      refine ⟨fun {Z W} f g h => ?_⟩
      have hb : (frTrPre G hiso).Base f.left.hom = (frTrPre G hiso).Base g.left.hom :=
        congrArg CommaMorphism.left h
      have hf : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      have hg : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w g)
      exact CommaMorphism.ext (WideSubcategory.hom_ext _
        (pullBack_uniq W.hom.property (hf.trans hg.symm) hb)) (Subsingleton.elim _ _)
    haveI hfu : (plBkOverFunctor (frTrPre G hiso) A).Full := by
      refine ⟨fun {Z W} u => ?_⟩
      obtain ⟨g, hg1, hg2⟩ := pullBack_lift W.hom.property Z.hom.hom u.left (Over.w u).symm
      have hgpb : IsPullBack (frTrPre G hiso) g :=
        isPullBack_of_comp_right (frTrPre G hiso) g W.hom.hom
          (by rw [hg1]; exact Z.hom.property) W.hom.property
      refine ⟨Over.homMk (⟨g, hgpb⟩ : Z.left ⟶ W.left) (WideSubcategory.hom_ext _ hg1), ?_⟩
      exact CommaMorphism.ext hg2 (Subsingleton.elim _ _)
    haveI hes : (plBkOverFunctor (frTrPre G hiso) A).EssSurj := by
      refine ⟨fun X => ?_⟩
      obtain ⟨A₀, π₀, ρ, hpb, hρ⟩ := exists_pullBack_over G A.1 X.hom
      have hA₀ : IsFrobeniusTrivial P A₀ := isFrobeniusTrivial_pullBack G hiso π₀ hpb A.2
      let A₀' : FrTrCat P := ⟨A₀, hA₀⟩
      let π₀' : A₀' ⟶ A := ⟨π₀, (G.core.pullBackLB π₀ hpb).1.2⟩
      have hπpb : IsPullBack (frTrPre G hiso) π₀' := frTr_isPullBack_of π₀' hpb
      refine ⟨Over.mk (show (⟨A₀'⟩ : PlBk (frTrPre G hiso)) ⟶ ⟨A⟩ from ⟨π₀', hπpb⟩), ⟨?_⟩⟩
      exact Over.isoMk ρ hρ
    exact ⟨hfa, hfu, hes⟩
  frobDegSurj A n := by
    obtain ⟨ζ, hdeg, hprop⟩ := frTr_isFrobeniusTrivial (G := G) (hiso := hiso) A
    exact ⟨A, ((ζ n : End A) : A ⟶ A), (hprop n).2, hdeg n⟩
  frobDegUniq A B E φ ψ hφ hψ hd := by
    obtain ⟨β, hβiso, hβ⟩ := G.core.frobDegUniq A.1 B.1 E.1 φ.1 ψ.1
      ⟨⟨prop_1_4_i P _ (fun Y _ => hiso Y), φ.2⟩, hφ.2⟩
      ⟨⟨prop_1_4_i P _ (fun Y _ => hiso Y), ψ.2⟩, hψ.2⟩ hd
    haveI := hβiso
    refine ⟨⟨β, isIsometric_of_isIso P _⟩, ⟨⟨inv β, isIsometric_of_isIso P _⟩,
      Subtype.ext (IsIso.hom_inv_id β), Subtype.ext (IsIso.inv_hom_id β)⟩, Subtype.ext hβ⟩
  coAngularComp _ _ _ _ := frTr_isCoAngular _
  coAngularOfPreStep _ _ _ _ := frTr_isCoAngular _
  otriFwd {A B} φ hc hs α hα := by
    haveI := frTr_isIso_of_preStep φ hs
    refine ⟨(inv φ ≫ (α : A ⟶ A) ≫ φ : End B), ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
    · show (frTrPre G hiso).Base (inv φ ≫ (α : A ⟶ A) ≫ φ) = (frTrPre G hiso).Base (𝟙 B)
      rw [(frTrPre G hiso).Base_comp, (frTrPre G hiso).Base_comp,
        show (frTrPre G hiso).Base ((α : A ⟶ A)) = (frTrPre G hiso).Base (𝟙 A) from hα.1,
        (frTrPre G hiso).Base_id, Category.id_comp, ← (frTrPre G hiso).Base_comp,
        IsIso.inv_hom_id]
    · show (frTrPre G hiso).degFr (inv φ ≫ (α : A ⟶ A) ≫ φ) = 1
      rw [(frTrPre G hiso).degFr_comp, (frTrPre G hiso).degFr_comp,
        show (frTrPre G hiso).degFr ((α : A ⟶ A)) = 1 from hα.2,
        degFr_of_isIso (frTrPre G hiso) (inv φ), hs.1, one_mul, mul_one]
    · rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
    · rintro β ⟨-, hβ⟩
      haveI : Epi φ := (frTrPre G hiso).totEpiC _ _ φ
      refine (cancel_epi φ).mp ?_
      rw [hβ, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  otriBwd {A B} φ hc hs β hβ := by
    haveI := frTr_isIso_of_preStep φ hs
    refine ⟨(φ ≫ (β : B ⟶ B) ≫ inv φ : End A), ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
    · show (frTrPre G hiso).Base (φ ≫ (β : B ⟶ B) ≫ inv φ) = (frTrPre G hiso).Base (𝟙 A)
      rw [(frTrPre G hiso).Base_comp, (frTrPre G hiso).Base_comp,
        show (frTrPre G hiso).Base ((β : B ⟶ B)) = (frTrPre G hiso).Base (𝟙 B) from hβ.1,
        (frTrPre G hiso).Base_id, Category.id_comp, ← (frTrPre G hiso).Base_comp,
        IsIso.hom_inv_id]
    · show (frTrPre G hiso).degFr (φ ≫ (β : B ⟶ B) ≫ inv φ) = 1
      rw [(frTrPre G hiso).degFr_comp, (frTrPre G hiso).degFr_comp,
        show (frTrPre G hiso).degFr ((β : B ⟶ B)) = 1 from hβ.2,
        degFr_of_isIso (frTrPre G hiso) (inv φ), hs.1, one_mul, mul_one]
    · rw [Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
    · rintro α' ⟨-, hα'⟩
      refine (?_ : (φ ≫ ((β : B ⟶ B)) ≫ inv φ) = ((α' : A ⟶ A))).symm
      rw [← Category.assoc, hα', Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  otriBase {A B} φ φ' hcφ hsφ hcφ' hsφ' hbe α hα β hβ heq := by
    have := G.core.otriBase φ.1 φ'.1
      (prop_1_4_i P _ (fun Y _ => hiso Y)) ⟨hsφ.1, hsφ.2⟩
      (prop_1_4_i P _ (fun Y _ => hiso Y)) ⟨hsφ'.1, hsφ'.2⟩ hbe
      (α : A ⟶ A).1 ⟨hα.1, hα.2⟩ (β : B ⟶ B).1 ⟨hβ.1, hβ.2⟩
      (congrArg Subtype.val heq)
    exact Subtype.ext this
  arbFactor {A B} φ := by
    obtain ⟨X₀, Y₀, γ₀, β₀, α₀, hfac, hγ, hβ, hα⟩ := G.core.arbFactor φ.1
    have hdα : P.Div α₀ = 0 := (G.core.pullBackLB α₀ hα).1.2
    have hdegα : P.degFr α₀ = 1 := (G.core.pullBackLB α₀ hα).2
    have hdγ : P.Div γ₀ = 0 := hγ.1.2
    haveI hbγ : IsIso (P.Base γ₀) := hγ.2
    have hdβ : P.Div β₀ = 0 := by
      have h := congrArg P.Div hfac
      rw [show P.Div φ.1 = 0 from φ.2, P.Div_comp, P.Div_comp, hdα, hdγ, map_zero,
        zero_add, smul_zero, add_zero, hdegα] at h
      have h2 : Φ.map (inv (P.Base γ₀)) (Φ.map (P.Base γ₀)
          (((1 : ℕ+) : ℕ) • P.Div β₀)) = Φ.map (inv (P.Base γ₀)) 0 := by rw [← h]
      rw [phi_map_inv_comp, map_zero, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul] at h2
      exact h2
    -- 中間対象が Frobenius-trivial であること
    obtain ⟨ζ, hdegζ, hpropζ⟩ := A.2
    obtain ⟨ρ, hρiso, hρ⟩ := G.core.frobDegUniq A.1 X₀ A.1 γ₀
      (((ζ (P.degFr γ₀) : End A.1) : A.1 ⟶ A.1)) hγ (hpropζ _).2 (hdegζ _).symm
    haveI := hρiso
    have hX : IsFrobeniusTrivial P X₀ := isFrobeniusTrivial_of_iso P G.core (asIso ρ).symm A.2
    haveI hβiso : IsIso β₀ := hiso X₀ Y₀ β₀ hdβ hβ
    have hY : IsFrobeniusTrivial P Y₀ := isFrobeniusTrivial_of_iso P G.core (asIso β₀) hX
    refine ⟨⟨X₀, hX⟩, ⟨Y₀, hY⟩, ⟨γ₀, hdγ⟩, ⟨β₀, hdβ⟩, ⟨α₀, hdα⟩, Subtype.ext hfac,
      ⟨⟨frTr_isCoAngular _, frTrPre_Div _⟩, hγ.2⟩, ⟨hβ.1, hβ.2⟩, frTr_isPullBack_of _ hα⟩
  arbFactorUniq X Y X' Y' γ β α γ' β' α' hfac hγ hβ hα hγ' hβ' hα' := by
    obtain ⟨δ, ε, h1, h2, h3⟩ := G.core.arbFactorUniq X.1 Y.1 X'.1 Y'.1 γ.1 β.1 α.1
      γ'.1 β'.1 α'.1 (congrArg Subtype.val hfac)
      ⟨⟨prop_1_4_i P _ (fun Z _ => hiso Z), γ.2⟩, hγ.2⟩ ⟨hβ.1, hβ.2⟩
      ((frTr_isPullBack_iff α).mp hα)
      ⟨⟨prop_1_4_i P _ (fun Z _ => hiso Z), γ'.2⟩, hγ'.2⟩ ⟨hβ'.1, hβ'.2⟩
      ((frTr_isPullBack_iff α').mp hα')
    refine ⟨⟨⟨δ.hom, isIsometric_of_isIso P _⟩, ⟨δ.inv, isIsometric_of_isIso P _⟩,
        Subtype.ext δ.hom_inv_id, Subtype.ext δ.inv_hom_id⟩,
      ⟨⟨ε.hom, isIsometric_of_isIso P _⟩, ⟨ε.inv, isIsometric_of_isIso P _⟩,
        Subtype.ext ε.hom_inv_id, Subtype.ext ε.inv_hom_id⟩,
      Subtype.ext h1, Subtype.ext h2, Subtype.ext h3⟩
  pullBackLB α hα := ⟨⟨frTr_isCoAngular _, frTrPre_Div _⟩,
    (G.core.pullBackLB α.1 ((frTr_isPullBack_iff α).mp hα)).2⟩
  preStepMono φ hφ := by
    haveI := frTr_isIso_of_preStep φ hφ
    infer_instance
  preStepFactor {A B} φ hφ :=
    ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, frTr_isCoAngular _, hφ,
      frTrPre_Div _, ⟨(frTrPre G hiso).degFr_id B, by
        show IsIso ((frTrPre G hiso).Base (𝟙 B))
        rw [(frTrPre G hiso).Base_id]; infer_instance⟩⟩
  preStepFactorUniq X X' β α β' α' heq _ hβ _ _ _ hβ' _ _ := by
    haveI := frTr_isIso_of_preStep β hβ
    haveI := frTr_isIso_of_preStep β' hβ'
    refine ⟨⟨inv β ≫ β', inv β' ≫ β, ?_, ?_⟩, ?_, ?_⟩
    · rw [Category.assoc, ← Category.assoc β', IsIso.hom_inv_id, Category.id_comp,
        IsIso.inv_hom_id]
    · rw [Category.assoc, ← Category.assoc β, IsIso.hom_inv_id, Category.id_comp,
        IsIso.inv_hom_id]
    · show α' = (inv β' ≫ β) ≫ α
      rw [Category.assoc, heq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    · show β' = β ≫ inv β ≫ β'
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  preStepFactor' {A B} φ hφ :=
    ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm, frTrPre_Div _,
      ⟨(frTrPre G hiso).degFr_id A, by
        show IsIso ((frTrPre G hiso).Base (𝟙 A))
        rw [(frTrPre G hiso).Base_id]; infer_instance⟩,
      frTr_isCoAngular _, hφ⟩
  preStepFactorUniq' X X' β α β' α' heq _ hβ _ _ _ hβ' _ _ := by
    haveI := frTr_isIso_of_preStep β hβ
    haveI := frTr_isIso_of_preStep β' hβ'
    refine ⟨⟨inv β ≫ β', inv β' ≫ β, ?_, ?_⟩, ?_, ?_⟩
    · rw [Category.assoc, ← Category.assoc β', IsIso.hom_inv_id, Category.id_comp,
        IsIso.inv_hom_id]
    · rw [Category.assoc, ← Category.assoc β, IsIso.hom_inv_id, Category.id_comp,
        IsIso.inv_hom_id]
    · show α' = (inv β' ≫ β) ≫ α
      rw [Category.assoc, heq, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
    · show β' = β ≫ inv β ≫ β'
      rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  faithfulUpToUnits {A B} φ ψ hbe _ _ hsφ _ hsψ := by
    haveI := frTr_isIso_of_preStep ψ hsψ
    haveI := frTr_isIso_of_preStep φ hsφ
    refine ⟨(inv ψ ≫ φ : End B), ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
    · show (frTrPre G hiso).Base (inv ψ ≫ φ) = (frTrPre G hiso).Base (𝟙 B)
      rw [(frTrPre G hiso).Base_comp, show (frTrPre G hiso).Base φ = (frTrPre G hiso).Base ψ
        from hbe, ← (frTrPre G hiso).Base_comp, IsIso.inv_hom_id]
    · show (frTrPre G hiso).degFr (inv ψ ≫ φ) = 1
      rw [(frTrPre G hiso).degFr_comp, degFr_of_isIso (frTrPre G hiso) (inv ψ),
        degFr_of_isIso (frTrPre G hiso) φ, one_mul]
    · exact (isUnit_iff_isIso ((inv ψ ≫ φ : End B))).mpr inferInstance
    · rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  isotropicHullExists A :=
    ⟨A, 𝟙 A, frTrPre_Div _,
      ⟨(frTrPre G hiso).degFr_id A, by
        show IsIso ((frTrPre G hiso).Base (𝟙 A))
        rw [(frTrPre G hiso).Base_id]; infer_instance⟩,
      frTr_isIsotropic A,
      fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
        simpa using hβ.symm⟩⟩
  isotropicClosed _ _ := frTr_isIsotropic _

end ABC3.Found.FrdI

/-! ## ★5. `Definition 1.3, (iii), (d)` の 2 つの圏同値 -/

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** —— `𝒞^Fr-tr` は **Frobenioid**。

★`(iii)(d)` の 2 つの圏同値は、`𝒞^Fr-tr` の co-angular pre-step が**すべて同型**で
あることから出る —— 余切片も切片も「可縮な亜群」になり、
`Order(0)`(1 対象 1 射)と同値になる。 -/
theorem frTr_frobenioid : Frobenioid (frTrPre G hiso) where
  core := frTr_frobenioidCore
  coaPreUnderEquiv := by
    letI := coaPreProp_isMultiplicative (frTrPre G hiso)
      (frTr_frobenioidCore (G := G) (hiso := hiso)).coAngularComp
    intro A
    refine ⟨⟨fun {Z W} {f g} _ => ?_⟩, ⟨fun {Z W} _ => ?_⟩, ⟨fun X => ?_⟩⟩
    · haveI : Epi Z.hom.hom := (frTrPre G hiso).totEpiC _ _ _
      have hf : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Under.w f)
      have hg : Z.hom.hom ≫ g.right.hom = W.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Under.w g)
      exact CommaMorphism.ext (Subsingleton.elim _ _)
        (WideSubcategory.hom_ext _ ((cancel_epi Z.hom.hom).mp (hf.trans hg.symm)))
    · haveI : IsIso Z.hom.hom := frTr_isIso_of_preStep _ Z.hom.property.2
      refine ⟨Under.homMk (show Z.right ⟶ W.right from
        ⟨inv Z.hom.hom ≫ W.hom.hom, frTr_isCoAngular _,
          ⟨by
            show (frTrPre G hiso).degFr (inv Z.hom.hom ≫ W.hom.hom) = 1
            rw [(frTrPre G hiso).degFr_comp, degFr_of_isIso (frTrPre G hiso) (inv Z.hom.hom),
              W.hom.property.2.1, one_mul],
           by
            show IsIso ((frTrPre G hiso).Base (inv Z.hom.hom ≫ W.hom.hom))
            rw [(frTrPre G hiso).Base_comp]
            haveI := W.hom.property.2.2
            haveI : IsIso ((frTrPre G hiso).Base (inv Z.hom.hom)) := by
              haveI := Z.hom.property.2.2
              have h : (frTrPre G hiso).Base (inv Z.hom.hom)
                  ≫ (frTrPre G hiso).Base Z.hom.hom = 𝟙 _ := by
                rw [← (frTrPre G hiso).Base_comp, IsIso.inv_hom_id]
                exact (frTrPre G hiso).Base_id _
              have h2 : (frTrPre G hiso).Base Z.hom.hom
                  ≫ (frTrPre G hiso).Base (inv Z.hom.hom) = 𝟙 _ := by
                rw [← (frTrPre G hiso).Base_comp, IsIso.hom_inv_id]
                exact (frTrPre G hiso).Base_id _
              exact ⟨_, h2, h⟩
            infer_instance⟩⟩)
        (WideSubcategory.hom_ext _ (by
          show Z.hom.hom ≫ inv Z.hom.hom ≫ W.hom.hom = W.hom.hom
          rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp])), Subsingleton.elim _ _⟩
    · exact ⟨Under.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp (frTrPre G hiso)))),
        ⟨eqToIso (Subsingleton.elim (α := PUnit.{w + 1}) _ _)⟩⟩
  coaPreOverEquiv := by
    letI := coaPreProp_isMultiplicative (frTrPre G hiso)
      (frTr_frobenioidCore (G := G) (hiso := hiso)).coAngularComp
    intro A
    refine ⟨⟨fun {Z W} {f g} _ => ?_⟩, ⟨fun {Z W} _ => ?_⟩, ⟨fun X => ?_⟩⟩
    · haveI : Mono W.hom.hom := (frTr_frobenioidCore (G := G) (hiso := hiso)).preStepMono
        _ W.hom.property.2
      have hf : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      have hg : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w g)
      exact CommaMorphism.ext
        (WideSubcategory.hom_ext _ ((cancel_mono W.hom.hom).mp (hf.trans hg.symm)))
        (Subsingleton.elim _ _)
    · haveI : IsIso W.hom.hom := frTr_isIso_of_preStep _ W.hom.property.2
      refine ⟨Over.homMk (show Z.left ⟶ W.left from
        ⟨Z.hom.hom ≫ inv W.hom.hom, frTr_isCoAngular _,
          ⟨by
            show (frTrPre G hiso).degFr (Z.hom.hom ≫ inv W.hom.hom) = 1
            rw [(frTrPre G hiso).degFr_comp, degFr_of_isIso (frTrPre G hiso) (inv W.hom.hom),
              Z.hom.property.2.1, mul_one],
           by
            show IsIso ((frTrPre G hiso).Base (Z.hom.hom ≫ inv W.hom.hom))
            rw [(frTrPre G hiso).Base_comp]
            haveI := Z.hom.property.2.2
            haveI : IsIso ((frTrPre G hiso).Base (inv W.hom.hom)) := by
              haveI := W.hom.property.2.2
              have h : (frTrPre G hiso).Base (inv W.hom.hom)
                  ≫ (frTrPre G hiso).Base W.hom.hom = 𝟙 _ := by
                rw [← (frTrPre G hiso).Base_comp, IsIso.inv_hom_id]
                exact (frTrPre G hiso).Base_id _
              have h2 : (frTrPre G hiso).Base W.hom.hom
                  ≫ (frTrPre G hiso).Base (inv W.hom.hom) = 𝟙 _ := by
                rw [← (frTrPre G hiso).Base_comp, IsIso.hom_inv_id]
                exact (frTrPre G hiso).Base_id _
              exact ⟨_, h2, h⟩
            infer_instance⟩⟩)
        (WideSubcategory.hom_ext _ (by
          show (Z.hom.hom ≫ inv W.hom.hom) ≫ W.hom.hom = Z.hom.hom
          rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id])), Subsingleton.elim _ _⟩
    · exact ⟨Over.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp (frTrPre G hiso)))),
        ⟨eqToIso (Subsingleton.elim (α := (OrderCat PUnit.{w + 1})ᵒᵖ) _ _)⟩⟩

/-! ## ★6. `Theorem 5.1, (iii)` の型 -/

/-- ★★★★**[FrdI] Theorem 5.1, (iii)** —— `𝒞^Fr-tr` は
**isotropic・group-like・base-trivial・Aut-ample 型**の Frobenioid。 -/
theorem thm_5_1_iii :
    Frobenioid (frTrPre G hiso) ∧
      (∀ A : FrTrCat P, IsIsotropic (frTrPre G hiso) A) ∧
      (∀ A : FrTrCat P, IsGroupLikeObj (frTrPre G hiso) A) ∧
      (∀ A : FrTrCat P, IsBaseTrivial (frTrPre G hiso) A) ∧
      (∀ A : FrTrCat P, IsAutAmple (frTrPre G hiso) A) :=
  ⟨frTr_frobenioid, frTr_isIsotropic, frTr_isGroupLikeObj, frTr_isBaseTrivial, frTr_isAutAmple⟩

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (iii)`。 -/
def thm_5_1_iii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 97, item := "Theorem 5.1, (iii)",
    sectionId := "frdi-thm-5-1" }
