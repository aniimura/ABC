/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Untr
import ABC3.Found.FrdI.Def28

/-!
# [FrdI] Corollary 5.7, (ii) —— unit-profinite 型は圏同値で移る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.108。

原文 (FrdI p.108):
> (ii) C1 is of unit-profinite type if and only if C2 is.

## ★★何が要るか

原文の証明は「定義を辿れば、`Ψ` が isotropic 対象・prime-Frobenius 射・pull-back 射・
birationalization・射影関手 `Cᵢ → Dᵢ`(したがって単元 `O^×(−)`)を保つことを示せば足りる」
と書き、それを `Theorem 3.4, (i), (iii)`; `Corollary 4.10`; `Corollary 4.11, (ii)` に帰す。

★**(ii) に効くのは「単元 `𝒪^×(−)` の保存」だけ**である。そしてそれは
`Cor411Untr.lean` の `otimes_map_of_baseSquare` が**底の 1-可換図式だけ**から与える。

## ★段取り

| 段 | 中身 |
|---|---|
| §1 | 底の四角形の**逆向き** —— `Base(Ψ f) = 𝟙` から `Base f = 𝟙`(`ΨBase` が忠実) |
| §2 | `𝒪^×(A) ≃* 𝒪^×(Ψ A)`(群同型) |
| §3 | 対象の同型に沿った `𝒪^×` の移送 |
| §4 | 位相を移して `IsOfUnitProfiniteType` を移す |

★★位相は `MulEquiv` に沿って `coinduced` で移す。全単射に沿った `coinduced` は
同相を与えるので、位相群・コンパクト・完全不連結がそのまま移る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section Cor57Unit

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-! ## ★1. 底の四角形の逆向き -/

/-- ★★**`ΨBase` が忠実なら base-identity は逆向きにも移る**。

★四角形の自然性 `ΨBase(Base φ) ≫ sq = sq ≫ Base(Ψ φ)` で、
右辺が `sq` になるので `ΨBase(Base φ) = 𝟙`、忠実性から `Base φ = 𝟙`。 -/
theorem isBaseIdentity_of_map_baseSquare (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂) [ΨBase.Faithful]
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj)
    {A : C₁} (φ : End A) (h : IsBaseIdentity P₂ (Ψ.map (φ : A ⟶ A))) :
    IsBaseIdentity P₁ φ := by
  have hnat : ΨBase.map (P₁.Base (φ : A ⟶ A)) ≫ sq.hom.app A
      = sq.hom.app A ≫ P₂.Base (Ψ.map (φ : A ⟶ A)) := sq.hom.naturality (φ : A ⟶ A)
  have hb : P₂.Base (Ψ.map (φ : A ⟶ A)) = P₂.Base (𝟙 (Ψ.obj A)) := h
  have hR : sq.hom.app A ≫ P₂.Base (Ψ.map (φ : A ⟶ A)) = sq.hom.app A := by
    rw [hb]
    exact (congrArg (fun t => sq.hom.app A ≫ t) (P₂.Base_id (Ψ.obj A))).trans
      (Category.comp_id (sq.hom.app A))
  have hid : ΨBase.map (P₁.Base (φ : A ⟶ A)) = 𝟙 _ :=
    (cancel_mono (sq.hom.app A)).mp
      (hnat.trans (hR.trans (Category.id_comp (sq.hom.app A)).symm))
  show P₁.Base (φ : A ⟶ A) = P₁.Base (𝟙 A)
  refine ΨBase.map_injective ?_
  rw [P₁.Base_id, ΨBase.map_id]
  exact hid

/-! ## ★2. `𝒪^×(A) ≃* 𝒪^×(Ψ A)` -/

/-- ★`Ψ` が誘導する `𝒪^×` の単系準同型。 -/
noncomputable def otimesMapHom (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj) (X : C₁) :
    OTimes P₁ X →* OTimes P₂ (Ψ.obj X) where
  toFun δ := ⟨Functor.mapEnd X Ψ (δ : End X),
    otimes_map_of_baseSquare Ψ ΨBase sq X (δ : End X) δ.2⟩
  map_one' := Subtype.ext (map_one (Functor.mapEnd X Ψ))
  map_mul' _ _ := Subtype.ext (map_mul (Functor.mapEnd X Ψ) _ _)

@[simp] theorem otimesMapHom_coe (Ψ : C₁ ⥤ C₂) (ΨBase : D₁ ⥤ D₂)
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj) (X : C₁) (δ : OTimes P₁ X) :
    ((otimesMapHom Ψ ΨBase sq X δ : End (Ψ.obj X)) : Ψ.obj X ⟶ Ψ.obj X)
      = Ψ.map ((δ : End X) : X ⟶ X) := rfl

/-- ★★★★**`Ψ` は `𝒪^×` を群同型に移す**。

★単射性は `Ψ` の忠実性、全射性は `Ψ` の充満性と §1(底の四角形の逆向き)。 -/
noncomputable def otimesEquivOfBaseSquare (Ψ : C₁ ⥤ C₂) [Ψ.Full] [Ψ.Faithful]
    (ΨBase : D₁ ⥤ D₂) [ΨBase.Faithful]
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj) (X : C₁) :
    OTimes P₁ X ≃* OTimes P₂ (Ψ.obj X) :=
  MulEquiv.ofBijective (otimesMapHom Ψ ΨBase sq X)
    ⟨by
      intro δ ε h
      refine Subtype.ext ?_
      exact Ψ.map_injective (congrArg
        (fun t : OTimes P₂ (Ψ.obj X) => ((t : End (Ψ.obj X)) : Ψ.obj X ⟶ Ψ.obj X)) h),
     by
      rintro ⟨ε, hε⟩
      refine ⟨⟨(Ψ.preimage ((ε : Ψ.obj X ⟶ Ψ.obj X)) : End X), ?_⟩, ?_⟩
      · refine (mem_otimes_iff P₁ _).mpr ⟨?_, ?_⟩
        · refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
          haveI : IsIso ((ε : Ψ.obj X ⟶ Ψ.obj X)) :=
            (CategoryTheory.isUnit_iff_isIso ε).mp ((mem_otimes_iff P₂ ε).mp hε).1
          haveI : IsIso (Ψ.map (Ψ.preimage ((ε : Ψ.obj X ⟶ Ψ.obj X)))) := by
            rw [Ψ.map_preimage]
            infer_instance
          show IsIso (Ψ.preimage ((ε : Ψ.obj X ⟶ Ψ.obj X)))
          exact isIso_of_reflects_iso _ Ψ
        · refine isBaseIdentity_of_map_baseSquare Ψ ΨBase sq _ ?_
          rw [Ψ.map_preimage]
          exact ((mem_otimes_iff P₂ ε).mp hε).2
      · refine Subtype.ext ?_
        show Ψ.map (Ψ.preimage ((ε : Ψ.obj X ⟶ Ψ.obj X))) = _
        exact Ψ.map_preimage _⟩

/-! ## ★3. 対象の同型に沿った移送 -/

/-- ★対象の同型に沿った `𝒪^×` の移送(共役)。 -/
noncomputable def otimesCongr {A B : C₂} (e : A ≅ B) :
    OTimes P₂ A ≃* OTimes P₂ B where
  toFun δ := ⟨(e.inv ≫ ((δ : End A) : A ⟶ A) ≫ e.hom : B ⟶ B), by
    refine (mem_otimes_iff P₂ _).mpr ⟨?_, ?_⟩
    · haveI : IsIso ((δ : End A) : A ⟶ A) :=
        (CategoryTheory.isUnit_iff_isIso (δ : End A)).mp δ.2.2
      exact (CategoryTheory.isUnit_iff_isIso
        (e.inv ≫ ((δ : End A) : A ⟶ A) ≫ e.hom : End B)).mpr inferInstance
    · have hδ : P₂.Base ((δ : End A) : A ⟶ A) = 𝟙 _ := by
        rw [show P₂.Base ((δ : End A) : A ⟶ A) = P₂.Base (𝟙 A) from δ.2.1.1, P₂.Base_id]
      show P₂.Base (e.inv ≫ ((δ : End A) : A ⟶ A) ≫ e.hom) = P₂.Base (𝟙 B)
      rw [P₂.Base_comp, P₂.Base_comp, hδ, Category.id_comp, ← P₂.Base_comp, e.inv_hom_id,
        P₂.Base_id]⟩
  invFun δ := ⟨(e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv : A ⟶ A), by
    refine (mem_otimes_iff P₂ _).mpr ⟨?_, ?_⟩
    · haveI : IsIso ((δ : End B) : B ⟶ B) :=
        (CategoryTheory.isUnit_iff_isIso (δ : End B)).mp δ.2.2
      exact (CategoryTheory.isUnit_iff_isIso
        (e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv : End A)).mpr inferInstance
    · have hδ : P₂.Base ((δ : End B) : B ⟶ B) = 𝟙 _ := by
        rw [show P₂.Base ((δ : End B) : B ⟶ B) = P₂.Base (𝟙 B) from δ.2.1.1, P₂.Base_id]
      show P₂.Base (e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv) = P₂.Base (𝟙 A)
      rw [P₂.Base_comp, P₂.Base_comp, hδ, Category.id_comp, ← P₂.Base_comp, e.hom_inv_id,
        P₂.Base_id]⟩
  left_inv δ := by
    refine Subtype.ext ?_
    show e.hom ≫ (e.inv ≫ ((δ : End A) : A ⟶ A) ≫ e.hom) ≫ e.inv = _
    simp
  right_inv δ := by
    refine Subtype.ext ?_
    show e.inv ≫ (e.hom ≫ ((δ : End B) : B ⟶ B) ≫ e.inv) ≫ e.hom = _
    simp
  map_mul' δ ε := by
    refine Subtype.ext ?_
    show e.inv ≫ (((ε : End A) : A ⟶ A) ≫ ((δ : End A) : A ⟶ A)) ≫ e.hom
      = (e.inv ≫ ((ε : End A) : A ⟶ A) ≫ e.hom) ≫ (e.inv ≫ ((δ : End A) : A ⟶ A) ≫ e.hom)
    simp only [Category.assoc, Iso.hom_inv_id_assoc]

/-! ## ★4. 位相の移送 -/

/-- ★★**群同型に沿って位相を移す** —— `coinduced` は全単射に沿って同相を与えるので、
位相群・コンパクト・完全不連結・位相的有限生成がそのまま移る。 -/
theorem exists_top_of_mulEquiv {G₁ : Type*} [CommGroup G₁] [TopologicalSpace G₁]
    [IsTopologicalGroup G₁] [CompactSpace G₁] [TotallyDisconnectedSpace G₁]
    {G₂ : Type*} [Group G₂] (E : G₁ ≃* G₂) (S : Finset G₁)
    (hS : closure ((Subgroup.closure (S : Set G₁) : Subgroup G₁) : Set G₁) = Set.univ) :
    ∃ t : TopologicalSpace G₂,
      (∀ a b : G₂, a * b = b * a) ∧
      (letI := t; IsTopologicalGroup G₂) ∧
      (letI := t; CompactSpace G₂) ∧
      (letI := t; TotallyDisconnectedSpace G₂) ∧
      ∃ T : Finset G₂,
        (letI := t;
          closure ((Subgroup.closure (T : Set G₂) : Subgroup G₂) : Set G₂) = Set.univ) := by
  classical
  letI t2 : TopologicalSpace G₂ := TopologicalSpace.coinduced (E : G₁ → G₂) inferInstance
  have hcont : Continuous (E : G₁ → G₂) := continuous_coinduced_rng
  have hopen : IsOpenMap (E : G₁ → G₂) := by
    intro U hU
    show IsOpen ((E : G₁ → G₂) '' U)
    rw [isOpen_coinduced, Set.preimage_image_eq U E.injective]
    exact hU
  let hh : G₁ ≃ₜ G₂ := E.toEquiv.toHomeomorphOfContinuousOpen hcont hopen
  haveI htg2 : IsTopologicalGroup G₂ := by
    refine Topology.IsInducing.topologicalGroup (E.symm : G₂ →* G₁) ?_
    exact hh.symm.isInducing
  haveI hcs2 : CompactSpace G₂ := hh.compactSpace
  haveI htd2 : TotallyDisconnectedSpace G₂ := hh.totallyDisconnectedSpace
  refine ⟨t2, ?_, htg2, hcs2, htd2, S.image E, ?_⟩
  · intro a b
    obtain ⟨a', rfl⟩ := E.surjective a
    obtain ⟨b', rfl⟩ := E.surjective b
    rw [← map_mul, ← map_mul, mul_comm]
  · have hcl : (Subgroup.closure ((S.image E : Finset G₂) : Set G₂) : Subgroup G₂)
        = Subgroup.map (E : G₁ →* G₂) (Subgroup.closure (S : Set G₁)) := by
      rw [Finset.coe_image, MonoidHom.map_closure]
      rfl
    rw [hcl, Subgroup.coe_map]
    show closure ((E : G₁ → G₂) '' ((Subgroup.closure (S : Set G₁) : Subgroup G₁) : Set G₁))
      = Set.univ
    have himg : closure ((E : G₁ → G₂) '' ((Subgroup.closure (S : Set G₁) : Subgroup G₁) :
        Set G₁)) = (E : G₁ → G₂) '' closure ((Subgroup.closure (S : Set G₁) : Subgroup G₁) :
        Set G₁) := (hh.image_closure _).symm
    rw [himg, hS, Set.image_univ]
    exact E.surjective.range_eq

/-! ## ★5. `Corollary 5.7, (ii)` 本体 -/

/-- ★★★★★**[FrdI] Corollary 5.7, (ii)** —— unit-profinite 型は圏同値で移る。

★★移すのは `Ψ` が単元を保つことだけで、それは**底の 1-可換図式**から出る
（`otimes_map_of_baseSquare`）。位相は `coinduced` で運ぶ。

原文 (FrdI p.108):
> (ii) C1 is of unit-profinite type if and only if C2 is. -/
theorem isOfUnitProfiniteType_map (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
    (ΨBase : D₁ ⥤ D₂) [ΨBase.Faithful]
    (sq : P₁.proj ⋙ ΨBase ≅ Ψ ⋙ P₂.proj)
    (h : IsOfUnitProfiniteType P₁) : IsOfUnitProfiniteType P₂ := by
  intro A₂
  obtain ⟨A₁, ⟨e⟩⟩ : ∃ A₁ : C₁, Nonempty (Ψ.obj A₁ ≅ A₂) :=
    ⟨Ψ.objPreimage A₂, ⟨Ψ.objObjPreimageIso A₂⟩⟩
  obtain ⟨t, hcomm, htg, hcpt, htd, S, hS⟩ := h A₁
  letI := t
  letI : CommGroup ↥(OTimes P₁ A₁) :=
    { (inferInstance : Group ↥(OTimes P₁ A₁)) with mul_comm := hcomm }
  haveI : IsTopologicalGroup ↥(OTimes P₁ A₁) := htg
  haveI : CompactSpace ↥(OTimes P₁ A₁) := hcpt
  haveI : TotallyDisconnectedSpace ↥(OTimes P₁ A₁) := htd
  exact exists_top_of_mulEquiv
    ((otimesEquivOfBaseSquare Ψ ΨBase sq A₁).trans (otimesCongr e)) S hS

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Corollary 5.7, (ii)`。 -/
def isOfUnitProfiniteType_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 108,
    item := "Corollary 5.7, (ii) — unit-profinite 型は圏同値で移る",
    sectionId := "frdi-cor-5-7" }

end Cor57Unit

end ABC3.Found.FrdI
