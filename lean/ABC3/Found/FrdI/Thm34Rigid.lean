import ABC3.Found.FrdI.Thm34

/-!
# [FrdI] Theorem 3.4, (i) の最終文 —— isotropification は rigid

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.62–p.63。

原文 (FrdI p.62):
> functors of this diagram is rigid.

## ★★原文の証明(測定、2026-08-17)

原文 (FrdI p.63):
> By Proposition 1.13, (ii), it suffices to show, for each A ∈Ob(Cistr) that the

★★**段は 3 つ**:

| 段 | 内容 | 道具 |
|---|---|---|
| A | `η` の成分は `𝒪^×` に属する | ★`Proposition 1.13, (ii)`(`𝒟` が slim) |
| B | `A` が **Frobenius-trivial** なら成分は恒等 | ★**Frobenius-normalized の関係式 ＋ 全射性** |
| C | 一般の `A` へ移す | ★`Definition 1.3, (i), (a), (b)` の**span** |

## ★★★段 A で原文が使う「`𝒞 → 𝔽_Φ` が rigid」の当て方

★★`isotropification` と `𝔽_Φ` への関手の合成は **`𝒞 → 𝔽_Φ` と同型**である ——
`hullMap A : A ⟶ A^istr` は **isometric pre-step** なので、その `𝔽_Φ` での像は
**同型**(底同型・零因子 0・次数 1)であり、四角形は `isotropification_square`。
★したがって `isRigidFunctor_of_iso` で rigidity が移る。

## ★★★段 C の span

★`Definition 1.3, (i), (a)`(`baseSurj`)で `Base A` を実現する
**Frobenius-trivial 対象 `B₀`** を取り、(i)(b)(`preStepSpan`)で
**pre-step の span `B₀ ← X → A`** を得る。

★★**左脚 `φ : X ⟶ B₀` は pre-step ゆえ mono**(`preStepMono`)なので
`η_X = 𝟙` が出る。★**右脚 `ψ : X ⟶ A` は epi**(totally epimorphic)なので
`η_A = 𝟙` が出る。★2 本の脚を**別々の理由**で消すのが要点である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★段 A の下ごしらえ —— hull の `𝔽_Φ` での像は同型 -/

include P in
/-- ★**isometric pre-step の `𝔽_Φ` での像は同型**。

★底同型(pre-step)・零因子 0(isometric)・次数 1(linear)の 3 つが揃う。 -/
theorem isIso_toElem_map_of_isometricPreStep {A B : C} (φ : A ⟶ B)
    (hi : IsIsometric P φ) (hs : IsPreStep P φ) : IsIso (P.toElem.map φ) := by
  refine (ElemFrobCat.isIso_iff (P.toElem.map φ)).mpr ⟨hs.2, ?_, hs.1⟩
  rw [show (P.toElem.map φ).div = P.Div φ from rfl, hi]
  exact isAddUnit_zero

include P in
theorem isIso_toElem_map_hullMap (F : FrobenioidCore P) (A : C) :
    IsIso (P.toElem.map (hullMap P F A)) :=
  isIso_toElem_map_of_isometricPreStep P (hullMap P F A)
    (hullMap_spec P F A).1 (hullMap_spec P F A).2.1

/-- ★hull の `𝔽_Φ` での像(同型として)。 -/
noncomputable def hullElemIso (F : FrobenioidCore P) (A : C) :
    P.toElem.obj A ≅ P.toElem.obj (hullObj P F A) :=
  @asIso _ _ _ _ _ (isIso_toElem_map_hullMap P F A)

@[simp] theorem hullElemIso_hom (F : FrobenioidCore P) (A : C) :
    (hullElemIso P F A).hom = P.toElem.map (hullMap P F A) := rfl

/-- ★★**`𝒞 → 𝔽_Φ` と `𝒞 → 𝒞^istr → 𝔽_Φ` は同型な関手**。

★成分は `hullMap` の像(上で同型)、四角形は `isotropification_square`。 -/
noncomputable def isotropificationElemIso (F : FrobenioidCore P) :
    P.toElem ≅ isotropification P F ⋙ (isotropicProp P).ι ⋙ P.toElem :=
  NatIso.ofComponents (hullElemIso P F)
    (fun {A B} f => by
      show P.toElem.map f ≫ P.toElem.map (hullMap P F B)
        = P.toElem.map (hullMap P F A) ≫ P.toElem.map (istrMap P F f)
      rw [← P.toElem.map_comp, ← P.toElem.map_comp, isotropification_square])

include P in
/-- ★★**段 A** —— `𝒞 → 𝒞^istr → 𝔽_Φ` は rigid。 -/
theorem isRigidFunctor_isotropification_comp_toElem (F : FrobenioidCore P)
    (hslim : IsSlimCat D) :
    IsRigidFunctor (isotropification P F ⋙ (isotropicProp P).ι ⋙ P.toElem) :=
  isRigidFunctor_of_iso (isotropificationElemIso P F).symm
    (prop_1_13_ii_global P F hslim)

include P in
/-- ★★**段 A の帰結** —— `isotropification` の自己同型の成分は `𝒪^×` に属する。 -/
theorem isotropification_aut_mem_otimes (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (η : isotropification P F ≅ isotropification P F) (A : C) :
    ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A) ∈ OTimes P (hullObj P F A) := by
  have h := isRigidFunctor_isotropification_comp_toElem P F hslim
    (Functor.isoWhiskerRight η ((isotropicProp P).ι ⋙ P.toElem))
  have hc : P.toElem.map ((η.hom.app A).hom) = 𝟙 (P.toElem.obj (hullObj P F A)) :=
    congrArg (fun t => t.hom.app A) h
  have hinv : ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A)
      ≫ ((η.inv.app A).hom : hullObj P F A ⟶ hullObj P F A) = 𝟙 _ :=
    congrArg (fun t : hullIstr P F A ⟶ hullIstr P F A => t.hom) (η.hom_inv_id_app A)
  have hinv' : ((η.inv.app A).hom : hullObj P F A ⟶ hullObj P F A)
      ≫ ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A) = 𝟙 _ :=
    congrArg (fun t : hullIstr P F A ⟶ hullIstr P F A => t.hom) (η.inv_hom_id_app A)
  refine ⟨⟨?_, ?_⟩, ?_⟩
  · show P.Base ((η.hom.app A).hom) = P.Base (𝟙 _)
    show (P.toElem.map ((η.hom.app A).hom)).base = (P.toElem.map (𝟙 _)).base
    rw [hc, P.toElem.map_id]
    rfl
  · show P.degFr ((η.hom.app A).hom) = 1
    show (P.toElem.map ((η.hom.app A).hom)).deg = 1
    rw [hc]; rfl
  · refine (CategoryTheory.isUnit_iff_isIso
      ((η.hom.app A).hom : End (hullObj P F A))).mpr ?_
    exact ⟨⟨((η.inv.app A).hom : hullObj P F A ⟶ hullObj P F A), hinv, hinv'⟩⟩

/-! ## ★段 B —— `A` が Frobenius-trivial なら成分は恒等

原文 (FrdI p.63):
> α follows from the functoriality of α with respect to base-identity endomorphisms

★★**自然性**が `F(ζ) ≫ η = η ≫ F(ζ)`、**Frobenius-normalized** が
`F(ζ) ≫ η^n = η ≫ F(ζ)` を与える。★両者を合わせ、`F(ζ)` が epi なので
`η^n = η`。★`n = 2` に取り、`η` が同型なので `η = 𝟙`。 -/

include P in
/-- ★★★**段 B** —— `A` が Frobenius-trivial かつ `A^istr` が Frobenius-normalized なら
`η` の `A` 成分は恒等。 -/
theorem isotropification_aut_eq_id_of_frobTrivial (F : FrobenioidCore P)
    (hslim : IsSlimCat D) (η : isotropification P F ≅ isotropification P F) (A : C)
    (hft : IsFrobeniusTrivial P A) (hfn : IsFrobeniusNormalized P (hullObj P F A)) :
    ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A) = 𝟙 (hullObj P F A) := by
  obtain ⟨ζ, hζd, hζb⟩ := hft
  set α : End (hullObj P F A) := (η.hom.app A).hom with hα
  have hmem : α ∈ OTimes P (hullObj P F A) := isotropification_aut_mem_otimes P F hslim η A
  set θ : hullObj P F A ⟶ hullObj P F A := istrMap P F ((ζ 2 : A ⟶ A)) with hθ
  have hθb : IsBaseIdentity P θ :=
    isotropification_baseIdentity P F (ζ 2 : A ⟶ A) (hζb 2).1
  have hθd : P.degFr θ = 2 := by
    rw [hθ, isotropification_degFr]; exact hζd 2
  have hnat : θ ≫ (α : hullObj P F A ⟶ hullObj P F A)
      = (α : hullObj P F A ⟶ hullObj P F A) ≫ θ :=
    congrArg (fun t : hullIstr P F A ⟶ hullIstr P F A => t.hom)
      (η.hom.naturality (ζ 2 : A ⟶ A))
  have hnorm := hfn θ hθb α (OTimes_le_OTri P _ hmem)
  rw [hθd, show ((2 : ℕ+) : ℕ) = 2 from rfl] at hnorm
  haveI : Epi θ := P.totEpiC _ _ θ
  have hpow : ((α ^ 2 : End (hullObj P F A)) : hullObj P F A ⟶ hullObj P F A)
      = (α : hullObj P F A ⟶ hullObj P F A) :=
    (cancel_epi θ).mp (hnorm.trans hnat.symm)
  have hsq : (α : hullObj P F A ⟶ hullObj P F A) ≫ (α : hullObj P F A ⟶ hullObj P F A)
      = (α : hullObj P F A ⟶ hullObj P F A) := by
    refine Eq.trans ?_ hpow
    show _ = ((α ^ 2 : End (hullObj P F A)) : hullObj P F A ⟶ hullObj P F A)
    rw [pow_two]
    rfl
  haveI : Epi (α : hullObj P F A ⟶ hullObj P F A) := P.totEpiC _ _ _
  exact (cancel_epi (α : hullObj P F A ⟶ hullObj P F A)).mp
    (hsq.trans (Category.comp_id _).symm)

/-! ## ★段 C —— span で一般の対象へ移す -/

include P in
/-- ★★**`Definition 1.3, (i), (a), (b)`** —— 任意の対象は
**Frobenius-trivial な対象と pre-step の span で結ばれる**。 -/
theorem exists_preStepSpan_frobTrivial (F : FrobenioidCore P) (A : C) :
    ∃ (B₀ X : C) (φ : X ⟶ B₀) (ψ : X ⟶ A), IsFrobeniusTrivial P B₀ ∧
      IsPreStep P φ ∧ IsPreStep P ψ := by
  obtain ⟨B₀, hft, ⟨e⟩⟩ := F.baseSurj ((P.toElem.obj A).base)
  obtain ⟨X, φ, ψ, hφ, hψ, -⟩ := F.preStepSpan B₀ A e.hom inferInstance
  exact ⟨B₀, X, φ, ψ, hft, hφ, hψ⟩

include P in
/-- ★★★★**[FrdI] Theorem 3.4, (i) の最終文(1 つの Frobenioid の場合)** ——
`𝒟` が slim で `𝒞` が Frobenius-normalized 型なら
**isotropification 関手は rigid** である。

原文 (FrdI p.62):
> functors of this diagram is rigid.

★段 A で成分を `𝒪^×` に落とし、段 B で Frobenius-trivial な対象で恒等にし、
★段 C で **pre-step の span**(左脚は mono、右脚は epi)で一般の対象へ移す。 -/
theorem isRigidFunctor_isotropification (F : FrobenioidCore P) (hslim : IsSlimCat D)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X) :
    IsRigidFunctor (isotropification P F) := by
  intro η
  refine Iso.ext (NatTrans.ext (funext fun A => ?_))
  refine InducedCategory.hom_ext ?_
  show ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A) = 𝟙 (hullObj P F A)
  obtain ⟨B₀, X, φ, ψ, hft, hφ, hψ⟩ := exists_preStepSpan_frobTrivial P F A
  have hB₀ : ((η.hom.app B₀).hom : hullObj P F B₀ ⟶ hullObj P F B₀)
      = 𝟙 (hullObj P F B₀) :=
    isotropification_aut_eq_id_of_frobTrivial P F hslim η B₀ hft (hfn _)
  haveI : Mono (istrMap P F φ) :=
    F.preStepMono (istrMap P F φ) ((isotropification_preStep_iff P F φ).mpr hφ)
  have hnatφ : istrMap P F φ ≫ ((η.hom.app B₀).hom : hullObj P F B₀ ⟶ hullObj P F B₀)
      = ((η.hom.app X).hom : hullObj P F X ⟶ hullObj P F X) ≫ istrMap P F φ :=
    congrArg (fun t : hullIstr P F X ⟶ hullIstr P F B₀ => t.hom) (η.hom.naturality φ)
  have hX : ((η.hom.app X).hom : hullObj P F X ⟶ hullObj P F X) = 𝟙 (hullObj P F X) := by
    have hrhs : istrMap P F φ ≫ ((η.hom.app B₀).hom : hullObj P F B₀ ⟶ hullObj P F B₀)
        = 𝟙 (hullObj P F X) ≫ istrMap P F φ := by
      rw [hB₀]
      exact (Category.comp_id _).trans (Category.id_comp _).symm
    exact (cancel_mono (istrMap P F φ)).mp (hnatφ.symm.trans hrhs)
  haveI : Epi (istrMap P F ψ) := P.totEpiC _ _ _
  have hnatψ : istrMap P F ψ ≫ ((η.hom.app A).hom : hullObj P F A ⟶ hullObj P F A)
      = ((η.hom.app X).hom : hullObj P F X ⟶ hullObj P F X) ≫ istrMap P F ψ :=
    congrArg (fun t : hullIstr P F X ⟶ hullIstr P F A => t.hom) (η.hom.naturality ψ)
  have hrhs : ((η.hom.app X).hom : hullObj P F X ⟶ hullObj P F X) ≫ istrMap P F ψ
      = istrMap P F ψ ≫ 𝟙 (hullObj P F A) := by
    rw [hX]
    exact (Category.id_comp _).trans (Category.comp_id _).symm
  exact (cancel_epi (istrMap P F ψ)).mp (hnatψ.trans hrhs)

/-! ## ★★★★2 つの Frobenioid —— 図式の 2 本の合成関手が rigid

原文 (FrdI p.62):
> functors of this diagram is rigid.

★★**片方は `isRigidFunctor_comp_of_isEquivalence`、もう片方は 1-可換性で移す。**
★これがまさに `Thm34.lean` の設計メモが言う「2 つの道具」である。 -/

section TwoFrob

universe v₁ u₁ v₂ u₂ v₃ u₃ w₃ v₄ u₄ w₄

variable {C₁ : Type u₁} [Category.{v₁} C₁] {C₂ : Type u₂} [Category.{v₂} C₂]
  (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
  {E₁ : Type u₃} [Category.{v₃} E₁] {Φ₁ : MonoidOn.{v₃, u₃, w₃} E₁}
  (P₁ : PreFrobenioid C₁ Φ₁)
  {E₂ : Type u₄} [Category.{v₄} E₂] {Φ₂ : MonoidOn.{v₄, u₄, w₄} E₂}
  (P₂ : PreFrobenioid C₂ Φ₂)

/-- ★★**図式の「右回り」の合成 `𝒞₁ → 𝒞₂ → 𝒞₂^istr` は rigid**。

★`Ψ` が圏同値なので `isRigidFunctor_comp_of_isEquivalence` がそのまま効く。 -/
theorem isRigidFunctor_psi_comp_isotropification (F₂ : FrobenioidCore P₂)
    (hslim₂ : IsSlimCat E₂) (hfn₂ : ∀ X : C₂, IsFrobeniusNormalized P₂ X) :
    IsRigidFunctor (Ψ ⋙ isotropification P₂ F₂) :=
  isRigidFunctor_comp_of_isEquivalence Ψ (isotropification P₂ F₂)
    (isRigidFunctor_isotropification P₂ F₂ hslim₂ hfn₂)

/-- ★★★★**[FrdI] Theorem 3.4, (i) の最終文** ——
`𝒟₁`・`𝒟₂` が slim で `𝒞₁`・`𝒞₂` が Frobenius-normalized 型なら、
**図式の 2 本の合成関手はいずれも rigid** である。

原文 (FrdI p.62):
> functors of this diagram is rigid.

★「左回り」`𝒞₁ → 𝒞₁^istr → 𝒞₂^istr` は 1-可換性(`isotropificationCommute`)で
「右回り」と同型なので、`isRigidFunctor_of_iso` で移る。 -/
theorem isRigidFunctor_isotropification_comp_psiIstr (F₁ : FrobenioidCore P₁)
    (F₂ : FrobenioidCore P₂) (h₁ : IsOfQuasiIsotropicType C₁ P₁)
    (h₂ : IsOfQuasiIsotropicType C₂ P₂) (hslim₂ : IsSlimCat E₂)
    (hfn₂ : ∀ X : C₂, IsFrobeniusNormalized P₂ X) :
    IsRigidFunctor (isotropification P₁ F₁ ⋙ psiIstr Ψ P₁ P₂ h₁ h₂) :=
  isRigidFunctor_of_iso (isotropificationCommute Ψ P₁ P₂ F₁ F₂ h₁ h₂)
    (isRigidFunctor_psi_comp_isotropification Ψ P₂ F₂ hslim₂ hfn₂)

end TwoFrob

/-! ## ★`Theorem 3.4, (i)` の 4 主張の対応

| 主張 | 実装 |
|---|---|
| `Ψ` は isotropic 対象を保つ | `thm_3_4_i_isotropic`(`Thm34.lean`) |
| `Ψ` は isotropic hull・isometric pre-step を保つ | `isotropicHull_map` / `isometricPreStep_map`(`Thm34.lean`) |
| 1-一意な `Ψ^istr` が 1-可換図式に入る | `psiIstr` / `isotropificationCommute` / `isotropificationCommute_uniq`(`Thm34.lean`) |
| **合成関手はいずれも rigid** | ★**本ファイル** |
-/

/-- ★★★★**[FrdI] Theorem 3.4, (i)** —— 4 主張が実装された。

★★(ii)–(v) は未実装なので**条つき**で記録する。 -/
def thm_3_4_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62, item := "Theorem 3.4, (i)",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
