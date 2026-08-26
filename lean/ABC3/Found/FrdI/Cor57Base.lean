/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor57Pair
import ABC3.Found.FrdI.Def27

/-!
# [FrdI] Corollary 5.7, (i) —— base-section は圏同値で移る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.108。

原文 (FrdI p.108):
> (i) Ψ maps base-sections (respectively, quasi-base-Frobenius pairs) of C1

## ★★何を作るか

`S : BaseSection P₁`(= `𝒫₁ ⊆ 𝒞₁`)から `BaseSection P₂` を作る。
★対象は **`Ψ` の像**で取る:

  `objP₂ A₂ ⟺ ∃ A₁, objP₁ A₁ ∧ Ψ A₁ = A₂`

★★**厳密な等号で取るのが要点**である。同型で取ると `𝒫₂` が skeletal でなくなる
(原文の (a)「`𝒫` は skeleton」が壊れる)。
等号で取れば `eqToHom` が入るが、証人を `obtain` して `rfl` で代入すれば消える。

## ★★★鍵は 1 本の補題

**`objP_eq_of_map_eq`** —— `𝒫₁` の 2 対象 `A`, `B` について `Ψ A = Ψ B` なら `A = B`。

★これが無いと `comp_mem`(中間対象が一致するか)も `skeletal` も書けない。
★証明は `𝒫₁ → 𝒟₁` の**充満性と忠実性**から: `Ψ` の充満忠実性で同型 `A ≅ B` を作り、
その底を `fullP` で `𝒫₁` に持ち上げ、合成が恒等であることを `faithfulP` で見て、
`skeletal` を当てる。

## ★仮定は原文が挙げるものそのまま

原文は「`Ψ` が isotropic 対象・prime-Frobenius 射・pull-back 射・birationalization・
射影関手を保つことを示せば足りる。それは `Theorem 3.4, (i), (iii)`;
`Corollary 4.10`; `Corollary 4.11, (ii)` から出る」と書く。
★そこで **pull-back 射の保存**と **Frobenius-trivial 性の保存**、そして
**底の 1-可換図式**(`BaseSquare`)を仮定として受け取る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section Cor57Base

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}
  {Ψ : C₁ ⥤ C₂} {ΨBase : D₁ ⥤ D₂}

/-! ## ★1. 対象の一致 -/

/-- ★★★**`𝒫₁` の 2 対象が `Ψ` で同じ像を持てば等しい**。

★`Ψ` の充満忠実性で同型 `A ≅ B` を作り、その底を `fullP` で `𝒫₁` に持ち上げ、
合成が恒等であることを `faithfulP` で見て、`skeletal` を当てる。 -/
theorem objP_eq_of_map_eq (S : BaseSection P₁) [Ψ.Full] [Ψ.Faithful]
    {A B : C₁} (hA : S.objP A) (hB : S.objP B) (h : Ψ.obj A = Ψ.obj B) : A = B := by
  set e : A ≅ B := Ψ.preimageIso (eqToIso h) with he
  obtain ⟨u, hu, hue⟩ := S.fullP hA hB (P₁.Base e.hom)
  obtain ⟨v, hv, hve⟩ := S.fullP hB hA (P₁.Base e.inv)
  have h1 : u ≫ v = 𝟙 A := by
    refine S.faithfulP (S.comp_mem hu hv) (S.id_mem hA) ?_
    rw [P₁.Base_comp, hue, hve, ← P₁.Base_comp, e.hom_inv_id]
  have h2 : v ≫ u = 𝟙 B := by
    refine S.faithfulP (S.comp_mem hv hu) (S.id_mem hB) ?_
    rw [P₁.Base_comp, hve, hue, ← P₁.Base_comp, e.inv_hom_id]
  exact S.skeletal hu hv h1 h2

/-! ## ★2. 輸送した述語 -/

/-- ★輸送した対象の述語 —— `F` の像。 -/
def objPmap (S : BaseSection P₁) (F : C₁ ⥤ C₂) (A₂ : C₂) : Prop :=
  ∃ A₁ : C₁, S.objP A₁ ∧ F.obj A₁ = A₂

/-- ★輸送した射の述語。 -/
def homPmap (S : BaseSection P₁) (F : C₁ ⥤ C₂) {A₂ B₂ : C₂} (g : A₂ ⟶ B₂) : Prop :=
  ∃ (A₁ B₁ : C₁) (f₁ : A₁ ⟶ B₁), S.objP A₁ ∧ S.objP B₁ ∧ S.homP f₁ ∧
    ∃ (hA : F.obj A₁ = A₂) (hB : F.obj B₁ = B₂),
      F.map f₁ = eqToHom hA ≫ g ≫ eqToHom hB.symm

theorem homPmap_of_map {S : BaseSection P₁} {F : C₁ ⥤ C₂} {A₁ B₁ : C₁} {f₁ : A₁ ⟶ B₁}
    (hA : S.objP A₁) (hB : S.objP B₁) (hf : S.homP f₁) :
    homPmap S F (F.map f₁) :=
  ⟨A₁, B₁, f₁, hA, hB, hf, rfl, rfl, by simp⟩

/-! ## ★3. 輸送した `BaseSection` -/

/-- ★★★★★★**[FrdI] Corollary 5.7, (i)** —— `Ψ` は base-section を base-section へ移す。

★仮定は原文が挙げるもの: 底の 1-可換図式、pull-back 射の保存、
Frobenius-trivial 性の保存。

原文 (FrdI p.108):
> (i) Ψ maps base-sections (respectively, quasi-base-Frobenius pairs) of C1 -/
noncomputable def BaseSection.map (S : BaseSection P₁)
    [Ψ.Full] [Ψ.Faithful] [ΨBase.Full] [ΨBase.Faithful] [ΨBase.EssSurj]
    (bs : BaseSquare Ψ ΨBase P₁ P₂)
    (hpb : ∀ {A B : C₁} (f : A ⟶ B), IsPullBack P₁ f → IsPullBack P₂ (Ψ.map f))
    (hft : ∀ {A : C₁}, S.objP A → IsFrobeniusTrivial P₂ (Ψ.obj A)) :
    BaseSection P₂ where
  objP := objPmap S Ψ
  homP := homPmap S Ψ
  id_mem := by
    rintro A₂ ⟨A₁, hA₁, rfl⟩
    have := homPmap_of_map (F := Ψ) hA₁ hA₁ (S.id_mem hA₁)
    rwa [Ψ.map_id] at this
  comp_mem := by
    rintro A₂ B₂ E₂ f g ⟨A₁, B₁, f₁, hA₁, hB₁, hf₁, rfl, hB, hfe⟩
      ⟨B₁', E₁, g₁, hB₁', hE₁, hg₁, hB', hE, hge⟩
    subst hB
    have hBB : B₁' = B₁ := objP_eq_of_map_eq S hB₁' hB₁ hB'
    subst hBB
    subst hE
    simp only [eqToHom_refl, Category.id_comp, Category.comp_id] at hfe hge
    refine ⟨A₁, E₁, f₁ ≫ g₁, hA₁, hE₁, S.comp_mem hf₁ hg₁, rfl, rfl, ?_⟩
    rw [Ψ.map_comp, hfe, hge]
    simp
  isPullBack := by
    rintro A₂ B₂ g ⟨A₁, B₁, f₁, hA₁, hB₁, hf₁, rfl, rfl, hfe⟩
    simp only [eqToHom_refl, Category.id_comp, Category.comp_id] at hfe
    rw [← hfe]
    exact hpb f₁ (S.isPullBack hf₁)
  skeletal := by
    rintro A₂ B₂ f g ⟨A₁, B₁, f₁, hA₁, hB₁, hf₁, rfl, hB, hfe⟩
      ⟨B₁', A₁', g₁, hB₁', hA₁', hg₁, hB', hA', hge⟩ hfg hgf
    subst hB
    have hBB : B₁ = B₁' := objP_eq_of_map_eq S hB₁ hB₁' hB'.symm
    subst hBB
    have hAA : A₁ = A₁' := objP_eq_of_map_eq S hA₁ hA₁' hA'.symm
    subst hAA
    simp only [eqToHom_refl, Category.id_comp, Category.comp_id] at hfe hge
    subst hfe
    subst hge
    have h1 : f₁ ≫ g₁ = 𝟙 A₁ := by
      refine Ψ.map_injective ?_
      rw [Ψ.map_comp, Ψ.map_id]
      exact hfg
    have h2 : g₁ ≫ f₁ = 𝟙 B₁ := by
      refine Ψ.map_injective ?_
      rw [Ψ.map_comp, Ψ.map_id]
      exact hgf
    exact congrArg Ψ.obj (S.skeletal hf₁ hg₁ h1 h2)
  frobTrivial := by
    rintro A₂ ⟨A₁, hA₁, rfl⟩
    exact hft hA₁
  essSurjP := by
    intro X₂
    obtain ⟨A₁, hA₁, ⟨eA⟩⟩ := S.essSurjP (ΨBase.objPreimage X₂)
    refine ⟨Ψ.obj A₁, ⟨A₁, hA₁, rfl⟩, ⟨?_⟩⟩
    exact (bs.app A₁).symm ≪≫ ΨBase.mapIso eA ≪≫ ΨBase.objObjPreimageIso X₂
  fullP := by
    rintro A₂ B₂ ⟨A₁, hA₁, rfl⟩ ⟨B₁, hB₁, rfl⟩ ψ₂
    obtain ⟨f₁, hf₁, hfe⟩ := S.fullP hA₁ hB₁
      (ΨBase.preimage ((bs.app A₁).hom ≫ ψ₂ ≫ (bs.app B₁).inv))
    refine ⟨Ψ.map f₁, homPmap_of_map hA₁ hB₁ hf₁, ?_⟩
    have hnat := bs.naturality f₁
    rw [hfe, ΨBase.map_preimage] at hnat
    have hL : ((bs.app A₁).hom ≫ ψ₂ ≫ (bs.app B₁).inv) ≫ (bs.app B₁).hom
        = (bs.app A₁).hom ≫ ψ₂ := by simp
    rw [hL] at hnat
    exact ((cancel_epi ((bs.app A₁).hom)).mp hnat).symm
  faithfulP := by
    rintro A₂ B₂ f g ⟨A₁, B₁, f₁, hA₁, hB₁, hf₁, rfl, hB, hfe⟩
      ⟨A₁', B₁', g₁, hA₁', hB₁', hg₁, hA', hB', hge⟩ hbase
    subst hB
    have hAA : A₁ = A₁' := objP_eq_of_map_eq S hA₁ hA₁' hA'.symm
    subst hAA
    have hBB : B₁ = B₁' := objP_eq_of_map_eq S hB₁ hB₁' hB'.symm
    subst hBB
    simp only [eqToHom_refl, Category.id_comp, Category.comp_id] at hfe hge
    subst hfe
    subst hge
    refine congrArg Ψ.map ?_
    refine S.faithfulP hf₁ hg₁ (ΨBase.map_injective ?_)
    have h1 := bs.naturality f₁
    have h2 := bs.naturality g₁
    rw [hbase] at h1
    rw [← h2] at h1
    exact (cancel_mono ((bs.app B₁).hom)).mp h1

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Corollary 5.7, (i)`。 -/
def BaseSection.map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 108,
    item := "Corollary 5.7, (i) — Ψ は base-section を base-section へ移す",
    sectionId := "frdi-cor-5-7" }

end Cor57Base

end ABC3.Found.FrdI
