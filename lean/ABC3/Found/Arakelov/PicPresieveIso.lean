import ABC3.Found.Arakelov.PicPresieveBij

/-!
# Arakelov (B1) 第 115 ブロック —— **生成族版の同型判定**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★実際に使える形になった

第 114 ブロックで「生成族だけで足りる」ことが出た。
★本ブロックはそれを **3 段**に組み上げる:

| 段 | 内容 |
|---|---|
| 1 | `Sheaf J Ab` 版 |
| 2 | `PresheafModulesOn X V` 版 |
| 3 | 切断から同型を作る版 |

★★これで**基本開集合だけを見れば** `IsLocallyTrivial` が出る。

## ★★本ブロックで取れるもの

| 定理・定義 | 内容 |
|---|---|
| `isIso_of_bijective_on_presieve` | ★★★★層の射の同型判定 |
| `isIso_of_bijective_on_presieve_mod` | ★★★★前層加群の射の同型判定 |
| `trivialIsoOfSectionPresieve` | ★★★★★★**切断から `𝟙_ ≅ P`** |
-/

universe v u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

section Sheaf

variable {C : Type u} [Category.{v} C] (J : GrothendieckTopology C)
  [HasWeakSheafify J AddCommGrpCat.{v}] [J.WEqualsLocallyBijective AddCommGrpCat.{v}]

/-- ★★★★**生成族の上で全単射なら層の射は同型である**。 -/
theorem isIso_of_bijective_on_presieve {F G : Sheaf J AddCommGrpCat.{v}} (f : F ⟶ G)
    (h : ∀ U : C, ∃ R : Presieve U, Sieve.generate R ∈ J U ∧
      ∀ ⦃W : C⦄ (i : W ⟶ U), R i → Function.Bijective (f.hom.app (op W))) :
    IsIso f := by
  haveI : Presheaf.IsLocallySurjective J f.hom := by
    refine isLocallySurjective_of_presieve J f.hom (fun U s => ?_)
    obtain ⟨R, hR, hb⟩ := h U
    exact ⟨R, hR, fun {W} i hi => (hb i hi).2 _⟩
  haveI : Presheaf.IsLocallyInjective J f.hom := by
    refine isLocallyInjective_of_presieve J f.hom (fun U x y hxy => ?_)
    obtain ⟨R, hR, hb⟩ := h U
    refine ⟨R, hR, fun {W} i hi => ?_⟩
    refine (hb i hi).1 ?_
    rw [NatTrans.naturality_apply f.hom i.op x, NatTrans.naturality_apply f.hom i.op y, hxy]
  exact (Sheaf.isLocallyBijective_iff_isIso f).1 ⟨inferInstance, inferInstance⟩

end Sheaf

section Mod

variable {X : Scheme.{u}} {V : X.Opens} {P Q : PresheafModulesOn X V}

/-- ★★★★**生成族の上で全単射なら前層加群の射は同型である**。 -/
theorem isIso_of_bijective_on_presieve_mod
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) P.presheaf)
    (hQ : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) Q.presheaf)
    (φ : P ⟶ Q)
    (h : ∀ W : Over V, ∃ R : Presieve W,
      Sieve.generate R ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), R i →
        Function.Bijective ((PresheafOfModules.toPresheaf _).map φ |>.app (op Z))) :
    IsIso φ := by
  haveI : IsIso ((⟨(PresheafOfModules.toPresheaf _).map φ⟩ :
      (⟨P.presheaf, hP⟩ : Sheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u})
        ⟶ ⟨Q.presheaf, hQ⟩)) :=
    isIso_of_bijective_on_presieve _ _ h
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map φ) := by
    have := (sheafToPresheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u}).map_isIso
      (⟨(PresheafOfModules.toPresheaf _).map φ⟩ :
      (⟨P.presheaf, hP⟩ : Sheaf ((Opens.grothendieckTopology X).over V) AddCommGrpCat.{u})
        ⟶ ⟨Q.presheaf, hQ⟩)
    exact this
  exact isIso_of_reflects_iso φ (PresheafOfModules.toPresheaf _)

/-- ★★★★★★**切断から `𝟙_ ≅ P`**(生成族版)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `IsLocallyTrivial` の同型そのものであり、
**基本開集合だけを見れば済む**形になっている。 -/
noncomputable def trivialIsoOfSectionPresieve (V : X.Opens) (P : PresheafModulesOn X V)
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) P.presheaf)
    (h1 : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      (𝟙_ (PresheafModulesOn X V)).presheaf)
    (s : P.obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, ∃ R : Presieve W,
      Sieve.generate R ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), R i →
        Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op Z) : Type u) =>
          c • restrictSec V P s Z)) :
    𝟙_ (PresheafModulesOn X V) ≅ P :=
  haveI : IsIso (unitHomOfSection V P s) := by
    refine isIso_of_bijective_on_presieve_mod h1 hP _ (fun W => ?_)
    obtain ⟨R, hR, hb⟩ := h W
    exact ⟨R, hR, fun {Z} i hi => bijective_unitHomOfSection_app V P s Z (hb i hi)⟩
  asIso (unitHomOfSection V P s)

end Mod

/-! ## ★出典の紐付け(`.src`) -/

def trivialIsoOfSectionPresieve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成族版の同型判定)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
