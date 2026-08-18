import ABC3.Found.Arakelov.PicSheafPullTensor

/-!
# Arakelov (B1) 第 42 ブロック —— **前層の段でも `f^* 𝒪 ≅ 𝒪`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★第 24 ブロックを `⊤` に当てる

    free (yoneda ⊤) の切断 = free (Hom(U, ⊤)) = free (1 点) ≅ R(U)
      ⟹ **free (yoneda ⊤) ≅ 𝟙_**(構造前層そのもの)

★`f⁻¹⊤ = ⊤` なので、第 24 ブロック(`f^*(free (yoneda V)) ≅ free (yoneda f⁻¹V)`)を
`V := ⊤` に当てると

    **f^*_pre 𝒪_Y ≅ 𝒪_X**(前層の段)

★★これが第 41 ブロックの仮定(局所階数 1)を示すための出発点である
——局所自明性は「局所的に `𝒪` と同型」だから。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeYonedaTopHom` | ★`free (yoneda ⊤) ⟶ 𝟙_` |
| `freeYonedaTopIso` | ★★★**`free (yoneda ⊤) ≅ 𝟙_`** |
| `opensMap_top` | ★`f⁻¹⊤ = ⊤` |
| `pullbackUnitPreIso` | ★★★★★**`f^*_pre 𝒪_Y ≅ 𝒪_X`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

/-! ## ★`free (yoneda ⊤)` は単位である -/

section Top

variable {α : Type u} [CompleteLattice α] {R : αᵒᵖ ⥤ CommRingCat.{u}}

/-- ★**`free (yoneda ⊤)` から単位への射**(生成元を `1` に送る)。 -/
noncomputable def freeYonedaTopHom :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj (⊤ : α))
      ⟶ 𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :=
  PresheafOfModules.freeObjDesc
    { app := fun U => ↾fun _ => show
        ((𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))).presheaf
          ⋙ forget Ab).obj U from
        (1 : ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u))
      naturality := by
        intro U V φ
        ext x
        exact (map_one ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).map φ).hom).symm }

/-- ★★切断ごとの同型——`Hom(U, ⊤)` は 1 点なので `free` は `R` そのものである。 -/
noncomputable def freeTopSectionIso (U : αᵒᵖ) :
    (ModuleCat.free ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u)).obj
        ((yoneda.obj (⊤ : α)).obj U)
      ≅ ModuleCat.of ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u)
        ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u) :=
  (ModuleCat.free _).mapIso (Equiv.toIso
    { toFun := fun _ => PUnit.unit
      invFun := fun _ => homOfLE le_top
      left_inv := fun _ => Subsingleton.elim (α := (unop U ⟶ (⊤ : α))) _ _
      right_inv := fun _ => rfl })
    ≪≫ (ModuleCat.FreeMonoidal.εIso _).symm

/-- ★★★**`free (yoneda ⊤) ≅ 𝟙_`**——構造前層そのものである。 -/
noncomputable def freeYonedaTopIso :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj (⊤ : α))
      ≅ 𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :=
  haveI happ : ∀ U, IsIso ((freeYonedaTopHom (R := R) (α := α)).app U) := by
    intro U
    have h : (freeYonedaTopHom (R := R) (α := α)).app U = (freeTopSectionIso (R := R) U).hom := by
      refine ModuleCat.free_hom_ext ?_
      intro x
      erw [ModuleCat.freeDesc_apply]
      erw [CategoryTheory.comp_apply, ModuleCat.free_map_apply,
        ModuleCat.FreeMonoidal.εIso_inv_freeMk]
      rfl
    rw [h]
    infer_instance
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (freeYonedaTopHom (R := R) (α := α))) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro U
    haveI := happ U
    exact inferInstanceAs
      (IsIso ((forget₂ _ AddCommGrpCat).map ((freeYonedaTopHom (R := R) (α := α)).app U)))
  haveI : IsIso (freeYonedaTopHom (R := R) (α := α)) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)
  asIso (freeYonedaTopHom (R := R) (α := α))

end Top

/-! ## ★★★★★前層の段の `f^* 𝒪 ≅ 𝒪` -/

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★**逆像は全体を全体に送る**。 -/
theorem opensMap_top : (Opens.map f.base).obj (⊤ : Y.Opens) = (⊤ : X.Opens) := by
  ext x
  simp

/-- ★★★★★★**前層の段でも `f^* 𝒪_Y ≅ 𝒪_X`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが第 41 ブロックの仮定(局所階数 1)を示すための出発点である。
★機構は第 24 ブロックを `V := ⊤` に当て、両端で `free (yoneda ⊤) ≅ 𝟙_` を使う。 -/
noncomputable def pullbackUnitPreIso :
    (pullbackPre f).obj (𝟙_ Y.PresheafOfModules) ≅ 𝟙_ X.PresheafOfModules :=
  (pullbackPre f).mapIso (freeYonedaTopIso (R := Y.presheaf) (α := Y.Opens)).symm
    ≪≫ pullbackFreeYonedaIso f ⊤
    ≪≫ eqToIso (by rw [opensMap_top]; rfl)
    ≪≫ freeYonedaTopIso (R := X.presheaf) (α := X.Opens)

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaTopIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——free (yoneda ⊤) が単位であること)",
    sectionId := "genell-def-1-1-i" }

def pullbackUnitPreIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——前層の段でも f^* 𝒪 ≅ 𝒪 であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
