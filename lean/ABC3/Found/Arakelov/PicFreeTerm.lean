import ABC3.Found.Arakelov.PicBCIso

/-!
# Arakelov (B1) 第 55 ブロック —— **`free (yoneda 終対象) ≅ 𝟙_`(一般形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★測ってから一般化した

第 42 ブロック(`free (yoneda ⊤) ≅ 𝟙_`)は **`CompleteLattice α`** で述べてある。
★§9-60 で「**`Over V` は圏としては束そのものではない**ので測ってから使うこと」
と旗を立てた。

★★**測った結果(2026-08-18)**:

| 問い | 結果 |
|---|---|
| `Over V` の終対象は `Over.mk (𝟙 V)` か | ★通った(`Hom` が一意) |
| 束でなく「終対象への `Hom` が一意」だけで足りるか | ★★通った |

★★★したがって **`∀ Z, Unique (Z ⟶ T)`** だけを仮定した形に一般化できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeYonedaTermHom` | ★`free (yoneda T) ⟶ 𝟙_` |
| `freeYonedaTermIso` | ★★★★**`free (yoneda T) ≅ 𝟙_`**(一般形) |
| `overTerminalUnique` | ★`Over V` の終対象への `Hom` は一意 |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

section Term

variable {C : Type u} [SmallCategory C] {R : Cᵒᵖ ⥤ CommRingCat.{u}}
  (T : C) (hT : ∀ Z : C, Unique (Z ⟶ T))

/-- ★**`free (yoneda T)` から単位への射**(生成元を `1` に送る)。 -/
noncomputable def freeYonedaTermHom :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj T)
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

/-- ★★切断ごとの同型——`Hom(Z, T)` は 1 点なので `free` は `R` そのものである。 -/
noncomputable def freeTermSectionIso (U : Cᵒᵖ) :
    (ModuleCat.free ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u)).obj
        ((yoneda.obj T).obj U)
      ≅ ModuleCat.of ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u)
        ((R ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj U : Type u) :=
  (ModuleCat.free _).mapIso (Equiv.toIso
    { toFun := fun _ => PUnit.unit
      invFun := fun _ => (hT U.unop).default
      left_inv := fun _ => ((hT U.unop).uniq _).symm
      right_inv := fun _ => rfl })
    ≪≫ (ModuleCat.FreeMonoidal.εIso _).symm

/-- ★★★★**`free (yoneda T) ≅ 𝟙_`**(一般形)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★仮定は「終対象への `Hom` が一意」だけである——束である必要は無い。 -/
noncomputable def freeYonedaTermIso :
    PresheafOfModules.freeObj (R := R ⋙ forget₂ CommRingCat.{u} RingCat.{u})
        (yoneda.obj T)
      ≅ 𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :=
  haveI happ : ∀ U, IsIso ((freeYonedaTermHom (R := R) T).app U) := by
    intro U
    have h : (freeYonedaTermHom (R := R) T).app U = (freeTermSectionIso (R := R) T hT U).hom := by
      refine ModuleCat.free_hom_ext ?_
      intro x
      erw [ModuleCat.freeDesc_apply]
      erw [CategoryTheory.comp_apply, ModuleCat.free_map_apply,
        ModuleCat.FreeMonoidal.εIso_inv_freeMk]
      rfl
    rw [h]
    infer_instance
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (freeYonedaTermHom (R := R) T)) := by
    rw [NatTrans.isIso_iff_isIso_app]
    intro U
    haveI := happ U
    exact inferInstanceAs
      (IsIso ((forget₂ _ AddCommGrpCat).map ((freeYonedaTermHom (R := R) T).app U)))
  haveI : IsIso (freeYonedaTermHom (R := R) T) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)
  asIso (freeYonedaTermHom (R := R) T)

end Term

/-! ## ★`Over V` の終対象 -/

variable {Y : Scheme.{u}} (V : Y.Opens)

/-- ★**`Over V` の終対象は `Over.mk (𝟙 V)` である**——`Hom` が一意。 -/
@[reducible] noncomputable def overTerminalUnique (Z : Over V) : Unique (Z ⟶ Over.mk (𝟙 V)) where
  default := Over.homMk Z.hom
  uniq _ := Subsingleton.elim _ _

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaTermIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——free (yoneda 終対象) が単位であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
