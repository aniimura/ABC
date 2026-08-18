import ABC3.Found.Arakelov.PicResFreeIso

/-!
# Arakelov (B1) 第 50 ブロック —— **`unit` の計算を `φ` について一般化する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★2 度目の出番

第 31–35 ブロックで、随伴関手定理で作られた**抽象な `unit`** を
生成元の上で計算し切った。★★それは**スキームの射 `f` について**述べてある。

★★★Beck–Chevalley の mate の計算には、**制限した `f|` についても同じもの**が要る。
★★★★どのブロックも証明は **`φ` に何も仮定していない**
——`CorepresentableBy` の一意性と `freeYonedaEquiv` だけである。

★したがって **`φ` について一般に述べ直せば、両方に効く。**

## ★★本ブロックで取れるもの(すべて第 24・31・33・34・35 の一般化)

| 定理 | 内容 |
|---|---|
| `pullbackCorepGen` | ★随伴からの余表現可能性 |
| `pullbackFreeIsoGen` | ★★`pullback φ (free (yoneda Z)) ≅ free (yoneda (F.obj Z))` |
| `corepHomEquivGen` | ★mathlib 側の `homEquiv` は `freeYonedaEquiv` の合成 |
| `unitAppFreeGen` | ★★★`unit` の書き下し |
| `freeYonedaEquivUnitGen` | ★★★`unit` の `freeYonedaEquiv` は `inv` のそれ |
| `isoHomUnitGenGen` | ★★★★★**`unit` は生成元を生成元に送る** |

★★これは**リファクタであって新しい数学ではない**。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory Opposite Limits

variable {C D : Type u} [SmallCategory C] [SmallCategory D]
  {F : C ⥤ D} {R : Dᵒᵖ ⥤ RingCat.{u}} {S : Cᵒᵖ ⥤ RingCat.{u}} (φ : S ⟶ F.op ⋙ R)

/-! ## ★随伴からの余表現可能性 -/

/-- ★**随伴は hom 集合を余表現する**(第 24 ブロックの一般化)。 -/
noncomputable def pullbackCorepGen (M : PresheafOfModules.{u} S) :
    (PresheafOfModules.pushforward φ ⋙ coyoneda.obj (op M)).CorepresentableBy
      ((PresheafOfModules.pullback φ).obj M) where
  homEquiv := (PresheafOfModules.pullbackPushforwardAdjunction φ).homEquiv _ _
  homEquiv_comp g h :=
    (PresheafOfModules.pullbackPushforwardAdjunction φ).homEquiv_naturality_right h g

/-- ★★**生成元の引き戻しは生成元である**(第 24 ブロックの一般化)。 -/
noncomputable def pullbackFreeIsoGen (Z : C) :
    (PresheafOfModules.pullback φ).obj ((PresheafOfModules.free S).obj (yoneda.obj Z))
      ≅ (PresheafOfModules.free R).obj (yoneda.obj (F.obj Z)) :=
  (pullbackCorepGen φ _).uniqueUpToIso
    (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ Z)

/-! ## ★mathlib 側の `homEquiv` -/

/-- ★**mathlib 側の余表現の `homEquiv` は `freeYonedaEquiv` の合成である**(`rfl`)。 -/
theorem corepHomEquivGen (Z : C) (N : PresheafOfModules.{u} R)
    (g : (PresheafOfModules.free R).obj (yoneda.obj (F.obj Z)) ⟶ N) :
    (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ Z).homEquiv g
      = (PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward φ).obj N) (X := Z)).symm
        (PresheafOfModules.freeYonedaEquiv (M := N) (X := F.obj Z) g) := rfl

/-! ## ★★★`unit` の書き下し -/

/-- ★★★**随伴の `unit` を生成元の上で書き下したもの**(第 31 ブロックの一般化)。 -/
theorem unitAppFreeGen (Z : C) :
    (PresheafOfModules.pullbackPushforwardAdjunction φ).unit.app
        ((PresheafOfModules.free S).obj (yoneda.obj Z))
      = (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ Z).homEquiv
        ((pullbackFreeIsoGen φ Z).inv) := by
  show _ = (PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ Z).homEquiv
    ((PresheafOfModules.pushforwardCompCoyonedaFreeYonedaCorepresentableBy φ Z).homEquiv.symm
      ((pullbackCorepGen φ _).homEquiv (𝟙 _)))
  rw [Equiv.apply_symm_apply]
  show _ = (PresheafOfModules.pullbackPushforwardAdjunction φ).homEquiv _ _ (𝟙 _)
  rw [Adjunction.homEquiv_unit]
  simp only [Functor.id_obj, CategoryTheory.Functor.map_id]
  exact (Category.comp_id _).symm

/-- ★★★**`unit` の `freeYonedaEquiv` は同型の `inv` のそれに等しい**(第 34 の一般化)。 -/
theorem freeYonedaEquivUnitGen (Z : C) :
    PresheafOfModules.freeYonedaEquiv
        ((PresheafOfModules.pullbackPushforwardAdjunction φ).unit.app
          ((PresheafOfModules.free S).obj (yoneda.obj Z)))
      = PresheafOfModules.freeYonedaEquiv
        (M := (PresheafOfModules.pullback φ).obj
          ((PresheafOfModules.free S).obj (yoneda.obj Z))) (X := F.obj Z)
        ((pullbackFreeIsoGen φ Z).inv) := by
  rw [unitAppFreeGen, corepHomEquivGen]
  exact Equiv.apply_symm_apply _ _

/-! ## ★★★★★`unit` は生成元を生成元に送る -/

/-- ★★★★★**`unit` は生成元を生成元に送る**(第 35 ブロックの一般化)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `δ`(第 40)と Beck–Chevalley の mate の**両方**に効く。 -/
theorem isoHomUnitGenGen (Z : C) :
    ((pullbackFreeIsoGen φ Z).hom).app (op (F.obj Z))
        (PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward φ).obj
            ((PresheafOfModules.pullback φ).obj
              ((PresheafOfModules.free S).obj (yoneda.obj Z)))) (X := Z)
          ((PresheafOfModules.pullbackPushforwardAdjunction φ).unit.app
            ((PresheafOfModules.free S).obj (yoneda.obj Z))))
      = ModuleCat.freeMk (𝟙 (F.obj Z)) := by
  erw [freeYonedaEquivUnitGen]
  erw [← PresheafOfModules.freeYonedaEquiv_comp, Iso.inv_hom_id, freeYonedaEquiv_id]

/-! ## ★出典の紐付け(`.src`) -/

def isoHomUnitGenGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unit が生成元を生成元に送ること(一般形))",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
