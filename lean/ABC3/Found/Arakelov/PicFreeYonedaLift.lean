import ABC3.Found.Arakelov.PicColimitLift
import Mathlib.Algebra.Category.ModuleCat.Presheaf.Generator

/-!
# Arakelov (B1) 第 29 ブロック —— **生成元で同型なら全対象で同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★これが「持ち上げ」である

mathlib は各前層加群 `M` を

    (kernel の余積) → (M の余積) → M

という**余積の余核**として表示する(`isColimitFreeYonedaCoproductsCokernelCofork`)。
★★したがって第 28 ブロックの器具を **2 回**当てればよい:

| 回 | 図式 | 得るもの |
|---|---|---|
| 1 | `Discrete (M.Elements)` | ★余積の上で同型 |
| 2 | `WalkingParallelPair` | ★★**余核の上で同型**——すなわち `M` の上で同型 |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isIso_app_freeYonedaCoproduct` | ★余積の上で同型 |
| `isIso_app_of_freeYoneda` | ★★★★★**生成元で同型なら全対象で同型** |

★★★これで `δ` の同型性は**生成元の 1 点だけ**に帰着した。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory Limits PresheafOfModules

variable {C : Type u} [SmallCategory C] {R : Cᵒᵖ ⥤ RingCat.{u}}
  {D : Type u} [Category.{u} D]
  {A B : PresheafOfModules.{u} R ⥤ D} (τ : A ⟶ B)
  [PreservesColimitsOfSize.{u, u} A] [PreservesColimitsOfSize.{u, u} B]

/-! ## ★余積の上で -/

/-- ★★**余積の上で同型**——第 28 ブロックを `Discrete` 図式に当てる。 -/
theorem isIso_app_freeYonedaCoproduct
    (h : ∀ (X : C), IsIso (τ.app ((PresheafOfModules.free R).obj (yoneda.obj X))))
    (M : PresheafOfModules.{u} R) : IsIso (τ.app M.freeYonedaCoproduct) :=
  isIso_app_of_colimit τ (Discrete.functor (Elements.freeYoneda (M := M)))
    (fun j => h j.as.1.unop)

/-! ## ★★★★★生成元から全対象へ -/

/-- ★★★★★★**生成元で同型なら全対象で同型**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `δ` の同型性を**生成元の 1 点**に帰着させる器具である。 -/
theorem isIso_app_of_freeYoneda
    (h : ∀ (X : C), IsIso (τ.app ((PresheafOfModules.free R).obj (yoneda.obj X))))
    (M : PresheafOfModules.{u} R) : IsIso (τ.app M) := by
  haveI := isIso_app_freeYonedaCoproduct τ h
  refine isIso_app_of_isColimit τ (parallelPair M.toFreeYonedaCoproduct 0) _
    M.isColimitFreeYonedaCoproductsCokernelCofork ?_
  rintro (_ | _)
  · exact isIso_app_freeYonedaCoproduct τ h _
  · exact isIso_app_freeYonedaCoproduct τ h M

/-! ## ★出典の紐付け(`.src`) -/

def isIso_app_of_freeYoneda.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元で同型なら全対象で同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
