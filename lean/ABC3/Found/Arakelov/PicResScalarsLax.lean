import ABC3.Found.Arakelov.PicSheafGroup
import Mathlib.Algebra.Category.ModuleCat.Monoidal.Adjunction

/-!
# Arakelov (B1) 第 18 ブロック —— **係数変換の lax monoidal(前層レベル)** (`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが引き戻し `f^*` への入口である

    PresheafOfModules.pullback φ := (pushforward φ).leftAdjoint
    pushforward φ = pushforward₀ F R ⋙ restrictScalars φ

★`pushforward₀` は mathlib で **strong monoidal**、
★★`restrictScalars` の lax monoidal は **ModuleCat レベルにしかない**。
★★★本ブロックはそれを**前層レベルへ持ち上げる**第一歩(`μ`)である。

## ★機構

成分ごとに `ModuleCat.restrictScalars` の `μ` を取り、自然性を純テンソルで確かめる。
★★★**`ModuleCat.restrictScalars_μ_tmul`(μ は純テンソル上で恒等)**が鍵。
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite

variable {C : Type u} [Category.{u} C] {R R' : Cᵒᵖ ⥤ CommRingCat.{u}} (α : R ⟶ R')

/-- ★係数変換(`CommRingCat` 版)。 -/
noncomputable abbrev resScalars :
    PresheafOfModules.{u} (R' ⋙ forget₂ CommRingCat.{u} RingCat.{u}) ⥤
      PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}) :=
  PresheafOfModules.restrictScalars
    (Functor.whiskerRight α (forget₂ CommRingCat.{u} RingCat.{u}))

/-- ★★★★**係数変換の lax monoidal の `μ`**(前層レベル)。

★★★成分ごとの `ModuleCat.restrictScalars` の `μ` を束ねる。
自然性は純テンソルで `ModuleCat.restrictScalars_μ_tmul`(μ は純テンソル上で恒等)で落ちる。 -/
noncomputable def resScalarsMu
    (M N : PresheafOfModules.{u} (R' ⋙ forget₂ CommRingCat.{u} RingCat.{u})) :
    (resScalars α).obj M ⊗ (resScalars α).obj N ⟶ (resScalars α).obj (M ⊗ N) where
  app X := Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app X).hom) (M.obj X) (N.obj X)
  naturality := by
    intro X Y f
    apply ModuleCat.hom_ext
    apply TensorProduct.ext'
    intro m n
    show Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Y).hom) (M.obj Y) (N.obj Y)
        ((M.map f m) ⊗ₜ (N.map f n))
      = ((resScalars α).obj (M ⊗ N)).map f
          (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app X).hom)
            (M.obj X) (N.obj X) (m ⊗ₜ n))
    rw [ModuleCat.restrictScalars_μ_tmul, ModuleCat.restrictScalars_μ_tmul]
    rfl


/-! ## ★出典の紐付け(`.src`) -/

def resScalarsMu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数変換の lax monoidal の μ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
