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


/-- ★★**係数変換の lax monoidal の `ε`**(前層レベル)。

★成分ごとに `ModuleCat.restrictScalars` の `ε`(= `α.app X` そのもの)を束ねる。 -/
noncomputable def resScalarsEps :
    𝟙_ (PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))
      ⟶ (resScalars α).obj (𝟙_ (PresheafOfModules.{u}
          (R' ⋙ forget₂ CommRingCat.{u} RingCat.{u}))) where
  app X := Functor.LaxMonoidal.ε (ModuleCat.restrictScalars (α.app X).hom)
  naturality := by
    intro X Y f
    apply ModuleCat.hom_ext
    apply LinearMap.ext
    intro r
    show Functor.LaxMonoidal.ε (ModuleCat.restrictScalars (α.app Y).hom)
        ((R.map f).hom r)
      = (R'.map f).hom
          (Functor.LaxMonoidal.ε (ModuleCat.restrictScalars (α.app X).hom) r)
    rw [ModuleCat.restrictScalars_η, ModuleCat.restrictScalars_η]
    have h : (R.map f ≫ α.app Y).hom r = (α.app X ≫ R'.map f).hom r :=
      congrArg (fun g : R.obj X ⟶ R'.obj Y => CommRingCat.Hom.hom g r) (α.naturality f)
    exact h

/-! ## ★★残り 5 条(2026-08-17 実測: `cat_disch` では落ちない)

    μ_natural_left / μ_natural_right / associativity /
    left_unitality / right_unitality

★★★いずれも `μ` / `ε` と**同じ手**で落ちる見込みである:

1. 成分に落とす(`PresheafOfModules` の射の等式 → 各 `X` での `ModuleCat` の射の等式)
2. `TensorProduct.ext'` で純テンソルへ
3. ★**`show` で両辺を明示的に書き下す**(片側だけでは駄目——本 turn で 8 回試して確定)
4. `ModuleCat.restrictScalars_μ_tmul` / `restrictScalars_η` で書き換え
5. `rfl`

★これが揃えば `Adjunction.leftAdjointOplaxMonoidal` で
`pullback` が oplax monoidal になり、`PicardData.pullback` が書ける。 -/

/-! ## ★出典の紐付け(`.src`) -/

def resScalarsEps.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数変換の lax monoidal の ε)",
    sectionId := "genell-def-1-1-i" }

def resScalarsMu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——係数変換の lax monoidal の μ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
