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

/-! ## ★★★★★★coherence 5 条

★★**recipe**(`μ` / `ε` で確立、5 条すべてに効く):

1. `PresheafOfModules.hom_ext` で成分に落とす
2. `ModuleCat.hom_ext` + `TensorProduct.ext'` で純テンソルへ
3. ★**`show` で両辺を明示的に書き下す**(片側だけでは駄目)
4. `ModuleCat.restrictScalars_μ_tmul` / `restrictScalars_η` で書き換え
5. `rfl`
-/

/-- ★★★★★★**係数変換は lax monoidal である**(前層レベル)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `pushforward φ` が lax monoidal になり、
`Adjunction.leftAdjointOplaxMonoidal` で **`pullback` が oplax monoidal** になる。 -/
noncomputable instance resScalarsLax : (resScalars α).LaxMonoidal where
  ε := resScalarsEps α
  μ M N := resScalarsMu α M N
  μ_natural_left := by
    intro M N f X'
    apply PresheafOfModules.hom_ext
    intro Z
    apply ModuleCat.hom_ext
    apply TensorProduct.ext'
    intro m n
    show Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) (N.obj Z) (X'.obj Z)
        ((f.app Z m) ⊗ₜ n)
      = ((resScalars α).map (f ▷ X')).app Z
          (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom)
            (M.obj Z) (X'.obj Z) (m ⊗ₜ n))
    rw [ModuleCat.restrictScalars_μ_tmul, ModuleCat.restrictScalars_μ_tmul]
    rfl
  μ_natural_right := by
    intro M N X' f
    apply PresheafOfModules.hom_ext
    intro Z
    apply ModuleCat.hom_ext
    apply TensorProduct.ext'
    intro m n
    show Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) (X'.obj Z) (N.obj Z)
        (m ⊗ₜ (f.app Z n))
      = ((resScalars α).map (X' ◁ f)).app Z
          (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom)
            (X'.obj Z) (M.obj Z) (m ⊗ₜ n))
    rw [ModuleCat.restrictScalars_μ_tmul, ModuleCat.restrictScalars_μ_tmul]
    rfl
  associativity := by
    intro M N P
    apply PresheafOfModules.hom_ext
    intro Z
    apply ModuleCat.hom_ext
    apply TensorProduct.ext_threefold
    intro m n p
    show ((resScalars α).map (α_ M N P).hom).app Z
        (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) _ (P.obj Z)
          ((Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom)
            (M.obj Z) (N.obj Z) (m ⊗ₜ n)) ⊗ₜ p))
      = Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) (M.obj Z) _
          (m ⊗ₜ (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom)
            (N.obj Z) (P.obj Z) (n ⊗ₜ p)))
    rw [ModuleCat.restrictScalars_μ_tmul, ModuleCat.restrictScalars_μ_tmul,
      ModuleCat.restrictScalars_μ_tmul, ModuleCat.restrictScalars_μ_tmul]
    rfl
  left_unitality := by
    intro M
    apply PresheafOfModules.hom_ext
    intro Z
    apply ModuleCat.hom_ext
    apply TensorProduct.ext'
    intro r m
    show r • m = ((resScalars α).map (λ_ M).hom).app Z
        (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) _ (M.obj Z)
          ((Functor.LaxMonoidal.ε (ModuleCat.restrictScalars (α.app Z).hom) r) ⊗ₜ m))
    rw [ModuleCat.restrictScalars_η, ModuleCat.restrictScalars_μ_tmul]
    rfl
  right_unitality := by
    intro M
    apply PresheafOfModules.hom_ext
    intro Z
    apply ModuleCat.hom_ext
    apply TensorProduct.ext'
    intro m r
    show r • m = ((resScalars α).map (ρ_ M).hom).app Z
        (Functor.LaxMonoidal.μ (ModuleCat.restrictScalars (α.app Z).hom) (M.obj Z) _
          (m ⊗ₜ (Functor.LaxMonoidal.ε (ModuleCat.restrictScalars (α.app Z).hom) r)))
    rw [ModuleCat.restrictScalars_η, ModuleCat.restrictScalars_μ_tmul]
    rfl

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
