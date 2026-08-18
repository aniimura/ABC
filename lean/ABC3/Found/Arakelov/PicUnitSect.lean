import ABC3.Found.Arakelov.PicUnitEnd

/-!
# Arakelov (B2) 第 165 ブロック —— **`unitEnd ∘ unitMul = id`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★`unitEnd` が全射であること

    unitEnd (unitMul c) = c

★これで `End(𝟙_) →+* Γ(X,U)` が全射である。★★単射も出れば環同型になり、
その逆で `Hom(F|_U, 𝟙_)` に `Γ(X,U)` 加群構造が入る。

## ★★逃げ道 3 つ(記録)

| 症状 | 逃げ道 |
|---|---|
| `NatTrans.comp_app` が使えない(`PresheafOfModules.Hom` は `NatTrans` でない) | ★**`PresheafOfModules.comp_app`**(mathlib に有る) |
| `freeYonedaTermIso.hom.app t (freeMk 𝟙) = 1` | ★`erw [ModuleCat.freeDesc_apply]` |
| `exact … _ _ c` が `whnf` で**タイムアウト** | ★★**明示引数**を与えると一瞬で通る |

★★★3 つ目が重要である——`_` を置くと単一化が爆発することがある。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (U : X.Opens)

set_option maxHeartbeats 1000000 in
theorem unitEnd_unitMul (c : (Γ(X, U) : Type u)) : unitEnd U (unitMul U c) = c := by
  show unitVal U (((unitMul U c).app (op (Over.mk (𝟙 U)))).hom (unitOne U)) = c
  unfold unitMul unitHomOfSection
  have h1 : ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom.app
      (op (Over.mk (𝟙 U)))).hom (ModuleCat.freeMk (𝟙 (Over.mk (𝟙 U)))) = unitOne U := by
    show (PresheafOfModules.freeObjDesc _).app _ _ = _
    erw [ModuleCat.freeDesc_apply]
    rfl
  have h2 : ((freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).inv.app
      (op (Over.mk (𝟙 U)))).hom (unitOne U) = ModuleCat.freeMk (𝟙 (Over.mk (𝟙 U))) := by
    rw [← h1, ← ConcreteCategory.comp_apply, ← PresheafOfModules.comp_app,
      (freeYonedaTermIso (Over.mk (𝟙 U)) (overTerminalUnique U)).hom_inv_id,
      PresheafOfModules.id_app]
    rfl
  rw [PresheafOfModules.comp_app, ConcreteCategory.comp_apply, h2]
  exact PresheafOfModules.freeYonedaEquiv_symm_app
    (𝟙_ (PresheafModulesOn X U)) (Over.mk (𝟙 U)) c


/-! ## ★出典の紐付け(`.src`) -/

def unitEnd_unitMul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unitEnd ∘ unitMul = id)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
