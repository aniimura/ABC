import ABC3.Found.Arakelov.PicUnitCalc

/-!
# Arakelov (B1) 第 35 ブロック —— **`unit` は生成元を生成元に送る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★抽象な `unit` が完全に具体化した

第 31・33・34 ブロックで `unit` を `freeYonedaEquiv` の言葉に翻訳した。
★本ブロックで、それを**生成元での値**として書き下す:

    (unit (free (yoneda V))).app (op V) (freeMk 𝟙)
      を第 24 ブロックの同型で送ると freeMk 𝟙 になる

★★すなわち **`unit` は生成元を生成元に送る**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeYonedaEquiv_apply` | ★`freeYonedaEquiv φ = φ.app _ (freeMk 𝟙)`(mathlib に無し) |
| `freeYonedaEquiv_id` | ★`freeYonedaEquiv 𝟙 = freeMk 𝟙` |
| `pullbackIso_hom_unit_gen` | ★★★★★**`unit` は生成元を生成元に送る** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

/-! ## ★`freeYonedaEquiv` は生成元での値である -/

section Generic

variable {C : Type u} [SmallCategory C] {R : Cᵒᵖ ⥤ RingCat.{u}}

/-- ★**`freeYonedaEquiv` は生成元での値である**(mathlib に無かったので書く)。

★機構は `freeYonedaEquiv_comp` に `𝟙` を入れるだけである。 -/
theorem freeYonedaEquiv_apply (M : PresheafOfModules.{u} R) (Z : C)
    (φ : (PresheafOfModules.free R).obj (yoneda.obj Z) ⟶ M) :
    PresheafOfModules.freeYonedaEquiv φ = φ.app (op Z) (ModuleCat.freeMk (𝟙 Z)) := by
  conv_lhs => rw [← Category.id_comp φ]
  rw [PresheafOfModules.freeYonedaEquiv_comp]
  congr 1

/-- ★恒等射の `freeYonedaEquiv` は生成元そのものである。 -/
theorem freeYonedaEquiv_id (Z : C) :
    PresheafOfModules.freeYonedaEquiv
        (𝟙 ((PresheafOfModules.free R).obj (yoneda.obj Z)))
      = ModuleCat.freeMk (𝟙 Z) := by
  rw [freeYonedaEquiv_apply]
  rfl

end Generic

/-! ## ★★★★★`unit` は生成元を生成元に送る -/

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★★★**`unit` は生成元を生成元に送る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★随伴関手定理で作られた**抽象な `unit`** が、これで完全に具体化した。

★機構は 3 段:
1. 第 34 ブロック: `freeYonedaEquiv (unit) = freeYonedaEquiv (同型の inv)`
2. `freeYonedaEquiv (inv ≫ hom) = hom.app (freeYonedaEquiv inv)`(合成則)
3. `inv ≫ hom = 𝟙` かつ `freeYonedaEquiv 𝟙 = freeMk 𝟙` -/
theorem pullbackIso_hom_unit_gen (V : Y.Opens) :
    ((pullbackFreeYonedaIso f V).hom).app (op ((Opens.map f.base).obj V))
        (PresheafOfModules.freeYonedaEquiv
          (M := (PresheafOfModules.pushforward (pullbackPhi f)).obj
            ((PresheafOfModules.pullback (pullbackPhi f)).obj (freeY V))) (X := V)
          ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (freeY V)))
      = ModuleCat.freeMk (𝟙 ((Opens.map f.base).obj V)) := by
  erw [freeYonedaEquiv_unit_free]
  erw [← PresheafOfModules.freeYonedaEquiv_comp, Iso.inv_hom_id, freeYonedaEquiv_id]

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——freeYonedaEquiv が生成元での値であること)",
    sectionId := "genell-def-1-1-i" }

def pullbackIso_hom_unit_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unit が生成元を生成元に送ること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
