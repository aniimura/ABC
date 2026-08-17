import ABC3.Found.Arakelov.PicGenVal

/-!
# Arakelov (B1) 第 36 ブロック —— **生成元の像は制限で決まる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★最終段の核

第 35 ブロックまでで、`unit` も `μ` も生成元の上で書き下せた。
★残るのは「両辺が一致すること」であり、その理由は **1 行**である:

> **両辺とも `freeYonedaEquiv (unit)` の `U₀` への制限である。**

★★本ブロックはその「制限である」を汎用の補題にする:

    φ.app (op Z) (freeMk i) = M.map i.op (freeYonedaEquiv φ)

★★★すなわち **`free (yoneda V)` からの射は、生成元での 1 つの値で決まり、
他の切断での値はその制限である**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeYoneda_map_gen` | ★生成元の制限は `freeMk i` |
| `free_app_freeMk` | ★★★★**`φ.app (freeMk i) = M.map i.op (freeYonedaEquiv φ)`** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite Limits

variable {C : Type u} [SmallCategory C] {R : Cᵒᵖ ⥤ RingCat.{u}}

/-! ## ★生成元の制限 -/

/-- ★**生成元 `freeMk 𝟙` を `i` で制限すると `freeMk i` になる**。 -/
theorem freeYoneda_map_gen {V Z : C} (i : Z ⟶ V) :
    ((PresheafOfModules.free R).obj (yoneda.obj V)).map i.op
        (ModuleCat.freeMk (𝟙 V))
      = ModuleCat.freeMk i := by
  show ModuleCat.freeDesc _ _ = _
  rw [ModuleCat.freeDesc_apply]
  show ModuleCat.freeMk (i ≫ 𝟙 V) = _
  rw [Category.comp_id]

/-! ## ★★★★生成元での値がすべてを決める -/

/-- ★★★★★**`free (yoneda V)` からの射は、生成元での 1 つの値で決まる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★他の切断での値は、その値の**制限**である。
★これが最終段(`δ` が生成元の上で同型であること)の核である——
`unit` の値も、第 24 ブロックの同型の `inv` の値も、
どちらも「生成元での値の制限」だから一致する。 -/
theorem free_app_freeMk {M : PresheafOfModules.{u} R} {V Z : C} (i : Z ⟶ V)
    (φ : (PresheafOfModules.free R).obj (yoneda.obj V) ⟶ M) :
    φ.app (op Z) (ModuleCat.freeMk i)
      = M.map i.op (PresheafOfModules.freeYonedaEquiv φ) := by
  rw [freeYonedaEquiv_apply, ← freeYoneda_map_gen i]
  exact NatTrans.naturality_apply
    ((PresheafOfModules.toPresheaf R).map φ) i.op (ModuleCat.freeMk (𝟙 V))

/-! ## ★出典の紐付け(`.src`) -/

def free_app_freeMk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——生成元での値が射を決めること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
