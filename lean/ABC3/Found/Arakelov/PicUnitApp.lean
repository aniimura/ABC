import ABC3.Found.Arakelov.PicGenMul

/-!
# Arakelov (B1) 第 104 ブロック —— **`unitHomOfSection` の `app` の値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★12 手詰まってから抜けた

§9-110〜§9-112 で **12 手**試して止めた所である。

★★決め手は **`freeYonedaEquiv.symm s = freeObjDesc (yonedaEquiv.symm s)` が `rfl`**
であること(2026-08-24 実測)。★★★そこまで `show` で降ろせば
`freeObjDesc_app`(`@[simps]`)と `ModuleCat.freeDesc_apply` が当たる。

★**`erw` が要る**——第 25 ブロックで環インスタンスの二経路を抜けたのと**同じ手**である。

## ★★詰まっていた理由(記録)

| 誤った道 | なぜ駄目か |
|---|---|
| `simp [freeYonedaEquiv, freeHomEquiv, yonedaEquiv]` | 展開が深すぎて閉じない |
| `NatTrans.naturality_apply` | ★**`PresheafOfModules.Hom` は `NatTrans` ではない** |
| `mulHom` を新しく作って自然性を示す | `ModuleCat.hom_ofHom` が剥がれない |

★★★**「新しく作る」より「既にあるものの値を計算する」方が短かった。**

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `freeYonedaEquiv_symm_eq_desc` | ★`freeYonedaEquiv.symm` は `freeObjDesc`(`rfl`) |
| `freeYonedaEquiv_symm_app_freeMk` | ★★★★**`app` の値は `P.map h.op s`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens) (P : PresheafModulesOn X V)

/-- ★**`freeYonedaEquiv.symm` は `freeObjDesc` である**(`rfl`)。 -/
theorem freeYonedaEquiv_symm_eq_desc (s : P.obj (op (Over.mk (𝟙 V)))) :
    PresheafOfModules.freeYonedaEquiv.symm s
      = PresheafOfModules.freeObjDesc (yonedaEquiv.symm s) := rfl

/-- ★★★★**`app` の値は `P.map h.op s` である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが局所自明性の同型を「乗法」として読むための鍵である。 -/
theorem freeYonedaEquiv_symm_app_freeMk (s : P.obj (op (Over.mk (𝟙 V))))
    (W : Over V) (h : W ⟶ Over.mk (𝟙 V)) :
    (PresheafOfModules.freeYonedaEquiv.symm s).app (op W) (ModuleCat.freeMk h)
      = P.map h.op s := by
  show (PresheafOfModules.freeObjDesc (yonedaEquiv.symm s)).app (op W) (ModuleCat.freeMk h) = _
  rw [PresheafOfModules.freeObjDesc_app]
  erw [ModuleCat.freeDesc_apply]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def freeYonedaEquiv_symm_app_freeMk.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unitHomOfSection の app の値)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
