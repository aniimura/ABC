import ABC3.Found.Arakelov.PicInfGen

/-!
# Arakelov (B1) 第 38 ブロック —— **`unit` の値と同型の `inv` の値は一致する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★最終段の成分等式

最終段の両辺は、成分ごとに

    左辺: unit_P.app (op Z) (freeMk i)          (随伴の unit)
    右辺: inv.app (op f⁻¹Z) (freeMk (f⁻¹i))     (第 24 ブロックの同型の inv)

となる。★★**この 2 つが等しい**というのが本ブロックである。

## ★★★理由は 2 行

1. どちらも「生成元での値の**制限**」である(第 36 ブロック)
2. 生成元での値は等しい(第 34 ブロック)

★押し出しの制限 `(f_* M).map i.op` は `M.map (f⁻¹i).op` そのものなので、
制限の仕方まで一致する。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `unit_app_eq_inv_app` | ★★★★★**`unit` の値 = 同型の `inv` の値** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★★★**`unit` の値と、第 24 ブロックの同型の `inv` の値は一致する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが最終段(`δ` が生成元の上で同型であること)の**成分等式**である。

★機構は 2 段:
1. どちらも「生成元での値の制限」(第 36 ブロック `free_app_freeMk`)
2. 生成元での値は等しい(第 34 ブロック `freeYonedaEquiv_unit_free`)

★押し出しの制限 `(f_* M).map i.op` は `M.map (f⁻¹i).op` **そのもの**なので、
制限の仕方まで一致する。 -/
theorem unit_app_eq_inv_app (V : Y.Opens) {Z : Y.Opens} (i : Z ⟶ V) :
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app
        (freeY V)).app (op Z) (ModuleCat.freeMk i)
      = ((pullbackFreeYonedaIso f V).inv).app (op ((Opens.map f.base).obj Z))
        (ModuleCat.freeMk ((Opens.map f.base).map i)) := by
  have h1 := free_app_freeMk
    (M := (PresheafOfModules.pushforward (pullbackPhi f)).obj
      ((PresheafOfModules.pullback (pullbackPhi f)).obj (freeY V))) i
    ((PresheafOfModules.pullbackPushforwardAdjunction (pullbackPhi f)).unit.app (freeY V))
  have h2 := free_app_freeMk
    (M := (PresheafOfModules.pullback (pullbackPhi f)).obj (freeY V))
    ((Opens.map f.base).map i) ((pullbackFreeYonedaIso f V).inv)
  erw [h1, h2, freeYonedaEquiv_unit_free]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def unit_app_eq_inv_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——unit の値と同型の inv の値が一致すること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
