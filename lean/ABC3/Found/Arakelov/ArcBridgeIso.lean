import ABC3.Found.Arakelov.ArcBridgeNat

/-!
# Arakelov (C3) 第 269 ブロック —— ★★★★★★★★**§9-297 の橋が完成した**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★`restrict F V.ι ≅ 𝒪_V`

`IsInvertibleSheaf`(`Interface` の可逆層の定義)が与える **`Over V` 上の自明化**から、
mathlib の語彙での自明化 **`restrict F V.ι ≅ unitModules V.toScheme`** を作った。

★★これが §9-297 で「元素レベルでは書けない、関手レベルで 5–15 ブロック」と
判定した橋である。★★★実際には**型付き恒等関数 5 本 + 3 法則 + 同型化**で済んだ。

| 段 | ブロック |
|---|---|
| 切断が `rfl` で一致 | 第 255 |
| 型付き恒等関数の橋 | ★第 265 |
| スカラー両立 | ★第 266 |
| `appIso` の自然性・射が `rfl` | ★第 267 |
| 自然性 | ★第 268 |
| ★★**同型化** | ★★★**本ブロック** |

## ★★同型化の摩擦

`SheafOfModules.forget` には `ReflectsIsomorphisms` の instance が**無い**。
★逆射を明示して `⟨⟨inv (bridgeHom …)⟩, Hom.ext …, Hom.ext …⟩` と書けば通る
——`(f ≫ g).val = f.val ≫ g.val` が `rfl` なので `Hom.ext` がそのまま効く。

★★`ConcreteCategory.isIso_iff_bijective` で切断ごとの全単射から `IsIso` を出し、
`NatTrans.isIso_iff_isIso_app` と `isIso_of_reflects_iso` で前層加群の同型に上げる
——第 55 ブロックと同じ 3 段である。

| 定義・定理 | 内容 |
|---|---|
| `bridgeInv_app` / `bridgeApp_bijective` | ★切断ごとに全単射 |
| `bridgeHom` / `bridgeSheafHom` | ★★★射 |
| `bridgeIso` | ★★★★★★★★**同型** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens) (F : X.Modules)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★もう一方の往復。 -/
theorem bridgeInv_app (W : (V.toScheme.Opens)ᵒᵖ)
    (y : ((unitModules V.toScheme).val.obj W : Type)) :
    bridgeApp F V e W (bridgeAppInv F V e W y) = y := by
  show (V.ι.appIso W.unop).hom.hom ((e.hom.app (op (overObj V W.unop))).hom
    ((e.inv.app (op (overObj V W.unop))).hom ((V.ι.appIso W.unop).inv.hom y))) = y
  have h1 : (e.hom.app (op (overObj V W.unop))).hom
      ((e.inv.app (op (overObj V W.unop))).hom ((V.ι.appIso W.unop).inv.hom y))
      = (V.ι.appIso W.unop).inv.hom y :=
    congrArg (fun (m : _ ⟶ _) => (m.app (op (overObj V W.unop))).hom
      ((V.ι.appIso W.unop).inv.hom y)) e.inv_hom_id
  rw [h1]
  exact congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m) y)
    (V.ι.appIso W.unop).inv_hom_id

/-- ★★切断ごとに全単射。 -/
theorem bridgeApp_bijective (W : (V.toScheme.Opens)ᵒᵖ) :
    Function.Bijective (bridgeApp F V e W) :=
  ⟨Function.LeftInverse.injective (bridgeApp_inv F V e W),
    Function.RightInverse.surjective (bridgeInv_app V F e W)⟩


/-- ★★★前層加群の射としての橋。 -/
noncomputable def bridgeHom :
    (Scheme.Modules.restrict F V.ι).val ⟶ (unitModules V.toScheme).val :=
  PresheafOfModules.homMk
    { app := fun W => AddCommGrpCat.ofHom
        (AddMonoidHom.mk' (bridgeApp F V e W) (fun x y => bridgeApp_add F V e W x y))
      naturality := by
        intro W W' i
        ext x
        exact bridgeApp_nat V F e i x }
    (fun W c s => bridgeApp_smul V F e W c s)

/-- ★★★★層加群の射。 -/
noncomputable def bridgeSheafHom :
    Scheme.Modules.restrict F V.ι ⟶ unitModules V.toScheme :=
  ⟨bridgeHom V F e⟩


/-- ★★★★★橋は同型である。 -/
noncomputable def bridgeIso :
    Scheme.Modules.restrict F V.ι ≅ unitModules V.toScheme := by
  haveI happ : ∀ W, IsIso (((PresheafOfModules.toPresheaf _).map (bridgeHom V F e)).app W) := by
    intro W
    rw [ConcreteCategory.isIso_iff_bijective]
    exact bridgeApp_bijective V F e W
  haveI : IsIso ((PresheafOfModules.toPresheaf _).map (bridgeHom V F e)) := by
    rw [NatTrans.isIso_iff_isIso_app]
    exact happ
  haveI : IsIso (bridgeHom V F e) :=
    isIso_of_reflects_iso _ (PresheafOfModules.toPresheaf _)
  haveI : IsIso (bridgeSheafHom V F e) :=
    ⟨⟨inv (bridgeHom V F e)⟩,
      SheafOfModules.Hom.ext (IsIso.hom_inv_id (bridgeHom V F e)),
      SheafOfModules.Hom.ext (IsIso.inv_hom_id (bridgeHom V F e))⟩
  exact asIso (bridgeSheafHom V F e)


/-! ## ★出典の紐付け(`.src`) -/

def bridgeIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——Over V の自明化から mathlib の自明化へ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
