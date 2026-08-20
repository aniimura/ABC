import ABC3.Found.Arakelov.ArcRestrNorm

/-!
# Arakelov (C3) 第 255 ブロック —— **`Over V` と `V.Opens` の橋(切断の段)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★残っていた 1 段の正体は「添字圏の違い」だけだった

`IsInvertibleSheaf`(`Interface` の可逆層の定義)が与える自明化は
**`Over V` 上の前層加群**の同型であり、第 254 が要求するのは
**`restrict F V.ι ≅ 𝒪_V`**(mathlib の語彙)である。

★★★実測すると、切断はどちらも `F.val.obj (op (V.ι ''ᵁ W))` で**`rfl` で一致する**:

    Γ(restrict F V.ι, W)                     = Γ(F, V.ι ''ᵁ W)      ★rfl
    ((restrictPresheafFunctor X V).obj F.val).obj (op (overObj W))  = 同上  ★rfl

★★したがって橋は**添字の付け替えだけ**であり、`V.ι.appIso` で係数環を合わせれば
切断ごとの同型が書ける。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `image_le` / `overObj` | ★`W ⊆ V` から `Over V` の対象へ |
| `restrict_sections` / `restrictPre_sections` | ★★**切断は `rfl` で一致** |
| `bridgeApp` / `bridgeAppInv` | ★★★切断ごとの写像(両向き) |
| `bridgeApp_inv` | ★★★★往復が恒等 |

★残るのは**自然性とスカラー両立**であり、それが済めば
`restrict F V.ι ≅ 𝒪_V` が出て第 254 に繋がる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (F : X.Modules) (V : X.Opens)

/-- ★像は `V` に含まれる。 -/
theorem image_le (W : V.toScheme.Opens) : V.ι ''ᵁ W ≤ V :=
  Scheme.Opens.ι_image_le V W

/-- ★対応する `Over V` の対象。 -/
def overObj (W : V.toScheme.Opens) : Over V :=
  Over.mk (homOfLE (image_le V W))

/-- ★★切断の橋: `Γ(restrict F ι, W) = Γ(F, ι ''ᵁ W)` は `rfl`。 -/
theorem restrict_sections (W : V.toScheme.Opens) :
    (((Scheme.Modules.restrict F V.ι).val.obj (op W)) : Type)
      = ((F.val.obj (op (V.ι ''ᵁ W))) : Type) :=
  rfl

/-- ★★制限した前層加群の切断も同じ。 -/
theorem restrictPre_sections (W : V.toScheme.Opens) :
    ((((restrictPresheafFunctor X V).obj F.val).obj (op (overObj V W))) : Type)
      = ((F.val.obj (op (V.ι ''ᵁ W))) : Type) :=
  rfl

variable (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★★★切断ごとの写像。 -/
noncomputable def bridgeApp (W : (V.toScheme.Opens)ᵒᵖ) :
    ((Scheme.Modules.restrict F V.ι).val.obj W : Type)
      → ((unitModules V.toScheme).val.obj W : Type) :=
  fun x => (V.ι.appIso W.unop).hom.hom
    ((e.hom.app (op (overObj V W.unop))).hom x)


/-- ★★★逆向きの写像。 -/
noncomputable def bridgeAppInv (W : (V.toScheme.Opens)ᵒᵖ) :
    ((unitModules V.toScheme).val.obj W : Type)
      → ((Scheme.Modules.restrict F V.ι).val.obj W : Type) :=
  fun y => (e.inv.app (op (overObj V W.unop))).hom ((V.ι.appIso W.unop).inv.hom y)

theorem bridgeApp_inv (W : (V.toScheme.Opens)ᵒᵖ)
    (x : ((Scheme.Modules.restrict F V.ι).val.obj W : Type)) :
    bridgeAppInv F V e W (bridgeApp F V e W x) = x := by
  show (e.inv.app (op (overObj V W.unop))).hom ((V.ι.appIso W.unop).inv.hom
      ((V.ι.appIso W.unop).hom.hom ((e.hom.app (op (overObj V W.unop))).hom x))) = x
  have h1 : (V.ι.appIso W.unop).inv.hom ((V.ι.appIso W.unop).hom.hom
      ((e.hom.app (op (overObj V W.unop))).hom x))
      = (e.hom.app (op (overObj V W.unop))).hom x :=
    congrArg (fun (m : _ ⟶ _) => (CommRingCat.Hom.hom m)
      ((e.hom.app (op (overObj V W.unop))).hom x)) (V.ι.appIso W.unop).hom_inv_id
  rw [h1]
  exact congrArg (fun (m : _ ⟶ _) => (m.app (op (overObj V W.unop))).hom x) e.hom_inv_id



/-- ★★加法性——`bridgeApp` は加法的である。 -/
theorem bridgeApp_add (W : (V.toScheme.Opens)ᵒᵖ)
    (x y : ((Scheme.Modules.restrict F V.ι).val.obj W : Type)) :
    bridgeApp F V e W (x + y) = bridgeApp F V e W x + bridgeApp F V e W y :=
  (congrArg (fun z => (V.ι.appIso W.unop).hom.hom z)
    (map_add (e.hom.app (op (overObj V W.unop))).hom x y)).trans (map_add _ _ _)

/-! ## ★出典の紐付け(`.src`) -/

def bridgeApp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——Over V と V.Opens の橋)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
