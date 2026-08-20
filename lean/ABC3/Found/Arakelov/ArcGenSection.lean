import ABC3.Found.Arakelov.ArcGlueZero

/-!
# Arakelov (C3) 第 264 ブロック —— ★★★**自明化から生成切断を取る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★非消滅の 2 つの半分のうち、片方

第 261 で仮定として受けた「生成切断の非消滅」は 2 つに割れる:

| 半分 | 内容 | 状態 |
|---|---|---|
| (B) `g` が `Γ(V,F)` を生成する | ★自明化 `e` から直ちに出る | ★★**本ブロック** |
| (A) `arcEvalOnTop` が零でない | ★`p^*` との比較(§9-297 の壁) | ★残り |

★★(B) は `g := e.inv(1)` と置けば `s = e.hom(s) · g` であり、
**`e.inv` の線形性と `c · 1 = c` だけ**で出る。

## ★摩擦 —— 単位対象の `1` は書けない

`(𝟙_ (PresheafModulesOn X V)).obj (op (Over.mk (𝟙 V)))` には `OfNat 1` が無い
(`Γ(X,V)` と定義的に等しいのに)。
★逃げ道は**結果の型を明示した `def`** で `1` を渡すこと(`unitOne`)——
期待型からの elaboration は定義的等しさを通す。

★★係数の側も同じで、`c • unitOne` の `c` は `forget₂` 側の綴りが要る。
★★★第 60 ブロックの `rVal`(係数の橋)がそのまま使えた——**在庫が効いた**。

| 定義・定理 | 内容 |
|---|---|
| `unitOne` | ★単位対象の切断としての `1` |
| `smul_unitOne` | ★★`c • 1 = c` |
| `genSection` | ★★自明化から取る生成切断 |
| `genSection_spans` | ★★★**生成すること** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} (F : X.Modules) (V : X.Opens)
  (e : (restrictPresheafFunctor X V).obj F.val ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★自明化から生成切断を取る。 -/
noncomputable def genSection : (F.val.obj (op V) : Type) :=
  (e.inv.app (op (Over.mk (𝟙 V)))).hom (unitOne V)

/-- ★★★生成切断はその名のとおり生成する。 -/
theorem genSection_spans (s : (F.val.obj (op V) : Type)) :
    ∃ c : ((X.presheaf.obj (op V)) : Type), s = c • genSection F V e := by
  refine ⟨(e.hom.app (op (Over.mk (𝟙 V)))).hom s, ?_⟩
  have hround : (e.inv.app (op (Over.mk (𝟙 V)))).hom
      ((e.hom.app (op (Over.mk (𝟙 V)))).hom s) = s :=
    congrArg (fun (m : _ ⟶ _) => (m.app (op (Over.mk (𝟙 V)))).hom s) e.hom_inv_id
  have hlin : (e.inv.app (op (Over.mk (𝟙 V)))).hom
      (((e.hom.app (op (Over.mk (𝟙 V)))).hom s) • unitOne V)
      = ((e.hom.app (op (Over.mk (𝟙 V)))).hom s) • genSection F V e :=
    map_smul _ _ _
  have h1 : (e.inv.app (op (Over.mk (𝟙 V)))).hom
      (((e.hom.app (op (Over.mk (𝟙 V)))).hom s) • unitOne V)
      = (e.inv.app (op (Over.mk (𝟙 V)))).hom
        ((e.hom.app (op (Over.mk (𝟙 V)))).hom s) :=
    congrArg (fun z => (e.inv.app (op (Over.mk (𝟙 V)))).hom z)
      (smul_unitOne V ((e.hom.app (op (Over.mk (𝟙 V)))).hom s))
  exact (hlin.symm.trans (h1.trans hround)).symm


/-! ## ★出典の紐付け(`.src`) -/

def genSection_spans.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——自明化から生成切断を取る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
