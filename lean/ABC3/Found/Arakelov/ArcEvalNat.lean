import ABC3.Found.Arakelov.ArcContCriterion

/-!
# Arakelov (C3) 第 250 ブロック —— **評価は層の射と可換である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★連続なノルムを作るための橋

計量を「双対の汎関数の絶対値の max」で作るとき、連続性は

    φ̄(arcEval p L s) = arcEval p 𝒪 (φ(s))     (φ : L ⟶ 𝒪 は層の射)

に帰着する。★右辺は**正則関数の点での値**なので `continuous_evalAffine` が効く。
★★本ブロックはこの可換性(引き戻しの単位射の自然性)を Lean に載せたものである。

## ★★★摩擦 3 つ —— どれも「構文照合 vs 定義的等しさ」

| # | 症状 | 逃げ道 |
|---|---|---|
| 1 | `(f ≫ g).val.app W` が割れない(`simp` が発火しない) | ★**明示束縛子**の `rfl` 補題 `hsplit` にする |
| 2 | `congrArg` の結果が **β 簡約されていない** | ★`have` の型を**書き下す** |
| 3 | `Eq.trans` が `(F ⋙ G).map φ` と `G.map (F.map φ)` を繋げない | ★中間項の**綴りを揃える** |

★★★#1 は「暗黙束縛子だとインスタンス透過性で落ちる」——**明示にすると `rfl` で通る**。
★★#3 は `Eq.trans` の単一化が **reducible 透過性**しか使わないためである
——`exact` は `default` 透過性なので通るが、`trans` は通らない。

| 定理 | 内容 |
|---|---|
| `hsplit` | ★合成の切断レベルの分解(`rfl`) |
| `arcEval_naturality` | ★★★★**評価は層の射と可換** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {X : Scheme.{0}}
theorem hsplit (A B C : X.Modules) (f : A ⟶ B) (g : B ⟶ C)
    (x : (A.val.obj (op ⊤) : Type)) :
    ((f ≫ g).val.app (op ⊤)).hom x = (g.val.app (op ⊤)).hom ((f.val.app (op ⊤)).hom x) :=
  rfl

theorem arcEval_naturality (p : Spec (CommRingCat.of ℂ) ⟶ X) {L M : X.Modules} (φ : L ⟶ M)
    (s : (L.val.obj (op ⊤) : Type)) :
    arcEval p M ((φ.val.app (op ⊤)).hom s)
      = ((((Scheme.Modules.pushforward p).map
          ((Scheme.Modules.pullback p).map φ)).val.app (op ⊤))).hom (arcEval p L s) := by
  have h : (((φ ≫ (Scheme.Modules.pullbackPushforwardAdjunction p).unit.app M).val.app
          (op ⊤))).hom s
      = (((((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app L)
        ≫ (Scheme.Modules.pullback p ⋙ Scheme.Modules.pushforward p).map
          φ).val.app (op ⊤))).hom s :=
    congrArg (fun (m : L ⟶ (Scheme.Modules.pushforward p).obj
        ((Scheme.Modules.pullback p).obj M)) => (m.val.app (op ⊤)).hom s)
      ((Scheme.Modules.pullbackPushforwardAdjunction p).unit.naturality φ)
  have h1 : (((((Scheme.Modules.pullbackPushforwardAdjunction p).unit.app L)
        ≫ (Scheme.Modules.pullback p ⋙ Scheme.Modules.pushforward p).map
          φ).val.app (op ⊤))).hom s
      = ((((Scheme.Modules.pushforward p).map
          ((Scheme.Modules.pullback p).map φ)).val.app (op ⊤))).hom (arcEval p L s) :=
    hsplit _ _ _ _ _ s
  have h2 : (((φ ≫ (Scheme.Modules.pullbackPushforwardAdjunction p).unit.app M).val.app
        (op ⊤))).hom s
      = arcEval p M ((φ.val.app (op ⊤)).hom s) :=
    hsplit _ _ _ _ _ s
  exact h2.symm.trans (h.trans h1)


/-! ## ★出典の紐付け(`.src`) -/

def arcEval_naturality.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——評価は層の射と可換)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
