import ABC3.Found.Arakelov.PicAppLEApply

/-!
# Arakelov (B2) 第 215 ブロック —— ★★★★**`appTop` を「制限 ≫ appIso」に分解**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-240 で特定した B2 の最後の 1 段

    g.appTop = (制限 ⊤ → g ''ᵁ ⊤) ≫ (g.appIso ⊤).hom     （g は開埋め込み)

★これで第 214(`appTop` の可換正方形)が、生成元の場合で要る
**`Γ(Y, B')` から始まる形**に繋がる。

## ★★機構

`appIso_hom'`(mathlib)が `(g.appIso U).hom = g.appLE (g ''ᵁ U) U _` を与える。
`appLE` を開いて `g.app` の**自然性**で `g.app ⊤` に寄せると、
残るのは `X.presheaf` の制限の合成が恒等であることだけ
——`Functor.map_id` で閉じる。

★★★摩擦は「`Opens` の射は Prop なので合成が自動で潰れない」点であった。
`congr 1` で `𝟙 = map (…)` の形に落としてから `map_id` を当てる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appTop_eq_res_appIso` | ★★★★**`appTop` の分解** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (g : X ⟶ Y) [IsOpenImmersion g]

/-- ★★★★**`appTop` は「制限 ≫ appIso」に分解する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any
★★★これが B2 の最後の 1 段である。 -/

theorem appTop_eq_res_appIso :
    g.appTop = Y.presheaf.map (homOfLE (le_top : g ''ᵁ (⊤ : X.Opens) ≤ ⊤)).op
      ≫ (g.appIso ⊤).hom := by
  rw [Scheme.Hom.appIso_hom' g ⊤, Scheme.Hom.appLE, ← Category.assoc,
    Scheme.Hom.naturality g (homOfLE (le_top : g ''ᵁ (⊤ : X.Opens) ≤ ⊤)).op]
  show Scheme.Hom.app g ⊤ = _
  rw [Category.assoc, ← Functor.map_comp, ← op_comp]
  exact (Category.comp_id _).symm.trans (by
    congr 1
    exact (CategoryTheory.Functor.map_id _ _).symm)

/-! ## ★出典の紐付け(`.src`) -/

def appTop_eq_res_appIso.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——appTop を制限と appIso に分解)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
