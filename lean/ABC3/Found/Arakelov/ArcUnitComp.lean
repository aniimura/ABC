import ABC3.Found.Arakelov.ArcXNormVia
import Mathlib.CategoryTheory.Adjunction.Mates

/-!
# Arakelov (C3) 第 278 ブロック —— **随伴の単位は合成と両立する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-320 の「残り(1)」を開ける鍵

C3 の `normSection_continuous` に残っていた唯一の穴は

    (arcFiberFactor j L p).hom (arcEval (p ≫ j) L s)
      = arcEval p (restrict L j) (restrictSection j L s)

であった。★★両辺は `p ≫ j` の随伴の単位と、`j` の単位・`p` の単位の**合成**であり、
`rfl` では出ない。

## ★★★mathlib に鍵が在った

| 補題 | 場所 | 内容 |
|---|---|---|
| `Adjunction.comp_unit_app` | `Adjunction/Basic.lean` | 合成随伴の単位の分解 |
| `unit_conjugateEquiv` | `Adjunction/Mates.lean` | ★★共役と単位の両立 |
| `Scheme.Modules.conjugateEquiv_pullbackComp_inv` | `Modules/Sheaf.lean` | ★★★`pullbackComp` の共役は `pushforwardComp` |

★★★3 つを繋ぐと、**射のレベルでの両立**が一発で出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `unit_conj_compat` | ★★共役版(mathlib をそのまま噛ませた形) |
| ★`Adjunction.comp_unit_app` | ★mathlib 側に在るので**証明の中で** `rw` する |

★これを `⊤` の切断に評価すれば §9-320 の穴が閉じる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable {U X : Scheme.{0}} (p : Spec (CommRingCat.of ℂ) ⟶ U) (j : U ⟶ X) (L : X.Modules)

/-- ★★**共役版**——`pullbackComp` の共役が `pushforwardComp` であることから直ちに出る。 -/
theorem unit_conj_compat :
    ((Scheme.Modules.pullbackPushforwardAdjunction j).comp
        (Scheme.Modules.pullbackPushforwardAdjunction p)).unit.app L
      ≫ (Scheme.Modules.pushforwardComp p j).hom.app
        ((Scheme.Modules.pullback j ⋙ Scheme.Modules.pullback p).obj L)
    = (Scheme.Modules.pullbackPushforwardAdjunction (p ≫ j)).unit.app L
      ≫ (Scheme.Modules.pushforward (p ≫ j)).map
        ((Scheme.Modules.pullbackComp p j).inv.app L) := by
  have h := unit_conjugateEquiv
    ((Scheme.Modules.pullbackPushforwardAdjunction j).comp
      (Scheme.Modules.pullbackPushforwardAdjunction p))
    (Scheme.Modules.pullbackPushforwardAdjunction (p ≫ j))
    (Scheme.Modules.pullbackComp p j).inv L
  rw [Scheme.Modules.conjugateEquiv_pullbackComp_inv] at h
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def unit_conj_compat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——随伴の単位が合成と両立すること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
