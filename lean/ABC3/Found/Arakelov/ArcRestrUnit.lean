import ABC3.Found.Arakelov.ArcGammaEval

/-!
# Arakelov (C3) 第 280 ブロック —— **制限の同型と単位の関係**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-320 の穴に残っていた**唯一の非自明な等式**

第 278–279 ブロックで、`arcEval` の整合性は次の 1 点に還元された:

    (restrictFunctorIsoPullback j).inv ∘ (j の随伴の単位) = 「制限写像」

★★これだけが `rfl` では出ない。

## ★★★鍵——`restrictFunctorIsoPullback` は `leftAdjointUniq` である

mathlib の定義は

    restrictFunctorIsoPullback j = (restrictAdjunction j).leftAdjointUniq
                                     (pullbackPushforwardAdjunction j)

であり、`Adjunction.unit_leftAdjointUniq_hom_app` が

    (restrictAdjunction j).unit ≫ (pushforward j).map (iso.hom) = (pullback の随伴).unit

を与える。★★★したがって `iso.inv` を掛けると `iso.hom ≫ iso.inv = 𝟙` で消え、
**残るのは `restrictAdjunction` の単位**——それは `restrictAdjunction_unit_app_app` により
**制限写像そのもの**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pushforwardMapTop` | ★押し出しは `⊤` で素通し(一般の射に対して) |
| `restrictAdjUnit_apply` | ★★`restrictAdjunction` の単位は制限写像(`rfl`) |
| `restrictIso_unit_apply` | ★★★★★**穴を埋める等式** |

## ★★★測定の記録

§9-320 で「`arcEval` の整合性は **3–6 ブロック**」と見積もった。
★★**実測 3 ブロック**(278・279・280)——見積もりが当たったのは久しぶりである。
★理由は明確で、**mathlib に鍵となる補題が既に在ることを先に確かめてから見積もった**からである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

/-- ★**押し出しは `⊤` の切断では素通しである**(一般の射に対して)。 -/
theorem pushforwardMapTop {A B : Scheme.{0}} (g : A ⟶ B) {M N : A.Modules} (φ : M ⟶ N)
    (v : (((Scheme.Modules.pushforward g).obj M).val.obj (op ⊤) : Type)) :
    (((Scheme.Modules.pushforward g).map φ).val.app (op ⊤)).hom v = (φ.val.app (op ⊤)).hom v := rfl

variable {U X : Scheme.{0}} (j : U ⟶ X) [IsOpenImmersion j] (L : X.Modules)

/-- ★★**`restrictAdjunction` の単位は、`⊤` の切断では制限写像そのものである**(`rfl`)。 -/
theorem restrictAdjUnit_apply (s : (L.val.obj (op ⊤) : Type)) :
    (((Scheme.Modules.restrictAdjunction j).unit.app L).val.app (op ⊤)).hom s
      = (L.val.map (homOfLE (le_top : (j ''ᵁ (⊤ : U.Opens)) ≤ ⊤)).op).hom s := rfl

/-- ★★★★★**§9-320 の穴を埋める等式**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`pullback` の随伴の単位を `restrictFunctorIsoPullback` で `restrict` 側へ移すと、
**ただの制限写像**になる。 -/
theorem restrictIso_unit_apply (s : (L.val.obj (op ⊤) : Type)) :
    ((((Scheme.Modules.restrictFunctorIsoPullback j).app L).inv.val.app (op ⊤)).hom
        ((((Scheme.Modules.pullbackPushforwardAdjunction j).unit.app L).val.app (op ⊤)).hom s))
      = (L.val.map (homOfLE (le_top : (j ''ᵁ (⊤ : U.Opens)) ≤ ⊤)).op).hom s := by
  rw [← Adjunction.unit_leftAdjointUniq_hom_app
    (Scheme.Modules.restrictAdjunction j) (Scheme.Modules.pullbackPushforwardAdjunction j) L]
  have hid : ∀ w : ((Scheme.Modules.restrict L j).val.obj (op ⊤) : Type),
      ((((Scheme.Modules.restrictFunctorIsoPullback j).app L).inv.val.app (op ⊤)).hom
        ((((Scheme.Modules.restrictFunctorIsoPullback j).app L).hom.val.app (op ⊤)).hom w)) = w := by
    intro w
    have h2 : ((((Scheme.Modules.restrictFunctorIsoPullback j).app L).hom
        ≫ ((Scheme.Modules.restrictFunctorIsoPullback j).app L).inv).val.app (op ⊤)).hom w = w := by
      rw [Iso.hom_inv_id]
      rfl
    exact h2
  exact hid ((((Scheme.Modules.restrictAdjunction j).unit.app L).val.app (op ⊤)).hom s)

/-! ## ★出典の紐付け(`.src`) -/

def restrictIso_unit_apply.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——制限の同型と随伴の単位の関係)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
