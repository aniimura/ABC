import ABC3.Found.Arakelov.PicModCoverIso

/-!
# Arakelov (B1) 第 110 ブロック —— **切断から局所自明性の同型を作る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★器具が全部繋がった

    切断 s があって、被覆の上で「s の倍」が全単射  ⟹  𝟙_ ≅ P

★★これが `IsLocallyTrivial` の同型そのものである。

## ★★機構 —— 積み上げた器具の合成

| 段 | 出典 |
|---|---|
| `unitHomOfSection` の `app` が全単射 | ★第 108 ブロック |
| 被覆で全単射なら同型(前層加群) | ★第 109 ブロック |

★**一発で通った**——第 99〜109 の 11 ブロックが噛み合った。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `overTermInst` | ★`Unique` を instance にする(経路を揃えるため) |
| `isIso_unitHomOfSection` | ★★★★★**切断から同型** |
| `trivialIsoOfSection` | ★★★★★★**`𝟙_ ≅ P`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★**`Unique` を instance にする**——`default` の経路を揃えるためである
([[inhabited-two-paths]])。 -/
noncomputable instance overTermInst (V : X.Opens) (Z : Over V) :
    Unique ((yoneda.obj (Over.mk (𝟙 V))).obj (op Z)) := overTerminalUnique V Z

/-- ★★★★★**切断から同型が出る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★被覆の上で「`s` の倍」が全単射なら、`unitHomOfSection s` は同型である。 -/
theorem isIso_unitHomOfSection (V : X.Opens) (P : PresheafModulesOn X V)
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) P.presheaf)
    (h1 : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      (𝟙_ (PresheafModulesOn X V)).presheaf)
    (s : P.obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, ∃ S : Sieve W, S ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), S.arrows i →
        Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op Z) : Type u) =>
          c • restrictSec V P s Z)) :
    IsIso (unitHomOfSection V P s) := by
  refine isIso_of_bijective_on_cover_mod h1 hP _ (fun W => ?_)
  obtain ⟨S, hS, hb⟩ := h W
  exact ⟨S, hS, fun {Z} i hi => bijective_unitHomOfSection_app V P s Z (hb i hi)⟩

/-- ★★★★★★**`𝟙_ ≅ P`**——局所自明性の同型そのもの。 -/
noncomputable def trivialIsoOfSection (V : X.Opens) (P : PresheafModulesOn X V)
    (hP : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V) P.presheaf)
    (h1 : Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      (𝟙_ (PresheafModulesOn X V)).presheaf)
    (s : P.obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, ∃ S : Sieve W, S ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), S.arrows i →
        Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op Z) : Type u) =>
          c • restrictSec V P s Z)) :
    𝟙_ (PresheafModulesOn X V) ≅ P :=
  haveI := isIso_unitHomOfSection V P hP h1 s h
  asIso (unitHomOfSection V P s)

/-! ## ★出典の紐付け(`.src`) -/

def trivialIsoOfSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断から局所自明性の同型を作る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
