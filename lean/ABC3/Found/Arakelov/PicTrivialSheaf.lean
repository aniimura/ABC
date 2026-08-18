import ABC3.Found.Arakelov.PicTrivialIso

/-!
# Arakelov (B1) 第 111 ブロック —— **層加群版の局所自明性の同型**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★層の仮定を埋める

第 110 ブロックの `trivialIsoOfSection` は 2 つの層の仮定を要求する。
★**層加群 `M : X.Modules` の制限**ならどちらも第 17 ブロックで出る:

| 仮定 | 材料 |
|---|---|
| `(restrict V).obj M.val` が層 | ★`isSheaf_restrict X V _ M.isSheaf` |
| `𝟙_` が層 | ★`isSheaf_restrict` + `sheafCompose` |

## ★★本ブロックで取れるもの

| 定理・定義 | 内容 |
|---|---|
| `isSheaf_restrictModules` | ★制限した層加群は層 |
| `isSheaf_unitOn` | ★単位は層 |
| `trivialIsoOfSectionSheaf` | ★★★★★★**層加群版の同型** |

## ★★★残り

これを `M := tilde M` に当て、切断を第 76 の生成元、
被覆を基本開集合(第 100)にすれば `IsLocallyTrivial` が出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens)

/-- ★**制限した層加群は層である**。 -/
theorem isSheaf_restrictModules (M : X.Modules) :
    Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      ((restrictPresheafFunctor X V).obj M.val).presheaf :=
  isSheaf_restrict X V _ M.isSheaf

/-- ★**制限した単位は層である**。 -/
theorem isSheaf_unitOn :
    Presheaf.IsSheaf ((Opens.grothendieckTopology X).over V)
      (𝟙_ (PresheafModulesOn X V)).presheaf :=
  isSheaf_restrict X V _
    ((sheafCompose (Opens.grothendieckTopology X)
      (forget₂ RingCat.{u} AddCommGrpCat.{u})).obj X.ringCatSheaf).property

/-- ★★★★★★**層加群版の局所自明性の同型**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★切断があって被覆の上で「その倍」が全単射なら、制限は単位と同型である。 -/
noncomputable def trivialIsoOfSectionSheaf (M : X.Modules)
    (s : ((restrictPresheafFunctor X V).obj M.val).obj (op (Over.mk (𝟙 V))))
    (h : ∀ W : Over V, ∃ S : Sieve W, S ∈ ((Opens.grothendieckTopology X).over V) W ∧
      ∀ ⦃Z : Over V⦄ (i : Z ⟶ W), S.arrows i →
        Function.Bijective (fun c : ((((Over.forget V).op ⋙ X.presheaf)
            ⋙ forget₂ CommRingCat.{u} RingCat.{u}).obj (op Z) : Type u) =>
          c • restrictSec V ((restrictPresheafFunctor X V).obj M.val) s Z)) :
    𝟙_ (PresheafModulesOn X V) ≅ (restrictPresheafFunctor X V).obj M.val :=
  trivialIsoOfSection V _ (isSheaf_restrictModules V M) (isSheaf_unitOn V) s h

/-! ## ★出典の紐付け(`.src`) -/

def trivialIsoOfSectionSheaf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——層加群版の局所自明性の同型)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
