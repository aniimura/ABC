import ABC3.Found.Arakelov.PicResSquare

/-!
# Arakelov (B1) 第 121 ブロック —— **2 つの前層は切断も制限も同じ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性へ渡す橋

第 120 ブロックの可換図式は **`modulesSpecToSheaf.obj (tilde M)`**(`R` 加群の前層)
について述べてある。★一方 `IsLocallyTrivial` は **`(tilde M).val`**
(`𝒪` 加群の前層)について述べる。

★★**両者は切断も制限射も `rfl` で一致する**(2026-08-24 実測)——
違うのは**加群構造だけ**である。

★★★したがって元のレベルでは自由に行き来できる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tildeSection_eq` | ★切断の台は同じ(`rfl`) |
| `tildeRes_eq` | ★★制限射も同じ(`rfl`) |

## ★★★これで第 120 の可換図式が `IsLocallyTrivial` の側で使える

    restrictSec(= (tilde M).val の制限)= modulesSpecToSheaf 側の制限
      = 局所化された生成元(第 120)
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (U : (Spec R).Opens)

/-- ★**切断の台は同じである**(`rfl`)。 -/
theorem tildeSection_eq : ((tilde M).val.obj (op U) : Type u)
    = ((modulesSpecToSheaf.obj (tilde M)).presheaf.obj (op U) : Type u) := rfl

/-- ★★**制限射も同じである**(`rfl`)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで第 120 ブロックの可換図式が `IsLocallyTrivial` の側で使える。 -/
theorem tildeRes_eq (W : (Spec R).Opens) (i : W ⟶ U)
    (x : ((tilde M).val.obj (op U) : Type u)) :
    ((tilde M).val.map i.op) x
      = ((modulesSpecToSheaf.obj (tilde M)).presheaf.map i.op).hom x := rfl

/-! ## ★出典の紐付け(`.src`) -/

def tildeRes_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——2 つの前層は切断も制限も同じ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
